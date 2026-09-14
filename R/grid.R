# The grid runner (spec 3, section 4; local/ARCHITECTURE.md decisions 8, 17,
# 19).
#
# One durable file per cell, written as soon as the cell finishes; a
# resume that re-reads existing files and never depends on a temporary
# one; a smoke mode in its own directory; a preflight before the loop;
# cells ordered so the first completed block spans the design.

#' The cells of a grid, replication-major
#' @noRd
grid_cells <- function(n_rows, reps) {
  tibble::tibble(
    row = rep(seq_len(n_rows), times = reps),
    rep = rep(seq_len(reps), each = n_rows)
  )
}

#' One seed per cell, derived from the master seed
#' @noRd
cell_seed <- function(seed, row, rep) {
  if (is.null(seed)) {
    return(NA_real_)
  }
  (as.double(seed) + 7919 * row + rep) %% .Machine$integer.max
}

#' The files of one cell
#' @noRd
cell_paths <- function(dir, row, rep) {
  stem <- file.path(dir, sprintf("cell-%d-rep-%d", row, rep))
  list(sim = paste0(stem, "-sim.rds"), fit = stem)
}

#' @noRd
check_grid <- function(grid, call = rlang::caller_env()) {
  if (!is.data.frame(grid)) {
    cli::cli_abort(
      "{.arg grid} must be a data frame, not {.obj_type_friendly {grid}}.",
      call = call
    )
  }
  missing <- setdiff(c("n_subjects", "n_trials"), names(grid))
  if (length(missing) > 0L) {
    cli::cli_abort(
      "{.arg grid} is missing the column{?s} {.val {missing}}.",
      call = call
    )
  }
  if (nrow(grid) == 0L) {
    cli::cli_abort("{.arg grid} has no rows.", call = call)
  }
  invisible(grid)
}

#' Apply a row's parameter and `sd_` columns
#' @noRd
override_pars <- function(pars, row) {
  for (p in intersect(names(row), names(pars))) {
    pars[[p]] <- row[[p]]
  }
  pars
}

#' @noRd
override_sds <- function(sds, row) {
  for (col in grep("^sd_", names(row), value = TRUE)) {
    if (is.null(sds)) sds <- numeric()
    sds[sub("^sd_", "", col)] <- row[[col]]
  }
  sds
}

#' The two terms of a `cor_<a>__<b>` column
#' @noRd
cor_column_terms <- function(col) {
  strsplit(sub("^cor_", "", col), "__", fixed = TRUE)[[1L]]
}

#' Apply a row's `cor_<a>__<b>` columns to a correlation matrix
#'
#' Starting from the identity over the named terms when `cors` is `NULL`;
#' terms `cors` lacks are added as uncorrelated dimensions first.
#'
#' @noRd
override_cors <- function(cors, row) {
  cols <- grep("^cor_", names(row), value = TRUE)
  if (length(cols) == 0L) {
    return(cors)
  }
  terms <- unique(unlist(lapply(cols, cor_column_terms)))
  have <- rownames(cors)
  dims <- c(have, setdiff(terms, have))
  out <- diag(length(dims))
  dimnames(out) <- list(dims, dims)
  if (!is.null(cors)) out[have, have] <- cors
  for (col in cols) {
    ab <- cor_column_terms(col)
    out[ab[[1L]], ab[[2L]]] <- row[[col]]
    out[ab[[2L]], ab[[1L]]] <- row[[col]]
  }
  out
}

#' Check the `cor_` columns of a grid against the model and covariates
#' @noRd
check_cor_columns <- function(grid, model, covariates,
                              call = rlang::caller_env()) {
  info <- model_parameters(model)
  known <- c(info$free, info$fixed, names(covariates))
  for (col in grep("^cor_", names(grid), value = TRUE)) {
    ab <- cor_column_terms(col)
    if (length(ab) != 2L || !all(ab %in% known) || ab[[1L]] == ab[[2L]]) {
      cli::cli_abort(
        c(
          "Grid column {.val {col}} does not name two known terms as \\
           {.code cor_<a>__<b>}.",
          i = "Known: {.val {known}}."
        ),
        call = call
      )
    }
  }
  invisible(grid)
}

#' Population values, SDs and correlations for one grid row
#'
#' A grid column named after a parameter overrides `pars` for that row;
#' a column `sd_<parameter>` overrides `sds`; a column `cor_<a>__<b>` sets
#' one correlation. A `function(row)` default becomes a zero-argument
#' function applying the same overrides, so that [simulate_recovery()]
#' evaluates it under the cell seed.
#'
#' @noRd
row_values <- function(row, pars, sds, cors = NULL) {
  wrap <- function(value, override) {
    if (is.function(value)) {
      force(value)
      return(function() override(value(row), row))
    }
    override(value, row)
  }
  list(
    pars = wrap(pars, override_pars),
    sds = wrap(sds, override_sds),
    cors = wrap(cors, override_cors)
  )
}

#' Fill in the truth tables a simulation written before 5.1 lacks
#'
#' Such a file has no SD, correlation or covariate tables; its draws were
#' uncorrelated and it had no covariates, so they are rebuilt from `sds`.
#'
#' @noRd
upgrade_simulation <- function(sim) {
  if (!is.null(sim$truth$sd)) {
    return(sim)
  }
  varying <- names(sim$sds)[sim$sds > 0]
  cors <- diag(length(varying))
  dimnames(cors) <- list(varying, varying)
  sim$truth$sd <- sd_table(sim$sds)
  sim$truth$cor <- cor_table(cors)
  sim$truth$covariates <- tibble::tibble(
    id = character(), term = character(), true_value = double()
  )
  sim["cors"] <- list(if (length(varying) < 2L) NULL else cors)
  sim["covariates"] <- list(NULL)
  sim
}

#' Simulate one cell, or read it back
#'
#' With `subjects = "fixed"`, later replications take replication 1's
#' realised values as numbers, so a function-valued truth is not
#' evaluated again.
#'
#' @noRd
cell_simulation <- function(paths, model, values, row, seed, subjects,
                            first_rep, generator, covariates) {
  if (file.exists(paths$sim)) {
    return(upgrade_simulation(readRDS(paths$sim)))
  }
  subject_pars <- NULL
  if (identical(subjects, "fixed") && !is.null(first_rep)) {
    subject_pars <- dplyr::bind_rows(
      first_rep$truth$subjects, first_rep$truth$covariates
    )
    values <- list(
      pars = first_rep$pars, sds = first_rep$sds, cors = first_rep$cors
    )
  }
  sim <- simulate_recovery(
    model, values$pars,
    n_subjects = row$n_subjects, n_trials = row$n_trials,
    sds = values$sds, cors = values$cors, covariates = covariates,
    subject_pars = subject_pars,
    generator = generator,
    seed = if (is.na(seed)) NULL else seed
  )
  write_atomic(paths$sim, function(tmp) saveRDS(sim, tmp))
  sim
}

#' The preflight: one short fit of the first cell before the loop
#' @noRd
run_preflight <- function(sim, model, formula, prior, dir, dots, fitter,
                          call = rlang::caller_env()) {
  short <- utils::modifyList(dots, list(chains = 1, iter = 200))
  if (!is.null(short$warmup) && short$warmup >= 200) short$warmup <- 100
  cache_args <- c(
    list(
      formula = formula, data = sim$data, model = model,
      file = file.path(dir, "preflight"), prior = prior,
      refit = "always", .fitter = fitter
    ),
    short
  )
  started <- Sys.time()
  tryCatch(
    rlang::exec(fit_cached, !!!cache_args),
    error = function(e) {
      cli::cli_abort(
        c(
          "The preflight fit of the first cell failed; no cell was run.",
          i = "Fix the model, formula, prior or data before the grid."
        ),
        parent = e,
        call = call
      )
    }
  )
  as.double(difftime(Sys.time(), started, units = "secs"))
}

#' Fit one cell and extract its estimates
#'
#' @return A list with `status`, `message`, `estimates`, `elapsed`,
#'   `converged`.
#' @noRd
run_cell <- function(sim, model, formula, prior, paths, seed, dots, fitter) {
  cache_args <- c(
    list(
      formula = formula, data = sim$data, model = model,
      file = paths$fit, prior = prior, .fitter = fitter
    ),
    if (!is.na(seed)) list(seed = seed),
    dots
  )
  started <- Sys.time()
  result <- tryCatch(
    {
      fit <- rlang::exec(fit_cached, !!!cache_args)
      estimates <- extract_estimates(
        fit,
        level = c("population", "subject")
      )
      list(status = "ok", message = NA_character_, estimates = estimates)
    },
    error = function(e) {
      list(status = "error", message = conditionMessage(e), estimates = NULL)
    }
  )
  result$elapsed <- as.double(difftime(Sys.time(), started, units = "secs"))
  result$converged <- NA
  if (!is.null(result$estimates)) {
    result$converged <- result$estimates$converged[[1L]]
  }
  result
}

#' Score every cell that has estimates
#' @noRd
score_cells <- function(runs, sims, cells, links, scale = "natural") {
  ok <- which(vapply(runs, function(r) r$status == "ok", logical(1)))
  if (length(ok) == 0L) {
    cli::cli_abort("No cell produced a fit; nothing to score.")
  }
  label <- function(i) sprintf("row-%d", cells$row[[i]])

  estimates <- dplyr::bind_rows(lapply(ok, function(i) {
    out <- runs[[i]]$estimates
    out$replication <- cells$rep[[i]]
    out$condition <- label(i)
    out
  }))
  truth_pop <- dplyr::bind_rows(lapply(ok, function(i) {
    out <- sims[[i]]$truth$population
    out$replication <- cells$rep[[i]]
    out$condition <- label(i)
    out
  }))
  truth_sub <- dplyr::bind_rows(lapply(ok, function(i) {
    out <- sims[[i]]$truth$subjects
    out$replication <- rep(cells$rep[[i]], nrow(out))
    out$condition <- rep(label(i), nrow(out))
    out
  }))

  population <- recover(
    estimates[estimates$level == "population", ],
    truth_pop,
    scale = scale, links = links
  )
  pieces <- list(population)
  if (nrow(truth_sub) > 0L) {
    pieces$subjects <- recover_subjects(
      estimates[estimates$level == "subject", ],
      truth_sub,
      scale = scale, links = links
    )
  }
  new_bmmtools_recovery(
    dplyr::bind_rows(lapply(pieces, tibble::as_tibble)),
    scale = scale,
    ci_level = attr(population, "ci_level"),
    call = attr(population, "call")
  )
}

#' Run a parameter-recovery grid
#'
#' The loop the validation scripts in bmm wrote by hand, with the
#' requirements their overnight runs taught (decision 19): one durable
#' file per cell written as soon as the cell finishes, a resume that
#' reads those files and needs no temporary state, a smoke mode in its
#' own directory, a preflight fit that catches a compile or init error
#' before any cell runs, and cells ordered so that the first completed
#' block spans the design rather than exhausting one level.
#'
#' @param model A `bmmodel`, or a function of one grid row (a one-row
#'   data frame) returning one, for models whose constructor depends on
#'   the design.
#' @param grid A data frame with the columns `n_subjects` and `n_trials`.
#'   A column named after a parameter gives that cell's population value
#'   on the link scale, overriding `pars`; a column `sd_<parameter>`
#'   overrides `sds`; a column `cor_<a>__<b>` sets the correlation of two
#'   parameters or covariates (the names in either order). A SimDesign
#'   design is a data frame and works as is.
#' @param pars,sds,cors Defaults for every cell, as in
#'   [simulate_recovery()]. Each may also be a `function(row)` of the
#'   one-row grid data frame, evaluated under the cell's seed, which draws
#'   new hyperparameters for every data set; the grid columns above are
#'   applied to its result.
#' @param covariates As in [simulate_recovery()], the same for every cell.
#' @param dir Directory for the per-cell files; created if missing.
#' @param reps Replications per cell.
#' @param formula A `bmmformula`; `NULL` means [recovery_formula()] of the
#'   cell's model.
#' @param prior Passed to [fit_cached()].
#' @param generator As in [simulate_recovery()].
#' @param seed Master seed. Each cell derives its own from it and the
#'   cell's row and replication, so a cell is reproducible on its own.
#'   `NULL` leaves everything unseeded and records `NA`.
#' @param subjects `"redraw"` draws new subject values in every
#'   replication; `"fixed"` draws them once per row and reuses them, so
#'   replications become a simulated retest of the same people. The
#'   realised population values, SDs, correlations and covariate values of
#'   replication 1 are reused too.
#' @param re_cor Passed to [recovery_formula()] when `formula` is `NULL`.
#' @param scale The scale recovery is scored on, as in [recover()];
#'   natural by default.
#' @param smoke `TRUE` runs the first two rows with two replications
#'   into `<dir>/smoke`, so a smoke run never overwrites a full one.
#' @param preflight Run the first cell once with one chain and 200
#'   iterations into `<dir>/preflight` before the loop. Skipped when the
#'   first cell already has a cached fit.
#' @param ... Passed to [fit_cached()] and so to the fitter: `chains`,
#'   `iter`, `warmup`, `backend`, `cores`, ...
#' @param .fitter As in [fit_cached()].
#'
#' @return A `bmmtools_recovery` object over every cell that produced a
#'   fit, at both the population and the subject level, with a
#'   `condition` column (`"row-<i>"`) naming the grid row. Two
#'   attributes: `cells`, a tibble with one row per cell (`condition`,
#'   `replication`, `n_subjects`, `n_trials`, `seed`, `file`, `status`,
#'   `elapsed` in seconds, `converged`), and `grid`. `summary()` groups
#'   by condition.
#'
#' @details
#' A cell whose fit errors is recorded with `status = "error"` and the
#' grid goes on; a warning at the end names the failed cells. A cell
#' whose data cannot be generated stops the grid, because that is a
#' design error rather than a sampler accident. Cell files are
#' `cell-<row>-rep-<rep>-sim.rds` (data and truth) and
#' `cell-<row>-rep-<rep>.rds` with its `.key` (the fit, through
#' [fit_cached()]); a resume reads an existing simulation file rather
#' than regenerating, so a cached fit and its truth always match.
#'
#' @examples
#' \dontrun{
#' out <- recovery_grid(
#'   bmm::mixture2p(resp_error = "y"),
#'   grid = expand.grid(n_subjects = c(20, 50), n_trials = c(30, 100)),
#'   pars = c(kappa = log(8), thetat = qlogis(0.75)),
#'   sds = c(kappa = 0.3, thetat = 0.5),
#'   dir = "fits/mixture2p", reps = 10, seed = 1,
#'   chains = 4, iter = 2000, backend = "cmdstanr"
#' )
#' summary(out)
#' }
#'
#' @export
recovery_grid <- function(model,
                          grid,
                          pars,
                          dir,
                          reps = 1,
                          sds = NULL,
                          cors = NULL,
                          covariates = NULL,
                          formula = NULL,
                          prior = NULL,
                          generator = NULL,
                          seed = NULL,
                          subjects = c("redraw", "fixed"),
                          re_cor = c("none", "all"),
                          scale = c("natural", "link"),
                          smoke = FALSE,
                          preflight = TRUE,
                          ...,
                          .fitter = NULL) {
  subjects <- rlang::arg_match(subjects)
  re_cor <- rlang::arg_match(re_cor)
  scale <- rlang::arg_match(scale)
  check_grid(grid)
  reps <- check_count(reps, "reps")
  dots <- rlang::list2(...)
  if (!rlang::is_bool(smoke) || !rlang::is_bool(preflight)) {
    cli::cli_abort(
      "{.arg smoke} and {.arg preflight} must be {.code TRUE} or \\
       {.code FALSE}."
    )
  }
  if (isTRUE(smoke)) {
    grid <- grid[seq_len(min(2L, nrow(grid))), , drop = FALSE]
    reps <- min(reps, 2L)
    dir <- file.path(dir, "smoke")
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  models <- vector("list", nrow(grid))
  model_for <- function(i) {
    if (is.null(models[[i]])) {
      models[[i]] <<- if (is.function(model)) {
        model(grid[i, , drop = FALSE])
      } else {
        model
      }
      check_model(models[[i]])
    }
    models[[i]]
  }
  formula_for <- function(i) {
    formula %||% recovery_formula(model_for(i), re_cor = re_cor)
  }

  check_cor_columns(grid, model_for(1L), covariates)

  cells <- grid_cells(nrow(grid), reps)
  cells$seed <- cell_seed(seed, cells$row, cells$rep)
  sims <- vector("list", nrow(cells))
  runs <- vector("list", nrow(cells))
  first_rep <- vector("list", nrow(grid))

  simulate_cell <- function(i) {
    row <- grid[cells$row[[i]], , drop = FALSE]
    values <- row_values(row, pars, sds, cors)
    sim <- cell_simulation(
      cell_paths(dir, cells$row[[i]], cells$rep[[i]]),
      model_for(cells$row[[i]]), values, row, cells$seed[[i]], subjects,
      first_rep[[cells$row[[i]]]], generator, covariates
    )
    if (cells$rep[[i]] == 1L) first_rep[[cells$row[[i]]]] <<- sim
    sim
  }

  first_paths <- cell_paths(dir, 1L, 1L)
  if (isTRUE(preflight) && !file.exists(paste0(first_paths$fit, ".rds"))) {
    sims[[1L]] <- simulate_cell(1L)
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    elapsed <- run_preflight(
      sims[[1L]], model_for(1L), formula_for(1L), prior, dir, dots, .fitter
    )
    cli::cli_inform(
      "Preflight passed in {round(elapsed, 1)} s; running \\
       {nrow(cells)} cell{?s}."
    )
  }

  for (i in seq_len(nrow(cells))) {
    if (is.null(sims[[i]])) sims[[i]] <- simulate_cell(i)
    runs[[i]] <- run_cell(
      sims[[i]], model_for(cells$row[[i]]), formula_for(cells$row[[i]]),
      prior, cell_paths(dir, cells$row[[i]], cells$rep[[i]]),
      cells$seed[[i]], dots, .fitter
    )
  }

  status <- vapply(runs, function(r) r$status, character(1))
  failed <- which(status == "error")
  if (length(failed) > 0L) {
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    labels <- sprintf("row-%d rep %d", cells$row[failed], cells$rep[failed])
    cli::cli_warn(c(
      "{length(failed)} cell{?s} failed to fit: {.val {labels}}.",
      i = "The first message: {runs[[failed[[1L]]]]$message}"
    ))
  }

  # one link table scores every cell, so a row-wise model must not
  # change it; a differing table would score some cells on the wrong
  # scale without any sign of it
  links <- model_for(1L)$links
  for (i in seq_len(nrow(grid))[-1L]) {
    if (!identical(model_for(i)$links, links)) {
      cli::cli_abort(c(
        "The model of grid row {i} has different links from row 1.",
        i = "A row-wise {.arg model} must keep one link table across rows."
      ))
    }
  }
  out <- score_cells(runs, sims, cells, unlist(links), scale)
  attr(out, "cells") <- tibble::tibble(
    condition = sprintf("row-%d", cells$row),
    replication = cells$rep,
    n_subjects = as.integer(grid$n_subjects[cells$row]),
    n_trials = as.integer(grid$n_trials[cells$row]),
    seed = cells$seed,
    file = paste0(vapply(seq_len(nrow(cells)), function(i) {
      cell_paths(dir, cells$row[[i]], cells$rep[[i]])$fit
    }, character(1)), ".rds"),
    status = status,
    elapsed = vapply(runs, function(r) r$elapsed, numeric(1)),
    converged = vapply(runs, function(r) as.logical(r$converged), logical(1))
  )
  attr(out, "grid") <- grid
  out
}
