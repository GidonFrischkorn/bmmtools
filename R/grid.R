# The grid runner.
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
  list(
    sim = paste0(stem, "-sim.rds"), fit = stem,
    est = paste0(stem, "-est.rds"), ml = paste0(stem, "-ml")
  )
}

#' The cells of a grid with their seeds
#' @noRd
grid_cell_table <- function(n_rows, reps, seed) {
  cells <- grid_cells(n_rows, reps)
  cells$seed <- cell_seed(seed, cells$row, cells$rep)
  cells
}

#' @noRd
check_grid <- function(grid, required = c("n_subjects", "n_trials"),
                       call = rlang::caller_env()) {
  if (!is.data.frame(grid)) {
    cli::cli_abort(
      "{.arg grid} must be a data frame, not {.obj_type_friendly {grid}}.",
      call = call
    )
  }
  missing <- setdiff(required, names(grid))
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

#' Check `reps`, `smoke` and `preflight`, apply smoke mode, create `dir`
#'
#' Shared by both forms of recovery_grid().
#'
#' @return A list with `grid`, `reps` and `dir` as run.
#' @noRd
grid_run_setup <- function(grid, reps, dir, smoke, preflight,
                           call = rlang::caller_env()) {
  reps <- check_count(reps, "reps", call = call)
  if (!rlang::is_bool(smoke) || !rlang::is_bool(preflight)) {
    cli::cli_abort(
      "{.arg smoke} and {.arg preflight} must be {.code TRUE} or \\
       {.code FALSE}.",
      call = call
    )
  }
  if (isTRUE(smoke)) {
    grid <- grid[seq_len(min(2L, nrow(grid))), , drop = FALSE]
    reps <- min(reps, 2L)
    dir <- file.path(dir, "smoke")
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  list(grid = grid, reps = reps, dir = dir)
}

#' Warn at the end of a grid about the cells that failed to fit
#' @param labels One label per failed cell (or cell and component).
#' @param message The first failure's message.
#' @noRd
warn_failed_cells <- function(labels, message, what = "cell") {
  if (length(labels) == 0L) {
    return(invisible(NULL))
  }
  # the noun is pluralised here: a string inside cli would set the quantity
  noun <- if (length(labels) == 1L) what else paste0(what, "s")
  cli::cli_warn(c(
    paste(length(labels), noun, "failed to fit: {.val {labels}}."),
    i = "The first message: {message}"
  ))
}

#' Apply a row's parameter columns
#'
#' Without tasks a column named after an entry of `pars` replaces it. With
#' tasks a column may also name a parameter whose task terms `pars` gives,
#' or a full task term: a bare column sets the bare value and every task
#' term `pars` already has, and a full-term column is applied after it, so
#' that within a row a full term still wins, as it does in
#' `simulate_recovery()`.
#'
#' @noRd
override_pars <- function(pars, row, tasks = NULL, task_col = NULL) {
  if (is.null(tasks)) {
    for (p in intersect(names(row), names(pars))) {
      pars[[p]] <- row[[p]]
    }
    return(pars)
  }
  override_task_values(pars, row, names(row), tasks, task_col)
}

#' @noRd
override_sds <- function(sds, row, tasks = NULL, task_col = NULL) {
  cols <- grep("^sd_", names(row), value = TRUE)
  if (length(cols) == 0L) {
    return(sds)
  }
  if (is.null(sds)) sds <- numeric()
  if (is.null(tasks)) {
    for (col in cols) {
      sds[sub("^sd_", "", col)] <- row[[col]]
    }
    return(sds)
  }
  # every sd_ column names a term, known to `sds` or not, as without tasks
  override_task_values(sds, row, cols, tasks, task_col,
    strip = "^sd_", add_all = TRUE
  )
}

#' The shared step of `override_pars()` and `override_sds()` with tasks
#'
#' @param cols The row columns to consider.
#' @param strip A pattern removed from a column name to give the term.
#' @param add_all Whether a column adds its term even when `x` has neither
#'   the term nor its task terms (the `sd_` columns do; a parameter column
#'   without a matching entry is some other grid column).
#' @noRd
override_task_values <- function(x, row, cols, tasks, task_col,
                                 strip = NULL, add_all = FALSE) {
  terms <- if (is.null(strip)) cols else sub(strip, "", cols)
  suffixes <- paste0("_", task_col, tasks)
  ends_in_task <- vapply(terms, function(t) {
    any(endsWith(t, suffixes))
  }, logical(1), USE.NAMES = FALSE)
  is_bare <- !ends_in_task
  for (k in which(is_bare)) {
    own <- intersect(paste0(terms[[k]], suffixes), names(x))
    if (terms[[k]] %in% names(x) || length(own) > 0L || add_all) {
      x[[terms[[k]]]] <- row[[cols[[k]]]]
      x[own] <- row[[cols[[k]]]]
    }
  }
  for (k in which(ends_in_task)) {
    term <- terms[[k]]
    suffix <- suffixes[endsWith(term, suffixes)][[1L]]
    par <- substr(term, 1L, nchar(term) - nchar(suffix))
    known <- term %in% names(x) || par %in% names(x) ||
      any(paste0(par, suffixes) %in% names(x))
    if (known || add_all) {
      x[[term]] <- row[[cols[[k]]]]
    }
  }
  x
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
check_cor_columns <- function(grid, model, covariates, tasks = NULL,
                              task_col = NULL, call = rlang::caller_env()) {
  info <- model_parameters(model)
  known <- c(
    task_terms(info$free, tasks, task_col), info$fixed, names(covariates)
  )
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
#' one correlation. With tasks, parameter and `sd_` columns may use full
#' task terms too. A `function(row)` default becomes a zero-argument
#' function applying the same overrides, so that [simulate_recovery()]
#' evaluates it under the cell seed.
#'
#' @noRd
row_values <- function(row, pars, sds, cors = NULL, tasks = NULL,
                       task_col = NULL) {
  wrap <- function(value, override) {
    if (is.function(value)) {
      force(value)
      return(function() override(value(row), row))
    }
    override(value, row)
  }
  list(
    pars = wrap(pars, function(x, r) override_pars(x, r, tasks, task_col)),
    sds = wrap(sds, function(x, r) override_sds(x, r, tasks, task_col)),
    cors = wrap(cors, override_cors)
  )
}

#' Fill in what a simulation written by an earlier version lacks
#'
#' A file from before 5.1 has no SD, correlation or covariate tables; its
#' draws were uncorrelated and it had no covariates, so they are rebuilt
#' from `sds`. A file from before 5.4 has no tasks.
#'
#' @noRd
upgrade_simulation <- function(sim) {
  for (field in setdiff(c("tasks", "task_col"), names(sim))) {
    sim[field] <- list(NULL)
  }
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
                            first_rep, generator, covariates, tasks = NULL,
                            task_col = "task") {
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
    tasks = tasks, task_col = task_col,
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
                          file = "preflight", call = rlang::caller_env()) {
  short <- utils::modifyList(dots, list(chains = 1, iter = 200))
  if (!is.null(short$warmup) && short$warmup >= 200) short$warmup <- 100
  cache_args <- c(
    list(
      formula = formula, data = sim$data, model = model,
      file = file.path(dir, file), prior = prior,
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

#' The version string a sidecar records
#'
#' A function of its own so that a test can pretend to be another version.
#'
#' @noRd
bmmtools_version <- function() {
  as.character(utils::packageVersion("bmmtools"))
}

#' What the grid extracts from every fit
#'
#' `ml` is `NULL`, or the record [ml_request()] makes of the `fit_ml()`
#' arguments: what a sidecar has to match before its ML rows are reused.
#'
#' @noRd
extraction_request <- function(levels = c("population", "subject"),
                               correlations = NULL,
                               cor_scale = "link",
                               links = NULL,
                               ml = NULL) {
  list(
    levels = levels, correlations = correlations, cor_scale = cor_scale,
    links = links, ml = ml
  )
}

#' The `fit_ml()` arguments the grid owns and refuses in `ml`
#' @noRd
grid_owned_ml_args <- function() {
  c("model", "data", "file", "refit", "seed", "by")
}

#' Validate `ml` of recovery_grid() into a list of `fit_ml()` arguments
#'
#' `FALSE` is no ML fit; `TRUE` is [fit_ml()] with its defaults; a named
#' list gives arguments to it, and a name [fit_ml()] has no formal for
#' goes through its `...` to the fitter (`refresh`, `silent`, ...). The
#' grid supplies `model`, `data`, `file`, `refit` and `seed` itself and
#' its subjects are always `id`, so those are refused; `method` and
#' `prior` are checked here, before any cell runs, because [fit_ml()]
#' would otherwise refuse them once per cell.
#'
#' @return `NULL`, or a named list of arguments for [fit_ml()].
#' @noRd
check_ml_arg <- function(ml, levels, call = rlang::caller_env()) {
  if (isFALSE(ml)) {
    return(NULL)
  }
  if (isTRUE(ml)) {
    ml <- list()
  }
  is_arg_list <- is.list(ml) && !is.data.frame(ml) &&
    (length(ml) == 0L || !is.null(names(ml)) && all(nzchar(names(ml))))
  if (!is_arg_list) {
    cli::cli_abort(
      "{.arg ml} must be {.code TRUE}, {.code FALSE} or a named list of \\
       arguments for {.fn fit_ml}, not {.obj_type_friendly {ml}}.",
      call = call
    )
  }
  owned <- intersect(names(ml), grid_owned_ml_args())
  if (length(owned) > 0L) {
    cli::cli_abort(
      c(
        "{.arg ml} cannot set {.arg {owned}}; the grid supplies \\
         {?it/them}.",
        i = "Each cell is fitted on its own simulated data, cached next \\
             to the cell's files, under the cell's seed, by {.val id}."
      ),
      call = call
    )
  }
  if (!is.null(ml$method)) {
    ml$method <- rlang::arg_match0(
      ml$method, c("stan", "optim"),
      error_call = call
    )
  }
  if (!is.null(ml$prior)) {
    ml$prior <- rlang::arg_match0(
      ml$prior, c("flat", "default"),
      error_call = call
    )
  }
  # each value is legal on its own; what the grid must catch is a pair
  # that is not, because fit_ml() refuses it once per cell --- after every
  # hierarchical fit of the design has already run
  wrong_route <- if (identical(ml$method, "optim")) {
    c(
      if (identical(ml$prior, "default")) "prior",
      intersect(names(ml), c("formula", ".fitter", "draws"))
    )
  } else {
    intersect(names(ml), "nll")
  }
  if (length(wrong_route) > 0L) {
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    other <- if (identical(ml$method, "optim")) "stan" else "optim"
    cli::cli_abort(
      c(
        "{.arg {wrong_route}} {?belongs/belong} to \\
         {.code method = \"{other}\"}, not to \\
         {.code method = \"{ml$method %||% \"stan\"}\"}.",
        i = "{cli::qty(length(wrong_route))}{.fn fit_ml} would refuse \\
             {?it/them} once per cell, after every hierarchical fit had \\
             run."
      ),
      call = call
    )
  }
  # nothing consumes the dots on the optim route, so an unknown name there
  # is a typo that would reach nothing and say nothing
  if (identical(ml$method, "optim")) {
    unknown <- setdiff(names(ml), names(formals(fit_ml)))
    if (length(unknown) > 0L) {
      cli::cli_abort(
        c(
          "{.arg ml} has no argument{?s} {.val {unknown}} for \\
           {.code method = \"optim\"}.",
          i = "Unknown names reach the Stan optimiser through \\
               {.fn fit_ml}'s dots; the optim route calls no fitter."
        ),
        call = call
      )
    }
  }
  if (!"subject" %in% levels) {
    cli::cli_abort(
      c(
        "{.arg ml} needs the {.val subject} level, but {.arg levels} \\
         does not include it.",
        i = "A subject-wise ML fit has estimates at the subject level \\
             only."
      ),
      call = call
    )
  }
  ml
}

#' The record of an ML request a sidecar stores and has to match
#'
#' Every argument that changes the ML rows, at its `fit_ml()` default
#' when not given, so that `ml = TRUE` and `ml = list(prior = "flat")`
#' make the same record. The fit itself is cached by its own key; this
#' record also covers what is applied after the fit (`ci_level`,
#' `max_abs_link`), so a change there re-extracts from the cached fit.
#' A formula enters deparsed, as in the cache key. `.fitter` is not part
#' of what the rows are.
#'
#' @return `NULL` when `ml` is, else a named list of plain values.
#' @noRd
ml_request <- function(ml) {
  if (is.null(ml)) {
    return(NULL)
  }
  record <- list(
    formula = NULL, method = "stan", prior = "flat", ci_level = 0.95,
    draws = 1000, max_abs_link = 20, start = NULL, nll = NULL
  )
  given <- intersect(names(ml), names(record))
  # `[<-` with a list keeps a NULL element where `$<-` would drop it
  record[given] <- ml[given]
  record["formula"] <- list(formula_key(record$formula))
  # as for `init` in the cache key: a closure differs in every environment
  # it is built in, so what enters is its text
  record["nll"] <- list(function_key(record$nll))
  record
}

#' Validate `levels`, `correlations` and `cor_scale` of recovery_grid()
#' @noRd
check_extraction_args <- function(levels, correlations, cor_scale,
                                  call = rlang::caller_env()) {
  if (!is.character(levels) || length(levels) == 0L) {
    cli::cli_abort(
      "{.arg levels} must be one or more of {.val {c(\"population\", \\
       \"subject\", \"sd\")}}.",
      call = call
    )
  }
  levels <- unique(rlang::arg_match(
    levels, c("population", "subject", "sd"),
    multiple = TRUE, error_call = call
  ))
  if (!is.null(correlations)) {
    bad <- !is.character(correlations) || length(correlations) == 0L ||
      anyNA(correlations)
    if (bad) {
      cli::cli_abort(
        "{.arg correlations} must be {.code NULL} or one or more of \\
         {.val {c(\"model\", \"draws\", \"point\")}}.",
        call = call
      )
    }
    correlations <- unique(rlang::arg_match(
      correlations, c("model", "draws", "point"),
      multiple = TRUE, error_call = call
    ))
  }
  cor_scale <- rlang::arg_match(
    cor_scale, c("link", "natural"),
    error_call = call
  )
  if ("model" %in% correlations && identical(cor_scale, "natural")) {
    cli::cli_abort(
      c(
        "The {.val model} estimator is on the link scale only.",
        i = "Use {.code cor_scale = \"link\"}, or leave {.val model} out of \\
             {.arg correlations}."
      ),
      call = call
    )
  }
  list(levels = levels, correlations = correlations, cor_scale = cor_scale)
}

#' Read a cell's sidecar when it matches the cell and the request
#'
#' It matches when its key is the cell's cache key, it was written by this
#' version of bmmtools, it holds every requested level and estimator
#' (and, for correlations, the requested scale), and, with an ML fit
#' requested, its ML record is the request's. What it holds beyond the
#' request is dropped, so that the result is what a fresh extraction
#' would give. An unreadable file counts as no sidecar.
#'
#' The set-level sidecar of the component form stores its keys, one per
#' component, under `keys` and no estimates; `key_field` names the field.
#'
#' @return The trimmed sidecar, or `NULL`.
#' @noRd
read_sidecar <- function(path, key, request, key_field = "key") {
  if (!file.exists(path)) {
    return(NULL)
  }
  stored <- tryCatch(readRDS(path), error = function(e) NULL)
  scale_ok <- is.null(request$correlations) ||
    identical(stored$cor_scale, request$cor_scale)
  matches <- is.list(stored) &&
    identical(stored[[key_field]], key) &&
    identical(stored$bmmtools_version, bmmtools_version()) &&
    all(request$levels %in% stored$levels) &&
    all(request$correlations %in% stored$correlations) &&
    scale_ok
  if (!matches) {
    return(NULL)
  }
  estimates <- stored$estimates
  if (!is.null(estimates)) {
    stored$estimates <- estimates[estimates$level %in% request$levels, ]
  }
  # `[["ml"]]`, not `$ml`: `$` on a list matches partially, so a sidecar
  # holding `ml_estimates` and no `ml` would answer with the rows
  if (is.null(request$ml)) {
    stored["ml_estimates"] <- list(NULL)
  } else if (!identical(stored[["ml"]], request$ml)) {
    # the ML record is stale or absent while the hierarchical extraction
    # is current. Dropping the whole sidecar would refit the hierarchical
    # model too --- and once the fit files have been deleted, as
    # ?recovery_grid says they may be, that is the whole design resampled
    # to change `max_abs_link`. The caller runs ML alone.
    stored["ml"] <- list(NULL)
    stored["ml_estimates"] <- list(NULL)
  }
  if (is.null(request$correlations)) {
    stored["cor_estimates"] <- list(NULL)
  } else if (!is.null(stored$cor_estimates)) {
    cors <- stored$cor_estimates
    cors <- cors[cors$estimator %in% request$correlations, ]
    stored$cor_estimates <- cors[
      order(match(cors$estimator, request$correlations)), ,
      drop = FALSE
    ]
  }
  stored
}

#' Write a sidecar: its key(s), what it was extracted with, the extraction
#'
#' @param keys A named list, `list(key = )` for a cell or a component,
#'   `list(keys = )` for the set-level sidecar.
#' @param levels Whether to record the requested levels (the set-level
#'   sidecar holds no estimates and records none).
#' @noRd
write_sidecar <- function(path, keys, request, extracted, levels = TRUE) {
  sidecar <- c(
    keys,
    list(bmmtools_version = bmmtools_version()),
    if (levels) list(levels = request$levels),
    list(
      correlations = request$correlations,
      cor_scale = request$cor_scale
    ),
    # only when ML rows are stored, so a sidecar without them keeps the
    # shape it had before there was an ML fit
    if (!is.null(request$ml)) list(ml = request$ml),
    extracted
  )
  write_atomic(path, function(tmp) saveRDS(sidecar, tmp))
}

#' Everything the grid keeps from one fit, while it is in memory
#' @noRd
extract_cell <- function(fit, sim, request) {
  estimates <- extract_estimates(fit, level = request$levels)
  converged <- NULL
  if (nrow(estimates) > 0L) converged <- as.logical(estimates$converged[[1L]])

  cor_estimates <- NULL
  if (!is.null(request$correlations)) {
    cov_names <- names(sim$covariates)
    covariates <- NULL
    if (length(cov_names) > 0L) {
      covariates <- as.data.frame(sim$data)[c("id", cov_names)]
    }
    cor_estimates <- extract_correlations(
      fit,
      estimator = request$correlations, covariates = covariates,
      scale = request$cor_scale, links = request$links,
      converged = converged
    )
  }

  subject_means <- empty_subject_means()
  if (nrow(sim$truth$subjects) > 0L) {
    subject_means <- subject_means_from_draws(
      extract_subject_draws(fit), sim$truth$covariates
    )
  }
  list(
    estimates = estimates, cor_estimates = cor_estimates,
    subject_means = subject_means
  )
}

#' The subject-wise ML fit of one cell, on the same simulated data
#'
#' One [fit_ml()] call per cell, cached at `path` under the cell's seed,
#' so a resume reuses the ML fit as it reuses the hierarchical one. A
#' failure is the cell's ML status, not the grid's: the hierarchical fit
#' of the cell is still scored.
#'
#' @return A list with `status`, `message`, `estimates` (a `bmmtools_ml`,
#'   or `NULL`) and `elapsed`.
#' @noRd
run_cell_ml <- function(sim, model, path, seed, ml) {
  args <- c(
    list(model = model, data = sim$data, file = path, refit = "on_change"),
    if (!is.na(seed)) list(seed = seed),
    ml
  )
  started <- Sys.time()
  result <- tryCatch(
    list(
      status = "ok", message = NA_character_,
      estimates = rlang::exec(fit_ml, !!!args)
    ),
    error = function(e) {
      list(status = "error", message = conditionMessage(e), estimates = NULL)
    }
  )
  result$elapsed <- as.double(difftime(Sys.time(), started, units = "secs"))
  result
}

#' Fit one cell and extract its estimates, or read them from its sidecar
#'
#' With `ml` (a list of [fit_ml()] arguments), the cell's ML fit follows
#' the hierarchical one and its rows go into the sidecar with the ML
#' record, so that a resume reads both or neither. When the ML fit fails
#' the sidecar is written without it and a resume tries the ML fit again.
#'
#' @return A list with `status`, `message`, `estimates`, `cor_estimates`,
#'   `subject_means`, `elapsed`, `converged`, and with `ml` also
#'   `ml_estimates`, `ml_status`, `ml_message` and `ml_elapsed`.
#' @noRd
run_cell <- function(sim, model, formula, prior, paths, seed, dots, fitter,
                     request = extraction_request(), ml = NULL) {
  fit_dots <- c(if (!is.na(seed)) list(seed = seed), dots)
  cache_args <- c(
    list(
      formula = formula, data = sim$data, model = model,
      file = paths$fit, prior = prior, .fitter = fitter
    ),
    fit_dots
  )
  # a cell whose hierarchical fit fails has no ML status: it is not scored
  ml_run <- list(
    status = NA_character_, message = NA_character_, estimates = NULL,
    elapsed = NA_real_
  )
  started <- Sys.time()
  result <- tryCatch(
    {
      key <- cache_key(formula, sim$data, model, prior, fit_dots)$key
      extracted <- read_sidecar(paths$est, key, request)
      if (is.null(extracted)) {
        fit <- rlang::exec(fit_cached, !!!cache_args)
        extracted <- extract_cell(fit, sim, request)
        written <- request
        if (!is.null(ml)) {
          ml_run <- run_cell_ml(sim, model, paths$ml, seed, ml)
          extracted["ml_estimates"] <- list(ml_run$estimates)
          if (!identical(ml_run$status, "ok")) written["ml"] <- list(NULL)
        }
        write_sidecar(paths$est, list(key = key), written, extracted)
      } else if (!is.null(ml) && is.null(extracted$ml_estimates)) {
        # the hierarchical rows are current and the ML ones are not: run
        # the ML fit alone and rewrite the sidecar around the rows already
        # there, rather than refitting a model nothing asked about
        ml_run <- run_cell_ml(sim, model, paths$ml, seed, ml)
        extracted["ml_estimates"] <- list(ml_run$estimates)
        written <- request
        if (!identical(ml_run$status, "ok")) written["ml"] <- list(NULL)
        write_sidecar(
          paths$est, list(key = key), written,
          extracted[c(
            "estimates", "cor_estimates", "subject_means", "ml_estimates"
          )]
        )
      } else if (!is.null(ml)) {
        # the sidecar matched the ML record, so its ML rows are current
        ml_run <- list(
          status = "ok", message = NA_character_,
          estimates = extracted$ml_estimates, elapsed = 0
        )
      }
      list(
        status = "ok", message = NA_character_,
        estimates = extracted$estimates,
        cor_estimates = extracted$cor_estimates,
        subject_means = extracted$subject_means
      )
    },
    error = function(e) {
      list(status = "error", message = conditionMessage(e), estimates = NULL)
    }
  )
  result$elapsed <- as.double(difftime(Sys.time(), started, units = "secs"))
  result$converged <- NA
  if (!is.null(result$estimates) && nrow(result$estimates) > 0L) {
    result$converged <- result$estimates$converged[[1L]]
  }
  if (!is.null(ml)) {
    result["ml_estimates"] <- list(ml_run$estimates)
    result$ml_status <- ml_run$status
    result$ml_message <- ml_run$message
    result$ml_elapsed <- ml_run$elapsed
  }
  result
}

#' Bind one table of every scored cell, labelled with its cell
#' @noRd
bind_cells <- function(ok, cells, table_of) {
  dplyr::bind_rows(lapply(ok, function(i) {
    out <- table_of(i)
    out$replication <- rep(cells$rep[[i]], nrow(out))
    out$condition <- rep(sprintf("row-%d", cells$row[[i]]), nrow(out))
    out
  }))
}

#' Score the correlations of every cell that has estimates
#' @noRd
score_cell_correlations <- function(runs, sims, cells, ok, request, links) {
  if (is.null(request$correlations)) {
    return(NULL)
  }
  estimates <- bind_cells(ok, cells, function(i) runs[[i]]$cor_estimates)
  truth <- lapply(
    stats::setNames(nm = c("cor", "subjects", "covariates")),
    function(name) bind_cells(ok, cells, function(i) sims[[i]]$truth[[name]])
  )
  if (nrow(estimates) == 0L || nrow(truth$cor) == 0L) {
    cli::cli_warn(
      "Correlations were requested, but no cell has a correlation pair to \\
       score."
    )
    return(NULL)
  }
  recover_correlations(
    estimates, truth,
    estimator = request$correlations, scale = request$cor_scale,
    links = links
  )
}

#' Score every cell that has estimates
#'
#' @return A `bmmtools_recovery` with the attributes `correlations` (a
#'   `bmmtools_cor_recovery`, or absent) and `subject_means`.
#' @noRd
score_cells <- function(runs, sims, cells, links, scale = "natural",
                        request = extraction_request()) {
  ok <- which(vapply(runs, function(r) r$status == "ok", logical(1)))
  if (length(ok) == 0L) {
    cli::cli_abort("No cell produced a fit; nothing to score.")
  }

  estimates <- bind_cells(ok, cells, function(i) runs[[i]]$estimates)
  truth_of <- function(name) {
    bind_cells(ok, cells, function(i) sims[[i]]$truth[[name]])
  }
  # the names truth_pop and truth_sub reach the `call` attribute
  truth_pop <- truth_of("population")
  truth_sub <- truth_of("subjects")

  pieces <- list()
  if ("population" %in% request$levels) {
    pieces$population <- recover(
      estimates[estimates$level == "population", ],
      truth_pop,
      scale = scale, links = links
    )
  }
  if ("subject" %in% request$levels && nrow(truth_sub) > 0L) {
    subjects <- estimates[estimates$level == "subject", ]
    if (!is.null(request$ml)) {
      # the ML rows are subject rows and nothing else: an ML fit has no
      # population or sd level. They are bound here, before scoring, so
      # that both estimators meet one truth and summary() groups them.
      # The natural scale needs `links`, which the ML rows do not carry.
      has_ml <- ok[vapply(ok, function(i) {
        !is.null(runs[[i]]$ml_estimates)
      }, logical(1))]
      ml_rows <- bind_cells(has_ml, cells, function(i) {
        tibble::as_tibble(runs[[i]]$ml_estimates)
      })
      subjects <- dplyr::bind_rows(subjects, ml_rows)
    }
    pieces$subjects <- recover_subjects(
      subjects, truth_sub,
      scale = scale, links = links
    )
  } else if (!is.null(request$ml)) {
    cli::cli_abort(c(
      "{.arg ml} was requested, but there is no subject level to score \\
       it at.",
      i = "The ML rows are subject rows; nothing varies between subjects \\
           in this grid."
    ))
  }
  truth_sd <- truth_of("sd")
  if ("sd" %in% request$levels && nrow(truth_sd) > 0L) {
    # SDs are scored on the link scale whatever `scale` says
    pieces$sd <- recover(
      estimates[estimates$level == "sd", ],
      truth_sd,
      level = "sd", scale = "link"
    )
  }
  if (length(pieces) == 0L) {
    cli::cli_abort(
      "Nothing to score at the level{?s} {.val {request$levels}}."
    )
  }
  out <- new_bmmtools_recovery(
    dplyr::bind_rows(lapply(pieces, tibble::as_tibble)),
    scale = scale,
    ci_level = attr(pieces[[1L]], "ci_level"),
    call = attr(pieces[[1L]], "call")
  )

  attr(out, "correlations") <- score_cell_correlations(
    runs, sims, cells, ok, request, links
  )
  subject_means <- dplyr::bind_rows(
    label_subject_means(
      empty_subject_means(), list(), character(), character(), integer()
    ),
    lapply(ok, function(i) {
      label_subject_means(
        runs[[i]]$subject_means, sims[[i]]$truth, names(sims[[i]]$covariates),
        condition = sprintf("row-%d", cells$row[[i]]),
        replication = cells$rep[[i]]
      )
    })
  )
  attr(subject_means, "links") <- links
  attr(out, "subject_means") <- subject_means
  out
}

#' Refuse the model estimator when the default formula correlates nothing
#'
#' Without it the grid would fail to find a model correlation only after
#' every cell ran. `"within"` without tasks falls back to `"none"` in
#' [recovery_formula()], so it is refused too.
#'
#' @noRd
check_model_correlations <- function(correlations, formula, re_cor, tasks,
                                     call = rlang::caller_env()) {
  uncorrelated <- identical(re_cor, "none") ||
    (identical(re_cor, "within") && is.null(tasks))
  if ("model" %in% correlations && is.null(formula) && uncorrelated) {
    hint <- if (is.null(tasks)) {
      "{.code re_cor = \"all\"}"
    } else {
      "{.code re_cor = \"within\"} or {.code re_cor = \"all\"}"
    }
    cli::cli_abort(
      c(
        "{.code correlations = \"model\"} needs correlated random effects, \\
         but the default formula has none.",
        i = paste0(
          "Use ", hint, ", a {.arg formula} with correlated terms such as \\
           {.code (1 | p | id)}, or the {.val draws} and {.val point} \\
           estimators."
        )
      ),
      call = call
    )
  }
  invisible(NULL)
}

#' The model of grid row 1, evaluated when `model` is a function
#'
#' Evaluated once, before anything else, because its form decides which
#' form of the grid runs. A grid that is not a data frame with rows is left
#' to `check_grid()`.
#'
#' @noRd
grid_first_model <- function(model, grid) {
  if (!is.function(model) || !is.data.frame(grid) || nrow(grid) == 0L) {
    return(model)
  }
  model(grid[1L, , drop = FALSE])
}

#' Is `x` a non-empty plain list of components?
#' @noRd
is_component_list <- function(x) {
  is.list(x) && !is.object(x) && length(x) > 0L &&
    all(vapply(x, inherits, logical(1), what = "bmmtools_component"))
}

#' The formula of one grid row
#'
#' `formula` is `NULL` (the default formula of the row's model), a
#' `bmmformula`, or a `function(row)` returning one.
#'
#' @noRd
grid_formula <- function(formula, row, i, model, re_cor, task_col,
                         call = rlang::caller_env()) {
  if (is.null(formula)) {
    return(recovery_formula(model, re_cor = re_cor, task_col = task_col))
  }
  if (!is.function(formula)) {
    return(formula)
  }
  out <- formula(row)
  if (!inherits(out, "bmmformula")) {
    cli::cli_abort(
      "{.arg formula} returned {.obj_type_friendly {out}} for grid row {i}; \\
       it must return a {.cls bmmformula}.",
      call = call
    )
  }
  out
}

#' Run a parameter-recovery grid
#'
#' The loop the validation scripts in bmm wrote by hand, with the
#' requirements their overnight runs taught: one durable
#' file per cell written as soon as the cell finishes, a resume that
#' reads those files and needs no temporary state, a smoke mode in its
#' own directory, a preflight fit that catches a compile or init error
#' before any cell runs, and cells ordered so that the first completed
#' block spans the design rather than exhausting one level.
#'
#' @param model A `bmmodel`, or a function of one grid row (a one-row
#'   data frame) returning one, for models whose constructor depends on
#'   the design. Or a list of [recovery_component()]s, or a function of the
#'   row returning such a list, to simulate several models for the same
#'   people and fit each on its own; see the section "Components".
#' @param grid A data frame with the columns `n_subjects` and `n_trials`
#'   (with components, `n_subjects` only; see "Components" for their
#'   columns).
#'   A column named after a parameter gives that cell's population value
#'   on the link scale, overriding `pars`; a column `sd_<parameter>`
#'   overrides `sds`; a column `cor_<a>__<b>` sets the correlation of two
#'   parameters or covariates (the names in either order). With `tasks`,
#'   these columns may also use full task terms (`kappa_task2`,
#'   `sd_kappa_task2`, `cor_kappa_task1__kappa_task2`); a bare parameter
#'   column sets every task, and a full-term column in the same row
#'   overrides it for its task. A SimDesign design is a data frame and
#'   works as is.
#' @param pars,sds,cors Defaults for every cell, as in
#'   [simulate_recovery()]. Each may also be a `function(row)` of the
#'   one-row grid data frame, evaluated under the cell's seed, which draws
#'   new hyperparameters for every data set; the grid columns above are
#'   applied to its result.
#' @param covariates As in [simulate_recovery()], the same for every cell.
#' @param tasks,task_col As in [simulate_recovery()], the same for every
#'   cell. With `tasks`, the default formula is
#'   `recovery_formula(model, re_cor = re_cor, task_col = task_col)`.
#' @param dir Directory for the per-cell files; created if missing.
#' @param reps Replications per cell.
#' @param formula A `bmmformula`; `NULL` means [recovery_formula()] of the
#'   cell's model; or a `function(row)` of the one-row grid data frame
#'   returning a `bmmformula`, called once per row, for designs whose
#'   formula depends on the row.
#' @param prior Passed to [fit_cached()]; with components, `NULL` or a list
#'   with a prior per component.
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
#' @param levels The estimate levels extracted from every fit and scored:
#'   one or more of `"population"`, `"subject"` and `"sd"` (the
#'   between-subject standard deviations, always scored on the link
#'   scale).
#' @param correlations `NULL` (the default) for no correlations, or one or
#'   more of the estimators of [extract_correlations()], `"model"`,
#'   `"draws"` and `"point"`, to extract and score the between-subject
#'   correlations of every cell, pairs with covariates included.
#' @param cor_scale The scale the correlations are extracted and scored
#'   on, `"link"` or `"natural"`. The `"model"` estimator exists on the
#'   link scale only.
#' @param ml `FALSE`, the default, fits every cell hierarchically only.
#'   `TRUE` also fits every cell subject by subject with [fit_ml()] on
#'   the same simulated data, and scores both estimators against one
#'   truth at the subject level; the ML rows carry `estimator = "ml"`. A
#'   named list gives arguments to [fit_ml()] (`prior`, `draws`,
#'   `max_abs_link`, `ci_level`, `formula`, `.fitter`, and through its
#'   `...` the fitter's own, such as `refresh = 0, silent = 2` to quiet
#'   the optimiser); `model`, `data`, `file`, `refit`, `seed` and `by`
#'   are the grid's. Needs `"subject"` in `levels`. See the section
#'   "Subject-wise ML".
#' @param smoke `TRUE` runs the first two rows with two replications
#'   into `<dir>/smoke`, so a smoke run never overwrites a full one.
#' @param preflight Run the first cell once with one chain and 200
#'   iterations into `<dir>/preflight` before the loop. Skipped when the
#'   first cell already has a cached fit or a sidecar that can be used.
#' @param ... Passed to [fit_cached()] and so to the fitter: `chains`,
#'   `iter`, `warmup`, `backend`, `cores`, ...
#' @param .fitter As in [fit_cached()].
#'
#' @return A `bmmtools_recovery` object over every cell that produced a
#'   fit, at the requested `levels`, with a `condition` column
#'   (`"row-<i>"`) naming the grid row. `summary()` groups by condition.
#'   Its attributes:
#'   * `cells`: a tibble with one row per cell (`condition`,
#'     `replication`, `n_subjects`, `n_trials`, `seed`, `file`, `status`,
#'     `elapsed` in seconds, `converged`).
#'   * `grid`: the grid as run.
#'   * `correlations`: with `correlations` requested, a
#'     `bmmtools_cor_recovery` from [recover_correlations()] over every
#'     cell, with `condition` and `replication`; otherwise `NULL`.
#'   * `subject_means`: a tibble with one row per cell, subject and term:
#'     `condition`, `replication`, `id`, `term`, `covariate` (whether the
#'     term is a covariate), the posterior `mean` and `median` on the link
#'     scale, and `true_value`, the generating value on the link scale (a
#'     covariate's value is its own mean, median and true value). Its
#'     attribute `links` is the model's link table. [subject_table()]
#'     turns it into one row per subject.
#'   * `ml_cells`: with `ml`, a tibble with one row per cell (`condition`,
#'     `replication`, `status`, `elapsed`, `n_subjects`, `n_converged`,
#'     `converged`, `file`) for the ML fit; see "Subject-wise ML".
#'     `cells` stays one row per cell and describes the hierarchical fit.
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
#' **The sidecar.** While a fit is in memory, everything scored from it is
#' written to `cell-<row>-rep-<rep>-est.rds`: the estimates at `levels`,
#' the correlations, the subject means, and what they were extracted
#' with (the cell's cache key, the bmmtools version, `levels`,
#' `correlations`, `cor_scale`). A resume uses the sidecar without
#' reading the fit when the key and the version match and it holds every
#' requested level and estimator, so the fit files may be deleted once a
#' grid has run. Otherwise the cell goes through [fit_cached()] again,
#' which reuses a cached fit, and the sidecar is rewritten.
#'
#' @section Subject-wise ML:
#' With `ml`, every cell is fitted twice on the same simulated data: the
#' hierarchical fit, and one [fit_ml()] call with no pooling, cached as
#' `cell-<row>-rep-<rep>-ml.rds` under the cell's seed. The ML rows are
#' bound to the hierarchical subject rows before scoring, so the result
#' has both estimators at the subject level, `"bayes"` and `"ml"`,
#' against one truth; `summary()` reports them as separate rows and
#' `plot_recovery(color_by = "estimator")` colours them. An ML fit has no
#' population or `sd` level, so those levels stay hierarchical only, and
#' its subject means are not added to `subject_means`. The `...` of the
#' grid (`chains`, `iter`, `backend`, `refresh`, ...) belong to the
#' hierarchical fit and do not reach [fit_ml()]; `ml = list(draws = ,
#' refresh = 0)` and the other [fit_ml()] arguments do.
#'
#' The sidecar records the ML request (its prior, draws, interval level,
#' link-scale range and formula) with the ML rows, and a resume reuses
#' them only when the request is the same; a change to `max_abs_link` or
#' `ci_level` re-extracts from the cached ML fit without refitting, and a
#' change to `prior`, `draws` or `formula` refits through [fit_cached()].
#' A cell whose ML fit fails keeps its hierarchical rows, is recorded
#' with `status = "error"` in `ml_cells`, and is named in a warning at
#' the end; its sidecar is written without ML rows, so a resume tries the
#' ML fit again. There is no ML preflight: a compile error of the
#' no-pooling model shows up as every cell's ML status.
#'
#' Within a cell, a subject whose ML estimate leaves the link-scale range
#' keeps its row with `estimate = NA` and `converged = FALSE`, as
#' [fit_ml()] documents; the recovery rows carry that verdict per subject
#' and term, and `ml_cells` counts the subjects per cell. Report `n` next
#' to any comparison of the two estimators.
#'
#' @section Components:
#' With `model` a list of [recovery_component()]s (or a `function(row)`
#' returning one, the same components by name in every row), each cell
#' simulates one set with [simulate_components()] under the cell's seed
#' and fits every component on its own through [fit_cached()]. `pars`,
#' `sds`, `generator`, `tasks`, `task_col`, `re_cor` and `formula` belong to
#' the components and are errors here; `cors`, `covariates`, `seed`,
#' `reps`, `scale`, `levels`, `correlations`, `cor_scale`, `smoke`,
#' `preflight` and `...` keep their meaning, and `prior` is `NULL` or a
#' list by component, as in [fit_components()]. `subjects = "fixed"` is not
#' supported for components yet. A component cannot be named `sim`, `est`,
#' `cor` or `sd`.
#'
#' Grid columns: `n_subjects`; `n_trials_<comp>` for a component's number
#' of trials (a plain `n_trials` column is an error); `<comp>_<par>` or
#' `<comp>_<par>_<task_col><level>` for population values, a bare parameter
#' setting every task and a full term overriding it; `sd_<comp>_<par>` for
#' SDs; `cor_<a>__<b>` with prefixed terms or covariates, in either order.
#' They apply to numeric `pars` and `sds` of a component and to the result
#' of a function-valued one. A missing value (`NA`) leaves the component's
#' value, so a column may name a task that only some rows have. A column
#' whose component or term is unknown is an error; so is a column whose
#' prefix is no component but whose remainder is a term of one, such as a
#' misspelt component name. Other columns are left for `model` and `cors`
#' functions.
#'
#' Files per cell: `cell-<row>-rep-<rep>-sim.rds` holds the set;
#' `cell-<row>-rep-<rep>-<comp>.rds` (with its `.key`) each fit; a sidecar
#' `cell-<row>-rep-<rep>-<comp>-est.rds` per component holds its estimates
#' and, with `"model"` requested, its model correlations; and
#' `cell-<row>-rep-<rep>-cor.rds` holds what needs every fit, the
#' correlations of [extract_correlations()] on the set of fits (covariates
#' included) and the subject means, with the components' keys. A resume in
#' which every sidecar matches reads no fit; otherwise the fits are obtained
#' through [fit_cached()], which reuses cached ones, and the stale sidecars
#' are rewritten. The preflight fits every component of the first cell
#' into `<dir>/preflight-<comp>`, and is skipped when each of them has a
#' cached fit or a usable sidecar.
#'
#' The result is one `bmmtools_recovery` with prefixed terms (`a_kappa`,
#' `b_c_task1`); `converged` is that component fit's verdict, and
#' `scale = "natural"` uses every component's links. The `cells` attribute
#' has one row per cell and component, with a `component` column; the
#' `correlations` and `subject_means` attributes use prefixed terms, so
#' [subject_table()] works on the result. A component whose fit fails is
#' recorded as `"error"` for that cell, the cell's other components are
#' still scored, its correlations and subject means are skipped, and the
#' grid goes on.
#'
#' Correlations between parameters of separately fitted components are
#' attenuated by the reliabilities of both estimates; see "Separate fits"
#' in [extract_correlations()]. The `"model"` estimator works within one
#' fit only and needs a component `formula` with correlated terms.
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
                          tasks = NULL,
                          task_col = "task",
                          formula = NULL,
                          prior = NULL,
                          generator = NULL,
                          seed = NULL,
                          subjects = c("redraw", "fixed"),
                          re_cor = c("none", "within", "all"),
                          scale = c("natural", "link"),
                          levels = c("population", "subject"),
                          correlations = NULL,
                          cor_scale = "link",
                          ml = FALSE,
                          smoke = FALSE,
                          preflight = TRUE,
                          ...,
                          .fitter = NULL) {
  first_model <- grid_first_model(model, grid)
  if (is_component_list(first_model)) {
    given <- c(
      pars = !missing(pars), sds = !missing(sds),
      tasks = !missing(tasks), task_col = !missing(task_col),
      formula = !missing(formula), generator = !missing(generator),
      re_cor = !missing(re_cor)
    )
    return(recovery_grid_components(
      model, first_model, grid,
      dir = dir, reps = reps, cors = cors, covariates = covariates,
      prior = prior, seed = seed, subjects = subjects, scale = scale,
      levels = levels, correlations = correlations, cor_scale = cor_scale,
      ml = ml, smoke = smoke, preflight = preflight,
      dots = rlang::list2(...), fitter = .fitter,
      given = names(given)[given]
    ))
  }
  subjects <- rlang::arg_match(subjects)
  re_cor <- rlang::arg_match(re_cor)
  scale <- rlang::arg_match(scale)
  extraction <- check_extraction_args(levels, correlations, cor_scale)
  check_model_correlations(extraction$correlations, formula, re_cor, tasks)
  ml_args <- check_ml_arg(ml, extraction$levels)
  check_grid(grid)
  dots <- rlang::list2(...)
  setup <- grid_run_setup(grid, reps, dir, smoke, preflight)
  grid <- setup$grid
  reps <- setup$reps
  dir <- setup$dir

  models <- vector("list", nrow(grid))
  model_for <- function(i) {
    if (is.null(models[[i]])) {
      models[[i]] <<- if (!is.function(model)) {
        model
      } else if (i == 1L) {
        first_model
      } else {
        model(grid[i, , drop = FALSE])
      }
      if (is_component_list(models[[i]])) {
        cli::cli_abort(c(
          "The model of grid row {i} is a list of components, but row 1's \\
           is a single model.",
          i = "Every row's {.arg model} must be of the same form."
        ))
      }
      check_model(models[[i]])
    }
    models[[i]]
  }
  formulas <- vector("list", nrow(grid))
  formula_for <- function(i) {
    if (is.null(formulas[[i]])) {
      formulas[[i]] <<- grid_formula(
        formula, grid[i, , drop = FALSE], i, model_for(i), re_cor,
        if (is.null(tasks)) NULL else task_col
      )
    }
    formulas[[i]]
  }

  design <- check_tasks(tasks, task_col, model_for(1L), covariates)
  check_cor_columns(
    grid, model_for(1L), covariates, design$tasks, design$task_col
  )

  cells <- grid_cell_table(nrow(grid), reps, seed)
  sims <- vector("list", nrow(cells))
  runs <- vector("list", nrow(cells))
  first_rep <- vector("list", nrow(grid))

  simulate_cell <- function(i) {
    row <- grid[cells$row[[i]], , drop = FALSE]
    values <- row_values(row, pars, sds, cors, design$tasks, design$task_col)
    sim <- cell_simulation(
      cell_paths(dir, cells$row[[i]], cells$rep[[i]]),
      model_for(cells$row[[i]]), values, row, cells$seed[[i]], subjects,
      first_rep[[cells$row[[i]]]], generator, covariates,
      tasks, task_col
    )
    if (cells$rep[[i]] == 1L) first_rep[[cells$row[[i]]]] <<- sim
    sim
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
  request <- extraction_request(
    extraction$levels, extraction$correlations, extraction$cor_scale,
    links = unlist(links), ml = ml_request(ml_args)
  )

  first_paths <- cell_paths(dir, 1L, 1L)
  run_first <- isTRUE(preflight) &&
    !file.exists(paste0(first_paths$fit, ".rds"))
  if (run_first && file.exists(first_paths$est)) {
    sims[[1L]] <- simulate_cell(1L)
    first_dots <- c(
      if (!is.na(cells$seed[[1L]])) list(seed = cells$seed[[1L]]), dots
    )
    first_key <- cache_key(
      formula_for(1L), sims[[1L]]$data, model_for(1L), prior, first_dots
    )$key
    run_first <- is.null(read_sidecar(first_paths$est, first_key, request))
  }
  if (run_first) {
    if (is.null(sims[[1L]])) sims[[1L]] <- simulate_cell(1L)
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
      cells$seed[[i]], dots, .fitter, request,
      ml = ml_args
    )
  }

  status <- vapply(runs, function(r) r$status, character(1))
  failed <- which(status == "error")
  warn_failed_cells(
    sprintf("row-%d rep %d", cells$row[failed], cells$rep[failed]),
    if (length(failed) > 0L) runs[[failed[[1L]]]]$message
  )
  if (!is.null(ml_args)) {
    ml_status <- vapply(runs, function(r) r$ml_status, character(1))
    ml_failed <- which(ml_status %in% "error")
    warn_failed_ml_cells(
      sprintf("row-%d rep %d", cells$row[ml_failed], cells$rep[ml_failed]),
      if (length(ml_failed) > 0L) runs[[ml_failed[[1L]]]]$ml_message
    )
  }

  out <- score_cells(runs, sims, cells, unlist(links), scale, request)
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
  if (!is.null(ml_args)) {
    attr(out, "ml_cells") <- ml_cell_table(runs, cells, dir, ml_args)
  }
  out
}

#' Warn at the end of a grid about the cells whose ML fit failed
#' @param labels One label per cell whose ML fit failed.
#' @param message The first failure's message.
#' @noRd
warn_failed_ml_cells <- function(labels, message) {
  if (length(labels) == 0L) {
    return(invisible(NULL))
  }
  cli::cli_warn(c(
    "{.fn fit_ml} failed in {length(labels)} cell{?s}: {.val {labels}}.",
    i = "The first message: {message}",
    i = "{cli::qty(length(labels))}The hierarchical rows of \\
         {?that cell/those cells} are scored; the ML rows are missing, \\
         and {.code attr(x, \"ml_cells\")} says where."
  ))
}

#' One row per cell for the ML fit, beside `cells`
#'
#' `n_subjects` and `n_converged` count the subjects of the cell's ML
#' fit, from its `ml_cells` attribute; a subject converged when every one
#' of its terms did. A cell whose ML fit failed, or was not attempted
#' because the hierarchical fit failed, has `NA` counts.
#'
#' `file` names the cell's cached ML fit, and is `NA` on the `optim`
#' route, which does not cache: a path to a file that will never exist is
#' worse than saying there is none.
#'
#' @noRd
ml_cell_table <- function(runs, cells, dir, ml = NULL) {
  count <- function(r, what) {
    rows <- attr(r$ml_estimates, "ml_cells")
    if (is.null(rows)) {
      return(NA_integer_)
    }
    if (identical(what, "n")) nrow(rows) else sum(rows$converged %in% TRUE)
  }
  n_subjects <- vapply(runs, count, integer(1), what = "n")
  n_converged <- vapply(runs, count, integer(1), what = "converged")
  tibble::tibble(
    condition = sprintf("row-%d", cells$row),
    replication = cells$rep,
    status = vapply(runs, function(r) r$ml_status, character(1)),
    elapsed = vapply(runs, function(r) r$ml_elapsed, numeric(1)),
    n_subjects = n_subjects,
    n_converged = n_converged,
    converged = n_converged == n_subjects,
    file = if (identical(ml$method, "optim")) {
      rep(NA_character_, nrow(cells))
    } else {
      paste0(vapply(seq_len(nrow(cells)), function(i) {
        cell_paths(dir, cells$row[[i]], cells$rep[[i]])$ml
      }, character(1)), ".rds")
    }
  )
}
