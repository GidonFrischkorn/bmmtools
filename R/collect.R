# Reading a finished grid back from its cell files, on any machine.

#' The cell files of one cell, read back into the shape score_cells() wants
#'
#' A cell with no sidecar is `status = "missing"`, which excludes it from
#' the scoring exactly as a failed fit is excluded, and says so in
#' `cells` rather than silently shortening the result.
#'
#' @noRd
collect_cell <- function(paths) {
  sim <- if (file.exists(paths$sim)) {
    tryCatch(upgrade_simulation(readRDS(paths$sim)), error = function(e) NULL)
  }
  stored <- if (file.exists(paths$est)) {
    tryCatch(readRDS(paths$est), error = function(e) NULL)
  }
  if (is.null(sim) || is.null(stored)) {
    return(list(sim = sim, run = list(
      status = "missing", message = NA_character_, estimates = NULL,
      elapsed = NA_real_, converged = NA,
      diagnostics = empty_diagnostics(), fit_seconds = NA_real_,
      sampler = cell_sampler(list()),
      ml_estimates = NULL, ml_status = NA_character_,
      ml_message = NA_character_, ml_elapsed = NA_real_
    )))
  }
  converged <- NA
  if (!is.null(stored$estimates) && nrow(stored$estimates) > 0L) {
    converged <- stored$estimates$converged[[1L]]
  }
  # a cell file written before the 50 % interval (decision 49) or the
  # posterior moments (Milestone 22) existed: what it lacks is unknown,
  # filled with NA, and collect_grid() says so
  tables <- c("estimates", "cor_estimates", "ml_estimates")
  lacks <- function(columns) {
    !all(vapply(
      tables, function(name) has_current_columns(stored[[name]], columns),
      logical(1)
    ))
  }
  no_inner <- lacks(inner_interval_columns())
  no_moments <- lacks(posterior_moment_columns())
  for (name in tables) {
    if (!is.null(stored[[name]])) {
      stored[[name]] <- fill_optional_estimates(stored[[name]])
    }
  }
  list(sim = sim, no_inner = no_inner, no_moments = no_moments, run = list(
    status = "ok", message = NA_character_,
    estimates = stored$estimates,
    cor_estimates = stored$cor_estimates,
    subject_means = stored$subject_means,
    diagnostics = stored$diagnostics %||% empty_diagnostics(),
    fit_seconds = stored$fit_seconds %||% NA_real_,
    sampler = stored$sampler %||% cell_sampler(list()),
    # a collected grid never ran a cell, so it has no wall time of its own
    elapsed = NA_real_,
    converged = converged,
    ml_estimates = stored$ml_estimates,
    ml_status = if (is.null(stored$ml_estimates)) NA_character_ else "ok",
    ml_message = NA_character_,
    ml_elapsed = NA_real_
  ))
}

#' Check that the sidecars hold what is being asked of them
#' @noRd
check_collected <- function(stored, levels, correlations,
                            call = rlang::caller_env()) {
  missing_levels <- setdiff(levels, stored$levels)
  if (length(missing_levels) > 0L) {
    cli::cli_abort(
      c(
        "The cells were not extracted at {.val {missing_levels}}.",
        i = "They hold {.val {stored$levels}}.",
        i = "Rerun the grid to extract another level; the fits are what \\
             a level is extracted from."
      ),
      call = call
    )
  }
  missing_cors <- setdiff(correlations, stored$correlations)
  if (length(missing_cors) > 0L) {
    cli::cli_abort(
      c(
        "The cells hold no {.val {missing_cors}} correlation estimator.",
        i = if (is.null(stored$correlations)) {
          "They hold none."
        } else {
          "They hold {.val {stored$correlations}}."
        }
      ),
      call = call
    )
  }
  invisible(NULL)
}

#' Read a finished recovery grid from its cell files
#'
#' What [recovery_grid()] returns, rebuilt from the files it left in
#' `dir`, without fitting anything and without reading a fit. The
#' per-cell sidecars hold everything that was scored, so a study can run
#' on one machine and be assembled on another.
#'
#' This is a reader, not a resume. It checks neither the cache key nor
#' the version of bmmtools that wrote the files, because on the machine
#' that collects the results neither can match: the key covers the
#' installed bmm, brms and Stan toolchain. A resume is what it has always
#' been --- rerunning the identical [recovery_grid()] call in the same
#' directory --- and it is the resume, not this function, that guarantees
#' the files describe the call.
#'
#' Cell files written before bmmtools stored the central 50 % interval
#' are read with `ci_low_50`, `ci_high_50` and `covered_50` set to `NA`,
#' so their `coverage_50` is `NA`, and a message says how many there
#' were. Rerunning the grid with the fits still in `dir` re-extracts them
#' without refitting. Cell files written before it stored the posterior
#' mean and SD are read with `post_mean_link` and `post_sd_link` set to
#' `NA`, and a second message counts them.
#'
#' @param dir The directory a [recovery_grid()] run wrote, holding
#'   `grid.rds` and the cell files. Written by bmmtools 0.2.0 or later;
#'   an earlier run left no `grid.rds` and cannot be collected.
#' @param scale The scale to score on, as in [recover()]. `NULL`, the
#'   default, uses the scale the grid was run on. Giving it re-scores the
#'   stored rows, so a grid run on the link scale can be read again on
#'   the natural one without refitting.
#' @param levels,correlations What to score, as in [recovery_grid()].
#'   `NULL` uses what the grid requested. Either may be narrowed to a
#'   subset of what the cells hold; asking for something they do not hold
#'   is an error, because a level is extracted from a fit and the fits
#'   are not read here.
#'
#' @return A `bmmtools_recovery` with the attributes [recovery_grid()]
#'   gives it. `cells$elapsed` is `NA`: no cell was run, so none took
#'   any time. `cells$fit_seconds` is the time each fit took when it was
#'   run, which is what a runtime table wants and what survives the trip
#'   between machines.
#'
#' @seealso [recovery_grid()], which writes what this reads.
#'
#' @examples
#' \dontrun{
#' # on the server
#' recovery_grid(model, grid, pars, dir = "runs/mixture2p", reps = 20)
#'
#' # on a laptop, with only runs/mixture2p/*-sim.rds and *-est.rds copied
#' out <- collect_grid("runs/mixture2p")
#' summary(out)
#' collect_grid("runs/mixture2p", scale = "link")
#' }
#'
#' @export
collect_grid <- function(dir, scale = NULL, levels = NULL,
                         correlations = NULL) {
  if (!is.character(dir) || length(dir) != 1L || is.na(dir)) {
    cli::cli_abort(
      "{.arg dir} must be a single path, not {.obj_type_friendly {dir}}."
    )
  }
  path <- file.path(dir, "grid.rds")
  if (!file.exists(path)) {
    cli::cli_abort(c(
      "{.file {path}} does not exist, so there is no grid to collect.",
      i = "A directory written before bmmtools 0.2.0 has no record of the \\
           call that made it; rerun the grid to write one.",
      i = "The component form of {.fn recovery_grid} is not collected."
    ))
  }
  record <- readRDS(path)
  if (!is.list(record) || is.null(record$grid)) {
    cli::cli_abort("{.file {path}} is not a grid record.")
  }

  scale <- scale %||% record$scale
  levels <- levels %||% record$levels
  correlations <- correlations %||% record$correlations
  check_collected(record, levels, correlations)

  grid <- record$grid
  cells <- grid_cell_table(nrow(grid), record$reps, record$seed)
  links <- unlist(record$links)
  request <- extraction_request(
    levels, correlations, record$cor_scale,
    links = links, ml = record$ml, convergence = record$convergence
  )

  collected <- lapply(seq_len(nrow(cells)), function(i) {
    collect_cell(cell_paths(dir, cells$row[[i]], cells$rep[[i]]))
  })
  sims <- lapply(collected, function(x) x$sim)
  runs <- lapply(collected, function(x) x$run)
  n_no_inner <- sum(vapply(
    collected, function(x) isTRUE(x$no_inner), logical(1)
  ))
  if (n_no_inner > 0L) {
    cli::cli_inform(c(
      "{n_no_inner} cell file{?s} {?has/have} no 50 % interval: written \\
       before bmmtools stored one.",
      i = "Their {.field coverage_50} is {.code NA}. Rerunning the grid \\
           with the fits in place re-extracts them without refitting."
    ))
  }
  n_no_moments <- sum(vapply(
    collected, function(x) isTRUE(x$no_moments), logical(1)
  ))
  if (n_no_moments > 0L) {
    cli::cli_inform(c(
      "{n_no_moments} cell file{?s} {?has/have} no posterior mean and SD: \\
       written before bmmtools stored them.",
      i = "Their {.field post_mean_link} and {.field post_sd_link} are \\
           {.code NA}. Rerunning the grid with the fits in place \\
           re-extracts them without refitting."
    ))
  }

  status <- vapply(runs, function(r) r$status, character(1))
  missing <- which(status == "missing")
  if (length(missing) > 0L) {
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    first <- sprintf(
      "row-%d rep %d", cells$row[missing[[1L]]], cells$rep[missing[[1L]]]
    )
    cli::cli_warn(c(
      "{length(missing)} cell{?s} {?has/have} no result on disk and \\
       {?is/are} not scored.",
      i = "The first is {.val {first}}.",
      i = "{.code attr(x, \"cells\")$status} names them all."
    ))
  }

  out <- score_cells(runs, sims, cells, links, scale, request)
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
    converged = vapply(runs, function(r) as.logical(r$converged), logical(1)),
    !!!cell_diagnostic_columns(runs)
  )
  attr(out, "grid") <- grid
  if (!is.null(record$ml)) {
    attr(out, "ml_cells") <- ml_cell_table(runs, cells, dir, record$ml)
  }
  out
}
