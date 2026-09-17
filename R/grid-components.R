# The component form of recovery_grid().
#
# `recovery_grid(model = <list of components>)` simulates one set per cell
# with simulate_components(), fits every component through fit_cached() and
# writes a sidecar per component plus one for the set, which holds what
# needs every fit at once: the cross-fit correlations and the subject
# means. The cell order, seeds, smoke mode, preflight, sidecar validation
# and scoring are the single-model grid's, shared through grid.R.

#' Component names the grid cannot use
#'
#' `sim`, `est` and `cor` would collide with the cell file suffixes
#' (`cell-1-rep-1-sim.rds` is the set, `cell-1-rep-1-<comp>.rds` a fit);
#' `sd` and `cor` would make `sd_<comp>_<term>` and `cor_<a>__<b>` grid
#' columns ambiguous.
#'
#' @noRd
reserved_component_names <- function() {
  c("sim", "est", "cor", "sd")
}

#' The files of one cell of a component grid
#' @noRd
component_cell_paths <- function(dir, row, rep, components) {
  stem <- file.path(dir, sprintf("cell-%d-rep-%d", row, rep))
  list(
    sim = paste0(stem, "-sim.rds"),
    fit = stats::setNames(paste0(stem, "-", components), components),
    est = stats::setNames(
      paste0(stem, "-", components, "-est.rds"), components
    ),
    cor = paste0(stem, "-cor.rds")
  )
}

#' Refuse the arguments that belong to the components, and `ml`
#' @noRd
check_component_grid_args <- function(given, subjects, ml = FALSE,
                                      call = rlang::caller_env()) {
  if (!isFALSE(ml)) {
    cli::cli_abort(
      c(
        "{.arg ml} is not supported for components yet.",
        i = "Use {.code ml = FALSE}, or a single model."
      ),
      call = call
    )
  }
  if (length(given) > 0L) {
    cli::cli_abort(
      c(
        "{.arg {given}} cannot be given to {.fn recovery_grid} together with \\
         components.",
        i = "Each component carries its own {.arg pars}, {.arg sds}, \\
             {.arg tasks}, {.arg task_col}, {.arg generator} and \\
             {.arg formula}; set them in {.fn recovery_component}. \\
             Correlated random effects need a component {.arg formula}."
      ),
      call = call
    )
  }
  subjects <- rlang::arg_match(
    subjects, c("redraw", "fixed"),
    error_call = call
  )
  if (identical(subjects, "fixed")) {
    cli::cli_abort(
      c(
        "{.code subjects = \"fixed\"} is not supported for components yet.",
        i = "Use {.code subjects = \"redraw\"}."
      ),
      call = call
    )
  }
  subjects
}

#' The components of every grid row, checked for one form and one set of
#' names
#'
#' @return A list over rows of named component lists.
#' @noRd
grid_component_rows <- function(model, first, grid,
                                call = rlang::caller_env()) {
  rows <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    specs <- if (i == 1L) {
      first
    } else if (is.function(model)) {
      model(grid[i, , drop = FALSE])
    } else {
      model
    }
    if (!is_component_list(specs)) {
      cli::cli_abort(
        c(
          "The model of grid row {i} is not a list of components, but row \\
           1's is.",
          i = "Every row's {.arg model} must be of the same form."
        ),
        call = call
      )
    }
    rows[[i]] <- check_components(specs, call = call)
    if (i > 1L && !identical(names(rows[[i]]), names(rows[[1L]]))) {
      cli::cli_abort(
        c(
          "The components of grid row {i} are {.val {names(rows[[i]])}}, \\
           but row 1's are {.val {names(rows[[1L]])}}.",
          i = "Every row must have the same components, in the same order."
        ),
        call = call
      )
    }
  }
  reserved <- intersect(names(rows[[1L]]), reserved_component_names())
  if (length(reserved) > 0L) {
    cli::cli_abort(
      c(
        "A component in a grid cannot be named {.val {reserved}}.",
        i = "The names {.val {reserved_component_names()}} are used by the \\
             grid's file names and columns; rename the component."
      ),
      call = call
    )
  }
  links <- component_links(rows[[1L]])
  for (i in seq_along(rows)[-1L]) {
    if (!identical(component_links(rows[[i]]), links)) {
      cli::cli_abort(
        c(
          "The components of grid row {i} have different links from row 1.",
          i = "A row-wise {.arg model} must keep one link table across rows."
        ),
        call = call
      )
    }
  }
  rows
}

#' Refuse the model estimator when every component uses the default formula
#' @noRd
check_component_model_cors <- function(correlations, rows,
                                       call = rlang::caller_env()) {
  if (!"model" %in% correlations) {
    return(invisible(NULL))
  }
  has_formula <- vapply(rows, function(specs) {
    any(!vapply(specs, function(s) is.null(s$formula), logical(1)))
  }, logical(1))
  if (!any(has_formula)) {
    cli::cli_abort(
      c(
        "{.code correlations = \"model\"} needs correlated random effects, \\
         but every component uses the default formula, which has none.",
        i = "Give a component a {.arg formula} with correlated terms such as \\
             {.code (1 | p | id)}, or use the {.val draws} and {.val point} \\
             estimators. Model correlations are within one fit only."
      ),
      call = call
    )
  }
  invisible(NULL)
}

# grid columns --------------------------------------------------------------

#' The terms a component's `pars` and `sds` may name, and its terms a
#' correlation may name, both unprefixed
#' @noRd
component_known_terms <- function(spec) {
  info <- model_parameters(spec$model)
  full <- task_terms(info$free, spec$tasks, spec$task_col)
  list(
    values = unique(c(info$free, full, info$fixed)),
    cor = c(full, info$fixed)
  )
}

#' Split `<comp>_<term>` at the first `_`
#' @return A list with `prefix` and `term` (`NA` without a `_`).
#' @noRd
split_component_column <- function(x) {
  has <- grepl("_", x, fixed = TRUE)
  list(
    prefix = ifelse(has, sub("_.*$", "", x), NA_character_),
    term = ifelse(has, sub("^[^_]*_", "", x), NA_character_)
  )
}

#' Check a component grid's columns against every row's components
#'
#' A column is checked in each row where its value is not missing: a
#' missing value leaves the component's own value, so a column may name a
#' task term that only some rows' components have.
#'
#' @noRd
check_component_columns <- function(grid, rows, covariates,
                                    call = rlang::caller_env()) {
  cols <- setdiff(names(grid), "n_subjects")
  if ("n_trials" %in% cols) {
    cli::cli_abort(
      c(
        "{.arg grid} has a column {.val n_trials}, but with components the \\
         number of trials is per component.",
        i = "Use {.val {paste0(\"n_trials_\", names(rows[[1L]]))}}."
      ),
      call = call
    )
  }
  bad <- function(col, i, why) {
    cli::cli_abort(
      c("Grid column {.val {col}} {why} in row {i}."),
      call = call
    )
  }
  for (i in seq_along(rows)) {
    specs <- rows[[i]]
    comps <- names(specs)
    known <- lapply(specs, component_known_terms)
    all_values <- unique(unlist(lapply(known, `[[`, "values")))
    cor_known <- c(
      unlist(lapply(comps, function(nm) prefix_terms(known[[nm]]$cor, nm))),
      names(covariates)
    )
    for (col in cols) {
      if (is.atomic(grid[[col]]) && is.na(grid[[col]][[i]])) next
      if (startsWith(col, "n_trials_")) {
        if (!sub("^n_trials_", "", col) %in% comps) {
          bad(col, i, "names no component")
        }
      } else if (startsWith(col, "cor_")) {
        ab <- cor_column_terms(col)
        is_pair <- length(ab) == 2L && all(ab %in% cor_known) &&
          ab[[1L]] != ab[[2L]]
        if (!is_pair) {
          bad(col, i, "does not name two known prefixed terms or covariates")
        }
      } else if (startsWith(col, "sd_")) {
        parts <- split_component_column(sub("^sd_", "", col))
        if (!parts$prefix %in% comps) bad(col, i, "names no component")
        if (!parts$term %in% known[[parts$prefix]]$values) {
          bad(col, i, "names no term of its component")
        }
      } else {
        parts <- split_component_column(col)
        if (is.na(parts$prefix)) next
        if (parts$prefix %in% comps) {
          if (!parts$term %in% known[[parts$prefix]]$values) {
            bad(col, i, "names no term of its component")
          }
        } else if (parts$term %in% all_values) {
          bad(col, i, "names no component")
        }
      }
    }
  }
  invisible(grid)
}

#' The non-missing values of a row's columns matching `pattern`, as a named
#' list keyed by what is left after `pattern`
#' @noRd
row_column_values <- function(row, pattern) {
  cols <- non_missing_columns(row, grep(pattern, names(row)))
  values <- lapply(cols, function(col) row[[col]])
  stats::setNames(values, sub(pattern, "", cols))
}

#' The names of the columns at `at` whose value in the row is not missing
#' @noRd
non_missing_columns <- function(row, at) {
  cols <- names(row)[at]
  missing <- vapply(cols, function(col) {
    value <- row[[col]]
    length(value) == 1L && is.na(value)
  }, logical(1))
  cols[!missing]
}

#' Apply a component's column values to its `pars` or `sds`
#'
#' Every column was checked to name a term of the component, so each one
#' sets its term; with tasks a bare parameter also sets the task terms
#' already given, and a full task term wins over it, as in the single-model
#' grid.
#'
#' @noRd
override_component_values <- function(x, values, tasks, task_col) {
  if (length(values) == 0L) {
    return(x)
  }
  if (is.null(x)) x <- numeric()
  if (is.null(tasks)) {
    for (term in names(values)) x[[term]] <- values[[term]]
    return(x)
  }
  override_task_values(x, values, names(values), tasks, task_col,
    add_all = TRUE
  )
}

#' One component with a grid row's columns applied
#'
#' A function of its own, so that each wrapped `pars` or `sds` closes over
#' its own component's values rather than a loop variable.
#'
#' @noRd
component_row_spec <- function(spec, row) {
  nm <- spec$name
  # `n_subjects` and `n_trials_*` are never parameter columns, even for a
  # component named `n`
  design <- names(row) == "n_subjects" | startsWith(names(row), "n_trials_")
  par_values <- row_column_values(row[!design], paste0("^", nm, "_"))
  sd_values <- row_column_values(row, paste0("^sd_", nm, "_"))
  n_trials <- row_column_values(row, paste0("^n_trials_", nm, "$"))
  wrap <- function(value, values) {
    if (is.function(value)) {
      force(value)
      return(function() {
        override_component_values(value(), values, spec$tasks, spec$task_col)
      })
    }
    override_component_values(value, values, spec$tasks, spec$task_col)
  }
  spec$pars <- wrap(spec$pars, par_values)
  spec$sds <- wrap(spec$sds, sd_values)
  if (length(n_trials) > 0L) {
    spec$n_trials <- check_count(n_trials[[1L]], paste0("n_trials_", nm))
  }
  spec
}

#' The components and correlations of one grid row
#'
#' Numeric `pars` and `sds` get the row's columns; a zero-argument function
#' becomes one that applies them to its result, so that
#' simulate_components() still evaluates it under the cell seed. `cors` may
#' be a `function(row)`, as in the single-model grid.
#'
#' @return A list with `specs` and `cors`.
#' @noRd
component_row_values <- function(specs, row, cors) {
  specs[] <- lapply(specs, component_row_spec, row = row)
  cor_row <- row[non_missing_columns(row, grep("^cor_", names(row)))]
  if (is.function(cors)) {
    cor_fun <- cors
    cors <- function() override_cors(cor_fun(row), cor_row)
  } else {
    cors <- override_cors(cors, cor_row)
  }
  list(specs = specs, cors = cors)
}

# cells ---------------------------------------------------------------------

#' Simulate one cell's set, or read it back
#' @noRd
component_cell_simulation <- function(path, values, n_subjects, covariates,
                                      seed) {
  if (file.exists(path)) {
    return(readRDS(path))
  }
  set <- simulate_components(
    values$specs,
    n_subjects = n_subjects, cors = values$cors, covariates = covariates,
    seed = if (is.na(seed)) NULL else seed
  )
  write_atomic(path, function(tmp) saveRDS(set, tmp))
  set
}

#' What a component sidecar keeps from one fit: its estimates and, with the
#' model estimator requested, its within-fit model correlations
#' @noRd
extract_component <- function(fit, request) {
  estimates <- extract_estimates(fit, level = request$levels)
  converged <- NULL
  if (nrow(estimates) > 0L) converged <- as.logical(estimates$converged[[1L]])
  cor_estimates <- NULL
  if ("model" %in% request$correlations) {
    cor_estimates <- extract_correlations(
      fit,
      estimator = "model", scale = "link", links = request$links,
      converged = converged
    )
  }
  list(estimates = estimates, cor_estimates = cor_estimates)
}

#' What the set sidecar keeps from all fits: the correlations across and
#' within them, and the subject means under prefixed terms
#'
#' @param verdicts A logical named by component, each fit's `converged`.
#' @noRd
extract_component_set <- function(fits, set, request, verdicts) {
  cor_estimates <- NULL
  if (!is.null(request$correlations)) {
    covariates <- NULL
    if (length(set$covariates) > 0L) {
      covariates <- as.data.frame(set$covariate_data)
    }
    cor_estimates <- extract_correlations(
      fits,
      estimator = request$correlations, covariates = covariates,
      scale = request$cor_scale, links = request$links,
      converged = verdicts
    )
  }
  subject_means <- empty_subject_means()
  if (nrow(set$truth$subjects) > 0L) {
    subject_means <- subject_means_from_draws(
      set_subject_draws(fits), set$truth$covariates
    )
  }
  list(cor_estimates = cor_estimates, subject_means = subject_means)
}

#' The cache key of one component's fit in a cell
#' @noRd
component_key <- function(set, nm, formulas, prior, fit_dots) {
  sim <- set$components[[nm]]
  cache_key(formulas[[nm]], sim$data, sim$model, prior[[nm]], fit_dots)$key
}

#' The extraction request of one component: the grid's, with its own links
#' @noRd
component_request <- function(request, set, nm) {
  request$links <- unlist(set$components[[nm]]$model$links)
  request
}

#' Fit every component of one cell and extract, or read the sidecars
#'
#' A component whose sidecar matches is not fitted and its fit is not read.
#' The set sidecar matches when its keys are the components' current keys;
#' otherwise every fit is obtained through fit_cached(), which reuses cached
#' fits, and the set is extracted again. A component that fails is
#' `"error"`; a failure of the set extraction marks every component, since
#' the cell is then incomplete.
#'
#' @return A list with `components` (per component: `status`, `message`,
#'   `estimates`, `elapsed`, `converged`), `cor_estimates` and
#'   `subject_means`.
#' @noRd
run_component_cell <- function(set, formulas, prior, paths, seed, dots,
                               fitter, request) {
  comps <- names(set$components)
  fit_dots <- c(if (!is.na(seed)) list(seed = seed), dots)
  fits <- list()
  fit_of <- function(nm) {
    if (is.null(fits[[nm]])) {
      sim <- set$components[[nm]]
      cache_args <- c(
        list(
          formula = formulas[[nm]], data = sim$data, model = sim$model,
          file = paths$fit[[nm]], prior = prior[[nm]], .fitter = fitter
        ),
        fit_dots
      )
      fits[[nm]] <<- rlang::exec(fit_cached, !!!cache_args)
    }
    fits[[nm]]
  }

  runs <- lapply(stats::setNames(nm = comps), function(nm) {
    started <- Sys.time()
    comp_request <- component_request(request, set, nm)
    run <- tryCatch(
      {
        key <- component_key(set, nm, formulas, prior, fit_dots)
        extracted <- read_sidecar(paths$est[[nm]], key, comp_request)
        if (is.null(extracted)) {
          extracted <- extract_component(fit_of(nm), comp_request)
          write_sidecar(
            paths$est[[nm]], list(key = key), comp_request, extracted
          )
        }
        list(
          status = "ok", message = NA_character_, key = key,
          estimates = extracted$estimates
        )
      },
      error = function(e) {
        list(status = "error", message = conditionMessage(e), estimates = NULL)
      }
    )
    run$elapsed <- as.double(difftime(Sys.time(), started, units = "secs"))
    run$converged <- NA
    if (!is.null(run$estimates) && nrow(run$estimates) > 0L) {
      run$converged <- as.logical(run$estimates$converged[[1L]])
    }
    run
  })

  out <- list(
    components = runs, cor_estimates = NULL,
    subject_means = empty_subject_means()
  )
  if (!all(vapply(runs, function(r) r$status == "ok", logical(1)))) {
    return(out)
  }
  keys <- vapply(runs, `[[`, character(1), "key")
  set_request <- request
  set_request$levels <- character()
  set_request$links <- set$links
  stored <- read_sidecar(paths$cor, keys, set_request, key_field = "keys")
  if (is.null(stored)) {
    for (nm in comps) {
      failure <- tryCatch(
        {
          fit_of(nm)
          NULL
        },
        error = conditionMessage
      )
      if (!is.null(failure)) {
        out$components[[nm]]$status <- "error"
        out$components[[nm]]$message <- failure
      }
    }
    if (any(vapply(out$components, `[[`, character(1), "status") == "error")) {
      return(out)
    }
    verdicts <- vapply(runs, function(r) as.logical(r$converged), logical(1))
    stored <- tryCatch(
      {
        fit_set <- structure(
          lapply(stats::setNames(nm = comps), fit_of),
          class = c("bmmtools_fit_set", "list")
        )
        extracted <- extract_component_set(fit_set, set, set_request, verdicts)
        write_sidecar(
          paths$cor, list(keys = keys), set_request, extracted,
          levels = FALSE
        )
        extracted
      },
      error = function(e) e
    )
    if (inherits(stored, "error")) {
      for (nm in comps) {
        out$components[[nm]]$status <- "error"
        out$components[[nm]]$message <- conditionMessage(stored)
      }
      return(out)
    }
  }
  out$cor_estimates <- stored$cor_estimates
  out$subject_means <- stored$subject_means
  out
}

#' One cell's run and truth in the form score_cells() reads
#'
#' Estimates of the components that produced a fit, with prefixed terms;
#' the truth reduced to those components; the correlations and subject
#' means only when every component produced one, since they need all fits.
#'
#' @return A list with `run` and `sim`.
#' @noRd
component_score_input <- function(cell, set) {
  status <- vapply(cell$components, `[[`, character(1), "status")
  ok <- names(status)[status == "ok"]
  estimates <- dplyr::bind_rows(lapply(ok, function(nm) {
    rows <- cell$components[[nm]]$estimates
    rows$term <- prefix_terms(rows$term, nm)
    rows
  }))
  truth <- set$truth
  for (name in c("population", "subjects", "sd")) {
    table <- truth[[name]]
    truth[[name]] <- table[sub("_.*$", "", table$term) %in% ok, ]
  }
  complete <- length(ok) == length(status)
  if (!complete) truth$cor <- truth$cor[0L, ]
  list(
    run = list(
      status = if (length(ok) > 0L) "ok" else "error",
      estimates = estimates,
      # bind_cells() labels every ok cell's table, so an empty one here
      cor_estimates = if (complete) {
        cell$cor_estimates
      } else {
        empty_correlation_rows()
      },
      subject_means = if (complete) {
        cell$subject_means
      } else {
        empty_subject_means()
      }
    ),
    sim = list(truth = truth, covariates = set$covariates)
  )
}

#' recovery_grid() with a list of components
#'
#' Called by recovery_grid() once row 1's model is a list of components;
#' `first` is that list. The arguments are recovery_grid()'s, with `given`
#' naming the component arguments the caller supplied.
#'
#' @noRd
recovery_grid_components <- function(model, first, grid, dir, reps, cors,
                                     covariates, prior, seed, subjects, scale,
                                     levels, correlations, cor_scale, ml,
                                     smoke, preflight, dots, fitter, given,
                                     call = rlang::caller_env()) {
  check_component_grid_args(given, subjects, ml, call = call)
  scale <- rlang::arg_match(scale, c("natural", "link"), error_call = call)
  extraction <- check_extraction_args(levels, correlations, cor_scale,
    call = call
  )
  check_grid(grid, required = "n_subjects", call = call)
  setup <- grid_run_setup(grid, reps, dir, smoke, preflight, call = call)
  grid <- setup$grid
  dir <- setup$dir

  rows <- grid_component_rows(model, first, grid, call = call)
  comps <- names(rows[[1L]])
  check_component_model_cors(extraction$correlations, rows, call = call)
  check_component_columns(grid, rows, covariates, call = call)
  prior <- check_component_prior(prior, comps, call = call)
  links <- component_links(rows[[1L]])
  request <- extraction_request(
    extraction$levels, extraction$correlations, extraction$cor_scale,
    links = links
  )

  formulas <- vector("list", nrow(grid))
  formulas_for <- function(r) {
    if (is.null(formulas[[r]])) {
      formulas[[r]] <<- lapply(rows[[r]], component_formula)
    }
    formulas[[r]]
  }
  cells <- grid_cell_table(nrow(grid), setup$reps, seed)
  sims <- vector("list", nrow(cells))
  simulate_cell <- function(i) {
    r <- cells$row[[i]]
    row <- grid[r, , drop = FALSE]
    component_cell_simulation(
      component_cell_paths(dir, r, cells$rep[[i]], comps)$sim,
      component_row_values(rows[[r]], row, cors),
      row$n_subjects,
      covariates, cells$seed[[i]]
    )
  }

  if (isTRUE(preflight)) {
    # `[<-` with list(): the preflight returns NULL when it simulated nothing
    sims[1L] <- list(component_preflight(
      sims[[1L]], simulate_cell, formulas_for(1L), prior,
      component_cell_paths(dir, 1L, 1L, comps), cells$seed[[1L]], request,
      dir, dots, fitter,
      n_cells = nrow(cells), call = call
    ))
  }

  runs <- vector("list", nrow(cells))
  for (i in seq_len(nrow(cells))) {
    if (is.null(sims[[i]])) sims[[i]] <- simulate_cell(i)
    runs[[i]] <- run_component_cell(
      sims[[i]], formulas_for(cells$row[[i]]), prior,
      component_cell_paths(dir, cells$row[[i]], cells$rep[[i]], comps),
      cells$seed[[i]], dots, fitter, request
    )
  }

  cell_table <- component_cell_table(runs, sims, cells, grid, dir, comps)
  failed <- which(cell_table$status == "error")
  first_failed <- NULL
  if (length(failed) > 0L) {
    first_failed <- runs[[cell_table$cell[[failed[[1L]]]]]]$components[[
      cell_table$component[[failed[[1L]]]]
    ]]$message
  }
  warn_failed_cells(
    sprintf(
      "row-%d rep %d (%s)", cells$row[cell_table$cell[failed]],
      cells$rep[cell_table$cell[failed]], cell_table$component[failed]
    ),
    first_failed,
    what = "component"
  )

  inputs <- lapply(seq_along(runs), function(i) {
    component_score_input(runs[[i]], sims[[i]])
  })
  out <- score_cells(
    lapply(inputs, `[[`, "run"), lapply(inputs, `[[`, "sim"), cells, links,
    scale, request
  )
  cell_table$cell <- NULL
  attr(out, "cells") <- cell_table
  attr(out, "grid") <- grid
  out
}

#' The preflight of a component grid: every component of the first cell
#' once, with one chain and 200 iterations
#'
#' Skipped when every component of the first cell has a cached fit or a
#' sidecar that can be used.
#'
#' @return The first cell's set when it had to be simulated, else `sim`.
#' @noRd
component_preflight <- function(sim, simulate_cell, formulas, prior, paths,
                                seed, request, dir, dots, fitter, n_cells,
                                call) {
  comps <- names(formulas)
  has_fit <- file.exists(paste0(paths$fit, ".rds"))
  if (all(has_fit)) {
    return(sim)
  }
  if (any(file.exists(paths$est[!has_fit]))) {
    sim <- sim %||% simulate_cell(1L)
    fit_dots <- c(if (!is.na(seed)) list(seed = seed), dots)
    usable <- vapply(comps[!has_fit], function(nm) {
      key <- component_key(sim, nm, formulas, prior, fit_dots)
      stored <- read_sidecar(
        paths$est[[nm]], key, component_request(request, sim, nm)
      )
      !is.null(stored)
    }, logical(1))
    if (all(usable)) {
      return(sim)
    }
  }
  sim <- sim %||% simulate_cell(1L)
  elapsed <- 0
  for (nm in comps) {
    part <- sim$components[[nm]]
    elapsed <- elapsed + in_component(nm,
      {
        run_preflight(
          part, part$model, formulas[[nm]], prior[[nm]], dir, dots, fitter,
          file = paste0("preflight-", nm), call = call
        )
      },
      call = call
    )
  }
  cli::cli_inform(
    "Preflight passed in {round(elapsed, 1)} s; running {n_cells} cell{?s}."
  )
  sim
}

#' The `cells` attribute of a component grid: one row per cell and
#' component, with the cell's index in `cell`
#' @noRd
component_cell_table <- function(runs, sims, cells, grid, dir, comps) {
  dplyr::bind_rows(lapply(seq_along(runs), function(i) {
    r <- cells$row[[i]]
    paths <- component_cell_paths(dir, r, cells$rep[[i]], comps)
    parts <- runs[[i]]$components[comps]
    tibble::tibble(
      cell = i,
      condition = sprintf("row-%d", r),
      replication = cells$rep[[i]],
      component = comps,
      n_subjects = as.integer(grid$n_subjects[[r]]),
      n_trials = vapply(comps, function(nm) {
        as.integer(sims[[i]]$components[[nm]]$n_trials)
      }, integer(1), USE.NAMES = FALSE),
      seed = cells$seed[[i]],
      file = paste0(unname(paths$fit), ".rds"),
      status = vapply(parts, `[[`, character(1), "status", USE.NAMES = FALSE),
      elapsed = vapply(parts, `[[`, numeric(1), "elapsed", USE.NAMES = FALSE),
      converged = vapply(parts, function(p) as.logical(p$converged),
        logical(1),
        USE.NAMES = FALSE
      )
    )
  }))
}
