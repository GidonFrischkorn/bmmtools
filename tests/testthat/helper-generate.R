# Helpers for the generate-layer and grid tests.
#
# The grid is tested end to end without a sampler: the mock fitter from
# helper-cache.R stands in for bmm::bmm(), and the method below turns
# what it returns into an estimates tibble, so extract_estimates(),
# recover() and summary() all run on it.

#' A mock fitter whose fits remember the model and the subject ids
#' @noRd
grid_mock_fitter <- function(fail_on = NULL) {
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  calls$log <- list()
  fitter <- function(formula, data, model, prior = NULL, ...) {
    calls$n <- calls$n + 1L
    dots <- list(...)
    calls$log[[calls$n]] <- list(
      n_rows = nrow(data), ids = levels(data$id), seed = dots$seed,
      chains = dots$chains, iter = dots$iter
    )
    if (!is.null(fail_on) && fail_on(calls$n, data, dots)) {
      stop("mock sampler failure")
    }
    structure(
      list(
        parameters = mock_terms(formula, data, model),
        ids = levels(data$id)
      ),
      class = "mockfit"
    )
  }
  list(fitter = fitter, calls = calls)
}

#' The terms a fit of `formula` would have: one per free parameter, or one
#' per parameter and task when the formula has cell means (`0 + task`)
#'
#' The task column is the variable of the first formula's right-hand side
#' other than the grouping and correlation ids, as `recovery_formula()`
#' writes it.
#'
#' @noRd
mock_terms <- function(formula, data, model) {
  p <- bmm::parameters(model)
  free <- p$parameter[!p$fixed]
  rhs <- if (length(formula) > 0L) all.vars(formula[[1L]][[3L]]) else NULL
  task_col <- setdiff(rhs, c("id", "p"))
  if (length(task_col) != 1L || !is.factor(data[[task_col]])) {
    return(free)
  }
  task_terms(free, levels(data[[task_col]]), task_col)
}

#' Estimates for a mock fit: 0 with a wide interval, converged
#'
#' The `"sd"` level repeats the population rows under that level, so a
#' grid can score SDs; `"cor"` gives no rows, as for an uncorrelated fit.
#'
#' @noRd
extract_estimates_mockfit <- function(fit,
                                      level = c("population", "subject"),
                                      group = NULL,
                                      ci_level = 0.95,
                                      ci_method = "eti",
                                      drop_constants = TRUE,
                                      converged = NULL,
                                      ...) {
  if (missing(level)) level <- "population"
  one <- function(term, id, lvl) {
    tibble::tibble(
      term = term, estimate = 0, ci_low = -10, ci_high = 10,
      ci_method = "eti", ci_level = ci_level, rhat = 1, ess_bulk = 1000,
      ess_tail = 1000, level = lvl, id = id, converged = TRUE
    )
  }
  pieces <- list()
  if ("population" %in% level) {
    pieces$population <- one(fit$parameters, NA_character_, "population")
  }
  if ("subject" %in% level) {
    pieces$subject <- dplyr::bind_rows(lapply(fit$ids, function(i) {
      one(fit$parameters, i, "subject")
    }))
  }
  if ("sd" %in% level) {
    pieces$sd <- one(fit$parameters, NA_character_, "sd")
  }
  dplyr::bind_rows(pieces)
}

registerS3method(
  "extract_estimates", "mockfit", extract_estimates_mockfit,
  envir = asNamespace("bmmtools")
)

#' Subject draws for a mock fit: deterministic, with spread in every margin
#'
#' Ten iterations and two chains per id and parameter. The values are a
#' smooth function of the indices, so a test can rely on them being the
#' same on every call and on every subject and draw differing.
#'
#' @noRd
extract_subject_draws_mockfit <- function(fit, group = NULL, ...) {
  ids <- as.character(fit$ids)
  terms <- fit$parameters
  n_iter <- 10L
  n_chain <- 2L
  out <- array(
    NA_real_,
    dim = c(n_iter, n_chain, length(ids), length(terms)),
    dimnames = list(
      iteration = as.character(seq_len(n_iter)),
      chain = as.character(seq_len(n_chain)),
      id = ids,
      term = terms
    )
  )
  for (t in seq_along(terms)) {
    for (i in seq_along(ids)) {
      draw <- outer(seq_len(n_iter), seq_len(n_chain), function(it, ch) {
        sin(i * t + 0.37 * it + 1.3 * ch) / 4
      })
      out[, , i, t] <- cos(1.7 * i * t) + draw
    }
  }
  structure(out, group = group %||% "id")
}

registerS3method(
  "extract_subject_draws", "mockfit", extract_subject_draws_mockfit,
  envir = asNamespace("bmmtools")
)

#' A long subject truth (`id`, `term`, `true_value`) as one column per term
#' @noRd
subjects_wide <- function(subjects) {
  terms <- unique(subjects$term)
  out <- lapply(terms, function(t) subjects$true_value[subjects$term == t])
  stats::setNames(as.data.frame(out), terms)
}

#' A two-row grid over subjects and trials
#' @noRd
small_grid <- function() {
  data.frame(n_subjects = c(3L, 4L), n_trials = c(10L, 20L))
}
