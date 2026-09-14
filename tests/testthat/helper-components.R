# Helpers for the component tests (spec 5, section 5.5a).
#
# `arrayfit` is a fit class that returns whatever subject-draws array and
# cor rows it was built with, so cross-fit correlations can be checked
# against hand-built arrays of any size, including unequal draw counts.

#' A fit that hands back a stored subject-draws array and cor rows
#' @noRd
array_fit <- function(draws, cor = NULL) {
  structure(list(draws = draws, cor = cor), class = "arrayfit")
}

#' @noRd
extract_subject_draws_arrayfit <- function(fit, group = NULL, ...) {
  fit$draws
}

#' @noRd
extract_estimates_arrayfit <- function(fit, level = "population", ...) {
  if (is.null(fit$cor)) {
    return(tibble::tibble(
      term = character(), estimate = double(), ci_low = double(),
      ci_high = double(), rhat = double(), ess_bulk = double(),
      ess_tail = double()
    ))
  }
  fit$cor
}

registerS3method(
  "extract_subject_draws", "arrayfit", extract_subject_draws_arrayfit,
  envir = asNamespace("bmmtools")
)
registerS3method(
  "extract_estimates", "arrayfit", extract_estimates_arrayfit,
  envir = asNamespace("bmmtools")
)

#' A generator that costs nothing: one constant response per trial
#' @noRd
zero_generator <- function(pars, n_trials, model) {
  data.frame(y = rep(0, n_trials))
}

#' A generator that uses the random number stream, so a shifted stream shows
#' @noRd
noisy_generator <- function(pars, n_trials, model) {
  data.frame(y = stats::rnorm(n_trials))
}
