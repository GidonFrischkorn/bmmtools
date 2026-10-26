# A mock fitter that draws from R's random number generator, for the
# seed tests of Milestone 34.
#
# bmm draws a fit's initial values from R's generator inside the fit:
# `new_initfun()` returns a closure that cmdstanr calls once per chain
# (measured on bmm 1.4.1). The wrapper draws in the same place, inside
# the call, so a fit is reproducible from its seed only if the call runs
# under it. Drawing before the call would hand the tests a fit that is
# reproducible for free.

#' Wrap a mock fitter so that every call draws once from R's generator
#'
#' The draw is stored on the fit as `draw` and logged in `draws$values`.
#' The fit gains the class `rngmockfit` in front of the inner fit's
#' classes, so the inner mock's methods still apply; the
#' `extract_estimates()` method below puts the draw into the estimates,
#' so that a grid's sidecars depend on it.
#'
#' @param inner A mock as built by [mock_fitter()], [grid_mock_fitter()] or
#'   [sbc_mock_fitter()].
#' @return A list with `fitter`, `calls` (the inner mock's) and `draws`,
#'   an environment holding `values`.
#' @noRd
rng_mock_fitter <- function(inner = mock_fitter()) {
  draws <- new.env(parent = emptyenv())
  draws$values <- numeric()
  inner_fitter <- inner$fitter
  fitter <- function(formula, data, model, prior = NULL, ...) {
    draw <- stats::runif(1L)
    draws$values <- c(draws$values, draw)
    fit <- inner_fitter(formula, data, model, prior, ...)
    fit$draw <- draw
    class(fit) <- c("rngmockfit", class(fit))
    fit
  }
  # SBC serialises the backend, and with it this closure; parented to the
  # namespace, as sbc_mock_fitter() is, the walk stops there
  environment(fitter) <- rlang::env(
    asNamespace("bmmtools"),
    draws = draws, inner_fitter = inner_fitter
  )
  list(fitter = fitter, calls = inner$calls, draws = draws)
}

#' Estimates for an RNG mock fit: the inner mock's, moved by the draw
#' @noRd
extract_estimates_rngmockfit <- function(fit, ...) {
  out <- extract_estimates_mockfit(fit, ...)
  out$estimate <- out$estimate + fit$draw
  out
}

registerS3method(
  "extract_estimates", "rngmockfit", extract_estimates_rngmockfit,
  envir = asNamespace("bmmtools")
)

#' A benchmark fitter that draws from R's generator inside the call
#' @noRd
rng_bench_fitter <- function() {
  draws <- new.env(parent = emptyenv())
  draws$values <- numeric()
  fitter <- function(data, ...) {
    draw <- stats::runif(1L)
    draws$values <- c(draws$values, draw)
    structure(
      list(label = "rng", rows = nrow(data), leapfrog = 1000, draw = draw),
      class = "mock_bench_fit"
    )
  }
  list(fitter = fitter, draws = draws)
}
