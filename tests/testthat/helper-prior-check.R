# Helpers for the prior-check tests (spec 4).
#
# The mock *backend* cannot stand in here: bmm::bmm(backend = "mock")
# returns a bmmfit whose $fit is the scalar mock_fit, so
# brms::posterior_predict() fails on it (measured 2026-09-08). The
# prediction step is therefore reached through the internal generic
# prior_predict_draws(), for which this file registers a method, in the
# same way helper-generate.R registers extract_estimates.mockfit.

#' A mock fitter whose fits carry prior-predictive draws
#'
#' @param yrep The matrix `prior_predict_draws()` will subset.
#' @param fail_on A function of the call number; `TRUE` makes that call
#'   error, so the failure path can be tested.
#' @return A list with `fitter` to inject and `calls`, an environment
#'   holding `n`, the call log, and the `ndraws` last asked for.
#' @noRd
prior_mock_fitter <- function(yrep = NULL, fail_on = NULL) {
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  calls$log <- list()
  calls$ndraws <- NA_integer_
  if (is.null(yrep)) {
    yrep <- matrix(seq_len(40) / 10, nrow = 8L, ncol = 5L)
  }
  fitter <- function(formula, data, model, prior = NULL, ...) {
    calls$n <- calls$n + 1L
    dots <- list(...)
    calls$log[[calls$n]] <- list(
      prior = prior,
      sample_prior = dots$sample_prior,
      seed = dots$seed,
      chains = dots$chains
    )
    if (!is.null(fail_on) && isTRUE(fail_on(calls$n))) {
      stop("mock sampler failure")
    }
    structure(
      list(yrep = yrep, calls = calls, call_number = calls$n),
      class = c("mockpriorfit", "mockfit")
    )
  }
  list(fitter = fitter, calls = calls)
}

#' Prior-predictive draws for a mock fit: the first `ndraws` rows
#' @noRd
mock_prior_draws <- function(fit, ndraws, seed = NULL, ...) {
  fit$calls$ndraws <- ndraws
  fit$yrep[seq_len(min(ndraws, nrow(fit$yrep))), , drop = FALSE]
}

registerS3method(
  "prior_predict_draws", "mockpriorfit", mock_prior_draws,
  envir = asNamespace("bmmtools")
)

#' A stand-in brmsprior
#'
#' `prior_sets()` must recognise a `brmsprior` before it recognises a
#' list, because a brmsprior is a data frame and a data frame is a list.
#' Building one by hand keeps the class tests free of brms.
#'
#' @noRd
fake_prior <- function(prior = "normal(0, 1)") {
  structure(
    data.frame(prior = prior, class = "b", stringsAsFactors = FALSE),
    class = c("brmsprior", "data.frame")
  )
}

#' Run prior_check() against the mock fitter with the fake model
#' @noRd
mock_prior_check <- function(mock, ...) {
  prior_check(
    fake_model(),
    fake_bmmformula(a = a ~ 1),
    fake_data(),
    .fitter = mock$fitter,
    ...
  )
}

#' A matrix whose statistics can be written down
#'
#' Twelve values 0 to 9 and two at 10: with `range = c(0, 10)` the floor
#' rate is 1/12, the ceiling rate 2/12, and the type-7 quantiles are
#' 5.5, 9.9 and 10.
#'
#' @noRd
known_yrep <- function() {
  matrix(c(0:9, 10, 10), nrow = 4L, ncol = 3L)
}
