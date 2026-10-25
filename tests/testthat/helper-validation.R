# Helpers for the run_validation() tests (Milestone 26).
#
# run_validation() calls the real stage functions -- prior_check(),
# sbc(), recovery_grid(), cross_check() -- through their `.fitter` seam,
# so one mock has to answer every reader those functions apply to a fit.
# Each reader is an S3 generic for which the helpers of the stage tests
# already register a method on a mock class; a fit of this mock carries
# every field those methods read and every class they dispatch on:
#
# * prior_predict_draws() (prior_check): `mockpriorfit`, reads `yrep`
# * prior_parameter_draws() and SBC_fit_to_draws_matrix() (sbc):
#   `sbcmockfit`, reads `prior_draws`
# * extract_estimates() (recovery_grid, cross_check, prior_sd = "fit"):
#   `valfit`, registered below, which adds the model correlations and,
#   on a prior-only fit, a wider posterior SD to the `mockfit` rows
#
# What the mock does not do is behave like a sampler: its estimates are
# 0 with an interval of -10 to 10, so coverage is 1 wherever the truth
# lies inside, and its draws are deterministic sequences, so SBC's Rhat
# is high on every fit. Tests assert on what the runner does with the
# stages' results, not on the numbers the stages compute.

#' Prior draws for every variable a correlated mixture2p fit has
#' @noRd
validation_prior_draws <- function(n_draws = 60L) {
  sbc_prior_draws(n_draws = n_draws)
}

#' A mock fitter for every stage of a validation
#'
#' @param fail_on `NULL`, or a function `(data, dots, prior)` returning
#'   `TRUE` for a call that should error. It is re-parented to the
#'   bmmtools namespace, as the fitter is, so it may use only its
#'   arguments and base or bmmtools functions.
#' @param prior_sd The link-scale SD a prior-only fit reports; the
#'   posterior SD of a fit is 5, so 10 gives a contraction of 0.75.
#' @return A list with `fitter` and `calls` (`n` and the call `log`).
#' @noRd
validation_mock_fitter <- function(fail_on = NULL, prior_sd = 10) {
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  calls$log <- list()
  draws <- validation_prior_draws()
  yrep <- matrix(seq(-3, 3, length.out = 40L), nrow = 8L)
  if (!is.null(fail_on)) environment(fail_on) <- asNamespace("bmmtools")
  fitter <- function(formula, data, model, prior = NULL, ...) {
    dots <- list(...)
    calls$n <- calls$n + 1L
    calls$log[[calls$n]] <- list(
      n_rows = nrow(data),
      ids = levels(factor(data$id)),
      sample_prior = dots$sample_prior,
      prior = prior,
      seed = dots$seed,
      chains = dots$chains,
      iter = dots$iter,
      backend = dots$backend,
      dot_names = names(dots)
    )
    if (!is.null(fail_on) && isTRUE(fail_on(data, dots, prior))) {
      stop("mock sampler failure")
    }
    prior_only <- identical(dots$sample_prior, "only")
    structure(
      list(
        parameters = model_parameters(model)$free,
        ids = levels(factor(data$id)),
        prior_draws = draws,
        yrep = yrep,
        post_sd = if (prior_only) prior_sd else 5
      ),
      class = c(
        "valfit", if (prior_only) "mockpriorfit", "sbcmockfit", "mockfit"
      )
    )
  }
  # SBC hashes the backend, which serialises this closure; parented to
  # the namespace, its environment chain never reaches the attached
  # package, which would make serialize() warn (see helper-sbc.R)
  environment(fitter) <- rlang::env(
    asNamespace("bmmtools"),
    calls = calls, draws = draws, yrep = yrep, fail_on = fail_on,
    prior_sd = prior_sd
  )
  list(fitter = fitter, calls = calls)
}

#' Estimates of a validation mock fit
#'
#' The `mockfit` rows with the fit's posterior SD, and at `level = "cor"`
#' one row per pair of free parameters, centred on 0, as brms names a
#' model correlation (`kappa__thetat`).
#' @noRd
extract_estimates_valfit <- function(fit, level = "population", ...) {
  rows <- extract_estimates_mockfit(fit, level = setdiff(level, "cor"), ...)
  rows$post_sd_link <- rep(fit$post_sd, nrow(rows))
  if ("cor" %in% level && length(fit$parameters) >= 2L) {
    pairs <- utils::combn(fit$parameters, 2L)
    cor <- tibble::tibble(
      term = paste(pairs[1L, ], pairs[2L, ], sep = "__"),
      estimate = 0, ci_low = -0.9, ci_high = 0.9, ci_method = "eti",
      ci_level = 0.95, rhat = 1, ess_bulk = 1000, ess_tail = 1000,
      level = "cor", id = NA_character_, converged = NA,
      ci_low_50 = -0.4, ci_high_50 = 0.4, post_mean_link = 0,
      post_sd_link = 0.4
    )
    rows <- dplyr::bind_rows(rows, cor)
  }
  rows
}

registerS3method(
  "extract_estimates", "valfit", extract_estimates_valfit,
  envir = asNamespace("bmmtools")
)

#' A plan small enough for the mock: two designs, two replications
#'
#' The prior check and the calibration get design points of their own
#' (3 and 4 subjects), so a `fail_on` can tell their fits apart from each
#' other and from the grid's.
#' @noRd
small_validation_plan <- function(model = NULL, ...) {
  model <- model %||% bmm::mixture2p(resp_error = "y")
  args <- list(
    grid = data.frame(n_subjects = c(5L, 6L), n_trials = c(10L, 20L)),
    reps = 2L,
    sbc = list(
      n_subjects = 4L, n_trials = 4L, n_sims = 5L, thin_ranks = 1,
      cores_per_fit = 1
    ),
    prior_check = list(n_subjects = 3L, n_trials = 4L, n_draws = 8L),
    seed = 7
  )
  given <- list(...)
  for (name in names(given)) args[name] <- list(given[[name]])
  rlang::exec(validation_plan, model, !!!args)
}

#' Evaluate `expr`, collecting its warnings instead of signalling them
#' @return A list with `value` and `warnings`, the messages in order.
#' @noRd
collect_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warnings)
}

#' Run a validation against the mock
#'
#' SBC's diagnostics warning fires on every mock run (the draws are
#' sequences, so every Rhat is high) and is muffled by its class; every
#' other warning reaches the test. Messages are dropped.
#' @noRd
mock_validation <- function(plan, dir, mock, ...) {
  suppressMessages(withCallingHandlers(
    run_validation(plan, dir, ..., .fitter = mock$fitter),
    bmmtools_sbc_diagnostics = function(w) invokeRestart("muffleWarning")
  ))
}

#' The fits a mock saw, by kind
#' @noRd
mock_calls_by_kind <- function(calls) {
  vapply(calls$log, function(entry) {
    if (identical(entry$sample_prior, "only")) "prior" else "fit"
  }, character(1))
}
