# Helpers for the SBC naming-layer tests (spec 6, section 1).
#
# bmm's mock backend returns a bmmfit whose $fit is the bare mock_fit
# number, so nothing that reads draws runs on it (measured 2026-09-08,
# [[bmm-mock-backend-no-draws]]). The prior draws are therefore reached
# through the internal generic prior_parameter_draws(), for which this
# file registers a method, in the same way helper-prior-check.R
# registers prior_predict_draws().

#' A `bmmodel`-shaped object with the free parameters a test needs
#'
#' `model_parameters()` reads the slots rather than calling
#' `bmm::parameters()`, so this stands in without bmm loaded: the free
#' set is the names of `links` minus the names of `fixed_parameters`.
#'
#' @noRd
sbc_model <- function(free = c("kappa", "thetat"), fixed = NULL) {
  names <- c(free, names(fixed))
  links <- as.list(stats::setNames(rep("identity", length(names)), names))
  structure(
    list(name = "toy", links = links, fixed_parameters = fixed),
    class = "bmmodel"
  )
}

#' The brms draw names the mock fit carries by default
#' @noRd
sbc_mock_variables <- function() {
  c("b_kappa_Intercept", "b_thetat_Intercept", "sd_id__kappa_Intercept")
}

#' A mock fit carrying prior draws under brms's own variable names
#'
#' The values are a deterministic function of the row and the column, so
#' a test can name any draw it wants to find again after subsetting.
#'
#' @noRd
sbc_mock_fit <- function(variables = sbc_mock_variables(), n_draws = 50L) {
  values <- vapply(
    seq_along(variables),
    function(k) seq_len(n_draws) + 100 * k,
    numeric(n_draws)
  )
  dimnames(values) <- list(NULL, variables)
  structure(
    list(prior_draws = posterior::as_draws_matrix(values)),
    class = c("sbcmockfit", "mockfit")
  )
}

#' Prior draws for a mock fit: the same subset-then-sample path
#'
#' Deliberately not a stub that ignores `variables`: the point of the
#' generic is that the selection and the draw count are exercised, and a
#' method that returned everything would hide exactly the mismatch 6.2
#' exists to catch. It goes through `subsample_prior_draws()`, the body
#' the default method uses, so the coverage check and the draw-count
#' error are the real ones rather than a second implementation of them.
#'
#' @noRd
mock_parameter_draws <- function(fit, n_sims, variables, seed = NULL, ...) {
  subsample_prior_draws(
    fit$prior_draws, check_count(n_sims, "n_sims"), variables, seed
  )
}

registerS3method(
  "prior_parameter_draws", "sbcmockfit", mock_parameter_draws,
  envir = asNamespace("bmmtools")
)

#' Prior draws a `mixture2p` simulation can actually be run from
#'
#' `sbc_mock_fit()`'s values are indices, which is all the naming layer
#' needs. `sbc()` hands them to `simulate_recovery()`, so they have to be
#' values the model can generate from: `kappa` on the log link and
#' `thetat` on the logit link, the SDs positive, the correlation inside
#' (-1, 1). Deterministic and free of the random number stream, so a test
#' can assert on a particular row.
#'
#' @param subjects The subject labels whose
#'   `r_<group>__<par>[<label>,Intercept]` draws the matrix carries, or
#'   `NULL` for none. Only the tests of the
#'   subject level ask for them: a mock that carried every possible name
#'   regardless of what was asked is the failure mode described under
#'   [sbc_mock_fitter()].
#' @noRd
sbc_prior_draws <- function(n_draws = 60L, group = "id", subjects = NULL) {
  k <- seq_len(n_draws)
  values <- cbind(
    b_kappa_Intercept = log(4) + 0.01 * k,
    b_thetat_Intercept = stats::qlogis(0.7) + 0.005 * k,
    sd_kappa = 0.2 + 0.002 * k,
    sd_thetat = 0.3 + 0.002 * k,
    cor_pair = 0.5 * cos(k / 7)
  )
  colnames(values) <- c(
    "b_kappa_Intercept", "b_thetat_Intercept",
    paste0("sd_", group, "__kappa_Intercept"),
    paste0("sd_", group, "__thetat_Intercept"),
    paste0("cor_", group, "__kappa_Intercept__thetat_Intercept")
  )
  for (term in c("kappa", "thetat")) {
    for (i in seq_along(subjects)) {
      # deviations around zero, distinct per subject and term
      column <- 0.1 * i * sin(k / 5 + i) * if (term == "kappa") 1 else -1
      values <- cbind(values, column)
      colnames(values)[ncol(values)] <- paste0(
        "r_", group, "__", term, "[", subjects[[i]], ",Intercept]"
      )
    }
  }
  posterior::as_draws_matrix(values)
}

#' A mock fitter for `sbc()`: one object serving as prior fit and as fit
#'
#' The prior fit is read through `prior_parameter_draws()` and every
#' dataset fit through SBC's `SBC_fit_to_draws_matrix()`, so the same
#' object answers both. What the test asks about afterwards is the log:
#' how often the fitter was called, and with what.
#'
#' @param variables Which columns of [sbc_prior_draws()] the fit carries.
#'   `NULL` carries every column **except** the `cor_` one, because that
#'   is what a real fit of the default `recovery_formula()` has: its
#'   group terms are uncorrelated, so brms writes no correlation. A mock
#'   that shipped every possible name regardless of the formula under
#'   test once let a HIGH through a green gate --- `sbc()` asked a real
#'   fit for a `cor_` draw it did not have, and every test "found" the
#'   one the mock happened to carry. Pass `variables` to build a fit
#'   whose draws deliberately do or do not cover what `sbc()` resolved.
#' @param subjects Passed to [sbc_prior_draws()]: the labels whose `r_`
#'   draws the fit carries. `NULL`, the default, carries none, so a test
#'   of the subject level has to say which labels its prior fit saw ---
#'   which is also what a real prior fit does, since it is fitted to the
#'   user's `data` and names its `r_` draws after that column's values.
#' @param draws A `draws_matrix` to carry instead of [sbc_prior_draws()]'s,
#'   for a fit whose names the default matrix cannot have --- a
#'   cell-means formula's `b_kappa_task1`, say. `variables` and
#'   `subjects` are then ignored: the caller has said exactly what the
#'   fit carries.
#' @return A list with `fitter` to inject and `calls`, an environment
#'   holding `n` and the call log.
#' @noRd
sbc_mock_fitter <- function(variables = NULL, n_draws = 60L, group = "id",
                            subjects = NULL, draws = NULL) {
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  calls$log <- list()
  if (is.null(draws)) {
    draws <- sbc_prior_draws(n_draws, group, subjects)
    variables <- variables %||%
      grep("^cor_", posterior::variables(draws), invert = TRUE, value = TRUE)
    draws <- posterior::subset_draws(draws, variable = variables)
  }
  fitter <- function(formula, data, model, prior = NULL, ...) {
    calls$n <- calls$n + 1L
    dots <- list(...)
    calls$log[[calls$n]] <- list(
      n_rows = nrow(data), columns = names(data),
      ids = unique(as.character(data[[1L]])),
      data = data, prior = prior,
      sample_prior = dots$sample_prior, cores = dots$cores,
      seed = dots$seed, chains = dots$chains, iter = dots$iter
    )
    structure(
      list(prior_draws = draws), class = c("sbcmockfit", "mockfit")
    )
  }
  # SBC hashes the backend with rlang::hash(), which serialises this
  # closure. Left in the test environment its parent chain reaches the
  # attached `package:bmmtools`, and serialize() then warns once per run
  # about it -- a warning that says nothing and would bury a real one.
  # Re-parenting to the namespace keeps everything the fitter needs and
  # stops the walk there.
  environment(fitter) <- rlang::env(
    asNamespace("bmmtools"),
    calls = calls, draws = draws
  )
  list(fitter = fitter, calls = calls)
}

#' The calls that fitted a data set, i.e. every call after the prior fit
#' @noRd
sbc_dataset_calls <- function(calls) {
  calls$log[-1L]
}

# SBC reads a fit through its own generic, so the mock fit needs a
# method there too. Registered into SBC's namespace, the way
# helper-generate.R registers extract_estimates.mockfit into bmmtools';
# guarded because SBC is in Suggests and the suite has to load without
# it.
if (requireNamespace("SBC", quietly = TRUE)) {
  registerS3method(
    "SBC_fit_to_draws_matrix", "sbcmockfit",
    function(fit, ...) fit$prior_draws,
    envir = asNamespace("SBC")
  )
}
