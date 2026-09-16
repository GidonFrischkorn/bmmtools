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
#' exists to catch.
#'
#' @noRd
mock_parameter_draws <- function(fit, n_sims, variables, seed = NULL, ...) {
  selected <- posterior::subset_draws(
    fit$prior_draws,
    variable = variables, regex = TRUE
  )
  posterior::subset_draws(selected, draw = seq_len(n_sims))
}

registerS3method(
  "prior_parameter_draws", "sbcmockfit", mock_parameter_draws,
  envir = asNamespace("bmmtools")
)
