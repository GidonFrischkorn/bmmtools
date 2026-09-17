# Skip helpers for parts of bmm that not every installed bmm has.
#
# `skip_if_not_installed("bmm")` is not a strong enough guard on its own.
# bmm's signal-detection stack --- the `sdt_yn` and `sdt_mafc` models, their
# `r*()` generators and their `d*()` densities --- exists only in the fork
# (1.4.1.9000); CRAN bmm 1.3.2 does not export any of it. Both GitHub
# workflows install `any::bmm`, so on a runner `skip_if_not_installed("bmm")`
# passes and `bmm::sdt_yn()` then errors with "not an exported object".
# Measured 2026-09-17 against CRAN bmm 1.3.2: 9 errors, 7 in test-generate.R
# and 2 in test-prior-check.R.
#
# The guard is on the symbol rather than on a version number because the
# version that will first carry the SDT stack is not yet decided upstream, so
# any number here would be a guess that silently keeps skipping if it is
# wrong.

#' Skip unless the installed bmm carries the signal-detection models
#' @noRd
skip_if_no_bmm_sdt <- function() {
  testthat::skip_if_not_installed("bmm")
  testthat::skip_if_not(
    exists("sdt_yn", envir = asNamespace("bmm"), inherits = FALSE),
    "installed bmm has no signal-detection models"
  )
}
