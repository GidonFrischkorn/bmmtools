# Skip helpers for parts of bmm that not every installed bmm has.
#
# `skip_if_not_installed("bmm")` is not a strong enough guard on its own:
# it answers "is bmm installed", never "does this build carry this model".
# bmm's models arrive one at a time, on `develop` before CRAN. CRAN bmm
# 1.4.1, the minimum bmmtools declares since 0.3.0, exports the whole
# signal-detection stack (`sdt_yn`, `sdt_mafc`, `sdt_ranking`, `sdt_rating`,
# `sdt_cdp`, their `r*()` and `d*()` functions) and `rmpt()`, so on it no
# test here skips. The guard stays for the CI develop leg and for whatever
# lands on `develop` next: a guard checked against one symbol would pass a
# run through to a model that build does not have yet, turning a skip into
# an error (measured 2026-09-17 against CRAN bmm 1.3.2: 9 errors from a
# guard that checked only `sdt_yn`).
#
# The guard is on the symbol rather than on a version number because which
# version carries which model is decided by bmm's releases, so a version
# number here would be a guess that silently keeps skipping once the guess
# is wrong.

#' Skip unless the installed bmm carries a model, its generator, and
#' (optionally) its density
#'
#' Checks the model's own constructor (`model`) and its generator
#' (`r<model>`) always; its density (`d<model>`) only when `density = TRUE`,
#' for a test that calls a density adapter.
#'
#' @param model The model's bmm constructor name, e.g. `"sdt_mafc"`.
#' @param density Also require `d<model>()` to be exported.
#' @noRd
skip_if_no_bmm_model <- function(model, density = FALSE) {
  testthat::skip_if_not_installed("bmm")
  wanted <- c(model, paste0("r", model))
  if (density) {
    wanted <- c(wanted, paste0("d", model))
  }
  ns <- asNamespace("bmm")
  missing <- wanted[!vapply(
    wanted, exists, logical(1),
    envir = ns, inherits = FALSE
  )]
  testthat::skip_if(
    length(missing) > 0,
    sprintf(
      "installed bmm does not export %s", paste(missing, collapse = ", ")
    )
  )
}
