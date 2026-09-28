# Skip helpers for parts of bmm that not every installed bmm has.
#
# `skip_if_not_installed("bmm")` is not a strong enough guard on its own:
# bmm's signal-detection stack --- the `sdt_yn` and `sdt_mafc` models, their
# `r*()` generators and their `d*()` densities --- arrives one model at a
# time. CRAN bmm 1.3.2 exports none of it; popov-lab's `develop` branch
# exports `sdt_yn` and, since 2026-09-28 (#361), `sdt_mafc`, but not yet
# `sdt_ranking`, `sdt_rating`, `sdt_cdp` or `rmpt`. A guard checked against
# one symbol would pass a `develop` run through to a model that branch does
# not have yet, turning a skip into an error --- measured 2026-09-17 against
# CRAN bmm 1.3.2: 9 errors, 7 in test-generate.R and 2 in test-prior-check.R,
# all from a guard that checked only `sdt_yn`.
#
# The guard is on the symbol rather than on a version number because no
# release has decided which version carries which model, so a version
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
