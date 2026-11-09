#' A model's or fit's parameter table, on every bmm build
#'
#' bmm `develop` deprecates `parameters()` (as of 1.4.0, removed in 1.6.0)
#' in favour of `parameter_info()`, which neither CRAN 1.3.2 nor the
#' 1.3.2.9000 build of 2026-09-30 has. This calls whichever the installed
#' bmm exports. Both are reached by name with `getExportedValue()`, so no
#' literal names a function one of the builds lacks.
#'
#' Measured 2026-10-05 on `develop` `1625e860`: `parameter_info()` is
#' `identical()` to `parameters()` for `mixture2p`, `sdm`, `mixture3p`,
#' `imm`, `m3` (`ss`, `custom`), `sdt_yn`, `ezdm`, `ddm` and the committed
#' `mixture2p` fit: a data frame `parameter`, `description`, `fixed`,
#' `value`, `link`, as `parameters()` returns on the other two builds.
#' `test-bmm-parameter-info.R` keeps that comparison where both exist.
#'
#' A test that needs only the free parameter names, and not bmm's own
#' answer as the thing compared against, uses `model_parameters()`.
#'
#' @noRd
bmm_parameter_info <- function(x, ...) {
  name <- if ("parameter_info" %in% getNamespaceExports("bmm")) {
    "parameter_info"
  } else {
    "parameters"
  }
  getExportedValue("bmm", name)(x, ...)
}

#' The name of a circular model's location: `mu1` on bmm 1.4.1, `mu` from
#' bmm 1.5.0 (#52)
#'
#' Read off the fixed parameters, where every bmm build keeps it (measured
#' 2026-11-09 on CRAN 1.4.1 and `develop` c13dd2c5), and not through the
#' package's own `location_parameter()`, so that a test of that function
#' does not grade itself.
#'
#' @noRd
location_of <- function(model) {
  intersect(c("mu", "mu1"), names(model$fixed_parameters))[[1L]]
}
