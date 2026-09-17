# The generator adapters.
#
# bmm's r<model>() functions do not take the model's parameter names on
# the model's link scale (measured 2026-09-07 and 2026-09-08), so each
# supported model gets a small function that maps one subject's
# natural-scale parameters onto the generator's arguments and the
# generator's output onto the model's column names. Column names come
# from the model object, never from literals. A user-supplied
# `generator =` always takes precedence over the built-in adapter. An
# upstream change letting r<model>() take the model's own parameter
# names would empty this table.

#' Resolve a bmm function by name, at call time
#'
#' The signal-detection stack (`sdt_yn`, `sdt_mafc` and their `r*()` and
#' `d*()` functions) is in no released bmm, so a literal `bmm::rsdt_yn()`
#' makes `R CMD check`'s "checking dependencies in R code" report a
#' missing object on any machine with a released bmm installed. That is a
#' WARNING, and `r-lib/actions` checks with `error-on = "warning"`, so it
#' fails the workflow. Resolving the name here keeps the check quiet and
#' turns an absent symbol into a named error rather than R's bare "not an
#' exported object".
#'
#' Used only for the four fork-only names; every other bmm call in the
#' package stays a literal `bmm::` call, which is what documents the
#' dependency.
#'
#' @param name The name of an exported bmm function.
#'
#' @return The function.
#' @noRd
bmm_fun <- function(name) {
  if (!requireNamespace("bmm", quietly = TRUE)) {
    cli::cli_abort("The {.pkg bmm} package is needed for {.fun {name}}.")
  }
  if (!name %in% getNamespaceExports("bmm")) {
    cli::cli_abort(c(
      "The installed {.pkg bmm} ({packageVersion('bmm')}) does not export \\
       {.fun {name}}.",
      i = "The signal-detection models are not in a released {.pkg bmm} yet."
    ))
  }
  getExportedValue("bmm", name)
}

#' The class name an adapter is registered under
#' @noRd
adapter_classes <- function() {
  c("sdt_yn", "sdt_mafc", "ezdm", "ddm", "cswald", "mixture2p", "sdm")
}

#' @noRd
adapter_name <- function(model) {
  for (cls in adapter_classes()) {
    if (inherits(model, cls)) {
      return(cls)
    }
  }
  NA_character_
}

#' Look an adapter up by the model's class
#'
#' `NULL` when the model has none, including the four-parameter EZ
#' variant, whose column layout the adapter does not know.
#'
#' @noRd
generator_for <- function(model) {
  if (inherits(model, "ezdm_4par")) {
    return(NULL)
  }
  switch(adapter_name(model),
    sdt_yn = generate_sdt_yn,
    sdt_mafc = generate_sdt_mafc,
    ezdm = generate_ezdm,
    ddm = generate_ddm,
    cswald = generate_cswald,
    mixture2p = generate_mixture2p,
    sdm = generate_sdm,
    NULL
  )
}

#' Name generated columns after the model's own column names
#' @noRd
name_columns <- function(data, names) {
  names(data) <- unlist(names, use.names = FALSE)
  data
}

#' Yes/no signal detection: one row per stimulus class
#' @noRd
generate_sdt_yn <- function(pars, n_trials, model) {
  stimulus <- c(1, 0)
  counts <- bmm_fun("rsdt_yn")(
    2L, n_trials, stimulus,
    d = pars$d, criterion = pars$criterion, sdratio = pars$sdratio,
    dist = model$other_vars$dist
  )
  name_columns(
    data.frame(counts, stimulus, n_trials),
    list(
      model$resp_vars$response, model$other_vars$stimulus,
      model$other_vars$n_trials
    )
  )
}

#' m-alternative forced choice: one row of correct counts
#' @noRd
generate_sdt_mafc <- function(pars, n_trials, model) {
  counts <- bmm_fun("rsdt_mafc")(
    1L, n_trials,
    m = model$other_vars$m, d = pars$d, dist = model$other_vars$dist
  )
  name_columns(
    data.frame(counts, n_trials),
    list(model$resp_vars$response, model$other_vars$n_trials)
  )
}

#' EZ diffusion: one row of summary statistics
#' @noRd
generate_ezdm <- function(pars, n_trials, model) {
  out <- bmm::rezdm(
    1L, n_trials,
    drift = pars$drift, bound = pars$bound, ndt = pars$ndt, s = pars$s
  )
  name_columns(
    out[c("mean_rt", "var_rt", "n_upper", "n_trials")],
    list(
      model$resp_vars$mean_rt, model$resp_vars$var_rt,
      model$resp_vars$n_upper, model$other_vars$n_trials
    )
  )
}

#' Diffusion decision model: one row per trial
#' @noRd
generate_ddm <- function(pars, n_trials, model) {
  out <- bmm::rddm(
    n_trials,
    drift = pars$drift, bound = pars$bound, ndt = pars$ndt, zr = pars$zr
  )
  name_columns(
    out[c("rt", "response")],
    list(model$resp_vars$rt, model$resp_vars$response)
  )
}

#' Censored-shifted Wald: one row per trial
#'
#' `bmm::rcswald()` takes no `version` argument --- measured 2026-09-17 on
#' the installed 1.4.1.9000, on the study's pinned `develop` and in the
#' 1.3.2 source. It draws from `rtdists::rdiffusion()` with `a = bound`
#' and `z = zr * bound`, which is the `crisk` parameterisation, so the two
#' versions differ here in the mapping and not in the generator.
#'
#' **`simple` doubles `bound`.** Its `bound` is the distance from an
#' unbiased starting point to the correct boundary, half the separation
#' `rcswald()` takes; bmm's own `?cswald` says to multiply by 2 to get the
#' full separation. Measured 2026-09-17 on 20,000 draws generated with
#' `bound = 2`: profiled over `bound`, the `simple` likelihood peaks at
#' 0.995 and the `crisk` likelihood at 1.977. Passing `bound` straight
#' through would have made every `simple` recovery report a bias of
#' `log(2)` on the log link and read as a defect in bmm rather than in
#' this adapter.
#'
#' `zr` is 0.5 for `simple`, which has no starting-point parameter and is
#' defined against an unbiased start. `crisk` fixes `zr` at 0.5 too but a
#' design may free it, so its own value is passed.
#'
#' `sndt` is passed only when the model carries it: the fork's cswald has
#' it and both the pinned `develop` and released 1.3.2 have no `sndt`
#' anywhere (measured 2026-09-17), and `rcswald()` there has no such
#' argument to take. Reading it off the model rather than off the
#' installed version keeps one adapter right on all three.
#' @noRd
generate_cswald <- function(pars, n_trials, model) {
  crisk <- inherits(model, "cswald_crisk")
  args <- list(
    n = n_trials,
    drift = pars$drift,
    bound = if (crisk) pars$bound else 2 * pars$bound,
    ndt = pars$ndt,
    zr = if (crisk) pars$zr %||% 0.5 else 0.5,
    s = pars$s %||% 1
  )
  if (!is.null(pars$sndt)) {
    args$sndt <- pars$sndt
  }
  out <- do.call(bmm::rcswald, args)
  name_columns(
    out[c("rt", "response")],
    list(model$resp_vars$rt, model$resp_vars$response)
  )
}

#' Two-parameter mixture model: one response error per trial
#' @noRd
generate_mixture2p <- function(pars, n_trials, model) {
  y <- bmm::rmixture2p(
    n_trials,
    mu = pars$mu1, kappa = pars$kappa, p_mem = pars$thetat
  )
  name_columns(data.frame(y), list(model$resp_vars$resp_error))
}

#' Signal discrimination model: one response error per trial
#' @noRd
generate_sdm <- function(pars, n_trials, model) {
  y <- bmm::rsdm(n_trials, mu = pars$mu, c = pars$c, kappa = pars$kappa)
  name_columns(data.frame(y), list(model$resp_vars$resp_error))
}
