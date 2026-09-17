# The density adapters.
#
# The mirror of R/adapters.R. Where a generator adapter maps one subject's
# natural-scale parameters onto r<model>()'s argument names, a density
# adapter maps the same parameters onto d<model>()'s, and reads the
# observations out of the columns the model names. Column names come from
# the model object, never from literals, exactly as in the generator table.
#
# These exist only for `fit_ml(method = "optim")`. Decision 34 keeps the
# Stan route free of any per-model code --- it optimises bmm's own generated
# likelihood --- but the optim route has to know the likelihood in R, so it
# is capped at the six models of adapter_classes() for decision 13's own
# reason: a general package cannot re-implement every model it validates.
#
# Two things measured 2026-09-17 before this file was written
# (local/dev/sim/ml-optim-probe.R):
#
# 1. The six d<model>() functions do NOT agree on their `log` default ---
#    dmixture2p(), dsdm(), dsdt_yn() and dsdt_mafc() default to FALSE while
#    dddm() and dezdm() default to TRUE. Every adapter passes log = TRUE
#    explicitly; relying on the default would silently sum probabilities for
#    four of the six.
# 2. sdt_yn and sdt_mafc exist only in the bmm fork. CRAN bmm 1.3.2 exports
#    neither the models nor their densities, so those two adapters cannot be
#    exercised on a runner and their tests guard with skip_if_no_bmm_sdt().

#' Look a density adapter up by the model's class
#'
#' `NULL` when the model has none, including the four-parameter EZ variant,
#' whose column layout the adapter does not know --- the same cap
#' [generator_for()] applies.
#'
#' @param model A `bmmodel`.
#'
#' @return A function of `(pars, data, model)` returning a vector of log
#'   densities, or `NULL`.
#' @noRd
density_for <- function(model) {
  if (inherits(model, "ezdm_4par")) {
    return(NULL)
  }
  switch(adapter_name(model),
    sdt_yn = density_sdt_yn,
    sdt_mafc = density_sdt_mafc,
    ezdm = density_ezdm,
    ddm = density_ddm,
    mixture2p = density_mixture2p,
    sdm = density_sdm,
    NULL
  )
}

#' Yes/no signal detection: one row per stimulus class
#' @noRd
density_sdt_yn <- function(pars, data, model) {
  bmm::dsdt_yn(
    data[[model$resp_vars$response]],
    data[[model$other_vars$n_trials]],
    data[[model$other_vars$stimulus]],
    d = pars$d, criterion = pars$criterion, sdratio = pars$sdratio,
    dist = model$other_vars$dist, log = TRUE
  )
}

#' m-alternative forced choice: one row of correct counts
#' @noRd
density_sdt_mafc <- function(pars, data, model) {
  bmm::dsdt_mafc(
    data[[model$resp_vars$response]],
    data[[model$other_vars$n_trials]],
    m = model$other_vars$m, d = pars$d,
    dist = model$other_vars$dist, log = TRUE
  )
}

#' EZ diffusion: one row of summary statistics
#' @noRd
density_ezdm <- function(pars, data, model) {
  bmm::dezdm(
    data[[model$resp_vars$mean_rt]],
    data[[model$resp_vars$var_rt]],
    data[[model$resp_vars$n_upper]],
    data[[model$other_vars$n_trials]],
    drift = pars$drift, bound = pars$bound, ndt = pars$ndt, s = pars$s,
    log = TRUE
  )
}

#' Diffusion decision model: one row per trial
#' @noRd
density_ddm <- function(pars, data, model) {
  bmm::dddm(
    data[[model$resp_vars$rt]],
    data[[model$resp_vars$response]],
    drift = pars$drift, bound = pars$bound, ndt = pars$ndt, zr = pars$zr,
    log = TRUE
  )
}

#' Two-parameter mixture model: one response error per trial
#' @noRd
density_mixture2p <- function(pars, data, model) {
  bmm::dmixture2p(
    data[[model$resp_vars$resp_error]],
    mu = pars$mu1, kappa = pars$kappa, p_mem = pars$thetat, log = TRUE
  )
}

#' Signal discrimination model: one response error per trial
#' @noRd
density_sdm <- function(pars, data, model) {
  bmm::dsdm(
    data[[model$resp_vars$resp_error]],
    mu = pars$mu, c = pars$c, kappa = pars$kappa, log = TRUE
  )
}
