# Subject-wise maximum likelihood, for comparison against a hierarchical fit.
#
# The point of this file is not inference. bmm discussion #263 records the
# wish that users not fit models this way; what it buys here is a *scored*
# contrast, so that the benefit of hierarchical estimation is a measured
# result rather than an assertion (ARCHITECTURE decision 34).
#
# Two things make it much smaller than it looks.
#
# 1. It is ONE fit, not one per subject. Under `p ~ 0 + id` every free
#    parameter is subject-specific and the rest are constant(), so the log
#    posterior is a sum of per-subject terms and the joint mode is exactly
#    the vector of per-subject modes. Measured over 40 subjects against
#    independent optimisation of bmm's own density, with no Stan involved:
#    the largest discrepancy was 3.0e-06 on the link scale and 2.5e-11 in
#    the per-subject log-likelihood (local/dev/sim/ml-factorization-check.R).
# 2. The per-subject estimates come back as *population* coefficients
#    (`b_kappa_id1`), because a no-pooling fit has no group-level effects at
#    all. So extract_estimates(level = "population") does the work --- and
#    drops bmm's constant() parameters on the way --- and this file only
#    renames the rows. extract_estimates(level = "subject") errors on such a
#    fit and is not the path.
#
# Decision 36: what comes back is an estimates tibble, not a fit, so
# extract_estimates() stays the only score-layer function touching a fit.

#' The no-pooling formula for a model
#'
#' `recovery_formula()` writes `p ~ 1 + (1 | id)`; this writes
#' `p ~ 0 + id`, the same parameters with the hierarchy taken out. It is
#' string construction on the same footing, and it recognises only the
#' shapes `recovery_formula()` itself writes (decision 38) --- anything
#' else is supplied ready-reduced through `formula`.
#'
#' @param model A `bmmodel`.
#' @param by The column defining one independent fit.
#' @param formula An already-reduced `bmmformula`, returned untouched.
#' @param call The calling environment, for error messages.
#'
#' @return A `bmmformula`.
#' @noRd
no_pooling_formula <- function(model, by = "id", formula = NULL,
                               call = rlang::caller_env()) {
  check_model(model, call = call)
  if (!is.null(formula)) {
    return(formula)
  }
  # `by = c("id", "condition")` has no Bayesian counterpart under the
  # current term naming: extract_estimates(level = "subject") names rows by
  # the plain group level and recover_subjects() joins on term + id, so an
  # ML id of "s01|A" would match nothing. Multi-column `by` belongs with a
  # cell-means formula whose *term* carries the condition, which is how the
  # task dimension already works.
  ok <- is.character(by) && length(by) == 1L && !is.na(by) && nzchar(by)
  if (!ok) {
    cli::cli_abort(
      c(
        "{.arg by} must be a single column name, \\
         not {.obj_type_friendly {by}}.",
        i = "Several columns need a cell-means formula whose term carries \\
             the condition; pass it through {.arg formula}."
      ),
      call = call
    )
  }
  rlang::check_installed("bmm", "to build a bmm formula.")
  free <- model_parameters(model)$free
  formulas <- lapply(free, function(p) {
    stats::as.formula(paste0(p, " ~ 0 + ", by), env = globalenv())
  })
  do.call(bmm::bmf, formulas)
}

#' A flat prior on every free parameter
#'
#' Decision 39: ML means flat priors. Under bmm's defaults the Laplace mode
#' is a MAP, and the comparison would become hierarchical shrinkage against
#' prior shrinkage --- a different and much weaker claim. An empty prior
#' string is what brms calls flat.
#'
#' Measured caveat, so the docs do not overstate it: on mixture2p at 100
#' trials a subject the default priors move the mode by a mean of 0.020 on
#' kappa and 0.069 on thetat. The default matters more at smaller trial
#' counts than at the scale the articles use.
#'
#' @noRd
flat_ml_prior <- function(model, call = rlang::caller_env()) {
  rlang::check_installed("brms", "to build a flat prior.")
  free <- model_parameters(model)$free
  priors <- lapply(free, function(p) {
    brms::set_prior("", class = "b", nlpar = p)
  })
  Reduce(`+`, priors)
}

#' Rename no-pooling population coefficients into subject rows
#'
#' The coefficients are named by `coefficient_terms()`: `<par>_<by><level>`
#' when the parameter has more than one coefficient, and the bare `<par>`
#' when it has one, which is what a single-subject fit gives.
#'
#' The map is built from the parameters and levels the caller already knows,
#' never by re-parsing the name. `split_coefficient()` splits on the *last*
#' underscore, so an id containing one (`"sub_01"` giving
#' `b_kappa_idsub_01`) would be read as the parameter `kappa_idsub`.
#'
#' @param est Population rows from [extract_estimates()].
#' @param free The model's free parameters.
#' @param levels The levels of `by`, as character.
#' @param by The grouping column's name.
#' @param max_abs_link The largest `abs(estimate)`, on the link scale, that
#'   still counts as converged.
#' @param call The calling environment, for error messages.
#'
#' @return A tibble of subject rows satisfying the estimates contract.
#' @noRd
ml_subject_rows <- function(est, free, levels, by, max_abs_link = 20,
                            call = rlang::caller_env()) {
  levels <- as.character(levels)
  map <- expand.grid(
    level = levels, par = free,
    stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE
  )
  map$name <- if (length(levels) == 1L) {
    map$par
  } else {
    paste0(map$par, "_", by, map$level)
  }

  hit <- match(est$term, map$name)
  unknown <- unique(est$term[is.na(hit)])
  if (length(unknown) > 0L) {
    cli::cli_abort(
      c(
        "{.fn fit_ml} did not recognise the coefficient{?s} {.val {unknown}}.",
        i = "A no-pooling fit of {.val {free}} by {.val {by}} should have \\
             one coefficient per parameter and level.",
        i = "Pass an already-reduced {.arg formula} for a design \\
             {.fn fit_ml} does not write itself."
      ),
      call = call
    )
  }
  # The map has to be total in both directions. Matching only the other way
  # leaves a coefficient the fit did not produce out of the result with no
  # error: the subject simply disappears, which is the outcome decision 40
  # exists to prevent --- and `check_estimator_balance()` would then be
  # comparing a table one subject short. A coefficient of exactly zero
  # posterior variance is dropped as a constant by `drop_constant_rows()`,
  # which a degenerate subject can trigger.
  absent <- setdiff(map$name, est$term)
  if (length(absent) > 0L) {
    cli::cli_abort(
      c(
        "{.fn fit_ml} has no estimate for the coefficient{?s} \\
         {.val {absent}}.",
        i = "A no-pooling fit of {.val {free}} by {.val {by}} has one \\
             coefficient per parameter and level.",
        i = "A coefficient with no posterior variance is dropped as a \\
             constant, which a degenerate subject can cause."
      ),
      call = call
    )
  }

  out <- est
  out$term <- map$par[hit]
  out$id <- map$level[hit]
  # brms keeps no optimiser return code on the fit --- `algorithm` is NULL
  # in the fit's stan_args --- so the range check is not an addition to
  # Stan's verdict, it is the only verdict available. It is also the more
  # reliable one: fitting one subject alone, Stan reported "relative
  # gradient magnitude is below tolerance" at a point 0.064 log-likelihood
  # short of the mode (local/dev/sim/ml-laplace-probe.out).
  ml_finalise(out, free, levels, max_abs_link, ci_method = "laplace")
}

#' The parts of an ML row that do not depend on the route
#'
#' Both routes return the same columns, the same class and the same
#' `estimator` (decision 37), so everything below the point estimate and
#' its interval is settled here rather than twice.
#'
#' `optim_ok` is the optimiser's own return code, which only the `optim`
#' route has: brms keeps none on a Laplace fit. Where it is absent
#' `converged` is the range check alone.
#'
#' @param out A data frame with at least `term`, `id`, `estimate`,
#'   `ci_low` and `ci_high`.
#' @param free The model's free parameters, for row order.
#' @param levels The levels of `by`, for row order.
#' @param max_abs_link The largest `abs(estimate)` on the link scale that
#'   still counts as converged.
#' @param ci_method `"laplace"` or `"wald"`.
#' @param optim_ok A logical vector, or `NULL` when the route has no
#'   optimiser flag.
#'
#' @return A tibble satisfying the estimates contract.
#' @noRd
ml_finalise <- function(out, free, levels, max_abs_link, ci_method,
                        optim_ok = NULL) {
  out$level <- "subject"
  out$ci_method <- ci_method
  out$estimator <- "ml"
  # Neither route produces MCMC draws. A Laplace fit is one chain of iid
  # draws from the normal approximation and the optim route has no draws at
  # all, so the values `posterior` would return describe the approximation
  # and not the estimate --- and check_convergence() would gate on them.
  out$rhat <- NA_real_
  out$ess_bulk <- NA_real_
  out$ess_tail <- NA_real_

  ok <- is.finite(out$estimate) & abs(out$estimate) <= max_abs_link
  if (!is.null(optim_ok)) {
    ok <- ok & optim_ok
  }
  out$converged <- ok
  # Decision 40: the row stays, so `n` differs visibly between estimators
  # and check_estimator_balance() warns. Dropping it would score ML on the
  # easiest subjects and the Bayesian estimator on all of them.
  out$estimate[!ok] <- NA_real_
  out$ci_low[!ok] <- NA_real_
  out$ci_high[!ok] <- NA_real_

  # ordered by the parameter and then the level as the data has them, not
  # by the character sort that would put id 10 before id 2
  order <- order(match(out$term, free), match(out$id, levels))
  out <- out[order, , drop = FALSE]
  tibble::as_tibble(out)
}

#' The objective one subject's optimisation minimises
#'
#' A density adapter returns log densities; what `optim()` needs is one
#' number. Wrapping here rather than asking the adapters for a sum keeps
#' the density table the exact mirror of the generator table, and makes
#' the user hook the more general of the two shapes: an `nll` that is not
#' a sum of per-row terms is still expressible.
#'
#' @param model A `bmmodel`.
#' @param nll A user objective, or `NULL` for the model's density adapter.
#' @param call The calling environment, for error messages.
#'
#' @return A function of `(pars, data, model)` returning one number.
#' @noRd
ml_objective <- function(model, nll = NULL, call = rlang::caller_env()) {
  if (!is.null(nll)) {
    if (!is.function(nll)) {
      cli::cli_abort(
        "{.arg nll} must be a function, not {.obj_type_friendly {nll}}.",
        call = call
      )
    }
    return(nll)
  }
  density <- density_for(model)
  if (is.null(density)) {
    cli::cli_abort(
      c(
        "{.fn fit_ml} has no density for {.cls {class(model)[1]}} with \\
         {.code method = \"optim\"}.",
        i = "It carries one for {.val {adapter_classes()}}.",
        i = "The {.code \"stan\"} method needs no density: it optimises \\
             bmm's own generated likelihood and works for any bmm model.",
        i = "Otherwise supply {.arg nll}."
      ),
      call = call
    )
  }
  function(pars, data, model) -sum(density(pars, data, model))
}

#' One subject's maximum-likelihood fit on the link scale
#'
#' The optimiser works on the link scale, which is the scale the estimates
#' contract stores, so the Hessian comes out on that scale too and the Wald
#' interval needs no transformation.
#'
#' A non-finite objective is returned as a large finite number rather than
#' `Inf`, because BFGS cannot take a step from `Inf` and would stop at the
#' start value while still reporting `convergence == 0`.
#'
#' @param y One subject's rows.
#' @param model A `bmmodel`.
#' @param objective As [ml_objective()] returns.
#' @param free The model's free parameters.
#' @param start Named start values on the link scale.
#'
#' @return The list [stats::optim()] returns, with `hessian`.
#' @noRd
ml_optim_one <- function(y, model, objective, free, start,
                         call = rlang::caller_env()) {
  eval_at <- function(theta) {
    objective(natural_pars(stats::setNames(theta, free), model), y, model)
  }
  # Evaluated once outside the guard below. A density that fails at an
  # extreme parameter is ordinary and must not stop the optimisation, but
  # an objective that cannot be called at all is a mistake in the call and
  # has to say so --- swallowed, it becomes a large finite number and the
  # fit "converges" at the start value.
  check_value <- function(v) {
    if (!is.numeric(v) || length(v) != 1L) {
      cli::cli_abort(
        c(
          "The objective must return a single number, \\
           not {.obj_type_friendly {v}}.",
          i = "{.arg nll} is called as {.code nll(pars, data, model)} and \\
               returns the negative log-likelihood of one subject."
        ),
        call = call
      )
    }
    v
  }
  check_value(eval_at(start))
  fn <- function(theta) {
    v <- tryCatch(eval_at(theta), error = function(e) NA_real_)
    # the type is checked at every step, not only at the start: an
    # objective that stops being a number once `optim()` moves away from
    # the start value would otherwise be indistinguishable from a density
    # that fails at an extreme parameter, which is what the sentinel is
    # for. `is.finite()` of a string is FALSE, so the two look alike.
    v <- check_value(v)
    if (!is.finite(v)) 1e100 else v
  }
  stats::optim(
    start, fn,
    method = "BFGS", control = list(reltol = 1e-14), hessian = TRUE
  )
}

#' Wald standard errors from an optimiser's Hessian
#'
#' The objective is a negative log-likelihood, so its Hessian at the mode
#' is the observed information and `solve()` of it is the covariance.
#'
#' A Hessian that cannot be inverted, or that has a non-positive variance
#' on the diagonal, yields `NA` standard errors. Gidon's decision
#' 2026-09-17: such a row keeps its estimate and its `converged` verdict
#' and loses only its interval, so it still counts in bias, RMSE and `r`
#' and drops out of `coverage` alone. Dropping the estimate would shrink
#' `n` on the ML side for a reason that is about the interval.
#'
#' @param hessian The Hessian [stats::optim()] returned.
#'
#' @return A numeric vector of standard errors, `NA` where unavailable.
#' @noRd
ml_wald_se <- function(hessian) {
  n <- nrow(hessian)
  v <- try(diag(solve(hessian)), silent = TRUE)
  if (inherits(v, "try-error")) {
    return(rep(NA_real_, n))
  }
  v[!is.finite(v) | v <= 0] <- NA_real_
  sqrt(v)
}

#' Subject-wise maximum likelihood by direct optimisation
#'
#' One `optim()` per subject on bmm's own R density, so no Stan code is
#' generated and no compiler is needed. This is the route that can run
#' anywhere, and the route whose convergence really is per subject: the
#' Stan route is one joint fit, in which a single boundary subject can
#' stall the whole optimisation and cannot be isolated.
#'
#' @inheritParams ml_optim_one
#' @param data The full data.
#' @param by The grouping column's name.
#' @param levels The levels of `by`, as character.
#' @param ci_level The interval level.
#' @param max_abs_link The convergence range on the link scale.
#'
#' @return A tibble satisfying the estimates contract.
#' @noRd
fit_ml_optim_rows <- function(model, data, by, levels, free, objective,
                              start, ci_level, max_abs_link) {
  # so that a missing density or a non-function `nll` aborts here rather
  # than inside the optimiser's error guard, which would turn it into a
  # large finite objective and a fit that "converged" at the start value
  force(objective)
  z <- stats::qnorm(1 - (1 - ci_level) / 2)
  n_par <- length(free)
  per_subject <- lapply(levels, function(lv) {
    y <- data[as.character(data[[by]]) == lv, , drop = FALSE]
    o <- ml_optim_one(y, model, objective, free, start)
    se <- ml_wald_se(o$hessian)
    data.frame(
      term = free, id = rep(lv, n_par), estimate = as.numeric(o$par),
      ci_low = as.numeric(o$par) - z * se,
      ci_high = as.numeric(o$par) + z * se,
      ci_level = rep(ci_level, n_par),
      optim_ok = rep(isTRUE(o$convergence == 0L), n_par),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, per_subject)
  optim_ok <- out$optim_ok
  out$optim_ok <- NULL
  ml_finalise(
    out, free, levels, max_abs_link,
    ci_method = "wald", optim_ok = optim_ok
  )
}

#' Start values for the optim route
#'
#' Zeros on the link scale. Measured over all 40 subjects of the 8.2 probe
#' (local/dev/sim/ml-optim-probe.R): a cold start at zero reaches the same
#' optimum as a start at the Stan mode to 4.3e-07 on both parameters, with
#' `optim()$convergence == 0` throughout, so there is nothing to warm up
#' from and `start` exists only for a model where that stops being true.
#'
#' @noRd
check_ml_start <- function(start, free, call = rlang::caller_env()) {
  if (is.null(start)) {
    return(stats::setNames(rep(0, length(free)), free))
  }
  bad <- !is.numeric(start) || is.null(names(start)) || anyNA(start)
  if (bad || !setequal(names(start), free)) {
    cli::cli_abort(
      c(
        "{.arg start} must be a named numeric vector of link-scale start \\
         values, one per free parameter.",
        i = "This model's free parameters are {.val {free}}."
      ),
      call = call
    )
  }
  start[free]
}

#' Per-cell status of an ML fit
#'
#' `count_converged()` counts *replications*, so for 50 subjects in one
#' replication with 3 failures it returns 1 --- true and useless. This is
#' the per-cell table the print method and the article read instead.
#'
#' @noRd
ml_cells_table <- function(rows) {
  dplyr::summarise(
    rows,
    n_terms = dplyr::n(),
    n_converged = sum(.data$converged %in% TRUE),
    converged = all(.data$converged %in% TRUE),
    .by = "id"
  )
}

#' Subject-wise maximum likelihood estimates of a bmm model
#'
#' Fits a model with no pooling --- one independent set of parameters per
#' subject --- and returns the estimates in the same tibble shape a
#' hierarchical fit gives, so that the two can be scored against one truth
#' and compared. It exists to measure what hierarchical estimation buys,
#' not to recommend maximum likelihood for inference.
#'
#' @section Two routes:
#' `method = "stan"` optimises bmm's own generated likelihood with
#' `algorithm = "laplace"`. It re-implements nothing, so it works for any
#' bmm model, and it needs a working CmdStan.
#'
#' `method = "optim"` maximises the likelihood directly with
#' [stats::optim()] on bmm's R density for the model. It needs no
#' compiler and is roughly twenty times faster --- 0.35 s against 8 s for
#' 40 subjects of `mixture2p` --- but it is capped at the models bmmtools
#' carries a density for, listed in the error it raises otherwise, and
#' `nll` is the escape hatch for anything else.
#'
#' The two routes return the same columns, the same class and the same
#' `estimator`, and differ only in `ci_method`. Measured on 40 subjects of
#' `mixture2p` at 100 trials, they agree to Monte Carlo error on the
#' estimate and their Wald and Laplace standard errors agree to a mean
#' ratio of 1.000 (range 0.95 to 1.05), giving the same `coverage` and
#' `calibration_slope` to three decimals. The route is not a moderator of
#' a scored number.
#'
#' @section One fit, not one per subject (`method = "stan"`):
#' Under `p ~ 0 + <by>` every free parameter belongs to one subject and
#' bmm's remaining parameters are constants, so the log posterior is a sum
#' of per-subject terms and the joint mode is the vector of per-subject
#' modes. The Stan route therefore runs a single fit for the whole data
#' set. This was checked rather than assumed: over 40 subjects, the joint
#' Stan mode and independent optimisation of bmm's own density agreed to
#' 3.0e-06 on the link scale.
#'
#' The cost is that convergence is no longer separable by subject --- one
#' subject at a boundary can stall the whole optimisation --- so the
#' per-cell verdict in `ml_cells` is the range check described below. The
#' `optim` route does not have this property: it is one independent
#' optimisation per subject.
#'
#' @section What counts as converged:
#' `converged` is `TRUE` where the estimate is finite and
#' `abs(estimate) <= max_abs_link` on the link scale, and, on the `optim`
#' route, where `optim()` also returned its success code. brms keeps no
#' optimiser return code on a Laplace fit, so on the Stan route the range
#' check is the only verdict available --- and not a weak one: fitting a
#' single subject, Stan reported convergence at a point 0.064
#' log-likelihood short of the mode.
#'
#' A row that fails is kept with `estimate = NA` rather than dropped,
#' because `metric_bias()` and `metric_rmse()` drop incomplete pairs
#' silently: dropping would compare the Bayesian estimator on every
#' subject against maximum likelihood on the easiest ones and report the
#' difference as the estimator's. Report `n` alongside any comparison.
#'
#' On the `optim` route a Hessian that cannot be inverted gives a row with
#' no interval: `ci_low` and `ci_high` are `NA` while the estimate and
#' `converged` stand, so the subject counts in `bias`, `rmse` and `r` and
#' drops out of `coverage` only.
#'
#' `rhat`, `ess_bulk` and `ess_tail` are always `NA`. Neither route
#' produces MCMC draws, so the values `posterior` would return describe
#' the approximation and not the estimate.
#'
#' @param model A `bmmodel`, as [bmm::mixture2p()] returns.
#' @param data The data to fit.
#' @param formula An already-reduced `bmmformula` with no group-level
#'   terms. `NULL`, the default, writes `p ~ 0 + <by>` for every free
#'   parameter. Supply one for a design `fit_ml()` does not write itself.
#'   `method = "optim"` builds no formula and does not accept one.
#' @param by The column defining one independent set of parameters. A
#'   single column in this version.
#' @param method `"stan"` optimises bmm's own generated likelihood through
#'   `algorithm = "laplace"`, re-implementing nothing, and needs CmdStan.
#'   `"optim"` maximises bmm's R density for the model directly with
#'   [stats::optim()] and needs no compiler. See the two-routes section.
#' @param prior `"flat"`, the default, puts an empty prior on every free
#'   parameter so the mode is the maximum-likelihood estimate.
#'   `"default"` keeps bmm's priors, which makes the mode a penalised
#'   estimate rather than an ML one. `method = "optim"` has no priors at
#'   all and accepts `"flat"` only.
#' @param ci_level The interval level. The interval is the Laplace
#'   approximation's on the Stan route (`ci_method = "laplace"`) and a
#'   Wald interval from the optimiser's Hessian on the other
#'   (`ci_method = "wald"`).
#' @param draws How many draws to take from the Laplace approximation.
#'   `method = "stan"` only.
#' @param max_abs_link The largest `abs(estimate)` on the link scale that
#'   still counts as converged.
#' @param start Named link-scale start values for `method = "optim"`, one
#'   per free parameter. `NULL`, the default, starts at zero, which was
#'   measured to reach the same optimum as a start at the Stan mode on
#'   every subject of the reference run.
#' @param nll A function of `(pars, data, model)` returning the negative
#'   log-likelihood of one subject's `data` as a single number, where
#'   `pars` is a named list on the natural scale with the model's fixed
#'   parameters included. Supply it for a model bmmtools carries no
#'   density for. `method = "optim"` only.
#' @param file,refit Passed to [fit_cached()]. `NULL` uses a temporary
#'   file, so nothing is cached between sessions. `method = "optim"` does
#'   not cache and ignores both: it produces no fit object, and at well
#'   under a second per call a cache would add an invalidation surface
#'   without saving anything.
#' @param seed Passed to the fitter. `method = "optim"` is deterministic
#'   and ignores it.
#' @param ... Passed to [fit_cached()] and on to the fitter.
#' @param .fitter As in [fit_cached()]. `method = "stan"` only.
#'
#' @return A `bmmtools_ml`: a tibble satisfying the estimates contract,
#'   one row per subject and term, with `level = "subject"` and
#'   `estimator = "ml"`. `attr(x, "ml_cells")` holds one row per subject
#'   with its term and convergence counts.
#'
#' @seealso [recover_subjects()] to score the result against a known
#'   truth, and [extract_estimates()] for the hierarchical side of the
#'   comparison. `bind_rows()` of the two is the join; there is no
#'   separate comparison function. [recovery_grid()] with `ml = TRUE`
#'   runs both fits on every cell of a design and does the join itself.
#'
#' @examples
#' \dontrun{
#' model <- bmm::mixture2p(resp_error = "y")
#' sim <- simulate_recovery(
#'   model,
#'   pars = c(kappa = 2, thetat = 1), sds = c(kappa = 0.3, thetat = 0.5),
#'   n_subjects = 40, n_trials = 100, seed = 1
#' )
#' ml <- fit_ml(model, sim$data)
#' summary(recover_subjects(ml, sim$truth$subjects))
#'
#' # the same thing without a compiler, and without Stan
#' ml2 <- fit_ml(model, sim$data, method = "optim")
#' }
#'
#' @export
fit_ml <- function(model,
                   data,
                   formula = NULL,
                   by = "id",
                   method = c("stan", "optim"),
                   prior = c("flat", "default"),
                   ci_level = 0.95,
                   draws = 1000,
                   max_abs_link = 20,
                   start = NULL,
                   nll = NULL,
                   file = NULL,
                   refit = c("on_change", "never", "always"),
                   seed = NULL,
                   ...,
                   .fitter = NULL) {
  check_model(model)
  method <- rlang::arg_match(method)
  prior <- rlang::arg_match(prior)
  refit <- rlang::arg_match(refit)
  if (!is.data.frame(data)) {
    cli::cli_abort(
      "{.arg data} must be a data frame, not {.obj_type_friendly {data}}."
    )
  }
  if (!is.numeric(max_abs_link) || length(max_abs_link) != 1L ||
        is.na(max_abs_link) || max_abs_link <= 0) {
    cli::cli_abort(
      "{.arg max_abs_link} must be a single positive number, \\
       not {.obj_type_friendly {max_abs_link}}."
    )
  }

  free <- model_parameters(model)$free
  optim <- identical(method, "optim")

  # Arguments that belong to the other route are named errors rather than
  # silent no-ops, because each of them would otherwise change nothing
  # while looking as though it had: a `prior` that is not applied, a
  # `formula` that is not read, a `.fitter` that is never called.
  if (optim) {
    if (!identical(prior, "flat")) {
      cli::cli_abort(c(
        "{.code method = \"optim\"} maximises the likelihood and has no \\
         priors, so {.code prior = \"default\"} has no meaning for it.",
        i = "Use {.code method = \"stan\"} for a penalised estimate under \\
             bmm's default priors."
      ))
    }
    if (!is.null(formula)) {
      cli::cli_abort(c(
        "{.code method = \"optim\"} fits each level of {.arg by} \\
         separately and builds no formula, so it cannot use {.arg formula}.",
        i = "Use {.code method = \"stan\"} for a design written as a formula."
      ))
    }
    if (!is.null(.fitter)) {
      cli::cli_abort(c(
        "{.arg .fitter} is the seam of {.code method = \"stan\"}; \\
         {.code method = \"optim\"} calls no fitter.",
        i = "Pass {.arg nll} to stand in for the model's density."
      ))
    }
    if (!missing(draws)) {
      cli::cli_abort(c(
        "{.arg draws} is the size of the Laplace sample of \\
         {.code method = \"stan\"}; {.code method = \"optim\"} draws \\
         nothing.",
        i = "Its intervals are Wald intervals from the Hessian."
      ))
    }
    # the optim route calls no fitter, so nothing consumes the dots. Left
    # unchecked, a mistyped argument reaches nothing and says nothing ---
    # `ci_levl = 0.8` would score 95% intervals in silence.
    rlang::check_dots_empty()
    group <- by
    if (!is.character(group) || length(group) != 1L || is.na(group) ||
          !nzchar(group)) {
      cli::cli_abort(
        "{.arg by} must be a single column name, \\
         not {.obj_type_friendly {by}}."
      )
    }
  } else {
    if (!is.null(nll)) {
      cli::cli_abort(c(
        "{.arg nll} is the density hook of {.code method = \"optim\"}.",
        i = "{.code method = \"stan\"} optimises bmm's own generated \\
             likelihood and never evaluates a density in R."
      ))
    }
    # `algorithm` and `backend` are not formals of fit_ml(), so a user's
    # value lands in the dots and is spliced into the fitter call beside
    # the one fit_ml() passes. The fitter then raises R's bare "matched by
    # multiple actual arguments", and cache_key() records the first of the
    # two --- a key for a fit nobody asked for.
    owned <- intersect(names(rlang::list2(...)), c("algorithm", "backend"))
    if (length(owned) > 0L) {
      cli::cli_abort(c(
        "{.fn fit_ml} sets {.arg {owned}} itself and {?it/they} cannot be \\
         passed through.",
        i = "A subject-wise ML fit is {.val laplace} on the \\
             {.val cmdstanr} backend; no other combination gives a mode."
      ))
    }
    f <- no_pooling_formula(model, by = by, formula = formula)
    group <- all.vars(f[[1L]][[3L]])[1L]
  }

  if (!group %in% names(data)) {
    cli::cli_abort("{.arg data} has no column {.val {group}}.")
  }
  levels <- levels(as.factor(data[[group]]))
  if (length(levels) == 0L) {
    cli::cli_abort(c(
      "{.arg data} has no non-missing value in {.val {group}}.",
      i = "{.fn fit_ml} fits one set of estimates per level of \\
           {.arg by}, and there is none."
    ))
  }

  if (optim) {
    rows <- fit_ml_optim_rows(
      model, data, group, levels, free,
      objective = ml_objective(model, nll),
      start = check_ml_start(start, free),
      ci_level = ci_level, max_abs_link = max_abs_link
    )
  } else {
    ml_prior <- if (identical(prior, "flat")) flat_ml_prior(model) else NULL
    fit <- fit_cached(
      formula = f,
      data = data,
      model = model,
      file = file %||% tempfile(pattern = "bmmtools-ml-"),
      prior = ml_prior,
      refit = refit,
      algorithm = "laplace",
      backend = "cmdstanr",
      draws = draws,
      seed = seed,
      ...,
      .fitter = .fitter
    )
    est <- extract_estimates(
      fit,
      level = "population",
      ci_level = ci_level,
      converged = TRUE,
      estimator = "ml"
    )
    rows <- ml_subject_rows(
      est, free, levels, group,
      max_abs_link = max_abs_link
    )
  }

  new_bmmtools_ml(
    rows,
    cells = ml_cells_table(rows),
    method = method,
    prior = prior,
    by = group
  )
}

#' @noRd
new_bmmtools_ml <- function(x, cells, method, prior, by) {
  x <- tibble::as_tibble(x)[names(estimates_contract())]
  structure(
    x,
    ml_cells = tibble::as_tibble(cells),
    ml_method = method,
    ml_prior = prior,
    ml_by = by,
    class = c("bmmtools_ml", class(x))
  )
}

#' @noRd
#' @importFrom dplyr dplyr_reconstruct
#' @exportS3Method dplyr::dplyr_reconstruct
dplyr_reconstruct.bmmtools_ml <- function(data, template) {
  if (!all(names(estimates_contract()) %in% names(data))) {
    return(tibble::as_tibble(data))
  }
  NextMethod()
}

#' Subset an ML estimates object
#'
#' As for [`[.bmmtools_recovery`]: `dplyr::select()` on a tibble subclass
#' subsets through `[` and never reaches `dplyr_reconstruct()`, so this
#' method is what drops the class once a contract column is gone.
#'
#' @param x A `bmmtools_ml` object.
#' @param ... Passed to the tibble method.
#' @return An ML object while the contract holds, a plain tibble once it
#'   does not.
#' @export
`[.bmmtools_ml` <- function(x, ...) {
  demote_if_incomplete(NextMethod(), names(estimates_contract()))
}

#' Format and print subject-wise ML estimates
#'
#' The header says how many cells converged without being asked, because
#' the difference between noticing three failed subjects and not noticing
#' them is the difference between a comparison and a wrong number.
#'
#' @param x A `bmmtools_ml` object.
#' @param ... Not used.
#'
#' @return `format()` returns a character vector; `print()` returns `x`
#'   invisibly.
#'
#' @export
format.bmmtools_ml <- function(x, ...) {
  cells <- attr(x, "ml_cells")
  n <- if (is.null(cells)) 0L else nrow(cells)
  ok <- if (is.null(cells)) 0L else sum(cells$converged %in% TRUE)
  terms <- unique(x$term)
  method <- attr(x, "ml_method") %||% "stan"
  # read off the rows rather than assumed from the method, so that a
  # bind_rows() of the two routes prints what it actually holds
  intervals <- unique(x$ci_method[!is.na(x$ci_method)])
  if (length(intervals) == 0L) {
    intervals <- if (identical(method, "optim")) "wald" else "laplace"
  }
  header <- c(
    "<bmmtools_ml>",
    paste0(
      n, " cell", if (n != 1L) "s", " by ", attr(x, "ml_by") %||% "id",
      ", ", ok, " converged."
    ),
    paste0(
      "Method: ", method, "/", paste(intervals, collapse = "/"), ", ",
      if (identical(method, "optim")) {
        "no priors (maximum likelihood)."
      } else {
        paste0(attr(x, "ml_prior") %||% "flat", " priors.")
      }
    ),
    paste0(
      length(terms), " parameter", if (length(terms) != 1L) "s", ": ",
      paste(terms, collapse = ", "), "."
    )
  )
  if (ok < n) {
    header <- c(header, paste0(
      "Failed cells keep their row with estimate = NA; report n when ",
      "comparing estimators."
    ))
  }
  c(
    header, "",
    utils::capture.output(print(tibble::as_tibble(x)))
  )
}

#' @rdname format.bmmtools_ml
#' @export
print.bmmtools_ml <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}
