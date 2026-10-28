# Validation defaults: the per-model table a validation plan fills from
# (Milestone 26, D79 and D81, amended 2026-10-24).
#
# No adapter carries generating values (measured 2026-10-05), so a
# "default plan per model" needs a table of its own: the truth, the
# standard grid, the sampler budget, the SBC design point and its prior,
# the prior-check layout, the reference and the owner. bmmtools ships an
# entry per card as cards are made; a model without one gets the generic
# grid and must be given `pars`.

#' The stages of a validation, in the order a card shows them
#'
#' One row per stage (D80): its name, its title on the card, whether a
#' card needs it to be called complete, and the issue that builds it
#' where it does not exist yet. `reads` names the stage whose stored
#' result a view stage is computed from.
#'
#' @noRd
validation_stage_table <- function() {
  tibble::tibble(
    stage = c(
      "run_record", "prior_predictive", "sbc", "recovery", "computation",
      "identifiability", "individual_differences", "data_requirements",
      "agreement", "model_recovery", "posterior_checks"
    ),
    order = seq_len(11L),
    title = c(
      "Run record", "Prior predictive", "Implementation (SBC)",
      "Recovery", "Computation", "Identifiability",
      "Individual differences", "Data requirements", "Agreement",
      "Model recovery", "Posterior checks"
    ),
    required = c(
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE
    ),
    built = c(
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE
    ),
    issue = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 6L, 20L),
    reads = c(
      NA, NA, NA, NA, "recovery", "recovery", "recovery", "recovery",
      "recovery", NA, NA
    )
  )
}

#' The stage names, in card order
#' @noRd
validation_stage_names <- function() {
  validation_stage_table()$stage
}

#' The standard grid of D81: 20 and 50 subjects by 30, 100 and 300 trials
#' @noRd
standard_validation_grid <- function() {
  data.frame(
    n_subjects = rep(c(20L, 50L), times = 3L),
    n_trials = rep(c(30L, 100L, 300L), each = 2L)
  )
}

#' The criteria a plan judges its stages by, unless it is given others
#'
#' The defaults of spec 26 § Q2 (D80). `tolerance` is what
#' [compare_validations()] uses for a number without a Monte Carlo
#' standard error: by metric, with `default` for the rest. A tolerance of
#' 0 flags every change, so a comparison over-reports rather than hides;
#' timings never count as moved.
#'
#' @noRd
default_validation_criteria <- function() {
  list(
    prior_bound_rate_max = 0.10,
    sbc_ecdf_prob = 0.95,
    convergence_rate_min = 0.95,
    coverage_level = 0.95,
    coverage_alpha = 0.05,
    contraction_min = 0,
    z_mean_max = 2,
    read_off_r_low = 0.80,
    read_off_n_subjects = 50,
    tolerance = c(default = 0, seconds = Inf, median_fit_seconds = Inf)
  )
}

#' The entry a model without one of its own is planned with
#'
#' D81's grid, replications and sampler budget, the SBC design point under
#' bmm's default priors, and no truth: `pars` has to be given.
#'
#' @noRd
generic_validation_defaults <- function() {
  list(
    model_args = NULL,
    pars = NULL,
    sds = NULL,
    cors = NULL,
    grid = standard_validation_grid(),
    reps = 10L,
    sampler = list(
      chains = 4L, iter = 2000L, warmup = 1000L, cores = 4L,
      backend = "cmdstanr"
    ),
    sbc = list(
      n_subjects = 20L, n_trials = 50L, n_sims = 100L, prior = NULL,
      prior_reason = NULL, level = c("population", "sd"),
      sampler = list(chains = 2L, iter = 500L, backend = "cmdstanr"),
      generator = NULL, cores_per_fit = NULL, thin_ranks = NULL
    ),
    prior_check = list(
      n_subjects = 20L, n_trials = 50L, n_draws = 500L, range = NULL,
      range_reason = NULL, summary = NULL,
      sampler = list(chains = 2L, iter = 1000L, backend = "cmdstanr")
    ),
    reference = NULL,
    owner = NA_character_,
    version = NA_character_,
    sections = list()
  )
}

#' The entries bmmtools ships
#'
#' A prior is written as a plain data frame with the columns of
#' `brms::set_prior()` (`prior`, `class`, `coef`, `nlpar` or `dpar`), so
#' the table can be read and printed without brms; the runner turns it
#' into a `brmsprior` when it fits.
#'
#' @noRd
validation_entries <- function() {
  generic <- generic_validation_defaults()
  mixture2p <- generic
  mixture2p$model_args <- list(resp_error = "y")
  mixture2p$pars <- c(kappa = log(8), thetat = stats::qlogis(0.75))
  mixture2p$sds <- c(kappa = 0.3, thetat = 0.5)
  mixture2p$sbc$prior <- data.frame(
    prior = c("normal(2, 0.5)", "normal(0, 0.3)", "normal(0, 0.5)"),
    class = c("b", "sd", "sd"),
    coef = c("Intercept", "", ""),
    nlpar = c("kappa", "kappa", "thetat")
  )
  mixture2p$sbc$prior_reason <- paste(
    "bmm's default SD prior on log kappa, a half student_t(3, 0, 2.5),",
    "drew values (9.5 on 2 of 20 draws, measured 2026-09-16) at which",
    "rmixture2p() overflows the node stack, so the calibration runs under",
    "a narrower prior on the kappa intercept and both SDs."
  )
  mixture2p$owner <- "Gidon Frischkorn"
  mixture2p$version <- "0.1.0"

  # the generating values of Milestone 31's standard grid (eval doc
  # 2026-10-24 § 3): accuracy 1 / (1 + exp(-1.5)) = 0.82, boundary 1,
  # non-decision time 0.3 s
  ezdm <- generic
  ezdm$model_args <- list(
    mean_rt = "mean_rt", var_rt = "var_rt", n_upper = "n_upper",
    n_trials = "n_trials"
  )
  ezdm$pars <- c(drift = 1.5, bound = 0, ndt = log(0.3))
  ezdm$sds <- c(drift = 0.3, bound = 0.2, ndt = 0.2)
  terms <- names(ezdm$pars)
  ezdm$cors <- diag(3L)
  dimnames(ezdm$cors) <- list(terms, terms)
  ezdm$cors["drift", "bound"] <- ezdm$cors["bound", "drift"] <- 0.4
  # The calibration's prior (26.3, #21). Under bmm 1.4.1's ezdm defaults,
  # 8 calibration fits at 20 x 50: 2 failed the card's convergence
  # criterion (Rhat 2.13 and 1.06), 4 had divergences (up to 109) and 5
  # hit the maximum treedepth (up to 499 of 500), at parameter draws such
  # as a drift of -31 and an SD of 2.03 on log bound. The scales put the
  # entry's SDs (0.3, 0.2, 0.2) near the prior medians and give those
  # tail values negligible prior mass; the bound and ndt intercepts keep
  # bmm's priors. At 100 simulations (plan seed 2603) the ECDF check
  # passes but the convergence criterion still fails on 5 fits (4 at Rhat
  # above 1.05, up to 1.63; 2 at a rank ESS ratio below 0.5; one fit
  # both): the two worst drew bound intercepts of -1.302 and -1.281
  # (boundaries 0.27 and 0.28, under bmm's intercept prior) and saturated
  # the treedepth, and three fail at Rhat 1.06 to 1.07 with no treedepth
  # hit. The marginal-Rhat failures are probably a budget question (no
  # larger-budget run has tested it); the small-boundary saturation is
  # open, and so is the remedy.
  ezdm$sbc$prior <- data.frame(
    prior = c(
      "normal(0, 2)", "normal(0, 0.5)", "normal(0, 0.3)", "normal(0, 0.3)"
    ),
    class = c("Intercept", "sd", "sd", "sd"),
    coef = "",
    dpar = c("drift", "drift", "bound", "ndt")
  )
  ezdm$sbc$prior_reason <- paste(
    "Under bmm's default priors, a cauchy(0, 1) drift intercept and",
    "exponential SD priors, 2 of 8 calibration fits failed the card's",
    "convergence criterion (Rhat up to 2.13) and 4 had divergences (up to",
    "109), at parameter draws such as a drift of -31 and an SD of 2.03 on",
    "log bound (measured 2026-10-27). The calibration runs under a normal",
    "drift intercept and half-normal SDs, scaled so that the entry's SDs",
    "sit near the prior medians and those values get negligible prior",
    "mass; the bound and non-decision time intercepts keep bmm's priors."
  )
  # bmm 1.4.1's posterior_predict() gives mean RT only, and prior_check()
  # knows no range for ezdm (ARCHITECTURE § 26.0, amendment 4)
  ezdm$prior_check$range <- c(0.15, 5)
  ezdm$prior_check$range_reason <- paste(
    "Mean RT is judged against a plausibility range of 0.15 to 5 s, the",
    "owner's convention and not a bound of the response: bmm's predictions",
    "for ezdm give mean RT only, not accuracy."
  )
  ezdm$owner <- "Gidon Frischkorn"
  # 0.2.0: the narrowed SBC prior (0.1.0 shipped with 0.3.0 without one)
  ezdm$version <- "0.2.0"

  list(mixture2p = mixture2p, ezdm = ezdm)
}

#' Validation defaults for a model
#'
#' `r lifecycle::badge("experimental")`
#'
#' The generating values, the standard grid and the budget a
#' [validation_plan()] fills from when an argument is left `NULL`. An
#' entry exists for a model once its validation card is made; a model
#' without one, including a model of your own, is planned with the generic
#' grid and needs `pars`.
#'
#' @param model `NULL` for every entry; a model name (`"mixture2p"`) or a
#'   `bmmodel` for that model's.
#'
#' @return With `model = NULL`, a named list of entries. Otherwise the
#'   entry, or `NULL` with no error for a model without one. An entry is a
#'   list: `model_args` (the arguments its constructor is called with when
#'   [validation_plan()] is given the name), `pars`, `sds` and `cors` (the
#'   truth on the link scale), `grid` (`n_subjects` by `n_trials`),
#'   `reps`, `sampler` (passed to the fitter), `sbc` (the calibration's
#'   design point, number of simulations, prior and the reason that prior
#'   is not bmm's default, levels and sampler), `prior_check` (its design
#'   point, draws and sampler), `reference`, `owner`, `version` and
#'   `sections` (model-specific card sections).
#'
#' @details
#' The standard grid is 20 and 50 subjects by 30, 100 and 300
#' trials per subject, ten replications each, fitted with 4 chains of 2000
#' iterations (1000 warmup) under cmdstanr; SBC runs 100 simulations at 20
#' subjects by 50 trials with 2 chains of 500. Ten replications give a
#' coverage a Monte Carlo standard error of about 0.07 to 0.1, which the
#' card prints next to every rate.
#'
#' Two entries ship:
#' * `mixture2p`: `kappa = log(8)`, `thetat = qlogis(0.75)`, SDs 0.3 and
#'   0.5. Its calibration runs under a narrower prior than bmm's, because
#'   bmm's default SD prior on log `kappa` draws values at which
#'   `bmm::rmixture2p()` overflows the stack; the entry records the reason
#'   and the card prints it.
#' * `ezdm` (three parameters): `drift = 1.5`, `bound = 0`,
#'   `ndt = log(0.3)`, SDs 0.3, 0.2 and 0.2, a correlation of 0.4 between
#'   drift and bound. Its calibration also runs under a narrower prior
#'   than bmm's, because bmm's `cauchy(0, 1)` drift intercept and
#'   exponential SD priors draw parameter values (a drift of -31, an SD
#'   of 2 on log bound) at which some fits fail the card's convergence
#'   criterion.
#'
#' Each entry names its `owner`, who changes its values by pull request;
#' the card prints the owner and the plan version.
#'
#' @seealso [validation_plan()], which fills a plan from an entry.
#'
#' @examples
#' names(validation_defaults())
#' validation_defaults("ezdm")$pars
#' validation_defaults("sdm")
#'
#' @export
validation_defaults <- function(model = NULL) {
  entries <- validation_entries()
  if (is.null(model)) {
    return(entries)
  }
  name <- validation_model_name(model)
  entries[[name]]
}

#' The name a model is planned under
#'
#' The adapter's name where there is one, else the most specific class.
#'
#' @noRd
validation_model_name <- function(model, call = rlang::caller_env()) {
  is_name <- is.character(model) && length(model) == 1L && !is.na(model) &&
    nzchar(model)
  if (is_name) {
    return(model)
  }
  if (!inherits(model, "bmmodel")) {
    cli::cli_abort(
      "{.arg model} must be a {.cls bmmodel} or a model name, \\
       not {.obj_type_friendly {model}}.",
      call = call
    )
  }
  name <- adapter_name(model)
  if (is.na(name)) name <- class(model)[[length(class(model))]]
  name
}
