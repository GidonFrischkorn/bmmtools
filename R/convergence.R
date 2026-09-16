# The convergence gate.
#
# Thresholds are arguments with the seeds' defaults. The numbers come
# from posterior::summarise_draws() and brms::nuts_params(), the only two
# sources of convergence numbers in the package, and the logic lives in
# convergence_from_parts() so that it can be tested on hand-built tables
# without a fit.

#' Validate the numeric thresholds of the gate
#' @noRd
check_thresholds <- function(rhat_max = 1.05,
                             ess_bulk_min = 400,
                             ess_tail_min = NULL,
                             divergent_max = 10,
                             call = rlang::caller_env()) {
  check_one <- function(value, name, allow_null = FALSE) {
    if (is.null(value) && allow_null) {
      return(invisible(NULL))
    }
    if (!is.numeric(value) || length(value) != 1L || is.na(value)) {
      cli::cli_abort(
        "{.arg {name}} must be a single number, \\
         not {.obj_type_friendly {value}}.",
        call = call
      )
    }
    invisible(NULL)
  }
  check_one(rhat_max, "rhat_max")
  check_one(ess_bulk_min, "ess_bulk_min")
  check_one(ess_tail_min, "ess_tail_min", allow_null = TRUE)
  check_one(divergent_max, "divergent_max")
  list(
    rhat_max = rhat_max,
    ess_bulk_min = ess_bulk_min,
    ess_tail_min = ess_tail_min,
    divergent_max = divergent_max
  )
}

#' The sampler's maximum tree depth, read from the fit
#'
#' Measured on the fixture (cmdstanr backend, brms 2.23.0):
#' `fit$fit@stan_args[[1]]$control$max_treedepth` is present. When it is
#' not, Stan's default of 10 is used.
#'
#' @noRd
fit_treedepth_max <- function(fit) {
  args <- tryCatch(fit$fit@stan_args[[1L]], error = function(e) NULL)
  value <- args$control$max_treedepth
  if (is.numeric(value) && length(value) == 1L && !is.na(value)) {
    return(as.double(value))
  }
  10
}

#' Count a NUTS diagnostic over the post-warmup iterations
#' @noRd
count_nuts <- function(nuts, parameter, test) {
  if (is.null(nuts)) {
    return(NA_integer_)
  }
  values <- nuts$Value[nuts$Parameter == parameter]
  if (length(values) == 0L) {
    return(NA_integer_)
  }
  as.integer(sum(test(values), na.rm = TRUE))
}

#' Apply the gate to a diagnostics table and a NUTS table
#'
#' `diagnostics` has the columns `variable`, `rhat`, `ess_bulk`,
#' `ess_tail` (what `posterior::summarise_draws()` returns with the
#' default convergence measures); `nuts` has `Parameter` and `Value`
#' (what `brms::nuts_params()` returns), or is `NULL` when the fit has no
#' sampler diagnostics.
#'
#' Rows with a missing rhat are dropped before anything is computed: a
#' constant has `NA` diagnostics and must neither fail the gate nor
#' count as assessed. Tree-depth hits are counted and not gated, because
#' hitting the maximum is an efficiency warning rather than a validity
#' failure, and the seeds never gated on it.
#'
#' @noRd
convergence_from_parts <- function(diagnostics, nuts, treedepth_max,
                                   thresholds) {
  keep <- !is.na(diagnostics$rhat)
  assessed <- diagnostics[keep, , drop = FALSE]
  n_variables <- nrow(assessed)

  max_or_na <- function(x) if (length(x) == 0L) NA_real_ else max(x)
  min_or_na <- function(x) if (length(x) == 0L) NA_real_ else min(x)

  max_rhat <- max_or_na(assessed$rhat)
  min_ess_bulk <- min_or_na(assessed$ess_bulk[!is.na(assessed$ess_bulk)])
  min_ess_tail <- min_or_na(assessed$ess_tail[!is.na(assessed$ess_tail)])
  n_divergent <- count_nuts(nuts, "divergent__", function(v) v > 0)
  n_max_treedepth <- count_nuts(
    nuts, "treedepth__", function(v) v >= treedepth_max
  )

  if (n_variables == 0L) {
    pass <- NA
    failed <- NA_character_
  } else {
    checks <- c(
      rhat = isTRUE(max_rhat > thresholds$rhat_max),
      ess_bulk = isTRUE(min_ess_bulk < thresholds$ess_bulk_min),
      ess_tail = !is.null(thresholds$ess_tail_min) &&
        isTRUE(min_ess_tail < thresholds$ess_tail_min),
      divergent = isTRUE(n_divergent > thresholds$divergent_max)
    )
    pass <- !any(checks)
    failed <- paste(names(checks)[checks], collapse = ", ")
  }

  out <- tibble::tibble(
    max_rhat = max_rhat,
    min_ess_bulk = min_ess_bulk,
    min_ess_tail = min_ess_tail,
    n_divergent = n_divergent,
    n_max_treedepth = n_max_treedepth,
    n_variables = as.integer(n_variables),
    pass = pass,
    failed = failed
  )
  attr(out, "thresholds") <- c(thresholds, list(treedepth_max = treedepth_max))
  out
}

#' Check whether a fit passes the convergence gate
#'
#' One row of worst-case diagnostics and a verdict. The thresholds are
#' arguments with the defaults the validation scripts in bmm used: rhat
#' at most 1.05, bulk ESS at least 400, at most ten divergent
#' transitions. A fit that fails is meant to be **scored but flagged**,
#' never dropped: [extract_estimates()] carries the verdict as its
#' `converged` column and `summary()` of a recovery object counts it as
#' `n_converged`.
#'
#' @param fit A `brmsfit`, and so also a `bmmfit`.
#' @param rhat_max Largest acceptable rhat over the variables considered.
#' @param ess_bulk_min Smallest acceptable bulk effective sample size.
#' @param ess_tail_min Smallest acceptable tail effective sample size.
#'   `NULL`, the default, reports it without gating on it.
#' @param divergent_max Largest acceptable number of divergent
#'   transitions after warmup, summed over chains.
#' @param treedepth_max The sampler's maximum tree depth, used to count
#'   the post-warmup iterations that hit it. `NULL` reads it from the fit
#'   and falls back to Stan's default of 10. Tree-depth hits are reported
#'   and never gated on: hitting the maximum is an efficiency warning,
#'   not a validity failure.
#' @param variables Names of the draws variables to consider. `NULL`
#'   means every variable `posterior::summarise_draws()` returns.
#'   Variables with a missing rhat, which is what a parameter fixed to a
#'   constant has, are dropped rather than failed on.
#' @param ... Passed to methods.
#'
#' @return A one-row tibble with `max_rhat`, `min_ess_bulk`,
#'   `min_ess_tail`, `n_divergent`, `n_max_treedepth`, `n_variables`,
#'   `pass` and `failed`, the last a comma-separated list of the
#'   criteria that failed (`""` when none did). `pass` is `NA`, not
#'   `TRUE`, when no variable had a finite rhat. The thresholds used are
#'   stored in the attribute `thresholds`.
#'
#' @details
#' The numbers come from `posterior::summarise_draws()` and
#' `brms::nuts_params()`. A fit without sampler diagnostics (a
#' variational fit, or one whose sampler parameters were not saved)
#' leaves `n_divergent` and `n_max_treedepth` as `NA` and does not fail
#' the divergence criterion, so a reader can see it was not checked.
#'
#' @examples
#' \dontrun{
#' fit <- bmm::bmm(bmm::bmf(kappa ~ 1, thetat ~ 1), data, model)
#' check_convergence(fit)
#' check_convergence(fit, rhat_max = 1.01, ess_tail_min = 400)
#' }
#'
#' @export
check_convergence <- function(fit, ...) {
  UseMethod("check_convergence")
}

#' @rdname check_convergence
#' @export
check_convergence.default <- function(fit, ...) {
  cli::cli_abort(
    "{.arg fit} must be a {.cls brmsfit}, not {.obj_type_friendly {fit}}."
  )
}

#' @rdname check_convergence
#' @export
check_convergence.brmsfit <- function(fit,
                                      rhat_max = 1.05,
                                      ess_bulk_min = 400,
                                      ess_tail_min = NULL,
                                      divergent_max = 10,
                                      treedepth_max = NULL,
                                      variables = NULL,
                                      ...) {
  rlang::check_installed("brms", "to check convergence of a fit.")
  rlang::check_dots_empty()
  thresholds <- check_thresholds(
    rhat_max, ess_bulk_min, ess_tail_min, divergent_max
  )
  if (is.null(treedepth_max)) {
    treedepth_max <- fit_treedepth_max(fit)
  }
  bad_treedepth <- !is.numeric(treedepth_max) ||
    length(treedepth_max) != 1L || is.na(treedepth_max)
  if (bad_treedepth) {
    cli::cli_abort(
      "{.arg treedepth_max} must be a single number or {.code NULL}, \\
       not {.obj_type_friendly {treedepth_max}}."
    )
  }

  convergence_from_fit(
    fit, posterior::as_draws_array(fit), thresholds, treedepth_max, variables
  )
}

#' The gate on draws already extracted from a fit
#'
#' Split out so that [extract_estimates()] extracts the draws once and
#' hands them to both the gate and the summaries.
#'
#' @noRd
convergence_from_fit <- function(fit, draws, thresholds, treedepth_max,
                                 variables = NULL) {
  if (!is.null(variables)) {
    if (!is.character(variables)) {
      cli::cli_abort(
        "{.arg variables} must be a character vector, \\
         not {.obj_type_friendly {variables}}."
      )
    }
    missing <- setdiff(variables, posterior::variables(draws))
    if (length(missing) > 0L) {
      cli::cli_abort(
        "Variable{?s} {.val {missing}} {?is/are} not in the fit."
      )
    }
    draws <- posterior::subset_draws(draws, variable = variables)
  }

  diagnostics <- posterior::summarise_draws(
    draws, posterior::default_convergence_measures()
  )
  nuts <- tryCatch(brms::nuts_params(fit), error = function(e) NULL)

  convergence_from_parts(diagnostics, nuts, treedepth_max, thresholds)
}
