# The prior SD of each term from its prior's closed form (Milestone 22,
# D88). The table it returns is one of the three shapes recover()'s
# `prior_sd` takes; a prior-only fit and a prior check are the others.

#' Prior standard deviations from a fit's priors
#'
#' The standard deviation on the link scale that each parameter's prior
#' implies, for the priors that have one in closed form. Pass the result
#' as `prior_sd` to [recover()] or [recovery_grid()] to get the
#' contraction of each posterior against its prior without fitting the
#' prior alone.
#'
#' @param x A fit (a `brmsfit` or `bmmfit`, read through
#'   [brms::prior_summary()]), or a prior table with the columns of one:
#'   `prior`, `class`, `coef`, `group`, and `nlpar` or `dpar`, as
#'   `bmm::default_prior(formula, data, model)` returns before anything is
#'   fitted.
#' @param group The grouping factor whose standard deviations are wanted,
#'   as in [extract_estimates()].
#' @param ... Not used.
#'
#' @return A tibble with the columns `term`, `level` (`"population"`,
#'   `"effect"` or `"sd"`), `prior` (the prior as brms writes it) and
#'   `prior_sd_link`, one row per estimated parameter, named as
#'   [extract_estimates()] names it.
#'
#' @details
#' A coefficient whose own prior is empty takes the prior of its class,
#' as brms does. The closed forms are, for coefficients (class `b`):
#' `normal(mu, sigma)` gives `sigma`, `logistic(mu, s)` gives
#' `s * pi / sqrt(3)`, and `student_t(nu, mu, sigma)` gives
#' `sigma * sqrt(nu / (nu - 2))` for `nu > 2`. For standard deviations
#' (class `sd`), which brms bounds at 0, the prior is the half of a
#' distribution centred on 0, except `exponential(lambda)`, bmm's
#' default there, which gives `1 / lambda`. Half-normal gives
#' `sigma * sqrt(1 - 2 / pi)`
#' and half-Student-t, for `nu > 2`, the square root of
#' `sigma^2 * nu / (nu - 2)` minus its squared mean
#' `2 * sigma * sqrt(nu) * gamma((nu + 1) / 2) /
#' (sqrt(pi) * gamma(nu / 2) * (nu - 1))`.
#'
#' Any other prior is `NA` --- a flat prior, a prior with bounds on a
#' coefficient, a Student t with `nu <= 2` (whose variance is infinite),
#' a half prior not centred on 0, another family --- and one message names
#' those terms. A prior-only fit covers them: pass it, or a
#' [prior_check()], as `prior_sd` instead. Fixed parameters
#' (`constant()`) have no row, as they have no estimate.
#'
#' A Student t prior on a standard deviation with `nu = 3`, brms's
#' default, has an infinite fourth moment, so the SD of a prior-only
#' fit's draws converges slowly; the closed form is exact.
#'
#' @seealso [recover()] for contraction, which this feeds.
#'
#' @examplesIf rlang::is_installed("bmm")
#' prior_sd_table(bmm::default_prior(
#'   bmm::bmf(kappa ~ 1 + (1 | id), thetat ~ 1 + (1 | id)),
#'   data = data.frame(
#'     y = c(0.1, -0.2, 0.3, 0.05), id = factor(c(1, 1, 2, 2))
#'   ),
#'   model = bmm::mixture2p(resp_error = "y")
#' ))
#'
#' @export
prior_sd_table <- function(x, group = "id", ...) {
  rlang::check_dots_empty()
  priors <- prior_rows(x)
  priors$prior <- inherited_priors(priors)
  sd_rows <- priors$class == "sd" & priors$group == group
  pieces <- list(
    coefficient_prior_rows(priors[priors$class == "b", , drop = FALSE]),
    sd_prior_rows(priors[sd_rows, , drop = FALSE])
  )
  out <- dplyr::bind_rows(empty_prior_sd_table(), pieces)
  missing <- out$term[is.na(out$prior_sd_link)]
  if (length(missing) > 0L) {
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    what <- paste0(missing, " (", out$level[is.na(out$prior_sd_link)], ")")
    cli::cli_inform(c(
      "No closed-form prior SD for {.val {what}}.",
      i = "Their contraction will be {.val {NA}}; a prior-only fit, or \\
           {.fn prior_check}, covers them."
    ))
  }
  out
}

#' @noRd
empty_prior_sd_table <- function() {
  tibble::tibble(
    term = character(), level = character(), prior = character(),
    prior_sd_link = double()
  )
}

#' The prior rows of a fit or a prior table, as plain character columns
#' @noRd
prior_rows <- function(x, call = rlang::caller_env()) {
  if (inherits(x, "brmsfit")) {
    rlang::check_installed("brms", "to read a fit's priors.")
    x <- brms::prior_summary(x)
  }
  if (!is.data.frame(x)) {
    cli::cli_abort(
      "{.arg x} must be a fit or a prior table, not \\
       {.obj_type_friendly {x}}.",
      call = call
    )
  }
  missing <- setdiff(c("prior", "class", "coef", "group"), names(x))
  if (!any(c("nlpar", "dpar") %in% names(x))) missing <- c(missing, "nlpar")
  if (length(missing) > 0L) {
    cli::cli_abort(
      c(
        "{.arg x} is not a prior table.",
        x = "It lacks the column{?s} {.val {missing}}.",
        i = "{.fn brms::prior_summary} of a fit, or \\
             {.fn bmm::default_prior}, gives one."
      ),
      call = call
    )
  }
  column <- function(name) {
    if (name %in% names(x)) {
      out <- as.character(x[[name]])
      out[is.na(out)] <- ""
      out
    } else {
      rep("", nrow(x))
    }
  }
  nlpar <- column("nlpar")
  dpar <- column("dpar")
  tibble::tibble(
    prior = trimws(column("prior")),
    class = column("class"),
    coef = column("coef"),
    group = column("group"),
    par = ifelse(nzchar(nlpar), nlpar, dpar),
    lb = column("lb"),
    ub = column("ub")
  )
}

#' Each row's prior after brms's inheritance
#'
#' An empty prior on a coefficient row is the prior of the class row of
#' its parameter and group, then of its parameter, then of the class.
#'
#' @noRd
inherited_priors <- function(priors) {
  own <- priors$prior
  find <- function(i, coef, group, par) {
    at <- priors$class == priors$class[[i]] & priors$coef == coef &
      priors$group == group & priors$par == par & nzchar(own)
    if (any(at)) own[which(at)[[1L]]] else ""
  }
  vapply(seq_len(nrow(priors)), function(i) {
    if (nzchar(own[[i]])) {
      return(own[[i]])
    }
    for (step in list(
      c("", priors$group[[i]], priors$par[[i]]),
      c("", "", priors$par[[i]]),
      c("", "", "")
    )) {
      found <- find(i, step[[1L]], step[[2L]], step[[3L]])
      if (nzchar(found)) {
        return(found)
      }
    }
    ""
  }, character(1))
}

#' Population and effect rows from the class `b` coefficient rows
#' @noRd
coefficient_prior_rows <- function(rows) {
  rows <- rows[nzchar(rows$coef) & !is_constant_prior(rows$prior), ,
    drop = FALSE
  ]
  if (nrow(rows) == 0L) {
    return(NULL)
  }
  bounded <- nzchar(rows$lb) | nzchar(rows$ub)
  sd <- vapply(rows$prior, closed_form_sd, double(1), half = FALSE)
  sd[bounded] <- NA_real_
  tibble::tibble(
    term = coefficient_terms(rows$par, rows$coef),
    level = coefficient_kinds(rows$par, rows$coef),
    prior = rows$prior,
    prior_sd_link = unname(sd)
  )
}

#' SD rows from the class `sd` coefficient rows of one group
#' @noRd
sd_prior_rows <- function(rows) {
  rows <- rows[nzchar(rows$coef), , drop = FALSE]
  if (nrow(rows) == 0L) {
    return(NULL)
  }
  sd <- vapply(rows$prior, closed_form_sd, double(1), half = TRUE)
  tibble::tibble(
    term = coefficient_terms(rows$par, rows$coef),
    level = "sd",
    prior = rows$prior,
    prior_sd_link = unname(sd)
  )
}

#' @noRd
is_constant_prior <- function(prior) {
  grepl("^constant\\(", prior)
}

#' The SD a prior string implies, or `NA` without a closed form
#'
#' `half` reads the prior as the half of it above 0, which is what brms's
#' lower bound on a standard deviation makes of it; only a distribution
#' centred on 0 halves into a known form.
#'
#' @noRd
closed_form_sd <- function(prior, half) {
  parsed <- regmatches(prior, regexec("^([A-Za-z_]+)\\((.*)\\)$", prior))[[1L]]
  if (length(parsed) != 3L) {
    return(NA_real_)
  }
  family <- parsed[[2L]]
  args <- suppressWarnings(as.double(trimws(strsplit(parsed[[3L]], ",")[[1L]])))
  if (anyNA(args)) {
    return(NA_real_)
  }
  if (identical(family, "normal") && length(args) == 2L && args[[2L]] > 0) {
    if (!half) {
      return(args[[2L]])
    }
    if (args[[1L]] != 0) {
      return(NA_real_)
    }
    return(args[[2L]] * sqrt(1 - 2 / pi))
  }
  # positive support already, so its half is itself; only on an SD
  exponential <- identical(family, "exponential") && length(args) == 1L &&
    half && args[[1L]] > 0
  if (exponential) {
    return(1 / args[[1L]])
  }
  logistic <- identical(family, "logistic") && length(args) == 2L &&
    !half && args[[2L]] > 0
  if (logistic) {
    return(args[[2L]] * pi / sqrt(3))
  }
  student <- identical(family, "student_t") && length(args) == 3L &&
    is.finite(args[[1L]]) && args[[1L]] > 2 && args[[3L]] > 0
  if (student) {
    nu <- args[[1L]]
    sigma <- args[[3L]]
    if (!half) {
      return(sigma * sqrt(nu / (nu - 2)))
    }
    if (args[[2L]] != 0) {
      return(NA_real_)
    }
    mean <- 2 * sigma * sqrt(nu) * gamma((nu + 1) / 2) /
      (sqrt(pi) * gamma(nu / 2) * (nu - 1))
    return(sqrt(sigma^2 * nu / (nu - 2) - mean^2))
  }
  NA_real_
}
