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
#' coefficient (bounds set on the class bound each of its coefficients,
#' as in brms), a Student t with `nu <= 2` (whose variance is infinite),
#' a half prior not centred on 0, another family --- and one message names
#' those terms. A prior-only fit covers them: pass it, or a
#' [prior_check()], as `prior_sd` instead. Fixed parameters
#' (`constant()`) have no row, as they have no estimate.
#'
#' A distributional parameter with an intercept and other population
#' coefficients (`kappa ~ 1 + x` under bmm's development version, which
#' fits `kappa` and `thetat` as distributional parameters) has its
#' intercept prior on brms's centred intercept `I`, while the intercept
#' reported is `I - sum(m_j * b_j)`, with `m_j` the mean of predictor `j`
#' in the design brms fits. Its prior SD is then
#' `sqrt(sd_I^2 + sum(m_j^2 * sd_j^2))`, with `sd_j` the closed-form SD
#' of `b_j`; a predictor with mean 0 (a sum contrast in a balanced
#' design) adds nothing. The means come from the data: a fit's own,
#' through `brms::standata()`, or the simulated data of
#' `recovery_grid(prior_sd = "analytic")`. A prior table such as
#' `bmm::default_prior()` returns carries no data, so such an intercept
#' is `NA` there and the message says it is centred; a fit whose data
#' cannot be read gives `NA` too, with a warning that says why. The
#' intercept is also `NA` where a slope with a non-zero mean has no
#' closed form, and the message names that slope.
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
  if (!inherits(x, "brmsfit")) {
    return(prior_sd_table_impl(prior_rows(x), group = group))
  }
  # the fit's own design; read only if an intercept is centred
  means <- function() {
    rlang::check_installed("brms", "to read a fit's design.")
    tryCatch(design_means(brms::standata(x)), error = function(e) {
      # nolint next: object_usage_linter. Used by cli's glue interpolation.
      reason <- conditionMessage(e)
      cli::cli_warn(c(
        "The fit's data could not be read for its predictor means.",
        x = "{reason}"
      ))
      NULL
    })
  }
  prior_sd_table_impl(prior_rows(x),
    group = group, means = means,
    no_means = "could not be read from the fit's data"
  )
}

#' `prior_sd_table()` on prior rows, with the predictor means if known
#'
#' `means` is `NULL`, a named list as [design_means()] returns, or a
#' function of no arguments returning one, called only when a parameter
#' has a centred intercept. `no_means` ends the message's sentence on a
#' centred intercept whose parameter has no means: why it has none.
#' @noRd
prior_sd_table_impl <- function(priors, group = "id", means = NULL,
                                no_means = no_data_reason()) {
  priors$prior <- inherited_priors(priors)
  priors[c("lb", "ub")] <- inherited_bounds(priors)
  sd_rows <- priors$class == "sd" & priors$group == group
  population <- population_prior_rows(priors)
  if (is.function(means)) {
    means <- if (any(population$centred)) means()
  }
  pieces <- list(
    coefficient_prior_rows(population, means),
    sd_prior_rows(priors[sd_rows, , drop = FALSE])
  )
  out <- dplyr::bind_rows(empty_prior_sd_table(), pieces)
  missing <- is.na(out$prior_sd_link)
  if (any(missing)) {
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    what <- paste0(out$term[missing], " (", out$level[missing], ")")
    own <- !is.na(coefficient_sd(population))
    unknown <- population$centred & own &
      !population$par %in% names(means)
    centred <- population$par[unknown]
    centred <- intersect(centred, out$term[missing])
    blocked <- centred_blocking_slopes(population, means)
    blocked <- blocked[names(blocked) %in% out$term[missing]]
    cli::cli_inform(c(
      "No closed-form prior SD for {.val {what}}.",
      i = if (length(centred) > 0L) {
        paste0(
          "The intercept{?s} of {.val {centred}} {?is/are} brms's centred ",
          "intercept: {?its/their} SD needs the means of the predictors, ",
          "which ", no_means, "."
        )
      },
      blocked_slope_bullets(blocked),
      i = "Contraction is {.val {NA}} where the prior SD is; a prior-only \\
           fit, or {.fn prior_check}, covers those terms."
    ))
  }
  out
}

#' Why a prior table has no predictor means, for the message
#' @noRd
no_data_reason <- function() "a prior table without data does not have"

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
coefficient_prior_rows <- function(rows, means = NULL) {
  all_rows <- rows
  rows <- rows[nzchar(rows$coef) & !is_constant_prior(rows$prior), ,
    drop = FALSE
  ]
  if (nrow(rows) == 0L) {
    return(NULL)
  }
  sd <- coefficient_sd(rows)
  centred <- which(rows[["centred"]])
  sd[centred] <- vapply(centred, function(i) {
    centred_intercept_sd(sd[[i]], all_rows, rows$par[[i]], means)
  }, double(1))
  tibble::tibble(
    term = coefficient_terms(rows$par, rows$coef),
    level = coefficient_kinds(rows$par, rows$coef),
    prior = rows$prior,
    prior_sd_link = unname(sd)
  )
}

#' The population coefficients' prior rows: class `b`, and the intercepts
#' brms puts on class `Intercept`
#'
#' bmm's `develop` (c13dd2c5, 2026-11-09) fits `mixture2p`'s `kappa` and
#' `thetat` as distributional parameters, and brms then gives an
#' intercept its own class, `Intercept`, with no `coef` (#52); bmm 1.4.1
#' fitted them as non-linear parameters, whose intercept is class `b`,
#' coef `Intercept`. Such a row becomes coef `Intercept` of class `b`
#' here. brms's `Intercept` class is the centred intercept, which is the
#' intercept reported only when the parameter has no other population
#' coefficient; with one, the row is kept, so that the other
#' coefficients keep their `<par>_<coef>` names, but marked `centred`:
#' [centred_intercept_sd()] gives it the prior SD of the reported
#' intercept where the predictor means are known. A constant prior (a
#' fixed parameter) is dropped by [coefficient_prior_rows()] as before.
#' @noRd
population_prior_rows <- function(priors) {
  b <- priors[priors$class == "b", , drop = FALSE]
  b$centred <- rep(FALSE, nrow(b))
  intercept <- priors[
    priors$class == "Intercept" & !nzchar(priors$coef) & nzchar(priors$par), ,
    drop = FALSE
  ]
  if (nrow(intercept) == 0L) {
    return(b)
  }
  sloped <- unique(b$par[nzchar(b$coef) & b$coef != "Intercept"])
  intercept$class <- "b"
  intercept$coef <- "Intercept"
  intercept$centred <- intercept$par %in% sloped
  rbind(b, intercept)
}

#' The closed-form SD of each class `b` row, `NA` under bounds
#' @noRd
coefficient_sd <- function(rows) {
  sd <- vapply(rows$prior, closed_form_sd, double(1), half = FALSE)
  sd[nzchar(rows$lb) | nzchar(rows$ub)] <- NA_real_
  unname(sd)
}

#' The prior SD of the reported intercept behind a centred one
#'
#' brms fits a distributional parameter with an intercept on centred
#' predictors and reports `b_Intercept = Intercept - sum(m_j * b_j)`,
#' with `m_j` the mean of predictor `j` in the design. Under independent
#' priors its SD is `sqrt(sd_I^2 + sum(m_j^2 * sd_j^2))`. A fixed slope
#' (`constant()`) adds nothing, nor does a predictor with mean 0; any
#' other slope without a closed form, or a parameter whose means are not
#' known, gives `NA`.
#' @noRd
centred_intercept_sd <- function(sd_intercept, rows, par, means) {
  m <- means[[par]]
  if (is.null(m) || is.na(sd_intercept)) {
    return(NA_real_)
  }
  m <- m[m != 0]
  sqrt(sd_intercept^2 + sum(m^2 * slope_sd(rows, par, names(m))^2))
}

#' The closed-form prior SD of each named slope of a parameter: 0 for a
#' fixed one, `NA` for one without a closed form or without a row
#' @noRd
slope_sd <- function(rows, par, coef) {
  slopes <- rows[rows$par == par & rows$coef %in% coef, , drop = FALSE]
  sd <- coefficient_sd(slopes)
  sd[is_constant_prior(slopes$prior)] <- 0
  sd[match(coef, slopes$coef)]
}

#' For each centred intercept with a closed-form prior of its own and
#' known means, the slopes with a non-zero mean and no closed-form SD
#' that leave it `NA`
#' @noRd
centred_blocking_slopes <- function(population, means) {
  own <- !is.na(coefficient_sd(population))
  pars <- population$par[population$centred & own]
  pars <- pars[pars %in% names(means)]
  out <- lapply(stats::setNames(nm = pars), function(par) {
    m <- means[[par]]
    m <- m[m != 0]
    names(m)[is.na(slope_sd(population, par, names(m)))]
  })
  out[lengths(out) > 0L]
}

#' One message bullet per centred intercept that its slopes leave `NA`,
#' interpolated in the caller's frame, which holds `blocked`
#' @noRd
blocked_slope_bullets <- function(blocked) {
  bullets <- vapply(seq_along(blocked), function(k) {
    sprintf(
      paste0(
        "The intercept of {.val {names(blocked)[[%d]]}} is brms's centred ",
        "intercept: its SD needs the prior SD of {?its slope/its slopes} ",
        "{.val {blocked[[%d]]}}, which {?has/have} no closed form."
      ),
      k, k
    )
  }, character(1))
  stats::setNames(bullets, rep("i", length(bullets)))
}

#' Each row's bounds after brms's inheritance
#'
#' brms takes `lb` and `ub` on a class row only (`coef` empty; measured
#' on brms 2.23 under bmm 1.4.2.9000: "Prior argument 'coef' may not be
#' specified when using boundaries") and bounds every coefficient of
#' that class and parameter by them, while [brms::prior_summary()]
#' leaves the coefficient rows' bounds empty: Stan declares
#' `vector<lower=0> b_kappa`. A coefficient row takes the bounds of the
#' class row of its parameter and group.
#' @noRd
inherited_bounds <- function(priors) {
  class_row <- match(
    paste(priors$class, priors$group, priors$par, sep = "\r"),
    ifelse(
      nzchar(priors$coef), NA_character_,
      paste(priors$class, priors$group, priors$par, sep = "\r")
    )
  )
  own <- nzchar(priors$lb) | nzchar(priors$ub) | !nzchar(priors$coef)
  take <- !own & !is.na(class_row)
  lb <- priors$lb
  ub <- priors$ub
  lb[take] <- priors$lb[class_row[take]]
  ub[take] <- priors$ub[class_row[take]]
  list(lb = lb, ub = ub)
}

#' The column means of each parameter's population design, without the
#' intercept, from a Stan data list (`bmm::standata()`,
#' `brms::standata()`): `X_<par>` is the design brms centres.
#' @noRd
design_means <- function(standata) {
  design <- grep("^X_", names(standata), value = TRUE)
  out <- lapply(standata[design], function(x) {
    x <- as.matrix(x)
    m <- colMeans(x)
    # a mean that is 0 but for rounding (poly(), a centred covariate) is
    # 0, so that it adds nothing even beside a flat slope
    scale <- apply(abs(x), 2L, max)
    m[abs(m) <= sqrt(.Machine$double.eps) * scale] <- 0
    m[names(m) != "Intercept"]
  })
  names(out) <- sub("^X_", "", design)
  out
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
