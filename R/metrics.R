# Recovery metrics (local/ARCHITECTURE.md decision 18).
#
# Every function here takes plain numeric vectors and returns a scalar
# (or a small named list). None of them touches a fit object, which is
# what lets the whole score layer be tested without brms.
#
# They are internal: summary.bmmtools_recovery() is the supported entry
# point. The exception is recovery_ccc(), exported for vectors that are not
# a recovery object (Amendment 2026-09-14); the rest is a 0.2 question.

#' Drop incomplete pairs
#'
#' Pairwise deletion, applied before any metric counts its sample size,
#' so that `n` always reports the pairs actually used.
#'
#' @return A list with the two filtered vectors and `n`.
#' @noRd
complete_pairs <- function(x, y, arg_x = "estimate", arg_y = "truth") {
  if (length(x) != length(y)) {
    cli::cli_abort(
      "{.arg {arg_x}} and {.arg {arg_y}} must be the same length, \\
       not {length(x)} and {length(y)}."
    )
  }
  keep <- !is.na(x) & !is.na(y)
  list(x = x[keep], y = y[keep], n = sum(keep))
}

#' Is a vector too degenerate for a correlation?
#'
#' Fewer than three pairs, or no spread on either side. Both cases return
#' `NA` rather than `0`: the seeds returned `0` here, which reads as "no
#' recovery" when the honest answer is "not estimable from this design".
#'
#' @noRd
correlation_undefined <- function(x, y) {
  length(x) < 3L ||
    !stats::var(x) > 0 ||
    !stats::var(y) > 0
}

#' Mean signed error
#' @noRd
metric_bias <- function(estimate, truth) {
  p <- complete_pairs(estimate, truth)
  if (p$n == 0L) {
    return(NA_real_)
  }
  mean(p$x - p$y)
}

#' Root mean squared error
#' @noRd
metric_rmse <- function(estimate, truth) {
  p <- complete_pairs(estimate, truth)
  if (p$n == 0L) {
    return(NA_real_)
  }
  sqrt(mean((p$x - p$y)^2))
}

#' Share of credible intervals containing the generating value
#'
#' The interval is closed: a bound exactly equal to truth counts as
#' covered.
#'
#' @noRd
metric_coverage <- function(truth, ci_low, ci_high) {
  if (length(truth) != length(ci_low) || length(truth) != length(ci_high)) {
    cli::cli_abort(
      "{.arg truth}, {.arg ci_low} and {.arg ci_high} must be the same \\
       length, not {length(truth)}, {length(ci_low)} and {length(ci_high)}."
    )
  }
  keep <- !is.na(truth) & !is.na(ci_low) & !is.na(ci_high)
  if (!any(keep)) {
    return(NA_real_)
  }
  mean(truth[keep] >= ci_low[keep] & truth[keep] <= ci_high[keep])
}

#' Mean credible-interval width
#' @noRd
metric_ci_width <- function(ci_low, ci_high) {
  p <- complete_pairs(ci_high, ci_low, "ci_high", "ci_low")
  if (p$n == 0L) {
    return(NA_real_)
  }
  mean(p$x - p$y)
}

#' Pearson correlation with a Fisher-z interval
#'
#' The interval is computed in base R rather than by calling
#' `stats::cor.test()`, so nothing at runtime depends on the shape of its
#' return object; a test asserts the two agree to 1e-8 across several
#' sample sizes and confidence levels.
#'
#' @return A named list: `r`, `r_low`, `r_high`, `n`.
#' @noRd
metric_r <- function(estimate, truth, ci_level = 0.95) {
  p <- complete_pairs(estimate, truth)
  out <- list(r = NA_real_, r_low = NA_real_, r_high = NA_real_, n = p$n)
  if (correlation_undefined(p$x, p$y)) {
    return(out)
  }

  out$r <- stats::cor(p$x, p$y, method = "pearson")

  # Fisher's z: the interval needs n > 3, so a 3-pair correlation gets a
  # point estimate and no interval rather than a division by zero.
  if (p$n > 3L && abs(out$r) < 1) {
    z <- atanh(out$r)
    se <- 1 / sqrt(p$n - 3)
    crit <- stats::qnorm(1 - (1 - ci_level) / 2)
    out$r_low <- tanh(z - crit * se)
    out$r_high <- tanh(z + crit * se)
  }
  out
}

#' Spearman rank correlation
#'
#' Reported alongside `metric_r()` because a model can order subjects
#' perfectly while a monotone link makes the Pearson correlation less
#' than one.
#'
#' @noRd
metric_rank_r <- function(estimate, truth) {
  p <- complete_pairs(estimate, truth)
  if (correlation_undefined(p$x, p$y)) {
    return(NA_real_)
  }
  stats::cor(p$x, p$y, method = "spearman")
}

#' Asymptotic variance of Lin's Z-transformed concordance
#'
#' The variance of `atanh(ccc)` for random bivariate normal pairs: Lin
#' (1989, p. 259, Eq. 2) as corrected in Lin (2000, Biometrics 56, p. 325,
#' item 3), where the constant of the second term is 2, not 4, and the 2
#' of the third term belongs in the denominator. `u` is the location
#' shift relative to the geometric mean of the standard deviations, and
#' `n - 2` replaces `n` as Lin recommends for less bias (p. 268).
#'
#' `NA` where Z has no finite variance or the formula divides by zero:
#' three or fewer pairs, a coefficient of +/-1, or a zero correlation. A
#' result that is not positive is also `NA`; it cannot arise from one data
#' set, whose `ccc`, `r` and `u` are tied together, only from inconsistent
#' input.
#'
#' @noRd
ccc_z_variance <- function(ccc, r, u, n) {
  if (anyNA(c(ccc, r, u, n)) || n <= 3 || abs(ccc) >= 1 || r == 0) {
    return(NA_real_)
  }
  one_minus <- 1 - ccc^2
  precision_term <- (1 - r^2) * ccc^2 / (one_minus * r^2)
  location_term <- 2 * ccc^3 * (1 - ccc) * u^2 / (r * one_minus^2)
  quartic_term <- ccc^4 * u^4 / (2 * r^2 * one_minus^2)
  out <- (precision_term + location_term - quartic_term) / (n - 2)
  if (!is.finite(out) || out <= 0) NA_real_ else out
}

#' Lin's concordance correlation coefficient
#'
#' Uses the population (divide by `n`) moments, which is Lin's original
#' definition (Lin, 1989, p. 258). `truth` is the reference throughout.
#'
#' Besides the coefficient, returns Lin's decomposition `ccc = r *
#' accuracy`, with `accuracy = 2 / (v + 1/v + u^2)`, the scale shift
#' `v = sd(estimate) / sd(truth)` and the location shift
#' `u = (mean(estimate) - mean(truth)) / sqrt(sd(estimate) sd(truth))`;
#' the calibration slope, r divided by `v`, which is the slope of truth
#' regressed on the estimate; and a 95% interval from Lin's Z transformation.
#'
#' A calibrated posterior mean has `v = r` and `calibration_slope = 1`,
#' not `v = 1`: shrinkage is not scale bias (local/dev/eval-ccc-2026-09-14.md).
#'
#' @return A named list: `ccc`, `ccc_low`, `ccc_high`, `accuracy`,
#'   `scale_shift`, `location_shift`, `calibration_slope`, `r`, `var_z`,
#'   `truth_sd`, `n`.
#' @noRd
metric_ccc <- function(estimate, truth) {
  p <- complete_pairs(estimate, truth)
  y <- p$x # estimate
  x <- p$y # truth
  n <- p$n

  out <- list(
    ccc = NA_real_, ccc_low = NA_real_, ccc_high = NA_real_,
    accuracy = NA_real_, scale_shift = NA_real_, location_shift = NA_real_,
    calibration_slope = NA_real_, r = NA_real_, var_z = NA_real_,
    truth_sd = NA_real_, n = n
  )
  if (n == 0L) {
    return(out)
  }
  mean_x <- mean(x)
  # population moments: Lin's original definition
  var_x <- sum((x - mean_x)^2) / n
  out$truth_sd <- sqrt(var_x)
  if (correlation_undefined(y, x)) {
    return(out)
  }

  mean_y <- mean(y)
  var_y <- sum((y - mean_y)^2) / n
  cov_xy <- sum((x - mean_x) * (y - mean_y)) / n
  sd_x <- sqrt(var_x)
  sd_y <- sqrt(var_y)

  out$ccc <- 2 * cov_xy / (var_x + var_y + (mean_x - mean_y)^2)
  out$r <- cov_xy / (sd_x * sd_y)
  out$scale_shift <- sd_y / sd_x
  out$location_shift <- (mean_y - mean_x) / sqrt(sd_x * sd_y)
  shifts <- out$scale_shift + 1 / out$scale_shift + out$location_shift^2
  out$accuracy <- 2 / shifts
  out$calibration_slope <- cov_xy / var_y

  out$var_z <- ccc_z_variance(out$ccc, out$r, out$location_shift, n)
  if (!is.na(out$var_z)) {
    crit <- stats::qnorm(0.975)
    z <- atanh(out$ccc)
    out$ccc_low <- tanh(z - crit * sqrt(out$var_z))
    out$ccc_high <- tanh(z + crit * sqrt(out$var_z))
  }
  out
}

#' Concordance between recovered and generating values
#'
#' Lin's concordance correlation coefficient with its 95% interval and
#' decomposition, for any pair of vectors: subject estimates against
#' their generating values, or recovered correlations and effects
#' against the ones that were simulated. `summary()` of a recovery
#' object reports the same columns per parameter.
#'
#' The coefficient is `1 - E[(estimate - truth)^2] / (sd_e^2 + sd_t^2 +
#' (mean_e - mean_t)^2)`, a rescaled mean squared deviation, and factors
#' into Pearson's `r` (precision) times `ccc_accuracy` (Lin, 1989, p. 258).
#' Accuracy depends on the scale shift `v = sd_e / sd_t` and the location
#' shift `u = (mean_e - mean_t) / sqrt(sd_e sd_t)`; all moments divide by
#' `n`.
#'
#' **Reading the scale shift.** Estimates from a hierarchical model are
#' shrunk towards the population mean, and a well-calibrated posterior
#' mean has `ccc_scale_shift` equal to `r`, not 1. `calibration_slope`
#' (r divided by `v`, the slope of truth regressed on the estimate) is 1 in that
#' case, above 1 when estimates are shrunk too much and below 1 when they
#' are shrunk too little. Because accuracy treats `v` and `1/v` alike,
#' `ccc` rewards too little shrinkage: read it together with the RMSE and
#' the calibration slope, and prefer `r` when only the rank order of
#' subjects matters.
#'
#' **The interval** uses Lin's Z transformation, `atanh(ccc)`, with the
#' asymptotic variance for random bivariate normal pairs, and is `NA`
#' with three or fewer pairs, at `ccc = +/-1` and when the Pearson
#' correlation is 0. When the
#' generating values are fixed by design rather than drawn, it is an
#' approximation. It is always a 95% interval.
#'
#' @param estimate,truth Numeric vectors of the same length. Incomplete
#'   pairs are dropped.
#'
#' @return A one-row tibble with `ccc`, `ccc_low`, `ccc_high`,
#'   `ccc_accuracy`, `ccc_scale_shift`, `ccc_location_shift`,
#'   `calibration_slope` and `n`, the number of complete pairs. The
#'   coefficient and its components are `NA` with fewer than three pairs
#'   or when either vector has no spread.
#'
#' @references Lin, L. I.-K. (1989). A concordance correlation
#'   coefficient to evaluate reproducibility. *Biometrics, 45*(1),
#'   255--268. \doi{10.2307/2532051}
#'
#'   Lin, L. I.-K. (2000). A note on the concordance correlation
#'   coefficient. *Biometrics, 56*(1), 324--325.
#'   \url{https://www.jstor.org/stable/2677159}
#'
#' @examples
#' truth <- c(0.2, 0.5, 0.9, 1.4, 1.8, 2.3)
#' # shrunk towards the mean and shifted up
#' estimate <- 0.7 * truth + 0.5
#' recovery_ccc(estimate, truth)
#'
#' @export
recovery_ccc <- function(estimate, truth) {
  if (!is.numeric(estimate)) {
    cli::cli_abort(
      "{.arg estimate} must be numeric, not {.obj_type_friendly {estimate}}."
    )
  }
  if (!is.numeric(truth)) {
    cli::cli_abort(
      "{.arg truth} must be numeric, not {.obj_type_friendly {truth}}."
    )
  }
  if (length(estimate) != length(truth)) {
    cli::cli_abort(
      "{.arg estimate} and {.arg truth} must be the same length, \\
       not {length(estimate)} and {length(truth)}."
    )
  }
  out <- metric_ccc(estimate, truth)
  tibble::tibble(
    ccc = out$ccc,
    ccc_low = out$ccc_low,
    ccc_high = out$ccc_high,
    ccc_accuracy = out$accuracy,
    ccc_scale_shift = out$scale_shift,
    ccc_location_shift = out$location_shift,
    calibration_slope = out$calibration_slope,
    n = as.integer(out$n)
  )
}
