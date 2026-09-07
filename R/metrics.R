# Recovery metrics (ARCHITECTURE.md decision 18).
#
# Every function here takes plain numeric vectors and returns a scalar
# (or a small named list). None of them touches a fit object, which is
# what lets the whole score layer be tested without brms.
#
# They are internal: summary.bmmtools_recovery() is the supported entry
# point. Exporting them is a 0.2 question.

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

#' Lin's concordance correlation coefficient
#'
#' Uses the population (divide by `n`) moments, which is Lin's original
#' definition. Reported next to `metric_r()` because a set of estimates
#' can correlate perfectly with truth while being shifted or stretched,
#' which is not recovery; `scale_bias` and `location_bias` say which of
#' the two happened.
#'
#' `truth` is the reference: `scale_bias` is the estimates' spread
#' relative to truth's, and `location_bias` is the estimates' offset.
#'
#' @return A named list: `ccc`, `scale_bias`, `location_bias`.
#' @noRd
metric_ccc <- function(estimate, truth) {
  p <- complete_pairs(estimate, truth)
  out <- list(ccc = NA_real_, scale_bias = NA_real_, location_bias = NA_real_)
  if (correlation_undefined(p$x, p$y)) {
    return(out)
  }

  y <- p$x # estimate
  x <- p$y # truth
  n <- p$n

  mean_x <- mean(x)
  mean_y <- mean(y)
  # population moments: Lin's original definition
  var_x <- sum((x - mean_x)^2) / n
  var_y <- sum((y - mean_y)^2) / n
  cov_xy <- sum((x - mean_x) * (y - mean_y)) / n

  out$ccc <- 2 * cov_xy / (var_x + var_y + (mean_x - mean_y)^2)
  out$scale_bias <- sqrt(var_y) / sqrt(var_x)
  out$location_bias <- (mean_y - mean_x) / sqrt(sqrt(var_x) * sqrt(var_y))
  out
}
