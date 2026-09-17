# Tests for the metric family, written against
# local/dev/spec-milestone-1-score-layer.md section 3.
#
# Every expected value here is either hand-computed (with the arithmetic
# shown in a comment) or produced by an oracle in the test itself
# (stats::cor.test, stats::cor). None is typed from memory.

# bias and rmse ----------------------------------------------------------

test_that("metric_bias is the mean signed difference", {
  # differences are 1, 1, 0 and 3; their mean is 5 / 4 = 1.25
  estimate <- c(2, 3, 4, 7)
  truth <- c(1, 2, 4, 4)
  expect_equal(metric_bias(estimate, truth), 1.25)
})

test_that("metric_bias is signed, not absolute", {
  # a model that under- and over-estimates symmetrically has zero bias
  expect_equal(metric_bias(c(0, 2), c(1, 1)), 0)
  expect_equal(metric_bias(c(1, 1), c(2, 4)), -2)
})

test_that("metric_rmse is the root mean squared difference", {
  # squared differences c(1, 1, 0, 9); mean = 11 / 4 = 2.75
  estimate <- c(2, 3, 4, 7)
  truth <- c(1, 2, 4, 4)
  expect_equal(metric_rmse(estimate, truth), sqrt(2.75))
})

test_that("metric_rmse is zero only for exact recovery", {
  expect_equal(metric_rmse(c(1, 2, 3), c(1, 2, 3)), 0)
  expect_gt(metric_rmse(c(1, 2, 3), c(1, 2, 3.001)), 0)
})

test_that("rmse punishes a single large miss more than bias does", {
  # the property that makes both worth reporting: offsetting errors
  # cancel in bias and accumulate in rmse
  estimate <- c(-10, 10)
  truth <- c(0, 0)
  expect_equal(metric_bias(estimate, truth), 0)
  expect_equal(metric_rmse(estimate, truth), 10)
})

test_that("metric_mae is the mean absolute difference", {
  # absolute differences c(1, 1, 0, 3); mean = 5 / 4 = 1.25
  estimate <- c(2, 3, 4, 7)
  truth <- c(1, 2, 4, 4)
  expect_equal(metric_mae(estimate, truth), 1.25)
})

test_that("mae does not cancel offsetting errors, unlike bias", {
  estimate <- c(-10, 10)
  truth <- c(0, 0)
  expect_equal(metric_bias(estimate, truth), 0)
  expect_equal(metric_mae(estimate, truth), 10)
})

# coverage and interval width --------------------------------------------

test_that("metric_coverage is the share of intervals containing truth", {
  truth <- c(1, 2, 3, 4)
  ci_low <- c(0, 5, 2, 3)
  ci_high <- c(2, 6, 4, 5)
  # covered: TRUE, FALSE, TRUE, TRUE -> 3 / 4
  expect_equal(metric_coverage(truth, ci_low, ci_high), 0.75)
})

test_that("coverage counts an interval whose bound equals truth", {
  # A closed interval. Stated as a test so a later reader knows this was
  # a choice and not an accident of the comparison operator.
  expect_equal(metric_coverage(1, 1, 2), 1)
  expect_equal(metric_coverage(2, 1, 2), 1)
  expect_equal(metric_coverage(2.0001, 1, 2), 0)
})

test_that("metric_ci_width is the mean interval width", {
  # widths c(2, 1, 6); mean = 3
  expect_equal(metric_ci_width(c(0, 1, 2), c(2, 2, 8)), 3)
})

test_that("coverage and width are NA when no complete row is left", {
  # Every row incomplete: NA, never 0 and never NaN from mean(logical(0)).
  # These two are the only metrics that survive a zero-variance truth, so
  # the all-missing case is the one place they can go wrong quietly.
  expect_identical(
    metric_coverage(c(NA_real_, 2), c(0, NA_real_), c(1, 3)),
    NA_real_
  )
  expect_identical(
    metric_ci_width(c(NA_real_, 1), c(2, NA_real_)),
    NA_real_
  )
})

# correlation ------------------------------------------------------------

test_that("metric_r matches stats::cor.test on point and interval", {
  withr::local_seed(1)
  truth <- stats::rnorm(30)
  estimate <- truth * 0.8 + stats::rnorm(30, sd = 0.5)

  oracle <- stats::cor.test(estimate, truth, method = "pearson")
  out <- metric_r(estimate, truth, ci_level = 0.95)

  expect_equal(out$r, unname(oracle$estimate))
  expect_equal(out$r_low, oracle$conf.int[1], tolerance = 1e-8)
  expect_equal(out$r_high, oracle$conf.int[2], tolerance = 1e-8)
  expect_equal(out$n, 30L)
})

test_that("metric_r matches cor.test across several n and ci_level", {
  withr::local_seed(2)
  for (n in c(5L, 12L, 60L)) {
    for (level in c(0.95, 0.89)) {
      truth <- stats::rnorm(n)
      estimate <- truth + stats::rnorm(n, sd = 0.7)
      oracle <- stats::cor.test(estimate, truth, conf.level = level)
      out <- metric_r(estimate, truth, ci_level = level)
      info <- paste("n =", n, "level =", level)
      expect_equal(out$r, unname(oracle$estimate), info = info)
      tol <- 1e-8
      expect_equal(out$r_low, oracle$conf.int[1], tolerance = tol, info = info)
      expect_equal(out$r_high, oracle$conf.int[2], tolerance = tol, info = info)
    }
  }
})

test_that("metric_r returns 1 for perfect recovery and -1 for reversal", {
  x <- c(1, 2, 3, 4, 5)
  expect_equal(metric_r(x, x)$r, 1)
  expect_equal(metric_r(-x, x)$r, -1)
})

test_that("metric_rank_r matches stats::cor with method spearman", {
  withr::local_seed(3)
  truth <- stats::rnorm(25)
  estimate <- exp(truth) + stats::rnorm(25, sd = 0.1)
  expect_equal(
    metric_rank_r(estimate, truth),
    stats::cor(estimate, truth, method = "spearman")
  )
})

test_that("rank_r is 1 under a monotone transform where r is not", {
  # the property that makes the rank version worth reporting alongside r
  truth <- c(1, 2, 3, 4, 5)
  estimate <- exp(truth)
  expect_equal(metric_rank_r(estimate, truth), 1)
  expect_lt(metric_r(estimate, truth)$r, 1)
})

# concordance ------------------------------------------------------------

test_that("metric_ccc matches a hand-computed value", {
  # truth x = 1:4, estimate y = x + 1
  # mean_x = 2.5, mean_y = 3.5
  # population variances: s_x^2 = s_y^2 = mean(2.25, .25, .25, 2.25) = 1.25
  # covariance s_xy = 1.25
  # so ccc is twice 1.25, over 1.25 plus 1.25 plus one squared:
  # 2.5 over 3.5, which is 5/7
  truth <- c(1, 2, 3, 4)
  estimate <- c(2, 3, 4, 5)
  out <- metric_ccc(estimate, truth)
  expect_equal(out$ccc, 5 / 7)
  # scale_shift is the estimates' SD over truth's, here equal
  expect_equal(out$scale_shift, 1)
  # location_shift is Lin's u = (mean_y - mean_x) / sqrt(sd_x * sd_y),
  # so 1 / sqrt(1.25) = 0.8944272, NOT 1 / 1.25. The denominator is the
  # geometric mean of the two standard deviations, not of the variances;
  # the decomposition test below is what pins this down.
  expect_equal(out$location_shift, 1 / sqrt(1.25))
})

test_that("ccc decomposes as r times Lin's bias-correction factor", {
  # ccc = r * C_b with C_b = 2 / (v + 1/v + u^2), v = scale_shift,
  # u = location_shift. This identity holds only for Lin's definitions of
  # u and v, so it is what fixes them: no reference implementation is
  # installed, and a hand-typed constant would only restate the code.
  withr::local_seed(11)
  for (i in 1:5) {
    truth <- stats::rnorm(20)
    estimate <- truth * stats::runif(1, 0.5, 2) +
      stats::rnorm(1) + stats::rnorm(20, sd = 0.3)

    out <- metric_ccc(estimate, truth)
    r <- metric_r(estimate, truth)$r
    v <- out$scale_shift
    u <- out$location_shift
    c_b <- 2 / (v + 1 / v + u^2)

    expect_equal(out$ccc, r * c_b, tolerance = 1e-10)
    expect_equal(out$accuracy, c_b, tolerance = 1e-10)
    expect_equal(out$ccc, stats::cor(estimate, truth) * out$accuracy,
      tolerance = 1e-10
    )
  }
})

test_that("ccc is 1 for exact agreement", {
  x <- c(2, 4, 6, 9)
  expect_equal(metric_ccc(x, x)$ccc, 1)
})

test_that("ccc penalises a location shift that r does not see", {
  # This is the reason CCC is reported next to r: a model whose estimates
  # are perfectly correlated with truth but shifted is not recovering it.
  truth <- c(1, 2, 3, 4)
  shifted <- truth + 5
  expect_equal(metric_r(shifted, truth)$r, 1)
  expect_lt(metric_ccc(shifted, truth)$ccc, 1)
})

test_that("ccc penalises a scale error that r does not see", {
  truth <- c(1, 2, 3, 4)
  stretched <- truth * 3
  expect_equal(metric_r(stretched, truth)$r, 1)
  expect_lt(metric_ccc(stretched, truth)$ccc, 1)
  expect_equal(metric_ccc(stretched, truth)$scale_shift, 3)
})

test_that("ccc is negative for a reversed relationship", {
  truth <- c(1, 2, 3, 4)
  expect_lt(metric_ccc(-truth, truth)$ccc, 0)
})

# concordance: interval and calibration ----------------------------------

# Population (divide-by-n) standardisation, so that two vectors share a
# mean and an SD exactly as metric_ccc() measures them.
pop_standardise <- function(x) {
  centred <- x - mean(x)
  centred / sqrt(mean(centred^2))
}

test_that("with equal moments ccc is r and its Z variance is 1/(n - 2)", {
  # With u = 0 and v = 1, C_b = 1 so ccc = r, and every term of Lin's
  # variance but the first vanishes: (1 - r^2) r^2 / ((1 - r^2) r^2),
  # divided by n - 2. The oracle needs no copy of the formula.
  withr::local_seed(101)
  n <- 30L
  truth <- stats::rnorm(n)
  estimate <- truth + stats::rnorm(n, sd = 0.6)
  truth <- pop_standardise(truth)
  estimate <- pop_standardise(estimate)

  out <- metric_ccc(estimate, truth)
  expect_equal(out$ccc, stats::cor(estimate, truth), tolerance = 1e-12)
  expect_equal(out$var_z, 1 / (n - 2), tolerance = 1e-12)
  crit <- stats::qnorm(0.975)
  expect_equal(out$ccc_low, tanh(atanh(out$ccc) - crit / sqrt(n - 2)),
    tolerance = 1e-12
  )
  expect_equal(out$ccc_high, tanh(atanh(out$ccc) + crit / sqrt(n - 2)),
    tolerance = 1e-12
  )
})

test_that("with no location shift the Z variance keeps only its first term", {
  withr::local_seed(102)
  n <- 40L
  truth <- pop_standardise(stats::rnorm(n))
  estimate <- 0.5 * pop_standardise(truth + stats::rnorm(n, sd = 0.8))

  out <- metric_ccc(estimate, truth)
  expect_equal(out$location_shift, 0, tolerance = 1e-12)
  r <- out$r
  p <- out$ccc
  expect_equal(out$var_z, (1 - r^2) * p^2 / ((1 - p^2) * r^2) / (n - 2),
    tolerance = 1e-10
  )
})

test_that("a location shift up or down gives the same interval", {
  # u enters the variance only as u^2 and u^4.
  withr::local_seed(103)
  truth <- stats::rnorm(25)
  estimate <- truth * 0.8 + stats::rnorm(25, sd = 0.5)
  estimate <- estimate - mean(estimate) + mean(truth)

  up <- metric_ccc(estimate + 0.7, truth)
  down <- metric_ccc(estimate - 0.7, truth)
  expect_equal(up$location_shift, -down$location_shift)
  expect_equal(up$ccc_low, down$ccc_low, tolerance = 1e-12)
  expect_equal(up$ccc_high, down$ccc_high, tolerance = 1e-12)
  expect_lt(up$ccc_low, up$ccc)
  expect_gt(up$ccc_high, up$ccc)
})

test_that("calibration_slope is the slope of truth regressed on estimate", {
  withr::local_seed(104)
  truth <- stats::rnorm(30)
  estimate <- 0.6 * truth + 0.3 + stats::rnorm(30, sd = 0.4)
  out <- metric_ccc(estimate, truth)
  oracle <- unname(stats::coef(stats::lm(truth ~ estimate))[[2L]])
  expect_equal(out$calibration_slope, oracle, tolerance = 1e-10)
  expect_equal(out$calibration_slope, out$r / out$scale_shift,
    tolerance = 1e-10
  )
})

test_that("swapping estimate and truth mirrors the components", {
  withr::local_seed(105)
  truth <- stats::rnorm(30)
  estimate <- 1.4 * truth - 0.2 + stats::rnorm(30, sd = 0.5)
  fwd <- metric_ccc(estimate, truth)
  rev <- metric_ccc(truth, estimate)

  expect_equal(rev$ccc, fwd$ccc)
  expect_equal(rev$ccc_low, fwd$ccc_low)
  expect_equal(rev$ccc_high, fwd$ccc_high)
  expect_equal(rev$accuracy, fwd$accuracy)
  expect_equal(rev$scale_shift, 1 / fwd$scale_shift)
  expect_equal(rev$location_shift, -fwd$location_shift)
  expect_equal(
    rev$calibration_slope,
    unname(stats::coef(stats::lm(estimate ~ truth))[[2L]]),
    tolerance = 1e-10
  )
})

test_that("truth_sd uses population moments and pairs with scale_shift", {
  truth <- c(1, 2, 3, 4, 10)
  estimate <- c(2, 2, 5, 3, 8)
  out <- metric_ccc(estimate, truth)
  expect_equal(out$truth_sd, sqrt(mean((truth - mean(truth))^2)))
  expect_equal(
    out$truth_sd * out$scale_shift,
    sqrt(mean((estimate - mean(estimate))^2))
  )
})

test_that("three pairs give a concordance but no interval", {
  out <- metric_ccc(c(1, 3, 2), c(1, 2, 3))
  expect_false(is.na(out$ccc))
  expect_false(is.na(out$accuracy))
  expect_true(is.na(out$ccc_low))
  expect_true(is.na(out$ccc_high))
  expect_true(is.na(out$var_z))
})

test_that("perfect agreement gives ccc 1 and no interval", {
  x <- c(2, 4, 6, 9, 11)
  out <- metric_ccc(x, x)
  expect_equal(out$ccc, 1)
  expect_true(is.na(out$ccc_low))
  expect_true(is.na(out$ccc_high))
  expect_equal(out$calibration_slope, 1)
})

test_that("zero covariance gives ccc 0, slope 0 and no interval", {
  # deviations of truth: -1.5, -.5, .5, 1.5; estimate: 1, -1, -1, 1;
  # their products sum to zero, so r = 0 and Lin's variance divides by 0
  out <- metric_ccc(c(1, -1, -1, 1), c(1, 2, 3, 4))
  expect_equal(out$ccc, 0)
  expect_equal(out$calibration_slope, 0)
  expect_true(is.na(out$ccc_low))
  expect_true(is.na(out$ccc_high))
})

test_that("the Z variance is NA rather than negative for impossible input", {
  # r = .9, ccc = .8, u = 3 cannot come from one data set (C_b would be
  # far below .8 / .9); the u^4 term then drives the sum negative.
  expect_true(is.na(ccc_z_variance(ccc = 0.8, r = 0.9, u = 3, n = 50)))
  expect_true(is.na(ccc_z_variance(ccc = 0.5, r = 0.6, u = 0, n = 3)))
  expect_true(is.na(ccc_z_variance(ccc = 1, r = 1, u = 0, n = 10)))
  expect_true(is.na(ccc_z_variance(ccc = 0, r = 0, u = 0.2, n = 10)))
})

test_that("no complete pair gives NA everywhere, truth_sd included", {
  out <- metric_ccc(c(NA, 1), c(2, NA))
  expect_identical(out$n, 0L)
  expect_true(all(is.na(unlist(out[setdiff(names(out), "n")]))))
})

test_that("the interval guards do not warn", {
  expect_silent(metric_ccc(c(1, 3, 2), c(1, 2, 3)))
  expect_silent(metric_ccc(c(1, -1, -1, 1), c(1, 2, 3, 4)))
  expect_silent(metric_ccc(1:5, 1:5))
  expect_silent(ccc_z_variance(ccc = 0.8, r = 0.9, u = 3, n = 50))
})

test_that("Lin's interval covers the population concordance", {
  # truth ~ N(0, 1), estimate = a + b * truth + N(0, s^2), so
  # ccc = 2b / (1 + b^2 + s^2 + a^2). 400 data sets of 50 pairs; the
  # band is the binomial 99% range around .95.
  withr::local_seed(106)
  a <- 0.4
  b <- 0.7
  s <- 0.5
  target <- 2 * b / (1 + b^2 + s^2 + a^2)
  hits <- vapply(seq_len(400), function(i) {
    truth <- stats::rnorm(50)
    estimate <- a + b * truth + stats::rnorm(50, sd = s)
    out <- metric_ccc(estimate, truth)
    out$ccc_low <= target && target <= out$ccc_high
  }, logical(1))
  expect_gte(mean(hits), 0.917)
  expect_lte(mean(hits), 0.983)
})

# guards -----------------------------------------------------------------

test_that("correlation metrics are NA with fewer than three pairs", {
  # NA, never 0: 0 reads as "no recovery" when the honest answer is
  # "not estimable from this design" (spec section 3, Guards).
  estimate <- c(1, 2)
  truth <- c(1, 3)
  expect_true(is.na(metric_r(estimate, truth)$r))
  expect_true(is.na(metric_rank_r(estimate, truth)))
  expect_true(is.na(metric_ccc(estimate, truth)$ccc))
})

test_that("bias, rmse, coverage and width survive fewer than three pairs", {
  estimate <- c(1, 2)
  truth <- c(1, 3)
  expect_equal(metric_bias(estimate, truth), -0.5)
  expect_false(is.na(metric_rmse(estimate, truth)))
  expect_false(is.na(metric_coverage(truth, c(0, 0), c(5, 5))))
  expect_false(is.na(metric_ci_width(c(0, 0), c(5, 5))))
})

test_that("zero spread in truth gives NA correlations, not zero", {
  # The seeds hit this on delta thresholds fixed by design.
  truth <- c(2, 2, 2, 2)
  estimate <- c(1, 2, 3, 4)
  expect_true(is.na(metric_r(estimate, truth)$r))
  expect_true(is.na(metric_rank_r(estimate, truth)))
  expect_true(is.na(metric_ccc(estimate, truth)$ccc))
  # and the ones that are still defined are still returned
  expect_equal(metric_bias(estimate, truth), 0.5)
})

test_that("zero spread in the estimate gives NA correlations", {
  truth <- c(1, 2, 3, 4)
  estimate <- c(2, 2, 2, 2)
  expect_true(is.na(metric_r(estimate, truth)$r))
  expect_true(is.na(metric_rank_r(estimate, truth)))
  expect_true(is.na(metric_ccc(estimate, truth)$ccc))
})

test_that("the guards do not warn", {
  # summary() is where a user learns why a cell is NA; a warning per
  # metric per parameter would drown the console on a real grid.
  truth <- c(2, 2, 2)
  estimate <- c(1, 2, 3)
  expect_silent(metric_r(estimate, truth))
  expect_silent(metric_rank_r(estimate, truth))
  expect_silent(metric_ccc(estimate, truth))
})

test_that("incomplete pairs are dropped pairwise and n reflects it", {
  estimate <- c(1, 2, NA, 4, 5)
  truth <- c(1, 2, 3, NA, 5)
  # complete pairs: positions 1, 2, 5 -> n = 3
  expect_equal(metric_r(estimate, truth)$n, 3L)
  expect_equal(metric_bias(estimate, truth), 0)
  expect_equal(metric_rmse(estimate, truth), 0)
})

test_that("dropping to fewer than three complete pairs triggers the guard", {
  estimate <- c(1, 2, NA, NA)
  truth <- c(1, 2, 3, 4)
  expect_true(is.na(metric_r(estimate, truth)$r))
  expect_equal(metric_r(estimate, truth)$n, 2L)
})

test_that("coverage ignores rows with a missing bound", {
  truth <- c(1, 2, 3)
  ci_low <- c(0, NA, 2)
  ci_high <- c(2, 3, 4)
  # usable rows 1 and 3, both covered
  expect_equal(metric_coverage(truth, ci_low, ci_high), 1)
})

test_that("all-missing input returns NA rather than NaN", {
  estimate <- c(NA_real_, NA_real_)
  truth <- c(NA_real_, NA_real_)
  expect_true(is.na(metric_bias(estimate, truth)))
  expect_true(is.na(metric_rmse(estimate, truth)))
  expect_true(is.na(metric_mae(estimate, truth)))
  expect_true(is.na(metric_r(estimate, truth)$r))
  expect_equal(metric_r(estimate, truth)$n, 0L)
})

test_that("metrics error on inputs of different lengths", {
  expect_error(metric_bias(c(1, 2), c(1, 2, 3)), "same length")
  expect_error(metric_rmse(c(1, 2), c(1, 2, 3)), "same length")
  expect_error(metric_r(c(1, 2), c(1, 2, 3)), "same length")
  expect_error(metric_coverage(c(1, 2), c(1, 2, 3), c(1, 2, 3)), "same length")
})
