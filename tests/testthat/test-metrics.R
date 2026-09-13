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
  # scale_bias is the estimates' SD over truth's, here equal
  expect_equal(out$scale_bias, 1)
  # location_bias is Lin's u = (mean_y - mean_x) / sqrt(sd_x * sd_y),
  # so 1 / sqrt(1.25) = 0.8944272, NOT 1 / 1.25. The denominator is the
  # geometric mean of the two standard deviations, not of the variances;
  # the decomposition test below is what pins this down.
  expect_equal(out$location_bias, 1 / sqrt(1.25))
})

test_that("ccc decomposes as r times Lin's bias-correction factor", {
  # ccc = r * C_b with C_b = 2 / (v + 1/v + u^2), v = scale_bias,
  # u = location_bias. This identity holds only for Lin's definitions of
  # u and v, so it is what fixes them: no reference implementation is
  # installed, and a hand-typed constant would only restate the code.
  withr::local_seed(11)
  for (i in 1:5) {
    truth <- stats::rnorm(20)
    estimate <- truth * stats::runif(1, 0.5, 2) +
      stats::rnorm(1) + stats::rnorm(20, sd = 0.3)

    out <- metric_ccc(estimate, truth)
    r <- metric_r(estimate, truth)$r
    v <- out$scale_bias
    u <- out$location_bias
    c_b <- 2 / (v + 1 / v + u^2)

    expect_equal(out$ccc, r * c_b, tolerance = 1e-10)
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
  expect_equal(metric_ccc(stretched, truth)$scale_bias, 3)
})

test_that("ccc is negative for a reversed relationship", {
  truth <- c(1, 2, 3, 4)
  expect_lt(metric_ccc(-truth, truth)$ccc, 0)
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
  expect_true(is.na(metric_r(estimate, truth)$r))
  expect_equal(metric_r(estimate, truth)$n, 0L)
})

test_that("metrics error on inputs of different lengths", {
  expect_error(metric_bias(c(1, 2), c(1, 2, 3)), "same length")
  expect_error(metric_rmse(c(1, 2), c(1, 2, 3)), "same length")
  expect_error(metric_r(c(1, 2), c(1, 2, 3)), "same length")
  expect_error(metric_coverage(c(1, 2), c(1, 2, 3), c(1, 2, 3)), "same length")
})
