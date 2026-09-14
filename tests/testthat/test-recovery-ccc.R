# Tests for recovery_ccc(), the exported concordance
# (local/dev/spec-milestone-1-score-layer.md, Amendment 2026-09-14, A4).
# The arithmetic is tested through metric_ccc() in test-metrics.R; these
# tests pin the public shape and the argument checks.

test_that("recovery_ccc returns one row with the documented columns", {
  withr::local_seed(41)
  truth <- stats::rnorm(30)
  estimate <- 0.8 * truth + stats::rnorm(30, sd = 0.4)
  out <- recovery_ccc(estimate, truth)

  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 1L)
  expect_named(out, c(
    "ccc", "ccc_low", "ccc_high", "ccc_accuracy", "ccc_scale_shift",
    "ccc_location_shift", "calibration_slope", "n"
  ))
  expect_identical(out$n, 30L)

  inner <- metric_ccc(estimate, truth)
  expect_equal(out$ccc, inner$ccc)
  expect_equal(out$ccc_low, inner$ccc_low)
  expect_equal(out$ccc_scale_shift, inner$scale_shift)
})

test_that("recovery_ccc drops incomplete pairs and counts the rest", {
  out <- recovery_ccc(c(1, 2, NA, 4, 5, 7), c(1, 3, 3, NA, 5, 6))
  expect_identical(out$n, 4L)
  expect_false(is.na(out$ccc))
})

test_that("recovery_ccc rejects non-numeric input, naming the argument", {
  expect_error(recovery_ccc(letters[1:4], 1:4), "estimate")
  expect_error(recovery_ccc(1:4, letters[1:4]), "truth")
})

test_that("recovery_ccc rejects vectors of different length", {
  expect_error(recovery_ccc(1:4, 1:5), "same length")
})
