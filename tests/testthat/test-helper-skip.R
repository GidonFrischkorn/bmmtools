# Tests for skip_if_no_bmm_model() (helper-skip.R) itself: a guard that
# is silently wrong is worse than no guard, per clean-runner-checks-before-
# release. testthat::skip() throws a condition of class "skip" rather than
# aborting the run, so it is caught rather than triggered.

run <- function(...) {
  tryCatch(
    {
      skip_if_no_bmm_model(...)
      "ran"
    },
    skip = function(e) e
  )
}

test_that("a model whose constructor and generator both exist does not skip", {
  skip_if_not_installed("bmm")
  # mixture2p ships in every bmm build measured so far (CRAN and develop).
  expect_identical(run("mixture2p"), "ran")
})

test_that("a name the installed bmm does not export skips, naming it", {
  skip_if_not_installed("bmm")
  cond <- run("no_such_bmm_model_xyz")
  expect_s3_class(cond, "skip")
  expect_match(conditionMessage(cond), "no_such_bmm_model_xyz")
})

test_that("density = TRUE also requires the density function", {
  skip_if_not_installed("bmm")
  cond <- run("no_such_bmm_model_xyz", density = TRUE)
  expect_s3_class(cond, "skip")
  expect_match(conditionMessage(cond), "dno_such_bmm_model_xyz")
})
