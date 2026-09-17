# Tests for the example objects in data/.
#
# They are built by data-raw/example-objects.R with Stan and saved, so a
# change to either class contract would leave them stale without failing
# anything else. These tests are that failure.

test_that("recovery_mixture2p satisfies the recovery contract", {
  x <- recovery_mixture2p
  expect_s3_class(x, "bmmtools_recovery")
  expect_true(all(recovery_contract_columns() %in% names(x)))
  expect_setequal(unique(x$level), c("population", "subject"))
  expect_false(anyNA(x$condition))
  expect_s3_class(attr(x, "cells"), "tbl_df")
})

test_that("summary() of recovery_mixture2p groups by condition", {
  s <- summary(recovery_mixture2p)
  expect_s3_class(s, "bmmtools_recovery_summary")
  expect_identical(names(s), c("condition", recovery_summary_columns()))
  expect_setequal(s$condition, unique(recovery_mixture2p$condition))
})

test_that("prior_check_sdt_yn satisfies the prior-check contract", {
  x <- prior_check_sdt_yn
  expect_s3_class(x, "bmmtools_prior_check")
  expect_true(all(names(prior_check_contract()) %in% names(x)))
  expect_identical(attr(x, "sets"), c("default", "narrow_sd"))
  s <- summary(x)
  expect_true("difference" %in% names(s))
})

test_that("the example objects plot", {
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot_recovery(recovery_mixture2p), "ggplot")
  expect_s3_class(plot_prior_check(prior_check_sdt_yn), "ggplot")
})

test_that("no example object carries a local path", {
  for (x in list(recovery_mixture2p, prior_check_sdt_yn)) {
    text <- paste(utils::capture.output(utils::str(attributes(x))),
      collapse = "\n"
    )
    expect_false(grepl("/Users/|/home/|[A-Z]:\\\\", text))
  }
})
