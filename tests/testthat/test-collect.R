# Tests for collect_grid(), written against
# local/dev/spec-milestone-9-study-support.md section 9.2.
#
# The point of the function is that it needs neither the fits nor the
# versions that made them, so every test deletes the fits first.

collect_run <- function(dir, ...) {
  suppressMessages(recovery_grid(
    bmm::mixture2p(resp_error = "y"),
    grid = small_grid(),
    pars = c(kappa = log(8), thetat = stats::qlogis(0.75)),
    dir = dir,
    reps = 2L,
    sds = c(kappa = 0.3, thetat = 0.4),
    seed = 100,
    chains = 2, iter = 400,
    ...,
    .fitter = grid_mock_fitter()$fitter
  ))
}

# the recovery rows alone, with the attributes that legitimately differ
# (the cells table's `elapsed`) left out of the comparison
strip_attrs <- function(x) {
  out <- as.data.frame(x)
  attributes(out) <- attributes(out)[c("names", "row.names", "class")]
  out
}

delete_fits <- function(dir) {
  fits <- setdiff(
    list.files(dir, pattern = "^cell-.*\\.rds$", full.names = TRUE),
    list.files(dir, pattern = "(-sim|-est)\\.rds$", full.names = TRUE)
  )
  unlink(fits)
  invisible(fits)
}

test_that("a grid is rebuilt from its cell files with the fits deleted", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  original <- collect_run(dir)
  # the fits and their meta files both go: nothing but -sim and -est is left
  expect_length(delete_fits(dir), 8L)

  out <- collect_grid(dir)

  expect_s3_class(out, "bmmtools_recovery")
  expect_identical(strip_attrs(out), strip_attrs(original))
  expect_equal(summary(out), summary(original))

  cells <- attr(out, "cells")
  before <- attr(original, "cells")
  expect_named(cells, names(before))
  expect_identical(cells$condition, before$condition)
  expect_identical(cells$status, before$status)
  expect_identical(cells$fit_seconds, before$fit_seconds)
  expect_identical(cells$chains, before$chains)
  # no cell was run, so no cell took any time
  expect_true(all(is.na(cells$elapsed)))
  expect_identical(attr(out, "grid"), attr(original, "grid"))
})

test_that("a cell with no result is reported, not silently dropped", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  collect_run(dir)
  delete_fits(dir)
  unlink(file.path(dir, "cell-1-rep-2-est.rds"))

  expect_warning(out <- collect_grid(dir), "row-1 rep 2")

  cells <- attr(out, "cells")
  expect_identical(sum(cells$status == "missing"), 1L)
  expect_identical(cells$status[cells$condition == "row-1"], c("ok", "missing"))
  expect_false(any(out$condition == "row-1" & out$replication == 2L))
})

test_that("a directory with no record says so and names the file", {
  dir <- withr::local_tempdir()
  err <- expect_error(collect_grid(dir))
  expect_match(conditionMessage(err), "grid.rds")
  expect_match(conditionMessage(err), "before bmmtools 0.2.0")
})

test_that("collecting re-scores on another scale without any fit", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  collect_run(dir, scale = "link")
  delete_fits(dir)

  link <- collect_grid(dir)
  natural <- collect_grid(dir, scale = "natural")

  expect_identical(unique(link$scale), "link")
  expect_identical(unique(natural$scale), "natural")
  expect_false(isTRUE(all.equal(link$estimate, natural$estimate)))
})

test_that("asking for what the cells do not hold is an error", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  collect_run(dir, levels = "population")
  delete_fits(dir)

  err <- expect_error(collect_grid(dir, levels = c("population", "subject")))
  expect_match(conditionMessage(err), "not extracted")
  err_cor <- expect_error(collect_grid(dir, correlations = "model"))
  expect_match(conditionMessage(err_cor), "correlation estimator")

  # narrowing to a subset of what is there is fine
  expect_s3_class(collect_grid(dir, levels = "population"), "bmmtools_recovery")
})

test_that("a bad dir is refused before anything is read", {
  expect_error(collect_grid(1), "single path")
  expect_error(collect_grid(c("a", "b")), "single path")
})
