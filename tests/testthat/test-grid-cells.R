# Tests for recovery_grid(cells = ) and grid_tasks(), written against
# local/dev/spec-milestone-25-parallel.md (D115, D116; invariants I1--I5,
# I7).
#
# The contract: a grid run as a setup call and then as disjoint slices,
# each by its own call (on a VM, its own process), writes the cells one
# sequential run writes. Every slice below gets its own mock, as a
# separate process would, and runs under its own caller seed, so that
# nothing a slice inherits from the session can make it agree.

cells_grid <- function() {
  data.frame(n_subjects = c(3L, 4L, 3L), n_trials = c(10L, 20L, 15L))
}

cells_run <- function(dir, mock, cells = NULL, grid = cells_grid(),
                      reps = 2L, seed = 100, ...) {
  recovery_grid(
    bmm::mixture2p(resp_error = "y"),
    grid = grid,
    pars = c(kappa = log(8), thetat = stats::qlogis(0.75)),
    dir = dir,
    reps = reps,
    sds = c(kappa = 0.3, thetat = 0.4),
    seed = seed,
    chains = 2, iter = 400,
    cells = cells,
    ...,
    .fitter = mock$fitter
  )
}

#' Every file of one kind in `dir`, read back, named by file
#' @noRd
cell_files <- function(dir, pattern, drop = character()) {
  files <- sort(list.files(dir, pattern))
  stats::setNames(lapply(file.path(dir, files), function(path) {
    x <- readRDS(path)
    x[setdiff(names(x), drop)]
  }), files)
}

#' The sidecars without the fit's wall time, which no two runs share
#' @noRd
cell_sidecars <- function(dir) {
  cell_files(dir, "-est\\.rds$", drop = "fit_seconds")
}

#' A recovery as a plain data frame, without its attributes
#' @noRd
cells_rows <- function(x) {
  out <- as.data.frame(x)
  attributes(out) <- attributes(out)[c("names", "row.names", "class")]
  out
}

#' The setup call, then every slice in the order given, each with its
#' own mock and its own caller seed
#' @noRd
run_sliced <- function(dir, make_mock, n, order = seq_len(n), ...) {
  suppressMessages(cells_run(dir, make_mock(), cells = integer(0), ...))
  slices <- grid_tasks(cells_grid(), 2L, n)
  for (k in order) {
    withr::with_seed(1000 + k, suppressMessages(
      cells_run(dir, make_mock(), cells = slices[[k]], ...)
    ))
  }
  invisible(slices)
}

# grid_tasks() -----------------------------------------------------------------

test_that("grid_tasks() splits the cells into contiguous runs", {
  grid <- cells_grid()
  expect_identical(grid_tasks(grid, 2, 3), list(1:2, 3:4, 5:6))
  expect_identical(grid_tasks(grid, 2, 4), list(1:2, 3:4, 5L, 6L))
  expect_identical(grid_tasks(grid, 2, 1), list(1:6))
  expect_identical(grid_tasks(grid, 1, 3), list(1L, 2L, 3L))
  # sizes differ by at most one, and together they are every cell once
  tasks <- grid_tasks(data.frame(n_subjects = 1:7, n_trials = 1L), 3, 4)
  expect_identical(lengths(tasks), c(6L, 5L, 5L, 5L))
  expect_identical(unlist(tasks), 1:21)
})

test_that("grid_tasks() refuses what is not a split of the grid", {
  grid <- cells_grid()
  expect_error(grid_tasks(grid, 2, 7), "6 cells")
  expect_error(grid_tasks(grid, 2, 0), "n")
  expect_error(grid_tasks(grid, 0, 2), "reps")
  expect_error(grid_tasks(grid, 2, 1.5), "n")
  expect_error(grid_tasks(list(1), 2, 2), "grid")
})

# the setup call ---------------------------------------------------------------

test_that("the setup call writes the record and every sim, and fits nothing", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()

  expect_message(
    out <- cells_run(dir, mock, cells = integer(0)),
    "6 cells"
  )

  # the preflight, and no cell
  expect_identical(mock$calls$n, 1L)
  expect_true(file.exists(file.path(dir, "grid.rds")))
  expect_length(list.files(dir, "-sim\\.rds$"), 6L)
  expect_length(list.files(dir, "-est\\.rds$"), 0L)
  expect_s3_class(out, "tbl_df")
  expect_identical(nrow(out), 6L)
  expect_identical(out$condition, sprintf("row-%d", c(1:3, 1:3)))
  expect_identical(out$replication, rep(1:2, each = 3L))
  expect_identical(out$seed, grid_cell_table(3L, 2L, 100)$seed)
})

test_that("the setup call without a preflight calls no fitter at all", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  suppressMessages(cells_run(dir, mock, cells = integer(0), preflight = FALSE))
  expect_identical(mock$calls$n, 0L)
  expect_length(list.files(dir, "-sim\\.rds$"), 6L)
})

# I1: identity -----------------------------------------------------------------

test_that("setup and three slices write the cells one sequential run writes", {
  skip_if_not_installed("bmm")
  make_mock <- function() rng_mock_fitter(grid_mock_fitter())
  sequential <- withr::local_tempdir()
  sliced <- withr::local_tempdir()

  out <- suppressMessages(cells_run(sequential, make_mock()))
  run_sliced(sliced, make_mock, 3L, order = c(3L, 1L, 2L))

  expect_identical(
    cell_files(sliced, "-sim\\.rds$"), cell_files(sequential, "-sim\\.rds$")
  )
  sidecars <- cell_sidecars(sliced)
  expect_length(sidecars, 6L)
  expect_identical(sidecars, cell_sidecars(sequential))

  collected <- suppressMessages(collect_grid(sliced))
  expect_identical(
    cells_rows(collected),
    cells_rows(suppressMessages(collect_grid(sequential)))
  )
  expect_identical(collected$estimate, out$estimate)
  expect_identical(collected$covered, out$covered)
  expect_true(all(attr(collected, "cells")$status == "ok"))
  # the fitter's draws reach the estimates, so the comparison can fail
  expect_gt(length(unique(out$estimate)), 1L)
})

test_that("a slice simulates a missing cell as the full grid would", {
  # subjects = "fixed": rep 2 takes its subjects from rep 1 of its row.
  # With the sims deleted after setup, the last slice (rep 2 of rows 2
  # and 3) holds no rep 1 of either row and must draw it first.
  skip_if_not_installed("bmm")
  sequential <- withr::local_tempdir()
  sliced <- withr::local_tempdir()
  suppressMessages(cells_run(sequential, grid_mock_fitter(),
    subjects = "fixed"
  ))

  suppressMessages(cells_run(sliced, grid_mock_fitter(),
    cells = integer(0), subjects = "fixed"
  ))
  unlink(list.files(sliced, "-sim\\.rds$", full.names = TRUE))
  slices <- grid_tasks(cells_grid(), 2L, 3L)
  for (k in 3:1) {
    suppressMessages(cells_run(sliced, grid_mock_fitter(),
      cells = slices[[k]], subjects = "fixed"
    ))
  }

  sims <- cell_files(sliced, "-sim\\.rds$")
  expect_length(sims, 6L)
  expect_identical(sims, cell_files(sequential, "-sim\\.rds$"))
  expect_identical(cell_sidecars(sliced), cell_sidecars(sequential))
  expect_identical(
    sims[["cell-2-rep-1-sim.rds"]]$truth$subjects,
    sims[["cell-2-rep-2-sim.rds"]]$truth$subjects
  )
})

test_that("a full run after the setup call is the sequential run", {
  skip_if_not_installed("bmm")
  sequential <- withr::local_tempdir()
  staged <- withr::local_tempdir()
  out <- suppressMessages(cells_run(sequential, grid_mock_fitter()))
  suppressMessages(cells_run(staged, grid_mock_fitter(), cells = integer(0)))
  again <- suppressMessages(cells_run(staged, grid_mock_fitter()))
  expect_identical(cells_rows(again), cells_rows(out))
  expect_identical(cell_sidecars(staged), cell_sidecars(sequential))
})

# I3, I4: the record and the preflight belong to the setup call ---------------

test_that("a slice fits its cells only and never writes the grid record", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  suppressMessages(cells_run(dir, grid_mock_fitter(), cells = integer(0)))
  record <- file.path(dir, "grid.rds")
  before <- unname(tools::md5sum(record))
  mtime <- file.mtime(record)

  mock <- grid_mock_fitter()
  out <- suppressMessages(cells_run(dir, mock, cells = 3:4))

  # two cells, no preflight
  expect_identical(mock$calls$n, 2L)
  expect_identical(unname(tools::md5sum(record)), before)
  expect_identical(file.mtime(record), mtime)
  expect_length(list.files(dir, "-est\\.rds$"), 2L)
  # what it returns is its own cells
  cells <- attr(out, "cells")
  expect_identical(cells$condition, c("row-3", "row-1"))
  expect_identical(cells$replication, c(1L, 2L))
  expect_identical(cells$seed, grid_cell_table(3L, 2L, 100)$seed[3:4])
  expect_setequal(out$condition, c("row-3", "row-1"))
})

test_that("a slice without a setup call stops before any fit", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  expect_error(
    suppressMessages(cells_run(dir, mock, cells = 1:2)),
    "cells = integer\\(0\\)"
  )
  expect_identical(mock$calls$n, 0L)
  expect_false(file.exists(file.path(dir, "grid.rds")))
})

test_that("a slice of another grid stops and leaves the record alone", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  suppressMessages(cells_run(dir, grid_mock_fitter(), cells = integer(0)))
  before <- unname(tools::md5sum(file.path(dir, "grid.rds")))
  mock <- grid_mock_fitter()

  expect_error(
    suppressMessages(cells_run(dir, mock, cells = 1:2, seed = 101)),
    "seed"
  )
  expect_error(
    suppressMessages(cells_run(dir, mock, cells = 1:2, warmup = 100)),
    "sampler"
  )
  expect_identical(mock$calls$n, 0L)
  expect_identical(unname(tools::md5sum(file.path(dir, "grid.rds"))), before)
  expect_length(list.files(dir, "-est\\.rds$"), 0L)
})

# I5: resume -------------------------------------------------------------------

test_that("a slice run again fits only the cells it has no sidecar for", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  make_mock <- function() rng_mock_fitter(grid_mock_fitter())
  suppressMessages(cells_run(dir, make_mock(), cells = integer(0)))
  suppressMessages(cells_run(dir, make_mock(), cells = 1:3))
  first <- cell_sidecars(dir)

  again <- make_mock()
  suppressMessages(cells_run(dir, again, cells = 1:3))
  expect_identical(again$calls$n, 0L)

  # a slice killed during cell 2: its fit and sidecar never got their names
  unlink(file.path(dir, c("cell-2-rep-1.rds", "cell-2-rep-1-est.rds")))
  resumed <- make_mock()
  withr::with_seed(5, suppressMessages(cells_run(dir, resumed, cells = 1:3)))
  expect_identical(resumed$calls$n, 1L)
  expect_identical(cell_sidecars(dir), first)
})

# prior_sd: fitted by the setup call, read by the slices -----------------------

#' The grid mock, but a prior-only call returns a prior fit and is counted
#' @noRd
cells_prior_fitter <- function() {
  mock <- grid_mock_fitter()
  prior_calls <- new.env(parent = emptyenv())
  prior_calls$n <- 0L
  fitter <- function(formula, data, model, prior = NULL, ...) {
    if (identical(list(...)$sample_prior, "only")) {
      prior_calls$n <- prior_calls$n + 1L
      return(structure(
        list(prior = structure(data.frame(), sample_prior = "only")),
        class = c("cellspriorfit", "brmsfit")
      ))
    }
    mock$fitter(formula, data, model, prior = prior, ...)
  }
  list(fitter = fitter, calls = mock$calls, prior_calls = prior_calls)
}

registerS3method(
  "extract_estimates", "cellspriorfit",
  function(fit, level = "population", ...) {
    tibble::tibble(
      term = c("kappa", "thetat"), level = "population",
      post_sd_link = c(10, 20)
    )
  },
  envir = asNamespace("bmmtools")
)

test_that("slices read the setup call's prior SDs and never fit or drop them", {
  skip_if_not_installed("bmm")
  sequential <- withr::local_tempdir()
  sliced <- withr::local_tempdir()
  out <- suppressMessages(cells_run(sequential, cells_prior_fitter(),
    prior_sd = "fit"
  ))

  setup <- cells_prior_fitter()
  suppressMessages(cells_run(sliced, setup,
    cells = integer(0), prior_sd = "fit"
  ))
  expect_identical(setup$prior_calls$n, 1L)
  record <- file.path(sliced, "prior-sd.rds")
  before <- unname(tools::md5sum(record))

  slices <- grid_tasks(cells_grid(), 2L, 2L)
  for (k in 2:1) {
    mock <- cells_prior_fitter()
    part <- suppressMessages(cells_run(sliced, mock,
      cells = slices[[k]], prior_sd = "fit"
    ))
    expect_identical(mock$prior_calls$n, 0L)
    expect_false(all(is.na(part$contraction[part$level == "population"])))
  }
  expect_identical(unname(tools::md5sum(record)), before)

  collected <- suppressMessages(collect_grid(sliced))
  expect_identical(collected$contraction, out$contraction)
})

test_that("a slice whose prior SDs are gone from disk stops", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  suppressMessages(cells_run(dir, cells_prior_fitter(),
    cells = integer(0), prior_sd = "fit"
  ))
  unlink(file.path(dir, "prior-sd.rds"))
  mock <- cells_prior_fitter()
  expect_error(
    suppressMessages(cells_run(dir, mock, cells = 1:2, prior_sd = "fit")),
    "prior-sd.rds"
  )
  expect_identical(mock$calls$n, 0L)
  expect_identical(mock$prior_calls$n, 0L)
  # a refused slice leaves no run record: the setup call's is the only one
  expect_length(list.files(file.path(dir, "run-info")), 1L)
})

test_that("a slice with another prior, control or prior_sd stops", {
  # what is not a field of the record still decides a fit: the setup call
  # records it as one key, and a slice is checked against that
  skip_if_not_installed("bmm")
  skip_if_not_installed("brms")
  dir <- withr::local_tempdir()
  suppressMessages(cells_run(dir, cells_prior_fitter(),
    cells = integer(0), prior_sd = "fit"
  ))
  mock <- cells_prior_fitter()
  run <- function(...) {
    suppressMessages(cells_run(dir, mock, cells = 1:2, ...))
  }
  # the text only a change of `fit` shows
  fit_changed <- "covers the prior"
  expect_error(
    run(prior_sd = "fit", prior = brms::prior(normal(0, 1), class = "b")),
    fit_changed
  )
  expect_error(
    run(prior_sd = "fit", control = list(adapt_delta = 0.99)), fit_changed
  )
  expect_error(run(), fit_changed)
  expect_identical(mock$calls$n, 0L)
  expect_identical(mock$prior_calls$n, 0L)
  # fewer cores change no draw, so a slice may name them
  run(prior_sd = "fit", cores = 1)
  expect_identical(mock$calls$n, 2L)
})

test_that("a record a call without `cells` wrote cannot be sliced", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  suppressMessages(cells_run(dir, grid_mock_fitter()))
  expect_error(
    suppressMessages(cells_run(dir, grid_mock_fitter(), cells = 1:2)),
    "not written by a setup call"
  )
  saveRDS("not a record", file.path(dir, "grid.rds"))
  expect_error(
    suppressMessages(cells_run(dir, grid_mock_fitter(), cells = 1:2)),
    "cannot be read"
  )
})

test_that("subject-level slices need no prior SDs", {
  # with no level that has a contraction, the setup call writes no
  # prior-sd.rds, and its slices must not ask for one
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  run <- function(cells) {
    suppressMessages(cells_run(dir, cells_prior_fitter(),
      cells = cells, prior_sd = "fit", levels = "subject"
    ))
  }
  run(integer(0))
  expect_false(file.exists(file.path(dir, "prior-sd.rds")))
  out <- run(1:3)
  expect_setequal(out$level, "subject")
})

test_that("cells in any order run the cells they name", {
  skip_if_not_installed("bmm")
  make_mock <- function() rng_mock_fitter(grid_mock_fitter())
  sequential <- withr::local_tempdir()
  sliced <- withr::local_tempdir()
  suppressMessages(cells_run(sequential, make_mock()))
  suppressMessages(cells_run(sliced, make_mock(), cells = integer(0)))
  out <- suppressMessages(cells_run(sliced, make_mock(), cells = c(5L, 2L)))
  expect_identical(attr(out, "cells")$condition, c("row-2", "row-2"))
  expect_identical(attr(out, "cells")$replication, c(2L, 1L))
  expect_identical(
    cell_sidecars(sliced),
    cell_sidecars(sequential)[c("cell-2-rep-1-est.rds", "cell-2-rep-2-est.rds")]
  )
})

# what `cells` accepts ---------------------------------------------------------

test_that("cells must name cells of the grid, each once", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  run <- function(cells, ...) {
    suppressMessages(cells_run(dir, mock, cells = cells, ...))
  }
  expect_error(run(0L), "between 1 and 6")
  expect_error(run(7L), "between 1 and 6")
  expect_error(run(c(1L, NA)), "between 1 and 6")
  expect_error(run(1.5), "whole")
  expect_error(run(c(1L, 1L)), "once")
  expect_error(run("1"), "whole")
  expect_error(run(1:2, smoke = TRUE), "smoke")
  expect_identical(mock$calls$n, 0L)
})
