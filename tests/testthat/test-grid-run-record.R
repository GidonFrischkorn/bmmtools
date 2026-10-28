# Tests for the run record a grid call leaves in its directory, written
# against local/dev/spec-milestone-25-parallel.md (D117; invariant I6).
#
# One record per call that ran, so that a grid split over processes says
# which machine and which job ran which cells. `run_info()` itself keeps
# its contract; the job columns live in the grid's record only.

record_run <- function(dir, mock = grid_mock_fitter(), cells = NULL, ...) {
  suppressMessages(recovery_grid(
    bmm::mixture2p(resp_error = "y"),
    grid = small_grid(),
    pars = c(kappa = log(8), thetat = stats::qlogis(0.75)),
    dir = dir,
    reps = 2L,
    sds = c(kappa = 0.3, thetat = 0.4),
    seed = 100,
    chains = 2, iter = 400,
    cells = cells,
    ...,
    .fitter = mock$fitter
  ))
}

no_scheduler <- function(env = parent.frame()) {
  withr::local_envvar(
    SLURM_JOB_ID = NA, SLURM_ARRAY_JOB_ID = NA, SLURM_ARRAY_TASK_ID = NA,
    SLURM_CPUS_PER_TASK = NA, PBS_JOBID = NA, PBS_ARRAY_INDEX = NA,
    PBS_ARRAYID = NA, NCPUS = NA, JOB_ID = NA, SGE_TASK_ID = NA,
    NSLOTS = NA, LSB_JOBID = NA, LSB_JOBINDEX = NA, LSB_DJOB_NUMPROC = NA,
    .local_envir = env
  )
}

# the job columns --------------------------------------------------------------

test_that("grid_job_info() reads no scheduler as NA", {
  no_scheduler()
  out <- grid_job_info()
  expect_s3_class(out, "tbl_df")
  expect_named(out, c(
    "scheduler", "job_id", "array_job_id", "array_task_id", "cores_granted"
  ))
  expect_identical(out$scheduler, NA_character_)
  expect_identical(out$job_id, NA_character_)
  expect_identical(out$cores_granted, NA_integer_)
})

test_that("grid_job_info() reads SLURM, PBS, SGE and LSF", {
  no_scheduler()
  withr::with_envvar(c(
    SLURM_JOB_ID = "812", SLURM_ARRAY_JOB_ID = "800",
    SLURM_ARRAY_TASK_ID = "12", SLURM_CPUS_PER_TASK = "8"
  ), {
    out <- grid_job_info()
    expect_identical(out$scheduler, "slurm")
    expect_identical(out$job_id, "812")
    expect_identical(out$array_job_id, "800")
    expect_identical(out$array_task_id, "12")
    expect_identical(out$cores_granted, 8L)
  })
  withr::with_envvar(c(
    PBS_JOBID = "77[3].server", PBS_ARRAY_INDEX = "3", NCPUS = "4"
  ), {
    out <- grid_job_info()
    expect_identical(out$scheduler, "pbs")
    expect_identical(out$array_task_id, "3")
    expect_identical(out$cores_granted, 4L)
  })
  withr::with_envvar(c(JOB_ID = "55", SGE_TASK_ID = "2", NSLOTS = "6"), {
    out <- grid_job_info()
    expect_identical(out$scheduler, "sge")
    expect_identical(out$array_task_id, "2")
    expect_identical(out$cores_granted, 6L)
  })
  withr::with_envvar(c(
    LSB_JOBID = "9", LSB_JOBINDEX = "0", LSB_DJOB_NUMPROC = "2"
  ), {
    out <- grid_job_info()
    expect_identical(out$scheduler, "lsf")
    # LSF's index is 0 outside an array
    expect_identical(out$array_task_id, NA_character_)
    expect_identical(out$cores_granted, 2L)
  })
  # an unreadable core count is NA, not an error
  withr::with_envvar(c(SLURM_JOB_ID = "1", SLURM_CPUS_PER_TASK = "lots"), {
    expect_identical(grid_job_info()$cores_granted, NA_integer_)
  })
})

test_that("run_info() keeps its columns", {
  # D117: the job columns are the grid record's, not run_info()'s
  expect_false(any(c("scheduler", "job_id") %in% names(run_info())))
})

# the record a call writes -----------------------------------------------------

test_that("a full run writes one record and returns it", {
  skip_if_not_installed("bmm")
  no_scheduler()
  dir <- withr::local_tempdir()
  out <- record_run(dir)

  files <- list.files(file.path(dir, "run-info"), "\\.rds$")
  expect_length(files, 1L)
  record <- attr(out, "run_info")
  expect_s3_class(record, "tbl_df")
  expect_identical(nrow(record), 1L)
  expect_true(all(names(run_info()) %in% names(record)))
  expect_identical(record$call, "all")
  expect_identical(record$cells, "1-4")
  expect_identical(record$n_cells, 4L)
  expect_identical(record$n_reused, 0L)
  expect_identical(record$n_failed, 0L)
  expect_identical(record$pid, Sys.getpid())
  expect_s3_class(record$started, "POSIXct")
  expect_false(is.na(record$ended))
  expect_true(record$ended >= record$started)
  expect_identical(record$scheduler, NA_character_)
})

test_that("setup and two slices leave three records, read back as one table", {
  skip_if_not_installed("bmm")
  no_scheduler()
  dir <- withr::local_tempdir()
  record_run(dir, cells = integer(0))
  withr::with_envvar(c(
    SLURM_JOB_ID = "11", SLURM_ARRAY_JOB_ID = "10",
    SLURM_ARRAY_TASK_ID = "2", SLURM_CPUS_PER_TASK = "4"
  ), record_run(dir, cells = 3:4))
  record_run(dir, cells = c(1L, 2L))

  collected <- suppressMessages(collect_grid(dir))
  record <- attr(collected, "run_info")
  expect_identical(nrow(record), 3L)
  # in the order the calls started
  expect_identical(record$call, c("setup", "slice", "slice"))
  expect_identical(record$cells, c("", "3-4", "1-2"))
  expect_identical(record$n_cells, c(0L, 2L, 2L))
  expect_identical(record$array_task_id, c(NA, "2", NA))
  expect_identical(record$cores_granted, c(NA, 4L, NA))
})

test_that("records are read in the order their calls started", {
  # file names sort by start time on one machine only; on several, the
  # host part decides, so the order comes from `started`
  skip_if_not_installed("bmm")
  no_scheduler()
  dir <- withr::local_tempdir()
  record_run(dir, cells = integer(0))
  record_run(dir, cells = 1:2)
  folder <- file.path(dir, "run-info")
  last <- sort(list.files(folder))[[2L]]
  file.rename(file.path(folder, last), file.path(folder, "0-first-by-name.rds"))
  # cells 3 and 4 have no results yet, and collect_grid() warns so
  collected <- suppressWarnings(suppressMessages(collect_grid(dir)))
  record <- attr(collected, "run_info")
  expect_identical(record$call, c("setup", "slice"))
})

test_that("a slice run again counts its reused cells", {
  skip_if_not_installed("bmm")
  no_scheduler()
  dir <- withr::local_tempdir()
  record_run(dir, cells = integer(0))
  record_run(dir, cells = 1:2)
  again <- record_run(dir, cells = 1:2)
  record <- attr(again, "run_info")
  expect_identical(nrow(record), 3L)
  expect_identical(record$n_reused, c(0L, 0L, 2L))
})

test_that("cells that fail are counted, and the record is still closed", {
  skip_if_not_installed("bmm")
  no_scheduler()
  dir <- withr::local_tempdir()
  record_run(dir, cells = integer(0))
  failing <- grid_mock_fitter(fail_on = function(n, data, dots) n <= 2L)
  # the grid warns that two cells failed, which is what this test makes
  out <- suppressWarnings(record_run(dir, failing, cells = 1:3))
  record <- attr(out, "run_info")
  last <- record[nrow(record), ]
  expect_identical(last$n_failed, 2L)
  expect_false(is.na(last$ended))
})

test_that("a call that stops before its cells writes no record", {
  skip_if_not_installed("bmm")
  no_scheduler()
  dir <- withr::local_tempdir()
  expect_error(record_run(dir, cells = 1:2), "integer\\(0\\)")
  expect_false(dir.exists(file.path(dir, "run-info")))
})

test_that("a record opened by a call that never finished reads as unfinished", {
  # the record is written when the cells start and again when they end,
  # so a process killed in between leaves one with `ended` NA
  skip_if_not_installed("bmm")
  no_scheduler()
  dir <- withr::local_tempdir()
  record_run(dir, cells = integer(0))
  record_run(dir, cells = 3:4)
  killer <- grid_mock_fitter(fail_on = function(n, data, dots) {
    rlang::interrupt()
  })
  tryCatch(record_run(dir, killer, cells = 1:2), interrupt = function(e) NULL)
  # the interrupted cells have no results, and collect_grid() warns so
  collected <- suppressWarnings(suppressMessages(collect_grid(dir)))
  record <- attr(collected, "run_info")
  expect_identical(nrow(record), 3L)
  expect_identical(record$cells, c("", "3-4", "1-2"))
  expect_false(is.na(record$ended[[2L]]))
  expect_true(is.na(record$ended[[3L]]))
  expect_true(is.na(record$n_failed[[3L]]))
})

test_that("an error outside a cell's fit leaves the record unfinished", {
  # data that cannot be generated stops the grid once its cells have
  # started; a failed fit would be counted instead
  skip_if_not_installed("bmm")
  no_scheduler()
  dir <- withr::local_tempdir()
  expect_error(suppressMessages(recovery_grid(
    bmm::mixture2p(resp_error = "y"),
    grid = small_grid(),
    pars = function(row) {
      if (row$n_subjects == 4L) stop("no values for row 2")
      c(kappa = log(8), thetat = stats::qlogis(0.75))
    },
    dir = dir, reps = 2L, sds = c(kappa = 0.3, thetat = 0.4), seed = 100,
    chains = 2, iter = 400, .fitter = grid_mock_fitter()$fitter
  )), "no values for row 2")
  record <- read_run_records(dir)
  expect_identical(nrow(record), 1L)
  expect_identical(record$call, "all")
  expect_true(is.na(record$ended))
  expect_true(is.na(record$n_failed))
  expect_true(is.na(record$n_reused))
})

test_that("a directory without records collects with no run_info", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  record_run(dir)
  unlink(file.path(dir, "run-info"), recursive = TRUE)
  expect_null(attr(suppressMessages(collect_grid(dir)), "run_info"))
})

test_that("a call whose preflight fails leaves no record", {
  # a call that errors before its cells leaves no record; only one stopped
  # once they had started leaves an unfinished one
  skip_if_not_installed("bmm")
  no_scheduler()
  failing <- grid_mock_fitter(fail_on = function(n, data, dots) TRUE)
  dir <- withr::local_tempdir()
  expect_error(record_run(dir, failing), "preflight")
  expect_false(dir.exists(file.path(dir, "run-info")))
  dir <- withr::local_tempdir()
  expect_error(record_run(dir, failing, cells = integer(0)), "preflight")
  expect_false(dir.exists(file.path(dir, "run-info")))
})

test_that("a record starts when its call starts, preflight included", {
  skip_if_not_installed("bmm")
  no_scheduler()
  mock <- grid_mock_fitter()
  seen <- new.env(parent = emptyenv())
  timed <- list(fitter = function(...) {
    if (is.null(seen$first)) seen$first <- Sys.time()
    mock$fitter(...)
  })
  dir <- withr::local_tempdir()
  record <- attr(record_run(dir, timed), "run_info")
  expect_true(record$started <= seen$first)
})

test_that("records that do not bind are left out, with a warning", {
  dir <- withr::local_tempdir()
  folder <- file.path(dir, "run-info")
  dir.create(folder)
  row <- function(started, cores) {
    tibble::tibble(call = "slice", started = started, cores = cores)
  }
  saveRDS(row(Sys.time(), 4L), file.path(folder, "a.rds"))
  saveRDS(row(Sys.time() + 1, "four"), file.path(folder, "b.rds"))
  expect_warning(out <- read_run_records(dir), "cannot be read as one table")
  expect_null(out)
})

test_that("a record's file name survives a host name of any encoding", {
  started <- as.POSIXct("2026-10-27 10:11:12.345", tz = "UTC")
  # a latin1 byte in a session that reads strings as UTF-8: invalid text
  host <- "caf\xe9 ../node"
  path <- run_record_path("grid", started, host, 42L)
  expect_identical(dirname(path), file.path("grid", "run-info"))
  expect_match(
    basename(path), "^20261027-101112\\.345-[A-Za-z0-9._-]+-42\\.rds$"
  )
})
