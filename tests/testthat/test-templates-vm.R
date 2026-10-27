# Tests for the VM launcher in inst/templates/vm (D118, D119).
#
# The templates are files a user copies, so the suite checks what can
# be checked without a grid: they ship, they parse, `grid-slice.R info`
# does the arithmetic the launcher relies on, and the launcher refuses
# more slices than the machine has cores for. A run of the launcher
# itself, on a real grid, is a scratch check (local/dev/STATE-milestone-25.md,
# 25.3); here it runs on a stub `grid-slice.R`, which needs no bmm.

vm_template <- function(file) {
  system.file("templates", "vm", file, package = "bmmtools", mustWork = TRUE)
}

#' A grid-args file that needs no package: what `info` reads
#' @noRd
write_args <- function(dir, chains = NULL, threads = NULL, reps = 2) {
  path <- file.path(dir, "grid-args.R")
  writeLines(c(
    "grid_args <- list(",
    "  grid = data.frame(n_subjects = c(20, 50, 80), n_trials = 30),",
    sprintf("  reps = %d,", reps),
    if (!is.null(chains)) sprintf("  chains = %d,", chains),
    if (!is.null(threads)) {
      sprintf("  threads = list(threads = %d),", threads)
    },
    "  dir = 'fits/grid'",
    ")"
  ), path)
  path
}

slice_info <- function(args) {
  out <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(shQuote(vm_template("grid-slice.R")), shQuote(args), "info"),
    stdout = TRUE, stderr = TRUE
  )
  strsplit(trimws(out[[length(out)]]), " ")[[1L]]
}

test_that("the three templates ship and the R ones parse", {
  for (file in c("grid-args.R", "grid-slice.R", "launch.sh")) {
    expect_true(file.exists(vm_template(file)))
  }
  expect_no_error(parse(vm_template("grid-args.R")))
  expect_no_error(parse(vm_template("grid-slice.R")))
})

test_that("launch.sh is valid bash", {
  skip_on_os("windows") # the launcher is a bash script for Linux and macOS
  skip_if(Sys.which("bash") == "", "no bash")
  status <- system2("bash", c("-n", shQuote(vm_template("launch.sh"))))
  expect_identical(status, 0L)
})

test_that("grid-slice.R info gives cores per slice, cells and dir", {
  dir <- withr::local_tempdir()
  # brms's default of 4 chains, one thread each
  expect_identical(slice_info(write_args(dir)), c("4", "6", "fits/grid"))
  expect_identical(
    slice_info(write_args(dir, chains = 2, threads = 3)),
    c("6", "6", "fits/grid")
  )
  # never 3e+05, which the launcher cannot compare
  expect_identical(
    slice_info(write_args(dir, reps = 100000L)), c("4", "300000", "fits/grid")
  )
})

test_that("launch.sh refuses more slices than the cores hold", {
  skip_on_os("windows") # the launcher is a bash script for Linux and macOS
  skip_if(Sys.which("bash") == "", "no bash")
  # R CMD check puts an Rscript on PATH that only says to give a path
  withr::local_path(R.home("bin"), action = "prefix")
  dir <- withr::local_tempdir()
  args <- write_args(dir, chains = 1)
  run <- function(...) {
    # system2() warns on a non-zero exit, which the tests read as `status`
    withr::with_dir(dir, suppressWarnings(system2(
      "bash", c(shQuote(vm_template("launch.sh")), ..., shQuote(args)),
      stdout = TRUE, stderr = TRUE
    )))
  }
  # more slices than any machine has cores
  out <- run("-k", 100000L)
  expect_identical(attr(out, "status"), 2L)
  expect_match(paste(out, collapse = "\n"), "is more than [0-9]+ cores")
  expect_match(paste(out, collapse = "\n"), "or --force")
  # nothing was started
  expect_false(dir.exists(file.path(dir, "fits", "grid", "logs")))
})

# the launcher on a stub ------------------------------------------------------

#' A launcher beside a stub grid-slice.R, in folders with spaces
#'
#' The stub prints a line before its `info` answer, as a grid-args.R
#' that prints would; its setup makes the grid folder; slice k writes
#' `cell-k-est.rds`, exits 3 when k is in `fail`, and waits while a file
#' `hold-k` exists.
#' @noRd
stub_launcher <- function(root, fail = integer(0)) {
  vm <- file.path(root, "my vm")
  dir.create(vm)
  file.copy(vm_template("launch.sh"), vm)
  writeLines(c(
    "args <- commandArgs(trailingOnly = TRUE)",
    "source(args[[1L]])",
    "dir <- grid_args$dir",
    "if (args[[2L]] == 'info') {",
    "  cat('a line from grid-args.R\\n')",
    "  cat(1L, 6L, dir, '\\n')",
    "  quit(save = 'no')",
    "}",
    "if (args[[2L]] == 'setup') {",
    "  dir.create(dir, recursive = TRUE, showWarnings = FALSE)",
    "  quit(save = 'no')",
    "}",
    "k <- as.integer(args[[2L]])",
    "while (file.exists(file.path(dir, paste0('hold-', k)))) Sys.sleep(0.1)",
    "if (k %in% grid_args$fail) quit(save = 'no', status = 3L)",
    "file.create(file.path(dir, sprintf('cell-%d-est.rds', k)))"
  ), file.path(vm, "grid-slice.R"))
  args <- file.path(vm, "grid-args.R")
  writeLines(sprintf(
    "grid_args <- list(dir = 'fits/my grid', fail = %s)", deparse(fail)
  ), args)
  list(
    launch = file.path(vm, "launch.sh"), args = args,
    dir = file.path(root, "fits", "my grid"),
    logs = file.path(root, "fits", "my grid", "logs")
  )
}

run_launcher <- function(root, launcher, ...) {
  # system2() warns on a non-zero exit, which the tests read as `status`
  withr::with_dir(root, suppressWarnings(system2(
    "bash", c(shQuote(launcher$launch), ..., shQuote(launcher$args)),
    stdout = TRUE, stderr = TRUE
  )))
}

#' The last line of a slice's log once it has an exit line, or NULL
#' @noRd
exit_line <- function(launcher, k, wait = 60) {
  log <- file.path(launcher$logs, sprintf("slice-%d.log", k))
  deadline <- Sys.time() + wait
  while (Sys.time() < deadline) {
    lines <- if (file.exists(log)) readLines(log, warn = FALSE) else ""
    end <- grep("^EXIT [0-9]+$", lines, value = TRUE)
    if (length(end) > 0L) {
      return(end[[length(end)]])
    }
    Sys.sleep(0.2)
  }
  NULL
}

#' The process id the launcher recorded for slice k, last start first
#' @noRd
pid_of_slice <- function(launcher, k) {
  pids <- read.table(file.path(launcher$logs, "pids"))
  utils::tail(pids[[2L]][pids[[1L]] == k], 1L)
}

skip_without_launcher <- function() {
  skip_on_os("windows") # the launcher is a bash script for Linux and macOS
  skip_if(Sys.which("bash") == "", "no bash")
}

test_that("launch.sh runs the setup call, then one logged slice per k", {
  skip_without_launcher()
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  launcher <- stub_launcher(root, fail = 2L)
  out <- run_launcher(root, launcher, "--force", "-k", "3")
  expect_null(attr(out, "status"))
  expect_true(any(grepl("started 3 slices over 6 cells", out)))
  expect_identical(
    utils::tail(readLines(file.path(launcher$logs, "setup.log")), 1L),
    "EXIT 0"
  )
  expect_identical(exit_line(launcher, 1L), "EXIT 0")
  expect_identical(exit_line(launcher, 2L), "EXIT 3")
  expect_identical(exit_line(launcher, 3L), "EXIT 0")
  pids <- read.table(file.path(launcher$logs, "pids"))
  expect_identical(pids[[1L]], 1:3)

  status <- run_launcher(root, launcher, "--status")
  expect_true("cells with a result: 2 of 6" %in% status)
  expect_true("slice-2.log: EXIT 3" %in% status)
  expect_identical(
    grep("^slice-", status, value = TRUE),
    c("slice-1.log: EXIT 0", "slice-2.log: EXIT 3", "slice-3.log: EXIT 0")
  )
})

test_that("launch.sh takes whole numbers for -k and --rerun only", {
  skip_without_launcher()
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  launcher <- stub_launcher(root)
  for (bad in c("abc", "0", "2.5", "-1")) {
    out <- run_launcher(root, launcher, "--force", "-k", bad)
    expect_identical(attr(out, "status"), 2L)
    expect_true(any(grepl("-k needs a whole number", out)))
  }
  out <- run_launcher(root, launcher, "--rerun", "x")
  expect_identical(attr(out, "status"), 2L)
  # nothing ran: no setup call, no grid folder
  expect_false(dir.exists(launcher$dir))
  # a rerun needs an earlier launch
  out <- run_launcher(root, launcher, "--rerun", "1")
  expect_identical(attr(out, "status"), 2L)
  expect_true(any(grepl("no earlier launch", out)))
})

test_that("a stopped slice says so, can be run again, and blocks no launch", {
  skip_without_launcher()
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  launcher <- stub_launcher(root)
  dir.create(launcher$dir, recursive = TRUE)
  hold <- file.path(launcher$dir, "hold-1")
  file.create(hold)
  withr::defer(unlink(hold))
  run_launcher(root, launcher, "--force", "-k", "2")
  expect_identical(exit_line(launcher, 2L), "EXIT 0")

  # slice 1 waits: it runs, and a second launch is refused
  status <- run_launcher(root, launcher, "--status")
  expect_true(any(grepl("^slice-1.log: running \\(pid [0-9]+\\)$", status)))
  again <- run_launcher(root, launcher, "--force", "-k", "2")
  expect_identical(attr(again, "status"), 2L)
  expect_true(any(grepl("still run", again)))
  rerun <- run_launcher(root, launcher, "--rerun", "1")
  expect_identical(attr(rerun, "status"), 2L)
  expect_true(any(grepl(
    sprintf("slice 1 still runs (pid %d)", pid_of_slice(launcher, 1L)),
    rerun,
    fixed = TRUE
  )))

  # a TERM to the recorded pid stops the R process too, and the log says so
  tools::pskill(pid_of_slice(launcher, 1L), tools::SIGTERM)
  expect_match(exit_line(launcher, 1L), "^EXIT [1-9][0-9]*$")

  unlink(hold)
  out <- run_launcher(root, launcher, "--rerun", "1")
  expect_null(attr(out, "status"))
  expect_identical(exit_line(launcher, 1L), "EXIT 0")
  status <- run_launcher(root, launcher, "--status")
  expect_true("cells with a result: 2 of 6" %in% status)
  expect_true("slice-1.log: EXIT 0" %in% status)
})
