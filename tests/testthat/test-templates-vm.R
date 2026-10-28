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

#' What `grid-slice.R info` answers: the fields of its marked line
#' @noRd
slice_info <- function(args) {
  out <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(shQuote(vm_template("grid-slice.R")), shQuote(args), "info"),
    stdout = TRUE, stderr = TRUE
  )
  line <- grep("^launch-info: ", out, value = TRUE)
  strsplit(trimws(sub("^launch-info: ", "", line)), " ")[[1L]]
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
  # a dir that does not exist stays as written
  withr::local_dir(dir)
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

test_that("grid-slice.R info expands a ~ in dir, as R does", {
  dir <- withr::local_tempdir()
  withr::local_envvar(HOME = dir)
  args <- write_args(dir)
  writeLines(sub("'fits/grid'", "'~/fits/grid'", readLines(args)), args)
  expect_identical(
    slice_info(args), c("4", "6", file.path(dir, "fits", "grid"))
  )
})

test_that("grid-slice.R info marks its line after a grid-args.R that prints", {
  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  args <- write_args(dir)
  # no newline: the info line would have run on from this one
  writeLines(c("cat('no newline here')", readLines(args)), args)
  expect_identical(slice_info(args), c("4", "6", "fits/grid"))
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
#' The stub answers `info` with the real grid-slice.R, after a
#' grid-args.R that prints without a newline. Its setup starts a child
#' with processx when `setup_child` is TRUE and writes its pid to
#' `setup-child.pid`, waits while a file `hold-setup` exists, exits 4
#' when a file `fail-setup` exists, and makes the grid folder. Slice k
#' starts a child process when `child[k]` is "shell" or "processx" and
#' writes its pid to `child-k.pid`, waits while a file `hold-k` exists,
#' exits 3 when k is in `fail`, and writes `cell-k-est.rds`. The hold
#' files are looked for in `hold_dir`, the grid folder by default.
#' @noRd
stub_launcher <- function(root, fail = integer(0), chains = 1L,
                          dir = "fits/my grid", child = character(0),
                          setup_child = FALSE, hold_dir = NULL, reps = 2L) {
  vm <- file.path(root, "my vm")
  dir.create(vm)
  file.copy(vm_template("launch.sh"), vm)
  writeLines(c(
    "args <- commandArgs(trailingOnly = TRUE)",
    sprintf(
      "if (args[[2L]] == 'info') source(%s)",
      deparse(vm_template("grid-slice.R"))
    ),
    "source(args[[1L]])",
    "dir <- path.expand(grid_args$dir)",
    "hold_dir <- grid_args$hold_dir",
    "if (is.null(hold_dir)) hold_dir <- dir",
    "hold <- function(name) {",
    "  while (file.exists(file.path(hold_dir, name))) Sys.sleep(0.1)",
    "}",
    "if (args[[2L]] == 'setup') {",
    "  if (isTRUE(grid_args$setup_child)) {",
    "    p <- processx::process$new('sleep', '300')",
    "    pid <- as.character(p$get_pid())",
    "    writeLines(pid, file.path(dir, 'setup-child.pid'))",
    "  }",
    "  hold('hold-setup')",
    "  if (file.exists(file.path(dir, 'fail-setup'))) {",
    "    quit(save = 'no', status = 4L)",
    "  }",
    "  dir.create(dir, recursive = TRUE, showWarnings = FALSE)",
    "  quit(save = 'no')",
    "}",
    "k <- as.integer(args[[2L]])",
    "pid_file <- file.path(dir, sprintf('child-%d.pid', k))",
    "kind <- grid_args$child[k]",
    "if (identical(kind, 'shell')) {",
    "  system(sprintf('sleep 300 & echo $! > %s', shQuote(pid_file)))",
    "} else if (identical(kind, 'processx')) {",
    "  p <- processx::process$new('sleep', '300')",
    "  writeLines(as.character(p$get_pid()), pid_file)",
    "}",
    "hold(paste0('hold-', k))",
    "if (k %in% grid_args$fail) quit(save = 'no', status = 3L)",
    "file.create(file.path(dir, sprintf('cell-%d-est.rds', k)))"
  ), file.path(vm, "grid-slice.R"))
  args <- file.path(vm, "grid-args.R")
  writeLines(c(
    "cat('a line from grid-args.R')",
    sprintf(
      paste(
        "grid_args <- list(grid = data.frame(x = 1:3), reps = %d,",
        "chains = %d, dir = %s, fail = %s, child = %s, setup_child = %s,",
        "hold_dir = %s)"
      ),
      as.integer(reps), chains, deparse(dir), deparse(fail), deparse(child),
      deparse(setup_child), deparse(hold_dir)
    )
  ), args)
  # path.expand() would not see a HOME the test has changed
  grid <- if (startsWith(dir, "~")) {
    sub("^~", Sys.getenv("HOME"), dir)
  } else {
    file.path(root, dir)
  }
  list(
    launch = file.path(vm, "launch.sh"), args = args,
    dir = grid, logs = file.path(grid, "logs")
  )
}

run_launcher <- function(root, launcher, ...) {
  # system2() warns on a non-zero exit, which the tests read as `status`;
  # the timeout turns a launch that waits on a held setup into a failure
  withr::with_dir(root, suppressWarnings(system2(
    "bash", c(shQuote(launcher$launch), ..., shQuote(launcher$args)),
    stdout = TRUE, stderr = TRUE, timeout = 120
  )))
}

#' Start the launcher in the background; its output goes to `out`
#' @noRd
launch_bg <- function(root, launcher, ..., out = file.path(root, "bg.out")) {
  args <- c(shQuote(launcher$launch), ..., shQuote(launcher$args))
  withr::with_dir(root, system2(
    "bash", args,
    stdout = out, stderr = paste0(out, ".err"), wait = FALSE
  ))
  out
}

#' What --status says of the setup call
#' @noRd
setup_line <- function(status) grep("^setup.log: ", status, value = TRUE)

#' The pids in a launch's lock once its setup call has started
#' @noRd
lock_pids <- function(launcher) {
  lock <- file.path(launcher$logs, "launch.lock")
  files <- file.path(lock, c("pid", "setup-pid"))
  ok <- wait_for(function() {
    all(file.exists(files)) && all(file.size(files) > 0)
  })
  if (!ok) {
    return(NULL)
  }
  c(
    launcher = as.integer(readLines(files[[1L]])),
    setup = as.integer(readLines(files[[2L]]))
  )
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

#' The launcher finds a slice's processes with ps; without a working ps
#' it counts any live process group, and it cannot see what R started
#' @noRd
skip_without_ps <- function() {
  status <- system2(
    "ps", c("-p", Sys.getpid()), stdout = FALSE, stderr = FALSE
  )
  skip_if(status != 0L, "no working ps")
}

#' Wait until `cond()` holds; FALSE when it does not within `wait` seconds
#' @noRd
wait_for <- function(cond, wait = 60) {
  deadline <- Sys.time() + wait
  while (Sys.time() < deadline) {
    if (isTRUE(cond())) {
      return(TRUE)
    }
    Sys.sleep(0.2)
  }
  FALSE
}

pid_alive <- function(pid) isTRUE(tools::pskill(pid, 0L))

#' Send a signal to a process group (the launcher records group ids)
#' @noRd
kill_group <- function(pgid, signal = "TERM") {
  cmd <- sprintf("kill -%s -- -%d 2>/dev/null", signal, pgid)
  invisible(system2("bash", c("-c", shQuote(cmd))))
}

#' An `nproc` on PATH that counts `cores`, or OMP_NUM_THREADS when set,
#' as GNU nproc does
#' @noRd
local_nproc <- function(root, cores, env = parent.frame()) {
  bin <- file.path(root, "bin")
  dir.create(bin)
  writeLines(c(
    "#!/bin/sh",
    sprintf('echo "${OMP_NUM_THREADS:-%d}"', cores)
  ), file.path(bin, "nproc"))
  Sys.chmod(file.path(bin, "nproc"), "755")
  withr::local_path(bin, action = "prefix", .local_envir = env)
}

#' A `ps` on PATH that fails, as where ps is refused
#' @noRd
local_broken_ps <- function(root, env = parent.frame()) {
  bin <- file.path(root, "bin")
  dir.create(bin)
  writeLines(c("#!/bin/sh", "exit 1"), file.path(bin, "ps"))
  Sys.chmod(file.path(bin, "ps"), "755")
  withr::local_path(bin, action = "prefix", .local_envir = env)
}

#' The pids whose parent is `pid`
#' @noRd
children_of <- function(pid) {
  rows <- system2("ps", c("-A", "-o", "pid=", "-o", "ppid="), stdout = TRUE)
  rows <- utils::read.table(text = rows)
  rows[[1L]][rows[[2L]] == pid]
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
  # the launch lock is gone once the slices are started
  expect_false(dir.exists(file.path(launcher$logs, "launch.lock")))

  status <- run_launcher(root, launcher, "--status")
  expect_identical(setup_line(status), "setup.log: EXIT 0")
  expect_true("cells with a result: 2 of 6" %in% status)
  expect_true("slice-2.log: EXIT 3" %in% status)
  expect_identical(
    grep("^slice-", status, value = TRUE),
    c("slice-1.log: EXIT 0", "slice-2.log: EXIT 3", "slice-3.log: EXIT 0")
  )
})

test_that("launch.sh sets K from the cores, not from OpenMP's settings", {
  skip_without_launcher()
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  local_nproc(root, 4L)
  withr::local_envvar(OMP_NUM_THREADS = "1")
  # 2 chains a slice: 4 cores hold 2 slices
  launcher <- stub_launcher(root, chains = 2L)
  out <- run_launcher(root, launcher)
  expect_null(attr(out, "status"))
  expect_true(any(grepl(
    "cores 4, cores per slice 2 (chains x threads): at most 2 slices",
    out,
    fixed = TRUE
  )))
  expect_true(any(grepl("started 2 slices over 6 cells", out)))
  expect_identical(exit_line(launcher, 2L), "EXIT 0")
})

test_that("launch.sh starts no more slices than the grid has cells", {
  skip_without_launcher()
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  launcher <- stub_launcher(root)
  out <- run_launcher(root, launcher, "--force", "-k", "10")
  expect_true(any(grepl("started 6 slices over 6 cells", out)))
  expect_identical(readLines(file.path(launcher$logs, "slices")), "6")
  expect_identical(exit_line(launcher, 6L), "EXIT 0")
})

test_that("launch.sh writes the logs where R writes a dir with ~", {
  skip_without_launcher()
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  withr::local_envvar(HOME = file.path(root, "home"))
  launcher <- stub_launcher(root, dir = "~/fits/g")
  expect_identical(launcher$dir, file.path(root, "home", "fits", "g"))
  out <- run_launcher(root, launcher, "--force", "-k", "2")
  expect_null(attr(out, "status"))
  expect_identical(exit_line(launcher, 2L), "EXIT 0")
  # no folder called ~ next to the launcher's working directory
  expect_false(file.exists(file.path(root, "~")))
  status <- run_launcher(root, launcher, "--status")
  expect_true("cells with a result: 2 of 6" %in% status)
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

test_that("launch.sh refuses an unknown option and a second file", {
  skip_without_launcher()
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  launcher <- stub_launcher(root)
  # a typo of --status before the file, and of --force after it
  out <- run_launcher(root, launcher, "--statsu")
  expect_identical(attr(out, "status"), 2L)
  expect_true(any(grepl("unknown option --statsu", out, fixed = TRUE)))
  out <- withr::with_dir(root, suppressWarnings(system2(
    "bash", c(shQuote(launcher$launch), shQuote(launcher$args), "--frce"),
    stdout = TRUE, stderr = TRUE
  )))
  expect_identical(attr(out, "status"), 2L)
  expect_true(any(grepl("unknown option --frce", out, fixed = TRUE)))
  out <- run_launcher(root, launcher, shQuote(launcher$args))
  expect_identical(attr(out, "status"), 2L)
  expect_true(any(grepl("give one grid-args.R file", out, fixed = TRUE)))
  # nothing ran
  expect_false(dir.exists(launcher$dir))
})

test_that("launch.sh names the cause when info gives no launch-info line", {
  skip_without_launcher()
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  launcher <- stub_launcher(root)
  writeLines("stop('grid-args.R broke')", launcher$args)
  out <- run_launcher(root, launcher)
  expect_identical(attr(out, "status"), 2L)
  expect_true(any(grepl("gave no launch-info line", out, fixed = TRUE)))
  expect_true(any(grepl("grid-args.R broke", out, fixed = TRUE)))
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

test_that("a TERM to a slice stops the processes its R started", {
  skip_without_launcher()
  skip_without_ps()
  skip_if_not_installed("processx")
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  # a child in R's process group, and one that processx starts in a
  # session of its own, as cmdstanr starts Stan's chains
  launcher <- stub_launcher(root, child = c("shell", "processx"))
  dir.create(launcher$dir, recursive = TRUE)
  holds <- file.path(launcher$dir, c("hold-1", "hold-2"))
  file.create(holds)
  withr::defer(unlink(holds))
  run_launcher(root, launcher, "--force", "-k", "2")
  pid_files <- file.path(launcher$dir, c("child-1.pid", "child-2.pid"))
  expect_true(wait_for(function() all(file.exists(pid_files))))
  Sys.sleep(0.5) # the pid is written once the file exists
  children <- vapply(pid_files, function(f) as.integer(readLines(f)), 1L)
  withr::defer(for (p in children) tools::pskill(p, tools::SIGKILL))
  expect_true(all(vapply(children, pid_alive, TRUE)))

  for (k in 1:2) tools::pskill(pid_of_slice(launcher, k), tools::SIGTERM)
  expect_match(exit_line(launcher, 1L), "^EXIT [1-9][0-9]*$")
  expect_match(exit_line(launcher, 2L), "^EXIT [1-9][0-9]*$")
  expect_true(wait_for(function() !any(vapply(children, pid_alive, TRUE)), 10))
})

test_that("a slice whose wrapper was killed with -9 still blocks a launch", {
  skip_without_launcher()
  skip_without_ps()
  skip_if_not_installed("processx")
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  # slice 1's R starts a child as cmdstanr starts a chain
  launcher <- stub_launcher(root, child = "processx")
  dir.create(launcher$dir, recursive = TRUE)
  hold <- file.path(launcher$dir, "hold-1")
  file.create(hold)
  withr::defer(unlink(hold))
  run_launcher(root, launcher, "--force", "-k", "2")
  expect_identical(exit_line(launcher, 2L), "EXIT 0")
  pid <- pid_of_slice(launcher, 1L)
  withr::defer(kill_group(pid, "KILL"))
  pid_file <- file.path(launcher$dir, "child-1.pid")
  expect_true(wait_for(function() file.exists(pid_file)))
  Sys.sleep(0.5) # the pid is written once the file exists
  child <- as.integer(readLines(pid_file))
  withr::defer(tools::pskill(child, tools::SIGKILL))

  # the wrapper goes, its R fits on
  tools::pskill(pid, tools::SIGKILL)
  expect_true(wait_for(function() !pid_alive(pid), 10))
  status <- run_launcher(root, launcher, "--status")
  line <- grep("^slice-1.log: ", status, value = TRUE)
  expect_match(line, sprintf(
    "^slice-1.log: running \\(pid %d is gone but its R, pid [0-9]+, is not;",
    pid
  ))
  again <- run_launcher(root, launcher, "--force", "-k", "2")
  expect_identical(attr(again, "status"), 2L)
  expect_true(any(grepl("still run", again)))

  # stopped as --status says, R and its child: dead without an exit
  # line, and no blocker
  listed <- sub(
    "^.*; kill ([0-9 ]+) stops R and what it started\\)$", "\\1", line
  )
  listed <- as.integer(strsplit(listed, " ")[[1L]])
  expect_true(child %in% listed)
  tools::pskill(listed, tools::SIGTERM)
  expect_true(wait_for(function() !pid_alive(child), 10))
  expect_true(wait_for(function() {
    any(grepl("slice-1.log: stopped without an exit line",
      run_launcher(root, launcher, "--status"),
      fixed = TRUE
    ))
  }, 10))
  unlink(hold)
  out <- run_launcher(root, launcher, "--force", "-k", "2")
  expect_null(attr(out, "status"))
  expect_identical(exit_line(launcher, 1L), "EXIT 0")
})

test_that("a recorded pid that another process has now does not count", {
  skip_without_launcher()
  skip_without_ps()
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  launcher <- stub_launcher(root)
  dir.create(launcher$logs, recursive = TRUE)
  # as after a reboot: slice 1's pid now leads a process group of an
  # unrelated process, as the slice's wrapper did
  other <- as.integer(system2(
    "bash", c("-c", shQuote("set -m; sleep 60 > /dev/null 2>&1 & echo $!")),
    stdout = TRUE
  ))
  withr::defer(tools::pskill(other, tools::SIGKILL))
  expect_true(pid_alive(other))
  writeLines(paste(1L, other), file.path(launcher$logs, "pids"))
  writeLines("1", file.path(launcher$logs, "slices"))
  status <- run_launcher(root, launcher, "--status")
  expect_true("slice-1.log: stopped without an exit line" %in% status)
  out <- run_launcher(root, launcher, "--force", "-k", "1")
  expect_null(attr(out, "status"))
  expect_identical(exit_line(launcher, 1L), "EXIT 0")
})

test_that("one launch at a time: a second is refused while one sets up", {
  skip_without_launcher()
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  launcher <- stub_launcher(root)
  # an earlier launch, so that there is a slice to rerun
  run_launcher(root, launcher, "--force", "-k", "2")
  expect_identical(exit_line(launcher, 2L), "EXIT 0")
  hold <- file.path(launcher$dir, "hold-setup")
  file.create(hold)
  withr::defer(unlink(hold))
  first <- launch_bg(root, launcher, "--force", "-k", "2")
  expect_false(is.null(lock_pids(launcher)))

  status <- run_launcher(root, launcher, "--status")
  expect_match(
    setup_line(status), "^setup.log: running \\(launch.sh pid [0-9]+\\)$"
  )
  second <- run_launcher(root, launcher, "--force", "-k", "3")
  expect_identical(attr(second, "status"), 2L)
  expect_true(any(grepl("another launch.sh (pid", second, fixed = TRUE)))
  # a rerun waits for the lock too
  rerun <- run_launcher(root, launcher, "--rerun", "1")
  expect_identical(attr(rerun, "status"), 2L)
  expect_true(any(grepl("another launch.sh (pid", rerun, fixed = TRUE)))

  # the first goes on and starts its slices alone
  unlink(hold)
  expect_true(wait_for(function() {
    any(grepl("^started 2 slices", readLines(first, warn = FALSE)))
  }))
  expect_identical(read.table(file.path(launcher$logs, "pids"))[[1L]], 1:2)
  expect_identical(exit_line(launcher, 2L), "EXIT 0")
})

test_that("a lock left by a killed launch is taken over", {
  skip_without_launcher()
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  launcher <- stub_launcher(root)
  lock <- file.path(launcher$logs, "launch.lock")
  dir.create(lock, recursive = TRUE)
  # the pid of a shell that has ended
  gone <- system2("bash", c("-c", shQuote("echo $$")), stdout = TRUE)
  writeLines(gone, file.path(lock, "pid"))
  status <- run_launcher(root, launcher, "--status")
  expect_identical(setup_line(status), "setup.log: none yet")
  out <- run_launcher(root, launcher, "--force", "-k", "2")
  expect_null(attr(out, "status"))
  expect_identical(exit_line(launcher, 2L), "EXIT 0")
  expect_false(dir.exists(lock))
  # the takeover lock goes with it, or the next stale lock is never taken
  expect_false(dir.exists(paste0(lock, ".takeover")))
})

test_that("a lock held by a slice's wrapper is taken over", {
  skip_without_launcher()
  skip_without_ps()
  withr::local_path(R.home("bin"), action = "prefix")
  root_a <- withr::local_tempdir()
  root_b <- withr::local_tempdir()
  a <- stub_launcher(root_a)
  b <- stub_launcher(root_b)
  dir.create(b$dir, recursive = TRUE)
  hold <- file.path(b$dir, "hold-1")
  file.create(hold)
  withr::defer(unlink(hold))
  run_launcher(root_b, b, "--force", "-k", "1")
  # bash .../launch.sh --run-slice ...: a launch.sh, but no launch
  other <- pid_of_slice(b, 1L)
  expect_true(any(grepl(
    "^slice-1.log: running", run_launcher(root_b, b, "--status")
  )))

  # as after a reboot: the pid in grid a's lock is grid b's live slice
  lock <- file.path(a$logs, "launch.lock")
  dir.create(lock, recursive = TRUE)
  writeLines(as.character(other), file.path(lock, "pid"))
  writeLines("a setup that stopped", file.path(a$logs, "setup.log"))
  status <- run_launcher(root_a, a, "--status")
  expect_identical(
    setup_line(status), "setup.log: stopped without an exit line"
  )
  out <- run_launcher(root_a, a, "--force", "-k", "2")
  expect_null(attr(out, "status"))
  expect_identical(exit_line(a, 2L), "EXIT 0")
  expect_false(dir.exists(lock))
})

test_that("a lock held by a process that is no launch.sh is taken over", {
  skip_without_launcher()
  skip_without_ps()
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  launcher <- stub_launcher(root)
  lock <- file.path(launcher$logs, "launch.lock")
  dir.create(lock, recursive = TRUE)
  # as after a reboot: the pid in the lock is a live, unrelated process
  writeLines(as.character(Sys.getpid()), file.path(lock, "pid"))
  out <- run_launcher(root, launcher, "--force", "-k", "2")
  expect_null(attr(out, "status"))
  expect_identical(exit_line(launcher, 2L), "EXIT 0")
})

test_that("of two launches that find a stale lock, one takes it over", {
  skip_without_launcher()
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  launcher <- stub_launcher(root)
  lock <- file.path(launcher$logs, "launch.lock")
  dir.create(lock, recursive = TRUE)
  gone <- system2("bash", c("-c", shQuote("echo $$")), stdout = TRUE)
  writeLines(gone, file.path(lock, "pid"))
  # the other launch holds the takeover lock: this one stands back
  dir.create(paste0(lock, ".takeover"))
  out <- run_launcher(root, launcher, "--force", "-k", "2")
  expect_identical(attr(out, "status"), 2L)
  expect_true(any(grepl("launch.lock.takeover", out, fixed = TRUE)))
  expect_false(file.exists(file.path(launcher$logs, "setup.log")))
})

test_that("a launch.sh killed with -9 in setup holds the lock while it runs", {
  skip_without_launcher()
  skip_without_ps()
  skip_if_not_installed("processx")
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  # the setup starts a child as cmdstanr starts a chain, in a session of
  # its own, which the kill that --status gives must name. The grid
  # folder does not exist yet, so `info` gives the relative dir as
  # written, and the launcher must mark the setup by the resolved one
  # that later calls see.
  launcher <- stub_launcher(root, setup_child = TRUE, hold_dir = root)
  hold <- file.path(root, "hold-setup")
  file.create(hold)
  withr::defer(unlink(hold))
  launch_bg(root, launcher, "--force", "-k", "2")
  pids <- lock_pids(launcher)
  withr::defer(tools::pskill(pids[["setup"]], tools::SIGKILL))
  pid_file <- file.path(launcher$dir, "setup-child.pid")
  expect_true(wait_for(function() file.exists(pid_file)))
  Sys.sleep(0.5) # the pid is written once the file exists
  child <- as.integer(readLines(pid_file))
  withr::defer(tools::pskill(child, tools::SIGKILL))
  expect_true(pid_alive(child))

  tools::pskill(pids[["launcher"]], tools::SIGKILL)
  expect_true(wait_for(function() !pid_alive(pids[["launcher"]]), 10))
  expect_true(pid_alive(pids[["setup"]]))
  status <- run_launcher(root, launcher, "--status")
  line <- setup_line(status)
  expect_match(line, sprintf(paste0(
    "^setup.log: running \\(pid %d\\), but its launch.sh is gone: no ",
    "slices will start; kill [0-9 ]+ stops it and what it started$"
  ), pids[["setup"]]))
  listed <- sub(
    "^.*; kill ([0-9 ]+) stops it and what it started$", "\\1", line
  )
  listed <- as.integer(strsplit(listed, " ")[[1L]])
  expect_true(pids[["setup"]] %in% listed)
  expect_true(child %in% listed)
  second <- run_launcher(root, launcher, "--force", "-k", "2")
  expect_identical(attr(second, "status"), 2L)
  expect_true(any(grepl(
    sprintf(
      "the setup call (pid %d) of a launch.sh that is gone", pids[["setup"]]
    ),
    second,
    fixed = TRUE
  )))

  # stopped as --status says: once the setup is gone, the lock is stale
  tools::pskill(listed, tools::SIGTERM)
  expect_true(wait_for(function() !pid_alive(pids[["setup"]]), 10))
  expect_true(wait_for(function() !pid_alive(child), 10))
  unlink(hold)
  status <- run_launcher(root, launcher, "--status")
  expect_identical(
    setup_line(status), "setup.log: stopped without an exit line"
  )
  # the next setup starts no child, which the test would have to stop
  args <- readLines(launcher$args)
  args <- sub("setup_child = TRUE", "setup_child = FALSE", args)
  writeLines(args, launcher$args)
  out <- run_launcher(root, launcher, "--force", "-k", "2")
  expect_null(attr(out, "status"))
  expect_identical(exit_line(launcher, 2L), "EXIT 0")
})

#' A grid folder with an U-umlaut in its name, and the C locale for the
#' launcher and what it starts: ps then prints the folder's path escaped
#' (macOS) or with ? (Linux), so the launcher cannot find its processes
#' by that path
#' @noRd
local_c_locale_dir <- function(env = parent.frame()) {
  withr::local_envvar(LC_ALL = "C", LANG = "C", .local_envir = env)
  "fits/\u00dcbung"
}

test_that("in the C locale, a setup call under a non-ASCII dir is found", {
  skip_without_launcher()
  skip_without_ps()
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  dir <- local_c_locale_dir()
  # the grid folder does not exist yet: the hold file is outside it
  launcher <- stub_launcher(root, dir = dir, hold_dir = root)
  hold <- file.path(root, "hold-setup")
  file.create(hold)
  withr::defer(unlink(hold))
  launch_bg(root, launcher, "--force", "-k", "2")
  pids <- lock_pids(launcher)
  expect_false(is.null(pids))
  withr::defer(tools::pskill(pids[["setup"]], tools::SIGKILL))

  tools::pskill(pids[["launcher"]], tools::SIGKILL)
  expect_true(wait_for(function() !pid_alive(pids[["launcher"]]), 10))
  status <- run_launcher(root, launcher, "--status")
  expect_match(setup_line(status), sprintf(
    "^setup.log: running \\(pid %d\\), but its launch.sh is gone",
    pids[["setup"]]
  ))
  second <- run_launcher(root, launcher, "--force", "-k", "2")
  expect_identical(attr(second, "status"), 2L)
  expect_true(any(grepl(
    sprintf(
      "the setup call (pid %d) of a launch.sh that is gone", pids[["setup"]]
    ),
    second,
    fixed = TRUE
  )))
  # no second setup call was started
  expect_identical(
    readLines(file.path(launcher$logs, "launch.lock", "setup-pid")),
    as.character(pids[["setup"]])
  )
})

test_that("in the C locale, a slice under a non-ASCII dir counts as running", {
  skip_without_launcher()
  skip_without_ps()
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  dir <- local_c_locale_dir()
  launcher <- stub_launcher(root, dir = dir, hold_dir = root)
  hold <- file.path(root, "hold-1")
  file.create(hold)
  withr::defer(unlink(hold))
  run_launcher(root, launcher, "--force", "-k", "2")
  expect_identical(exit_line(launcher, 2L), "EXIT 0")
  pid <- pid_of_slice(launcher, 1L)
  withr::defer(kill_group(pid, "KILL"))

  status <- run_launcher(root, launcher, "--status")
  expect_true(
    sprintf("slice-1.log: running (pid %d)", pid) %in% status
  )
  rerun <- run_launcher(root, launcher, "--rerun", "1")
  expect_identical(attr(rerun, "status"), 2L)
  expect_true(any(grepl(
    sprintf("slice 1 still runs (pid %d)", pid), rerun, fixed = TRUE
  )))
  again <- run_launcher(root, launcher, "--force", "-k", "2")
  expect_identical(attr(again, "status"), 2L)
  expect_true(any(grepl("still run", again)))
  # one copy of slice 1, the first
  expect_identical(pid_of_slice(launcher, 1L), pid)
  unlink(hold)
  expect_identical(exit_line(launcher, 1L), "EXIT 0")
})

test_that("slice 1 is not slice 10: the mark is matched whole", {
  skip_without_launcher()
  skip_without_ps()
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  launcher <- stub_launcher(root, reps = 4L, hold_dir = root)
  hold <- file.path(root, "hold-10")
  file.create(hold)
  withr::defer(unlink(hold))
  run_launcher(root, launcher, "--force", "-k", "10")
  for (k in 1:9) expect_identical(exit_line(launcher, k), "EXIT 0")
  ten <- pid_of_slice(launcher, 10L)
  withr::defer(kill_group(ten, "KILL"))

  # as after a reboot: slice 1's recorded pid is slice 10's, which runs
  cat("1", ten, "\n", file = file.path(launcher$logs, "pids"), append = TRUE)
  writeLines("a slice that stopped", file.path(launcher$logs, "slice-1.log"))
  status <- run_launcher(root, launcher, "--status")
  expect_true("slice-1.log: stopped without an exit line" %in% status)
  expect_true(sprintf("slice-10.log: running (pid %d)", ten) %in% status)
  again <- run_launcher(root, launcher, "--force", "-k", "10")
  expect_identical(attr(again, "status"), 2L)
  expect_true(any(grepl(
    sprintf("still run: slice 10 (pid %d).", ten), again, fixed = TRUE
  )))
  unlink(hold)
  expect_identical(exit_line(launcher, 10L), "EXIT 0")
})

test_that("a TERM to launch.sh in setup stops the setup and frees the lock", {
  skip_without_launcher()
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  launcher <- stub_launcher(root)
  dir.create(launcher$dir, recursive = TRUE)
  hold <- file.path(launcher$dir, "hold-setup")
  file.create(hold)
  withr::defer(unlink(hold))
  launch_bg(root, launcher, "--force", "-k", "2")
  pids <- lock_pids(launcher)
  withr::defer(tools::pskill(pids[["setup"]], tools::SIGKILL))

  tools::pskill(pids[["launcher"]], tools::SIGTERM)
  expect_true(wait_for(function() !pid_alive(pids[["setup"]]), 10))
  expect_true(wait_for(function() {
    !dir.exists(file.path(launcher$logs, "launch.lock"))
  }, 10))
  expect_false(file.exists(file.path(launcher$logs, "slices")))
})

signals <- c(TERM = tools::SIGTERM, HUP = tools::SIGHUP)
for (signal in names(signals)) {
  test_that(sprintf(
    "a %s to launch.sh in setup also stops what the setup started", signal
  ), {
    skip_without_launcher()
    skip_without_ps()
    skip_if_not_installed("processx")
    if (signal == "HUP") {
      # a HUP that is ignored when bash starts (nohup) cannot be trapped
      caught <- system2("bash", c(
        "-c", shQuote("trap 'echo caught' HUP; kill -HUP $$")
      ), stdout = TRUE)
      skip_if(!identical(caught, "caught"), "bash cannot trap HUP here")
    }
    withr::local_path(R.home("bin"), action = "prefix")
    root <- withr::local_tempdir()
    # the setup starts a child as cmdstanr starts a chain, in a session of
    # its own, which a signal to the process group misses
    launcher <- stub_launcher(root, setup_child = TRUE)
    dir.create(launcher$dir, recursive = TRUE)
    hold <- file.path(launcher$dir, "hold-setup")
    file.create(hold)
    withr::defer(unlink(hold))
    launch_bg(root, launcher, "--force", "-k", "2")
    pids <- lock_pids(launcher)
    withr::defer(tools::pskill(pids[["setup"]], tools::SIGKILL))
    pid_file <- file.path(launcher$dir, "setup-child.pid")
    expect_true(wait_for(function() file.exists(pid_file)))
    Sys.sleep(0.5) # the pid is written once the file exists
    child <- as.integer(readLines(pid_file))
    withr::defer(tools::pskill(child, tools::SIGKILL))
    expect_true(pid_alive(child))

    tools::pskill(pids[["launcher"]], signals[[signal]])
    expect_true(wait_for(function() !pid_alive(pids[["setup"]]), 10))
    expect_true(wait_for(function() !pid_alive(child), 10))
    expect_true(wait_for(function() {
      !dir.exists(file.path(launcher$logs, "launch.lock"))
    }, 10))
    expect_false(file.exists(file.path(launcher$logs, "slices")))
  })
}

test_that("a failed setup call starts no slices and frees the lock", {
  skip_without_launcher()
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  launcher <- stub_launcher(root)
  dir.create(launcher$dir, recursive = TRUE)
  file.create(file.path(launcher$dir, "fail-setup"))
  out <- run_launcher(root, launcher, "--force", "-k", "2")
  expect_identical(attr(out, "status"), 1L)
  expect_true(any(grepl("the setup call failed", out, fixed = TRUE)))
  expect_identical(
    utils::tail(readLines(file.path(launcher$logs, "setup.log")), 1L),
    "EXIT 4"
  )
  expect_false(file.exists(file.path(launcher$logs, "slices")))
  expect_false(file.exists(file.path(launcher$logs, "pids")))
  expect_false(dir.exists(file.path(launcher$logs, "launch.lock")))
  status <- run_launcher(root, launcher, "--status")
  expect_identical(setup_line(status), "setup.log: EXIT 4")
})

test_that("a dir with a backslash in it: a running slice counts as running", {
  skip_without_launcher()
  skip_without_ps()
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  # a backslash and a t, which awk -v would read as a tab
  launcher <- stub_launcher(root, dir = "fits/a\\tb")
  dir.create(launcher$dir, recursive = TRUE)
  hold <- file.path(launcher$dir, "hold-1")
  file.create(hold)
  withr::defer(unlink(hold))
  run_launcher(root, launcher, "--force", "-k", "2")
  expect_identical(exit_line(launcher, 2L), "EXIT 0")
  status <- run_launcher(root, launcher, "--status")
  expect_true(any(grepl("^slice-1.log: running \\(pid [0-9]+\\)$", status)))
  again <- run_launcher(root, launcher, "--force", "-k", "2")
  expect_identical(attr(again, "status"), 2L)
  expect_true(any(grepl("still run", again)))
  unlink(hold)
  expect_identical(exit_line(launcher, 1L), "EXIT 0")
})

test_that("a recorded pid that another grid's slice has now does not count", {
  skip_without_launcher()
  skip_without_ps()
  withr::local_path(R.home("bin"), action = "prefix")
  root_a <- withr::local_tempdir()
  root_b <- withr::local_tempdir()
  a <- stub_launcher(root_a)
  b <- stub_launcher(root_b)
  dir.create(b$dir, recursive = TRUE)
  hold <- file.path(b$dir, "hold-1")
  file.create(hold)
  withr::defer(unlink(hold))
  run_launcher(root_b, b, "--force", "-k", "1")
  other <- pid_of_slice(b, 1L)
  expect_true(any(grepl(
    "^slice-1.log: running", run_launcher(root_b, b, "--status")
  )))

  # grid a's slice 1 is recorded with grid b's live slice 1
  dir.create(a$logs, recursive = TRUE)
  writeLines(paste(1L, other), file.path(a$logs, "pids"))
  writeLines("1", file.path(a$logs, "slices"))
  status <- run_launcher(root_a, a, "--status")
  expect_true("slice-1.log: stopped without an exit line" %in% status)
})

#' A lock of grid `a` as a launch.sh killed with -9 leaves it, whose
#' setup-pid is now `setup_pid`; returns the lock's path
#' @noRd
stale_setup_lock <- function(a, setup_pid) {
  lock <- file.path(a$logs, "launch.lock")
  dir.create(lock, recursive = TRUE)
  gone <- system2("bash", c("-c", shQuote("echo $$")), stdout = TRUE)
  writeLines(gone, file.path(lock, "pid"))
  writeLines(as.character(setup_pid), file.path(lock, "setup-pid"))
  writeLines("a setup that stopped", file.path(a$logs, "setup.log"))
  lock
}

test_that("a setup-pid that another grid's slice has now does not count", {
  skip_without_launcher()
  skip_without_ps()
  skip_if_not_installed("processx")
  withr::local_path(R.home("bin"), action = "prefix")
  root_a <- withr::local_tempdir()
  root_b <- withr::local_tempdir()
  a <- stub_launcher(root_a)
  # grid b's slice R starts a child as cmdstanr starts a chain
  b <- stub_launcher(root_b, child = "processx")
  dir.create(b$dir, recursive = TRUE)
  hold <- file.path(b$dir, "hold-1")
  file.create(hold)
  withr::defer(unlink(hold))
  run_launcher(root_b, b, "--force", "-k", "1")
  pid_file <- file.path(b$dir, "child-1.pid")
  expect_true(wait_for(function() file.exists(pid_file)))
  Sys.sleep(0.5) # the pid is written once the file exists
  child <- as.integer(readLines(pid_file))
  withr::defer(tools::pskill(child, tools::SIGKILL))
  # the slice's R, which runs grid-slice.R as a setup call does
  r <- children_of(pid_of_slice(b, 1L))
  expect_length(r, 1L)

  # as after a reboot: grid a's setup-pid is grid b's slice R
  lock <- stale_setup_lock(a, r)
  status <- run_launcher(root_a, a, "--status")
  expect_identical(
    setup_line(status), "setup.log: stopped without an exit line"
  )
  # no advice to kill grid b's slice and its chain
  expect_false(any(grepl("kill", status, fixed = TRUE)))
  # while another launch takes the stale lock over, the refusal names
  # the launch.sh, not grid b's slice as a setup call
  dir.create(paste0(lock, ".takeover"))
  refused <- run_launcher(root_a, a, "--force", "-k", "2")
  expect_identical(attr(refused, "status"), 2L)
  expect_true(any(grepl("another launch.sh (pid", refused, fixed = TRUE)))
  expect_false(any(grepl("the setup call", refused, fixed = TRUE)))
  unlink(paste0(lock, ".takeover"), recursive = TRUE)
  out <- run_launcher(root_a, a, "--force", "-k", "2")
  expect_null(attr(out, "status"))
  expect_identical(exit_line(a, 2L), "EXIT 0")
  expect_false(dir.exists(lock))

  # grid b's slice was left alone and ends as it would have
  expect_true(pid_alive(r))
  expect_true(pid_alive(child))
  unlink(hold)
  expect_identical(exit_line(b, 1L), "EXIT 0")
})

test_that("a setup-pid that another grid's setup call has now does not count", {
  skip_without_launcher()
  skip_without_ps()
  withr::local_path(R.home("bin"), action = "prefix")
  root_a <- withr::local_tempdir()
  root_b <- withr::local_tempdir()
  a <- stub_launcher(root_a)
  b <- stub_launcher(root_b)
  dir.create(b$dir, recursive = TRUE)
  hold <- file.path(b$dir, "hold-setup")
  file.create(hold)
  withr::defer(unlink(hold))
  b_out <- launch_bg(root_b, b, "--force", "-k", "1")
  b_pids <- lock_pids(b)
  expect_false(is.null(b_pids))

  # as after a reboot: grid a's setup-pid is grid b's setup call
  lock <- stale_setup_lock(a, b_pids[["setup"]])
  status <- run_launcher(root_a, a, "--status")
  expect_identical(
    setup_line(status), "setup.log: stopped without an exit line"
  )
  expect_false(any(grepl("kill", status, fixed = TRUE)))
  out <- run_launcher(root_a, a, "--force", "-k", "2")
  expect_null(attr(out, "status"))
  expect_identical(exit_line(a, 2L), "EXIT 0")
  expect_false(dir.exists(lock))

  # grid b's launch goes on and starts its slice
  expect_true(pid_alive(b_pids[["setup"]]))
  unlink(hold)
  expect_true(wait_for(function() {
    any(grepl("^started 1 slices", readLines(b_out, warn = FALSE)))
  }))
  expect_identical(exit_line(b, 1L), "EXIT 0")
})

test_that("without a working ps the launcher says what it cannot check", {
  skip_without_launcher()
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  local_broken_ps(root)
  launcher <- stub_launcher(root)
  out <- run_launcher(root, launcher, "--force", "-k", "2")
  expect_null(attr(out, "status"))
  expect_identical(
    sum(grepl("ps does not work here", out, fixed = TRUE)), 1L
  )
  expect_identical(exit_line(launcher, 2L), "EXIT 0")
  status <- run_launcher(root, launcher, "--status")
  expect_true("slice-2.log: EXIT 0" %in% status)
})

test_that("without a working ps, --status claims no more of the setup's kill", {
  skip_without_launcher()
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  local_broken_ps(root)
  launcher <- stub_launcher(root)
  dir.create(launcher$dir, recursive = TRUE)
  hold <- file.path(launcher$dir, "hold-setup")
  file.create(hold)
  withr::defer(unlink(hold))
  launch_bg(root, launcher, "--force", "-k", "2")
  pids <- lock_pids(launcher)
  withr::defer(tools::pskill(pids[["setup"]], tools::SIGKILL))

  tools::pskill(pids[["launcher"]], tools::SIGKILL)
  expect_true(wait_for(function() !pid_alive(pids[["launcher"]]), 10))
  status <- run_launcher(root, launcher, "--status")
  # what the setup started in a session of its own cannot be found, and
  # outlives the kill: the line names the setup's pid and claims no more
  expect_identical(setup_line(status), sprintf(paste0(
    "setup.log: running (pid %d), but its launch.sh is gone: no slices ",
    "will start; kill %d stops it"
  ), pids[["setup"]], pids[["setup"]]))
  tools::pskill(pids[["setup"]], tools::SIGTERM)
  expect_true(wait_for(function() !pid_alive(pids[["setup"]]), 10))
})

test_that("without a working ps, a slice whose wrapper was killed with -9", {
  skip_without_launcher()
  withr::local_path(R.home("bin"), action = "prefix")
  root <- withr::local_tempdir()
  local_broken_ps(root)
  launcher <- stub_launcher(root, hold_dir = root)
  hold <- file.path(root, "hold-1")
  file.create(hold)
  withr::defer(unlink(hold))
  run_launcher(root, launcher, "--force", "-k", "2")
  expect_identical(exit_line(launcher, 2L), "EXIT 0")
  pid <- pid_of_slice(launcher, 1L)
  withr::defer(kill_group(pid, "KILL"))

  # the wrapper goes, its R fits on in the wrapper's process group
  tools::pskill(pid, tools::SIGKILL)
  expect_true(wait_for(function() !pid_alive(pid), 10))
  status <- run_launcher(root, launcher, "--status")
  # the advice is to signal the group, not the pid that is gone
  expect_true(sprintf(paste0(
    "slice-1.log: running (pid %d is gone but processes of its group are ",
    "not; kill -- -%d stops them)"
  ), pid, pid) %in% status)
  kill_group(pid, "TERM")
  expect_true(wait_for(function() {
    "slice-1.log: stopped without an exit line" %in%
      run_launcher(root, launcher, "--status")
  }, 10))
})
