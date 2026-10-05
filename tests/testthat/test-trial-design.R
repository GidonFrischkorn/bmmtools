# Tests for `trial_design`, the per-trial design input of D52-D55, written
# against local/dev/spec-milestone-21-design-and-adapters.md section 6
# (rows B1-B16). Most run on a model stub without bmm: the engine only
# reads `links` and `resp_vars`, and a user generator does the rest.

td_model <- function() {
  structure(
    list(
      name = "toy",
      links = list(a = "identity", b = "identity"),
      resp_vars = list(y = "y")
    ),
    class = "bmmodel"
  )
}

#' A generator that records the design it was handed and uses it
td_recorder <- function() {
  seen <- new.env(parent = emptyenv())
  seen$designs <- list()
  generator <- function(pars, n_trials, model, trial_design) {
    seen$designs[[length(seen$designs) + 1L]] <- trial_design
    data.frame(y = pars$a + trial_design$ss + stats::runif(n_trials))
  }
  list(generator = generator, seen = seen)
}

td_sim <- function(..., n_subjects = 3, n_trials = 4, seed = 21) {
  simulate_recovery(
    td_model(), c(a = 0.5, b = 2),
    n_subjects = n_subjects, n_trials = n_trials,
    sds = c(a = 0.25), seed = seed, ...
  )
}

# B1: no design is the 0.2.0 engine -------------------------------------

td_old_generator <- function(pars, n_trials, model) {
  data.frame(y = pars$a + pars$b * stats::runif(n_trials))
}

td_old_run <- function(...) {
  simulate_recovery(
    td_model(), c(a = 0.5, b = 2),
    n_subjects = 4, n_trials = 6,
    sds = c(a = 0.25), generator = td_old_generator, seed = 21, ...
  )
}

test_that("without a design the data are those of 0.2.0", {
  # Saved on macOS arm64 from these calls, after checking that they give
  # the hashes recorded on b3cfb06, before trial_design existed (next
  # test). Compared with a tolerance, not bit for bit: the subject values
  # come from rnorm(), whose qnorm() arithmetic differs in the last ulp on
  # the x86_64 Linux runner (CI on PR #25: same R 4.6.1, other hashes).
  # A changed random number stream still fails this loudly.
  old <- readRDS(test_path("fixtures", "simulation-b3cfb06.rds"))
  expect_equal(td_old_run()$data, old$plain)
  expect_equal(
    td_old_run(
      tasks = c("1", "2"), covariates = list(g = c(mean = 0, sd = 1))
    )$data,
    old$tasks_covariates
  )
})

test_that("without a design the data are byte-identical to 0.2.0", {
  # Hashes recorded 2026-10-04 on b3cfb06 (local/dev/STATE-milestone-21.md).
  # cache_key() hashes `data` as one component, so equal data is an equal
  # key, but only on the platform the key was written on: the hashes hold
  # on macOS arm64 and not on the Linux runner (previous test).
  skip_if_not(
    Sys.info()[["sysname"]] == "Darwin" && R.version$arch == "aarch64",
    "hashes recorded on macOS arm64; rnorm() differs by an ulp elsewhere"
  )
  expect_identical(
    rlang::hash(td_old_run()$data), "f9c487efbb40577b2095c232743d4dc2"
  )
  expect_identical(
    rlang::hash(td_old_run(
      tasks = c("1", "2"), covariates = list(g = c(mean = 0, sd = 1))
    )$data),
    "9838bcdfe9f7c233a55f704b9219a686"
  )
})

test_that("a 0.2.0 generator is called with three arguments", {
  called_with <- NULL
  old <- function(pars, n_trials, model) {
    called_with <<- nargs()
    data.frame(y = stats::runif(n_trials))
  }
  sim <- td_sim(generator = old)
  expect_identical(called_with, 3L)
  expect_null(sim$trial_design)
  expect_identical(sim$trial_design_columns, character())
})

# B2-B4: the design reaches the generator and the data unchanged --------

test_that("a shared data frame reaches every generator call and the data", {
  rec <- td_recorder()
  design <- data.frame(ss = c(1, 2, 3, 4), cue = c("a", "b", "a", "b"))
  rownames(design) <- c("w", "x", "y", "z")
  sim <- td_sim(trial_design = design, generator = rec$generator)

  expect_length(rec$seen$designs, 3L)
  expected <- design
  rownames(expected) <- NULL
  for (d in rec$seen$designs) expect_identical(d, expected)

  expect_named(sim$data, c("id", "ss", "cue", "y"))
  for (i in 1:3) {
    rows <- as.data.frame(sim$data[sim$data$id == i, c("ss", "cue")])
    expect_identical(rows, expected)
  }
  expect_identical(sim$trial_design, expected)
  expect_identical(sim$trial_design_columns, c("ss", "cue"))
})

test_that("a data frame with id gives each subject its own rows", {
  rec <- td_recorder()
  design <- data.frame(
    id = rep(c(2, 1, 3), each = 4),
    ss = c(21:24, 11:14, 31:34)
  )
  sim <- td_sim(trial_design = design, generator = rec$generator)

  expect_identical(rec$seen$designs[[1L]], data.frame(ss = 11:14))
  expect_identical(rec$seen$designs[[2L]], data.frame(ss = 21:24))
  expect_identical(rec$seen$designs[[3L]], data.frame(ss = 31:34))
  expect_identical(sim$data$ss, c(11:14, 21:24, 31:34))
  expect_named(sim$data, c("id", "ss", "y"))
})

test_that("a function is called once per subject and task under the seed", {
  calls <- 0L
  design <- function(n_trials) {
    calls <<- calls + 1L
    data.frame(ss = stats::runif(n_trials))
  }
  rec <- td_recorder()
  sim <- td_sim(
    trial_design = design, generator = rec$generator, tasks = c("1", "2")
  )
  expect_identical(calls, 6L)
  expect_named(sim$data, c("id", "task", "ss", "y"))
  # fresh draws per call, reproducible under the seed
  expect_false(identical(rec$seen$designs[[1L]], rec$seen$designs[[2L]]))
  again <- td_sim(
    trial_design = design, generator = rec$generator, tasks = c("1", "2")
  )
  expect_identical(again$data, sim$data)
  expect_type(sim$trial_design, "character")
  expect_identical(sim$trial_design_columns, "ss")
})

test_that("the design columns follow id, the covariates and the task", {
  rec <- td_recorder()
  sim <- td_sim(
    trial_design = data.frame(ss = 1:4), generator = rec$generator,
    tasks = c("1", "2"), covariates = list(g = c(mean = 0, sd = 1))
  )
  expect_named(sim$data, c("id", "g", "task", "ss", "y"))
})

test_that("a generator with ... receives the design by name", {
  got <- NULL
  generator <- function(pars, n_trials, model, ...) {
    got <<- list(...)$trial_design
    data.frame(y = stats::runif(n_trials))
  }
  td_sim(trial_design = data.frame(ss = 1:4), generator = generator)
  expect_identical(got, data.frame(ss = 1:4))
})

# B5-B6: a design that conflicts with n_trials --------------------------

test_that("a design whose rows differ from n_trials errors naming both", {
  rec <- td_recorder()
  expect_error(
    td_sim(trial_design = data.frame(ss = 1:3), generator = rec$generator),
    "3 rows.*`n_trials` is 4"
  )
  per_id <- data.frame(id = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3), ss = 0)
  expect_error(
    td_sim(trial_design = per_id, generator = rec$generator),
    "subject 2.*3 rows.*`n_trials` is 4"
  )
  short <- function(n_trials) data.frame(ss = seq_len(n_trials - 1L))
  expect_error(
    td_sim(trial_design = short, generator = rec$generator),
    "3 rows.*`n_trials` is 4"
  )
  expect_length(rec$seen$designs, 0L)
})

test_that("an id-keyed design must name exactly the subjects", {
  rec <- td_recorder()
  design <- data.frame(id = rep(c(1, 2, 5), each = 4), ss = 0)
  expect_error(
    td_sim(trial_design = design, generator = rec$generator),
    "1 to 3"
  )
})

test_that("a grid whose n_trials disagrees with a design stops before any cell", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  generator <- function(pars, n_trials, model, trial_design) {
    data.frame(y = stats::runif(n_trials, -1, 1))
  }
  expect_error(
    recovery_grid(
      bmm::mixture2p(resp_error = "y"),
      grid = small_grid(),
      pars = c(kappa = log(8), thetat = stats::qlogis(0.75)),
      dir = dir, trial_design = data.frame(ss = 1:10),
      generator = generator, seed = 1, .fitter = mock$fitter
    ),
    "row 2.*10 rows.*`n_trials` is 20"
  )
  expect_identical(mock$calls$n, 0L)
})

# B7-B11: what is refused -----------------------------------------------

test_that("a generator without a trial_design formal is refused with a design", {
  old <- function(pars, n_trials, model) data.frame(y = 0)
  expect_error(
    td_sim(trial_design = data.frame(ss = 1:4), generator = old),
    "fourth argument"
  )
})

test_that("an adapter that reads no design refuses one", {
  skip_if_not_installed("bmm")
  expect_error(
    simulate_recovery(
      bmm::mixture2p(resp_error = "y"),
      c(kappa = log(8), thetat = 1),
      n_subjects = 2, n_trials = 3,
      trial_design = data.frame(ss = 1:3)
    ),
    "mixture2p.*does not use a trial design"
  )
})

test_that("trial_design must be NULL, a data frame or a function", {
  rec <- td_recorder()
  expect_error(
    td_sim(trial_design = 1:4, generator = rec$generator),
    "must be `NULL`, a data frame or a function"
  )
  expect_error(
    td_sim(
      trial_design = data.frame(id = rep(1:3, each = 4)),
      generator = rec$generator
    ),
    "at least one column"
  )
  with_id <- function(n_trials) data.frame(id = 1, ss = seq_len(n_trials))
  expect_error(
    td_sim(trial_design = with_id, generator = rec$generator),
    "returned a column \"id\""
  )
  not_df <- function(n_trials) seq_len(n_trials)
  expect_error(
    td_sim(trial_design = not_df, generator = rec$generator),
    "must return a data frame"
  )
})

test_that("a design column may not clash with the data's own columns", {
  rec <- td_recorder()
  expect_error(
    td_sim(
      trial_design = data.frame(ss = 1:4, y = 0),
      generator = rec$generator
    ),
    "\"y\""
  )
  expect_error(
    td_sim(
      trial_design = data.frame(ss = 1:4, task = 0),
      generator = rec$generator, tasks = c("1", "2")
    ),
    "\"task\""
  )
  expect_error(
    td_sim(
      trial_design = data.frame(ss = 1:4, g = 0), generator = rec$generator,
      covariates = list(g = c(mean = 0, sd = 1))
    ),
    "\"g\""
  )
})

test_that("a generator may not return a design column", {
  generator <- function(pars, n_trials, model, trial_design) {
    data.frame(y = 0, ss = trial_design$ss)
  }
  expect_error(
    td_sim(trial_design = data.frame(ss = 1:4), generator = generator),
    "\"ss\".*trial design"
  )
})

test_that("a generator must return one row per design row", {
  generator <- function(pars, n_trials, model, trial_design) {
    data.frame(y = c(0, 0))
  }
  expect_error(
    td_sim(trial_design = data.frame(ss = 1:4), generator = generator),
    "2 rows.*4"
  )
})

# B12: the key follows the data ------------------------------------------

test_that("designs that differ in one value give different keys", {
  rec <- td_recorder()
  a <- td_sim(trial_design = data.frame(ss = 1:4), generator = rec$generator)
  b <- td_sim(
    trial_design = data.frame(ss = c(1:3, 5)), generator = rec$generator
  )
  key <- function(sim) {
    cache_key(fake_bmmformula(y ~ 1), sim$data, td_model(), NULL, list())
  }
  expect_false(identical(key(a)$key, key(b)$key))
  differs <- names(key(a)$components)[
    key(a)$components != key(b)$components
  ]
  expect_identical(differs, "data")
})

# B13: the grid record ---------------------------------------------------

test_that("a grid without a design writes the record it wrote before", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  run <- function(...) {
    recovery_grid(
      bmm::mixture2p(resp_error = "y"),
      grid = small_grid(), reps = 1L,
      pars = c(kappa = log(8), thetat = stats::qlogis(0.75)),
      dir = dir, seed = 1, .fitter = mock$fitter, ...
    )
  }
  suppressMessages(run())
  record <- readRDS(file.path(dir, "grid.rds"))
  expect_false("trial_design" %in% names(record))
  # a record written before 21.1 has no field either, and is not rewritten
  expect_no_message(run(), message = "Rewriting")
  expect_identical(mock$calls$n, 3L)
})

test_that("a grid passes the design to every cell and records it", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  generator <- function(pars, n_trials, model, trial_design) {
    data.frame(y = stats::runif(n_trials, -1, 1))
  }
  design <- function(n_trials) data.frame(ss = sample(2:6, n_trials, TRUE))
  run <- function(trial_design) {
    recovery_grid(
      bmm::mixture2p(resp_error = "y"),
      grid = small_grid(), reps = 1L,
      pars = c(kappa = log(8), thetat = stats::qlogis(0.75)),
      dir = dir, seed = 1, trial_design = trial_design,
      generator = generator, .fitter = mock$fitter
    )
  }
  suppressMessages(run(design))
  for (row in 1:2) {
    sim <- readRDS(file.path(dir, sprintf("cell-%d-rep-1-sim.rds", row)))
    expect_named(sim$data, c("id", "ss", "y"))
    expect_identical(sim$trial_design_columns, "ss")
  }
  record <- readRDS(file.path(dir, "grid.rds"))
  expect_identical(record$trial_design, function_key(design))
  # the same function built again is the same record
  design2 <- function(n_trials) data.frame(ss = sample(2:6, n_trials, TRUE))
  expect_no_message(run(design2), message = "Rewriting")
})

# B15-B16: print and stored simulations ---------------------------------

test_that("print names the design columns", {
  rec <- td_recorder()
  sim <- td_sim(trial_design = data.frame(ss = 1:4), generator = rec$generator)
  expect_output(print(sim), "trial design: ss")
})

test_that("a simulation stored before trial_design is upgraded", {
  rec <- td_recorder()
  sim <- unclass(td_sim(generator = function(pars, n_trials, model) {
    data.frame(y = 0)
  }))
  sim$trial_design <- NULL
  sim$trial_design_columns <- NULL
  old <- structure(sim, class = "bmmtools_simulation")
  up <- upgrade_simulation(old)
  expect_true("trial_design" %in% names(up))
  expect_null(up$trial_design)
  expect_identical(up$trial_design_columns, character())
})

# components --------------------------------------------------------------

test_that("a component passes its trial design to its generator", {
  rec <- td_recorder()
  comp <- recovery_component(
    td_model(), c(a = 0.5, b = 2),
    n_trials = 4,
    sds = c(a = 0.25), trial_design = data.frame(ss = 1:4),
    generator = rec$generator, name = "m"
  )
  set <- simulate_components(list(comp), n_subjects = 2, seed = 3)
  expect_named(set$components$m$data, c("id", "ss", "y"))
  expect_length(rec$seen$designs, 2L)
})

test_that("a component refuses a design that disagrees with its n_trials", {
  rec <- td_recorder()
  expect_error(
    recovery_component(
      td_model(), c(a = 0.5, b = 2),
      n_trials = 4,
      trial_design = data.frame(ss = 1:3),
      generator = rec$generator, name = "m"
    ),
    "3 rows.*`n_trials` is 4"
  )
})

# trial_design_columns() ---------------------------------------------------

test_that("no built-in adapter reads a design yet", {
  skip_if_not_installed("bmm")
  for (m in list(
    bmm::mixture2p(resp_error = "y"), bmm::sdm(resp_error = "y"),
    bmm::ddm(rt = "rt", response = "r")
  )) {
    expect_identical(trial_design_columns(m), character())
  }
})

# an adapter that reads a design (none does until 21.2) ----------------------

test_that("an adapter that reads a design needs one, with its columns", {
  skip_if_not_installed("bmm")
  design_adapter <- function(pars, n_trials, model, trial_design) {
    data.frame(y = stats::runif(n_trials, -1, 1))
  }
  local_mocked_bindings(
    generator_for = function(model) design_adapter,
    trial_design_columns = function(model) "ss"
  )
  model <- bmm::mixture2p(resp_error = "y")
  run <- function(trial_design) {
    simulate_recovery(
      model, c(kappa = log(8), thetat = 1),
      n_subjects = 2, n_trials = 3, trial_design = trial_design, seed = 1
    )
  }
  expect_error(run(NULL), "needs `trial_design`.*\"ss\"")
  expect_error(run(data.frame(cue = 1:3)), "lacks the column \"ss\"")
  expect_error(
    run(function(n_trials) data.frame(cue = seq_len(n_trials))),
    "lacks the column \"ss\""
  )
  sim <- run(data.frame(ss = 1:3))
  expect_identical(sim$generator, "adapter:mixture2p")
  expect_named(sim$data, c("id", "ss", "y"))
})

test_that("a design function must return the same columns on every call", {
  rec <- td_recorder()
  calls <- 0L
  design <- function(n_trials) {
    calls <<- calls + 1L
    out <- data.frame(ss = seq_len(n_trials))
    if (calls > 1L) out$extra <- 0
    out
  }
  expect_error(
    td_sim(trial_design = design, generator = rec$generator),
    "on one call"
  )
})

# review fixes (21.1) ---------------------------------------------------------

td_grid <- function(dir, mock, ..., generator = NULL) {
  recovery_grid(
    bmm::mixture2p(resp_error = "y"),
    grid = small_grid(), reps = 1L,
    pars = c(kappa = log(8), thetat = stats::qlogis(0.75)),
    dir = dir, seed = 1, generator = generator, .fitter = mock$fitter, ...
  )
}

td_grid_generator <- function(pars, n_trials, model, trial_design) {
  data.frame(y = stats::runif(n_trials, -1, 1))
}

test_that("dropping a design on a rerun is named in the grid record", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  design <- function(n_trials) data.frame(ss = seq_len(n_trials))
  suppressMessages(
    td_grid(dir, mock, trial_design = design, generator = td_grid_generator)
  )
  expect_message(
    td_grid(dir, mock, generator = td_grid_generator),
    "Rewriting.*trial_design"
  )
  expect_false("trial_design" %in% names(readRDS(file.path(dir, "grid.rds"))))
})

test_that("a grid refuses a generator without the design before writing", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  old <- function(pars, n_trials, model) data.frame(y = rep(0, n_trials))
  design <- function(n_trials) data.frame(ss = seq_len(n_trials))
  expect_error(
    td_grid(dir, mock, trial_design = design, generator = old),
    "fourth argument"
  )
  expect_false(file.exists(file.path(dir, "grid.rds")))
  expect_identical(mock$calls$n, 0L)

  # a resume whose first cell is already on disk still asks
  suppressMessages(td_grid(dir, mock, generator = old))
  expect_error(
    td_grid(dir, mock, trial_design = design, generator = old),
    "fourth argument"
  )
  # and so does an adapter that reads no design
  expect_error(
    td_grid(dir, mock, trial_design = design),
    "does not use a trial design"
  )
})

test_that("a component grid checks each row's n_trials before any cell", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  comp <- recovery_component(
    bmm::mixture2p(resp_error = "y"), c(kappa = 2, thetat = 0),
    n_trials = 3, trial_design = data.frame(ss = 1:3),
    generator = td_grid_generator, name = "a"
  )
  expect_error(
    recovery_grid(
      list(comp),
      data.frame(n_subjects = c(2L, 2L), n_trials_a = c(3L, 4L)),
      dir = dir, seed = 1, .fitter = mock$fitter
    ),
    "row 2, component a.*3 rows.*`n_trials` is 4"
  )
  expect_identical(mock$calls$n, 0L)
})

test_that("an id-keyed design with interleaved rows and a factor id", {
  rec <- td_recorder()
  design <- data.frame(
    id = factor(rep(c(3, 1, 2), times = 4), levels = c(3, 2, 1)),
    ss = c(9, 4, 7, 1, 8, 2, 6, 5, 3, 0, 11, 10)
  )
  sim <- td_sim(trial_design = design, generator = rec$generator)
  for (i in 1:3) {
    expect_identical(
      rec$seen$designs[[i]]$ss,
      design$ss[as.character(design$id) == as.character(i)]
    )
    expect_identical(
      sim$data$ss[sim$data$id == i],
      design$ss[as.character(design$id) == as.character(i)]
    )
  }
})

test_that("a component saved before trial_design existed still simulates", {
  comp <- recovery_component(
    td_model(), c(a = 0.5, b = 2),
    n_trials = 2,
    generator = function(pars, n_trials, model) {
      data.frame(y = rep(0, n_trials))
    },
    name = "m"
  )
  old <- unclass(comp)
  old$trial_design <- NULL
  old <- structure(old, class = "bmmtools_component")
  set <- simulate_components(list(old), n_subjects = 2, seed = 1)
  expect_named(set$components$m$data, c("id", "y"))
})
