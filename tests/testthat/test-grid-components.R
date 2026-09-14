# Tests for recovery_grid() with a list of components, written against
# local/dev/spec-milestone-5-correlation-recovery.md section 5.5 as amended
# for 5.5b. Fits come from the grid's mock fitter (helper-generate.R);
# nothing compiles Stan.

# The grid saves each set with its components, generators included. A
# generator from a helper file encloses the load_all() package environment,
# which saveRDS() warns it may not find when loading; this one encloses the
# base environment only.
cg_generator <- local(
  function(pars, n_trials, model) data.frame(y = rep(0, n_trials)),
  envir = new.env(parent = baseenv())
)

cg_a <- function(name = "a", ...) {
  recovery_component(
    bmm::mixture2p(resp_error = "y"), c(kappa = 2, thetat = 0),
    n_trials = 2, sds = c(kappa = 0.3, thetat = 0.2),
    generator = cg_generator, name = name, ...
  )
}

cg_b <- function(tasks = c("1", "2"), ...) {
  recovery_component(
    bmm::sdm(resp_error = "y"), c(c = 1, kappa = 1), n_trials = 2,
    sds = c(c_task1 = 0.5), generator = cg_generator, name = "b",
    tasks = tasks, ...
  )
}

cg_grid <- function() {
  data.frame(
    n_subjects = c(5L, 6L), n_trials_a = c(2L, 3L),
    cor_b_c_task1__a_kappa = c(0.2, 0.5)
  )
}

cg_run <- function(dir, mock, grid = cg_grid(), reps = 2L,
                   model = list(cg_a(), cg_b()),
                   covariates = list(G = c(mean = 0, sd = 1)), ...) {
  recovery_grid(
    model, grid,
    dir = dir, reps = reps, seed = 11, covariates = covariates,
    chains = 2, iter = 400, ...,
    .fitter = mock$fitter
  )
}

# everything but the timings, which differ between two runs
cg_parts <- function(out) {
  cells <- attr(out, "cells")
  cells$elapsed <- NULL
  list(
    rows = unclass(out)[names(out)],
    cells = cells,
    correlations = attr(out, "correlations"),
    subject_means = attr(out, "subject_means")
  )
}

cg_fit_files <- function(dir) {
  list.files(dir, pattern = "^cell-[0-9]+-rep-[0-9]+-(a|b)\\.rds$",
    full.names = TRUE
  )
}

cg_sim <- function(dir, row, rep = 1L) {
  readRDS(file.path(dir, sprintf("cell-%d-rep-%d-sim.rds", row, rep)))
}

test_that("a component grid scores prefixed terms per cell and component", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()

  out <- suppressMessages(cg_run(
    dir, mock,
    correlations = c("draws", "point")
  ))

  expect_s3_class(out, "bmmtools_recovery")
  expect_true(all(grepl("^(a|b)_", out$term)))
  expect_true(all(c("a_kappa", "b_c_task1", "b_kappa_task2") %in% out$term))
  expect_setequal(out$condition, c("row-1", "row-2"))
  expect_setequal(out$replication, 1:2)
  # a preflight per component, then 4 cells of 2 components
  expect_identical(mock$calls$n, 10L)

  cells <- attr(out, "cells")
  expect_named(cells, c(
    "condition", "replication", "component", "n_subjects", "n_trials",
    "seed", "file", "status", "elapsed", "converged"
  ))
  expect_equal(nrow(cells), 8L)
  expect_equal(cells$component, rep(c("a", "b"), 4))
  expect_equal(cells$n_trials, c(2L, 2L, 3L, 2L, 2L, 2L, 3L, 2L))
  expect_true(all(cells$status == "ok"))
  expect_true(all(file.exists(cells$file)))
  expect_true(all(file.exists(file.path(
    dir, c("cell-1-rep-1-a-est.rds", "cell-1-rep-1-b-est.rds",
           "cell-1-rep-1-cor.rds", "cell-1-rep-1-sim.rds")
  ))))
  expect_equal(attr(out, "grid"), cg_grid())

  cors <- attr(out, "correlations")
  expect_s3_class(cors, "bmmtools_cor_recovery")
  expect_true(all(c("condition", "replication") %in% names(cors)))
  expect_true("a_kappa__b_c_task1" %in% cors$term)
  expect_true("G__b_c_task1" %in% cors$term)
  truth_of <- function(row) {
    rows <- cors[cors$term == "a_kappa__b_c_task1" &
                   cors$condition == row, ]
    unique(rows$true_value)
  }
  expect_equal(truth_of("row-1"), 0.2)
  expect_equal(truth_of("row-2"), 0.5)

  means <- attr(out, "subject_means")
  expect_equal(attr(means, "links"), cg_sim(dir, 1L)$links)
  table <- subject_table(out)
  expect_true(all(c(
    "true_a_kappa", "est_a_kappa", "true_b_c_task1", "est_b_c_task1", "G"
  ) %in% names(table)))
  expect_equal(nrow(table), (5 + 6) * 2)

  sidecar <- readRDS(file.path(dir, "cell-1-rep-1-cor.rds"))
  expect_named(sidecar, c(
    "keys", "bmmtools_version", "correlations", "cor_scale",
    "cor_estimates", "subject_means"
  ))
  expect_named(sidecar$keys, c("a", "b"))
  expect_named(readRDS(file.path(dir, "cell-1-rep-1-b-est.rds")), c(
    "key", "bmmtools_version", "levels", "correlations", "cor_scale",
    "estimates", "cor_estimates"
  ))
})

# resume --------------------------------------------------------------------

# counts the calls of fit_cached(), which reads (or fits) a cell's fit
local_fit_counter <- function(env = parent.frame()) {
  counter <- new.env(parent = emptyenv())
  counter$n <- 0L
  original <- fit_cached
  local_mocked_bindings(
    fit_cached = function(...) {
      counter$n <- counter$n + 1L
      original(...)
    },
    .env = env
  )
  counter
}

test_that("a resume with the fits deleted reads sidecars only", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  first <- suppressMessages(cg_run(
    dir, grid_mock_fitter(),
    correlations = c("draws", "point")
  ))
  unlink(cg_fit_files(dir))
  expect_length(cg_fit_files(dir), 0L)

  mock <- grid_mock_fitter()
  counter <- local_fit_counter()
  expect_no_message(
    second <- cg_run(dir, mock, correlations = c("draws", "point")),
    message = "Preflight"
  )

  expect_identical(mock$calls$n, 0L)
  expect_identical(counter$n, 0L)
  expect_length(cg_fit_files(dir), 0L)
  expect_identical(cg_parts(second), cg_parts(first))
})

test_that("a different bmmtools version re-extracts through the cached fits", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  suppressMessages(cg_run(dir, grid_mock_fitter(), correlations = "draws"))

  local_mocked_bindings(bmmtools_version = function() "99.0.0")
  mock <- grid_mock_fitter()
  out <- suppressMessages(cg_run(dir, mock, correlations = "draws"))

  # every fit is cached, so the preflight is skipped and nothing is refitted
  expect_identical(mock$calls$n, 0L)
  expect_equal(
    readRDS(file.path(dir, "cell-2-rep-2-b-est.rds"))$bmmtools_version,
    "99.0.0"
  )
  expect_equal(
    readRDS(file.path(dir, "cell-2-rep-2-cor.rds"))$bmmtools_version,
    "99.0.0"
  )
  expect_s3_class(attr(out, "correlations"), "bmmtools_cor_recovery")
})

test_that("stale sidecars without fits run the preflight and refit", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  suppressMessages(cg_run(dir, grid_mock_fitter(), reps = 1L))
  unlink(cg_fit_files(dir))

  local_mocked_bindings(bmmtools_version = function() "99.0.0")
  mock <- grid_mock_fitter()
  messages <- testthat::capture_messages(
    suppressWarnings(cg_run(dir, mock, reps = 1L))
  )
  expect_true(any(grepl("Preflight passed", messages)))
  expect_identical(mock$calls$n, 2L + 4L)
})

test_that("deleting only the set sidecar reads the fits again, fitting none", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  first <- suppressMessages(cg_run(
    dir, grid_mock_fitter(),
    correlations = "point"
  ))
  unlink(list.files(dir, pattern = "-cor\\.rds$", full.names = TRUE))

  mock <- grid_mock_fitter()
  counter <- local_fit_counter()
  second <- suppressMessages(cg_run(dir, mock, correlations = "point"))

  expect_identical(mock$calls$n, 0L)
  # both fits of each of the four cells, through the cache
  expect_identical(counter$n, 8L)
  expect_length(list.files(dir, pattern = "-cor\\.rds$"), 4L)
  expect_identical(cg_parts(second), cg_parts(first))
})

test_that("a fit that cannot be obtained for the set marks its component", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  suppressMessages(cg_run(dir, grid_mock_fitter(), reps = 1L))
  unlink(cg_fit_files(dir))
  unlink(file.path(dir, "cell-1-rep-1-cor.rds"))

  mock <- grid_mock_fitter(fail_on = function(n, data, dots) {
    "task" %in% names(data)
  })
  expect_warning(
    out <- suppressMessages(cg_run(dir, mock, reps = 1L)),
    "1 component failed to fit"
  )
  cells <- attr(out, "cells")
  expect_equal(cells$status, c("ok", "error", "ok", "ok"))
  means <- attr(out, "subject_means")
  expect_equal(unique(means$condition), "row-2")
})

test_that("a failing set extraction marks every component of the cell", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  local_mocked_bindings(extract_component_set = function(...) {
    stop("set extraction failed")
  })
  expect_error(
    suppressWarnings(suppressMessages(cg_run(
      dir, grid_mock_fitter(),
      reps = 1L, preflight = FALSE
    ))),
    "No cell produced a fit"
  )
  expect_false(file.exists(file.path(dir, "cell-1-rep-1-cor.rds")))
  expect_true(file.exists(file.path(dir, "cell-1-rep-1-a-est.rds")))
})

# grid columns --------------------------------------------------------------

test_that("component grid columns change each row's set", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  grid <- data.frame(
    n_subjects = 4L, n_trials_b = c(2L, 5L), a_kappa = c(NA, 3),
    b_c = c(NA, 2), b_c_task2 = c(NA, 4), sd_a_thetat = c(NA, 0.7),
    sd_b_c_task2 = c(NA, 0.6), cor_a_kappa__b_c_task1 = c(NA, 0.3),
    cor_G__a_thetat = c(NA, -0.4), ntasks = 1:2, n_tasks = 1:2
  )
  suppressMessages(cg_run(
    dir, grid_mock_fitter(),
    grid = grid, reps = 1L, preflight = FALSE
  ))
  one <- cg_sim(dir, 1L)
  two <- cg_sim(dir, 2L)

  expect_equal(one$components$b$n_trials, 2L)
  expect_equal(two$components$b$n_trials, 5L)
  expect_equal(one$pars[["a_kappa"]], 2)
  expect_equal(two$pars[["a_kappa"]], 3)
  # a bare column sets every task, a full term overrides it
  expect_equal(two$pars[c("b_c_task1", "b_c_task2")],
    c(b_c_task1 = 2, b_c_task2 = 4)
  )
  expect_equal(one$pars[["b_c_task2"]], 1)
  expect_equal(two$sds[["a_thetat"]], 0.7)
  expect_equal(two$sds[["b_c_task2"]], 0.6)
  cor_of <- function(sim, term) {
    sim$truth$cor$true_value[sim$truth$cor$term == term]
  }
  expect_equal(cor_of(one, "a_kappa__b_c_task1"), 0)
  expect_equal(cor_of(two, "a_kappa__b_c_task1"), 0.3)
  expect_equal(cor_of(one, "G__a_thetat"), 0)
  expect_equal(cor_of(two, "G__a_thetat"), -0.4)
})

test_that("columns apply to function-valued pars and a function of the row", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  a <- recovery_component(
    bmm::mixture2p(resp_error = "y"),
    pars = function() c(kappa = stats::rnorm(1, 2, 0.1), thetat = 0),
    n_trials = 2, sds = function() c(kappa = 0.3),
    generator = cg_generator, name = "a"
  )
  cors <- function(row) {
    m <- diag(2)
    dimnames(m) <- rep(list(c("a_kappa", "b_c_task1")), 2)
    m[1, 2] <- m[2, 1] <- row$rho
    m
  }
  grid <- data.frame(
    n_subjects = 4L, rho = c(0.1, 0.6), a_thetat = c(1, 2),
    sd_a_kappa = c(0.2, 0.5)
  )
  suppressMessages(cg_run(
    dir, grid_mock_fitter(),
    grid = grid, model = list(a, cg_b()), reps = 1L, cors = cors,
    covariates = NULL, preflight = FALSE
  ))
  one <- cg_sim(dir, 1L)
  two <- cg_sim(dir, 2L)
  expect_equal(two$pars[["a_thetat"]], 2)
  expect_equal(one$pars[["a_thetat"]], 1)
  expect_equal(two$sds[["a_kappa"]], 0.5)
  expect_equal(
    two$truth$cor$true_value[two$truth$cor$term == "a_kappa__b_c_task1"],
    0.6
  )
})

test_that("model may be a function of the row with a varying number of tasks", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  model <- function(row) {
    list(cg_a(), cg_b(tasks = as.character(seq_len(row$ntasks))))
  }
  grid <- data.frame(
    n_subjects = 4L, ntasks = c(2L, 3L), b_c_task3 = c(NA, 2.5)
  )
  out <- suppressMessages(cg_run(
    dir, grid_mock_fitter(),
    grid = grid, model = model, reps = 1L, covariates = NULL
  ))
  terms_of <- function(row) unique(out$term[out$condition == row])
  expect_false("b_c_task3" %in% terms_of("row-1"))
  expect_true("b_c_task3" %in% terms_of("row-2"))
  expect_equal(cg_sim(dir, 2L)$pars[["b_c_task3"]], 2.5)
})

test_that("a component grid refuses what does not fit the component form", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  run <- function(...) cg_run(dir, mock, reps = 1L, preflight = FALSE, ...)
  with_column <- function(...) {
    grid <- cg_grid()
    extra <- list(...)
    for (nm in names(extra)) grid[[nm]] <- extra[[nm]]
    grid
  }

  expect_error(run(grid = with_column(x_kappa = 1)), "x_kappa")
  expect_error(run(grid = with_column(a_foo = 1)), "a_foo")
  expect_error(run(grid = with_column(sd_a_foo = 1)), "sd_a_foo")
  expect_error(run(grid = with_column(sd_x_kappa = 1)), "sd_x_kappa")
  expect_error(run(grid = with_column(n_trials_x = 1)), "n_trials_x")
  expect_error(run(grid = with_column(cor_a_kappa__x_c = 0.1)), "cor_a_kappa")
  expect_error(run(grid = with_column(n_trials = 3)), "n_trials_a")

  expect_error(run(pars = c(kappa = 1)), "pars")
  expect_error(run(sds = c(kappa = 1)), "sds")
  expect_error(run(formula = NULL), "formula")
  expect_error(run(tasks = c("1", "2")), "tasks")
  expect_error(run(task_col = "cond"), "task_col")
  expect_error(run(re_cor = "all"), "re_cor")
  expect_error(run(generator = zero_generator), "generator")
  expect_error(run(subjects = "fixed"), "not supported for components")

  mixed <- function(row) {
    if (row$n_subjects == 5L) list(cg_a(), cg_b()) else cg_a()$model
  }
  expect_error(run(model = mixed), "row 2")
  renamed <- function(row) {
    if (row$n_subjects == 5L) {
      list(cg_a(), cg_b())
    } else {
      list(cg_a(), cg_a(name = "c"))
    }
  }
  expect_error(run(model = renamed), "row 2")
  expect_error(run(model = list(cg_a(name = "sim"), cg_b())), "sim")
  expect_error(run(model = list(cg_a(name = "cor"), cg_b())), "cor")

  single_then_components <- function(row) {
    if (row$n_subjects == 5L) cg_a()$model else list(cg_a(), cg_b())
  }
  expect_error(
    recovery_grid(
      single_then_components,
      grid = data.frame(n_subjects = c(5L, 6L), n_trials = 2L),
      pars = c(kappa = 1, thetat = 0), dir = dir, preflight = FALSE,
      generator = zero_generator, .fitter = mock$fitter
    ),
    "row 2"
  )
  expect_identical(mock$calls$n, 0L)
})

# failures, estimators, preflight -------------------------------------------

test_that("a failed component fit is recorded and the grid goes on", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  # component b of row 2 (6 subjects, with a task column) fails
  mock <- grid_mock_fitter(fail_on = function(n, data, dots) {
    "task" %in% names(data) && length(levels(data$id)) == 6L
  })
  expect_warning(
    out <- suppressMessages(cg_run(
      dir, mock,
      preflight = FALSE, correlations = "draws"
    )),
    "2 components failed to fit"
  )
  cells <- attr(out, "cells")
  expect_equal(
    cells$status[cells$condition == "row-2"],
    c("ok", "error", "ok", "error")
  )
  expect_true(all(cells$status[cells$condition == "row-1"] == "ok"))
  expect_true(all(is.na(cells$converged[cells$status == "error"])))
  row_2_terms <- unique(out$term[out$condition == "row-2"])
  expect_true(length(row_2_terms) > 0L)
  expect_true(all(startsWith(row_2_terms, "a_")))
  expect_equal(unique(attr(out, "correlations")$condition), "row-1")
  expect_equal(unique(attr(out, "subject_means")$condition), "row-1")
})

test_that("the model estimator needs a component formula", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  expect_error(
    cg_run(dir, mock, correlations = "model"),
    "every component uses the default formula"
  )
  expect_identical(mock$calls$n, 0L)

  model <- bmm::mixture2p(resp_error = "y")
  a <- cg_a(formula = recovery_formula(model, re_cor = "all"))
  expect_warning(
    out <- suppressMessages(cg_run(
      dir, mock,
      model = list(a, cg_b()), reps = 1L, correlations = "model"
    )),
    "no cell has a correlation pair"
  )
  expect_s3_class(out, "bmmtools_recovery")
  sidecar <- readRDS(file.path(dir, "cell-1-rep-1-a-est.rds"))
  expect_equal(sidecar$correlations, "model")
  expect_s3_class(sidecar$cor_estimates, "tbl_df")
})

test_that("a failing preflight names its component and runs no cell", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter(fail_on = function(n, data, dots) {
    "task" %in% names(data)
  })
  expect_error(
    suppressMessages(cg_run(dir, mock)),
    "In component"
  )
  expect_length(cg_fit_files(dir), 0L)
})

test_that("smoke mode runs a component grid in its own directory", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  grid <- data.frame(n_subjects = c(4L, 5L, 6L))
  out <- suppressMessages(cg_run(
    dir, grid_mock_fitter(),
    grid = grid, reps = 3L, smoke = TRUE, preflight = FALSE
  ))
  expect_equal(nrow(attr(out, "cells")), 2L * 2L * 2L)
  expect_true(file.exists(file.path(dir, "smoke", "cell-2-rep-2-cor.rds")))
})

test_that("rows must keep one link table and sd_ columns add missing SDs", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  relinked <- function(row) {
    if (row$n_subjects == 5L) {
      list(cg_a(), cg_b())
    } else {
      other <- recovery_component(
        bmm::sdm(resp_error = "y"), c(c = 1, kappa = 1), n_trials = 2,
        generator = cg_generator, name = "a"
      )
      list(other, cg_b())
    }
  }
  expect_error(
    cg_run(dir, mock, model = relinked, reps = 1L, preflight = FALSE),
    "different links"
  )

  a <- recovery_component(
    bmm::mixture2p(resp_error = "y"), c(kappa = 2, thetat = 0),
    n_trials = 2, generator = cg_generator, name = "a"
  )
  grid <- data.frame(n_subjects = 4L, sd_a_kappa = 0.4)
  suppressMessages(cg_run(
    dir, mock,
    grid = grid, model = list(a, cg_b()), reps = 1L, covariates = NULL,
    preflight = FALSE
  ))
  expect_equal(cg_sim(dir, 1L)$sds[["a_kappa"]], 0.4)
})
