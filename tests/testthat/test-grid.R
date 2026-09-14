# Tests for recovery_grid(), written against
# local/dev/spec-milestone-3-generate-layer.md section 4. The mock fitter and
# the registered extract_estimates() method for it (helper-generate.R)
# run the pipeline end to end; bmm is needed for the generate step, a
# sampler never is.

grid_run <- function(dir, mock, grid = small_grid(), reps = 2L, ...) {
  recovery_grid(
    bmm::mixture2p(resp_error = "y"),
    grid = grid,
    pars = c(kappa = log(8), thetat = stats::qlogis(0.75)),
    dir = dir,
    reps = reps,
    sds = c(kappa = 0.3, thetat = 0.4),
    seed = 100,
    chains = 2, iter = 400,
    ...,
    .fitter = mock$fitter
  )
}

test_that("a small grid runs end to end and scores both levels", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()

  out <- suppressMessages(grid_run(dir, mock))

  expect_s3_class(out, "bmmtools_recovery")
  expect_setequal(out$condition, c("row-1", "row-2"))
  expect_setequal(out$replication, c(1L, 2L))
  expect_setequal(out$level, c("population", "subject"))
  expect_true(all(out$covered))
  expect_true(all(out$converged))

  cells <- attr(out, "cells")
  expect_equal(nrow(cells), 4L)
  expect_named(cells, c(
    "condition", "replication", "n_subjects", "n_trials", "seed", "file",
    "status", "elapsed", "converged"
  ))
  expect_true(all(cells$status == "ok"))
  expect_true(all(file.exists(cells$file)))
  expect_equal(
    sort(list.files(dir, pattern = "sim\\.rds$")),
    sort(sprintf("cell-%d-rep-%d-sim.rds", c(1, 2, 1, 2), c(1, 1, 2, 2)))
  )
  # the preflight fit plus four cells
  expect_identical(mock$calls$n, 5L)
  expect_equal(attr(out, "grid"), small_grid())
})

test_that("cells run replication-major and get distinct seeds", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()

  out <- suppressMessages(grid_run(dir, mock))

  # the first call is the preflight; then rows 1, 2 at rep 1, rows 1, 2
  # at rep 2, identified by their subject counts 3, 4, 3, 4
  subjects <- vapply(mock$calls$log[-1], function(x) length(x$ids), 1L)
  expect_equal(subjects, c(3L, 4L, 3L, 4L))
  cells <- attr(out, "cells")
  expect_equal(cells$condition, c("row-1", "row-2", "row-1", "row-2"))
  expect_equal(cells$replication, c(1L, 1L, 2L, 2L))
  expect_equal(length(unique(cells$seed)), 4L)
  expect_true(is.na(cell_seed(NULL, 1L, 1L)))
  expect_equal(cell_seed(100, 1L, 1L), cell_seed(100, 1L, 1L))
  expect_false(cell_seed(100, 1L, 1L) == cell_seed(100, 2L, 1L))
})

test_that("the same seed reproduces the simulation files", {
  skip_if_not_installed("bmm")
  dir_a <- withr::local_tempdir()
  dir_b <- withr::local_tempdir()
  suppressMessages(grid_run(dir_a, grid_mock_fitter(), reps = 1L))
  suppressMessages(grid_run(dir_b, grid_mock_fitter(), reps = 1L))
  a <- readRDS(file.path(dir_a, "cell-1-rep-1-sim.rds"))
  b <- readRDS(file.path(dir_b, "cell-1-rep-1-sim.rds"))
  expect_identical(a$data, b$data)
  expect_identical(a$truth, b$truth)
})

test_that("a second call resumes without fitting", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  first <- suppressMessages(grid_run(dir, mock))
  n_first <- mock$calls$n

  second <- suppressMessages(grid_run(dir, mock))

  expect_identical(mock$calls$n, n_first)
  expect_equal(as.list(second), as.list(first), ignore_attr = TRUE)
  cells <- attr(second, "cells")
  expect_true(all(cells$status == "ok"))
})

test_that("subjects = 'fixed' reuses subject truths across replications", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  suppressMessages(grid_run(dir, grid_mock_fitter(), subjects = "fixed"))
  s11 <- readRDS(file.path(dir, "cell-1-rep-1-sim.rds"))$truth$subjects
  s12 <- readRDS(file.path(dir, "cell-1-rep-2-sim.rds"))$truth$subjects
  s21 <- readRDS(file.path(dir, "cell-2-rep-1-sim.rds"))$truth$subjects
  expect_identical(s11, s12)
  expect_false(identical(s11, s21))

  dir2 <- withr::local_tempdir()
  suppressMessages(grid_run(dir2, grid_mock_fitter()))
  r11 <- readRDS(file.path(dir2, "cell-1-rep-1-sim.rds"))$truth$subjects
  r12 <- readRDS(file.path(dir2, "cell-1-rep-2-sim.rds"))$truth$subjects
  expect_false(identical(r11, r12))
})

test_that("grid columns override the population values and the sds", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  grid <- data.frame(
    n_subjects = c(3L, 3L), n_trials = c(10L, 10L),
    kappa = c(log(4), log(16)), sd_kappa = c(0, 0.5)
  )
  out <- suppressMessages(grid_run(dir, grid_mock_fitter(),
    grid = grid,
    reps = 1L, scale = "link"
  ))
  pop <- out[out$level == "population" & out$term == "kappa", ]
  expect_equal(pop$true_value[pop$condition == "row-1"], log(4))
  expect_equal(pop$true_value[pop$condition == "row-2"], log(16))
  s1 <- readRDS(file.path(dir, "cell-1-rep-1-sim.rds"))
  s2 <- readRDS(file.path(dir, "cell-2-rep-1-sim.rds"))
  expect_false("kappa" %in% s1$truth$subjects$term)
  expect_true("kappa" %in% s2$truth$subjects$term)
})

test_that("smoke = TRUE runs into its own directory", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  grid <- data.frame(n_subjects = c(3L, 4L, 5L), n_trials = 10L)
  mock <- grid_mock_fitter()

  out <- suppressMessages(grid_run(dir, mock,
    grid = grid, reps = 5L,
    smoke = TRUE
  ))

  expect_equal(nrow(attr(out, "cells")), 4L)
  expect_setequal(out$condition, c("row-1", "row-2"))
  expect_true(dir.exists(file.path(dir, "smoke")))
  expect_length(list.files(dir, pattern = "\\.rds$"), 0L)
  expect_true(length(list.files(file.path(dir, "smoke"))) > 0L)
})

test_that("the preflight runs first with a short chain and can stop the grid", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  suppressMessages(grid_run(dir, mock, reps = 1L))
  first <- mock$calls$log[[1L]]
  expect_equal(first$chains, 1)
  expect_equal(first$iter, 200)
  expect_true(file.exists(file.path(dir, "preflight.rds")))

  dir2 <- withr::local_tempdir()
  failing <- grid_mock_fitter(fail_on = function(n, data, dots) n == 1L)
  expect_error(
    suppressMessages(grid_run(dir2, failing, reps = 1L)),
    "preflight"
  )
  expect_length(list.files(dir2, pattern = "^cell-.*-rep-[0-9]+\\.rds$"), 0L)

  dir3 <- withr::local_tempdir()
  mock3 <- grid_mock_fitter()
  suppressMessages(grid_run(dir3, mock3, reps = 1L, preflight = FALSE))
  expect_false(file.exists(file.path(dir3, "preflight.rds")))
  expect_identical(mock3$calls$n, 2L)
})

test_that("a failing cell is recorded and the rest is scored", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  # fail the second cell fit (call 3: preflight, cell 1, cell 2)
  failing <- grid_mock_fitter(fail_on = function(n, data, dots) n == 3L)

  expect_warning(
    out <- suppressMessages(grid_run(dir, failing, reps = 1L)),
    "row-2"
  )
  cells <- attr(out, "cells")
  expect_equal(cells$status, c("ok", "error"))
  expect_setequal(out$condition, "row-1")
})

test_that("model may be a function of the grid row", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  seen <- new.env()
  seen$n <- 0L
  builder <- function(row) {
    seen$n <- seen$n + 1L
    bmm::mixture2p(resp_error = "y")
  }
  out <- suppressMessages(recovery_grid(
    builder,
    grid = small_grid(),
    pars = c(kappa = log(8), thetat = stats::qlogis(0.75)),
    dir = dir, reps = 1L, seed = 1, preflight = FALSE,
    .fitter = grid_mock_fitter()$fitter
  ))
  expect_equal(seen$n, 2L)
  expect_s3_class(out, "bmmtools_recovery")
})

test_that("summary groups by condition and the plot facets by it", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  out <- suppressMessages(grid_run(dir, grid_mock_fitter()))

  s <- summary(out)
  expect_equal(names(s)[[1L]], "condition")
  expect_equal(nrow(s), 2L * 2L * 2L)
  expect_setequal(s$n_converged, 2L)

  skip_if_not_installed("ggplot2")
  p <- plot_recovery(out, facet_by = "condition")
  expect_s3_class(p, "ggplot")
})

test_that("a grid without its required columns errors", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  expect_error(
    grid_run(dir, grid_mock_fitter(), grid = data.frame(n_trials = 10L)),
    "n_subjects"
  )
})

test_that("recover() still fills condition with NA when absent", {
  estimates <- fake_estimates("a", estimate = 1)
  truth <- fake_truth("a", true_value = 1)
  out <- recover(estimates, truth, scale = "link")
  expect_true("condition" %in% names(out))
  expect_identical(out$condition, NA_character_)
  expect_false("condition" %in% names(summary(out)))
})

# correlated truths (spec 5, section 5.1) ------------------------------

test_that("cor_ columns set a correlation per row, in either name order", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  grid <- data.frame(
    n_subjects = c(3L, 3L), n_trials = c(10L, 10L),
    cor_kappa__thetat = c(0, 0.5)
  )
  suppressMessages(grid_run(dir, grid_mock_fitter(), grid = grid, reps = 1L))
  s1 <- readRDS(file.path(dir, "cell-1-rep-1-sim.rds"))
  s2 <- readRDS(file.path(dir, "cell-2-rep-1-sim.rds"))
  expect_equal(s1$truth$cor$true_value, 0)
  expect_equal(s2$truth$cor$true_value, 0.5)

  dir2 <- withr::local_tempdir()
  names(grid)[3] <- "cor_thetat__kappa"
  suppressMessages(grid_run(dir2, grid_mock_fitter(), grid = grid, reps = 1L))
  r2 <- readRDS(file.path(dir2, "cell-2-rep-1-sim.rds"))
  expect_equal(r2$truth$cor$term, "kappa__thetat")
  expect_equal(r2$truth$cor$true_value, 0.5)

  dir3 <- withr::local_tempdir()
  names(grid)[3] <- "cor_kappa__zeta"
  expect_error(
    suppressMessages(grid_run(dir3, grid_mock_fitter(), grid = grid)),
    "cor_kappa__zeta"
  )
})

test_that("covariates reach every cell's data and truth", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  grid <- data.frame(n_subjects = 3L, n_trials = 10L, cor_G__kappa = 0.3)
  mock <- grid_mock_fitter()
  base <- diag(2)
  dimnames(base) <- rep(list(c("kappa", "thetat")), 2)
  base["kappa", "thetat"] <- base["thetat", "kappa"] <- 0.4
  suppressMessages(grid_run(
    dir, mock,
    grid = grid, reps = 1L, cors = base,
    covariates = list(G = c(mean = 0, sd = 1))
  ))
  sim <- readRDS(file.path(dir, "cell-1-rep-1-sim.rds"))
  expect_true("G" %in% names(sim$data))
  cors <- sim$truth$cor
  cor_of <- function(term) cors$true_value[cors$term == term]
  # the column adds a pair to the default matrix and keeps its entries
  expect_equal(cor_of("G__kappa"), 0.3)
  expect_equal(cor_of("kappa__thetat"), 0.4)
})

test_that("function-valued truths use the cell seed and reproduce", {
  skip_if_not_installed("bmm")
  run <- function(dir) {
    suppressMessages(recovery_grid(
      bmm::mixture2p(resp_error = "y"),
      grid = small_grid(),
      pars = function(row) {
        c(kappa = stats::rnorm(1, log(8), 0.2), thetat = 1)
      },
      sds = c(kappa = 0.3, thetat = 0.4),
      cors = function(row) {
        r <- stats::runif(1, 0.2, 0.8)
        matrix(c(1, r, r, 1), 2, dimnames = rep(list(c("kappa", "thetat")), 2))
      },
      dir = dir, reps = 2L, seed = 5, preflight = FALSE,
      .fitter = grid_mock_fitter()$fitter
    ))
  }
  dir_a <- withr::local_tempdir()
  dir_b <- withr::local_tempdir()
  run(dir_a)
  run(dir_b)
  a11 <- readRDS(file.path(dir_a, "cell-1-rep-1-sim.rds"))
  b11 <- readRDS(file.path(dir_b, "cell-1-rep-1-sim.rds"))
  a21 <- readRDS(file.path(dir_a, "cell-2-rep-1-sim.rds"))
  expect_identical(a11$truth, b11$truth)
  expect_false(identical(a11$cors, a21$cors))
  expect_false(identical(a11$pars, a21$pars))
})

test_that("subjects = 'fixed' reuses rep 1's realised truths", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  suppressMessages(recovery_grid(
    bmm::mixture2p(resp_error = "y"),
    grid = small_grid(),
    pars = function(row) c(kappa = stats::rnorm(1, 2, 0.2), thetat = 1),
    sds = c(kappa = 0.3, thetat = 0.4),
    cors = function(row) {
      r <- stats::runif(1, 0.2, 0.8)
      matrix(c(1, r, r, 1), 2, dimnames = rep(list(c("kappa", "thetat")), 2))
    },
    covariates = list(G = c(mean = 0, sd = 1)),
    subjects = "fixed",
    dir = dir, reps = 2L, seed = 5, preflight = FALSE,
    .fitter = grid_mock_fitter()$fitter
  ))
  s11 <- readRDS(file.path(dir, "cell-1-rep-1-sim.rds"))
  s12 <- readRDS(file.path(dir, "cell-1-rep-2-sim.rds"))
  expect_identical(s12$pars, s11$pars)
  expect_identical(s12$cors, s11$cors)
  expect_identical(s12$truth$subjects, s11$truth$subjects)
  expect_identical(s12$truth$covariates, s11$truth$covariates)
  expect_identical(s12$truth$cor, s11$truth$cor)
})

test_that("re_cor reaches the default formula", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  seen <- new.env()
  fitter <- function(formula, data, model, prior = NULL, ...) {
    seen$formula <- formula
    mock$fitter(formula, data, model, prior, ...)
  }
  suppressMessages(recovery_grid(
    bmm::mixture2p(resp_error = "y"),
    grid = small_grid()[1, ], pars = c(kappa = 2, thetat = 1),
    sds = c(kappa = 0.3, thetat = 0.4), re_cor = "all",
    dir = dir, reps = 1L, seed = 1, preflight = FALSE, .fitter = fitter
  ))
  expect_equal(deparse(seen$formula$kappa), "kappa ~ 1 + (1 | p | id)")
})

test_that("a simulation file from before 5.1 is upgraded on resume", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  suppressMessages(grid_run(dir, mock, reps = 1L))
  path <- file.path(dir, "cell-1-rep-1-sim.rds")
  old <- readRDS(path)
  old$truth <- old$truth[c("population", "subjects")]
  old$cors <- NULL
  old$covariates <- NULL
  saveRDS(old, path)

  sim <- upgrade_simulation(readRDS(path))
  expect_equal(sim$truth$sd$term, c("kappa", "thetat"))
  expect_equal(sim$truth$cor$term, "kappa__thetat")
  expect_equal(sim$truth$cor$true_value, 0)
  expect_equal(nrow(sim$truth$covariates), 0L)
  expect_equal(unname(sim$cors), diag(2))

  out <- suppressMessages(grid_run(dir, mock, reps = 1L))
  expect_s3_class(out, "bmmtools_recovery")
})

test_that("a row-wise model must keep one link table across rows", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  builder <- function(row) {
    model <- bmm::mixture2p(resp_error = "y")
    if (row$n_subjects == 4L) model$links$kappa <- "identity"
    model
  }
  expect_error(
    suppressMessages(recovery_grid(
      builder,
      grid = small_grid(),
      pars = c(kappa = log(8), thetat = stats::qlogis(0.75)),
      dir = dir, reps = 1L, seed = 1, preflight = FALSE,
      .fitter = grid_mock_fitter()$fitter
    )),
    "links"
  )
})

# the sidecar, sd level and correlations (spec 5, section 5.3, D31) -----

grid_fit_files <- function(dir) {
  list.files(dir, pattern = "^(cell-[0-9]+-rep-[0-9]+|preflight)\\.rds$",
    full.names = TRUE
  )
}

# everything but the timings, which differ between two runs
grid_result_parts <- function(out) {
  cells <- attr(out, "cells")
  cells$elapsed <- NULL
  list(
    # the columns only; the attributes are compared one by one below
    rows = unclass(out)[names(out)],
    cells = cells,
    correlations = attr(out, "correlations"),
    subject_means = attr(out, "subject_means")
  )
}

test_that("the defaults extract no correlations and keep today's levels", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  out <- suppressMessages(grid_run(dir, grid_mock_fitter()))
  expect_null(attr(out, "correlations"))
  expect_setequal(out$level, c("population", "subject"))
  sidecar <- readRDS(file.path(dir, "cell-1-rep-1-est.rds"))
  expect_equal(sidecar$levels, c("population", "subject"))
  expect_null(sidecar$correlations)
  expect_null(sidecar$cor_estimates)
})

test_that("each cell writes a sidecar with the documented fields", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  suppressMessages(grid_run(
    dir, grid_mock_fitter(),
    correlations = c("draws", "point")
  ))

  expect_equal(
    sort(list.files(dir, pattern = "-est\\.rds$")),
    sort(sprintf("cell-%d-rep-%d-est.rds", c(1, 2, 1, 2), c(1, 1, 2, 2)))
  )
  sidecar <- readRDS(file.path(dir, "cell-2-rep-1-est.rds"))
  expect_named(sidecar, c(
    "key", "bmmtools_version", "levels", "correlations", "cor_scale",
    "estimates", "cor_estimates", "subject_means"
  ))

  sim <- readRDS(file.path(dir, "cell-2-rep-1-sim.rds"))
  key <- cache_key(
    recovery_formula(sim$model), sim$data, sim$model, NULL,
    list(seed = cell_seed(100, 2L, 1L), chains = 2, iter = 400)
  )$key
  expect_identical(sidecar$key, key)
  expect_identical(
    sidecar$bmmtools_version,
    as.character(utils::packageVersion("bmmtools"))
  )
  expect_equal(sidecar$levels, c("population", "subject"))
  expect_equal(sidecar$correlations, c("draws", "point"))
  expect_equal(sidecar$cor_scale, "link")
  expect_setequal(sidecar$estimates$level, c("population", "subject"))
  expect_setequal(sidecar$cor_estimates$estimator, c("draws", "point"))
  expect_equal(unique(sidecar$cor_estimates$scale), "link")

  means <- sidecar$subject_means
  expect_named(means, c("id", "term", "mean", "median"))
  expect_equal(nrow(means), 4L * 2L)
  mock_fit <- structure(
    list(parameters = c("kappa", "thetat"), ids = as.character(1:4)),
    class = "mockfit"
  )
  draws <- extract_subject_draws(mock_fit)
  at <- means$id == "3" & means$term == "thetat"
  expect_equal(means$mean[at], mean(draws[, , "3", "thetat"]))
  expect_equal(means$median[at], stats::median(draws[, , "3", "thetat"]))
})

test_that("a resume with the fits deleted reads the sidecars, fitting none", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  first <- suppressMessages(grid_run(
    dir, grid_mock_fitter(),
    correlations = c("draws", "point"),
    covariates = list(G = c(mean = 0, sd = 1))
  ))
  unlink(grid_fit_files(dir))
  expect_length(grid_fit_files(dir), 0L)

  mock <- grid_mock_fitter()
  second <- suppressMessages(grid_run(
    dir, mock,
    correlations = c("draws", "point"),
    covariates = list(G = c(mean = 0, sd = 1))
  ))

  expect_identical(mock$calls$n, 0L)
  expect_length(grid_fit_files(dir), 0L)
  expect_identical(grid_result_parts(second), grid_result_parts(first))
  expect_true(all(attr(second, "cells")$status == "ok"))
})

test_that("a different bmmtools version re-extracts through the fit", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  suppressMessages(grid_run(dir, grid_mock_fitter(), correlations = "draws"))
  unlink(grid_fit_files(dir))

  local_mocked_bindings(bmmtools_version = function() "99.0.0")
  mock <- grid_mock_fitter()
  out <- suppressMessages(grid_run(dir, mock, correlations = "draws"))

  # the preflight (the first sidecar is stale) plus four cells
  expect_identical(mock$calls$n, 5L)
  expect_equal(
    readRDS(file.path(dir, "cell-1-rep-1-est.rds"))$bmmtools_version,
    "99.0.0"
  )
  expect_s3_class(attr(out, "correlations"), "bmmtools_cor_recovery")
})

test_that("a sidecar lacking a requested level, estimator or scale is redone", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  suppressMessages(grid_run(dir, grid_mock_fitter(), correlations = "draws"))

  # an estimator the sidecar lacks
  unlink(grid_fit_files(dir))
  mock <- grid_mock_fitter()
  out <- suppressMessages(grid_run(
    dir, mock,
    correlations = c("draws", "point")
  ))
  expect_identical(mock$calls$n, 5L)
  expect_setequal(attr(out, "correlations")$estimator, c("draws", "point"))

  # a level the sidecar lacks
  unlink(grid_fit_files(dir))
  mock <- grid_mock_fitter()
  suppressMessages(grid_run(
    dir, mock,
    levels = c("population", "subject", "sd"),
    correlations = c("draws", "point")
  ))
  expect_identical(mock$calls$n, 5L)

  # another correlation scale
  unlink(grid_fit_files(dir))
  mock <- grid_mock_fitter()
  out <- suppressMessages(grid_run(
    dir, mock,
    correlations = "point", cor_scale = "natural"
  ))
  expect_identical(mock$calls$n, 5L)
  expect_equal(attr(out, "scale", exact = TRUE), "natural")
  expect_equal(unique(attr(out, "correlations")$scale), "natural")

  # a request the sidecar covers, and one it holds more than
  unlink(grid_fit_files(dir))
  mock <- grid_mock_fitter()
  subset <- suppressMessages(grid_run(dir, mock, levels = "population"))
  expect_identical(mock$calls$n, 0L)
  expect_setequal(subset$level, "population")
  expect_null(attr(subset, "correlations"))
})

test_that("a sidecar that cannot be read is replaced", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  suppressMessages(grid_run(dir, grid_mock_fitter(), reps = 1L))
  path <- file.path(dir, "cell-2-rep-1-est.rds")
  writeLines("not an rds file", path)
  mock <- grid_mock_fitter()
  out <- suppressMessages(grid_run(dir, mock, reps = 1L))
  # the fit is still cached, so fit_cached() reads it without the fitter
  expect_identical(mock$calls$n, 0L)
  expect_named(readRDS(path)[1:2], c("key", "bmmtools_version"))
  expect_true(all(attr(out, "cells")$status == "ok"))
})

test_that("levels may include 'sd', scored on the link scale", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  out <- suppressMessages(grid_run(
    dir, grid_mock_fitter(),
    levels = c("population", "subject", "sd")
  ))
  expect_setequal(out$level, c("population", "subject", "sd"))
  sd_rows <- out[out$level == "sd", ]
  expect_equal(nrow(sd_rows), 2L * 2L * 2L)
  expect_true(all(sd_rows$scale == "link"))
  expect_equal(
    sd_rows$true_value[sd_rows$term == "kappa" & sd_rows$condition == "row-1"],
    c(0.3, 0.3)
  )
  expect_equal(attr(out, "scale"), "natural")
  expect_equal(unique(out$scale[out$level == "population"]), "natural")

  s <- summary(out)
  expect_true("sd" %in% s$level)
})

test_that("score_cells attaches correlations and subject means", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  out <- suppressMessages(grid_run(
    dir, grid_mock_fitter(),
    correlations = c("draws", "point"),
    covariates = list(G = c(mean = 0, sd = 1))
  ))

  cors <- attr(out, "correlations")
  expect_s3_class(cors, "bmmtools_cor_recovery")
  expect_setequal(cors$condition, c("row-1", "row-2"))
  expect_setequal(cors$replication, c(1L, 2L))
  expect_setequal(cors$term, c("G__kappa", "G__thetat", "kappa__thetat"))
  expect_equal(nrow(cors), 2L * 2L * 3L * 2L)

  # sample_value is the in-sample correlation of that cell's true values
  sim <- readRDS(file.path(dir, "cell-2-rep-2-sim.rds"))
  wide <- subjects_wide(sim$truth$subjects)
  at <- cors$condition == "row-2" & cors$replication == 2L &
    cors$term == "kappa__thetat"
  expect_equal(
    unique(cors$sample_value[at]),
    stats::cor(wide$kappa, wide$thetat)
  )
  expect_equal(unique(cors$true_value[at]), 0)

  means <- attr(out, "subject_means")
  expect_named(means, c(
    "condition", "replication", "id", "term", "covariate", "mean",
    "median", "true_value"
  ))
  # (3 + 4) subjects, two replications, two parameters and one covariate
  expect_equal(nrow(means), 7L * 2L * 3L)
  at_g <- means$term == "G" & means$condition == "row-2" &
    means$replication == 2L
  g <- means[at_g, ]
  expect_true(all(g$covariate))
  expect_equal(g$mean, g$true_value)
  expect_equal(
    g$true_value,
    sim$truth$covariates$true_value[match(g$id, sim$truth$covariates$id)]
  )
  expect_equal(attr(means, "links"), unlist(sim$model$links))

  # the cells and grid attributes are still there
  expect_equal(nrow(attr(out, "cells")), 4L)
  expect_equal(attr(out, "grid"), small_grid())
})

test_that("recovery_grid validates levels, correlations and cor_scale", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  expect_error(grid_run(dir, mock, levels = "bogus"), "levels")
  expect_error(grid_run(dir, mock, levels = character(0)), "levels")
  expect_error(grid_run(dir, mock, correlations = "bogus"), "correlations")
  expect_error(
    grid_run(dir, mock, correlations = character(0)),
    "correlations"
  )
  expect_error(grid_run(dir, mock, correlations = NA), "correlations")
  expect_error(
    grid_run(dir, mock, correlations = "draws", cor_scale = "bogus"),
    "cor_scale"
  )
  expect_error(
    grid_run(dir, mock, correlations = "model", cor_scale = "natural"),
    "link scale"
  )
  # the default uncorrelated formula cannot give model correlations
  expect_error(
    grid_run(dir, mock, correlations = "model"),
    "re_cor"
  )
  expect_identical(mock$calls$n, 0L)
})

test_that("a grid with nothing to score at its levels says so", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  one_sd <- function(dir, ...) {
    suppressMessages(recovery_grid(
      bmm::mixture2p(resp_error = "y"),
      grid = small_grid()[1, ],
      pars = c(kappa = log(8), thetat = stats::qlogis(0.75)),
      dir = dir, reps = 1L, seed = 1, preflight = FALSE, ...,
      .fitter = grid_mock_fitter()$fitter
    ))
  }

  # one varying parameter and no covariate: no correlation pair
  expect_warning(
    out <- one_sd(dir, sds = c(kappa = 0.3), correlations = "draws"),
    "no cell has a correlation pair"
  )
  expect_null(attr(out, "correlations"))
  expect_setequal(out$level, c("population", "subject"))

  # nothing varies, so there is no subject level to score
  dir2 <- withr::local_tempdir()
  expect_error(one_sd(dir2, levels = "subject"), "Nothing to score")
})
