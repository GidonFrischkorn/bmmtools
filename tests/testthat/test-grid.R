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
  list.files(dir,
    pattern = "^(cell-[0-9]+-rep-[0-9]+|preflight)\\.rds$",
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

# the task dimension (spec 5, section 5.4) --------------------------------

task_grid_run <- function(dir, mock = grid_mock_fitter(), grid = small_grid(),
                          reps = 1L, ...) {
  suppressMessages(recovery_grid(
    bmm::mixture2p(resp_error = "y"),
    grid = grid,
    pars = c(kappa = log(8), kappa_task2 = log(4), thetat = 0.5),
    dir = dir,
    reps = reps,
    sds = c(kappa = 0.3, thetat_task1 = 0.4),
    tasks = c("1", "2"),
    seed = 200,
    ...,
    .fitter = mock$fitter
  ))
}

test_that("a grid with tasks scores task terms at every level", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  seen <- new.env()
  mock <- grid_mock_fitter()
  fitter <- function(formula, data, model, prior = NULL, ...) {
    seen$formula <- formula
    mock$fitter(formula, data, model, prior, ...)
  }
  out <- task_grid_run(
    dir, list(fitter = fitter),
    re_cor = "within",
    levels = c("population", "subject", "sd"),
    correlations = c("draws", "point")
  )

  expect_equal(
    deparse(seen$formula$kappa), "kappa ~ 0 + task + (0 + task | id)"
  )
  terms <- c("kappa_task1", "kappa_task2", "thetat_task1", "thetat_task2")
  expect_setequal(out$term[out$level == "population"], terms)
  expect_setequal(
    out$term[out$level == "subject"],
    c("kappa_task1", "kappa_task2", "thetat_task1")
  )
  expect_setequal(
    out$term[out$level == "sd"],
    c("kappa_task1", "kappa_task2", "thetat_task1")
  )
  pop <- out[out$level == "population" & out$condition == "row-1", ]
  expect_equal(
    pop$true_value[match(terms, pop$term)],
    c(8, 4, stats::plogis(0.5), stats::plogis(0.5))
  )

  cors <- attr(out, "correlations")
  expect_s3_class(cors, "bmmtools_cor_recovery")
  expect_setequal(cors$estimator, c("draws", "point"))
  expect_true("kappa_task1__kappa_task2" %in% cors$term)
  sim <- readRDS(file.path(dir, "cell-1-rep-1-sim.rds"))
  expect_identical(sim$tasks, c("1", "2"))
  expect_named(sim$data, c("id", "task", "y"))

  table <- subject_table(out)
  expect_true(all(c(
    "true_kappa_task1", "est_kappa_task1", "true_thetat_task2",
    "est_thetat_task2"
  ) %in% names(table)))
  expect_equal(nrow(table), 3L + 4L)
})

test_that("grid columns may set bare and full task terms", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  grid <- data.frame(
    n_subjects = 3L, n_trials = 2L,
    kappa = c(1, 2), kappa_task1 = c(1, 5),
    sd_thetat = c(0.2, 0.2), sd_kappa_task2 = c(0.1, 0.6),
    cor_kappa_task2__thetat_task1 = c(0, 0.4)
  )
  task_grid_run(dir, grid = grid, preflight = FALSE)
  s1 <- readRDS(file.path(dir, "cell-1-rep-1-sim.rds"))
  s2 <- readRDS(file.path(dir, "cell-2-rep-1-sim.rds"))

  # a bare column sets every task, a full-term column wins for its task
  kappas <- c("kappa_task1", "kappa_task2")
  expect_equal(s1$pars[kappas], c(kappa_task1 = 1, kappa_task2 = 1))
  expect_equal(s2$pars[kappas], c(kappa_task1 = 5, kappa_task2 = 2))
  expect_equal(
    s2$sds,
    c(
      kappa_task1 = 0.3, kappa_task2 = 0.6,
      thetat_task1 = 0.2, thetat_task2 = 0.2
    )
  )
  expect_equal(
    s2$truth$cor$true_value[s2$truth$cor$term == "kappa_task2__thetat_task1"],
    0.4
  )

  # bare or unknown terms in a cor_ column are refused with tasks
  bad <- data.frame(n_subjects = 3L, n_trials = 2L, cor_kappa__thetat = 0.2)
  expect_error(
    task_grid_run(withr::local_tempdir(), grid = bad, preflight = FALSE),
    "cor_kappa__thetat"
  )
})

test_that("subjects = 'fixed' with tasks reuses rep 1's task values", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  suppressMessages(recovery_grid(
    bmm::mixture2p(resp_error = "y"),
    grid = small_grid(),
    pars = function(row) c(kappa = stats::rnorm(1, 2, 0.2), thetat = 1),
    sds = c(kappa = 0.3),
    tasks = c("1", "2"), subjects = "fixed",
    dir = dir, reps = 2L, seed = 5, preflight = FALSE,
    .fitter = grid_mock_fitter()$fitter
  ))
  s11 <- readRDS(file.path(dir, "cell-1-rep-1-sim.rds"))
  s12 <- readRDS(file.path(dir, "cell-1-rep-2-sim.rds"))
  expect_named(s11$pars, c(
    "kappa_task1", "kappa_task2", "thetat_task1", "thetat_task2"
  ))
  expect_identical(s12$pars, s11$pars)
  expect_identical(s12$sds, s11$sds)
  expect_identical(s12$truth$subjects, s11$truth$subjects)
  expect_false(identical(s12$data, s11$data))
})

test_that("formula may be a function of the grid row, called once per row", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  rows <- new.env()
  rows$seen <- list()
  formula <- function(row) {
    rows$seen[[length(rows$seen) + 1L]] <- row
    recovery_formula(
      bmm::mixture2p(resp_error = "y"),
      re_cor = if (row$n_subjects == 3L) "none" else "all",
      task_col = "task"
    )
  }
  calls <- new.env()
  calls$formulas <- list()
  mock <- grid_mock_fitter()
  fitter <- function(formula, data, model, prior = NULL, ...) {
    calls$formulas[[length(calls$formulas) + 1L]] <- formula
    mock$fitter(formula, data, model, prior, ...)
  }
  task_grid_run(dir, list(fitter = fitter), reps = 2L, formula = formula)

  expect_length(rows$seen, 2L)
  expect_equal(rows$seen[[1L]], small_grid()[1, , drop = FALSE])
  deparsed <- vapply(calls$formulas, function(f) deparse(f$kappa), "")
  # preflight (row 1), then rows 1, 2 in each replication
  expect_equal(deparsed, c(
    "kappa ~ 0 + task + (0 + task || id)",
    rep(c(
      "kappa ~ 0 + task + (0 + task || id)",
      "kappa ~ 0 + task + (0 + task | p | id)"
    ), 2L)
  ))

  # a fixed bmmformula is used as given for every row
  fixed <- recovery_formula(
    bmm::mixture2p(resp_error = "y"),
    re_cor = "within", task_col = "task"
  )
  calls$formulas <- list()
  task_grid_run(
    withr::local_tempdir(), list(fitter = fitter),
    formula = fixed, preflight = FALSE
  )
  expect_true(all(vapply(calls$formulas, identical, TRUE, fixed)))

  not_formula <- function(row) "kappa ~ 1"
  err <- expect_error(
    task_grid_run(withr::local_tempdir(), formula = not_formula),
    "bmmformula"
  )
  expect_match(conditionMessage(err), "row 1")
})

test_that("tasks change the data and so the cache key", {
  skip_if_not_installed("bmm")
  model <- bmm::mixture2p(resp_error = "y")
  sim_for <- function(tasks) {
    simulate_recovery(
      model, c(kappa = 1, thetat = 0),
      n_subjects = 2, n_trials = 3, tasks = tasks, seed = 1
    )
  }
  a <- sim_for(c("1", "2"))
  b <- sim_for(c("1", "3"))
  formula <- recovery_formula(model, task_col = "task")
  expect_false(identical(
    cache_key(formula, a$data, model, NULL, list())$key,
    cache_key(formula, b$data, model, NULL, list())$key
  ))
})

test_that("grid argument checks know about tasks and re_cor = 'within'", {
  skip_if_not_installed("bmm")
  mock <- grid_mock_fitter()
  err <- expect_error(
    task_grid_run(withr::local_tempdir(), mock, correlations = "model"),
    "re_cor"
  )
  expect_match(conditionMessage(err), "within")
  err <- expect_error(
    grid_run(
      withr::local_tempdir(), mock,
      correlations = "model", re_cor = "within"
    ),
    "re_cor"
  )
  expect_no_match(conditionMessage(err), "within")
  expect_error(
    task_grid_run(withr::local_tempdir(), mock, tasks = 1:2),
    "tasks"
  )
  expect_error(
    task_grid_run(withr::local_tempdir(), mock, task_col = "id"),
    "task_col"
  )
  expect_identical(mock$calls$n, 0L)
})

test_that("an older simulation file gains NULL tasks on resume", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  suppressMessages(grid_run(dir, grid_mock_fitter(), reps = 1L))
  path <- file.path(dir, "cell-1-rep-1-sim.rds")
  old <- unclass(readRDS(path))
  old <- structure(
    old[setdiff(names(old), c("tasks", "task_col"))],
    class = "bmmtools_simulation"
  )
  sim <- upgrade_simulation(old)
  expect_true(all(c("tasks", "task_col") %in% names(sim)))
  expect_null(sim$tasks)
  expect_identical(sim$truth, old$truth)
})

# subject-wise ML (milestone 8, stage 8.4) ---------------------------------

# the Bayesian mock is `mock`; the ML mock goes in through `ml = list(.fitter)`
# because the grid's `.fitter` returns fits whose coefficients are not the
# `<par>_id<level>` shape a no-pooling fit has
ml_grid_run <- function(dir, mock = grid_mock_fitter(),
                        ml_mock = ml_mock_fitter(), ml_args = list(),
                        reps = 2L, ...) {
  quietly(suppressMessages(grid_run(
    dir, mock,
    reps = reps,
    ml = c(list(.fitter = ml_mock$fitter), ml_args),
    ...
  )))
}

#' Muffle the balance warning, and only it
#'
#' These fixtures fit at 10 and 20 trials a subject, where some MLEs run to
#' the boundary, so two estimators really are scored on different subjects
#' and the warning is the right answer --- it is asserted on its own below.
#' Every other warning still surfaces.
#'
#' @noRd
quietly <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("same subjects", conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}

test_that("the grid runs the optim route with no ML fitter at all", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  # no `.fitter` in the ml list: the optim route calls no fitter, so this
  # is the one ML path a runner can exercise end to end. The hierarchical
  # side is still the mock.
  out <- quietly(suppressMessages(
    grid_run(dir, mock, ml = list(method = "optim"))
  ))

  expect_setequal(out$estimator, c("bayes", "ml"))
  expect_setequal(out$estimator[out$level == "population"], "bayes")
  cells <- attr(out, "ml_cells")
  expect_equal(nrow(cells), 4L)
  # the fit ran in every cell; whether every subject converged is another
  # question, and at this grid's 10 and 20 trials the answer is no
  expect_true(all(cells$status == "ok"))
  expect_equal(cells$n_subjects, c(3L, 4L, 3L, 4L))
  expect_true(all(cells$n_converged <= cells$n_subjects))

  # decision 40 through the grid: a subject whose MLE runs to the boundary
  # keeps its row with estimate NA, so the two estimators are scored on the
  # same subjects and the difference between them is not selection
  sub <- out[out$level == "subject", ]
  expect_equal(sum(sub$estimator == "ml"), sum(sub$estimator == "bayes"))
  failed <- sub[sub$estimator == "ml" & !sub$converged, ]
  expect_true(all(is.na(failed$estimate)))
})

test_that("changing the ML route re-fits rather than reusing the sidecar", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  first <- quietly(suppressMessages(
    grid_run(dir, grid_mock_fitter(), ml = list(method = "optim"))
  ))
  # the sidecar records `method`, so the same directory asked for the other
  # route must not answer with the rows it already has
  again <- quietly(suppressMessages(
    grid_run(dir, grid_mock_fitter(), ml = list(method = "optim"))
  ))
  # the rows, not the attributes: `elapsed` is 0 on the resume precisely
  # because nothing was fitted, so as.data.frame() is not enough --- a
  # tibble subclass carries its attributes through it
  plain <- function(x) as.data.frame(lapply(x, identity))
  expect_equal(plain(first), plain(again))
  expect_true(all(attr(again, "ml_cells")$elapsed == 0))

  ml_mock <- ml_mock_fitter()
  switched <- quietly(suppressMessages(grid_run(
    dir, grid_mock_fitter(),
    ml = list(method = "stan", .fitter = ml_mock$fitter)
  )))
  expect_true(ml_mock$calls$n > 0L)
  expect_false(identical(
    first$estimate[first$estimator == "ml"],
    switched$estimate[switched$estimator == "ml"]
  ))
})

test_that("ml = TRUE scores both estimators at the subject level", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  ml_mock <- ml_mock_fitter(estimates = c(kappa_id1 = 0.4, thetat_id2 = -0.3))

  out <- ml_grid_run(dir, mock, ml_mock)

  expect_s3_class(out, "bmmtools_recovery")
  expect_setequal(out$estimator, c("bayes", "ml"))
  # an ML fit has no population level: those rows stay hierarchical only
  expect_setequal(out$estimator[out$level == "population"], "bayes")
  expect_setequal(out$estimator[out$level == "subject"], c("bayes", "ml"))
  # every cell has as many ML subject rows as Bayesian ones
  sub <- out[out$level == "subject", ]
  counts <- table(sub$condition, sub$replication, sub$estimator)
  expect_equal(unname(counts[, , "ml"]), unname(counts[, , "bayes"]))
  # the ML rows are what fit_ml() produced, scored against the same truth
  ml <- sub[sub$estimator == "ml", ]
  expect_true(all(ml$ci_method == "laplace"))
  expect_true(all(is.na(ml$rhat)))
  bayes <- sub[sub$estimator == "bayes", ]
  key <- c("condition", "replication", "id", "term")
  merged <- merge(
    ml[c(key, "true_value")], bayes[c(key, "true_value")],
    by = key
  )
  expect_equal(merged$true_value.x, merged$true_value.y)
  # one ML fit per cell, none for the preflight
  expect_identical(ml_mock$calls$n, 4L)
  expect_identical(mock$calls$n, 5L)
  # the ML fit is cached beside the cell's other files
  expect_true(all(file.exists(
    file.path(dir, sprintf(
      "cell-%d-rep-%d-ml.rds", c(1, 2, 1, 2), c(1, 1, 2, 2)
    ))
  )))

  s <- summary(out)
  expect_true("estimator" %in% names(s))
  expect_setequal(s$estimator[s$level == "subject"], c("bayes", "ml"))
  expect_setequal(s$estimator[s$level == "population"], "bayes")
})

test_that("ml = TRUE leaves the cells attribute as it was", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  out <- ml_grid_run(dir)
  cells <- attr(out, "cells")
  expect_equal(nrow(cells), 4L)
  expect_named(cells, c(
    "condition", "replication", "n_subjects", "n_trials", "seed", "file",
    "status", "elapsed", "converged"
  ))
})

test_that("ml_cells has one row per cell with the subject counts", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  # subject 2's kappa leaves the link range in every cell
  out <- ml_grid_run(
    dir,
    ml_mock = ml_mock_fitter(estimates = c(kappa_id2 = 50))
  )

  ml_cells <- attr(out, "ml_cells")
  expect_s3_class(ml_cells, "tbl_df")
  expect_named(ml_cells, c(
    "condition", "replication", "status", "elapsed", "n_subjects",
    "n_converged", "converged", "file"
  ))
  expect_equal(ml_cells$condition, c("row-1", "row-2", "row-1", "row-2"))
  expect_equal(ml_cells$replication, c(1L, 1L, 2L, 2L))
  expect_true(all(ml_cells$status == "ok"))
  expect_equal(ml_cells$n_subjects, c(3L, 4L, 3L, 4L))
  expect_equal(ml_cells$n_converged, c(2L, 3L, 2L, 3L))
  expect_false(any(ml_cells$converged))
  expect_true(all(file.exists(ml_cells$file)))

  # the failed subject keeps its row, unscored, so n differs visibly
  ml <- out[out$estimator == "ml" & out$term == "kappa", ]
  expect_true(all(is.na(ml$estimate[ml$id == "2"])))
  expect_false(any(ml$converged[ml$id == "2"]))
  expect_true(all(ml$converged[ml$id != "2"]))

  # without ml there is no such attribute
  out2 <- suppressMessages(grid_run(withr::local_tempdir(), grid_mock_fitter()))
  expect_null(attr(out2, "ml_cells"))
})

test_that("the sidecar stores the ML record and rows; a resume reads them", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  first <- ml_grid_run(dir, ml_args = list(draws = 200))

  sidecar <- readRDS(file.path(dir, "cell-2-rep-1-est.rds"))
  expect_named(sidecar, c(
    "key", "bmmtools_version", "levels", "correlations", "cor_scale", "ml",
    "estimates", "cor_estimates", "subject_means", "ml_estimates"
  ))
  # every fit_ml() argument that changes the rows, at its default when not
  # given --- `start` and `nll` belong to the optim route and are NULL here
  expect_equal(sidecar$ml, list(
    formula = NULL, method = "stan", prior = "flat", ci_level = 0.95,
    draws = 200, max_abs_link = 20, start = NULL, nll = NULL
  ))
  expect_s3_class(sidecar$ml_estimates, "bmmtools_ml")
  expect_equal(nrow(sidecar$ml_estimates), 4L * 2L)

  # with every fit deleted, both mocks stay idle and the result is the same
  unlink(grid_fit_files(dir))
  unlink(list.files(dir, pattern = "-ml\\.(rds|key)$", full.names = TRUE))
  mock <- grid_mock_fitter()
  ml_mock <- ml_mock_fitter()
  second <- ml_grid_run(dir, mock, ml_mock, ml_args = list(draws = 200))
  expect_identical(mock$calls$n, 0L)
  expect_identical(ml_mock$calls$n, 0L)
  expect_identical(grid_result_parts(second), grid_result_parts(first))
  ml_cells <- attr(second, "ml_cells")
  expect_true(all(ml_cells$status == "ok"))
  expect_true(all(ml_cells$elapsed == 0))
})

test_that("ml added to a finished grid fits ML only; dropped, it is trimmed", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  suppressMessages(grid_run(dir, grid_mock_fitter()))
  expect_false("ml" %in% names(readRDS(file.path(dir, "cell-1-rep-1-est.rds"))))

  # the hierarchical fits are cached, so only the ML mock works
  mock <- grid_mock_fitter()
  ml_mock <- ml_mock_fitter()
  with_ml <- ml_grid_run(dir, mock, ml_mock)
  expect_identical(mock$calls$n, 0L)
  expect_identical(ml_mock$calls$n, 4L)
  expect_setequal(with_ml$estimator, c("bayes", "ml"))
  expect_true(
    "ml_estimates" %in%
      names(readRDS(file.path(dir, "cell-1-rep-1-est.rds")))
  )

  # a run without ml reads the same sidecars and drops the ML rows
  mock <- grid_mock_fitter()
  without <- suppressMessages(grid_run(dir, mock))
  expect_identical(mock$calls$n, 0L)
  expect_setequal(without$estimator, "bayes")
  expect_null(attr(without, "ml_cells"))
})

test_that("a changed ML request re-extracts from the cached ML fit", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  ml_mock <- ml_mock_fitter(estimates = c(kappa_id2 = 50))
  strict <- ml_grid_run(dir, ml_mock = ml_mock)
  expect_equal(attr(strict, "ml_cells")$n_converged, c(2L, 3L, 2L, 3L))

  # max_abs_link is applied after the fit: the fit is reused, the rows redone
  lenient <- ml_grid_run(
    dir,
    ml_mock = ml_mock, ml_args = list(max_abs_link = 100)
  )
  expect_identical(ml_mock$calls$n, 4L)
  expect_equal(attr(lenient, "ml_cells")$n_converged, c(3L, 4L, 3L, 4L))
  at <- lenient$estimator == "ml" & lenient$term == "kappa" & lenient$id == "2"
  ml <- lenient[at, ]
  expect_true(all(is.finite(ml$estimate)))
  sidecar <- readRDS(file.path(dir, "cell-1-rep-1-est.rds"))
  expect_equal(sidecar$ml$max_abs_link, 100)

  # draws enters the ML fit's own cache key, so it refits
  ml_grid_run(dir, ml_mock = ml_mock, ml_args = list(draws = 50))
  expect_identical(ml_mock$calls$n, 8L)
})

test_that("an ML failure is recorded per cell and retried on resume", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  calls <- new.env()
  calls$n <- 0L
  inner <- ml_mock_fitter()
  # the second ML fit (row-2 rep 1) fails
  flaky <- function(formula, data, model, prior = NULL, ...) {
    calls$n <- calls$n + 1L
    if (calls$n == 2L) stop("mock laplace failure")
    inner$fitter(formula, data, model, prior, ...)
  }
  mock <- grid_mock_fitter()

  # two warnings, both right: the cell whose ML fit failed, and --- because
  # that cell then has no ML rows at all while its ids appear in every
  # other cell --- the balance check, which only sees this because it is
  # keyed by the cell rather than pooled over ids
  expect_warning(
    expect_warning(
      out <- suppressMessages(grid_run(dir, mock, ml = list(.fitter = flaky))),
      "row-2 rep 1"
    ),
    "same subjects"
  )
  ml_cells <- attr(out, "ml_cells")
  expect_equal(ml_cells$status, c("ok", "error", "ok", "ok"))
  expect_true(is.na(ml_cells$n_subjects[[2L]]))
  # the hierarchical rows of that cell are scored; its ML rows are absent
  expect_true(all(attr(out, "cells")$status == "ok"))
  at <- out$condition == "row-2" & out$replication == 1L
  expect_setequal(out$estimator[at], "bayes")
  expect_setequal(out$estimator[!at], c("bayes", "ml"))
  # its sidecar has no ML record, so a resume tries again, once
  expect_false("ml" %in% names(readRDS(file.path(dir, "cell-2-rep-1-est.rds"))))
  expect_true("ml" %in% names(readRDS(file.path(dir, "cell-1-rep-1-est.rds"))))

  mock2 <- grid_mock_fitter()
  again <- suppressMessages(grid_run(dir, mock2, ml = list(.fitter = flaky)))
  expect_identical(calls$n, 5L)
  expect_identical(mock2$calls$n, 0L)
  expect_true(all(attr(again, "ml_cells")$status == "ok"))
  at <- again$condition == "row-2" & again$replication == 1L
  expect_setequal(again$estimator[at], c("bayes", "ml"))
})

test_that("a cell whose hierarchical fit fails has no ML status", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  failing <- grid_mock_fitter(fail_on = function(n, data, dots) n == 3L)
  ml_mock <- ml_mock_fitter()
  expect_warning(
    out <- ml_grid_run(dir, failing, ml_mock, reps = 1L),
    "row-2"
  )
  expect_identical(ml_mock$calls$n, 1L)
  ml_cells <- attr(out, "ml_cells")
  expect_equal(ml_cells$status, c("ok", NA_character_))
  expect_true(is.na(ml_cells$elapsed[[2L]]))
})

test_that("the ML fit gets the cell's data and seed, and flat priors", {
  skip_if_not_installed("bmm")
  skip_if_not_installed("brms")
  dir <- withr::local_tempdir()
  ml_mock <- ml_mock_fitter()
  ml_grid_run(dir, ml_mock = ml_mock, reps = 1L, ml_args = list(refresh = 0))
  last <- ml_mock$calls$last
  # the last cell is row-2 rep 1: four subjects at twenty trials
  expect_equal(nrow(last$data), 4L * 20L)
  expect_equal(last$dots$seed, cell_seed(100, 2L, 1L))
  expect_equal(last$dots$algorithm, "laplace")
  # the grid's sampler arguments do not reach the ML fit; the list's do
  expect_null(last$dots$chains)
  expect_null(last$dots$iter)
  expect_equal(last$dots$refresh, 0)
  expect_s3_class(last$prior, "brmsprior")
  expect_true(all(!nzchar(last$prior$prior)))

  ml_grid_run(
    withr::local_tempdir(),
    ml_mock = ml_mock, reps = 1L, ml_args = list(prior = "default")
  )
  expect_null(ml_mock$calls$last$prior)
})

test_that("ml is validated before any cell runs", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  run <- function(...) grid_run(dir, mock, ...)

  expect_error(run(ml = "yes"), "ml")
  expect_error(run(ml = list(1, 2)), "named list")
  expect_error(run(ml = list(model = 1)), "grid supplies")
  expect_error(run(ml = list(by = "id")), "grid supplies")
  expect_error(run(ml = list(method = "laplace")), "method")
  expect_error(run(ml = list(prior = "steep")), "prior")
  expect_error(run(ml = TRUE, levels = "population"), "subject")
  expect_identical(mock$calls$n, 0L)
})

test_that("ml with nothing varying between subjects is an error, not silence", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  expect_error(
    suppressMessages(recovery_grid(
      bmm::mixture2p(resp_error = "y"),
      grid = small_grid()[1, ],
      pars = c(kappa = log(8), thetat = stats::qlogis(0.75)),
      dir = dir, reps = 1L, seed = 1, preflight = FALSE,
      ml = list(.fitter = ml_mock_fitter()$fitter),
      .fitter = grid_mock_fitter()$fitter
    )),
    "subject level"
  )
})

test_that("ml is refused with components", {
  skip_if_not_installed("bmm")
  comps <- list(
    a = recovery_component(
      bmm::mixture2p(resp_error = "y"), c(kappa = 2, thetat = 1),
      n_trials = 5L, sds = c(kappa = 0.3), name = "a"
    )
  )
  expect_error(
    recovery_grid(
      comps,
      grid = data.frame(n_subjects = 3L),
      dir = withr::local_tempdir(), ml = TRUE,
      .fitter = grid_mock_fitter()$fitter
    ),
    "components"
  )
})

# what the review of milestone 8 found ------------------------------------

test_that("an ml route clash is refused before any cell is fitted", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  run <- function(ml) {
    suppressMessages(grid_run(dir, mock, ml = ml))
  }
  # check_ml_arg() matches `method` and `prior` one at a time, so the pair
  # clears the gate and fit_ml() refuses it once per cell --- after every
  # hierarchical fit of the design has run
  expect_error(run(list(method = "optim", prior = "default")), "prior")
  expect_identical(mock$calls$n, 0L)
  expect_error(run(list(method = "optim", draws = 50)), "draws")
  expect_error(run(list(nll = function(pars, data, model) 0)), "nll")
  expect_identical(mock$calls$n, 0L)
})

test_that("ml_cells names no file on the optim route, which does not cache", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  out <- quietly(suppressMessages(
    grid_run(dir, grid_mock_fitter(), ml = list(method = "optim"))
  ))
  cells <- attr(out, "ml_cells")
  expect_true(all(is.na(cells$file)))
})

test_that("a changed ml record keeps the hierarchical extraction", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  ml_mock <- ml_mock_fitter()
  ml_grid_run(dir, ml_mock = ml_mock)
  # ?recovery_grid says the fit files may be deleted once a grid has run
  unlink(list.files(
    dir,
    pattern = "^cell-.*-rep-[0-9]+\\.rds$", full.names = TRUE
  ))

  # only the ML record changes, so the hierarchical rows in the sidecar are
  # still current and nothing should be sampled again
  mock <- grid_mock_fitter()
  again <- ml_grid_run(
    dir, mock,
    ml_mock = ml_mock_fitter(), ml_args = list(max_abs_link = 100)
  )
  expect_identical(mock$calls$n, 0L)
  expect_setequal(again$estimator, c("bayes", "ml"))
})

test_that("a grid whose ML fit fails on some subjects says so", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  # 10 and 20 trials a subject is far too little for a stable mixture2p
  # MLE, so some subjects come back NA (decision 40 keeps their rows) and
  # the two estimators are not scored on the same people. That is the one
  # thing the user has to be told, because the metrics drop the pairs in
  # silence.
  expect_warning(
    suppressMessages(
      grid_run(dir, grid_mock_fitter(), ml = list(method = "optim"))
    ),
    "same subjects"
  )
})

test_that("an unknown ml name is refused on the optim route", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  # on the stan route an unknown name is deliberate: it reaches the Stan
  # optimiser through fit_ml()'s dots. The optim route calls no fitter, so
  # the same name reaches nothing and must not pass in silence.
  expect_error(
    suppressMessages(
      grid_run(dir, mock, ml = list(method = "optim", ci_levl = 0.8))
    ),
    "ci_levl"
  )
  expect_identical(mock$calls$n, 0L)
})
