# Alternative parameterizations (Milestone 30.1, #30; D91, D92).
#
# bmm 1.4.1 lists `criterion` and `bias` as alternatives on `sdt_yn()` and
# `sdt_rating()`: both are in `$links`, neither in `$fixed_parameters`,
# and a fit estimates exactly one. One member of each group is free: the
# one `pars` (or a formula) names, else the group's first, which is bmm's
# own default; naming both is an error. The other member appears
# nowhere. No test here compiles Stan.

skip_if_no_alternatives <- function() {
  skip_if_no_bmm_model("sdt_yn")
  skip_if(
    is.null(bmm::sdt_yn("y", "s", "n")$alternative_parameters),
    "installed bmm has no alternative parameters"
  )
}

alt_yn <- function() {
  bmm::sdt_yn(response = "hits", stimulus = "stim", n_trials = "n")
}

mock_bmm_alt <- function(formula, data, model) {
  bmm::bmm(
    formula, data, model,
    backend = "mock", mock_fit = 1, rename = FALSE, silent = 2
  )
}

# model_parameters -------------------------------------------------------------

test_that("one member of an alternative group is free: the one named", {
  skip_if_no_alternatives()
  model <- alt_yn()

  default <- model_parameters(model)
  expect_identical(default$free, c("d", "criterion"))
  expect_identical(default$alternatives, list(c("criterion", "bias")))
  expect_identical(default$fixed, "sdratio")

  expect_identical(
    model_parameters(model, named = c("d", "criterion"))$free,
    c("d", "criterion")
  )
  expect_identical(
    model_parameters(model, named = c("d", "bias"))$free, c("d", "bias")
  )
  # a task term names its parameter
  expect_identical(
    model_parameters(model, named = c("d_task1", "bias_task2"))$free,
    c("d", "bias")
  )
  # neither named: bmm's own default, the group's first member
  expect_identical(
    model_parameters(model, named = "d")$free, c("d", "criterion")
  )
})

test_that("a model without alternatives is unchanged", {
  skip_if_not_installed("bmm")
  model <- bmm::mixture2p(resp_error = "y")
  info <- model_parameters(model)
  expect_identical(info$free, c("kappa", "thetat"))
  expect_identical(info$alternatives, list())
  expect_identical(model_parameters(model, named = "kappa"), info)
})

test_that("sdt_rating's default member is criterion", {
  skip_if_no_bmm_model("sdt_rating")
  model <- bmm::sdt_rating(paste0("r", 1:4), "stim")
  skip_if(is.null(model$alternative_parameters), "no alternative parameters")
  free <- model_parameters(model)$free
  expect_true("criterion" %in% free)
  expect_false("bias" %in% free)
  expect_identical(
    names(recovery_formula(model)),
    setdiff(names(model$links), c("bias", names(model$fixed_parameters)))
  )
})

# simulate_recovery ------------------------------------------------------------

test_that("pars chooses the member, and the truth names only it", {
  skip_if_no_alternatives()
  model <- alt_yn()

  crit <- simulate_recovery(
    model, c(d = 1, criterion = 0.5),
    n_subjects = 4, n_trials = 50, sds = c(criterion = 0.2), seed = 3
  )
  expect_identical(crit$truth$population$term, c("d", "criterion"))
  expect_identical(unique(crit$truth$subjects$term), "criterion")
  expect_false("bias" %in% crit$truth$sd$term)

  bias <- simulate_recovery(
    model, c(d = 1, bias = 0),
    n_subjects = 4, n_trials = 50, seed = 3
  )
  expect_identical(bias$truth$population$term, c("d", "bias"))
  expect_named(bias$pars, c("d", "bias"))

  # with d = 1, bias 0 is criterion 0.5 (bias = criterion - d / 2), so the
  # same seed draws the same counts through either member
  same <- simulate_recovery(
    model, c(d = 1, criterion = 0.5),
    n_subjects = 4, n_trials = 50, seed = 3
  )
  expect_identical(bias$data, same$data)
})

test_that("the generator receives only the chosen member", {
  skip_if_no_alternatives()
  seen <- new.env()
  recorder <- function(pars, n_trials, model) {
    seen$pars <- pars
    data.frame(hits = c(5, 2), stim = c(1, 0), n = 10)
  }
  simulate_recovery(
    alt_yn(), c(d = 1.5, bias = 0.1),
    n_subjects = 1, n_trials = 10, generator = recorder
  )
  expect_identical(sort(names(seen$pars)), c("bias", "d", "sdratio"))
})

test_that("naming both members is an error that quotes bmm's rule", {
  skip_if_no_alternatives()
  expect_error(
    simulate_recovery(
      alt_yn(), c(d = 1, criterion = 0.3, bias = 0),
      n_subjects = 2, n_trials = 10, seed = 1
    ),
    "alternative parameterizations of the same quantity"
  )
  expect_error(
    recovery_formula(alt_yn(), pars = c(d = 1, criterion = 0, bias = 0)),
    "alternative parameterizations"
  )
})

test_that("naming neither asks for the default and names its alternative", {
  skip_if_no_alternatives()
  err <- expect_error(
    simulate_recovery(alt_yn(), c(d = 1), n_subjects = 2, n_trials = 10),
    "missing \"criterion\""
  )
  expect_match(conditionMessage(err), "bias", fixed = TRUE)
})

test_that("an sds or cors entry for the unchosen member is refused", {
  skip_if_no_alternatives()
  expect_error(
    simulate_recovery(
      alt_yn(), c(d = 1, criterion = 0.3),
      n_subjects = 2, n_trials = 10, sds = c(bias = 0.2)
    ),
    "bias.*criterion"
  )
  cors <- diag(2)
  dimnames(cors) <- rep(list(c("d", "bias")), 2)
  cors[1, 2] <- cors[2, 1] <- 0.3
  expect_error(
    simulate_recovery(
      alt_yn(), c(d = 1, criterion = 0.3),
      n_subjects = 2, n_trials = 10, sds = c(d = 0.2, criterion = 0.2),
      cors = cors
    ),
    "bias.*criterion"
  )
})

test_that("tasks expand the chosen member", {
  skip_if_no_alternatives()
  sim <- simulate_recovery(
    alt_yn(), c(d = 1, bias = 0, bias_task2 = 0.2),
    n_subjects = 2, n_trials = 10, tasks = c("1", "2"), seed = 1
  )
  expect_identical(
    sim$truth$population$term,
    c("d_task1", "d_task2", "bias_task1", "bias_task2")
  )
  formula <- recovery_formula(sim$model, task_col = "task", pars = sim$pars)
  expect_named(formula, c("d", "bias"))
})

# recovery_formula -------------------------------------------------------------

test_that("recovery_formula() writes the chosen member, and bmm accepts it", {
  skip_if_no_alternatives()
  skip_if_not_installed("brms")
  model <- alt_yn()
  expect_named(recovery_formula(model), c("d", "criterion"))
  expect_named(
    recovery_formula(model, pars = c(d = 1, bias = 0)), c("d", "bias")
  )
  for (pars in list(c(d = 1, criterion = 0.3), c(d = 1, bias = -0.2))) {
    sim <- simulate_recovery(
      model, pars,
      n_subjects = 3, n_trials = 30, sds = pars * 0 + 0.2, seed = 5
    )
    formula <- recovery_formula(model, re_cor = "all", pars = sim$pars)
    expect_named(formula, names(pars))
    fit <- mock_bmm_alt(formula, sim$data, model)
    expect_s3_class(fit, "bmmfit")
  }
})

test_that("recovery_formula() refuses a pars without names", {
  skip_if_no_alternatives()
  expect_error(recovery_formula(alt_yn(), pars = c(1, 2)), "names")
})

# recovery_grid ----------------------------------------------------------------

formula_recorder <- function() {
  inner <- grid_mock_fitter()
  seen <- new.env(parent = emptyenv())
  seen$formulas <- list()
  fitter <- function(formula, data, model, prior = NULL, ...) {
    seen$formulas[[length(seen$formulas) + 1L]] <- names(formula)
    inner$fitter(formula, data, model, prior, ...)
  }
  list(fitter = fitter, seen = seen, calls = inner$calls)
}

test_that("recovery_grid() fits each cell with the member pars chose", {
  skip_if_no_alternatives()
  mock <- formula_recorder()
  out <- suppressMessages(recovery_grid(
    alt_yn(), small_grid(),
    pars = c(d = 1, bias = 0.1), sds = c(d = 0.2, bias = 0.2),
    dir = withr::local_tempdir(), seed = 4, .fitter = mock$fitter
  ))
  expect_true(all(vapply(
    mock$seen$formulas, identical, logical(1), c("d", "bias")
  )))
  expect_setequal(out$term, c("d", "bias"))
  expect_true(all(out$covered))
})

test_that("recovery_grid() refuses a formula that names the other member", {
  skip_if_no_alternatives()
  mock <- formula_recorder()
  err <- expect_error(suppressMessages(recovery_grid(
    alt_yn(), small_grid(),
    pars = c(d = 1, criterion = 0.3),
    formula = bmm::bmf(d ~ 1 + (1 | id), bias ~ 1 + (1 | id)),
    dir = withr::local_tempdir(), seed = 4, .fitter = mock$fitter
  )))
  expect_match(conditionMessage(err), "bias", fixed = TRUE)
  expect_match(conditionMessage(err), "criterion", fixed = TRUE)
  expect_identical(mock$calls$n, 0L)

  # the same with a function `pars`, refused once it has been evaluated,
  # still before any fit
  mock <- formula_recorder()
  expect_error(suppressMessages(recovery_grid(
    alt_yn(), small_grid(),
    pars = function(row) c(d = 1, criterion = 0.3),
    formula = bmm::bmf(d ~ 1 + (1 | id), bias ~ 1 + (1 | id)),
    dir = withr::local_tempdir(), seed = 4, .fitter = mock$fitter
  )), "bias")
  expect_identical(mock$calls$n, 0L)
})

test_that("slices fit the member pars chose, as one sequential run does", {
  # D115: a function `pars` names the member only through a simulated
  # cell, which a slice holding replication 2 alone must read first
  skip_if_no_alternatives()
  skip_if_no_bmm_model("sdt_yn", density = TRUE)
  run <- function(dir, mock, cells = NULL) {
    # at 10 and 20 trials some MLEs run to the boundary, so the balance
    # warning is the right answer here; only it is muffled
    withCallingHandlers(
      suppressMessages(recovery_grid(
        alt_yn(), small_grid(),
        pars = function(row) c(d = 1, bias = 0.1),
        sds = c(d = 0.2, bias = 0.2),
        dir = dir, seed = 4, reps = 2L, cells = cells,
        ml = list(method = "optim"), .fitter = mock$fitter
      )),
      warning = function(w) {
        if (grepl("same subjects", conditionMessage(w), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }
  sidecars <- function(dir) {
    files <- sort(list.files(dir, "-est\\.rds$"))
    lapply(file.path(dir, files), function(path) {
      x <- readRDS(path)
      x[setdiff(names(x), "fit_seconds")]
    })
  }
  sequential <- withr::local_tempdir()
  sliced <- withr::local_tempdir()
  run(sequential, formula_recorder())
  run(sliced, formula_recorder(), cells = integer(0))
  slices <- grid_tasks(small_grid(), 2L, 2L)
  mocks <- list()
  for (k in 2:1) {
    mocks[[k]] <- formula_recorder()
    run(sliced, mocks[[k]], cells = slices[[k]])
  }

  for (mock in mocks) {
    expect_identical(mock$calls$n, 2L)
    expect_true(all(vapply(
      mock$seen$formulas, identical, logical(1), c("d", "bias")
    )))
  }
  expect_length(sidecars(sliced), 4L)
  expect_identical(sidecars(sliced), sidecars(sequential))
})

test_that("a slice fits a row with the member its first replication chose", {
  # D91: the realised `pars` of a row's first replication choose its
  # member for every replication. Here `pars` draws the member under each
  # cell's seed, and under seed 6 both rows chose `bias` in replication 1
  # and `criterion` in replication 2, so the slice holding replication 2
  # alone must read replication 1 before it builds a formula
  skip_if_no_alternatives()
  run <- function(dir, mock, cells = NULL) {
    # replication 2's truth names `criterion`, which a fit of `bias` does
    # not estimate; that warning is the documented consequence of D91
    withCallingHandlers(
      suppressMessages(recovery_grid(
        alt_yn(), small_grid(),
        pars = function(row) {
          if (stats::runif(1) < 0.5) {
            c(d = 1, bias = 0.1)
          } else {
            c(d = 1, criterion = 0.3)
          }
        },
        sds = c(d = 0.2),
        dir = dir, seed = 6, reps = 2L, cells = cells, .fitter = mock$fitter
      )),
      warning = function(w) {
        if (grepl("no matching estimate", conditionMessage(w), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }
  members <- function(dir) {
    files <- sprintf(
      "cell-%d-rep-%d-sim.rds", c(1L, 2L, 1L, 2L), c(1L, 1L, 2L, 2L)
    )
    vapply(file.path(dir, files), function(path) {
      setdiff(names(readRDS(path)$pars), "d")
    }, character(1), USE.NAMES = FALSE)
  }
  sidecars <- function(dir) {
    files <- sort(list.files(dir, "-est\\.rds$"))
    lapply(file.path(dir, files), function(path) {
      x <- readRDS(path)
      x[setdiff(names(x), "fit_seconds")]
    })
  }
  sequential <- withr::local_tempdir()
  sliced <- withr::local_tempdir()
  run(sequential, formula_recorder())
  run(sliced, formula_recorder(), cells = integer(0))
  # the premise: the two replications of each row chose different members
  expect_identical(members(sliced), c("bias", "bias", "criterion", "criterion"))
  slices <- grid_tasks(small_grid(), 2L, 2L)
  mocks <- list()
  for (k in 2:1) {
    mocks[[k]] <- formula_recorder()
    run(sliced, mocks[[k]], cells = slices[[k]])
  }

  # replication 2 is fitted with replication 1's member, as in one run
  for (mock in mocks) {
    expect_identical(mock$calls$n, 2L)
    expect_true(all(vapply(
      mock$seen$formulas, identical, logical(1), c("d", "bias")
    )))
  }
  expect_identical(sidecars(sliced), sidecars(sequential))
})

test_that("a grid's ML fit estimates the member pars chose", {
  skip_if_no_alternatives()
  skip_if_no_bmm_model("sdt_yn", density = TRUE)
  out <- suppressMessages(recovery_grid(
    alt_yn(), data.frame(n_subjects = 3L, n_trials = 40L),
    pars = c(d = 1, bias = 0.1), sds = c(d = 0.2, bias = 0.2),
    dir = withr::local_tempdir(), seed = 4,
    ml = list(method = "optim"), .fitter = grid_mock_fitter()$fitter
  ))
  # subject rows are scored against the varying terms only, so an ML fit
  # of `criterion` would leave the `bias` rows missing rather than wrong
  ml_rows <- out[out$estimator == "ml", ]
  expect_setequal(ml_rows$term, c("d", "bias"))
  expect_equal(sum(ml_rows$term == "bias"), 3L)
  expect_false(any(grepl("criterion", out$term)))
})

# fit_ml -----------------------------------------------------------------------

test_that("fit_ml() follows its formula or its start values", {
  skip_if_no_alternatives()
  skip_if_no_bmm_model("sdt_yn", density = TRUE)
  model <- alt_yn()
  sim <- simulate_recovery(
    model, c(d = 1, bias = 0.1),
    n_subjects = 3, n_trials = 200, sds = c(d = 0.2), seed = 9
  )
  optim_default <- fit_ml(model, sim$data, method = "optim")
  expect_setequal(optim_default$term, c("d", "criterion"))
  optim_bias <- fit_ml(
    model, sim$data,
    method = "optim", start = c(d = 0, bias = 0)
  )
  expect_setequal(optim_bias$term, c("d", "bias"))
  # criterion = bias + d / 2, subject by subject
  pick <- function(x, term) {
    rows <- x[x$term == term, ]
    stats::setNames(rows$estimate, as.character(rows$id))[c("1", "2", "3")]
  }
  expect_equal(
    pick(optim_default, "criterion"),
    pick(optim_bias, "bias") + pick(optim_bias, "d") / 2,
    tolerance = 1e-3
  )

  expect_error(
    fit_ml(
      model, sim$data,
      method = "optim", start = c(d = 0, criterion = 0, bias = 0)
    ),
    "alternative parameterizations"
  )

  ml_mock <- ml_mock_fitter()
  stan_bias <- fit_ml(
    model, sim$data,
    formula = bmm::bmf(d ~ 0 + id, bias ~ 0 + id), .fitter = ml_mock$fitter
  )
  expect_setequal(stan_bias$term, c("d", "bias"))
})

test_that("fit_ml()'s default formula and flat prior name the default member", {
  skip_if_no_alternatives()
  skip_if_not_installed("brms")
  model <- alt_yn()
  expect_named(no_pooling_formula(model), c("d", "criterion"))
  expect_named(
    no_pooling_formula(model, named = c("d", "bias")), c("d", "bias")
  )
  prior <- flat_ml_prior(model, named = c("d", "bias"))
  expect_setequal(prior$nlpar, c("d", "bias"))
})

# components -------------------------------------------------------------------

test_that("a component's formula must name the member its pars chose", {
  skip_if_no_alternatives()
  expect_error(
    recovery_component(
      alt_yn(),
      pars = c(d = 1, criterion = 0.3), n_trials = 20,
      formula = bmm::bmf(d ~ 1 + (1 | id), bias ~ 1 + (1 | id)),
      name = "yn"
    ),
    "bias"
  )
  comp <- recovery_component(
    alt_yn(),
    pars = c(d = 1, bias = 0.1), n_trials = 20, sds = c(bias = 0.2),
    name = "yn"
  )
  set <- simulate_components(list(comp), n_subjects = 3, seed = 1)
  expect_identical(set$components$yn$truth$population$term, c("d", "bias"))
  expect_named(
    component_formula(set$specs$yn, set$components$yn), c("d", "bias")
  )
})

# sbc --------------------------------------------------------------------------

test_that("sbc()'s draw patterns follow the member its formula names", {
  skip_if_no_alternatives()
  patterns <- sbc_variables(
    alt_yn(), "id", c("population", "sd"),
    named = c("d", "bias")
  )
  expect_true(any(grepl("bias", patterns)))
  expect_false(any(grepl("criterion", patterns)))
  default <- sbc_variables(alt_yn(), "id", "population")
  expect_true(any(grepl("criterion", default)))
  expect_false(any(grepl("bias", default)))
})

test_that("sbc() simulates sdt_yn under a bias formula", {
  skip_if_not_installed("SBC")
  skip_on_cran()
  skip_if_no_alternatives()
  data <- data.frame(
    id = rep(c("a", "b", "c"), each = 2), stim = rep(c(1, 0), 3),
    n = 50, hits = 0
  )
  draws <- posterior::as_draws_matrix(cbind(
    b_d_Intercept = 1 + 0.001 * seq_len(60),
    b_bias_Intercept = 0.1 + 0.001 * seq_len(60),
    sd_id__d_Intercept = 0.1, sd_id__bias_Intercept = 0.1
  ))
  mock <- sbc_mock_fitter(draws = draws)
  formula <- bmm::bmf(d ~ 1 + (1 | id), bias ~ 1 + (1 | id))
  suppressMessages(withCallingHandlers(
    sbc(
      alt_yn(), formula, data,
      n_sims = 2L, thin_ranks = 1, cores_per_fit = 1, level = "population",
      .fitter = mock$fitter
    ),
    bmmtools_sbc_diagnostics = function(w) invokeRestart("muffleWarning")
  ))
  calls <- sbc_dataset_calls(mock$calls)
  expect_length(calls, 2L)
  # rsdt_yn() needs exactly one of criterion and bias, so a generator
  # handed the wrong member would have stopped before any data set
  for (call in calls) {
    expect_identical(call$data$stim, data$stim)
    expect_true(all(call$data$hits <= call$data$n))
  }
  expect_error(
    sbc(
      alt_yn(),
      bmm::bmf(d ~ 1 + (1 | id), criterion ~ 1, bias ~ 1 + (1 | id)),
      data,
      n_sims = 2L, .fitter = mock$fitter
    ),
    "alternative parameterizations"
  )
})

# review follow-ups ------------------------------------------------------------

test_that("a missing term of the chosen member does not ask for the other", {
  skip_if_no_alternatives()
  err <- expect_error(
    simulate_recovery(
      alt_yn(), c(d = 1, bias_task1 = 0),
      n_subjects = 2, n_trials = 10, tasks = c("1", "2")
    ),
    "bias_task2"
  )
  expect_no_match(conditionMessage(err), "has the alternative")
})

test_that("a function pars chooses each cell's default formula member", {
  skip_if_no_alternatives()
  mock <- formula_recorder()
  out <- suppressMessages(recovery_grid(
    alt_yn(), small_grid(),
    pars = function(row) c(d = 1, bias = stats::rnorm(1, 0, 0.1)),
    sds = c(d = 0.2, bias = 0.2), reps = 2,
    dir = withr::local_tempdir(), seed = 4, .fitter = mock$fitter
  ))
  # the preflight and four cells, each fitted for bias
  expect_length(mock$seen$formulas, 5L)
  expect_true(all(vapply(
    mock$seen$formulas, identical, logical(1), c("d", "bias")
  )))
  expect_setequal(out$term, c("d", "bias"))
})

test_that("a grid's ml formula or start must estimate the member pars chose", {
  skip_if_no_alternatives()
  mock <- formula_recorder()
  expect_error(suppressMessages(recovery_grid(
    alt_yn(), small_grid(),
    pars = c(d = 1, bias = 0.1),
    ml = list(method = "optim", start = c(d = 0, criterion = 0)),
    dir = withr::local_tempdir(), seed = 4, .fitter = mock$fitter
  )), "start")
  expect_error(suppressMessages(recovery_grid(
    alt_yn(), small_grid(),
    pars = c(d = 1, bias = 0.1),
    ml = list(formula = bmm::bmf(d ~ 0 + id, criterion ~ 0 + id)),
    dir = withr::local_tempdir(), seed = 4, .fitter = mock$fitter
  )), "criterion")
  expect_identical(mock$calls$n, 0L)
})

test_that("a grid column for the member pars did not choose is refused", {
  skip_if_no_alternatives()
  mock <- formula_recorder()
  for (grid in list(
    data.frame(n_subjects = 3L, n_trials = 20L, bias = c(0, 1)),
    data.frame(n_subjects = 3L, n_trials = 20L, sd_bias = c(0.1, 0.2))
  )) {
    expect_error(suppressMessages(recovery_grid(
      alt_yn(), grid,
      pars = c(d = 1, criterion = 0.3),
      dir = withr::local_tempdir(), seed = 4, .fitter = mock$fitter
    )), "bias")
  }
  expect_identical(mock$calls$n, 0L)
})

test_that("a grid of components fits each component's chosen member", {
  skip_if_no_alternatives()
  mock <- formula_recorder()
  comp <- recovery_component(
    alt_yn(),
    pars = c(d = 1, bias = 0.1), n_trials = 20, sds = c(bias = 0.2),
    name = "yn"
  )
  out <- suppressMessages(recovery_grid(
    list(comp), data.frame(n_subjects = c(3L, 4L)),
    dir = withr::local_tempdir(), reps = 1, seed = 4,
    .fitter = mock$fitter
  ))
  expect_gt(length(mock$seen$formulas), 0L)
  expect_true(all(vapply(
    mock$seen$formulas, identical, logical(1), c("d", "bias")
  )))
  expect_false(any(grepl("criterion", out$term)))
})

test_that("a missing name never removes a parameter from the free set", {
  skip_if_no_alternatives()
  expect_identical(
    model_parameters(alt_yn(), named = c("d", NA))$free, c("d", "criterion")
  )
})
