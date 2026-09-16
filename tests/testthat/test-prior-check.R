# Tests for prior_check() and plot_prior_check(), written against the
# spec local/dev/spec-milestone-4-prior-check.md in full.
#
# Three routes, and the split is deliberate. The pure functions are
# tested on hand-built matrices with values that can be written down.
# The wiring is tested with the mock fitter of Milestone 2 plus the
# prior_predict_draws() method helper-prior-check.R registers, because
# the mock *backend* has no draws (measured 2026-09-08). The prediction
# step itself is tested on the committed fixture, which is a posterior
# fit, not a prior-only one: it exercises the same
# fit -> posterior_predict() -> matrix path, and nothing here claims it
# came from a prior.

# 1. the default summary ------------------------------------------------

test_that("default_prior_summary computes the statistics of the spec", {
  out <- default_prior_summary(c(0, 10), fake_model())(known_yrep(), NULL)

  expect_s3_class(out, "tbl_df")
  expect_true(all(c("statistic", "value") %in% names(out)))
  expect_equal(
    out$statistic,
    c("floor_rate", "ceiling_rate", "q50", "q90", "q95")
  )
  expect_equal(out$value[out$statistic == "floor_rate"], 1 / 12)
  expect_equal(out$value[out$statistic == "ceiling_rate"], 2 / 12)
  expect_equal(out$value[out$statistic == "q50"], 5.5)
  expect_equal(out$value[out$statistic == "q90"], 9.9)
  expect_equal(out$value[out$statistic == "q95"], 10)
  expect_equal(unique(out$response), "response")
})

test_that("a vector ceiling is compared column-wise", {
  ceilings <- c(1, 2, 3, 4)
  at_ceiling <- matrix(rep(ceilings, each = 3L), nrow = 3L)
  summarise <- default_prior_summary(list(floor = 0, ceiling = ceilings),
    model = fake_model()
  )

  full <- summarise(at_ceiling, NULL)
  expect_equal(full$value[full$statistic == "ceiling_rate"], 1)

  below <- at_ceiling
  below[, 4L] <- 0
  part <- summarise(below, NULL)
  expect_equal(
    part$value[part$statistic == "ceiling_rate"],
    1 - 1 / ncol(at_ceiling)
  )
})

test_that("an unknown range gives NA rates and finite quantiles", {
  out <- default_prior_summary(c(NA, NA), fake_model())(known_yrep(), NULL)

  rates <- out$value[out$statistic %in% c("floor_rate", "ceiling_rate")]
  expect_true(all(is.na(rates)))
  expect_true(all(is.finite(out$value[startsWith(out$statistic, "q")])))
})

# 2. the response range -------------------------------------------------

test_that("response_range knows the count models by their own columns", {
  skip_if_not_installed("bmm")
  model <- bmm::sdt_mafc(response = "k", n_trials = "n", m = 4L)
  data <- data.frame(k = c(1, 2), n = c(10, 20))

  out <- response_range(model, data)

  expect_equal(out$floor, 0)
  expect_equal(out$ceiling, c(10, 20))
})

test_that("response_range knows the circular models", {
  skip_if_not_installed("bmm")
  for (model in list(
    bmm::mixture2p(resp_error = "y"),
    bmm::sdm(resp_error = "y")
  )) {
    out <- response_range(model, data.frame(y = 0))
    expect_equal(out$floor, -pi)
    expect_equal(out$ceiling, pi)
  }
})

test_that("response_range says NA where it has no boundary to give", {
  skip_if_not_installed("bmm")
  model <- bmm::ddm(rt = "rt", response = "resp")

  expect_message(out <- response_range(model, data.frame(rt = 1, resp = 1)))
  expect_true(is.na(out$floor))
  expect_true(is.na(out$ceiling))

  expect_message(unknown <- response_range(fake_model(), fake_data()))
  expect_true(is.na(unknown$ceiling))
})

test_that("a count model without its n_trials column degrades to NA", {
  skip_if_not_installed("bmm")
  model <- bmm::sdt_mafc(response = "k", n_trials = "n", m = 4L)

  expect_message(out <- response_range(model, data.frame(k = 1)))
  expect_true(is.na(out$ceiling))
})

# 3. shaping the draws --------------------------------------------------

test_that("yrep_list names a matrix after the model's response", {
  skip_if_not_installed("bmm")
  model <- bmm::mixture2p(resp_error = "y")

  out <- yrep_list(known_yrep(), model)

  expect_named(out, "y")
  expect_identical(out$y, known_yrep())
})

test_that("yrep_list splits an array and falls back on names", {
  named <- array(
    seq_len(12), dim = c(2L, 3L, 2L),
    dimnames = list(NULL, NULL, c("rt", "resp"))
  )
  expect_named(yrep_list(named, fake_model()), c("rt", "resp"))
  expect_equal(dim(yrep_list(named, fake_model())$rt), c(2L, 3L))

  bare <- array(seq_len(12), dim = c(2L, 3L, 2L))
  expect_named(yrep_list(bare, fake_model()), c("response_1", "response_2"))
})

test_that("yrep_list refuses what it cannot shape", {
  expect_error(yrep_list(list(1, 2), fake_model()), "list")
})

# 4. prior sets ---------------------------------------------------------

test_that("prior_sets names the sets", {
  expect_identical(prior_sets(NULL), list(default = NULL))
  expect_named(prior_sets(fake_prior()), "user")
  expect_named(
    prior_sets(list(loose = fake_prior(), tight = fake_prior())),
    c("loose", "tight")
  )
  expect_named(prior_sets(list(default = NULL, tight = fake_prior())),
    c("default", "tight")
  )
})

test_that("a list of prior sets must be named", {
  expect_error(
    prior_sets(list(fake_prior(), fake_prior())),
    "name"
  )
})

# 5. wiring -------------------------------------------------------------

test_that("one prior set fits once from the prior and records the seed", {
  mock <- prior_mock_fitter()

  out <- suppressMessages(mock_prior_check(mock, seed = 7))

  expect_s3_class(out, "bmmtools_prior_check")
  expect_named(out, c("prior", "response", "statistic", "value"))
  expect_identical(unique(out$prior), "default")
  expect_identical(mock$calls$n, 1L)
  expect_identical(mock$calls$log[[1L]]$sample_prior, "only")
  expect_identical(mock$calls$log[[1L]]$seed, 7)
  expect_equal(attr(out, "seed"), 7)
  expect_equal(nrow(out), 5L)
})

test_that("no seed is recorded as NA and none is passed on", {
  mock <- prior_mock_fitter()

  out <- suppressMessages(mock_prior_check(mock))

  expect_true(is.na(attr(out, "seed")))
  expect_null(mock$calls$log[[1L]]$seed)
})

test_that("two prior sets are fitted, labelled and compared", {
  mock <- prior_mock_fitter()
  p <- fake_prior("normal(0, 5)")

  out <- suppressMessages(mock_prior_check(
    mock,
    prior = list(default = NULL, tight = p)
  ))

  expect_identical(mock$calls$n, 2L)
  expect_null(mock$calls$log[[1L]]$prior)
  expect_identical(mock$calls$log[[2L]]$prior, p)
  expect_identical(unique(out$prior), c("default", "tight"))

  table <- summary(out)
  expect_s3_class(table, "bmmtools_prior_check_summary")
  expect_true(all(c("default", "tight", "difference") %in% names(table)))
  expect_equal(table$difference, table$tight - table$default)
})

test_that("three prior sets get no difference column", {
  mock <- prior_mock_fitter()

  out <- suppressMessages(mock_prior_check(mock, prior = list(
    a = fake_prior(), b = fake_prior(), c = fake_prior()
  )))

  expect_false("difference" %in% names(summary(out)))
  expect_identical(mock$calls$n, 3L)
})

test_that("a file caches the fits and a second call refits nothing", {
  dir <- withr::local_tempdir()
  mock <- prior_mock_fitter()
  args <- list(prior = list(default = NULL, tight = fake_prior()))

  first <- suppressMessages(do.call(
    mock_prior_check, c(list(mock), args, file = file.path(dir, "pc"))
  ))
  expect_identical(mock$calls$n, 2L)
  expect_true(all(file.exists(file.path(
    dir, c("pc-default.rds", "pc-default.key", "pc-tight.rds", "pc-tight.key")
  ))))

  second <- suppressMessages(do.call(
    mock_prior_check, c(list(mock), args, file = file.path(dir, "pc"))
  ))
  expect_identical(mock$calls$n, 2L)
  expect_equal(second$value, first$value)
})

test_that("file = NULL leaves nothing behind in the working directory", {
  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  mock <- prior_mock_fitter()

  suppressMessages(mock_prior_check(mock))

  expect_identical(list.files(dir), character(0))
})

test_that("n_draws reaches the prediction step", {
  mock <- prior_mock_fitter()

  suppressMessages(mock_prior_check(mock, n_draws = 3L))

  expect_identical(mock$calls$ndraws, 3L)
})

test_that("a user summary is called with the draws and the data", {
  mock <- prior_mock_fitter()
  seen <- NULL
  user_summary <- function(yrep, data) {
    seen <<- list(yrep = yrep, data = data)
    tibble::tibble(statistic = "mean", value = mean(yrep))
  }

  out <- suppressMessages(mock_prior_check(mock, summary = user_summary))

  expect_equal(dim(seen$yrep), c(8L, 5L))
  expect_equal(seen$data, fake_data())
  expect_identical(out$statistic, "mean")
  expect_equal(out$value, mean(seen$yrep))
  # a summary that gives no response gets the model's single response name
  expect_identical(out$response, "response")
})

test_that("a user summary keeps a response column it supplies", {
  mock <- prior_mock_fitter()
  user_summary <- function(yrep, data) {
    tibble::tibble(
      response = c("a", "b"), statistic = "mean", value = c(1, 2)
    )
  }

  out <- suppressMessages(mock_prior_check(mock, summary = user_summary))

  expect_identical(out$response, c("a", "b"))
})

test_that("a failing fitter names the set and keeps the cause", {
  mock <- prior_mock_fitter(fail_on = function(n) n == 2L)

  expect_error(
    suppressMessages(mock_prior_check(mock, prior = list(
      ok = fake_prior(), bad = fake_prior()
    ))),
    "bad"
  )
  err <- tryCatch(
    suppressMessages(mock_prior_check(
      prior_mock_fitter(fail_on = function(n) TRUE),
      prior = list(bad = fake_prior())
    )),
    error = function(e) e
  )
  expect_match(conditionMessage(err$parent), "mock sampler failure")
})

# 6. validation ---------------------------------------------------------

test_that("prior_check validates its arguments", {
  mock <- prior_mock_fitter()
  form <- fake_bmmformula(a = a ~ 1)

  expect_error(
    prior_check("nope", form, fake_data(), .fitter = mock$fitter),
    "bmmodel"
  )
  expect_error(
    prior_check(fake_model(), form, "nope", .fitter = mock$fitter),
    "data frame"
  )
  expect_error(
    prior_check(fake_model(), form, fake_data()[0, ], .fitter = mock$fitter),
    "row"
  )
  expect_error(mock_prior_check(mock, n_draws = 0), "n_draws")
  expect_error(mock_prior_check(mock, n_draws = c(1, 2)), "n_draws")
  expect_error(mock_prior_check(mock, range = 1), "range")
  expect_error(mock_prior_check(mock, summary = "nope"), "summary")
  expect_error(mock_prior_check(mock, summary = function(x) x), "summary")
  expect_error(mock_prior_check(mock, seed = c(1, 2)), "seed")
  expect_error(mock_prior_check(mock, sample_prior = "yes"), "sample_prior")
})

test_that("a summary that breaks its contract is named", {
  mock <- prior_mock_fitter()

  expect_error(
    suppressMessages(mock_prior_check(mock, summary = function(yrep, data) 1)),
    "data frame"
  )
  expect_error(
    suppressMessages(mock_prior_check(
      mock,
      summary = function(yrep, data) tibble::tibble(stat = "x", value = 1)
    )),
    "statistic"
  )
})

test_that("n_draws above what the fit holds is reported, not refused", {
  mock <- prior_mock_fitter()

  # fit_cached() speaks first, so the messages are collected rather than
  # matched one at a time
  messages <- testthat::capture_messages(
    out <- mock_prior_check(mock, n_draws = 500L)
  )

  expect_true(any(grepl("500", messages)))
  expect_s3_class(out, "bmmtools_prior_check")
})

test_that("a flat population-level slope is warned about", {
  skip_if_not_installed("bmm")
  mock <- prior_mock_fitter()
  model <- bmm::mixture2p(resp_error = "y")
  data <- data.frame(
    id = factor(rep(1:3, each = 4L)),
    cond = factor(rep(c("a", "b"), 6L)),
    y = stats::runif(12L, -pi, pi)
  )

  expect_warning(
    suppressMessages(prior_check(
      model, bmm::bmf(kappa ~ 1 + cond + (1 | id), thetat ~ 1 + (1 | id)),
      data,
      .fitter = mock$fitter
    )),
    "condb"
  )
  expect_no_warning(
    suppressMessages(prior_check(
      model, bmm::bmf(kappa ~ 1 + (1 | id), thetat ~ 1 + (1 | id)),
      data,
      .fitter = mock$fitter
    ))
  )
})

# 7. the class ----------------------------------------------------------

test_that("dplyr verbs keep the class only while the contract holds", {
  mock <- prior_mock_fitter()
  out <- suppressMessages(mock_prior_check(mock))

  kept <- dplyr::filter(out, .data$statistic == "q50")
  expect_s3_class(kept, "bmmtools_prior_check")

  dropped <- dplyr::select(out, "value")
  expect_false(inherits(dropped, "bmmtools_prior_check"))
  expect_s3_class(dropped, "tbl_df")
})

test_that("the recovery class still round-trips through dplyr", {
  # demote_if_incomplete() gained an argument; Milestone 1 must not move
  recovery <- recover(
    fake_estimates(c("kappa", "thetat"), estimate = c(1, 2)),
    fake_truth(c("kappa", "thetat"), true_value = c(1, 2)),
    scale = "link"
  )
  expect_s3_class(
    dplyr::filter(recovery, .data$term == "kappa"), "bmmtools_recovery"
  )
  expect_false(
    inherits(dplyr::select(recovery, "term"), "bmmtools_recovery")
  )
})

test_that("print names the model and the sets", {
  mock <- prior_mock_fitter()
  out <- suppressMessages(mock_prior_check(mock, prior = list(
    loose = fake_prior(), tight = fake_prior()
  )))

  text <- paste(utils::capture.output(print(out)), collapse = "\n")

  expect_match(text, "bmmtools_prior_check")
  expect_match(text, "loose")
  expect_match(text, "tight")
  # the fake model has no known response range, so the reason is printed
  expect_match(text, "range")
})

# 8. the prediction step, on the fixture ---------------------------------

test_that("prior_predict_draws reaches brms::posterior_predict", {
  skip_if_not_installed("brms")
  fit <- readRDS(test_path("fixtures", "mixture2p-fit.rds"))

  yrep <- prior_predict_draws(fit, ndraws = 20L, seed = 1)

  expect_true(is.matrix(yrep))
  expect_equal(dim(yrep), c(20L, 240L))
  expect_type(yrep, "double")
})

test_that("the same seed gives the same subsample", {
  skip_if_not_installed("brms")
  fit <- readRDS(test_path("fixtures", "mixture2p-fit.rds"))

  expect_equal(
    prior_predict_draws(fit, ndraws = 5L, seed = 3),
    prior_predict_draws(fit, ndraws = 5L, seed = 3)
  )
  expect_false(identical(
    prior_predict_draws(fit, ndraws = 5L, seed = 3),
    prior_predict_draws(fit, ndraws = 5L, seed = 4)
  ))
})

test_that("the default summary reads a real predictive matrix", {
  skip_if_not_installed("brms")
  skip_if_not_installed("bmm")
  fit <- readRDS(test_path("fixtures", "mixture2p-fit.rds"))
  yrep <- prior_predict_draws(fit, ndraws = 20L, seed = 1)
  model <- bmm::mixture2p(resp_error = "y")

  out <- default_prior_summary(c(-pi, pi), model)(yrep, fit$data)

  expect_equal(out$value[out$statistic == "floor_rate"], 0)
  expect_equal(out$value[out$statistic == "ceiling_rate"], 0)
  expect_equal(
    out$value[out$statistic %in% c("q50", "q90", "q95")],
    unname(stats::quantile(yrep, c(0.5, 0.9, 0.95)))
  )
  expect_identical(unique(out$response), "y")
})

# 9. the plot -----------------------------------------------------------

test_that("plot_prior_check returns a ggplot for every type", {
  skip_if_not_installed("ggplot2")
  mock <- prior_mock_fitter()
  out <- suppressMessages(mock_prior_check(mock))

  for (type in c("density", "histogram", "statistic")) {
    expect_s3_class(plot_prior_check(out, type = type), "ggplot")
  }
})

test_that("plot_prior_check validates what it is given", {
  skip_if_not_installed("ggplot2")
  mock <- prior_mock_fitter()
  out <- suppressMessages(mock_prior_check(mock))

  expect_error(plot_prior_check(fake_data()), "bmmtools_prior_check")
  expect_error(plot_prior_check(out, facet_by = "nope"), "nope")
  expect_error(plot_prior_check(out, type = "spaghetti"), "type")
})

test_that("the observed overlay is one layer that can be turned off", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("bmm")
  mock <- prior_mock_fitter(yrep = matrix(stats::runif(60, -pi, pi), nrow = 6L))
  data <- data.frame(
    id = factor(rep(1:2, each = 5L)), y = stats::runif(10L, -pi, pi)
  )
  out <- suppressMessages(prior_check(
    bmm::mixture2p(resp_error = "y"),
    bmm::bmf(kappa ~ 1, thetat ~ 1),
    data,
    .fitter = mock$fitter
  ))

  with_obs <- plot_prior_check(out, observed = TRUE)
  without <- plot_prior_check(out, observed = FALSE)

  expect_equal(length(with_obs$layers), length(without$layers) + 1L)
})

test_that("draws above what the object holds are capped", {
  skip_if_not_installed("ggplot2")
  mock <- prior_mock_fitter()
  out <- suppressMessages(mock_prior_check(mock))

  expect_s3_class(plot_prior_check(out, draws = 1000L), "ggplot")
})

test_that("two prior sets give a colour scale with two levels", {
  skip_if_not_installed("ggplot2")
  mock <- prior_mock_fitter()
  out <- suppressMessages(mock_prior_check(mock, prior = list(
    loose = fake_prior(), tight = fake_prior()
  )))

  p <- plot_prior_check(out, type = "density")
  built <- ggplot2::ggplot_build(p)

  expect_length(unique(built$data[[1L]]$colour), 2L)
})

# check_improper_priors() ------------------------------------------------

#' A formula with cell means on one parameter, so brms writes coefficient
#' rows that a class-level prior can cover.
task_prior_data <- function() {
  data.frame(
    id = rep(1:3, each = 4L),
    y = rep(c(-0.5, 0.1, 0.4, -0.2), times = 3L),
    task = rep(c("task1", "task2"), length.out = 12L)
  )
}

test_that("a class-level prior covers its coefficient rows", {
  skip_if_not_installed("bmm")

  expect_no_warning(
    check_improper_priors(
      bmm::bmf(kappa ~ 0 + task + (1 | id), thetat ~ 1 + (1 | id)),
      task_prior_data(),
      bmm::mixture2p(resp_error = "y"),
      list(covered = brms::set_prior(
        "normal(0, 1)", class = "b", nlpar = "kappa"
      ))
    )
  )
})

test_that("a coefficient with no prior at either level is still flat", {
  skip_if_not_installed("bmm")

  expect_warning(
    check_improper_priors(
      bmm::bmf(kappa ~ 0 + task + (1 | id), thetat ~ 1 + (1 | id)),
      task_prior_data(),
      bmm::mixture2p(resp_error = "y"),
      list(none = NULL)
    ),
    "tasktask1"
  )
})

test_that("a class-level prior covers its own parameter only", {
  skip_if_not_installed("bmm")

  expect_warning(
    check_improper_priors(
      bmm::bmf(kappa ~ 0 + task + (1 | id), thetat ~ 0 + task + (1 | id)),
      task_prior_data(),
      bmm::mixture2p(resp_error = "y"),
      list(half = brms::set_prior(
        "normal(0, 1)", class = "b", nlpar = "kappa"
      ))
    ),
    "thetat"
  )
})

test_that("per-coefficient priors are covered as before", {
  skip_if_not_installed("bmm")

  expect_no_warning(
    check_improper_priors(
      bmm::bmf(kappa ~ 0 + task + (1 | id), thetat ~ 1 + (1 | id)),
      task_prior_data(),
      bmm::mixture2p(resp_error = "y"),
      list(per_coef = brms::set_prior(
        "normal(0, 1)", class = "b", coef = "tasktask1", nlpar = "kappa"
      ) + brms::set_prior(
        "normal(0, 1)", class = "b", coef = "tasktask2", nlpar = "kappa"
      ))
    )
  )
})
