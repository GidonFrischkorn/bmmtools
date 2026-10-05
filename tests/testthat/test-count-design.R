# Count models with a trial design (Milestone 29.1, #26).
#
# One row of a count model holds many trials. A `trial_design` gives each
# row its own number of trials (and, for `sdt_yn`, its stimulus class), so
# the simulated data have the rows the design describes. Without a design
# the adapters behave as in 0.2.0.

skip_if_no_sdt <- function() {
  skip_if_not_installed("bmm")
  skip_if_not(
    exists("sdt_yn", envir = asNamespace("bmm"), inherits = FALSE),
    "installed bmm has no signal-detection models"
  )
}

yn_model <- function() {
  bmm::sdt_yn(response = "hits", stimulus = "stim", n_trials = "n")
}

mafc_model <- function() {
  bmm::sdt_mafc(response = "k", n_trials = "n", m = 4)
}

ez_model <- function() {
  bmm::ezdm(mean_rt = "mrt", var_rt = "vrt", n_upper = "nu", n_trials = "n")
}

m3_count_model <- function() {
  bmm::m3(
    resp_cats = c("corr", "other", "npl"), num_options = c(1, 4, 5),
    choice_rule = "simple", version = "ss"
  )
}

# simulate_recovery ------------------------------------------------------------

test_that("count adapters name their optional design columns", {
  skip_if_not_installed("bmm")
  expect_identical(count_design_columns(ez_model()), "n")
  expect_identical(count_design_columns(m3_count_model()), "n_trials")
  expect_identical(count_design_columns(bmm::mixture2p("y")), character())
  # a model no adapter knows
  expect_identical(
    count_design_columns(structure(list(), class = "bmmodel")), character()
  )
  # an m3 with option columns is a per-trial model
  m3_cols <- bmm::m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c("n_corr", "n_other", "n_npl"),
    choice_rule = "simple", version = "ss"
  )
  expect_identical(count_design_columns(m3_cols), character())
  skip_if_no_sdt()
  expect_identical(count_design_columns(yn_model()), c("stim", "n"))
  expect_identical(count_design_columns(mafc_model()), "n")
})

test_that("sdt_yn takes each row's stimulus and trials from the design", {
  skip_if_no_sdt()
  design <- data.frame(stim = c(0, 1), n = c(50, 100))
  sim <- simulate_recovery(
    yn_model(), c(d = 1, criterion = 0),
    n_subjects = 3, n_trials = 2, trial_design = design, seed = 1
  )
  expect_named(sim$data, c("id", "stim", "n", "hits"))
  expect_identical(sim$data$stim, rep(c(0, 1), 3))
  expect_identical(sim$data$n, rep(c(50, 100), 3))
  expect_true(all(sim$data$hits >= 0 & sim$data$hits <= sim$data$n))
  expect_identical(sim$trial_design_columns, c("stim", "n"))
})

test_that("sdt_mafc takes a per-subject number of trials from the design", {
  skip_if_no_sdt()
  design <- data.frame(id = 1:3, n = c(40, 80, 120))
  sim <- simulate_recovery(
    mafc_model(), c(d = 1),
    n_subjects = 3, n_trials = 1, trial_design = design, seed = 1
  )
  expect_identical(sim$data$n, c(40, 80, 120))
  expect_true(all(sim$data$k <= sim$data$n))
})

test_that("ezdm simulates each row from its own number of trials", {
  skip_if_not_installed("bmm")
  design <- data.frame(n = c(100, 3))
  sim <- simulate_recovery(
    ez_model(), c(drift = 1, bound = log(1.2), ndt = log(0.3)),
    n_subjects = 2, n_trials = 2, trial_design = design, seed = 1
  )
  expect_named(sim$data, c("id", "n", "mrt", "vrt", "nu"))
  expect_identical(sim$data$n, c(100, 3, 100, 3))
  expect_true(all(sim$data$nu <= sim$data$n))
})

test_that("m3 with numbers of options takes each row's total from n_trials", {
  skip_if_not_installed("bmm")
  design <- data.frame(n_trials = c(100, 60))
  sim <- simulate_recovery(
    m3_count_model(), c(c = log(3), a = log(0.5)),
    n_subjects = 2, n_trials = 2, trial_design = design, seed = 1
  )
  expect_named(sim$data, c("id", "n_trials", "corr", "other", "npl"))
  totals <- rowSums(as.data.frame(sim$data)[c("corr", "other", "npl")])
  expect_equal(totals, c(100, 60, 100, 60))
})

test_that("a count design from a function is read per call", {
  skip_if_not_installed("bmm")
  sim <- simulate_recovery(
    ez_model(), c(drift = 1, bound = log(1.2), ndt = log(0.3)),
    n_subjects = 2, n_trials = 1,
    trial_design = function(n_trials) data.frame(n = 50),
    seed = 1
  )
  expect_identical(sim$data$n, c(50, 50))
})

test_that("without a design the count adapters are unchanged", {
  skip_if_not_installed("bmm")
  sim <- simulate_recovery(
    ez_model(), c(drift = 1, bound = log(1.2), ndt = log(0.3)),
    n_subjects = 2, n_trials = 40, seed = 1
  )
  expect_named(sim$data, c("id", "mrt", "vrt", "nu", "n"))
  expect_equal(sim$data$n, c(40, 40))
  expect_null(sim$trial_design)
})

test_that("a count design without its count column is an error", {
  skip_if_not_installed("bmm")
  expect_error(
    simulate_recovery(
      ez_model(), c(drift = 1, bound = log(1.2), ndt = log(0.3)),
      n_subjects = 2, n_trials = 1, trial_design = data.frame(other = 1)
    ),
    "lacks the column"
  )
})

test_that("trial counts that are not whole and positive are refused by row", {
  skip_if_not_installed("bmm")
  pars <- c(drift = 1, bound = log(1.2), ndt = log(0.3))
  expect_error(
    simulate_recovery(
      ez_model(), pars,
      n_subjects = 1, n_trials = 3, trial_design = data.frame(n = c(10, 2.5, 0))
    ),
    "2 and 3"
  )
  expect_error(
    simulate_recovery(
      m3_count_model(), c(c = 1, a = 0),
      n_subjects = 1, n_trials = 1, trial_design = data.frame(n_trials = NA)
    ),
    "whole number"
  )
})

test_that("ezdm refuses fewer than 3 trials in a row, before any draw", {
  skip_if_not_installed("bmm")
  expect_error(
    simulate_recovery(
      ez_model(), c(drift = 1, bound = log(1.2), ndt = log(0.3)),
      n_subjects = 1, n_trials = 2, trial_design = data.frame(n = c(10, 2))
    ),
    "at least 3"
  )
})

test_that("sdt_yn refuses a stimulus that is not 0 or 1", {
  skip_if_no_sdt()
  expect_error(
    simulate_recovery(
      yn_model(), c(d = 1, criterion = 0),
      n_subjects = 1, n_trials = 2,
      trial_design = data.frame(stim = c(1, 2), n = 10)
    ),
    "0 or 1"
  )
})

test_that("a per-trial adapter still refuses a design", {
  skip_if_not_installed("bmm")
  expect_error(
    simulate_recovery(
      bmm::mixture2p("y"), c(kappa = 1, thetat = 0),
      n_subjects = 1, n_trials = 2, trial_design = data.frame(n = 1:2)
    ),
    "does not use a trial design"
  )
})

# sbc --------------------------------------------------------------------------

# population draws as given, and an SD for each, which the `(1 | id)` of
# recovery_formula() implies
count_draws <- function(...) {
  k <- seq_len(60)
  values <- list(...)
  sds <- stats::setNames(
    rep(list(0.1), length(values)),
    sub("^b_(.*)_Intercept$", "sd_id__\\1_Intercept", names(values))
  )
  posterior::as_draws_matrix(do.call(cbind, lapply(c(values, sds), function(v) {
    v + 0.001 * k
  })))
}

run_count_sbc <- function(model, data, draws) {
  mock <- sbc_mock_fitter(draws = draws)
  suppressMessages(withCallingHandlers(
    sbc(
      model, recovery_formula(model), data,
      n_sims = 2L, thin_ranks = 1, cores_per_fit = 1, level = "population",
      .fitter = mock$fitter
    ),
    bmmtools_sbc_diagnostics = function(w) invokeRestart("muffleWarning")
  ))
  sbc_dataset_calls(mock$calls)
}

test_that("sbc() simulates sdt_yn with the trials and stimuli of data", {
  skip_if_not_installed("SBC")
  skip_on_cran()
  skip_if_no_sdt()
  data <- data.frame(
    id = rep(c("b", "a", "c"), each = 2), stim = c(1, 0, 1, 0, 0, 1),
    n = c(100, 50, 80, 80, 30, 60), hits = 0
  )
  calls <- run_count_sbc(
    yn_model(), data,
    count_draws(b_d_Intercept = 1, b_criterion_Intercept = 0)
  )
  expect_length(calls, 2L)
  for (call in calls) {
    expect_identical(call$data$stim, data$stim)
    expect_identical(call$data$n, data$n)
    expect_true(all(call$data$hits <= call$data$n))
  }
})

test_that("sbc() simulates sdt_mafc with the trials of data", {
  skip_if_not_installed("SBC")
  skip_on_cran()
  skip_if_no_sdt()
  data <- data.frame(id = c("a", "b", "c"), n = c(100, 40, 70), k = 0)
  calls <- run_count_sbc(mafc_model(), data, count_draws(b_d_Intercept = 1))
  for (call in calls) {
    expect_identical(call$data$n, data$n)
  }
})

test_that("sbc() runs ezdm with the trials of data", {
  skip_if_not_installed("SBC")
  skip_on_cran()
  skip_if_not_installed("bmm")
  data <- data.frame(id = c("a", "b"), n = c(100, 40), mrt = 0, vrt = 0, nu = 0)
  calls <- run_count_sbc(
    ez_model(), data,
    count_draws(
      b_drift_Intercept = 1, b_bound_Intercept = log(1.2),
      b_ndt_Intercept = log(0.3)
    )
  )
  expect_length(calls, 2L)
  for (call in calls) {
    expect_identical(call$data$n, data$n)
    expect_true(all(call$data$nu <= call$data$n))
  }
})

test_that("sbc() simulates m3 with the row totals of data", {
  skip_if_not_installed("SBC")
  skip_on_cran()
  skip_if_not_installed("bmm")
  data <- data.frame(
    id = c("a", "b"), corr = c(60, 20), other = c(30, 10), npl = c(10, 0)
  )
  calls <- run_count_sbc(
    m3_count_model(), data,
    count_draws(b_c_Intercept = log(3), b_a_Intercept = log(0.5))
  )
  for (call in calls) {
    totals <- rowSums(as.data.frame(call$data)[c("corr", "other", "npl")])
    expect_equal(unname(totals), c(100, 30))
  }
})

test_that("sbc() refuses a design error before the prior fit", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  mock <- sbc_mock_fitter(draws = count_draws(b_drift_Intercept = 1))
  data <- data.frame(id = c("a", "b"), n = c(100, 2), mrt = 0, vrt = 0, nu = 0)
  err <- expect_error(
    sbc(
      ez_model(), recovery_formula(ez_model()), data,
      n_sims = 2L, .fitter = mock$fitter
    ),
    "at least 3"
  )
  expect_no_match(conditionMessage(err), "prior")
  expect_identical(mock$calls$n, 0L)
})

# assertions a wrong implementation would fail (29.1 review, MEDIUM 2) ---------

test_that("each sdt_mafc row is drawn from its own number of trials", {
  skip_if_no_sdt()
  # with the first row's 5 trials for every row, row 2 could not exceed 5
  sim <- simulate_recovery(
    mafc_model(), c(d = 4),
    n_subjects = 1, n_trials = 2,
    trial_design = data.frame(n = c(5, 1000)), seed = 1
  )
  expect_gt(sim$data$k[[2]], 500)
})

test_that("each sdt_yn row is drawn from its own stimulus class", {
  skip_if_no_sdt()
  sim <- simulate_recovery(
    yn_model(), c(d = 3, criterion = 0),
    n_subjects = 1, n_trials = 2,
    trial_design = data.frame(stim = c(0, 1), n = 1000), seed = 1
  )
  # a hit rate far above the false-alarm rate, in the rows the design says
  expect_gt(sim$data$hits[[2]] - sim$data$hits[[1]], 500)
})

test_that("repeated m3 totals land in their own rows", {
  skip_if_not_installed("bmm")
  sizes <- c(100, 10, 100, 10, 50)
  sim <- simulate_recovery(
    m3_count_model(), c(c = log(3), a = log(0.5)),
    n_subjects = 2, n_trials = 5,
    trial_design = data.frame(n_trials = sizes), seed = 1
  )
  totals <- rowSums(as.data.frame(sim$data)[c("corr", "other", "npl")])
  expect_equal(unname(totals), rep(sizes, 2))
})

test_that("sbc() keeps each subject's counts when data interleave subjects", {
  skip_if_not_installed("SBC")
  skip_on_cran()
  skip_if_not_installed("bmm")
  data <- data.frame(
    id = c("a", "b", "a", "b"), n = c(100, 40, 60, 20),
    mrt = 0, vrt = 0, nu = 0
  )
  calls <- run_count_sbc(
    ez_model(), data,
    count_draws(
      b_drift_Intercept = 1, b_bound_Intercept = log(1.2),
      b_ndt_Intercept = log(0.3)
    )
  )
  for (call in calls) {
    by_id <- split(call$data$n, as.character(call$data$id))
    expect_equal(by_id, list(a = c(100, 60), b = c(40, 20)))
  }
})

test_that("sbc()'s design errors name data, not trial_design", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  mock <- sbc_mock_fitter(draws = count_draws(b_drift_Intercept = 1))
  expect_error(
    sbc(
      ez_model(), recovery_formula(ez_model()),
      data.frame(id = c("a", "b"), n = c(100, 2), mrt = 0, vrt = 0, nu = 0),
      n_sims = 2L, .fitter = mock$fitter
    ),
    "Column \"n\" of `data`"
  )
  # m3: the count is the row total, and data without the categories says so
  m3 <- m3_count_model()
  expect_error(
    sbc(
      m3, recovery_formula(m3),
      data.frame(id = c("a", "b"), corr = c(10, 0), other = 0, npl = 0),
      n_sims = 2L, .fitter = mock$fitter
    ),
    "total of each row"
  )
  expect_error(
    sbc(
      m3, recovery_formula(m3), data.frame(id = c("a", "b"), x = 1),
      n_sims = 2L, .fitter = mock$fitter
    ),
    "lacks the columns"
  )
  expect_error(
    sbc(
      m3, recovery_formula(m3),
      data.frame(id = c("a", "b"), corr = "1", other = 1, npl = 1),
      n_sims = 2L, .fitter = mock$fitter
    ),
    "must be numeric"
  )
  expect_identical(mock$calls$n, 0L)
})

test_that("a factor stimulus is refused, not passed to the generator", {
  skip_if_no_sdt()
  expect_error(
    simulate_recovery(
      yn_model(), c(d = 1, criterion = 0),
      n_subjects = 1, n_trials = 2,
      trial_design = data.frame(stim = factor(c(1, 0)), n = 10)
    ),
    "as numbers"
  )
})

test_that("the m3 total message renders without a stray backslash", {
  skip_if_not_installed("bmm")
  err <- expect_error(
    check_count_design(
      data.frame(n_trials = c(10, 0)), m3_count_model(),
      in_data = TRUE
    )
  )
  expect_match(
    conditionMessage(err),
    "The total of each row of `data` over \"corr\", \"other\", and \"npl\"",
    fixed = TRUE
  )
  expect_no_match(conditionMessage(err), "\\", fixed = TRUE)
})

test_that("sbc() reads a factor or character 0/1 stimulus as bmm does", {
  skip_if_not_installed("SBC")
  skip_on_cran()
  skip_if_no_sdt()
  draws <- count_draws(b_d_Intercept = 1, b_criterion_Intercept = 0)
  for (stim in list(factor(c(1, 0)), c("1", "0"))) {
    data <- data.frame(id = c("a", "a"), n = 10, hits = 1)
    data$stim <- stim
    calls <- run_count_sbc(yn_model(), data, draws)
    for (call in calls) {
      expect_identical(call$data$stim, c(1, 0))
    }
  }
  # labels bmm would not read are refused before the prior fit
  data <- data.frame(id = "a", stim = c("signal", "noise"), n = 10, hits = 1)
  expect_error(run_count_sbc(yn_model(), data, draws), "0 or 1")
})
