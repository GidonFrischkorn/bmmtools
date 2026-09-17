test_that("no_pooling_formula() writes one cell-means formula per parameter", {
  skip_if_not_installed("bmm")
  f <- no_pooling_formula(fake_ml_model(), by = "id")
  expect_s3_class(f, "bmmformula")
  expect_named(f, c("kappa", "thetat"))
  expect_equal(
    vapply(f, function(x) deparse1(x), character(1)),
    c(kappa = "kappa ~ 0 + id", thetat = "thetat ~ 0 + id")
  )
})

test_that("no_pooling_formula() honours another grouping column", {
  skip_if_not_installed("bmm")
  f <- no_pooling_formula(fake_ml_model(free = "kappa"), by = "subject")
  expect_equal(deparse1(f[[1L]]), "kappa ~ 0 + subject")
})

test_that("no_pooling_formula() has no group-level term", {
  skip_if_not_installed("bmm")
  f <- no_pooling_formula(fake_ml_model(), by = "id")
  # the reduction is the whole point: a `(1 | id)` left in place would be a
  # partially pooled fit labelled "ml"
  expect_false(any(grepl("|", vapply(f, deparse1, character(1)), fixed = TRUE)))
})

test_that("a user formula passes through untouched", {
  skip_if_not_installed("bmm")
  mine <- bmm::bmf(kappa ~ 0 + id)
  expect_identical(no_pooling_formula(fake_ml_model(), "id", formula = mine), mine)
})

test_that("by must be a single column in this version", {
  skip_if_not_installed("bmm")
  expect_error(
    no_pooling_formula(fake_ml_model(), by = c("id", "condition")),
    "single column"
  )
  expect_error(no_pooling_formula(fake_ml_model(), by = 1), "single column")
})

test_that("flat_ml_prior() blanks the class-level prior of every parameter", {
  skip_if_not_installed("brms")
  p <- flat_ml_prior(fake_ml_model())
  expect_s3_class(p, "brmsprior")
  expect_setequal(p$nlpar, c("kappa", "thetat"))
  expect_true(all(p$class == "b"))
  # an empty prior string is what brms calls flat
  expect_true(all(!nzchar(p$prior)))
})

# the rename ---------------------------------------------------------------

ml_population_rows <- function(terms, estimate = 0, sd = 0.5) {
  tibble::tibble(
    term = terms,
    estimate = rep_len(estimate, length(terms)),
    ci_low = rep_len(estimate, length(terms)) - sd,
    ci_high = rep_len(estimate, length(terms)) + sd,
    ci_method = "eti", ci_level = 0.95,
    rhat = 1, ess_bulk = 1000, ess_tail = 1000,
    level = "population", id = NA_character_,
    converged = TRUE, estimator = "bayes"
  )
}

test_that("ml_subject_rows() maps <par>_<by><level> onto term and id", {
  est <- ml_population_rows(c("kappa_id1", "kappa_id2", "thetat_id1", "thetat_id2"))
  out <- ml_subject_rows(est, c("kappa", "thetat"), c("1", "2"), "id")
  expect_equal(out$term, c("kappa", "kappa", "thetat", "thetat"))
  expect_equal(out$id, c("1", "2", "1", "2"))
  expect_true(all(out$level == "subject"))
})

test_that("ml_subject_rows() marks the rows as a laplace ML estimate", {
  est <- ml_population_rows(c("kappa_id1", "kappa_id2"))
  out <- ml_subject_rows(est, "kappa", c("1", "2"), "id")
  expect_true(all(out$ci_method == "laplace"))
  expect_true(all(out$estimator == "ml"))
})

test_that("ml_subject_rows() blanks rhat and ESS", {
  # a Laplace fit is one chain of iid draws: posterior returns a finite
  # split-Rhat and ESS that mean nothing, and check_convergence() would
  # otherwise gate on them
  est <- ml_population_rows(c("kappa_id1", "kappa_id2"))
  out <- ml_subject_rows(est, "kappa", c("1", "2"), "id")
  expect_true(all(is.na(out$rhat)))
  expect_true(all(is.na(out$ess_bulk)))
  expect_true(all(is.na(out$ess_tail)))
})

test_that("an estimate outside the link range is reported, not dropped", {
  est <- ml_population_rows(c("kappa_id1", "kappa_id2"), estimate = c(1, 40))
  out <- ml_subject_rows(est, "kappa", c("1", "2"), "id", max_abs_link = 20)
  expect_equal(nrow(out), 2L)
  expect_equal(out$converged, c(TRUE, FALSE))
  expect_equal(out$estimate, c(1, NA_real_))
  expect_true(is.na(out$ci_low[2L]))
  expect_true(is.na(out$ci_high[2L]))
})

test_that("a non-finite estimate fails the convergence rule", {
  est <- ml_population_rows(c("kappa_id1", "kappa_id2"), estimate = c(1, Inf))
  out <- ml_subject_rows(est, "kappa", c("1", "2"), "id")
  expect_equal(out$converged, c(TRUE, FALSE))
  expect_true(is.na(out$estimate[2L]))
})

test_that("a single level keeps the bare parameter name", {
  # with one coefficient coefficient_terms() does not append it, so the
  # term arrives as "kappa" rather than "kappa_id1"
  est <- ml_population_rows(c("kappa", "thetat"))
  out <- ml_subject_rows(est, c("kappa", "thetat"), "7", "id")
  expect_equal(out$term, c("kappa", "thetat"))
  expect_equal(out$id, c("7", "7"))
})

test_that("an id containing an underscore still maps correctly", {
  # split_coefficient() splits on the LAST underscore, so re-parsing the
  # name would read "kappa_idsub" as the parameter
  est <- ml_population_rows(c("kappa_idsub_01", "kappa_idsub_02"))
  out <- ml_subject_rows(est, "kappa", c("sub_01", "sub_02"), "id")
  expect_equal(out$term, c("kappa", "kappa"))
  expect_equal(out$id, c("sub_01", "sub_02"))
})

test_that("an unrecognised coefficient is an error that names it", {
  est <- ml_population_rows(c("kappa_id1", "kappa_task2"))
  expect_error(
    ml_subject_rows(est, "kappa", "1", "id"),
    "kappa_task2"
  )
})

# fit_ml() -----------------------------------------------------------------

test_that("fit_ml() runs ONE fit for all subjects, on the full data", {
  skip_if_not_installed("bmm")
  mock <- ml_mock_fitter()
  out <- fit_ml(
    fake_ml_model(), fake_ml_data(n = 6L),
    by = "id",
    .fitter = mock$fitter, file = withr::local_tempfile()
  )
  # the joint no-pooling mode IS the vector of per-subject modes, so a fit
  # per subject would be six times the work for the same answer
  expect_equal(mock$calls$n, 1L)
  expect_equal(nrow(mock$calls$last$data), 30L)
  expect_setequal(levels(mock$calls$last$data$id), as.character(1:6))
})

test_that("fit_ml() sends a no-pooling formula to the fitter", {
  skip_if_not_installed("bmm")
  mock <- ml_mock_fitter()
  fit_ml(
    fake_ml_model(), fake_ml_data(),
    .fitter = mock$fitter,
    file = withr::local_tempfile()
  )
  written <- vapply(mock$calls$last$formula, deparse1, character(1))
  expect_true(all(grepl("~ 0 \\+ id$", written)))
  expect_false(any(grepl("|", written, fixed = TRUE)))
})

test_that("fit_ml() asks for a laplace fit on the cmdstanr backend", {
  skip_if_not_installed("bmm")
  mock <- ml_mock_fitter()
  fit_ml(
    fake_ml_model(), fake_ml_data(),
    .fitter = mock$fitter,
    file = withr::local_tempfile()
  )
  expect_equal(mock$calls$last$dots$algorithm, "laplace")
  expect_equal(mock$calls$last$dots$backend, "cmdstanr")
})

test_that("fit_ml() returns subject rows satisfying the estimates contract", {
  skip_if_not_installed("bmm")
  mock <- ml_mock_fitter()
  out <- fit_ml(
    fake_ml_model(), fake_ml_data(n = 3L),
    .fitter = mock$fitter,
    file = withr::local_tempfile()
  )
  expect_s3_class(out, "bmmtools_ml")
  expect_named(out, names(estimates_contract()))
  expect_equal(nrow(out), 6L)
  expect_setequal(out$term, c("kappa", "thetat"))
  expect_setequal(out$id, as.character(1:3))
  expect_true(all(out$level == "subject"))
  expect_true(all(out$estimator == "ml"))
})

test_that("fit_ml() records per-cell status in ml_cells", {
  skip_if_not_installed("bmm")
  mock <- ml_mock_fitter(estimates = c(kappa_id2 = 50))
  out <- fit_ml(
    fake_ml_model(free = "kappa"), fake_ml_data(n = 3L),
    .fitter = mock$fitter, file = withr::local_tempfile()
  )
  cells <- attr(out, "ml_cells")
  expect_s3_class(cells, "tbl_df")
  expect_equal(nrow(cells), 3L)
  expect_named(cells, c("id", "n_terms", "n_converged", "converged"))
  expect_equal(cells$converged, c(TRUE, FALSE, TRUE))
})

test_that("fit_ml() reports a failed cell rather than dropping it", {
  skip_if_not_installed("bmm")
  mock <- ml_mock_fitter(estimates = c(kappa_id2 = 50))
  out <- fit_ml(
    fake_ml_model(free = "kappa"), fake_ml_data(n = 3L),
    .fitter = mock$fitter, file = withr::local_tempfile()
  )
  # dropping would score ML on the two easiest subjects and the Bayesian
  # estimator on all three, and report the difference as the estimator's
  expect_equal(nrow(out), 3L)
  expect_equal(sum(is.na(out$estimate)), 1L)
  expect_equal(out$converged, c(TRUE, FALSE, TRUE))
})

test_that("fit_ml() passes flat priors by default and none on request", {
  skip_if_not_installed("bmm")
  skip_if_not_installed("brms")
  mock <- ml_mock_fitter()
  fit_ml(
    fake_ml_model(), fake_ml_data(),
    .fitter = mock$fitter,
    file = withr::local_tempfile()
  )
  expect_s3_class(mock$calls$last$prior, "brmsprior")
  expect_true(all(!nzchar(mock$calls$last$prior$prior)))

  fit_ml(
    fake_ml_model(), fake_ml_data(),
    prior = "default",
    .fitter = mock$fitter, file = withr::local_tempfile()
  )
  expect_null(mock$calls$last$prior)
})

test_that("fit_ml() reuses a cached fit", {
  skip_if_not_installed("bmm")
  mock <- ml_mock_fitter()
  file <- withr::local_tempfile()
  fit_ml(fake_ml_model(), fake_ml_data(), .fitter = mock$fitter, file = file)
  suppressMessages(
    fit_ml(fake_ml_model(), fake_ml_data(), .fitter = mock$fitter, file = file)
  )
  expect_equal(mock$calls$n, 1L)
})


# the optim route ------------------------------------------------------
#
# Every test below runs without Stan and without a compiler, which is the
# point of the route: it is the only one a runner can exercise end to end.
# Most use `nll` with an objective whose minimum and curvature are known in
# closed form, so the expected estimate and the expected standard error are
# written out rather than recomputed by the same code under test.
#
# With the toy model's `log` link, pars$kappa is exp(theta), so an objective
# written in log(pars$kappa) is an objective in theta itself.

test_that("the optim route finds the analytic optimum and its Wald interval", {
  skip_if_not_installed("bmm")
  dat <- fake_ml_data(n = 3L, n_trials = 5L)
  out <- fit_ml(
    fake_ml_model(free = "kappa"), dat,
    method = "optim", nll = ml_square_nll()
  )
  # the minimum of sum((theta - y)^2)/2 is mean(y); the Hessian is n, so
  # the standard error is 1/sqrt(n)
  expected <- vapply(split(dat$y, dat$id), mean, numeric(1))
  se <- 1 / sqrt(5)
  expect_equal(out$estimate, unname(expected), tolerance = 1e-6)
  expect_equal(
    out$ci_low, unname(expected) - stats::qnorm(0.975) * se,
    tolerance = 1e-5
  )
  expect_equal(
    out$ci_high, unname(expected) + stats::qnorm(0.975) * se,
    tolerance = 1e-5
  )
  expect_true(all(out$converged))
})

test_that("the optim route fills the contract the Stan route fills", {
  skip_if_not_installed("bmm")
  out <- fit_ml(
    fake_ml_model(free = "kappa"), fake_ml_data(n = 3L),
    method = "optim", nll = ml_square_nll()
  )
  expect_s3_class(out, "bmmtools_ml")
  expect_named(out, names(estimates_contract()))
  expect_equal(unique(out$ci_method), "wald")
  expect_equal(unique(out$estimator), "ml")
  expect_equal(unique(out$level), "subject")
  # no draws on either route, so no draw diagnostics
  expect_true(all(is.na(out$rhat)))
  expect_true(all(is.na(out$ess_bulk)))
  expect_true(all(is.na(out$ess_tail)))
})

test_that("both routes return the same columns and types", {
  skip_if_not_installed("bmm")
  mock <- ml_mock_fitter(estimates = c(kappa_id1 = 1.2))
  stan <- fit_ml(
    fake_ml_model(free = "kappa"), fake_ml_data(n = 3L),
    .fitter = mock$fitter, file = withr::local_tempfile()
  )
  optim <- fit_ml(
    fake_ml_model(free = "kappa"), fake_ml_data(n = 3L),
    method = "optim", nll = ml_square_nll()
  )
  # decision 37: the routes differ in ci_method and in nothing else
  expect_identical(names(stan), names(optim))
  expect_identical(
    vapply(stan, function(x) class(x)[1], character(1)),
    vapply(optim, function(x) class(x)[1], character(1))
  )
  expect_s3_class(dplyr::bind_rows(stan, optim), "tbl_df")
})

test_that("ci_level sets the width of the Wald interval", {
  skip_if_not_installed("bmm")
  args <- list(
    fake_ml_model(free = "kappa"), fake_ml_data(n = 2L),
    method = "optim", nll = ml_square_nll()
  )
  wide <- do.call(fit_ml, c(args, list(ci_level = 0.95)))
  narrow <- do.call(fit_ml, c(args, list(ci_level = 0.5)))
  expect_true(all((narrow$ci_high - narrow$ci_low) <
    (wide$ci_high - wide$ci_low)))
  expect_equal(unique(narrow$ci_level), 0.5)
})

test_that("a singular Hessian costs the interval and not the estimate", {
  skip_if_not_installed("bmm")
  # an objective with no curvature: optim stops at the start value and the
  # Hessian is zero, so solve() fails. Gidon's decision 2026-09-17 is that
  # such a row keeps its estimate and its verdict and loses the interval,
  # so the subject still counts in bias, rmse and r.
  flat <- function(pars, data, model) 0
  out <- fit_ml(
    fake_ml_model(free = "kappa"), fake_ml_data(n = 2L),
    method = "optim", nll = flat
  )
  expect_true(all(is.finite(out$estimate)))
  expect_true(all(out$converged))
  expect_true(all(is.na(out$ci_low)))
  expect_true(all(is.na(out$ci_high)))
})

test_that("an estimate outside max_abs_link is not converged", {
  skip_if_not_installed("bmm")
  far <- function(pars, data, model) (log(pars$kappa) - 50)^2
  out <- fit_ml(
    fake_ml_model(free = "kappa"), fake_ml_data(n = 2L),
    method = "optim", nll = far, max_abs_link = 20
  )
  # decision 40: the row stays so that n is visibly equal across estimators
  expect_equal(nrow(out), 2L)
  expect_false(any(out$converged))
  expect_true(all(is.na(out$estimate)))
})

test_that("start is validated and is where the optimiser begins", {
  skip_if_not_installed("bmm")
  model <- fake_ml_model(free = "kappa")
  dat <- fake_ml_data(n = 2L)
  flat <- function(pars, data, model) 0
  expect_error(
    fit_ml(model, dat, method = "optim", nll = flat, start = c(3)),
    "named numeric vector"
  )
  expect_error(
    fit_ml(model, dat, method = "optim", nll = flat, start = c(sigma = 3)),
    "free parameters"
  )
  # with a flat objective the optimiser cannot move, so the estimate is the
  # start value and nothing else
  out <- fit_ml(model, dat,
    method = "optim", nll = flat,
    start = c(kappa = 3)
  )
  expect_equal(out$estimate, rep(3, 2L))
})

test_that("each route refuses the other's arguments by name", {
  skip_if_not_installed("bmm")
  model <- fake_ml_model(free = "kappa")
  dat <- fake_ml_data(n = 2L)
  nll <- ml_square_nll()
  expect_error(
    fit_ml(model, dat, method = "optim", nll = nll, prior = "default"),
    "no priors"
  )
  expect_error(
    fit_ml(model, dat,
      method = "optim", nll = nll,
      formula = bmm::bmf(kappa ~ 0 + id)
    ),
    "builds no formula"
  )
  expect_error(
    fit_ml(model, dat, method = "optim", nll = nll, .fitter = function(...) 1),
    "calls no fitter"
  )
  expect_error(
    fit_ml(model, dat, method = "stan", nll = nll),
    "density hook"
  )
})

test_that("nll must be a function and is required without a density", {
  skip_if_not_installed("bmm")
  model <- fake_ml_model(free = "kappa")
  dat <- fake_ml_data(n = 2L)
  expect_error(
    fit_ml(model, dat, method = "optim", nll = 3),
    "must be a function"
  )
  # the toy model is not one of the six, so there is no density to fall
  # back on and the error has to say what the alternatives are
  expect_error(fit_ml(model, dat, method = "optim"), "has no density")
  expect_error(fit_ml(model, dat, method = "optim"), "stan")
})

test_that("ml_cells counts subjects on the optim route too", {
  skip_if_not_installed("bmm")
  out <- fit_ml(
    fake_ml_model(free = "kappa"), fake_ml_data(n = 4L),
    method = "optim", nll = ml_square_nll()
  )
  cells <- attr(out, "ml_cells")
  expect_equal(nrow(cells), 4L)
  expect_equal(cells$n_converged, rep(1L, 4L))
  expect_true(all(cells$converged))
})

test_that("rows follow the data's level order, not a character sort", {
  skip_if_not_installed("bmm")
  dat <- fake_ml_data(levels = as.character(c(1, 2, 10, 11)))
  out <- fit_ml(
    fake_ml_model(free = "kappa"), dat,
    method = "optim", nll = ml_square_nll()
  )
  expect_equal(out$id, as.character(c(1, 2, 10, 11)))
})

test_that("print names the Wald interval and says there are no priors", {
  skip_if_not_installed("bmm")
  out <- fit_ml(
    fake_ml_model(free = "kappa"), fake_ml_data(n = 3L),
    method = "optim", nll = ml_square_nll()
  )
  txt <- paste(format(out), collapse = "\n")
  expect_match(txt, "optim/wald")
  expect_match(txt, "no priors")
  expect_false(grepl("flat priors", txt, fixed = TRUE))
})

# the density adapters -------------------------------------------------

test_that("density_for() knows six models and nothing else", {
  skip_if_not_installed("bmm")
  skip_if_no_bmm_sdt()
  expect_equal(adapter_classes(), c(
    "sdt_yn", "sdt_mafc", "ezdm", "ddm", "mixture2p", "sdm"
  ))
  known <- list(
    bmm::sdt_yn(response = "hits", stimulus = "stim", n_trials = "n"),
    bmm::sdt_mafc(response = "k", n_trials = "n", m = 4),
    bmm::ezdm(mean_rt = "m", var_rt = "v", n_upper = "u", n_trials = "n"),
    bmm::ddm(rt = "rt", response = "resp"),
    bmm::mixture2p(resp_error = "y"),
    bmm::sdm(resp_error = "y")
  )
  for (m in known) {
    expect_true(is.function(density_for(m)))
  }
  # the same cap generator_for() applies, for the same reason
  expect_null(density_for(fake_ml_model()))
})

test_that("every density adapter scores data its own generator made", {
  skip_if_not_installed("bmm")
  # not guarded on the SDT stack: ml_density_cases() drops those two when
  # the installed bmm has no signal-detection models, so a runner still
  # exercises the other four rather than skipping the lot
  cases <- ml_density_cases()
  expect_gte(length(cases), 4L)
  for (case in cases) {
    pars <- natural_pars(case$link, case$model)
    dat <- generator_for(case$model)(pars, 20L, case$model)
    ll <- density_for(case$model)(pars, dat, case$model)
    expect_true(is.numeric(ll), info = case$label)
    expect_true(all(is.finite(ll)), info = case$label)
    # every adapter is on the log scale, checked against the same bmm call
    # made independently here. A log density is NOT bounded above by zero
    # for a continuous model --- dddm() and dezdm() return positive values
    # wherever the density exceeds 1 --- so the comparison is with the
    # density itself, not with a sign.
    expect_equal(ll, log(case$raw(pars, dat, case$model)), info = case$label)
  }
})

test_that("the optim route recovers a real bmm model without Stan", {
  skip_if_not_installed("bmm")
  model <- bmm::mixture2p(resp_error = "y")
  sim <- simulate_recovery(
    model,
    pars = c(kappa = 2, thetat = 1), sds = c(kappa = 0.3, thetat = 0.5),
    n_subjects = 12L, n_trials = 300L, seed = 85
  )
  out <- fit_ml(model, sim$data, method = "optim")
  expect_s3_class(out, "bmmtools_ml")
  expect_equal(unique(out$ci_method), "wald")
  expect_true(all(out$converged))
  s <- summary(suppressMessages(recover_subjects(out, sim$truth$subjects)))
  # a loose gate: this asserts the route is wired to the right likelihood,
  # not how well ML does. The measured contrast lives in local/dev/sim/.
  expect_true(all(s$r > 0.5))
  expect_true(all(abs(s$bias) < 0.5))
})

test_that("fit_ml() output scores through recover_subjects()", {
  skip_if_not_installed("bmm")
  mock <- ml_mock_fitter(estimates = c(
    kappa_id1 = 1.2, kappa_id2 = 0.8,
    kappa_id3 = 1.6
  ))
  out <- fit_ml(
    fake_ml_model(free = "kappa"), fake_ml_data(n = 3L),
    .fitter = mock$fitter, file = withr::local_tempfile()
  )
  truth <- tibble::tibble(
    id = as.character(1:3), term = "kappa", true_value = c(1, 1, 2)
  )
  rec <- suppressMessages(recover_subjects(out, truth))
  s <- summary(rec)
  expect_equal(s$estimator, "ml")
  expect_equal(s$n, 3L)
})

test_that("print.bmmtools_ml() names the cells, the method and the prior", {
  skip_if_not_installed("bmm")
  mock <- ml_mock_fitter(estimates = c(kappa_id2 = 50))
  out <- fit_ml(
    fake_ml_model(free = "kappa"), fake_ml_data(n = 3L),
    .fitter = mock$fitter, file = withr::local_tempfile()
  )
  txt <- paste(format(out), collapse = "\n")
  expect_match(txt, "bmmtools_ml")
  expect_match(txt, "3 cells")
  expect_match(txt, "2 converged")
  expect_match(txt, "laplace")
  expect_match(txt, "flat")
})

test_that("a bmmtools_ml drops its class when a contract column goes", {
  skip_if_not_installed("bmm")
  mock <- ml_mock_fitter()
  out <- fit_ml(
    fake_ml_model(), fake_ml_data(),
    .fitter = mock$fitter,
    file = withr::local_tempfile()
  )
  expect_s3_class(out[, setdiff(names(out), "estimator")], "tbl_df")
  expect_false(inherits(out[, setdiff(names(out), "estimator")], "bmmtools_ml"))
  expect_s3_class(dplyr::filter(out, .data$term == "kappa"), "bmmtools_ml")
})

# what the review of milestone 8 found ------------------------------------

test_that("a coefficient missing from the fit is an error, not a lost row", {
  est <- fake_ml_population(c("kappa_id1", "kappa_id3"))
  # the map is a cross product of `free` and the levels, so level 2 is
  # expected and absent. Dropping it silently would score the two
  # estimators on different subjects, which decision 40 exists to prevent.
  expect_error(
    ml_subject_rows(est, "kappa", c("1", "2", "3"), "id"),
    "kappa_id2"
  )
  full <- fake_ml_population(c("kappa_id1", "kappa_id2", "kappa_id3"))
  expect_equal(nrow(ml_subject_rows(full, "kappa", c("1", "2", "3"), "id")), 3L)
})

test_that("the arguments fit_ml() sets itself are refused in the dots", {
  skip_if_not_installed("bmm")
  mock <- ml_mock_fitter()
  run <- function(...) {
    fit_ml(
      fake_ml_model(), fake_ml_data(),
      .fitter = mock$fitter,
      file = withr::local_tempfile(), ...
    )
  }
  # both are spliced into the fitter call beside the value fit_ml() passes,
  # so without a guard the user sees R's bare "matched by multiple actual
  # arguments" and the cache key records the value they did not ask for
  expect_error(run(algorithm = "sampling"), "algorithm")
  expect_error(run(backend = "rstan"), "backend")
})

test_that("an objective that stops being a number away from the start errors", {
  skip_if_not_installed("bmm")
  # the start-value check catches an objective that is wrong at zero; this
  # one is a number at the start and a string one step away, which the
  # sentinel would turn into an ordinary "the density failed here"
  bad <- function(pars, data, model) {
    if (!isTRUE(all.equal(log(pars$kappa), 0))) {
      return("oops")
    }
    0
  }
  expect_error(
    fit_ml(
      fake_ml_model(free = "kappa"), fake_ml_data(n = 1L),
      method = "optim", nll = bad
    ),
    "single number"
  )
})

test_that("arguments of the other route are refused, not silently ignored", {
  skip_if_not_installed("bmm")
  run <- function(...) {
    fit_ml(
      fake_ml_model(free = "kappa"), fake_ml_data(n = 2L),
      method = "optim", nll = ml_square_nll(), ...
    )
  }
  # `draws` is the size of the Laplace sample and does nothing here; an
  # unknown name has nothing to reach, because the optim route calls no
  # fitter and so consumes no dots
  expect_error(run(draws = 50), "draws")
  expect_error(run(ci_levl = 0.8), "ci_levl")
})

test_that("a grouping column with no level is a named error", {
  skip_if_not_installed("bmm")
  data <- fake_ml_data(n = 2L)
  data$id <- NA_character_
  expect_error(
    fit_ml(
      fake_ml_model(free = "kappa"), data,
      method = "optim", nll = ml_square_nll()
    ),
    "id"
  )
})

test_that("the remaining argument guards of fit_ml() have their errors", {
  skip_if_not_installed("bmm")
  run <- function(...) {
    fit_ml(
      fake_ml_model(free = "kappa"), fake_ml_data(n = 2L),
      method = "optim", nll = ml_square_nll(), ...
    )
  }
  expect_error(run(by = 1), "by")
  expect_error(run(by = "subject"), "no column")
  expect_error(run(max_abs_link = -1), "max_abs_link")
  expect_error(
    fit_ml(fake_ml_model(), "not a data frame", method = "optim"),
    "data frame"
  )
})

test_that("a bmmtools_ml prints what format() builds", {
  skip_if_not_installed("bmm")
  out <- fit_ml(
    fake_ml_model(free = "kappa"), fake_ml_data(n = 2L),
    method = "optim", nll = ml_square_nll()
  )
  expect_output(print(out), "bmmtools_ml")
  expect_identical(print(out), out)
  # a tibble with the contract but no rows still names its interval method
  expect_match(paste(format(out[0, ]), collapse = "\n"), "wald")
})
