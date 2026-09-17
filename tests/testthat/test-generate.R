# Tests for simulate_recovery(), recovery_formula() and the adapter
# table, written against local/dev/spec-milestone-3-generate-layer.md sections
# 1 to 3. They need bmm for the model objects and the r*() generators;
# nothing here compiles Stan. Generated data are validated against bmm's
# own data checks through brms's mock backend.

mock_bmm <- function(formula, data, model) {
  bmm::bmm(
    formula, data, model,
    backend = "mock", mock_fit = 1, rename = FALSE, silent = 2
  )
}

m2p_pars <- c(kappa = log(8), thetat = stats::qlogis(0.75))

# simulate_recovery -----------------------------------------------------

test_that("mixture2p simulates the right shape and truth tables", {
  skip_if_not_installed("bmm")
  model <- bmm::mixture2p(resp_error = "y")

  sim <- simulate_recovery(
    model, m2p_pars,
    n_subjects = 5, n_trials = 20,
    sds = c(kappa = 0.3, thetat = 0.5), seed = 1
  )

  expect_s3_class(sim, "bmmtools_simulation")
  expect_equal(nrow(sim$data), 100L)
  expect_named(sim$data, c("id", "y"))
  expect_s3_class(sim$data$id, "factor")
  expect_setequal(sim$truth$population$term, c("kappa", "thetat"))
  expect_equal(sim$truth$population$true_value[
    sim$truth$population$term == "kappa"
  ], log(8))
  expect_equal(nrow(sim$truth$subjects), 10L)
  expect_named(sim$truth$subjects, c("id", "term", "true_value"))
  expect_equal(sim$seed, 1)
  expect_equal(sim$generator, "adapter:mixture2p")

  only_kappa <- simulate_recovery(
    model, m2p_pars,
    n_subjects = 5, n_trials = 20,
    sds = c(kappa = 0.3), seed = 1
  )
  expect_equal(nrow(only_kappa$truth$subjects), 5L)
  expect_setequal(only_kappa$truth$subjects$term, "kappa")
  expect_equal(only_kappa$sds, c(kappa = 0.3, thetat = 0))
})

test_that("no sds means no subject truth and identical subjects", {
  skip_if_not_installed("bmm")
  model <- bmm::mixture2p(resp_error = "y")
  sim <- simulate_recovery(model, m2p_pars, n_subjects = 3, n_trials = 5)
  expect_equal(nrow(sim$truth$subjects), 0L)
  expect_named(sim$truth$subjects, c("id", "term", "true_value"))
  expect_true(is.na(sim$seed))
})

test_that("the seed reproduces data and truth", {
  skip_if_not_installed("bmm")
  model <- bmm::mixture2p(resp_error = "y")
  run <- function(seed) {
    simulate_recovery(
      model, m2p_pars,
      n_subjects = 4, n_trials = 10,
      sds = c(kappa = 0.3), seed = seed
    )
  }
  a <- run(11)
  b <- run(11)
  c <- run(12)
  expect_identical(a$data, b$data)
  expect_identical(a$truth, b$truth)
  expect_false(identical(a$data, c$data))
  expect_false(identical(a$truth$subjects, c$truth$subjects))
})

test_that("subject_pars are used verbatim instead of drawn", {
  skip_if_not_installed("bmm")
  model <- bmm::mixture2p(resp_error = "y")
  given <- tibble::tibble(
    id = as.character(1:3),
    term = "kappa",
    true_value = c(1.5, 2.0, 2.5)
  )
  sim <- simulate_recovery(
    model, m2p_pars,
    n_subjects = 3, n_trials = 5,
    sds = c(kappa = 0.3), subject_pars = given, seed = 1
  )
  expect_equal(sim$truth$subjects, given)

  # passing the recorded subjects back reproduces the data under the
  # same seed: the draw step is skipped and the generator sees the same
  # values
  again <- simulate_recovery(
    model, m2p_pars,
    n_subjects = 3, n_trials = 5,
    sds = c(kappa = 0.3), subject_pars = sim$truth$subjects, seed = 1
  )
  expect_identical(again$data, sim$data)
})

test_that("the generator receives natural-scale values, fixed ones too", {
  skip_if_not_installed("bmm")
  skip_if_no_bmm_sdt()
  seen <- new.env()
  recorder <- function(pars, n_trials, model) {
    seen$pars <- pars
    seen$n_trials <- n_trials
    data.frame(y = rep(0, n_trials))
  }
  model <- bmm::mixture2p(resp_error = "y")
  sim <- simulate_recovery(
    model, m2p_pars,
    n_subjects = 1, n_trials = 7, generator = recorder
  )
  expect_equal(seen$pars$kappa, 8)
  expect_equal(seen$pars$thetat, 0.75)
  expect_equal(seen$pars$mu1, 2 * atan(0))
  expect_equal(seen$n_trials, 7)
  expect_equal(sim$generator, "user")

  yn <- bmm::sdt_yn(response = "hits", stimulus = "stim", n_trials = "n")
  simulate_recovery(
    yn, c(d = 1.5, criterion = 0.2),
    n_subjects = 1, n_trials = 10,
    generator = function(pars, n_trials, model) {
      seen$yn <- pars
      data.frame(hits = c(5, 2), stim = c(1, 0), n = 10)
    }
  )
  expect_equal(seen$yn$sdratio, 1)
  expect_equal(seen$yn$d, 1.5)
})

test_that("every adapter's data passes bmm's checks (mock backend)", {
  skip_if_not_installed("bmm")
  skip_if_no_bmm_sdt()
  skip_if_not_installed("brms")
  cases <- list(
    list(
      model = bmm::sdt_yn(response = "hits", stimulus = "stim", n_trials = "n"),
      pars = c(d = 1.5, criterion = 0.2)
    ),
    list(
      model = bmm::sdt_mafc(response = "k", n_trials = "n", m = 4),
      pars = c(d = 1)
    ),
    list(
      model = bmm::ezdm(
        mean_rt = "mrt", var_rt = "vrt", n_upper = "nu", n_trials = "n"
      ),
      pars = c(drift = 2, bound = log(1.2), ndt = log(0.3))
    ),
    list(
      model = bmm::ddm(rt = "rt", response = "resp"),
      pars = c(drift = 2, bound = log(1.2), ndt = log(0.3))
    ),
    list(model = bmm::mixture2p(resp_error = "y"), pars = m2p_pars),
    list(
      model = bmm::sdm(resp_error = "y"),
      pars = c(c = log(4), kappa = log(3))
    )
  )
  for (case in cases) {
    label <- class(case$model)[[length(class(case$model))]]
    sim <- simulate_recovery(
      case$model, case$pars,
      n_subjects = 3, n_trials = 30, seed = 5
    )
    expect_true(startsWith(sim$generator, "adapter:"), label = label)
    fit <- mock_bmm(recovery_formula(case$model), sim$data, case$model)
    expect_s3_class(fit, "bmmfit")
  }
})

test_that("adapters lay trials out the way each model expects", {
  skip_if_not_installed("bmm")
  skip_if_no_bmm_sdt()
  yn <- simulate_recovery(
    bmm::sdt_yn(response = "hits", stimulus = "stim", n_trials = "n"),
    c(d = 1.5, criterion = 0.2),
    n_subjects = 2, n_trials = 25, seed = 1
  )
  expect_equal(nrow(yn$data), 4L)
  expect_equal(yn$data$stim, c(1, 0, 1, 0))
  expect_true(all(yn$data$n == 25))
  expect_true(all(yn$data$hits >= 0 & yn$data$hits <= 25))

  mafc <- simulate_recovery(
    bmm::sdt_mafc(response = "k", n_trials = "n", m = 4), c(d = 1),
    n_subjects = 2, n_trials = 25, seed = 1
  )
  expect_equal(nrow(mafc$data), 2L)
  expect_named(mafc$data, c("id", "k", "n"))

  ez <- simulate_recovery(
    bmm::ezdm(mean_rt = "mrt", var_rt = "vrt", n_upper = "nu", n_trials = "n"),
    c(drift = 2, bound = log(1.2), ndt = log(0.3)),
    n_subjects = 2, n_trials = 40, seed = 1
  )
  expect_equal(nrow(ez$data), 2L)
  expect_named(ez$data, c("id", "mrt", "vrt", "nu", "n"))

  ddm <- simulate_recovery(
    bmm::ddm(rt = "rt", response = "resp"),
    c(drift = 2, bound = log(1.2), ndt = log(0.3)),
    n_subjects = 2, n_trials = 15, seed = 1
  )
  expect_equal(nrow(ddm$data), 30L)
  expect_named(ddm$data, c("id", "rt", "resp"))

  sdm <- simulate_recovery(
    bmm::sdm(resp_error = "err"), c(c = log(4), kappa = log(3)),
    n_subjects = 2, n_trials = 15, seed = 1
  )
  expect_equal(nrow(sdm$data), 30L)
  expect_named(sdm$data, c("id", "err"))
})

test_that("bad inputs are refused with the offending names", {
  skip_if_not_installed("bmm")
  model <- bmm::mixture2p(resp_error = "y")
  expect_error(simulate_recovery(1, m2p_pars, 2, 2), "bmmodel")
  expect_error(
    simulate_recovery(model, c(kappa = 1), 2, 2),
    "thetat"
  )
  expect_error(
    simulate_recovery(model, c(m2p_pars, zeta = 1), 2, 2),
    "zeta"
  )
  expect_error(
    simulate_recovery(model, m2p_pars, 2, 2, sds = c(zeta = 1)),
    "zeta"
  )
  expect_error(
    simulate_recovery(model, m2p_pars, 2, 2, sds = c(mu1 = 1)),
    "mu1"
  )
  expect_error(
    simulate_recovery(model, m2p_pars, 2, 2, sds = c(kappa = -1)),
    "kappa"
  )
  expect_error(
    simulate_recovery(model, m2p_pars, 0, 2),
    "n_subjects"
  )
  expect_error(
    simulate_recovery(
      model, m2p_pars, 2, 2,
      generator = function(pars, n_trials, model) data.frame(z = 1)
    ),
    "y"
  )
  no_adapter <- structure(
    list(links = list(a = "identity"), fixed_parameters = list()),
    class = c("bmmodel", "nothing_here")
  )
  expect_error(
    simulate_recovery(no_adapter, c(a = 1), 2, 2),
    "generator"
  )
  expect_error(
    simulate_recovery(
      model, m2p_pars, 2, 2,
      sds = c(kappa = 0.1),
      subject_pars = tibble::tibble(id = "1", term = "kappa", true_value = 1)
    ),
    "subject_pars"
  )
})

test_that("print gives one line with the model and the subject count", {
  skip_if_not_installed("bmm")
  sim <- simulate_recovery(
    bmm::mixture2p(resp_error = "y"), m2p_pars,
    n_subjects = 6, n_trials = 3, sds = c(kappa = 0.2)
  )
  expect_output(print(sim), "6 subjects")
  expect_output(print(sim), "kappa")
})

# correlated truths (spec 5, section 5.1) ------------------------------

m2p_cors <- function(r) {
  m <- matrix(c(1, r, r, 1), 2, 2)
  dimnames(m) <- rep(list(c("kappa", "thetat")), 2)
  m
}

test_that("NULL and identity cors reproduce the Milestone 3 draws", {
  skip_if_not_installed("bmm")
  model <- bmm::mixture2p(resp_error = "y")
  sds <- c(kappa = 0.3, thetat = 0.5)
  run <- function(cors) {
    simulate_recovery(
      model, m2p_pars,
      n_subjects = 6, n_trials = 10, sds = sds, cors = cors, seed = 9
    )
  }
  none <- run(NULL)
  identity <- run(m2p_cors(0))
  expect_identical(identity$data, none$data)
  expect_identical(identity$truth$subjects, none$truth$subjects)

  oracle <- withr::with_seed(9, vapply(names(m2p_pars), function(p) {
    m2p_pars[[p]] + stats::rnorm(6, 0, sds[[p]])
  }, numeric(6)))
  expect_identical(
    none$truth$subjects$true_value,
    as.double(oracle)
  )
})

test_that("a covariate adds a data column and leaves the parameters alone", {
  skip_if_not_installed("bmm")
  model <- bmm::mixture2p(resp_error = "y")
  sds <- c(kappa = 0.3, thetat = 0.5)
  base <- simulate_recovery(
    model, m2p_pars,
    n_subjects = 5, n_trials = 4, sds = sds, cors = m2p_cors(0.5), seed = 2
  )
  cors <- diag(3)
  dimnames(cors) <- rep(list(c("kappa", "thetat", "G")), 2)
  cors["kappa", "thetat"] <- cors["thetat", "kappa"] <- 0.5
  cors["G", "kappa"] <- cors["kappa", "G"] <- 0.4
  with_g <- simulate_recovery(
    model, m2p_pars,
    n_subjects = 5, n_trials = 4, sds = sds, cors = cors,
    covariates = list(G = c(mean = 100, sd = 15)), seed = 2
  )

  expect_identical(with_g$truth$subjects, base$truth$subjects)
  expect_identical(with_g$data$y, base$data$y)
  expect_named(with_g$data, c("id", "G", "y"))
  expect_type(with_g$data$G, "double")
  # subject-constant, and the same values as the covariate truth
  per_id <- tapply(with_g$data$G, with_g$data$id, function(g) length(unique(g)))
  expect_true(all(per_id == 1L))
  expect_named(with_g$truth$covariates, c("id", "term", "true_value"))
  expect_equal(nrow(with_g$truth$covariates), 5L)
  expect_equal(
    with_g$truth$covariates$true_value,
    as.double(tapply(with_g$data$G, with_g$data$id, unique))
  )
  expect_equal(with_g$covariates, list(G = c(mean = 100, sd = 15)))
})

test_that("the truth gains sd, cor and covariate tables", {
  skip_if_not_installed("bmm")
  skip_if_no_bmm_sdt()
  model <- bmm::sdt_yn(response = "hits", stimulus = "stim", n_trials = "n")
  sim <- simulate_recovery(
    model, c(d = 1, criterion = 0.3),
    n_subjects = 4, n_trials = 10,
    sds = c(d = 0.4, criterion = 0.2),
    covariates = list(G = c(mean = 0, sd = 1)),
    seed = 1
  )
  expect_equal(sim$truth$sd, tibble::tibble(
    term = c("d", "criterion"), true_value = c(0.4, 0.2)
  ))
  expect_named(sim$truth$cor, c("term", "var1", "var2", "true_value"))
  # all three pairs, names sorted within the pair in the C locale (upper
  # case first), zeros included
  expect_setequal(
    sim$truth$cor$term, c("criterion__d", "G__d", "G__criterion")
  )
  pairs <- mapply(
    function(a, b) identical(c(a, b), sort(c(a, b), method = "radix")),
    sim$truth$cor$var1, sim$truth$cor$var2
  )
  expect_true(all(pairs))
  expect_true(all(sim$truth$cor$true_value == 0))
  expect_equal(dim(sim$cors), c(3L, 3L))

  one_varying <- simulate_recovery(
    model, c(d = 1, criterion = 0.3),
    n_subjects = 4, n_trials = 10, sds = c(d = 0.4), seed = 1
  )
  expect_equal(one_varying$truth$sd$term, "d")
  expect_equal(nrow(one_varying$truth$cor), 0L)
  expect_named(one_varying$truth$cor, c("term", "var1", "var2", "true_value"))
  expect_null(one_varying$cors)
  expect_equal(nrow(one_varying$truth$covariates), 0L)
  expect_null(one_varying$covariates)
})

test_that("large samples reach the requested correlation", {
  skip_if_not_installed("bmm")
  sim <- simulate_recovery(
    bmm::mixture2p(resp_error = "y"), m2p_pars,
    n_subjects = 4000, n_trials = 1,
    sds = c(kappa = 0.3, thetat = 0.5), cors = m2p_cors(0.5),
    covariates = list(G = c(mean = 0, sd = 1)),
    generator = function(pars, n_trials, model) data.frame(y = 0),
    seed = 11
  )
  wide <- subjects_wide(sim$truth$subjects)
  expect_lt(abs(stats::cor(wide$kappa, wide$thetat) - 0.5), 0.05)
  expect_equal(
    sim$truth$cor$true_value[sim$truth$cor$term == "kappa__thetat"], 0.5
  )
})

test_that("function-valued truths are evaluated under the seed and stored", {
  skip_if_not_installed("bmm")
  model <- bmm::mixture2p(resp_error = "y")
  run <- function(seed) {
    simulate_recovery(
      model,
      pars = function() c(kappa = stats::rnorm(1, 2, 0.1), thetat = 1),
      sds = function() c(kappa = stats::runif(1, 0.2, 0.4), thetat = 0.5),
      cors = function() m2p_cors(stats::runif(1, 0.2, 0.8)),
      n_subjects = 5, n_trials = 3, seed = seed
    )
  }
  a <- run(21)
  b <- run(21)
  c <- run(22)
  expect_identical(a$truth, b$truth)
  expect_false(identical(a$pars, c$pars))
  expect_false(identical(a$cors, c$cors))
  expect_type(a$pars, "double")
  expect_true(is.matrix(a$cors))
  expect_equal(
    a$truth$cor$true_value, a$cors["kappa", "thetat"]
  )
})

test_that("subject_pars skip the draw, and must carry the covariates", {
  skip_if_not_installed("bmm")
  model <- bmm::mixture2p(resp_error = "y")
  covariates <- list(G = c(mean = 0, sd = 1))
  first <- simulate_recovery(
    model, m2p_pars,
    n_subjects = 3, n_trials = 4, sds = c(kappa = 0.3),
    covariates = covariates, seed = 5
  )
  given <- dplyr::bind_rows(first$truth$subjects, first$truth$covariates)
  again <- simulate_recovery(
    model, m2p_pars,
    n_subjects = 3, n_trials = 4, sds = c(kappa = 0.3),
    covariates = covariates, subject_pars = given, seed = 5
  )
  expect_identical(again$truth$covariates, first$truth$covariates)
  expect_identical(again$truth$subjects, first$truth$subjects)
  expect_identical(again$data$G, first$data$G)

  expect_error(
    simulate_recovery(
      model, m2p_pars,
      n_subjects = 3, n_trials = 4, sds = c(kappa = 0.3),
      covariates = covariates, subject_pars = first$truth$subjects
    ),
    "G"
  )
})

test_that("bad covariates and cors are refused with the offending name", {
  skip_if_not_installed("bmm")
  model <- bmm::mixture2p(resp_error = "y")
  sim <- function(...) {
    simulate_recovery(
      model, m2p_pars,
      n_subjects = 3, n_trials = 2, sds = c(kappa = 0.3, thetat = 0.2), ...
    )
  }
  expect_error(sim(covariates = c(G = 1)), "covariates")
  expect_error(sim(covariates = list(G = c(mean = 0))), "sd")
  expect_error(sim(covariates = list(G = c(mean = 0, sd = 0))), "G")
  expect_error(sim(covariates = list(y = c(mean = 0, sd = 1))), "y")
  expect_error(sim(covariates = list(kappa = c(mean = 0, sd = 1))), "kappa")
  expect_error(sim(covariates = list(id = c(mean = 0, sd = 1))), "id")
  expect_error(sim(covariates = list(g_1 = c(mean = 0, sd = 1))), "g_1")
  expect_error(sim(cors = m2p_cors(2)), "cors")
  expect_error(sim(cors = function() "x"), "function")
  expect_error(sim(pars = function() "x"), "function")
  expect_error(sim(cors = function(row) diag(2)), "no arguments")
  collide <- function(pars, n_trials, model) data.frame(y = 0, G = 1)
  expect_error(
    sim(covariates = list(G = c(mean = 0, sd = 1)), generator = collide),
    "G"
  )
})

test_that("print names covariates and nonzero correlations", {
  skip_if_not_installed("bmm")
  sim <- simulate_recovery(
    bmm::mixture2p(resp_error = "y"), m2p_pars,
    n_subjects = 3, n_trials = 2, sds = c(kappa = 0.3, thetat = 0.2),
    cors = m2p_cors(0.5), covariates = list(G = c(mean = 0, sd = 1))
  )
  expect_output(print(sim), "covariates: G")
  expect_output(print(sim), "kappa__thetat = 0.5")
})

# recovery_formula ------------------------------------------------------

test_that("re_cor = 'all' correlates every random intercept", {
  skip_if_not_installed("bmm")
  skip_if_no_bmm_sdt()
  skip_if_not_installed("brms")
  model <- bmm::mixture2p(resp_error = "y")
  f <- recovery_formula(model, re_cor = "all")
  expect_equal(deparse(f$kappa), "kappa ~ 1 + (1 | p | id)")
  expect_equal(deparse(f$thetat), "thetat ~ 1 + (1 | p | id)")
  expect_identical(
    deparse(recovery_formula(model, re_cor = "none")$kappa),
    "kappa ~ 1 + (1 | id)"
  )
  sim <- simulate_recovery(
    model, m2p_pars,
    n_subjects = 4, n_trials = 20, sds = c(kappa = 0.3, thetat = 0.2),
    cors = m2p_cors(0.5), covariates = list(G = c(mean = 0, sd = 1)), seed = 1
  )
  fit <- mock_bmm(f, sim$data, model)
  expect_true(all(fit$ranef$cor))

  mafc <- bmm::sdt_mafc(response = "k", n_trials = "n", m = 4)
  expect_message(
    single <- recovery_formula(mafc, re_cor = "all"),
    "one"
  )
  expect_equal(deparse(single$d), "d ~ 1 + (1 | id)")
})

test_that("recovery_formula gives every free parameter a random intercept", {
  skip_if_not_installed("bmm")
  f <- recovery_formula(bmm::mixture2p(resp_error = "y"))
  expect_s3_class(f, "bmmformula")
  expect_named(f, c("kappa", "thetat"))
  expect_equal(deparse(f$kappa), "kappa ~ 1 + (1 | id)")
  g <- recovery_formula(bmm::mixture2p(resp_error = "y"), group = "subject")
  expect_equal(deparse(g$thetat), "thetat ~ 1 + (1 | subject)")
})

# the adapter table -----------------------------------------------------

test_that("generator_for knows six models and nothing else", {
  skip_if_not_installed("bmm")
  skip_if_no_bmm_sdt()
  models <- list(
    bmm::sdt_yn(response = "hits", stimulus = "stim", n_trials = "n"),
    bmm::sdt_mafc(response = "k", n_trials = "n", m = 4),
    bmm::ezdm(mean_rt = "mrt", var_rt = "vrt", n_upper = "nu", n_trials = "n"),
    bmm::ddm(rt = "rt", response = "resp"),
    bmm::mixture2p(resp_error = "y"),
    bmm::sdm(resp_error = "y")
  )
  for (m in models) expect_type(generator_for(m), "closure")
  expect_null(generator_for(structure(list(), class = c("bmmodel", "mpt"))))
})

test_that("subject_pars may only give parameters that vary", {
  skip_if_not_installed("bmm")
  model <- bmm::mixture2p(resp_error = "y")
  given <- tibble::tibble(
    id = as.character(1:2), term = "thetat", true_value = c(0.5, 1.5)
  )
  expect_error(
    simulate_recovery(
      model, m2p_pars, 2, 5,
      sds = c(kappa = 0.3), subject_pars = given
    ),
    "thetat"
  )
})

# the task dimension (spec 5, section 5.4) ------------------------------

# A cheap generator: one row per call, carrying the natural-scale values it
# was given, so a test can read them back from the data.
echo_generator <- function(pars, n_trials, model) {
  data.frame(y = 0, kappa_seen = pars$kappa, thetat_seen = pars$thetat)
}

task_pars <- c(kappa = log(8), kappa_task2 = log(4), thetat = 0.5)

test_that("tasks = NULL reproduces the simulations of 273cc47", {
  skip_if_not_installed("bmm")
  # Built by running exactly these calls at commit 273cc47, before the
  # task dimension existed, and saving the fields compared here.
  # Compared with a tolerance, not bit for bit: the fixture was built on
  # macOS and the correlated draw goes through a different LAPACK on the
  # Linux runner, where the values differ in the last one or two ulp.
  # testthat's default tolerance is ~1.5e-8, eight orders of magnitude
  # above that noise, so a changed RNG stream or a changed simulation
  # still fails this loudly --- which is the regression it exists for.
  head <- readRDS(test_path("fixtures", "simulation-273cc47.rds"))
  model <- bmm::mixture2p(resp_error = "y")
  gen <- function(pars, n_trials, model) {
    data.frame(y = stats::rnorm(n_trials, pars$kappa, pars$thetat))
  }
  cors <- diag(3)
  dimnames(cors) <- rep(list(c("kappa", "thetat", "G")), 2)
  cors["kappa", "thetat"] <- cors["thetat", "kappa"] <- 0.5
  cors["G", "kappa"] <- cors["kappa", "G"] <- 0.3
  correlated <- simulate_recovery(
    model, c(kappa = log(8), thetat = stats::qlogis(0.75)),
    n_subjects = 5, n_trials = 4, sds = c(kappa = 0.3, thetat = 0.5),
    cors = cors, covariates = list(G = c(mean = 0, sd = 1)),
    generator = gen, seed = 42
  )
  plain <- simulate_recovery(
    model, c(kappa = log(8), thetat = stats::qlogis(0.75), mu1 = 0),
    n_subjects = 3, n_trials = 6, sds = c(kappa = 0.3), generator = gen,
    seed = 7
  )
  fields <- names(head$correlated)
  expect_equal(unclass(correlated)[fields], head$correlated)
  expect_equal(unclass(plain)[fields], head$plain)
  expect_null(correlated$tasks)
  expect_null(correlated$task_col)
  expect_true(all(c("tasks", "task_col") %in% names(correlated)))
})

test_that("tasks expand pars and sds to full terms, full terms winning", {
  skip_if_not_installed("bmm")
  sim <- simulate_recovery(
    bmm::mixture2p(resp_error = "y"), task_pars,
    n_subjects = 3, n_trials = 1,
    sds = c(kappa = 0.2, kappa_task1 = 0.4, thetat_task2 = 0.1),
    tasks = c("1", "2"), generator = echo_generator, seed = 1
  )
  expect_equal(sim$pars, c(
    kappa_task1 = log(8), kappa_task2 = log(4),
    thetat_task1 = 0.5, thetat_task2 = 0.5
  ))
  expect_equal(sim$sds, c(
    kappa_task1 = 0.4, kappa_task2 = 0.2,
    thetat_task1 = 0, thetat_task2 = 0.1
  ))
  expect_equal(sim$truth$population$term, names(sim$pars))
  expect_equal(
    sim$truth$sd$term, c("kappa_task1", "kappa_task2", "thetat_task2")
  )
  expect_setequal(
    unique(sim$truth$subjects$term),
    c("kappa_task1", "kappa_task2", "thetat_task2")
  )
  expect_setequal(sim$truth$cor$term, c(
    "kappa_task1__kappa_task2", "kappa_task1__thetat_task2",
    "kappa_task2__thetat_task2"
  ))
  expect_identical(sim$tasks, c("1", "2"))
  expect_identical(sim$task_col, "task")
  expect_output(print(sim), "; tasks: 1, 2")
  expect_output(print(sim), "kappa_task1")

  # a fixed parameter stays bare and reaches every task
  with_fixed <- simulate_recovery(
    bmm::mixture2p(resp_error = "y"), c(task_pars, mu1 = 1),
    n_subjects = 1, n_trials = 1, tasks = c("1", "2"),
    generator = function(pars, n_trials, model) data.frame(y = pars$mu1)
  )
  expect_equal(names(with_fixed$pars)[[5L]], "mu1")
  expect_equal(with_fixed$data$y, rep(inverse_link(1, "tan_half"), 2L))
})

test_that("the generator runs per subject and task with that task's values", {
  skip_if_not_installed("bmm")
  calls <- new.env()
  calls$pars <- list()
  recorder <- function(pars, n_trials, model) {
    calls$pars[[length(calls$pars) + 1L]] <- pars
    data.frame(y = rep(0, n_trials))
  }
  sim <- simulate_recovery(
    bmm::mixture2p(resp_error = "y"),
    c(kappa = log(8), kappa_condB = log(4), thetat = 0.5),
    n_subjects = 3, n_trials = 2, sds = c(kappa_condA = 0.3),
    tasks = c("A", "B"), task_col = "cond", generator = recorder, seed = 2
  )
  expect_length(calls$pars, 6L)
  expect_named(sim$data, c("id", "cond", "y"))
  expect_s3_class(sim$data$cond, "factor")
  expect_identical(levels(sim$data$cond), c("A", "B"))
  expect_equal(nrow(sim$data), 3L * 2L * 2L)
  # subject-major, task-minor
  expect_equal(as.character(sim$data$id), rep(c("1", "2", "3"), each = 4L))
  expect_equal(
    as.character(sim$data$cond), rep(rep(c("A", "B"), each = 2L), 3L)
  )

  subjects <- sim$truth$subjects
  for (i in 1:3) {
    a <- calls$pars[[2L * i - 1L]]
    b <- calls$pars[[2L * i]]
    kappa_a <- subjects$true_value[subjects$id == i]
    expect_equal(a$kappa, exp(kappa_a))
    expect_equal(b$kappa, 4)
    expect_equal(a$thetat, stats::plogis(0.5))
    expect_equal(b$mu1, 0)
    expect_false(any(grepl("_", names(a))))
  }
})

test_that("covariates sit between id and the task column", {
  skip_if_not_installed("bmm")
  sim <- simulate_recovery(
    bmm::mixture2p(resp_error = "y"), task_pars,
    n_subjects = 2, n_trials = 1, sds = c(kappa = 0.3),
    covariates = list(G = c(mean = 0, sd = 1)),
    tasks = c("1", "2"), generator = echo_generator, seed = 3
  )
  expect_named(sim$data, c("id", "G", "task", "y", "kappa_seen", "thetat_seen"))
  per_id <- tapply(sim$data$G, sim$data$id, function(g) length(unique(g)))
  expect_true(all(per_id == 1L))
  expect_true("G__kappa_task1" %in% sim$truth$cor$term)
})

test_that("correlated task terms reach their target in large samples", {
  skip_if_not_installed("bmm")
  cors <- diag(3)
  terms <- c("kappa_task1", "kappa_task2", "thetat_task1")
  dimnames(cors) <- list(terms, terms)
  cors[1, 2] <- cors[2, 1] <- 0.5
  cors[1, 3] <- cors[3, 1] <- -0.3
  sim <- simulate_recovery(
    bmm::mixture2p(resp_error = "y"), task_pars,
    n_subjects = 4000, n_trials = 1,
    sds = c(kappa = 0.3, thetat_task1 = 0.5), cors = cors,
    tasks = c("1", "2"),
    generator = function(pars, n_trials, model) data.frame(y = 0),
    seed = 4
  )
  wide <- subjects_wide(sim$truth$subjects)
  expect_lt(abs(stats::cor(wide$kappa_task1, wide$kappa_task2) - 0.5), 0.05)
  expect_lt(abs(stats::cor(wide$kappa_task1, wide$thetat_task1) + 0.3), 0.05)
  expect_equal(
    sim$truth$cor$true_value[sim$truth$cor$term == "kappa_task1__kappa_task2"],
    0.5
  )
})

test_that("realised task values passed back reproduce the simulation", {
  skip_if_not_installed("bmm")
  model <- bmm::mixture2p(resp_error = "y")
  first <- simulate_recovery(
    model, task_pars,
    n_subjects = 3, n_trials = 5, sds = c(kappa = 0.3),
    tasks = c("1", "2"), seed = 6
  )
  replay <- function() {
    simulate_recovery(
      model, first$pars,
      n_subjects = 3, n_trials = 5, sds = first$sds,
      subject_pars = first$truth$subjects,
      tasks = c("1", "2"), seed = 6
    )
  }
  again <- replay()
  expect_identical(again$pars, first$pars)
  expect_identical(again$sds, first$sds)
  expect_identical(again$truth, first$truth)
  expect_identical(replay()$data, again$data)

  # tasks shift the random stream relative to a task-free simulation: the
  # first varying term takes the same numbers, the responses do not
  none <- simulate_recovery(
    model, c(kappa = log(8), thetat = 0.5),
    n_subjects = 3, n_trials = 5, sds = c(kappa = 0.3), seed = 6
  )
  task1 <- first$truth$subjects$term == "kappa_task1"
  expect_identical(
    none$truth$subjects$true_value, first$truth$subjects$true_value[task1]
  )
  expect_false(identical(none$data$y, first$data$y[first$data$task == "1"]))
})

test_that("bad tasks, task_col and task terms are refused by name", {
  skip_if_not_installed("bmm")
  model <- bmm::mixture2p(resp_error = "y")
  sim <- function(pars = task_pars, ...) {
    simulate_recovery(
      model, pars,
      n_subjects = 2, n_trials = 1, generator = echo_generator, ...
    )
  }
  expect_error(sim(tasks = c(1, 2)), "tasks")
  expect_error(sim(tasks = "1"), "at least")
  expect_error(sim(tasks = c("1", NA)), "tasks")
  expect_error(sim(tasks = c("a_1", "b")), "a_1")
  expect_error(sim(tasks = c("x", "x")), "x")
  expect_error(sim(tasks = c("1", "2"), task_col = "1task"), "task_col")
  expect_error(sim(tasks = c("1", "2"), task_col = c("a", "b")), "task_col")
  expect_error(sim(tasks = c("1", "2"), task_col = "id"), "id")
  expect_error(sim(tasks = c("1", "2"), task_col = "y"), "y")
  expect_error(sim(tasks = c("1", "2"), task_col = "kappa"), "kappa")
  expect_error(sim(tasks = c("1", "2"), task_col = "mu1"), "mu1")
  expect_error(
    sim(
      tasks = c("1", "2"), task_col = "G",
      covariates = list(G = c(mean = 0, sd = 1))
    ),
    "G"
  )
  # task_col is not checked without tasks
  expect_no_error(sim(pars = c(kappa = 1, thetat = 0), task_col = "id"))

  expect_error(
    sim(pars = c(task_pars, kappa_task3 = 1), tasks = c("1", "2")),
    "kappa_task3"
  )
  expect_error(
    sim(pars = c(task_pars, mu1_task1 = 1), tasks = c("1", "2")),
    "mu1_task1"
  )
  expect_error(
    sim(tasks = c("1", "2"), sds = c(thetat_task9 = 0.1)),
    "thetat_task9"
  )
  expect_error(sim(tasks = c("1", "2"), sds = c(mu1 = 0.1)), "mu1")
  expect_error(
    sim(pars = c(kappa_task1 = 1, thetat = 0), tasks = c("1", "2")),
    "kappa_task2"
  )
  bare_cors <- matrix(c(1, 0.5, 0.5, 1), 2,
    dimnames = rep(list(c("kappa", "thetat")), 2)
  )
  err <- expect_error(
    sim(
      tasks = c("1", "2"), sds = c(kappa = 0.3, thetat = 0.2),
      cors = bare_cors
    ),
    "kappa"
  )
  expect_match(conditionMessage(err), "full terms")
  collide <- function(pars, n_trials, model) data.frame(y = 0, task = 1)
  expect_error(
    simulate_recovery(
      model, task_pars,
      n_subjects = 2, n_trials = 1, tasks = c("1", "2"), generator = collide
    ),
    "task column"
  )
})

# recovery_formula with tasks -------------------------------------------

test_that("task formulas are cell means with three random-effect structures", {
  skip_if_not_installed("bmm")
  skip_if_no_bmm_sdt()
  model <- bmm::mixture2p(resp_error = "y")
  f <- function(re_cor) {
    recovery_formula(model, re_cor = re_cor, task_col = "task")
  }
  expect_equal(
    deparse(f("none")$kappa), "kappa ~ 0 + task + (0 + task || id)"
  )
  expect_equal(
    deparse(f("within")$thetat), "thetat ~ 0 + task + (0 + task | id)"
  )
  expect_equal(
    deparse(f("all")$kappa), "kappa ~ 0 + task + (0 + task | p | id)"
  )

  expect_message(
    within <- recovery_formula(model, re_cor = "within"),
    "nothing to correlate within a parameter"
  )
  expect_equal(deparse(within$kappa), "kappa ~ 1 + (1 | id)")

  mafc <- bmm::sdt_mafc(response = "k", n_trials = "n", m = 4)
  expect_message(
    single <- recovery_formula(mafc, re_cor = "all", task_col = "cond"),
    "within"
  )
  expect_equal(deparse(single$d), "d ~ 0 + cond + (0 + cond | id)")

  expect_error(recovery_formula(model, task_col = "my_task"), "task_col")
  expect_error(recovery_formula(model, task_col = "my.task"), "task_col")
  expect_error(recovery_formula(model, task_col = 1), "task_col")
  expect_error(recovery_formula(model, task_col = c("a", "b")), "task_col")
})

test_that("contrast coding writes an intercept and a task slope", {
  skip_if_not_installed("bmm")
  skip_if_no_bmm_sdt()
  model <- bmm::mixture2p(resp_error = "y")
  f <- function(re_cor) {
    recovery_formula(
      model,
      re_cor = re_cor, task_col = "task", coding = "contrast"
    )
  }
  expect_equal(
    deparse(f("none")$kappa), "kappa ~ 1 + task + (1 + task || id)"
  )
  expect_equal(
    deparse(f("within")$thetat), "thetat ~ 1 + task + (1 + task | id)"
  )
  expect_equal(
    deparse(f("all")$kappa), "kappa ~ 1 + task + (1 + task | p | id)"
  )

  # cell means stay the default, untouched by the new argument
  expect_equal(
    deparse(recovery_formula(model, task_col = "task")$kappa),
    "kappa ~ 0 + task + (0 + task || id)"
  )

  expect_error(
    recovery_formula(model, coding = "contrast"), "task_col"
  )
  expect_error(
    recovery_formula(model, task_col = "task", coding = "sum"), "coding"
  )
})

test_that("task formulas pass the mock backend and match the truth terms", {
  skip_if_not_installed("bmm")
  skip_if_not_installed("brms")
  model <- bmm::mixture2p(resp_error = "y")
  sim <- simulate_recovery(
    model, task_pars,
    n_subjects = 4, n_trials = 20, sds = c(kappa = 0.3, thetat = 0.2),
    tasks = c("1", "2"), seed = 1
  )
  for (re_cor in c("none", "within", "all")) {
    formula <- recovery_formula(model, re_cor = re_cor, task_col = "task")
    fit <- mock_bmm(formula, sim$data, model)
    expect_s3_class(fit, "bmmfit")
    ranef <- as.data.frame(fit$ranef)
    expect_equal(
      paste0(ranef$nlpar, "_", ranef$coef),
      c("kappa_task1", "kappa_task2", "thetat_task1", "thetat_task2"),
      label = re_cor
    )
    expect_equal(all(ranef$cor), re_cor != "none", label = re_cor)
  }

  prior <- as.data.frame(bmm::default_prior(formula, sim$data, model))
  b <- prior[prior$class == "b" & nzchar(prior$coef), ]
  expect_setequal(paste0(b$nlpar, "_", b$coef), sim$truth$population$term)
})

# resolving bmm's fork-only functions ------------------------------------

test_that("bmm_fun() resolves an export and names one the build lacks", {
  skip_if_not_installed("bmm")
  # The signal-detection stack is in no released bmm, so the adapters that
  # need it resolve the name at call time. A literal bmm::rsdt_yn() makes
  # "checking dependencies in R code" report a missing object on every
  # machine with a released bmm --- a WARNING, which fails the check
  # workflow, because r-lib/actions checks with error-on = "warning".
  expect_true(is.function(bmm_fun("rmixture2p")))
  expect_error(bmm_fun("no_such_bmm_function"), "no_such_bmm_function")
})
