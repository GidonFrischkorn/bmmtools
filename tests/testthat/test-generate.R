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
