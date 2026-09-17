# Tests for components and cross-model correlations with separate fits,
# written against local/dev/spec-milestone-5-correlation-recovery.md
# section 5.5 as amended for 5.5a. Nothing here compiles Stan: fits come
# from the grid's mock fitter (helper-generate.R) or from `array_fit()`
# (helper-components.R).

m2p <- function() bmm::mixture2p(resp_error = "y")
sdm_model <- function() bmm::sdm(resp_error = "y")

m2p_component <- function(name = "a", ...) {
  recovery_component(
    m2p(),
    pars = c(kappa = log(8), thetat = stats::qlogis(0.75)),
    n_trials = 5, sds = c(kappa = 0.3, thetat = 0.5), name = name, ...
  )
}

sdm_component <- function(name = "b", ...) {
  recovery_component(
    sdm_model(),
    pars = c(c = log(3), kappa = log(5)),
    n_trials = 4, sds = c(c = 0.4), name = name, ...
  )
}

# recovery_component() ------------------------------------------------------

test_that("a component stores what it was given and prints one line", {
  skip_if_not_installed("bmm")
  comp <- m2p_component(generator = zero_generator)
  expect_s3_class(comp, "bmmtools_component")
  expect_equal(comp$name, "a")
  expect_equal(comp$n_trials, 5L)
  expect_null(comp$tasks)
  expect_null(comp$formula)
  expect_identical(comp$generator, zero_generator)
  out <- capture.output(print(comp))
  expect_length(out, 1L)
  expect_match(out, "<bmmtools_component> a: mixture2p", fixed = TRUE)
  expect_invisible(print(comp))

  tasked <- recovery_component(
    m2p(), c(kappa = 2, thetat = 0),
    n_trials = 3,
    tasks = c("1", "2"), task_col = "cond", name = "t"
  )
  expect_equal(tasked$tasks, c("1", "2"))
  expect_equal(tasked$task_col, "cond")
  expect_match(capture.output(print(tasked)), "tasks: 1, 2", fixed = TRUE)
})

test_that("component names must be syntactic without _ or .", {
  skip_if_not_installed("bmm")
  for (bad in list("a_b", "a.b", "", "1a", c("a", "b"), 1, NA_character_)) {
    expect_error(m2p_component(name = bad), "name", label = deparse(bad))
  }
  expect_error(m2p_component(name = "a_b"), "ab", fixed = TRUE)
  expect_error(
    recovery_component(m2p(), c(kappa = 1, thetat = 0), n_trials = 2),
    "name"
  )
  expect_s3_class(m2p_component(name = "m3"), "bmmtools_component")
})

test_that("a component validates what needs no draw", {
  skip_if_not_installed("bmm")
  pars <- c(kappa = 1, thetat = 0)
  expect_error(
    recovery_component(list(), pars, n_trials = 2, name = "a"),
    "bmmodel"
  )
  expect_error(
    recovery_component(m2p(), pars, n_trials = 0, name = "a"),
    "n_trials"
  )
  expect_error(
    recovery_component(m2p(), pars, n_trials = 2, tasks = "1", name = "a"),
    "tasks"
  )
  expect_error(
    recovery_component(
      structure(list(), class = c("bmmodel", "mpt")), c(a = 1),
      n_trials = 2, name = "a"
    ),
    "No generator"
  )
  expect_error(
    recovery_component(m2p(), pars, n_trials = 2, generator = 1, name = "a"),
    "generator"
  )
  expect_error(
    recovery_component(m2p(), pars, n_trials = 2, formula = y ~ 1, name = "a"),
    "bmmformula"
  )
  expect_error(
    recovery_component(m2p(), c(kappa = 1), n_trials = 2, name = "a"),
    "missing"
  )
  expect_error(
    recovery_component(
      m2p(), pars,
      n_trials = 2, sds = c(bogus = 1), name = "a"
    ),
    "bogus"
  )
  expect_error(
    recovery_component(
      m2p(), function() pars,
      n_trials = 2, sds = c(mu1 = 1), name = "a"
    ),
    "fixes"
  )
  expect_error(
    recovery_component(m2p(), function(x) pars, n_trials = 2, name = "a"),
    "no arguments"
  )
  # a task term that does not exist is caught before any draw
  expect_error(
    recovery_component(
      m2p(), c(pars, kappa_task3 = 1),
      n_trials = 2,
      tasks = c("1", "2"), name = "a"
    ),
    "kappa_task3"
  )
  expect_s3_class(
    recovery_component(
      m2p(), function() pars,
      n_trials = 2, sds = function() c(kappa = 1),
      formula = bmm::bmf(kappa ~ 1, thetat ~ 1), name = "a"
    ),
    "bmmtools_component"
  )
})

# simulate_components() -----------------------------------------------------

test_that("the joint truth carries prefixed terms, tasks included", {
  skip_if_not_installed("bmm")
  a <- recovery_component(
    m2p(), c(kappa = 2, thetat = 0),
    n_trials = 2,
    sds = c(kappa = 0.3), tasks = c("1", "2"), name = "a",
    generator = zero_generator
  )
  b <- sdm_component(generator = zero_generator)
  cors <- diag(3)
  dimnames(cors) <- rep(list(c("a_kappa_task1", "b_c", "G")), 2)
  cors["a_kappa_task1", "b_c"] <- cors["b_c", "a_kappa_task1"] <- 0.4
  set <- simulate_components(
    list(a, b),
    n_subjects = 6, cors = cors,
    covariates = list(G = c(mean = 0, sd = 1)), seed = 3
  )

  expect_s3_class(set, "bmmtools_simulation_set")
  expect_named(set$components, c("a", "b"))
  expect_s3_class(set$components$a, "bmmtools_simulation")
  expect_equal(
    set$truth$population$term,
    c(
      "a_kappa_task1", "a_kappa_task2", "a_thetat_task1", "a_thetat_task2",
      "b_c", "b_kappa"
    )
  )
  varying <- c("a_kappa_task1", "a_kappa_task2", "b_c")
  expect_equal(unique(set$truth$subjects$term), varying)
  expect_equal(set$truth$sd$term, varying)
  expect_equal(names(set$sds)[set$sds > 0], varying)
  expect_equal(unique(set$truth$covariates$term), "G")

  # every pair of varying terms and covariates, zeros included, D27 names
  expect_equal(nrow(set$truth$cor), choose(4, 2))
  expect_true(all(c(
    "a_kappa_task1__b_c", "G__a_kappa_task1", "G__b_c",
    "a_kappa_task1__a_kappa_task2"
  ) %in% set$truth$cor$term))
  expect_equal(
    set$truth$cor$true_value[set$truth$cor$term == "a_kappa_task1__b_c"], 0.4
  )
  expect_equal(dimnames(set$cors)[[1L]], c(varying, "G"))

  # the component keeps its own terms, unprefixed, and no covariates
  expect_equal(unique(set$components$a$truth$subjects$term), c(
    "kappa_task1", "kappa_task2"
  ))
  expect_null(set$components$a$covariates)
  expect_equal(set$components$a$tasks, c("1", "2"))
  expect_null(set$components$b$cors)

  # the covariate data are one row per subject, id a factor
  expect_named(set$covariate_data, c("id", "G"))
  expect_s3_class(set$covariate_data$id, "factor")
  expect_type(set$covariate_data$G, "double")
  expect_equal(nrow(set$covariate_data), 6L)
  expect_equal(
    set$covariate_data$G,
    set$truth$covariates$true_value
  )

  # links: every model's table, prefixed, resolved by the longest prefix
  expect_equal(set$links[["a_kappa"]], "log")
  expect_equal(set$links[["b_c"]], "log")
  expect_equal(
    link_of(c("a_kappa_task1", "a_thetat_task2", "b_c", "G"), set$links),
    c("log", "logit", "log", "identity")
  )
  expect_equal(set$n_subjects, 6L)
  expect_equal(set$seed, 3)

  out <- capture.output(print(set))
  expect_length(out, 1L)
  expect_match(out, "a (mixture2p, 2 trials)", fixed = TRUE)
  expect_match(out, "covariates: G", fixed = TRUE)
  expect_match(out, "a_kappa_task1__b_c = 0.4", fixed = TRUE)
})

test_that("components share one subject draw across models", {
  skip_if_not_installed("bmm")
  a <- recovery_component(
    m2p(), c(kappa = 2, thetat = 0),
    n_trials = 1,
    sds = c(kappa = 0.3), generator = zero_generator, name = "a"
  )
  b <- recovery_component(
    sdm_model(), c(c = 1, kappa = 1),
    n_trials = 1,
    sds = c(c = 0.5), generator = zero_generator, name = "b"
  )
  cors <- diag(3)
  dimnames(cors) <- rep(list(c("a_kappa", "b_c", "G")), 2)
  cors["a_kappa", "b_c"] <- cors["b_c", "a_kappa"] <- 0.6
  cors["G", "b_c"] <- cors["b_c", "G"] <- 0.3
  set <- simulate_components(
    list(a = a, b = b),
    n_subjects = 4000, cors = cors,
    covariates = list(G = c(mean = 0, sd = 1)), seed = 12
  )
  wide <- subjects_wide(set$truth$subjects)
  expect_lt(abs(stats::cor(wide$a_kappa, wide$b_c) - 0.6), 0.05)
  expect_lt(
    abs(stats::cor(set$covariate_data$G, wide$b_c) - 0.3), 0.05
  )
  # the components' own truths are the same numbers
  expect_identical(
    set$components$a$truth$subjects$true_value, wide$a_kappa
  )
  expect_identical(
    set$components$b$truth$subjects$true_value, wide$b_c
  )
})

test_that("one component and no covariates equals simulate_recovery()", {
  skip_if_not_installed("bmm")
  model <- m2p()
  pars <- function() c(kappa = stats::rnorm(1, 2, 0.1), thetat = 0.5)
  sds <- c(kappa = 0.3, thetat = 0.5)
  cors <- matrix(c(1, 0.4, 0.4, 1), 2, dimnames = rep(list(c(
    "kappa", "thetat"
  )), 2))
  direct <- simulate_recovery(
    model, pars,
    n_subjects = 7, n_trials = 3, sds = sds, cors = cors,
    generator = noisy_generator, seed = 5
  )
  component <- recovery_component(
    model, pars,
    n_trials = 3, sds = sds, generator = noisy_generator, name = "m"
  )
  prefixed <- cors
  dimnames(prefixed) <- rep(list(c("m_kappa", "m_thetat")), 2)
  set <- simulate_components(
    list(component),
    n_subjects = 7, cors = prefixed, seed = 5
  )
  inner <- set$components$m
  expect_true(is.na(inner$seed))
  expect_equal(direct$seed, 5)
  expect_identical(inner$generator, direct$generator)
  fields <- setdiff(names(unclass(direct)), "seed")
  expect_identical(unclass(inner)[fields], unclass(direct)[fields])

  # the adapter path, with tasks and no correlations
  adapter <- simulate_recovery(
    model, c(kappa = 2, thetat = 0.5),
    n_subjects = 3, n_trials = 4, sds = c(kappa = 0.2),
    tasks = c("x", "y"), seed = 8
  )
  set2 <- simulate_components(
    list(recovery_component(
      model, c(kappa = 2, thetat = 0.5),
      n_trials = 4, sds = c(kappa = 0.2), tasks = c("x", "y"), name = "m"
    )),
    n_subjects = 3, seed = 8
  )
  expect_identical(
    unclass(set2$components$m)[fields], unclass(adapter)[fields]
  )
  expect_equal(set2$components$m$generator, "adapter:mixture2p")
})

test_that("a second component leaves the first one's subject values", {
  skip_if_not_installed("bmm")
  a <- recovery_component(
    m2p(), c(kappa = 2, thetat = 0),
    n_trials = 2,
    sds = c(kappa = 0.3, thetat = 0.2), generator = noisy_generator,
    name = "a"
  )
  b <- recovery_component(
    sdm_model(), c(c = 1, kappa = 1),
    n_trials = 2,
    sds = c(c = 0.5, kappa = 0.1), generator = noisy_generator, name = "b"
  )
  one_cors <- matrix(c(1, 0.5, 0.5, 1), 2,
    dimnames = rep(list(c("a_kappa", "a_thetat")), 2)
  )
  two_cors <- diag(4)
  dimnames(two_cors) <- rep(list(c("a_kappa", "a_thetat", "b_c", "b_kappa")), 2)
  two_cors[1:2, 1:2] <- one_cors
  alone <- simulate_components(list(a), 20, cors = one_cors, seed = 9)
  both <- simulate_components(list(a, b), 20, cors = two_cors, seed = 9)
  # the Cholesky factor is upper triangular, so appended terms leave the
  # earlier columns of the draw alone
  expect_equal(
    both$components$a$truth$subjects, alone$components$a$truth$subjects
  )
  # the generator runs after the whole draw, whose extra columns take
  # random numbers first, so the responses do not survive
  expect_false(identical(both$components$a$data$y, alone$components$a$data$y))
})

test_that("covariates are drawn after every generator", {
  skip_if_not_installed("bmm")
  a <- recovery_component(
    m2p(), c(kappa = 2, thetat = 0),
    n_trials = 3,
    sds = c(kappa = 0.3), generator = noisy_generator, name = "a"
  )
  b <- recovery_component(
    sdm_model(), c(c = 1, kappa = 1),
    n_trials = 3,
    sds = c(c = 0.5), generator = noisy_generator, name = "b"
  )
  cors <- diag(3)
  dimnames(cors) <- rep(list(c("a_kappa", "b_c", "G")), 2)
  cors["G", "a_kappa"] <- cors["a_kappa", "G"] <- 0.5
  plain <- simulate_components(list(a, b), 10, seed = 2)
  with_g <- simulate_components(
    list(a, b), 10,
    cors = cors, covariates = list(G = c(mean = 1, sd = 2)), seed = 2
  )
  expect_identical(with_g$components$b$data, plain$components$b$data)
  expect_identical(with_g$truth$subjects, plain$truth$subjects)
  expect_null(plain$covariates)
  expect_named(plain$covariate_data, "id")
  expect_equal(nrow(plain$truth$covariates), 0L)
})

test_that("simulate_components refuses bad sets", {
  skip_if_not_installed("bmm")
  a <- m2p_component(generator = zero_generator)
  b <- sdm_component(generator = zero_generator)
  expect_error(simulate_components(a, 5), "list")
  expect_error(simulate_components(list(), 5), "list")
  expect_error(simulate_components(list(a, "b"), 5), "bmmtools_component")
  expect_error(simulate_components(list(a, a), 5), "more than once")
  expect_error(
    simulate_components(list(x = a, b = b), 5),
    "names"
  )
  expect_error(simulate_components(list(a, b), 0), "n_subjects")
  expect_error(simulate_components(list(a, b), 5, seed = "x"), "seed")
  expect_error(
    simulate_components(
      list(a, b), 5,
      covariates = list(a = c(mean = 0, sd = 1))
    ),
    "component name"
  )
  expect_error(
    simulate_components(
      list(a, b), 5,
      covariates = list(kappa = c(mean = 0, sd = 1))
    ),
    "cannot be used"
  )
  tasked <- recovery_component(
    m2p(), c(kappa = 2, thetat = 0),
    n_trials = 2, sds = c(kappa = 0.3),
    tasks = c("1", "2"), task_col = "G", name = "t"
  )
  expect_error(
    simulate_components(
      list(tasked, b), 5,
      covariates = list(G = c(mean = 0, sd = 1))
    ),
    "task_col"
  )
  bare <- diag(2)
  dimnames(bare) <- rep(list(c("t_kappa", "b_c")), 2)
  expect_error(
    simulate_components(list(tasked, b), 5, cors = bare),
    "t_kappa_G1"
  )
  unprefixed <- diag(2)
  dimnames(unprefixed) <- rep(list(c("kappa", "c")), 2)
  expect_error(simulate_components(list(a, b), 5, cors = unprefixed), "a_kappa")
  # an error inside a component names the component
  broken <- recovery_component(
    m2p(), function() c(kappa = 1),
    n_trials = 2, name = "broken"
  )
  expect_error(simulate_components(list(broken), 5), "broken")
})

# fit_components() ----------------------------------------------------------

two_component_set <- function(n_subjects = 8, seed = 4, ...) {
  a <- recovery_component(
    m2p(), c(kappa = 2, thetat = 0),
    n_trials = 2,
    sds = c(kappa = 0.3, thetat = 0.2), generator = zero_generator,
    name = "a"
  )
  b <- recovery_component(
    sdm_model(), c(c = 1, kappa = 1),
    n_trials = 2,
    sds = c(c_task1 = 0.5), generator = zero_generator, name = "b",
    tasks = c("1", "2")
  )
  cors <- diag(4)
  dimnames(cors) <- rep(list(c("a_kappa", "a_thetat", "b_c_task1", "G")), 2)
  cors["a_kappa", "b_c_task1"] <- cors["b_c_task1", "a_kappa"] <- 0.5
  simulate_components(
    list(a, b), n_subjects,
    cors = cors, covariates = list(G = c(mean = 0, sd = 1)), seed = seed, ...
  )
}

test_that("fit_components caches one fit per component and reuses them", {
  skip_if_not_installed("bmm")
  set <- two_component_set()
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  prior <- list(b = "a prior for b")

  fits <- suppressMessages(fit_components(
    set, dir,
    prior = prior, chains = 2, .fitter = mock$fitter
  ))
  expect_s3_class(fits, "bmmtools_fit_set")
  expect_named(fits, c("a", "b"))
  expect_equal(mock$calls$n, 2L)
  expect_true(file.exists(file.path(dir, "a.rds")))
  expect_true(file.exists(file.path(dir, "b.rds")))
  expect_true(file.exists(file.path(dir, "b.key")))
  # the set's seed and the dots reach the fitter
  expect_equal(mock$calls$log[[1L]]$seed, 4)
  expect_equal(mock$calls$log[[2L]]$chains, 2)
  # the default formula of a tasked component has cell means
  expect_equal(
    fits$b$parameters,
    c("c_task1", "c_task2", "kappa_task1", "kappa_task2")
  )
  expect_equal(fits$a$parameters, c("kappa", "thetat"))

  again <- suppressMessages(fit_components(
    set, dir,
    prior = prior, chains = 2, .fitter = mock$fitter
  ))
  expect_equal(mock$calls$n, 2L)
  expect_true(attr(again$a, "bmmtools_cache")$reused)

  out <- capture.output(print(fits))
  expect_length(out, 3L)
  expect_match(out[[2L]], "a:", fixed = TRUE)
  expect_invisible(print(fits))
})

test_that("per-component prior, seed and refit reach fit_cached", {
  skip_if_not_installed("bmm")
  set <- two_component_set()
  dir <- withr::local_tempdir()
  seen <- new.env()
  seen$calls <- list()
  fitter <- function(formula, data, model, prior = NULL, ...) {
    seen$calls[[length(seen$calls) + 1L]] <- list(
      prior = prior, dots = list(...)
    )
    structure(list(parameters = "x", ids = levels(data$id)), class = "mockfit")
  }
  suppressMessages(fit_components(
    set, dir,
    prior = list(a = "pa"), seed = 99, .fitter = fitter
  ))
  expect_equal(seen$calls[[1L]]$prior, "pa")
  expect_null(seen$calls[[2L]]$prior)
  expect_equal(seen$calls[[1L]]$dots$seed, 99)
  suppressMessages(fit_components(
    set, dir,
    refit = "always", seed = 99, .fitter = fitter
  ))
  expect_length(seen$calls, 4L)

  unseeded <- set
  unseeded$seed <- NA_real_
  suppressMessages(fit_components(
    unseeded, withr::local_tempdir(),
    .fitter = fitter
  ))
  expect_null(seen$calls[[5L]]$dots$seed)

  custom <- recovery_component(
    m2p(), c(kappa = 2, thetat = 0),
    n_trials = 2, sds = c(kappa = 0.3),
    generator = zero_generator, name = "c",
    formula = bmm::bmf(kappa ~ 1 + (1 | id), thetat ~ 1)
  )
  formulas <- new.env()
  formulas$seen <- list()
  capture <- function(formula, data, model, prior = NULL, ...) {
    formulas$seen[[length(formulas$seen) + 1L]] <- formula
    structure(list(parameters = "x", ids = levels(data$id)), class = "mockfit")
  }
  one <- simulate_components(list(custom), 4, seed = 1)
  suppressMessages(fit_components(
    one, withr::local_tempdir(),
    .fitter = capture
  ))
  expect_identical(formulas$seen[[1L]], custom$formula)
})

test_that("fit_components validates its arguments", {
  skip_if_not_installed("bmm")
  set <- two_component_set()
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  expect_error(fit_components(list(), dir), "bmmtools_simulation_set")
  expect_error(fit_components(set, 1), "dir")
  expect_error(
    fit_components(set, dir, prior = data.frame(prior = "x")),
    "named list"
  )
  expect_error(
    fit_components(set, dir, prior = structure(list(), class = "brmsprior")),
    "named list"
  )
  expect_error(fit_components(set, dir, prior = list("x")), "named list")
  expect_error(
    fit_components(set, dir, prior = list(zz = "x")),
    "zz"
  )
  expect_error(
    fit_components(set, dir, file = "x", .fitter = mock$fitter),
    "file"
  )
  expect_equal(mock$calls$n, 0L)
})

# extract_correlations() on a set --------------------------------------------

test_that("cross-fit draws and point match a hand-bound oracle", {
  a <- fake_subject_draws(
    n_subjects = 7, terms = c("kappa", "thetat"), seed = 1
  )
  b <- fake_subject_draws(n_subjects = 7, terms = "c", seed = 2)
  fits <- list(m1 = array_fit(a), m2 = array_fit(b))
  out <- extract_correlations(fits, estimator = c("draws", "point"))

  expect_setequal(out$term, c(
    "m1_kappa__m1_thetat", "m1_kappa__m2_c", "m1_thetat__m2_c"
  ))
  bound <- array(
    c(a, b),
    dim = c(40, 2, 7, 3),
    dimnames = list(NULL, NULL, NULL, c("m1_kappa", "m1_thetat", "m2_c"))
  )
  per_draw <- apply(bound, 1:2, function(s) {
    stats::cor(s[, "m1_kappa"], s[, "m2_c"])
  })
  row <- out[out$estimator == "draws" & out$term == "m1_kappa__m2_c", ]
  expect_equal(row$estimate, stats::median(per_draw))
  expect_equal(
    row$ci_low,
    stats::quantile(per_draw, 0.025, names = FALSE)
  )
  expect_equal(row$var1, "m1_kappa")
  expect_equal(row$n, 7L)

  means <- colMeans(bound, dims = 2L)
  test <- stats::cor.test(means[, "m1_thetat"], means[, "m2_c"])
  point <- out[out$estimator == "point" & out$term == "m1_thetat__m2_c", ]
  expect_equal(point$estimate, unname(test$estimate))
  expect_equal(point$ci_low, test$conf.int[[1L]])
  expect_equal(point$ci_high, test$conf.int[[2L]])

  # the within-fit pair is what the fit alone gives
  alone <- extract_correlations(array_fit(a), estimator = "draws")
  within <- out[out$estimator == "draws" & out$term == "m1_kappa__m1_thetat", ]
  expect_equal(within$estimate, alone$estimate)

  # a fit set is read the same way
  set <- structure(fits, class = c("bmmtools_fit_set", "list"))
  expect_equal(extract_correlations(set, estimator = "point"), out[
    out$estimator == "point",
  ], ignore_attr = TRUE)
})

test_that("unequal draw counts are trimmed with one message", {
  a <- fake_subject_draws(
    n_subjects = 5, terms = "kappa", n_iter = 40, n_chain = 2
  )
  b <- fake_subject_draws(
    n_subjects = 5, terms = "c", n_iter = 30, n_chain = 3, seed = 2
  )
  fits <- list(a = array_fit(a), b = array_fit(b))
  msgs <- testthat::capture_messages(
    out <- extract_correlations(fits, estimator = "draws")
  )
  expect_length(msgs, 1L)
  expect_match(msgs, "30 iterations")
  expect_match(msgs, "2 chains")
  bound <- array(
    c(a[1:30, 1:2, , ], b[1:30, 1:2, , ]),
    dim = c(30, 2, 5, 2)
  )
  per_draw <- apply(bound, 1:2, function(s) stats::cor(s[, 1], s[, 2]))
  expect_equal(out$estimate, stats::median(per_draw))

  # equal counts say nothing
  expect_silent(extract_correlations(
    list(a = array_fit(a), b = array_fit(a[, , , 1, drop = FALSE])),
    estimator = "point"
  ))
})

test_that("subject ids are aligned across fits and must agree", {
  a <- fake_subject_draws(n_subjects = 5, terms = "kappa")
  b <- fake_subject_draws(n_subjects = 5, terms = "c", seed = 3)
  shuffled <- b[, , c(5, 3, 1, 2, 4), , drop = FALSE]
  aligned <- extract_correlations(
    list(a = array_fit(a), b = array_fit(b)),
    estimator = "point"
  )
  reordered <- extract_correlations(
    list(a = array_fit(a), b = array_fit(shuffled)),
    estimator = "point"
  )
  expect_equal(reordered$estimate, aligned$estimate)

  other <- fake_subject_draws(n_subjects = 6, terms = "c")
  expect_error(
    extract_correlations(
      list(a = array_fit(a), b = array_fit(other)),
      estimator = "point"
    ),
    "a.*b"
  )
})

test_that("model rows are per fit, prefixed, and converged joins the fits", {
  cor_rows <- tibble::tibble(
    term = "thetat__kappa", estimate = 0.3, ci_low = 0.1, ci_high = 0.5,
    rhat = 1, ess_bulk = 500, ess_tail = 500
  )
  a <- fake_subject_draws(n_subjects = 5, terms = c("kappa", "thetat"))
  b <- fake_subject_draws(n_subjects = 5, terms = "c", seed = 2)
  fits <- list(a = array_fit(a, cor = cor_rows), b = array_fit(b))
  covariates <- data.frame(id = as.character(1:5), G = c(1, 3, 2, 5, 4))

  out <- extract_correlations(
    fits,
    covariates = covariates,
    converged = c(a = TRUE, b = FALSE)
  )
  model <- out[out$estimator == "model", ]
  expect_equal(model$term, "a_kappa__a_thetat")
  expect_equal(model$estimate, 0.3)
  expect_equal(model$n, 5L)
  conv <- stats::setNames(out$converged, paste(out$estimator, out$term))
  expect_true(conv[["model a_kappa__a_thetat"]])
  expect_false(conv[["draws a_kappa__b_c"]])
  expect_true(conv[["point G__a_kappa"]])
  expect_false(conv[["point G__b_c"]])
  expect_true(all(c("G__a_kappa", "G__b_c") %in% out$term))

  one <- extract_correlations(fits, estimator = "point", converged = FALSE)
  expect_false(any(one$converged))
  default <- extract_correlations(fits, estimator = "point")
  expect_true(all(is.na(default$converged)))

  expect_error(
    extract_correlations(fits, converged = c(a = TRUE)),
    "converged"
  )
  expect_error(
    extract_correlations(fits, converged = c(TRUE, FALSE)),
    "converged"
  )
  expect_error(
    extract_correlations(fits, covariates = "G"),
    "data frame"
  )
})

test_that("a set's links are every fit's table, prefixed", {
  a <- fake_subject_draws(n_subjects = 5, terms = "kappa")
  b <- fake_subject_draws(n_subjects = 5, terms = "c", seed = 2)
  fa <- array_fit(a)
  fb <- array_fit(b)
  natural <- extract_correlations(
    list(a = fa, b = fb),
    estimator = "point", scale = "natural",
    links = c(a_kappa = "log", b_c = "logit")
  )
  means <- colMeans(a, dims = 2L)[, "kappa"]
  means_b <- colMeans(b, dims = 2L)[, "c"]
  expect_equal(
    natural$estimate,
    stats::cor(exp(means), stats::plogis(means_b))
  )
  expect_message(
    fallback <- extract_correlations(
      list(a = fa, b = fb),
      estimator = "point", scale = "natural"
    ),
    "link scale"
  )
  expect_equal(fallback$scale, "link")

  # a bmm fit's table is read and prefixed
  class(fa) <- c("brmsfit", "arrayfit")
  fa$bmm <- list(model = list(links = list(kappa = "log")))
  expect_equal(set_model_links(list(a = fa, b = fb)), c(a_kappa = "log"))
  expect_null(set_model_links(list(b = fb)))
})

test_that("a plain list of fits needs valid names to be a set", {
  a <- array_fit(fake_subject_draws(n_subjects = 4, terms = "kappa"))
  expect_error(extract_correlations(list(a, a)), "brmsfit")
  expect_error(extract_correlations(list(x_1 = a, b = a)), "component")
  expect_error(extract_correlations(list(a = a, a = a)), "component")
})

# recover_correlations() and subject_table() on a set ------------------------

test_that("simulate, fit and score a set end to end with the mock", {
  skip_if_not_installed("bmm")
  set <- two_component_set()
  fits <- suppressMessages(fit_components(
    set, withr::local_tempdir(),
    .fitter = grid_mock_fitter()$fitter
  ))
  out <- recover_correlations(fits, set, estimator = c("draws", "point"))
  expect_s3_class(out, "bmmtools_cor_recovery")
  expect_true("a_kappa__b_c_task1" %in% out$term)
  expect_true("G__b_c_task1" %in% out$term)
  expect_equal(unique(out$replication), 1L)
  row <- out[out$term == "a_kappa__b_c_task1" & out$estimator == "point", ]
  expect_equal(row$true_value, 0.5)
  wide <- subjects_wide(set$truth$subjects)
  expect_equal(row$sample_value, stats::cor(wide$a_kappa, wide$b_c_task1))
  g_row <- out[out$term == "G__b_c_task1" & out$estimator == "point", ]
  means <- colMeans(extract_subject_draws(fits$b), dims = 2L)
  expect_equal(
    g_row$estimate,
    stats::cor(set$covariate_data$G, means[, "c_task1"])
  )
  summ <- summary(out)
  expect_s3_class(summ, "bmmtools_cor_recovery_summary")
  expect_equal(nrow(summ), nrow(out))

  # the same with a plain named list, and with the natural scale from the
  # set's links
  plain <- recover_correlations(unclass(fits), set, estimator = "point")
  expect_equal(plain$estimate, out$estimate[out$estimator == "point"])
  natural <- recover_correlations(
    fits, set,
    estimator = "point", scale = "natural"
  )
  expect_equal(unique(natural$scale), "natural")

  # an extracted tibble scores against the set
  extracted <- extract_correlations(
    fits,
    estimator = "point",
    covariates = as.data.frame(set$covariate_data)
  )
  from_tibble <- recover_correlations(extracted, set)
  expect_equal(from_tibble$sample_value, plain$sample_value)

  table <- subject_table(set, fits)
  expect_equal(nrow(table), 8L)
  terms <- c(
    "a_kappa", "a_thetat", "b_c_task1", "b_c_task2", "b_kappa_task1",
    "b_kappa_task2"
  )
  expect_named(table, c(
    "condition", "replication", "id",
    as.vector(rbind(paste0("true_", terms), paste0("est_", terms))),
    "G"
  ))
  expect_identical(table$condition, rep(NA_character_, 8L))
  expect_identical(table$replication, rep(1L, 8L))
  expect_equal(table$true_b_c_task1, wide$b_c_task1)
  expect_equal(
    table$est_a_kappa,
    unname(colMeans(extract_subject_draws(fits$a)[, , , "kappa"], dims = 2L))
  )
  # c_task2 does not vary: its truth is the population value
  expect_equal(table$true_b_c_task2, rep(1, 8L))
  expect_equal(table$G, set$covariate_data$G)
})

test_that("replications of sets pair by position", {
  skip_if_not_installed("bmm")
  sets <- list(two_component_set(seed = 1), two_component_set(seed = 2))
  mock <- grid_mock_fitter()
  fits <- lapply(sets, function(s) {
    suppressMessages(fit_components(
      s, withr::local_tempdir(),
      .fitter = mock$fitter
    ))
  })
  out <- recover_correlations(fits, sets, estimator = "point")
  expect_equal(sort(unique(out$replication)), 1:2)
  first <- recover_correlations(fits[[1L]], sets[[1L]], estimator = "point")
  expect_equal(
    out$sample_value[out$replication == 1L],
    first$sample_value
  )
  expect_error(
    recover_correlations(fits, sets[1L], estimator = "point"),
    "2"
  )
})

test_that("sets and single simulations do not mix", {
  skip_if_not_installed("bmm")
  set <- two_component_set()
  fits <- suppressMessages(fit_components(
    set, withr::local_tempdir(),
    .fitter = grid_mock_fitter()$fitter
  ))
  expect_error(
    recover_correlations(fits, set$components$a),
    "bmmtools_simulation_set"
  )
  expect_error(
    recover_correlations(fits$a, set),
    "set"
  )
  renamed <- fits
  names(renamed) <- c("a", "z")
  expect_error(recover_correlations(unclass(renamed), set), "z")
  expect_error(subject_table(set), "fit")
  expect_error(subject_table(set, fits$a), "set")
  expect_error(subject_table(set, unclass(renamed)), "z")
})
