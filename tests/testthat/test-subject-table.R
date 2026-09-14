# Tests for subject_table(), written against
# local/dev/spec-milestone-5-correlation-recovery.md section 5.3. The mock
# fit and its extract_subject_draws() method come from helper-generate.R.

# A simulation carrying only what subject_table() reads: five subjects,
# kappa varying, thetat fixed at its population value, one covariate.
table_simulation <- function() {
  ids <- as.character(1:5)
  structure(
    list(
      truth = list(
        population = tibble::tibble(
          term = c("kappa", "thetat"), true_value = c(2, 0.5)
        ),
        subjects = tibble::tibble(
          id = ids, term = "kappa", true_value = c(1.5, 2.2, 1.9, 2.6, 1.1)
        ),
        covariates = tibble::tibble(
          id = ids, term = "G", true_value = c(-1, 0.3, 0.8, -0.2, 1.4)
        )
      ),
      covariates = list(G = c(mean = 0, sd = 1)),
      model = list(links = list(kappa = "log", thetat = "logit"))
    ),
    class = "bmmtools_simulation"
  )
}

table_fit <- function(ids = 1:5) {
  structure(
    list(parameters = c("kappa", "thetat"), ids = as.character(ids)),
    class = "mockfit"
  )
}

test_that("a simulation and its fit give one row per subject", {
  sim <- table_simulation()
  fit <- table_fit()
  out <- subject_table(sim, fit)

  expect_s3_class(out, "tbl_df")
  expect_named(out, c(
    "condition", "replication", "id", "true_kappa", "est_kappa",
    "true_thetat", "est_thetat", "G"
  ))
  expect_equal(nrow(out), 5L)
  expect_identical(out$condition, rep(NA_character_, 5L))
  expect_identical(out$replication, rep(1L, 5L))
  expect_identical(out$id, as.character(1:5))

  draws <- extract_subject_draws(fit)
  expect_equal(out$est_kappa, unname(colMeans(draws[, , , "kappa"], dims = 2L)))
  expect_equal(
    out$est_thetat,
    unname(colMeans(draws[, , , "thetat"], dims = 2L))
  )
  expect_equal(out$true_kappa, sim$truth$subjects$true_value)
  # a parameter that does not vary has its population value for everyone
  expect_equal(out$true_thetat, rep(0.5, 5L))
  expect_equal(out$G, sim$truth$covariates$true_value)
})

test_that("point = 'median' takes the posterior median", {
  fit <- table_fit()
  out <- subject_table(table_simulation(), fit, point = "median")
  draws <- extract_subject_draws(fit)
  expect_equal(
    out$est_kappa,
    unname(apply(draws[, , , "kappa"], 3L, stats::median))
  )
})

test_that("the natural scale back-transforms parameters, not covariates", {
  sim <- table_simulation()
  fit <- table_fit()
  link <- subject_table(sim, fit)
  natural <- subject_table(sim, fit, scale = "natural")

  expect_equal(natural$est_kappa, exp(link$est_kappa))
  expect_equal(natural$true_kappa, exp(link$true_kappa))
  expect_equal(natural$est_thetat, stats::plogis(link$est_thetat))
  expect_equal(natural$true_thetat, stats::plogis(link$true_thetat))
  expect_equal(natural$G, link$G)

  # links given explicitly win over the model's table
  identity <- subject_table(
    sim, fit,
    scale = "natural", links = c(kappa = "identity", thetat = "identity")
  )
  expect_equal(identity$est_kappa, link$est_kappa)

  # no links anywhere: the link scale, with a message
  sim$model <- NULL
  expect_message(
    fallback <- subject_table(sim, fit, scale = "natural"),
    "link scale"
  )
  expect_equal(fallback, link)
})

test_that("a simulation without covariates has no covariate columns", {
  sim <- table_simulation()
  sim$covariates <- NULL
  sim$truth$covariates <- sim$truth$covariates[0, ]
  out <- subject_table(sim, table_fit())
  expect_named(out, c(
    "condition", "replication", "id", "true_kappa", "est_kappa",
    "true_thetat", "est_thetat"
  ))
})

test_that("a grid result gives one row per condition, replication and id", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  out <- suppressMessages(recovery_grid(
    bmm::mixture2p(resp_error = "y"),
    grid = small_grid(),
    pars = c(kappa = log(8), thetat = stats::qlogis(0.75)),
    dir = dir, reps = 2L,
    sds = c(kappa = 0.3, thetat = 0.4),
    covariates = list(G = c(mean = 0, sd = 1)),
    seed = 100, preflight = FALSE,
    .fitter = grid_mock_fitter()$fitter
  ))

  tab <- subject_table(out)
  expect_named(tab, c(
    "condition", "replication", "id", "true_kappa", "est_kappa",
    "true_thetat", "est_thetat", "G"
  ))
  # (3 + 4) subjects in each of two replications
  expect_equal(nrow(tab), 14L)
  expect_equal(
    nrow(dplyr::distinct(tab[c("condition", "replication", "id")])),
    14L
  )

  sim <- readRDS(file.path(dir, "cell-2-rep-2-sim.rds"))
  rows <- tab[tab$condition == "row-2" & tab$replication == 2L, ]
  kappa <- sim$truth$subjects[sim$truth$subjects$term == "kappa", ]
  expect_equal(rows$true_kappa, kappa$true_value[match(rows$id, kappa$id)])
  expect_equal(
    rows$G,
    sim$truth$covariates$true_value[match(rows$id, sim$truth$covariates$id)]
  )
  draws <- extract_subject_draws(table_fit(1:4))
  expect_equal(
    rows$est_thetat,
    unname(colMeans(draws[, , rows$id, "thetat"], dims = 2L))
  )

  natural <- subject_table(out, scale = "natural")
  expect_equal(natural$est_kappa, exp(tab$est_kappa))
  expect_equal(natural$true_thetat, stats::plogis(tab$true_thetat))
  expect_equal(natural$G, tab$G)
})

test_that("subject_table refuses what it cannot read", {
  sim <- table_simulation()
  fit <- table_fit()
  expect_error(subject_table(sim), "fit")
  expect_error(subject_table(1:3), "bmmtools_simulation")
  expect_error(subject_table(sim, fit, point = "bogus"), "point")
  expect_error(subject_table(sim, fit, scale = "bogus"), "scale")
  expect_error(
    subject_table(sim, fit, scale = "natural", links = c("log", "logit")),
    "links"
  )
  expect_error(subject_table(sim, 1:3), "extract_subject_draws")

  recovery <- recover(
    fake_estimates("a", estimate = 1), fake_truth("a", true_value = 1),
    scale = "link"
  )
  expect_error(subject_table(recovery), "subject means")

  with_means <- recovery
  attr(with_means, "subject_means") <- tibble::tibble(
    condition = "row-1", replication = 1L, id = "1", term = "kappa",
    covariate = FALSE, mean = 1, median = 1, true_value = 1
  )
  expect_error(subject_table(with_means, fit), "fit")
  attr(with_means, "subject_means") <- attr(with_means, "subject_means")[0, ]
  expect_error(subject_table(with_means), "no subject means")
})

test_that("task terms give one column pair per parameter and task", {
  skip_if_not_installed("bmm")
  sim <- simulate_recovery(
    bmm::mixture2p(resp_error = "y"),
    c(kappa = 2, thetat = 0.5),
    n_subjects = 4, n_trials = 3, sds = c(kappa = 0.3),
    tasks = c("1", "2"), seed = 1
  )
  fit <- grid_mock_fitter()$fitter(
    recovery_formula(sim$model, task_col = "task"), sim$data, sim$model
  )
  expect_equal(
    dimnames(extract_subject_draws(fit))$term,
    c("kappa_task1", "kappa_task2", "thetat_task1", "thetat_task2")
  )
  out <- subject_table(sim, fit)
  expect_named(out, c(
    "condition", "replication", "id",
    "true_kappa_task1", "est_kappa_task1", "true_kappa_task2",
    "est_kappa_task2", "true_thetat_task1", "est_thetat_task1",
    "true_thetat_task2", "est_thetat_task2"
  ))
  subjects <- sim$truth$subjects
  expect_equal(
    out$true_kappa_task2,
    subjects$true_value[subjects$term == "kappa_task2"]
  )
  expect_equal(out$true_thetat_task1, rep(0.5, 4L))

  natural <- subject_table(sim, fit, scale = "natural")
  expect_equal(natural$true_kappa_task1, exp(out$true_kappa_task1))

  cors <- extract_correlations(fit, estimator = "draws")
  expect_true("kappa_task1__kappa_task2" %in% cors$term)
})
