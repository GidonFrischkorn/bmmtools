# The combined design (Milestone 12.3, spec section 12.3): population
# intercepts, a contrast-coded effect (D44), correlated subject parameters
# and their correlations, all from one recovery_grid() call, reported
# split by estimand family. The mock fitter stands in for the sampler;
# whether a real fit's draw names match the mock's is a scratch check
# outside the suite (local/dev/STATE-milestone-12.md).

combined_grid_run <- function(dir, fitter, reps = 2L) {
  suppressMessages(recovery_grid(
    bmm::mixture2p(resp_error = "y"),
    grid = small_grid(),
    pars = c(kappa = log(8), kappa_task2 = log(4), thetat = 0.5),
    dir = dir,
    reps = reps,
    sds = c(kappa = 0.3, thetat = 0.4),
    tasks = c("1", "2"),
    coding = "contrast",
    re_cor = "all",
    levels = c("population", "effect", "sd", "subject"),
    correlations = c("draws", "point"),
    seed = 300,
    .fitter = fitter
  ))
}

test_that("one grid scores intercepts, an effect, SDs, subjects and correlations", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  seen <- new.env()
  mock <- grid_mock_fitter()
  fitter <- function(formula, data, model, prior = NULL, ...) {
    seen$formula <- formula
    mock$fitter(formula, data, model, prior, ...)
  }
  reps <- 2L
  out <- combined_grid_run(dir, fitter, reps = reps)

  # one correlated contrast model per cell: the intercept, the effect and
  # both their subject deviations share one correlation block
  expect_equal(
    deparse1(seen$formula$kappa), "kappa ~ 1 + task + (1 + task | p | id)"
  )
  expect_true(all(attr(out, "cells")$status == "ok"))

  n_cells <- nrow(small_grid())
  subjects <- sum(small_grid()$n_subjects)
  all_terms <- c("kappa", "kappa_task1", "thetat", "thetat_task1")
  expected <- list(
    population = list(terms = c("kappa", "thetat"), per_rep = 2L * n_cells),
    effect = list(
      terms = c("kappa_task1", "thetat_task1"), per_rep = 2L * n_cells
    ),
    sd = list(terms = all_terms, per_rep = 4L * n_cells),
    subject = list(terms = all_terms, per_rep = 4L * subjects)
  )

  families <- split(tibble::as_tibble(out), out$level)
  expect_setequal(names(families), names(expected))
  # every row lands in exactly one family
  expect_identical(sum(vapply(families, nrow, 0L)), nrow(out))
  for (family in names(expected)) {
    rows <- families[[family]]
    expect_setequal(rows$term, expected[[family]]$terms)
    expect_identical(
      nrow(rows), expected[[family]]$per_rep * reps,
      info = family
    )
    key <- rows[c("condition", "replication", "term", "id")]
    expect_identical(nrow(unique(key)), nrow(rows), info = family)
  }

  # effects and SDs on the link scale, intercepts and subjects on the
  # grid's scale
  expect_equal(unique(families$effect$scale), "link")
  expect_equal(unique(families$sd$scale), "link")
  expect_equal(unique(families$population$scale), "natural")
  expect_equal(unique(families$subject$scale), "natural")

  # the truths of the contrast (D44): the effect is the log ratio of the
  # cells, its SD that of a difference of two independent cells
  effect <- families$effect
  expect_equal(
    unique(effect$true_value[effect$term == "kappa_task1"]), log(4) - log(8)
  )
  sd_rows <- families$sd
  expect_equal(
    unique(sd_rows$true_value[sd_rows$term == "kappa_task1"]), 0.3 * sqrt(2)
  )

  # the fifth family: every pair of the four correlated terms, for each
  # estimator, cell and replication
  cors <- attr(out, "correlations")
  expect_s3_class(cors, "bmmtools_cor_recovery")
  expect_identical(
    nrow(cors), as.integer(choose(length(all_terms), 2)) * 2L * n_cells * reps
  )
  expect_identical(
    nrow(unique(cors[c("condition", "replication", "term", "estimator")])),
    nrow(cors)
  )

  # no row lost between the cell files and the scored result
  sidecars <- lapply(
    list.files(dir, pattern = "-est\\.rds$", full.names = TRUE), readRDS
  )
  expect_length(sidecars, n_cells * reps)
  expect_identical(
    sum(vapply(sidecars, function(s) nrow(s$estimates), 0L)), nrow(out)
  )
  expect_identical(
    sum(vapply(sidecars, function(s) nrow(s$cor_estimates), 0L)), nrow(cors)
  )
})

test_that("each family summarises to its own rows of the combined summary", {
  skip_if_not_installed("bmm")
  out <- combined_grid_run(withr::local_tempdir(), grid_mock_fitter()$fitter)
  whole <- summary(out)
  expect_identical(
    nrow(whole),
    nrow(unique(tibble::as_tibble(out)[c("condition", "level", "term")]))
  )
  for (family in unique(out$level)) {
    part <- summary(out[out$level == family, ])
    expect_equal(
      tibble::as_tibble(part),
      tibble::as_tibble(whole[whole$level == family, ]),
      info = family
    )
  }
  cor_whole <- summary(attr(out, "correlations"))
  expect_setequal(cor_whole$estimator, c("draws", "point"))
})
