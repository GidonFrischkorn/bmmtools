# Tests for the `estimator` column, written against
# local/dev/spec-milestone-8-ml-comparison.md decision 35.
#
# The column exists so that two estimators of the same parameter ---
# a hierarchical posterior and a subject-wise ML fit --- can be scored
# against one truth without being pooled into a single bias and RMSE.
# Everything here runs on hand-built tibbles: no fit, no bmm, no Stan.

# the fill ---------------------------------------------------------------

test_that("an estimates tibble without an estimator is filled 'bayes'", {
  estimates <- fake_estimates(c("a", "b"), estimate = c(2, -1))
  truth <- fake_truth(c("a", "b"), true_value = c(1.5, -3))

  out <- recover(estimates, truth, scale = "link")

  expect_true("estimator" %in% names(out))
  expect_equal(unique(out$estimator), "bayes")
})

test_that("an estimator given on the estimates tibble survives recover", {
  estimates <- fake_estimates(c("a", "b"), estimate = c(2, -1))
  estimates$estimator <- "ml"
  truth <- fake_truth(c("a", "b"), true_value = c(1.5, -3))

  out <- recover(estimates, truth, scale = "link")

  expect_equal(unique(out$estimator), "ml")
})

test_that("the contract carries estimator and recover returns it in order", {
  expect_true("estimator" %in% names(recovery_contract()))
  expect_true("estimator" %in% names(estimates_contract()))
  # optional on input, exactly as converged is
  expect_false("estimator" %in% estimates_required_columns())

  estimates <- fake_estimates("a", estimate = 1)
  out <- recover(estimates, fake_truth("a", 1), scale = "link")
  expect_named(out, recovery_contract_columns())
})

# the reason the column exists -------------------------------------------

test_that("summary does not pool two estimators", {
  # one estimator is biased +1, the other -1; pooled they would cancel
  bayes <- fake_estimates(c("a", "a"), estimate = c(3, 3), replication = 1L)
  bayes$estimator <- "bayes"
  ml <- fake_estimates(c("a", "a"), estimate = c(1, 1), replication = 1L)
  ml$estimator <- "ml"
  truth <- fake_truth("a", true_value = 2)

  out <- recover(dplyr::bind_rows(bayes, ml), truth, scale = "link")
  s <- summary(out)

  expect_equal(nrow(s), 2L)
  expect_setequal(s$estimator, c("bayes", "ml"))
  expect_equal(s$bias[s$estimator == "bayes"], 1)
  expect_equal(s$bias[s$estimator == "ml"], -1)
  # the pooled answer, which is what a missing grouping key would give
  expect_false(any(s$bias == 0))
})

test_that("the summary column list carries estimator", {
  expect_true("estimator" %in% recovery_summary_columns())

  estimates <- fake_estimates("a", estimate = 1)
  s <- summary(recover(estimates, fake_truth("a", 1), scale = "link"))
  expect_named(s, recovery_summary_columns())
})

test_that("summary keeps its row order when a second estimator appears", {
  one <- fake_estimates(c("b", "a"), estimate = c(1, 2))
  one$estimator <- "bayes"
  two <- fake_estimates(c("b", "a"), estimate = c(1, 2))
  two$estimator <- "ml"
  truth <- fake_truth(c("a", "b"), true_value = c(2, 1))

  s <- summary(recover(dplyr::bind_rows(one, two), truth, scale = "link"))

  # terms appear in the order the estimates brought them, not alphabetically
  expect_equal(s$term[s$estimator == "bayes"], c("b", "a"))
})

# the class --------------------------------------------------------------

test_that("dropping estimator demotes the recovery class", {
  estimates <- fake_estimates("a", estimate = 1)
  out <- recover(estimates, fake_truth("a", 1), scale = "link")

  expect_s3_class(out, "bmmtools_recovery")
  expect_false(inherits(dplyr::select(out, -"estimator"), "bmmtools_recovery"))
})

# the shared-helper trap lives in test-correlations.R, where the
# three_reps() helper that builds a correlation object is defined.

# printing ---------------------------------------------------------------

test_that("format names the estimators only when there is more than one", {
  one <- fake_estimates("a", estimate = 1)
  single <- recover(one, fake_truth("a", 1), scale = "link")
  expect_false(any(grepl("Estimator", format(single))))

  bayes <- fake_estimates("a", estimate = 1)
  bayes$estimator <- "bayes"
  ml <- fake_estimates("a", estimate = 1)
  ml$estimator <- "ml"
  both <- recover(
    dplyr::bind_rows(bayes, ml), fake_truth("a", 1),
    scale = "link"
  )
  expect_true(any(grepl("Estimators: bayes, ml", format(both), fixed = TRUE)))
})

test_that("an NA estimator is filled rather than becoming its own group", {
  # bind_rows of a labelled and an unlabelled tibble leaves NA
  unlabelled <- fake_estimates("a", estimate = 1)
  labelled <- fake_estimates("a", estimate = 1)
  labelled$estimator <- "ml"

  out <- recover(
    dplyr::bind_rows(unlabelled, labelled), fake_truth("a", 1),
    scale = "link"
  )

  expect_false(anyNA(out$estimator))
  expect_setequal(out$estimator, c("bayes", "ml"))
})

# annotating a two-estimator panel ---------------------------------------

test_that("a panel split by estimator annotates one line per estimator", {
  skip_if_not_installed("ggplot2")
  terms <- c("a", "b")
  bayes <- fake_estimates(terms, estimate = c(1, 2))
  bayes$estimator <- "bayes"
  ml <- fake_estimates(terms, estimate = c(1.5, 2.5))
  ml$estimator <- "ml"
  truth <- fake_truth(terms, true_value = c(1, 2))
  out <- recover(dplyr::bind_rows(bayes, ml), truth, scale = "link")

  labels <- recovery_panel_labels(out, "term", "estimator")

  expect_equal(nrow(labels), 2L)
  expect_true(all(grepl("bayes:", labels$label)))
  expect_true(all(grepl("ml:", labels$label)))

  expect_s3_class(
    plot_recovery(out, color_by = "estimator", annotate = TRUE),
    "ggplot"
  )
})

test_that("a panel mixing estimators still errors without color_by", {
  skip_if_not_installed("ggplot2")
  bayes <- fake_estimates("a", estimate = 1)
  bayes$estimator <- "bayes"
  ml <- fake_estimates("a", estimate = 1.5)
  ml$estimator <- "ml"
  out <- recover(
    dplyr::bind_rows(bayes, ml), fake_truth("a", 1),
    scale = "link"
  )

  expect_error(
    plot_recovery(out, annotate = TRUE),
    "2 summary rows"
  )
})

# the balance warning ----------------------------------------------------

test_that("scoring warns when two estimators cover different subjects", {
  bayes <- fake_estimates(
    rep("a", 3),
    estimate = c(1, 2, 3), level = "subject", id = c("s1", "s2", "s3")
  )
  bayes$estimator <- "bayes"
  ml <- fake_estimates(
    rep("a", 2),
    estimate = c(1, 2), level = "subject", id = c("s1", "s2")
  )
  ml$estimator <- "ml"
  truth <- fake_truth(
    rep("a", 3),
    true_value = c(1, 2, 3), id = c("s1", "s2", "s3")
  )

  expect_warning(
    recover_subjects(dplyr::bind_rows(bayes, ml), truth, scale = "link"),
    "s3"
  )
})

test_that("a failed ML row warns, although decision 40 keeps the row", {
  # the shape the package actually produces: fit_ml() keeps a failed
  # subject with estimate = NA rather than dropping it, so comparing which
  # ids are present finds nothing. What the metrics drop is the pair.
  bayes <- fake_estimates(
    rep("a", 3),
    estimate = c(1, 2, 3), level = "subject", id = c("s1", "s2", "s3")
  )
  bayes$estimator <- "bayes"
  ml <- fake_estimates(
    rep("a", 3),
    estimate = c(1.1, 2.1, NA_real_), level = "subject",
    id = c("s1", "s2", "s3")
  )
  ml$ci_low[[3L]] <- NA_real_
  ml$ci_high[[3L]] <- NA_real_
  ml$estimator <- "ml"
  truth <- fake_truth(
    rep("a", 3),
    true_value = c(1, 2, 3), id = c("s1", "s2", "s3")
  )

  expect_warning(
    recover_subjects(dplyr::bind_rows(bayes, ml), truth, scale = "link"),
    "s3"
  )
})

test_that("the balance warning is per cell, not pooled over conditions", {
  # subject ids repeat in every cell of a grid, so a whole cell's failed ML
  # fit intersects away against the same ids in another cell unless the
  # comparison is keyed by the cell too
  cell <- function(estimator, condition, ids, estimate) {
    out <- fake_estimates(
      rep("a", length(ids)),
      estimate = estimate, level = "subject", id = ids
    )
    out$estimator <- estimator
    out$condition <- condition
    out
  }
  rows <- dplyr::bind_rows(
    cell("bayes", "row-1", c("s1", "s2"), c(1, 2)),
    cell("bayes", "row-2", c("s1", "s2"), c(1, 2)),
    cell("ml", "row-1", c("s1", "s2"), c(1.1, NA_real_)),
    cell("ml", "row-2", c("s1", "s2"), c(1.1, 2.1))
  )
  truth <- fake_truth(rep("a", 2), true_value = c(1, 2), id = c("s1", "s2"))

  expect_warning(
    recover_subjects(rows, truth, scale = "link"),
    "s2"
  )
})
