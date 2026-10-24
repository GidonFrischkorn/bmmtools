# Grid design columns, the ci_level column and smallest_design()
# (Milestone 24.3, D97 and D98 with the 2026-10-24 amendments). Grid
# summaries here are hand-made tibbles; every expected design is read off
# the rule by hand. The Wilson bound of a perfect rate is n / (n + z^2).

z95 <- stats::qnorm(0.975)

#' The certified design for one target
#' @noRd
smallest_at <- function(x, ...) smallest_design(x, c(...))$smallest

#' The Wilson score lower bound, written out from its definition
#' @noRd
wilson_lower <- function(p, n, z = z95) {
  centre <- p + z^2 / (2 * n)
  spread <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  (centre - spread) / (1 + z^2 / n)
}

#' A hand-made grid summary: one row per design for one term and level
#' @noRd
design_summary <- function(n_trials, ..., n_subjects = 50, term = "kappa",
                           level = "subject", n = n_subjects,
                           n_replications = 10) {
  tibble::tibble(
    condition = paste0("row-", seq_along(n_trials)),
    n_subjects = n_subjects, n_trials = n_trials,
    term = term, estimator = "bayes", level = level, scale = "natural",
    n = n, n_replications = n_replications, ...
  )
}

# the summary gains ci_level and the grid's columns -------------------------

test_that("the summary carries ci_level after sign_recovery_mcse", {
  columns <- recovery_summary_columns()
  at <- match("sign_recovery_mcse", columns)
  expect_identical(columns[at + 1L], "ci_level")
  expect_type(empty_recovery_summary()$ci_level, "double")
  s <- summary(recovery_mixture2p)
  expect_true(all(s$ci_level == attr(recovery_mixture2p, "ci_level")))
})

test_that("a group that mixes ci_levels has ci_level NA", {
  estimates <- fake_estimates("kappa", rep(1, 4),
    replication = 1:4, ci_level = c(0.95, 0.95, 0.9, 0.9)
  )
  out <- summary(recover(estimates, fake_truth("kappa", 1), scale = "link"))
  expect_identical(out$ci_level, NA_real_)
})

test_that("a grid summary carries the grid's columns after condition", {
  s <- summary(recovery_mixture2p)
  grid <- attr(recovery_mixture2p, "grid")
  expect_identical(names(s)[1:3], c("condition", "n_subjects", "n_trials"))
  row <- as.integer(sub("row-", "", s$condition))
  expect_equal(s$n_subjects, grid$n_subjects[row])
  expect_equal(s$n_trials, grid$n_trials[row])
  expect_identical(
    names(s),
    c("condition", "n_subjects", "n_trials", recovery_summary_columns())
  )
})

test_that("the design columns survive filter() before and after summary()", {
  cut <- dplyr::filter(recovery_mixture2p, condition == "row-4")
  s <- summary(cut)
  expect_identical(unique(s$n_subjects), 50)
  expect_identical(unique(s$n_trials), 100)
  both <- dplyr::bind_rows(
    dplyr::filter(summary(recovery_mixture2p), condition == "row-1"),
    s
  )
  expect_identical(unique(both$n_trials), c(30, 100))
})

test_that("a grid column named like a summary column joins as grid_<name>", {
  x <- recovery_mixture2p
  grid <- attr(x, "grid")
  grid$bias <- c(0.1, 0.2, 0.3, 0.4)
  grid$kappa <- 2
  attr(x, "grid") <- grid
  s <- summary(x)
  expect_identical(
    names(s)[1:5],
    c("condition", "n_subjects", "n_trials", "grid_bias", "kappa")
  )
  row <- as.integer(sub("row-", "", s$condition))
  expect_equal(s$grid_bias, grid$bias[row])
  # the metric keeps its name and its meaning
  expect_false(isTRUE(all.equal(s$bias, s$grid_bias)))
})

test_that("a recovery with conditions but no grid has no design columns", {
  rows <- dplyr::mutate(fake_estimates("kappa", c(1, 2)),
    condition = c("a", "b"), replication = 1
  )
  s <- summary(recover(rows, fake_truth("kappa", 1), scale = "link"))
  expect_identical(names(s), c("condition", recovery_summary_columns()))
})

# smallest_design(): the rule ------------------------------------------------

test_that("a bound column is read as it is; the smallest design is found", {
  x <- design_summary(c(30, 100, 300), r_low = c(0.6, 0.82, 0.9))
  out <- smallest_design(x, c(r_low = 0.8))
  expect_s3_class(out, "bmmtools_smallest_design")
  expect_identical(out$smallest, 100)
  expect_identical(out$smallest_point, 100)
  expect_identical(out$status, "met")
})

test_that("noise that passes at 30 and fails at 100 does not give 30", {
  x <- design_summary(c(30, 100, 300), r_low = c(0.85, 0.7, 0.9))
  out <- smallest_design(x, c(r_low = 0.8))
  expect_identical(out$smallest, 300)
})

test_that("met at the grid's smallest design says the minimum may be lower", {
  x <- design_summary(c(30, 100), r_low = c(0.85, 0.9))
  out <- smallest_design(x, c(r_low = 0.8))
  expect_identical(out$smallest, 30)
  expect_identical(out$status, "met at the grid's smallest design")
})

test_that("never met: not reached in this grid, both columns NA", {
  x <- design_summary(c(30, 100), r_low = c(0.5, 0.7))
  out <- smallest_design(x, c(r_low = 0.8))
  expect_identical(c(out$smallest, out$smallest_point), c(NA_real_, NA_real_))
  expect_identical(out$status, "not reached in this grid")
})

test_that("a perfect rate below 16 replications is undecided", {
  # 10/10: the Wilson bound is 10 / (10 + z^2) = 0.722 < 0.8
  x <- design_summary(c(30, 100),
    level = "population", n = 10, detected = c(0.6, 1), detected_mcse = c(
      sqrt(0.6 * 0.4 / 10), 0
    )
  )
  out <- smallest_design(x, c(detected = 0.8))
  expect_identical(out$smallest, NA_real_)
  expect_identical(out$smallest_point, 100)
  expect_identical(out$status, "undecided: too few replications")
  expect_identical(out$replications_needed, 16L)
  expect_match(
    paste(format(out), collapse = "\n"),
    "at least 16 replications"
  )
})

test_that("16 of 16 certifies a rate of 0.80", {
  expect_gte(16 / (16 + z95^2), 0.8)
  expect_lt(15 / (15 + z95^2), 0.8)
  x <- design_summary(100,
    level = "population", n = 16, n_replications = 16,
    detected = 1, detected_mcse = 0
  )
  expect_identical(smallest_design(x, c(detected = 0.8))$smallest, 100)
})

test_that("a rate's bound is Wilson on the effective n", {
  # population: the effective n of a binomial MCSE is the count, 9 here
  p <- 7 / 9
  x <- design_summary(c(30, 100),
    level = "population", n = 9, n_replications = 9,
    sign_recovery = c(0.3, p),
    sign_recovery_mcse = c(0.1, sqrt(p * (1 - p) / 9))
  )
  bound <- wilson_lower(p, 9)
  expect_identical(
    smallest_design(x, c(sign_recovery = bound - 1e-9))$smallest, 100
  )
  expect_identical(
    smallest_design(x, c(sign_recovery = bound + 1e-9))$status,
    "undecided: too few replications"
  )
})

test_that("subject precision_met: effective n, capped at every interval", {
  # 50 subjects x 5 replications; MCSE between replications
  p <- 0.852
  mcse <- 0.0496
  n_eff <- p * (1 - p) / mcse^2
  x <- design_summary(100,
    n_replications = 5, precision_met = p, precision_met_mcse = mcse
  )
  bound <- wilson_lower(p, n_eff)
  expect_identical(smallest_at(x, precision_met = bound - 1e-9), 100)
  expect_identical(smallest_at(x, precision_met = bound + 1e-9), NA_real_)
  # an MCSE so small that n_eff would exceed the 250 intervals: capped
  x$precision_met_mcse <- 1e-4
  capped <- wilson_lower(p, 250)
  expect_identical(smallest_at(x, precision_met = capped - 1e-9), 100)
  expect_identical(smallest_at(x, precision_met = capped + 1e-9), NA_real_)
})

test_that("a rate of 1 or a zero MCSE falls back to the replications", {
  x <- design_summary(100,
    n_replications = 5, precision_met = 0.6, precision_met_mcse = 0
  )
  bound <- wilson_lower(0.6, 5)
  expect_identical(smallest_at(x, precision_met = bound - 1e-9), 100)
  expect_identical(smallest_at(x, precision_met = bound + 1e-9), NA_real_)
})

test_that("a mean's bound is the estimate minus z MCSE, z from level", {
  x <- design_summary(c(30, 100),
    reliability_model = c(0.6, 0.85), reliability_model_mcse = c(0.05, 0.02)
  )
  at95 <- 0.85 - z95 * 0.02
  at90 <- 0.85 - stats::qnorm(0.95) * 0.02
  expect_identical(
    smallest_design(x, c(reliability_model = at95 - 1e-9))$smallest, 100
  )
  expect_identical(
    smallest_design(x, c(reliability_model = at95 + 1e-9))$smallest, NA_real_
  )
  expect_identical(
    smallest_design(x, c(reliability_model = at95 + 1e-9),
      level = 0.9
    )$smallest,
    100
  )
  expect_lt(at95 + 1e-9, at90)
})

test_that("one row per by group, term, level and target", {
  x <- dplyr::bind_rows(
    design_summary(c(30, 100), n_subjects = 20, r_low = c(0.5, 0.7)),
    design_summary(c(30, 100), n_subjects = 50, r_low = c(0.7, 0.82)),
    design_summary(c(30, 100),
      n_subjects = 50, term = "thetat",
      r_low = c(0.81, 0.9)
    )
  )
  out <- smallest_design(x, c(r_low = 0.8))
  expect_identical(nrow(out), 3L)
  expect_identical(out$n_subjects, c(20, 50, 50))
  expect_identical(out$smallest, c(NA, 100, 30))
  expect_true(all(c("term", "level", "target", "threshold") %in% names(out)))
})

test_that("the conjunction of targets or terms is the largest single answer", {
  x <- design_summary(c(30, 100, 300),
    r_low = c(0.85, 0.7, 0.9),
    reliability_model = c(0.9, 0.95, 0.99),
    reliability_model_mcse = 0.01
  )
  out <- smallest_design(x, c(r_low = 0.8, reliability_model = 0.8))
  expect_identical(out$smallest, c(300, 30))
  # both hold at a design and at every larger one only from 300 on
  both <- x$r_low >= 0.8 & x$reliability_model - z95 * 0.01 >= 0.8
  from <- min(which(rev(cumprod(rev(both))) == 1))
  expect_identical(x$n_trials[[from]], max(out$smallest))
})

test_that("a group with the target column NA throughout is left out", {
  x <- dplyr::bind_rows(
    design_summary(c(30, 100),
      reliability_model = c(0.85, 0.9),
      reliability_model_mcse = 0.01
    ),
    design_summary(c(30, 100),
      level = "population", n = 10,
      reliability_model = NA_real_, reliability_model_mcse = NA_real_
    )
  )
  out <- smallest_design(x, c(reliability_model = 0.8))
  expect_identical(out$level, "subject")
})

# smallest_design(): D83's coverage leg ---------------------------------------

test_that("coverage_check requires population coverage within one MCSE", {
  # subject r_low passes from 30; population coverage is off at 30
  coverage_mcse <- sqrt(0.95 * 0.05 / 10)
  x <- dplyr::bind_rows(
    design_summary(c(30, 100), r_low = c(0.85, 0.9)),
    design_summary(c(30, 100),
      level = "population", n = 10,
      coverage = c(0.95 - 2 * coverage_mcse, 0.95 - 0.5 * coverage_mcse),
      coverage_mcse = coverage_mcse, ci_level = 0.95
    )
  )
  x$ci_level <- 0.95
  x$coverage_mcse <- coverage_mcse
  plain <- smallest_design(dplyr::filter(x, level == "subject"), c(r_low = 0.8))
  expect_identical(plain$smallest, 30)
  checked <- smallest_design(x, c(r_low = 0.8), coverage_check = TRUE)
  checked <- checked[checked$level == "subject", ]
  expect_identical(checked$smallest, 100)
})

test_that("coverage_check without population rows is refused", {
  x <- design_summary(c(30, 100),
    r_low = c(0.85, 0.9), coverage = 0.9,
    coverage_mcse = 0.05, ci_level = 0.95
  )
  expect_error(
    smallest_design(x, c(r_low = 0.8), coverage_check = TRUE),
    "population"
  )
})

# smallest_design() on the shipped example (spec 24, section 1d) -------------

test_that("the shipped example gives the spec's worked read-off", {
  # 50 subjects: thetat reaches r_low >= .80 at 100 trials (0.820),
  # kappa does not within the grid (0.690 at 100)
  out <- smallest_design(recovery_mixture2p, c(r_low = 0.8))
  at50 <- out[out$n_subjects == 50 & out$level == "subject", ]
  expect_identical(at50$smallest[at50$term == "thetat"], 100)
  expect_identical(
    at50$status[at50$term == "kappa"], "not reached in this grid"
  )
})

# refusals ------------------------------------------------------------------

test_that("targets without a Monte Carlo bound are refused, naming one", {
  x <- design_summary(c(30, 100), r_low = c(0.5, 0.9))
  expect_snapshot(error = TRUE, {
    smallest_design(x, c(r = 0.8))
    smallest_design(x, c(ccc = 0.8))
    smallest_design(x, c(ci_width = 1))
    smallest_design(x, c(coverage = 0.95))
    smallest_design(x, c(bias = 0))
    smallest_design(x, c(precision_met = 0.8))
    smallest_design(x, 0.8)
    smallest_design(x, c(r_low = NA))
    smallest_design(x, c(r_low = 0.8), along = "n_items")
    smallest_design(x, c(r_low = 0.8), level = 1.2)
    smallest_design(list(), c(r_low = 0.8))
    smallest_design(x, c(r_low = 0.8), along = c("n_trials", "n_subjects"))
    smallest_design(x, c(r_low = 0.8), by = 1)
    smallest_design(dplyr::mutate(x, term = NA), c(r_low = 0.8))
    smallest_design(dplyr::mutate(x, r_low = NA_real_), c(r_low = 0.8))
  })
})

# review of 24.3 ---------------------------------------------------------

test_that("the Wilson bound matches the spec's table (section 1c)", {
  # independent numbers, read off the spec rather than this file's helper
  expect_equal(round(wilson_lower_bound(1, 10, 0.95), 3), 0.722)
  expect_equal(round(wilson_lower_bound(0.9, 10, 0.95), 3), 0.596)
  expect_equal(round(wilson_lower_bound(1, 20, 0.95), 3), 0.839)
  expect_equal(round(wilson_lower_bound(0.9, 20, 0.95), 3), 0.699)
})

test_that("two rows at one design are refused, naming the fix", {
  # two numbers of subjects pooled by `by = NULL`
  x <- dplyr::bind_rows(
    design_summary(c(30, 100), n_subjects = 20, r_low = c(0.9, 0.9)),
    design_summary(c(30, 100), n_subjects = 50, r_low = c(0.5, 0.9))
  )
  expect_error(smallest_design(x, c(r_low = 0.8), by = NULL), "add the grid")
  # a third grid column left out of `by`
  x$bias <- rep(0:1, each = 2)
  x$n_subjects <- 50
  expect_error(smallest_design(x, c(r_low = 0.8)), "more than one row")
  out <- smallest_design(x, c(r_low = 0.8), by = c("n_subjects", "bias"))
  expect_identical(out$bias, 0:1)
  expect_identical(out$smallest, c(30, 100))
})

test_that("a missing or non-numeric along is refused", {
  x <- design_summary(c(30, NA), r_low = c(0.9, 0.9))
  expect_error(smallest_design(x, c(r_low = 0.8)), "no missing values")
  x$n_trials <- factor(c(30, 100))
  expect_error(smallest_design(x, c(r_low = 0.8)), "numeric")
})

test_that("by = NULL reads one group; along keeps its type", {
  x <- design_summary(c(30L, 100L), r_low = c(0.7, 0.9))
  out <- smallest_design(x, c(r_low = 0.8), by = NULL)
  expect_false("n_subjects" %in% names(out))
  expect_identical(out$smallest, 100L)
})

test_that("a component's trials can be the design column", {
  x <- design_summary(c(30, 100), r_low = c(0.7, 0.9))
  names(x)[names(x) == "n_trials"] <- "n_trials_a"
  out <- smallest_design(x, c(r_low = 0.8), along = "n_trials_a")
  expect_identical(out$smallest, 100)
  expect_match(format(out)[[2L]], "n_trials_a = 100")
})

test_that("two estimators get one answer each, told apart in the print", {
  x <- dplyr::bind_rows(
    design_summary(c(30, 100), r_low = c(0.85, 0.9)),
    dplyr::mutate(design_summary(c(30, 100), r_low = c(0.7, 0.9)),
      estimator = "ml"
    )
  )
  out <- smallest_design(x, c(r_low = 0.8))
  expect_identical(out$estimator, c("bayes", "ml"))
  expect_identical(out$smallest, c(30, 100))
  expect_match(format(out)[[3L]], "kappa \\[ml\\] \\(subject\\)")
})

test_that("coverage_check needs the coverage columns", {
  x <- dplyr::bind_rows(
    design_summary(c(30, 100), r_low = c(0.85, 0.9)),
    design_summary(c(30, 100), level = "population", n = 10, coverage = 0.5)
  )
  expect_error(
    smallest_design(x, c(r_low = 0.8), coverage_check = TRUE),
    "coverage_mcse"
  )
})

test_that("the replications note names a perfect rate only when it decides", {
  # p = 1 at 30 trials, p = 0.85 with n_eff 12.75 at 100: the deciding
  # design is not a perfect rate
  p <- 0.85
  x <- design_summary(c(30, 100),
    level = "population", n = 10,
    detected = c(1, p), detected_mcse = c(0, 0.1)
  )
  out <- smallest_design(x, c(detected = 0.8))
  expect_identical(out$status, "undecided: too few replications")
  expect_identical(out$replications_needed, NA_integer_)
  expect_no_match(format(out)[[2L]], "perfect rate")
})

test_that("a rate target of 1 or more, or a target named twice, is refused", {
  x <- design_summary(100,
    level = "population", n = 10, detected = 1, detected_mcse = 0
  )
  expect_error(smallest_design(x, c(detected = 1)), "strictly between")
  expect_error(
    smallest_design(x, c(detected = 0.8, detected = 0.7)),
    "more than once"
  )
})

test_that("the print states each status in a sentence", {
  x <- dplyr::bind_rows(
    design_summary(c(30, 100), r_low = c(0.85, 0.9)),
    design_summary(c(30, 100), term = "thetat", r_low = c(0.7, 0.9)),
    design_summary(c(30, 100), term = "c", r_low = c(0.5, 0.6)),
    design_summary(c(30, 100),
      level = "population", n = 10, term = "mu",
      detected = c(0.5, 1), detected_mcse = c(0.15, 0)
    )
  )
  expect_snapshot(smallest_design(x, c(r_low = 0.8, detected = 0.8)))
})

test_that("grid columns named condition or already prefixed stay apart", {
  x <- recovery_mixture2p
  grid <- attr(x, "grid")
  grid$condition <- letters[1:4]
  grid$bias <- 1:4
  grid$grid_bias <- 5:8
  attr(x, "grid") <- grid
  s <- summary(x)
  renamed <- c("grid_condition", "grid_grid_bias", "grid_bias")
  expect_true(all(renamed %in% names(s)))
  row <- as.integer(sub("row-", "", s$condition))
  expect_identical(s$grid_grid_bias, grid$bias[row])
  expect_identical(s$grid_bias, grid$grid_bias[row])
})
