# summary(precision = ): the share of intervals no wider than a target
# (Milestone 24.2, D95). Widths are hand-built, so every share is a count.

#' Population rows of one term whose `ci_level` intervals have given widths
#' @noRd
width_rows <- function(width, term = "kappa", level = "population",
                       replication = seq_along(width), id = NA_character_) {
  fake_estimates(term,
    estimate = 1, ci_low = 1 - width / 2, ci_high = 1 + width / 2,
    level = level, id = id, replication = replication
  )
}

precision_columns <- c("precision_target", "precision_met", "precision_met_mcse")

test_that("without precision the summary has no precision columns", {
  out <- summary(recover(width_rows(c(1, 2)), fake_truth("kappa", 1),
    scale = "link"
  ))
  expect_false(any(precision_columns %in% names(out)))
  expect_named(out, recovery_summary_columns())
})

test_that("a named vector is a population target; the share is inclusive", {
  # widths 0.5, 1, 1, 2, 3 against a target of 1: three of five
  rows <- width_rows(c(0.5, 1, 1, 2, 3))
  out <- summary(recover(rows, fake_truth("kappa", 1), scale = "link"),
    precision = c(kappa = 1)
  )
  expect_identical(names(out)[length(names(out)) - 2:0], precision_columns)
  expect_equal(out$precision_target, 1)
  expect_equal(out$precision_met, 3 / 5)
  expect_equal(out$precision_met_mcse, sqrt(0.6 * 0.4 / 5))
})

test_that("a width equal to the target in floating point still counts", {
  # 1.1 - 0.8 is 0.30000000000000004, one ulp above 0.3
  rows <- fake_estimates("kappa", 1, ci_low = 0.8, ci_high = 1.1)
  out <- summary(recover(rows, fake_truth("kappa", 1), scale = "link"),
    precision = c(kappa = 0.3)
  )
  expect_equal(out$precision_met, 1)
})

test_that("a width a little above the target is not met", {
  # the slack covers rounding of ci_high - ci_low, not a 1e-6 excess
  rows <- fake_estimates("kappa", 1, ci_low = 0.5, ci_high = 1.5 + 1e-6)
  out <- summary(recover(rows, fake_truth("kappa", 1), scale = "link"),
    precision = c(kappa = 1)
  )
  expect_equal(out$precision_met, 0)
})

test_that("a table sets targets per term and level; others are NA", {
  rows <- dplyr::bind_rows(
    width_rows(c(1, 3), "kappa"),
    width_rows(c(0.1, 0.4), "thetat"),
    width_rows(c(0.2, 0.2), "kappa", level = "sd")
  )
  truth <- list(
    population = fake_truth(c("kappa", "thetat"), 1),
    sd = fake_truth("kappa", 1)
  )
  recovery <- recover(rows, truth,
    level = c("population", "sd"),
    scale = "link"
  )
  target <- tibble::tibble(
    term = c("kappa", "kappa"), level = c("population", "sd"),
    width = c(2, 0.1)
  )
  out <- summary(recovery, precision = target)
  key <- paste(out$term, out$level)
  expect_equal(
    out$precision_target[match(c(
      "kappa population", "kappa sd",
      "thetat population"
    ), key)],
    c(2, 0.1, NA)
  )
  expect_equal(
    out$precision_met[match(c(
      "kappa population", "kappa sd",
      "thetat population"
    ), key)],
    c(0.5, 0, NA)
  )
  expect_identical(out$precision_met_mcse[key == "thetat population"], NA_real_)
})

test_that("the target is read on the row's own scale", {
  # natural kappa intervals are exp() of the link ones, so the same target
  # is met on one scale and not on the other
  rows <- fake_estimates("kappa", 1, ci_low = 0.5, ci_high = 1.5)
  truth <- fake_truth("kappa", 1)
  linked <- summary(recover(rows, truth, scale = "link"),
    precision = c(kappa = 1)
  )
  natural <- summary(
    recover(rows, truth, scale = "natural", links = c(kappa = "log")),
    precision = c(kappa = 1)
  )
  expect_identical(c(linked$scale, natural$scale), c("link", "natural"))
  expect_equal(linked$precision_met, 1)
  expect_equal(natural$precision_met, 0)
})

test_that("subject rows share every interval; the MCSE is between replications", {
  # replication 1: 3 of 4 subjects within 1; replication 2: 1 of 4
  width <- c(0.5, 1, 0.8, 2, 3, 0.9, 2, 4)
  rows <- width_rows(width,
    level = "subject", id = rep(as.character(1:4), 2),
    replication = rep(1:2, each = 4)
  )
  truth <- fake_truth("kappa", 1,
    id = rep(as.character(1:4), 2), replication = rep(1:2, each = 4)
  )
  target <- tibble::tibble(term = "kappa", level = "subject", width = 1)
  out <- summary(recover_subjects(rows, truth, scale = "link"),
    precision = target
  )
  expect_equal(out$precision_met, 4 / 8)
  expect_equal(out$precision_met_mcse, stats::sd(c(0.75, 0.25)) / sqrt(2))
  # never the binomial over the eight intervals
  expect_false(isTRUE(all.equal(out$precision_met_mcse, sqrt(0.25 / 8))))
})

test_that("unbalanced replications: the share counts intervals", {
  rows <- width_rows(c(0.5, 0.5, 0.5, 2, 2, 2, 0.5),
    level = "subject", id = as.character(c(1:3, 1:4)),
    replication = c(1, 1, 1, 2, 2, 2, 2)
  )
  truth <- fake_truth("kappa", 1,
    id = as.character(c(1:3, 1:4)), replication = c(1, 1, 1, 2, 2, 2, 2)
  )
  out <- summary(recover_subjects(rows, truth, scale = "link"),
    precision = data.frame(term = "kappa", level = "subject", width = 1)
  )
  expect_equal(out$precision_met, 4 / 7)
  expect_equal(out$precision_met_mcse, stats::sd(c(1, 0.25)) / sqrt(2))
})

test_that("a missing interval is left out of the share", {
  rows <- width_rows(c(0.5, 2, 0.5))
  rows$ci_high[[3]] <- NA
  out <- summary(recover(rows, fake_truth("kappa", 1), scale = "link"),
    precision = c(kappa = 1)
  )
  expect_equal(out$precision_met, 0.5)
  expect_equal(out$precision_met_mcse, sqrt(0.25 / 2))
})

test_that("a table with a condition column sets targets per condition", {
  rows <- dplyr::bind_rows(
    dplyr::mutate(width_rows(c(1, 3)), condition = "row-1"),
    dplyr::mutate(width_rows(c(1, 3)), condition = "row-2")
  )
  recovery <- recover(rows, fake_truth("kappa", 1), scale = "link")
  expect_identical(unique(recovery$condition), c("row-1", "row-2"))
  target <- tibble::tibble(
    condition = c("row-1", "row-2"), term = "kappa",
    level = "population", width = c(1, 3)
  )
  out <- summary(recovery, precision = target)
  expect_equal(out$precision_target, c(1, 3))
  expect_equal(out$precision_met, c(0.5, 1))
})

test_that("without a condition column a target applies to every condition", {
  s <- summary(recovery_mixture2p, precision = c(kappa = 2))
  population <- s$level == "population" & s$term == "kappa"
  expect_true(all(s$precision_target[population] == 2))
  expect_true(all(is.na(s$precision_target[!population])))
  # by hand from the rows: the share of population kappa intervals <= 2
  rows <- recovery_mixture2p[recovery_mixture2p$level == "population" &
    recovery_mixture2p$term == "kappa", ]
  by_hand <- tapply(rows$ci_high - rows$ci_low <= 2, rows$condition, mean)
  expect_equal(
    s$precision_met[population],
    as.vector(by_hand[s$condition[population]])
  )
})

test_that("a key may repeat across conditions, not within one", {
  rows <- dplyr::bind_rows(
    dplyr::mutate(width_rows(c(1, 3)), condition = "row-1"),
    dplyr::mutate(width_rows(c(1, 3)), condition = "row-2")
  )
  recovery <- recover(rows, fake_truth("kappa", 1), scale = "link")
  twice <- tibble::tibble(
    condition = c("row-1", "row-1"), term = "kappa",
    level = "population", width = c(1, 2)
  )
  expect_error(summary(recovery, precision = twice), "more than once")
  one <- tibble::tibble(
    condition = "row-1", term = "kappa", level = "population", width = 1
  )
  out <- summary(recovery, precision = one)
  expect_equal(out$precision_target, c(1, NA))
})

test_that("a target applies to every estimator of its term and level", {
  rows <- dplyr::bind_rows(
    dplyr::mutate(width_rows(c(1, 3)), estimator = "bayes"),
    dplyr::mutate(width_rows(c(0.5, 0.5)), estimator = "ml")
  )
  out <- summary(recover(rows, fake_truth("kappa", 1), scale = "link"),
    precision = c(kappa = 1)
  )
  expect_identical(out$estimator, c("bayes", "ml"))
  expect_equal(out$precision_met, c(0.5, 1))
})

test_that("one replication of subjects has a share and no MCSE", {
  rows <- width_rows(c(0.5, 2, 0.5),
    level = "subject", id = as.character(1:3), replication = 1
  )
  truth <- fake_truth("kappa", 1, id = as.character(1:3), replication = 1)
  out <- summary(recover_subjects(rows, truth, scale = "link"),
    precision = data.frame(term = "kappa", level = "subject", width = 1)
  )
  expect_equal(out$precision_met, 2 / 3)
  expect_identical(out$precision_met_mcse, NA_real_)
})

test_that("a group without intervals echoes its target and has no share", {
  rows <- width_rows(c(1, 2))
  rows$ci_low <- NA_real_
  out <- summary(recover(rows, fake_truth("kappa", 1), scale = "link"),
    precision = c(kappa = 1)
  )
  expect_equal(out$precision_target, 1)
  expect_identical(c(out$precision_met, out$precision_met_mcse), c(NA_real_, NA))
})

test_that("a group that mixes ci_levels has no share", {
  rows <- fake_estimates("kappa", c(1, 1),
    ci_low = 0.5, ci_high = 1.5, replication = 1:2, ci_level = c(0.95, 0.9)
  )
  out <- summary(recover(rows, fake_truth("kappa", 1), scale = "link"),
    precision = c(kappa = 2)
  )
  expect_identical(out$precision_met, NA_real_)
})

test_that("a summary of a summary refuses precision", {
  s <- summary(recover(width_rows(c(1, 2)), fake_truth("kappa", 1),
    scale = "link"
  ))
  expect_error(summary(s, precision = c(kappa = 1)), class = "rlib_error_dots")
})

test_that("an empty summary with precision has the columns", {
  empty <- recover(width_rows(1), fake_truth("kappa", 1), scale = "link")[0, ]
  out <- summary(empty, precision = c(kappa = 1))
  expect_identical(names(out)[length(names(out)) - 2:0], precision_columns)
  for (column in precision_columns) expect_type(out[[column]], "double")
})

test_that("the print shows the target and the share, not the MCSE", {
  out <- summary(recover(width_rows(c(1, 2)), fake_truth("kappa", 1),
    scale = "link"
  ), precision = c(kappa = 1))
  shown <- summary_print_columns(out)
  expect_true(all(c("precision_target", "precision_met") %in% shown))
  expect_false("precision_met_mcse" %in% shown)
})

# refusals ----------------------------------------------------------------

test_that("precision is refused in a shape it cannot use", {
  recovery <- recover(width_rows(c(1, 2)), fake_truth("kappa", 1),
    scale = "link"
  )
  expect_snapshot(error = TRUE, {
    summary(recovery, precision = 1)
    summary(recovery, precision = c(kappa = 1, 2))
    summary(recovery, precision = c(kappa = 0))
    summary(recovery, precision = c(kappa = -1))
    summary(recovery, precision = c(kappa = Inf))
    summary(recovery, precision = c(kappa = NA))
    summary(recovery, precision = "1")
    summary(recovery, precision = c(thetat = 1))
    summary(recovery, precision = data.frame(term = "kappa", width = 1))
    summary(recovery, precision = data.frame(
      term = "kappa", level = "subject", width = 1
    ))
    summary(recovery, precision = data.frame(
      term = c("kappa", "kappa"), level = "population", width = 1:2
    ))
    summary(recovery, precision = data.frame(
      condition = "row-1", term = "kappa", level = "population", width = 1
    ))
    summary(recovery_mixture2p, precision = data.frame(
      condition = "rwo-1", term = "kappa", level = "population", width = 1
    ))
  })
})

test_that("a misspelt argument is not silently ignored", {
  recovery <- recover(width_rows(c(1, 2)), fake_truth("kappa", 1),
    scale = "link"
  )
  expect_error(summary(recovery, precison = c(kappa = 1)), class = "rlib_error_dots")
})
