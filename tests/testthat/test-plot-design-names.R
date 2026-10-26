# Tests for how plot_design() finds and names the design columns: a
# design column named like a column of the plot data, a duplicate that no
# design column explains, and a summary whose columns have been moved.

skip_if_not_installed("ggplot2")

# the fixture with its grid column `from` renamed `to`
rename_design <- function(x, from, to) {
  grid <- attr(x, "grid")
  names(grid)[names(grid) == from] <- to
  attr(x, "grid") <- grid
  x
}

# the fixture with its four designs twice, as conditions row-5 to row-8
doubled_designs <- function(x) {
  shift <- function(condition) {
    sprintf("row-%d", as.integer(sub("row-", "", condition)) + 4L)
  }
  rows <- tibble::as_tibble(x)
  more <- rows
  more$condition <- shift(rows$condition)
  y <- new_bmmtools_recovery(
    dplyr::bind_rows(rows, more),
    scale = "link", ci_level = 0.95
  )
  cells <- attr(x, "cells")
  more_cells <- cells
  more_cells$condition <- shift(cells$condition)
  attr(y, "cells") <- dplyr::bind_rows(cells, more_cells)
  grid <- attr(x, "grid")
  attr(y, "grid") <- rbind(grid, grid)
  y
}

# design columns named like the plot data -------------------------------

test_that("a design column named like a profile's own column is refused", {
  x <- rename_design(design_fixture(), "n_subjects", "value")
  # bias alone has no ties, so this drew a colour per bias value
  expect_error(
    plot_design(x, metric = "bias", color_by = "value"),
    "design column value"
  )
  expect_error(plot_design(x, color_by = "value"), "design column value")
  expect_error(
    plot_design(x, color_by = "value"),
    "Rename it in the grid"
  )
  # also when it is on neither axis, as it is still copied in
  y <- design_fixture()
  grid <- attr(y, "grid")
  grid$lower <- c(1, 2, 3, 4)
  attr(y, "grid") <- grid
  expect_error(plot_design(y, metric = "bias"), "design column lower")
  # the cell profile copies the design columns as well
  expect_error(
    plot_design(x, metric = "pass_rate", color_by = "value"),
    "design column value"
  )
  # and so does a grid summary
  expect_error(
    plot_design(summary(x), color_by = "value"), "design column value"
  )
})

test_that("a design column named like a tile's own column is refused", {
  x <- rename_design(design_fixture(), "n_subjects", "label")
  expect_error(
    plot_design(x, type = "convergence", color_by = "label"),
    "design column label"
  )
  y <- rename_design(design_fixture(), "n_subjects", "pass")
  expect_error(
    plot_design(y, type = "convergence", color_by = "pass"),
    "design column pass"
  )
})

test_that("a design column that is no plot column's name still draws", {
  x <- design_fixture()
  p <- plot_design(x, metric = "bias")
  expect_s3_class(p, "ggplot")
  expect_identical(levels(p$data$colour_group), c("10", "20"))
  # a tile's own column does not clash with a profile, nor pooled coverage
  y <- rename_design(x, "n_subjects", "label")
  p <- plot_design(y, metric = "bias", color_by = "label")
  expect_identical(levels(p$data$colour_group), c("10", "20"))
  expect_s3_class(plot_design(y, type = "coverage"), "ggplot")
  z <- rename_design(x, "n_subjects", "value")
  expect_s3_class(plot_design(z, type = "coverage"), "ggplot")
  tiles <- plot_design(x, type = "convergence")$data
  expect_identical(nrow(tiles), 4L)
})

# duplicates no design column explains ----------------------------------

test_that("a duplicate no design column explains names no empty list", {
  data <- data.frame(
    n_trials = c(30L, 30L), n_subjects = c(10L, 10L), term = "a",
    metric = "bias", value = c(0.1, 0.2)
  )
  keys <- c("n_trials", "n_subjects", "term", "metric")
  err <- expect_error(
    check_design_unique(data, keys, c("n_subjects", "n_trials"), NULL),
    "No design column tells them apart"
  )
  expect_no_match(conditionMessage(err), "differ in\\s*,")

  # the same designs twice in a grid: tiles and cell profiles
  y <- doubled_designs(design_fixture())
  err <- expect_error(
    plot_design(y, type = "convergence"), "No design column"
  )
  expect_no_match(conditionMessage(err), "differ in\\s*,")
  expect_error(plot_design(y, metric = "pass_rate"), "No design column")
})

test_that("a duplicate a design column explains still names that column", {
  data <- data.frame(
    n_trials = c(30L, 30L), n_subjects = c(10L, 10L), prior_sd = c(1, 2)
  )
  expect_error(
    check_design_unique(
      data, c("n_trials", "n_subjects"),
      c("n_subjects", "n_trials", "prior_sd"), NULL
    ),
    "They differ in .*prior_sd"
  )
})

# summaries whose columns were moved -------------------------------------

test_that("a summary finds its design columns by name, not position", {
  x <- design_fixture()
  s <- summary(x)
  moved <- dplyr::relocate(s, n_trials, .after = term)
  expect_s3_class(moved, "bmmtools_recovery_summary")
  metrics <- c("bias", "rmse", "coverage")
  from_summary <- plot_design(s, metric = metrics)$data
  from_moved <- plot_design(moved, metric = metrics)$data
  key <- function(d) paste(d$metric, d$term, d$n_subjects, d$n_trials)
  from_moved <- from_moved[match(key(from_summary), key(from_moved)), ]
  expect_equal(from_moved$value, from_summary$value)
  expect_equal(from_moved$lower, from_summary$lower)
  expect_equal(from_moved$upper, from_summary$upper)
  # the axes still read the moved column
  p <- plot_design(moved, x_axis = "n_subjects", color_by = "n_trials")
  expect_identical(levels(p$data$colour_group), c("30", "100"))
})

test_that("a summary with precision columns draws them as metrics", {
  x <- design_fixture()
  s <- summary(x, precision = c(a = 1, b = 1))
  p <- plot_design(s, metric = "precision_met")
  expect_s3_class(p, "ggplot")
  population <- s[s$level == "population", ]
  d <- p$data
  at <- match(
    paste(population$term, population$n_subjects, population$n_trials),
    paste(d$term, d$n_subjects, d$n_trials)
  )
  expect_equal(d$value[at], population$precision_met)
})

test_that("a column added to a summary per term is a metric, not a design", {
  x <- design_fixture()
  s <- dplyr::mutate(summary(x), half_bias = bias / 2)
  expect_s3_class(s, "bmmtools_recovery_summary")
  p <- plot_design(s, metric = "half_bias")
  rows <- s[s$level == "population", ]
  key <- function(d) paste(d$term, d$n_subjects, d$n_trials)
  expect_equal(
    p$data$value, rows$half_bias[match(key(p$data), key(rows))]
  )
})

test_that("a metric with one value per condition is drawn as a metric", {
  x <- design_fixture()
  s <- summary(x)
  s <- s[s$level == "population" & s$term == "a", ]
  s$twice_bias <- 2 * s$bias
  p <- plot_design(s, metric = "twice_bias")
  key <- function(d) paste(d$n_subjects, d$n_trials)
  expect_equal(p$data$value, s$twice_bias[match(key(p$data), key(s))])
  # unasked for, the same column is read as a design column, which a
  # profile along the other two still draws
  expect_s3_class(plot_design(s, metric = "bias"), "ggplot")
})

test_that("the tiles refuse a design column named dark", {
  x <- rename_design(design_fixture(), "n_subjects", "dark")
  expect_error(
    plot_design(x, type = "convergence", color_by = "dark"), "dark"
  )
})
