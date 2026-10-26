# plot_design(): which cells and correlations a figure reads: designs
# whose fits errored, filters on replications, the cell join, pooled
# estimators and correlations on two scales.
# Every expectation is computed from the rows or the cells themselves.

skip_if_not_installed("ggplot2")

# a copy of `x` with its rows replaced, keeping the grid's attributes
rebind_rows <- function(x, rows) {
  y <- new_bmmtools_recovery(rows, scale = "link", ci_level = 0.95)
  for (name in c("cells", "grid", "correlations")) {
    attr(y, name) <- attr(x, name)
  }
  y
}

# the tile, or the pass-rate point, of one design
tile_of <- function(data, n_subjects, n_trials) {
  data[data$n_subjects == n_subjects & data$n_trials == n_trials, ]
}

pass_point <- function(data, n_subjects, n_trials) {
  tile_of(data[data$metric == "pass_rate", ], n_subjects, n_trials)
}

# designs whose every fit errored ----------------------------------------

test_that("a design whose every fit errored keeps its tile and pass rate", {
  x <- design_fixture()
  # what recovery_grid() writes for errored fits: no rows, an error cell
  x <- dplyr::filter(x, condition != "row-4")
  cells <- attr(x, "cells")
  at <- cells$condition == "row-4"
  cells$status[at] <- "error"
  cells$converged[at] <- NA
  cells$failed[at] <- NA_character_
  attr(x, "cells") <- cells
  expect_false("row-4" %in% x$condition)

  n_fits <- sum(at)
  tile <- tile_of(plot_design(x, type = "convergence")$data, 20L, 100L)
  expect_identical(nrow(tile), 1L)
  expect_match(tile$label, paste0("0/", n_fits, " pass"), fixed = TRUE)
  expect_match(tile$label, paste0("errors: ", n_fits), fixed = TRUE)
  expect_identical(tile$pass, 0)

  point <- pass_point(plot_design(x, metric = "pass_rate")$data, 20L, 100L)
  expect_identical(nrow(point), 1L)
  expect_identical(point$value, 0)
  expect_equal(c(point$lower, point$upper), wilson_by_hand(0, n_fits))
})

test_that("a grid whose fits of one design all errored draws that design", {
  skip_if_not_installed("bmm")
  mock <- grid_mock_fitter(
    fail_on = function(n, data, dots) nrow(data) == 6L * 40L
  )
  out <- suppressWarnings(suppressMessages(recovery_grid(
    bmm::mixture2p(resp_error = "y"),
    grid = expand.grid(n_subjects = c(4L, 6L), n_trials = c(20L, 40L)),
    pars = c(kappa = log(8), thetat = stats::qlogis(0.75)),
    dir = withr::local_tempdir(), reps = 2L,
    sds = c(kappa = 0.3, thetat = 0.4), seed = 100, chains = 2, iter = 400,
    .fitter = mock$fitter
  )))
  cells <- attr(out, "cells")
  at <- cells$n_subjects == 6L & cells$n_trials == 40L
  # the premise: every fit of 6 x 40 errored and left no rows
  expect_true(all(cells$status[at] != "ok"))
  expect_false(any(unique(cells$condition[at]) %in% out$condition))

  tiles <- plot_design(out, type = "convergence")$data
  expect_identical(nrow(tiles), length(unique(cells$condition)))
  tile <- tile_of(tiles, 6L, 40L)
  expect_match(tile$label, paste0("0/", sum(at), " pass"), fixed = TRUE)
  expect_match(tile$label, paste0("errors: ", sum(at)), fixed = TRUE)

  points <- plot_design(out, metric = "pass_rate")$data
  expect_identical(nrow(points), length(unique(cells$condition)))
  expect_identical(pass_point(points, 6L, 40L)$value, 0)
})

# a filter on replications -----------------------------------------------

test_that("a filter on replications reaches the correlations, not the cells", {
  x <- design_fixture()
  f <- dplyr::filter(x, replication <= 2)
  cells <- attr(x, "cells")

  # an errored fit has no rows, so the cells count every replication
  expect_message(
    tiles <- plot_design(f, type = "convergence")$data, "replication"
  )
  expect_message(
    points <- plot_design(f, metric = "pass_rate")$data, "replication"
  )
  for (i in seq_len(nrow(attr(x, "grid")))) {
    design <- attr(x, "grid")[i, ]
    d <- cells[cells$condition == sprintf("row-%d", i), ]
    k <- sum(d$converged %in% TRUE)
    expect_identical(nrow(d), 4L)
    tile <- tile_of(tiles, design$n_subjects, design$n_trials)
    expect_match(tile$label, paste0("^", k, "/", nrow(d), " pass"))
    point <- pass_point(points, design$n_subjects, design$n_trials)
    expect_equal(point$value, k / nrow(d))
    expect_equal(c(point$lower, point$upper), wilson_by_hand(k, nrow(d)))
  }

  cors <- tibble::as_tibble(attr(x, "correlations"))
  at <- cors$condition == "row-3" & cors$replication <= 2
  p <- plot_design(f, level = "correlation", metric = c("bias", "coverage"))
  bias <- tile_of(p$data[p$data$metric == "bias", ], 10L, 100L)
  expect_equal(bias$value, mean(cors$estimate[at] - 0.4))
  coverage <- tile_of(p$data[p$data$metric == "coverage", ], 10L, 100L)
  expect_equal(
    c(coverage$lower, coverage$upper),
    wilson_by_hand(sum(cors$covered[at]), sum(at))
  )
})

# grid columns the cells also use ----------------------------------------

test_that("a design column named like a cell column is refused by name", {
  for (name in c("status", "converged", "failed", "fit_seconds")) {
    x <- design_fixture()
    grid <- attr(x, "grid")
    grid[[name]] <- 1
    attr(x, "grid") <- grid
    pattern <- paste0("design column named .*", name)
    expect_error(plot_design(x, type = "convergence"), pattern)
    expect_error(plot_design(x, metric = "pass_rate"), pattern)
    # the summary profiles do not read the cells
    expect_s3_class(plot_design(x, metric = "bias"), "ggplot")
  }
})

# pooled coverage of one estimator ---------------------------------------

test_that("pooled coverage refuses rows of two estimators", {
  x <- design_fixture()
  rows <- tibble::as_tibble(x)
  ml <- rows[rows$level == "population", ]
  ml$estimator <- "ml"
  y <- rebind_rows(x, dplyr::bind_rows(rows, ml))
  expect_error(plot_design(y, type = "coverage"), "estimators")
  # one estimator pools as before
  expect_s3_class(plot_design(x, type = "coverage"), "ggplot")
})

# correlations on two scales ---------------------------------------------

test_that("correlations on two scales draw the attribute's scale", {
  x <- design_fixture()
  cors <- tibble::as_tibble(attr(x, "correlations"))
  natural <- cors
  natural$scale <- "natural"
  natural$estimate <- natural$estimate + 10
  attr(x, "correlations") <- new_bmmtools_cor_recovery(
    dplyr::bind_rows(cors, natural), "link", 0.95
  )
  expect_message(
    p <- plot_design(x, level = "correlation", metric = "bias"), "natural"
  )
  at <- cors$condition == "row-3"
  got <- p$data[p$data$n_subjects == 10L & p$data$n_trials == 100L, ]
  expect_equal(got$value, mean(cors$estimate[at] - 0.4))

  # without a scale to choose, the second one is refused by name
  attr(attr(x, "correlations"), "scale") <- NULL
  expect_error(
    plot_design(x, level = "correlation", metric = "bias"),
    "scales"
  )
})

test_that("a replication whose every fit errored keeps its error cells", {
  x <- design_fixture()
  # replication 4 errored in every design, so x holds no row of it
  x <- dplyr::filter(x, replication != 4L)
  cells <- attr(x, "cells")
  at <- cells$replication == 4L
  cells$status[at] <- "error"
  cells$converged[at] <- NA
  cells$failed[at] <- NA_character_
  attr(x, "cells") <- cells
  tiles <- plot_design(x, type = "convergence")$data
  for (i in seq_len(nrow(tiles))) {
    expect_match(tiles$label[[i]], "/4 pass", fixed = TRUE)
    expect_match(tiles$label[[i]], "errors: 1", fixed = TRUE)
  }
  passing <- tapply(
    cells$converged %in% TRUE, cells$condition, sum
  )[c("row-1", "row-2", "row-3", "row-4")]
  expect_equal(sort(tiles$pass), sort(as.vector(passing) / 4))
})

test_that("a design column named like a correlation column is refused", {
  x <- design_fixture()
  grid <- attr(x, "grid")
  names(grid)[names(grid) == "n_trials"] <- "true_value"
  attr(x, "grid") <- grid
  expect_error(
    plot_design(x, level = "correlation", x_axis = "true_value"),
    "true_value"
  )
  # the other levels have no such column
  expect_s3_class(
    plot_design(x, metric = "bias", x_axis = "true_value"), "ggplot"
  )
})

# an errored fit in `x`, as recovery_grid() writes it: an error cell and
# no rows for that condition and replication
with_errored_fit <- function(x, condition, replication) {
  drop <- x$condition == condition & x$replication == replication
  y <- rebind_rows(x, tibble::as_tibble(x)[!drop, ])
  cells <- attr(y, "cells")
  at <- cells$condition == condition & cells$replication == replication
  cells$status[at] <- "error"
  cells$converged[at] <- NA
  cells$failed[at] <- NA_character_
  attr(y, "cells") <- cells
  y
}

test_that("a filter on replications keeps every fit of the design", {
  # the errored fit has no rows, so no filter on `x` can tell its
  # replication apart from one filtered out
  x <- with_errored_fit(design_fixture(), "row-4", 2L)
  two <- dplyr::filter(x, replication <= 2L)
  expect_message(
    tiles <- plot_design(two, type = "convergence")$data, "replication"
  )
  tile <- tile_of(tiles, 20L, 100L)
  expect_identical(nrow(tile), 1L)
  expect_match(tile$label, "3/4 pass", fixed = TRUE)
  expect_match(tile$label, "errors: 1", fixed = TRUE)

  # filtered to the errored replication, the design has no rows left;
  # its fits ran in other replications only, so it is drawn and named
  x <- with_errored_fit(design_fixture(), "row-4", 1L)
  one <- dplyr::filter(x, replication == 1L)
  expect_message(
    tiles <- plot_design(one, type = "convergence")$data, "row-4"
  )
  expect_identical(nrow(tiles), 4L)
  tile <- tile_of(tiles, 20L, 100L)
  expect_match(tile$label, "3/4 pass", fixed = TRUE)
  expect_match(tile$label, "errors: 1", fixed = TRUE)

  x <- with_errored_fit(design_fixture(), "row-2", 4L)
  kept <- dplyr::filter(x, replication <= 2L)
  expect_message(
    tile <- tile_of(plot_design(kept, type = "convergence")$data, 20L, 30L),
    "replication"
  )
  expect_match(tile$label, "3/4 pass", fixed = TRUE)
  expect_match(tile$label, "errors: 1", fixed = TRUE)
})

test_that("a filter on one design keeps that design's errored fits", {
  x <- with_errored_fit(design_fixture(), "row-2", 3L)
  one <- dplyr::filter(x, condition == "row-2")
  expect_no_message(tiles <- plot_design(one, type = "convergence")$data)
  expect_identical(nrow(tiles), 1L)
  expect_match(tiles$label, "3/4 pass", fixed = TRUE)
  expect_match(tiles$label, "errors: 1", fixed = TRUE)
  point <- plot_design(one, metric = "pass_rate")$data
  expect_equal(point$value, 3 / 4)
})

test_that("an errored design a filter cannot reach is named when left out", {
  x <- design_fixture()
  for (r in 1:4) x <- with_errored_fit(x, "row-4", r)
  # row-4 (20 x 100) has no rows, so asking for it selects nothing
  two <- dplyr::filter(x, condition %in% c("row-2", "row-4"))
  expect_message(
    tiles <- plot_design(two, type = "convergence")$data, "row-4"
  )
  expect_identical(nrow(tiles), 1L)
  # unfiltered, it is drawn
  tiles <- plot_design(x, type = "convergence")$data
  expect_match(tile_of(tiles, 20L, 100L)$label, "errors: 4", fixed = TRUE)
})

test_that("a filter on converged rows does not reach the cells", {
  x <- design_fixture()
  rows <- tibble::as_tibble(x)
  cells <- attr(x, "cells")
  # the rows carry their fit's verdict, as recovery_grid() writes them
  verdict <- cells$converged[match(
    paste(rows$condition, rows$replication),
    paste(cells$condition, cells$replication)
  )]
  rows$converged <- verdict
  y <- rebind_rows(x, rows)
  passed <- dplyr::filter(y, converged %in% TRUE)
  tile <- tile_of(plot_design(passed, type = "convergence")$data, 10L, 30L)
  at <- cells$condition == "row-1"
  expect_match(
    tile$label, paste0(sum(cells$converged[at]), "/", sum(at), " pass"),
    fixed = TRUE
  )
  expect_equal(tile$pass, mean(cells$converged[at]))
})

test_that("a filter on conditions drops errored designs it filtered out", {
  x <- design_fixture()
  rows <- tibble::as_tibble(x)
  shifted <- rows
  shifted$condition <- sprintf(
    "row-%d", as.integer(sub("row-", "", rows$condition)) + 4L
  )
  y <- rebind_rows(x, dplyr::bind_rows(rows, shifted))
  grid <- attr(x, "grid")
  attr(y, "grid") <- rbind(
    cbind(grid, prior_sd = 0.5), cbind(grid, prior_sd = 1)
  )
  cells <- attr(x, "cells")
  more <- cells
  more$condition <- sprintf(
    "row-%d", as.integer(sub("row-", "", cells$condition)) + 4L
  )
  attr(y, "cells") <- dplyr::bind_rows(cells, more)
  # every fit of one prior_sd = 1 design errored
  for (r in 1:4) y <- with_errored_fit(y, "row-8", r)
  kept <- dplyr::filter(y, condition %in% sprintf("row-%d", 1:4))
  expect_message(
    tiles <- plot_design(kept, type = "convergence")$data, "row-8"
  )
  expect_identical(nrow(tiles), 4L)
  expect_identical(
    nrow(suppressMessages(plot_design(kept, metric = "pass_rate"))$data), 4L
  )
  # the other half keeps its errored design
  other <- dplyr::filter(y, condition %in% sprintf("row-%d", 5:8))
  tile <- tile_of(plot_design(other, type = "convergence")$data, 20L, 100L)
  expect_match(tile$label, "errors: 4", fixed = TRUE)
})

test_that("a filter on replications is not taken for one on conditions", {
  # row-4 errored in every replication, row-2 in replication 2: filtered
  # to replication 2, row-2 has no rows, yet no design was filtered out
  x <- design_fixture()
  for (r in 1:4) x <- with_errored_fit(x, "row-4", r)
  x <- with_errored_fit(x, "row-2", 2L)
  two <- dplyr::filter(x, replication == 2L)
  expect_message(
    tiles <- plot_design(two, type = "convergence")$data, "row-2"
  )
  expect_identical(nrow(tiles), 4L)
  expect_match(tile_of(tiles, 20L, 30L)$label, "errors: 1", fixed = TRUE)
  expect_match(tile_of(tiles, 20L, 100L)$label, "errors: 4", fixed = TRUE)
})
