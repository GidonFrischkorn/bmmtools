# Tests for plot_design() (Milestone 33.1, D102-D103). What each figure
# must show is worked out from the rows of the fixture, by the formulas
# the documentation states: Morris, White & Crowther (2019, Table 6) for
# the MCSE of a bias and of an MSE, the Wilson score interval for rates.
# Numbers are read from the plot's data; the drawing from ggplot_build().

skip_if_not_installed("ggplot2")

# one row per design and term from the fixture's rows at one level
by_design <- function(x, level) {
  rows <- tibble::as_tibble(x)[x$level == level, ]
  grid <- attr(x, "grid")
  grid$condition <- sprintf("row-%d", seq_len(nrow(grid)))
  merge(rows, grid, by = "condition")
}

pick_row <- function(data, metric, term, n_subjects, n_trials) {
  keep <- data$metric == metric & data$term == term &
    data$n_subjects == n_subjects & data$n_trials == n_trials
  out <- data[keep, ]
  stopifnot(nrow(out) == 1L)
  out
}

# the summary row of one level, term and design
summary_row <- function(x, level, term, n_subjects, n_trials) {
  s <- summary(x)
  keep <- s$level == level & s$term == term & s$n_subjects == n_subjects
  out <- s[keep & s$n_trials == n_trials, ]
  stopifnot(nrow(out) == 1L)
  out
}

# the plot data of one cell metric and design
cell_row <- function(data, metric, n_subjects, n_trials) {
  keep <- data$metric == metric & data$n_subjects == n_subjects
  data[keep & data$n_trials == n_trials, ]
}

# the profile ------------------------------------------------------------

test_that("a profile: n_trials on a log axis, a line per n_subjects", {
  x <- design_fixture()
  p <- plot_design(x)
  expect_s3_class(p, "ggplot")

  built <- ggplot2::ggplot_build(p)
  x_scale <- built$layout$panel_scales_x[[1]]
  expect_identical(x_scale$get_transformation()$name, "log-10")
  expect_equal(sort(x_scale$get_breaks()), log10(c(30, 100)))

  lines <- design_layer(p, "GeomLine")
  # every panel holds one line per subject count, each through both
  # trial counts
  per_panel <- tapply(lines$group, lines$PANEL, function(g) length(unique(g)))
  expect_true(all(per_panel == 2L))
  expect_true(all(table(lines$PANEL, lines$group) %in% c(0L, 2L)))
  # Okabe-Ito blue and orange, in the order of the subject counts
  expect_setequal(unique(lines$colour), c("#0072B2", "#D55E00"))
  expect_identical(
    unique(lines$colour[order(lines$group)]), c("#0072B2", "#D55E00")
  )
})

test_that("panels are metrics by terms, terms in the order of the object", {
  x <- design_fixture()
  p <- plot_design(x, metric = c("coverage", "bias"))
  layout <- ggplot2::ggplot_build(p)$layout$layout
  expect_identical(levels(layout$metric), c("coverage", "bias"))
  expect_identical(levels(layout$term), c("a", "b"))
  expect_identical(nrow(layout), 4L)

  # the population default: bias, RMSE, coverage and CrI width
  default <- plot_design(x)
  expect_identical(
    levels(default$data$metric), c("bias", "rmse", "coverage", "ci_width")
  )
  expect_identical(
    levels(plot_design(x, level = "subject")$data$metric),
    c("r", "reliability_model", "coverage")
  )
})

test_that("bias is the mean error with its Table 6 MCSE", {
  x <- design_fixture()
  p <- plot_design(x, metric = "bias")
  rows <- by_design(x, "population")
  for (term in c("a", "b")) {
    for (n_trials in c(30L, 100L)) {
      e <- rows$estimate - rows$true_value
      at <- rows$term == term & rows$n_trials == n_trials &
        rows$n_subjects == 20L
      se <- sqrt(sum((e[at] - mean(e[at]))^2) / (sum(at) * (sum(at) - 1)))
      got <- pick_row(p$data, "bias", term, 20L, n_trials)
      expect_equal(got$value, mean(e[at]))
      expect_equal(c(got$lower, got$upper), mean(e[at]) + c(-1, 1) * se)
    }
  }
})

test_that("RMSE carries the MCSE of the MSE through the delta method", {
  x <- design_fixture()
  p <- plot_design(x, level = "sd", metric = "rmse")
  rows <- by_design(x, "sd")
  at <- rows$term == "b" & rows$n_subjects == 10L & rows$n_trials == 100L
  e2 <- (rows$estimate[at] - rows$true_value[at])^2
  mse <- mean(e2)
  n <- sum(at)
  mse_se <- sqrt(sum((e2 - mse)^2) / (n * (n - 1)))
  got <- pick_row(p$data, "rmse", "b", 10L, 100L)
  expect_equal(got$value, sqrt(mse))
  expect_equal(got$upper - got$value, mse_se / (2 * sqrt(mse)))
  expect_equal(got$value - got$lower, mse_se / (2 * sqrt(mse)))
})

test_that("CrI width and contraction carry the SD over replications", {
  x <- design_fixture()
  p <- plot_design(x, metric = c("ci_width", "contraction"))
  rows <- by_design(x, "population")
  at <- rows$term == "a" & rows$n_subjects == 10L & rows$n_trials == 30L
  width <- rows$ci_high[at] - rows$ci_low[at]
  got <- pick_row(p$data, "ci_width", "a", 10L, 30L)
  expect_equal(got$value, mean(width))
  # the fixture's widths are equal within a design, so the bar is zero
  expect_equal(got$upper - got$lower, 2 * stats::sd(width) / sqrt(sum(at)))
  got <- pick_row(p$data, "contraction", "a", 10L, 30L)
  ctr <- rows$contraction[at]
  expect_equal(got$value, mean(ctr))
  expect_equal(got$upper - got$value, stats::sd(ctr) / sqrt(sum(at)))
})

test_that("design coverage has a Wilson interval and the nominal line", {
  x <- design_fixture()
  p <- plot_design(x, metric = c("bias", "coverage"))
  rows <- by_design(x, "population")
  for (term in c("a", "b")) {
    at <- rows$term == term & rows$n_subjects == 20L & rows$n_trials == 30L
    got <- pick_row(p$data, "coverage", term, 20L, 30L)
    expect_equal(got$value, mean(rows$covered[at]))
    expect_equal(
      c(got$lower, got$upper), wilson_by_hand(sum(rows$covered[at]), sum(at))
    )
  }

  ref <- design_layer(p, "GeomHline")
  coverage_line <- unique(ref$yintercept[ref$metric == "coverage"])
  bias_line <- unique(ref$yintercept[ref$metric == "bias"])
  expect_identical(coverage_line, 0.95)
  expect_identical(bias_line, 0)
})

test_that("a design where every interval covered keeps a lower bound below 1", {
  x <- design_fixture()
  rows <- tibble::as_tibble(x)
  rows$covered[rows$level == "population"] <- TRUE
  attrs <- attributes(x)[c("cells", "grid", "correlations")]
  y <- new_bmmtools_recovery(rows, scale = "link", ci_level = 0.95)
  for (name in names(attrs)) attr(y, name) <- attrs[[name]]
  got <- pick_row(
    plot_design(y, metric = "coverage")$data,
    "coverage", "a", 10L, 30L
  )
  z <- stats::qnorm(0.975)
  expect_equal(c(got$value, got$upper), c(1, 1))
  expect_equal(got$lower, 4 / (4 + z^2))
})

test_that("subject r is drawn with its pooled interval from summary()", {
  x <- design_fixture()
  p <- plot_design(x, level = "subject", metric = "r")
  s <- summary_row(x, "subject", "b", 20L, 100L)
  got <- pick_row(p$data, "r", "b", 20L, 100L)
  expect_equal(c(got$value, got$lower, got$upper), c(s$r, s$r_low, s$r_high))
})

test_that("subject coverage takes the Wilson interval at the effective n", {
  x <- design_fixture()
  p <- plot_design(x, level = "subject", metric = "coverage")
  s <- summary_row(x, "subject", "a", 10L, 30L)
  # the binomial count behind the between-replication MCSE, capped at
  # every interval there is (D97)
  n <- min(
    s$coverage * (1 - s$coverage) / s$coverage_mcse^2,
    s$n * s$n_replications
  )
  got <- pick_row(p$data, "coverage", "a", 10L, 30L)
  expect_equal(
    c(got$lower, got$upper), wilson_by_hand(s$coverage * n, n)
  )
})

test_that("reliability is drawn with its summary MCSE", {
  x <- design_fixture()
  p <- plot_design(x, level = "subject", metric = "reliability_model")
  s <- summary_row(x, "subject", "a", 20L, 30L)
  got <- pick_row(p$data, "reliability_model", "a", 20L, 30L)
  expect_equal(
    c(got$value, got$lower, got$upper),
    s$reliability_model + c(0, -1, 1) * s$reliability_model_mcse
  )
})

test_that("pass rate and seconds per fit are profiles of the cells", {
  x <- design_fixture()
  cells <- attr(x, "cells")
  p <- plot_design(x, metric = c("pass_rate", "fit_seconds"))
  expect_false("term" %in% names(p$data))

  at <- cells$n_subjects == 10L & cells$n_trials == 30L
  got <- cell_row(p$data, "pass_rate", 10L, 30L)
  expect_equal(got$value, 0.5)
  expect_equal(c(got$lower, got$upper), wilson_by_hand(2, 4))

  got <- cell_row(p$data, "fit_seconds", 10L, 30L)
  secs <- cells$fit_seconds[at]
  expect_equal(got$value, mean(secs))
  expect_equal(got$upper - got$value, stats::sd(secs) / sqrt(sum(at)))
})

test_that("model correlations are a level of their own", {
  x <- design_fixture()
  p <- plot_design(x, level = "correlation", metric = c("bias", "coverage"))
  cors <- tibble::as_tibble(attr(x, "correlations"))
  at <- cors$condition == "row-3"
  got <- pick_row(p$data, "bias", "a__b", 10L, 100L)
  expect_equal(got$value, mean(cors$estimate[at] - 0.4))
  got <- pick_row(p$data, "coverage", "a__b", 10L, 100L)
  expect_equal(
    c(got$lower, got$upper), wilson_by_hand(sum(cors$covered[at]), sum(at))
  )
})

test_that("a grid with another design column can be drawn along it", {
  x <- design_fixture()
  p <- plot_design(x,
    metric = "bias", x_axis = "n_subjects",
    color_by = "n_trials"
  )
  built <- ggplot2::ggplot_build(p)
  expect_equal(
    sort(built$layout$panel_scales_x[[1]]$get_breaks()), log10(c(10, 20))
  )
  expect_identical(levels(p$data$colour_group), c("30", "100"))
})

# summaries -------------------------------------------------------------

test_that("a grid summary draws the same values, bars where it has them", {
  x <- design_fixture()
  from_rows <- plot_design(x, metric = c("bias", "rmse", "coverage"))$data
  from_summary <- plot_design(
    summary(x),
    metric = c("bias", "rmse", "coverage")
  )$data
  key <- function(d) paste(d$metric, d$term, d$n_subjects, d$n_trials)
  from_summary <- from_summary[match(key(from_rows), key(from_summary)), ]
  expect_equal(from_summary$value, from_rows$value)
  keep <- from_rows$metric != "rmse"
  expect_equal(from_summary$lower[keep], from_rows$lower[keep])
  # the MSE's MCSE needs the replications, which a summary has dropped
  expect_true(all(is.na(from_summary$lower[!keep])))
})

test_that("a summary without design columns is refused", {
  x <- design_fixture()
  s <- summary(x)
  s$n_subjects <- NULL
  s$n_trials <- NULL
  expect_error(plot_design(s), "design")
  expect_error(plot_design(summary(x), type = "convergence"), "recovery")
})

# the convergence tile and pooled coverage --------------------------------

test_that("the convergence tile shows the share of fits passing per design", {
  x <- design_fixture()
  p <- plot_design(x, type = "convergence")
  expect_identical(nrow(p$data), 4L)
  smallest <- p$data[p$data$n_subjects == 10L & p$data$n_trials == 30L, ]
  expect_equal(smallest$pass, 0.5)
  expect_match(smallest$label, "2/4", fixed = TRUE)
  expect_match(smallest$label, "bulk ESS", fixed = TRUE)
  expect_match(smallest$label, "divergen", fixed = TRUE)
  others <- p$data[!(p$data$n_subjects == 10L & p$data$n_trials == 30L), ]
  expect_true(all(others$pass == 1))
  expect_true(all(grepl("4/4", others$label, fixed = TRUE)))
  tiles <- design_layer(p, "GeomTile")
  expect_identical(nrow(tiles), 4L)
})

test_that("pooled coverage pools every design per level and term", {
  x <- design_fixture()
  p <- plot_design(x, type = "coverage")
  rows <- tibble::as_tibble(x)
  for (level in c("population", "sd")) {
    at <- rows$level == level & rows$term == "a"
    at_a <- p$data$level == level & p$data$term == "a"
    got <- p$data[at_a & p$data$nominal == 0.95, ]
    expect_equal(got$value, mean(rows$covered[at]))
    expect_equal(
      c(got$lower, got$upper), wilson_by_hand(sum(rows$covered[at]), sum(at))
    )
    got50 <- p$data[at_a & p$data$nominal == 0.5, ]
    expect_equal(got50$value, mean(rows$covered_50[at]))
  }
  expect_false("subject" %in% p$data$level)
  ref <- design_layer(p, "GeomVline")
  expect_setequal(unique(ref$xintercept), c(0.5, 0.95))
})

# the shipped grid --------------------------------------------------------

test_that("the shipped example grid draws at both of its levels", {
  p <- plot_design(recovery_mixture2p)
  expect_identical(levels(p$data$term), c("kappa", "thetat"))
  expect_identical(nrow(p$data), 4L * 2L * 4L)
  p <- plot_design(recovery_mixture2p, level = "subject")
  expect_true(all(!is.na(p$data$value[p$data$metric == "r"])))
  expect_s3_class(
    plot_design(recovery_mixture2p, type = "convergence"),
    "ggplot"
  )
  expect_error(
    plot_design(recovery_mixture2p, level = "sd"), "does not have"
  )
})

# errors ------------------------------------------------------------------

test_that("plot_design() refuses what it cannot draw", {
  x <- design_fixture()
  plain <- x
  attr(plain, "grid") <- NULL
  expect_error(plot_design(plain), "grid")
  expect_error(plot_design(x, x_axis = "seed"), "x_axis")
  expect_error(plot_design(x, color_by = "n_trials"), "color_by")
  expect_error(plot_design(x, metric = "nonsense"), "nonsense")
  expect_error(plot_design(x, metric = c("bias", "pass_rate")), "cells")
  expect_error(plot_design(x, level = "effect"), "effect")
  expect_error(plot_design(x, type = "tiles"), "type")
  no_cells <- x
  attr(no_cells, "cells") <- NULL
  expect_error(plot_design(no_cells, type = "convergence"), "cells")
  no_cors <- x
  attr(no_cors, "correlations") <- NULL
  expect_error(plot_design(no_cors, level = "correlation"), "correlation")
  expect_error(plot_design(1), "bmmtools_recovery")
  expect_error(plot_design(x, foo = 1), "foo")
})

test_that("a grid run through the mock fitter draws end to end", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  out <- suppressMessages(recovery_grid(
    bmm::mixture2p(resp_error = "y"),
    grid = expand.grid(n_subjects = c(4L, 6L), n_trials = c(20L, 40L)),
    pars = c(kappa = log(8), thetat = stats::qlogis(0.75)),
    dir = dir, reps = 2L, sds = c(kappa = 0.3, thetat = 0.4), seed = 100,
    chains = 2, iter = 400, .fitter = grid_mock_fitter()$fitter
  ))
  p <- plot_design(out)
  expect_setequal(unique(p$data$n_trials), c(20L, 40L))
  expect_setequal(unique(p$data$n_subjects), c(4L, 6L))
  expect_no_error(ggplot2::ggplot_build(p))
  expect_no_error(
    ggplot2::ggplot_build(plot_design(out, type = "convergence"))
  )
})

# what the review of 33.1 found ------------------------------------------

# a copy of `x` with its rows replaced, keeping the grid's attributes
with_rows <- function(x, rows) {
  y <- new_bmmtools_recovery(rows, scale = "link", ci_level = 0.95)
  for (name in c("cells", "grid", "correlations")) {
    attr(y, name) <- attr(x, name)
  }
  y
}

test_that("designs that share both axes are refused, not drawn as one line", {
  x <- design_fixture()
  # the four designs crossed with a third design column
  rows <- tibble::as_tibble(x)
  shifted <- rows
  shifted$condition <- sprintf(
    "row-%d", as.integer(sub("row-", "", rows$condition)) + 4L
  )
  y <- with_rows(x, dplyr::bind_rows(rows, shifted))
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
  expect_error(plot_design(y, metric = "bias"), "prior_sd")
  expect_error(plot_design(y, type = "convergence"), "prior_sd")
  expect_error(plot_design(y, metric = "pass_rate"), "prior_sd")
  # one value of the third column draws
  kept <- dplyr::filter(y, condition %in% sprintf("row-%d", 1:4))
  expect_s3_class(plot_design(kept, metric = "bias"), "ggplot")
  expect_identical(
    nrow(plot_design(kept, type = "convergence")$data), 4L
  )
})

test_that("cells without a failed column or seconds still draw their tiles", {
  x <- design_fixture()
  cells <- attr(x, "cells")
  cells$failed <- NULL
  cells$fit_seconds <- NULL
  attr(x, "cells") <- cells
  p <- plot_design(x, type = "convergence")
  smallest <- p$data[p$data$n_subjects == 10L & p$data$n_trials == 30L, ]
  expect_identical(smallest$label, "2/4 pass")
  expect_error(plot_design(x, metric = "fit_seconds"), "fit_seconds")
  expect_equal(
    plot_design(x, metric = "pass_rate")$data$value[[1L]], 0.5
  )
})

test_that("a design whose fits failed shows them as errors, and none passing", {
  x <- design_fixture()
  cells <- attr(x, "cells")
  at <- cells$condition == "row-4"
  cells$status[at] <- "error"
  cells$converged[at] <- NA
  cells$failed[at] <- NA_character_
  # a fit that ran, but the gate assessed no variable
  cells$converged[cells$condition == "row-3" & cells$replication == 1L] <- NA
  attr(x, "cells") <- cells
  p <- plot_design(x, type = "convergence")
  failed <- p$data[p$data$n_subjects == 20L & p$data$n_trials == 100L, ]
  expect_identical(failed$pass, 0)
  expect_match(failed$label, "errors: 4", fixed = TRUE)
  unassessed <- p$data[p$data$n_subjects == 10L & p$data$n_trials == 100L, ]
  expect_match(unassessed$label, "3/4 pass", fixed = TRUE)
  expect_match(unassessed$label, "not assessed: 1", fixed = TRUE)
})

test_that("pooled coverage needs no trial column; the profiles do", {
  x <- design_fixture()
  grid <- attr(x, "grid")
  names(grid)[names(grid) == "n_trials"] <- "n_trials_a"
  attr(x, "grid") <- grid
  expect_s3_class(plot_design(x, type = "coverage"), "ggplot")
  expect_error(plot_design(x), "x_axis")
  expect_s3_class(
    plot_design(x, type = "convergence", x_axis = "n_trials_a"), "ggplot"
  )
})

test_that("a grid column named like a summary column is renamed", {
  x <- design_fixture()
  grid <- attr(x, "grid")
  names(grid)[names(grid) == "n_subjects"] <- "bias"
  attr(x, "grid") <- grid
  expect_error(plot_design(x), "grid_bias")
  p <- plot_design(x, metric = "bias", color_by = "grid_bias")
  expect_identical(levels(p$data$colour_group), c("10", "20"))
})

test_that("a rate's Wilson count is the rows that entered it", {
  x <- design_fixture()
  rows <- tibble::as_tibble(x)
  first <- rows$condition == "row-1" & rows$term == "a"
  pop <- which(first & rows$level == "population")
  # one replication without an interval, one with a truth of zero
  rows$ci_low[pop[[1L]]] <- NA_real_
  rows$ci_high[pop[[1L]]] <- NA_real_
  rows$covered[pop[[1L]]] <- NA
  rows$true_value[pop[[2L]]] <- 0
  rows$bias[pop[[2L]]] <- rows$estimate[pop[[2L]]]
  rows$covered[pop[[2L]]] <- rows$ci_low[pop[[2L]]] <= 0 &
    rows$ci_high[pop[[2L]]] >= 0
  y <- with_rows(x, rows)
  p <- plot_design(y, metric = c("coverage", "detected", "sign_recovery"))
  d <- rows[pop, ]
  got <- pick_row(p$data, "coverage", "a", 10L, 30L)
  expect_equal(
    c(got$lower, got$upper),
    wilson_by_hand(sum(d$covered, na.rm = TRUE), 3)
  )
  got <- pick_row(p$data, "detected", "a", 10L, 30L)
  expect_equal(c(got$lower, got$upper), wilson_by_hand(got$value * 3, 3))
  got <- pick_row(p$data, "sign_recovery", "a", 10L, 30L)
  expect_equal(c(got$lower, got$upper), wilson_by_hand(got$value * 3, 3))
})

test_that("correlations of several estimators draw the model's own", {
  x <- design_fixture()
  cors <- tibble::as_tibble(attr(x, "correlations"))
  draws <- cors
  draws$estimator <- "draws"
  draws$estimate <- draws$estimate + 10
  attr(x, "correlations") <- new_bmmtools_cor_recovery(
    dplyr::bind_rows(cors, draws), "link", 0.95
  )
  expect_message(
    p <- plot_design(x, level = "correlation", metric = "bias"), "draws"
  )
  got <- pick_row(p$data, "bias", "a__b", 10L, 100L)
  at <- cors$condition == "row-3"
  expect_equal(got$value, mean(cors$estimate[at] - 0.4))
})

test_that("metric names each metric once", {
  x <- design_fixture()
  expect_error(plot_design(x, metric = c("bias", "bias")), "once")
  expect_error(plot_design(x, metric = character()), "once")
})

test_that("defaults, labels and the bars a level does not have", {
  x <- design_fixture()
  cor_default <- plot_design(x, level = "correlation")
  expect_identical(
    levels(cor_default$data$metric),
    c("bias", "coverage", "ci_width", "rejection_rate")
  )
  # the errors of one replication's subjects are not independent
  subject_bias <- plot_design(x, level = "subject", metric = "bias")$data
  expect_true(all(is.na(subject_bias$lower) & is.na(subject_bias$upper)))
  # strips read as words, and a column bmmtools does not name keeps its own
  labeller <- plot_design(x)$facet$params$labeller
  strips <- labeller(data.frame(metric = c("bias", "ci_width", "mae")))
  expect_identical(
    unlist(strips, use.names = FALSE), c("Bias", "CrI width", "MAE")
  )
  strips <- labeller(data.frame(metric = "calibration_slope"))
  expect_identical(unlist(strips, use.names = FALSE), "calibration_slope")
})

test_that("a design column that is not a positive number gets a plain axis", {
  x <- design_fixture()
  grid <- attr(x, "grid")
  grid$n_trials <- ifelse(grid$n_trials == 30L, "short", "long")
  attr(x, "grid") <- grid
  p <- plot_design(x, metric = "bias")
  x_scale <- ggplot2::ggplot_build(p)$layout$panel_scales_x[[1]]
  expect_true(x_scale$is_discrete())
})

test_that("the summary method and estimators refuse what they cannot draw", {
  x <- design_fixture()
  s <- summary(x)
  expect_error(plot_design(s, metric = "pass_rate"), "cells")
  expect_error(plot_design(s, level = "effect"), "does not have")

  rows <- tibble::as_tibble(x)
  ml <- rows[rows$level == "population", ]
  ml$estimator <- "ml"
  y <- with_rows(x, dplyr::bind_rows(rows, ml))
  expect_error(plot_design(y), "estimator")

  cors <- tibble::as_tibble(attr(x, "correlations"))
  cors$estimator <- "draws"
  point <- cors
  point$estimator <- "point"
  attr(x, "correlations") <- new_bmmtools_cor_recovery(
    dplyr::bind_rows(cors, point), "link", 0.95
  )
  expect_error(plot_design(x, level = "correlation"), "no \"model\"|no model")

  subjects_only <- dplyr::filter(design_fixture(), level == "subject")
  expect_error(plot_design(subjects_only, type = "coverage"), "SD rows")
})
