# A hand-built recovery grid for the plot_design() tests (Milestone 33,
# D103): two subject counts by two trial counts, four replications, two
# terms, every level, model correlations and a cells table. The values are
# random but seeded, so a test computes what a figure must show from the
# rows themselves rather than from the function under test.

design_fixture <- function(reps = 4L) {
  withr::local_seed(33)
  grid <- data.frame(
    n_subjects = c(10L, 20L, 10L, 20L),
    n_trials = c(30L, 30L, 100L, 100L)
  )
  truth <- c(a = 1, b = -0.5)
  one_cell <- function(i, rep) {
    condition <- sprintf("row-%d", i)
    # precision grows with trials, so the profiles have a direction
    noise <- 0.4 / sqrt(grid$n_trials[[i]] / 30)
    group <- function(level, ids, true_value) {
      n_rows <- length(true_value)
      estimate <- true_value + stats::rnorm(n_rows, 0, noise)
      # every interval its own width, so that the widths' bars over the
      # replications differ from a bar over every row
      post_sd <- noise * stats::runif(n_rows, 0.7, 1.3)
      hw <- 1.96 * post_sd
      # outside tibble(), whose `level` column would mask the argument
      contraction <- if (level == "subject") {
        stats::runif(n_rows, 0.3, 0.8)
      } else {
        stats::runif(n_rows, 0.9, 0.99)
      }
      tibble::tibble(
        term = names(true_value), estimate = estimate,
        ci_low = estimate - hw, ci_high = estimate + hw, ci_method = "eti",
        ci_level = 0.95, rhat = 1, ess_bulk = 1000, ess_tail = 1000,
        true_value = unname(true_value), bias = estimate - true_value,
        covered = true_value >= estimate - hw & true_value <= estimate + hw,
        scale = "link", level = level, id = ids, converged = TRUE,
        condition = condition, estimator = "bayes",
        ci_low_50 = estimate - hw / 2.9, ci_high_50 = estimate + hw / 2.9,
        covered_50 = abs(estimate - true_value) <= hw / 2.9,
        post_mean_link = estimate, post_sd_link = post_sd,
        z = (estimate - true_value) / post_sd,
        contraction = contraction,
        replication = rep
      )
    }
    n <- grid$n_subjects[[i]]
    subject_truth <- c(
      stats::setNames(stats::rnorm(n, 1, 0.5), rep("a", n)),
      stats::setNames(stats::rnorm(n, -0.5, 0.5), rep("b", n))
    )
    dplyr::bind_rows(
      group("population", NA_character_, truth),
      group("sd", NA_character_, c(a = 0.5, b = 0.5)),
      group("subject", rep(as.character(seq_len(n)), 2L), subject_truth)
    )
  }
  rows <- dplyr::bind_rows(lapply(seq_len(nrow(grid)), function(i) {
    dplyr::bind_rows(lapply(seq_len(reps), function(rep) one_cell(i, rep)))
  }))
  x <- new_bmmtools_recovery(rows, scale = "link", ci_level = 0.95)

  cells <- expand.grid(i = seq_len(nrow(grid)), replication = seq_len(reps))
  cells <- tibble::tibble(
    condition = sprintf("row-%d", cells$i),
    replication = cells$replication,
    n_subjects = grid$n_subjects[cells$i],
    n_trials = grid$n_trials[cells$i],
    status = "ok",
    fit_seconds = stats::runif(nrow(cells), 1, 2) * grid$n_subjects[cells$i]
  )
  # the smallest design fails twice, once on each check
  cells$converged <- !(cells$condition == "row-1" & cells$replication <= 2L)
  cells$failed <- ifelse(
    cells$converged, "",
    c("ess_bulk", "ess_bulk, divergent")[pmin(cells$replication, 2L)]
  )
  attr(x, "cells") <- cells
  attr(x, "grid") <- grid
  attr(x, "correlations") <- design_cor_fixture(grid, reps)
  x
}

#' Model correlations for the fixture: one pair, generating value 0.4
#' @noRd
design_cor_fixture <- function(grid, reps) {
  rows <- dplyr::bind_rows(lapply(seq_len(nrow(grid)), function(i) {
    noise <- 0.3 / sqrt(grid$n_subjects[[i]] / 10)
    estimate <- 0.4 + stats::rnorm(reps, 0, noise)
    hw <- 0.5
    tibble::tibble(
      term = "a__b", var1 = "a", var2 = "b", estimator = "model",
      estimate = estimate, ci_low = estimate - hw, ci_high = estimate + hw,
      ci_method = "eti", ci_level = 0.95, rhat = 1, ess_bulk = 1000,
      ess_tail = 1000, true_value = 0.4, sample_value = 0.4,
      bias = estimate - 0.4, bias_sample = estimate - 0.4,
      covered = abs(estimate - 0.4) <= hw,
      covered_sample = abs(estimate - 0.4) <= hw,
      excludes_zero = estimate - hw > 0, scale = "link",
      n = grid$n_subjects[[i]], converged = TRUE,
      condition = sprintf("row-%d", i), ci_low_50 = estimate - hw / 2.9,
      ci_high_50 = estimate + hw / 2.9,
      covered_50 = abs(estimate - 0.4) <= hw / 2.9,
      post_mean_link = estimate, post_sd_link = hw / 1.96,
      z = (estimate - 0.4) / (hw / 1.96), replication = seq_len(reps)
    )
  }))
  new_bmmtools_cor_recovery(rows, scale = "link", ci_level = 0.95)
}

#' The 95 % Wilson score interval, written out here so the tests do not
#' check the package against itself
#' @noRd
wilson_by_hand <- function(k, n, z = stats::qnorm(0.975)) {
  p <- k / n
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  spread <- z / (1 + z^2 / n) * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  c(centre - spread, centre + spread)
}

#' The data of the layer drawn with `geom`, with the panel's facet values
#' @noRd
design_layer <- function(p, geom) {
  built <- ggplot2::ggplot_build(p)
  which <- vapply(
    p$layers, function(l) inherits(l$geom, geom), logical(1)
  )
  stopifnot(sum(which) == 1L)
  data <- built$data[[which(which)]]
  layout <- built$layout$layout
  merge(data, layout, by = "PANEL")
}
