# Tests for plot_recovery(), written against
# dev/spec-milestone-1-score-layer.md section 7.
#
# These assert structure, not pixels: which layers a plot has, how many
# rows reach them and what the axes are labelled. A pixel comparison
# (vdiffr) would break on every ggplot2 release and tell us nothing about
# whether the plot shows the right numbers.

plot_example <- function(n_replications = 6L, scale = "link") {
  terms <- c("kappa", "thetat")
  estimates <- list()
  truth <- list()
  for (rep in seq_len(n_replications)) {
    true_value <- c(1, 2) + rep / 10
    estimates[[rep]] <- fake_estimates(
      terms,
      estimate = true_value + c(0.05, -0.05), replication = rep
    )
    truth[[rep]] <- fake_truth(terms,
      true_value = true_value,
      replication = rep
    )
  }
  links <- c(kappa = "log", thetat = "logit")
  recover(
    dplyr::bind_rows(estimates), dplyr::bind_rows(truth),
    scale = scale, links = if (scale == "natural") links else NULL
  )
}

geom_classes <- function(p) {
  vapply(p$layers, function(layer) class(layer$geom)[[1L]], character(1))
}

test_that("plot_recovery returns a ggplot", {
  skip_if_not_installed("ggplot2")
  p <- plot_recovery(plot_example())
  expect_s3_class(p, "ggplot")
})

test_that("every scored row reaches the point layer", {
  skip_if_not_installed("ggplot2")
  x <- plot_example()
  p <- plot_recovery(x)
  built <- ggplot2::ggplot_build(p)

  point_layer <- which(geom_classes(p) == "GeomPoint")
  expect_length(point_layer, 1L)
  expect_equal(nrow(built$data[[point_layer]]), nrow(x))
})

test_that("truth is on x and the estimate on y", {
  skip_if_not_installed("ggplot2")
  x <- plot_example()
  p <- plot_recovery(x)
  built <- ggplot2::ggplot_build(p)
  point_layer <- which(geom_classes(p) == "GeomPoint")

  expect_setequal(built$data[[point_layer]]$x, x$true_value)
  expect_setequal(built$data[[point_layer]]$y, x$estimate)
})

test_that("intervals = FALSE removes the interval layer", {
  skip_if_not_installed("ggplot2")
  x <- plot_example()
  with_intervals <- plot_recovery(x, intervals = TRUE)
  without <- plot_recovery(x, intervals = FALSE)

  expect_equal(length(without$layers), length(with_intervals$layers) - 1L)
  expect_true("GeomLinerange" %in% geom_classes(with_intervals))
  expect_false("GeomLinerange" %in% geom_classes(without))
})

test_that("identity_line = FALSE removes the identity layer", {
  skip_if_not_installed("ggplot2")
  x <- plot_example()
  with_line <- plot_recovery(x, identity_line = TRUE)
  without <- plot_recovery(x, identity_line = FALSE)

  expect_equal(length(without$layers), length(with_line$layers) - 1L)
  expect_true("GeomAbline" %in% geom_classes(with_line))
  expect_false("GeomAbline" %in% geom_classes(without))
})

test_that("facet_by splits the panels and NULL keeps one", {
  skip_if_not_installed("ggplot2")
  x <- plot_example()

  faceted <- ggplot2::ggplot_build(plot_recovery(x, facet_by = "term"))
  single <- ggplot2::ggplot_build(plot_recovery(x, facet_by = NULL))

  expect_equal(length(unique(faceted$data[[1L]]$PANEL)), 2L)
  expect_equal(length(unique(single$data[[1L]]$PANEL)), 1L)
})

test_that("color_by maps a column to colour", {
  skip_if_not_installed("ggplot2")
  x <- plot_example()
  p <- plot_recovery(x, facet_by = NULL, color_by = "term")
  built <- ggplot2::ggplot_build(p)
  point_layer <- which(geom_classes(p) == "GeomPoint")

  expect_equal(length(unique(built$data[[point_layer]]$colour)), 2L)
})

test_that("the axis labels name the scale the object was scored on", {
  skip_if_not_installed("ggplot2")
  on_link <- plot_recovery(plot_example(scale = "link"))
  on_natural <- plot_recovery(plot_example(scale = "natural"))

  expect_match(on_link$labels$x, "link")
  expect_match(on_link$labels$y, "link")
  expect_match(on_natural$labels$x, "natural")
})

test_that("an unknown facet or colour column is an error naming it", {
  skip_if_not_installed("ggplot2")
  x <- plot_example()
  expect_error(plot_recovery(x, facet_by = "condition"), "condition")
  expect_error(plot_recovery(x, color_by = "condition"), "condition")
})

test_that("plot_recovery rejects anything that is not a recovery object", {
  skip_if_not_installed("ggplot2")
  expect_error(plot_recovery(tibble::tibble(a = 1)), "bmmtools_recovery")
  expect_error(plot_recovery(1:10), "bmmtools_recovery")
})
