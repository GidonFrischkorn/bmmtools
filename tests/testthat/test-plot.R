# Tests for plot_recovery(), written against
# local/dev/spec-milestone-1-score-layer.md section 7.
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
  expect_error(plot_recovery(x, facet_by = "nonexistent"), "nonexistent")
  expect_error(plot_recovery(x, color_by = "nonexistent"), "nonexistent")
})

# annotation -------------------------------------------------------------

text_layers <- function(p) which(geom_classes(p) == "GeomText")

test_that("the default plot carries no annotation", {
  skip_if_not_installed("ggplot2")
  p <- plot_recovery(plot_example())
  expect_length(text_layers(p), 0L)
})

test_that("annotate = TRUE labels each panel with r and CCC", {
  skip_if_not_installed("ggplot2")
  x <- plot_example()
  p <- plot_recovery(x, annotate = TRUE)
  layer <- text_layers(p)
  expect_length(layer, 1L)

  labels <- ggplot2::layer_data(p, layer)$label
  expect_length(labels, 2L)
  summarised <- summary(x)
  for (term in summarised$term) {
    row <- summarised[summarised$term == term, ]
    ccc <- sub("^(-?)0\\.", "\\1.", sprintf("%.2f", row$ccc))
    expect_true(any(grepl(paste0("CCC = ", ccc), labels, fixed = TRUE)))
  }
})

test_that("annotate works with a colour mapping", {
  skip_if_not_installed("ggplot2")
  p <- plot_recovery(plot_example(), color_by = "term", annotate = TRUE)
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})

test_that("a single panel is annotated when facet_by is NULL", {
  skip_if_not_installed("ggplot2")
  x <- dplyr::filter(plot_example(), term == "kappa")
  p <- plot_recovery(x, facet_by = NULL, annotate = TRUE)
  labels <- ggplot2::layer_data(p, text_layers(p))$label
  expect_length(labels, 1L)

  # two terms in one panel are two summary rows
  expect_error(
    plot_recovery(plot_example(), facet_by = NULL, annotate = TRUE),
    "mixes terms"
  )
})

test_that("a metric that is not estimable is labelled NA", {
  skip_if_not_installed("ggplot2")
  # two replications: too few pairs for r or CCC
  p <- plot_recovery(plot_example(n_replications = 2L), annotate = TRUE)
  labels <- ggplot2::layer_data(p, text_layers(p))$label
  expect_true(all(grepl("r = NA", labels, fixed = TRUE)))
  expect_true(all(grepl("CCC = NA", labels, fixed = TRUE)))
})

test_that("a panel whose facet value is NA gets its label too", {
  skip_if_not_installed("ggplot2")
  # condition is NA for every row of a recovery that is not a grid
  x <- plot_example()
  p <- plot_recovery(dplyr::filter(x, term == "kappa"),
    facet_by = "condition", annotate = TRUE
  )
  expect_length(ggplot2::layer_data(p, text_layers(p))$label, 1L)

  # and when only some rows are NA, every drawn panel is labelled
  x$condition <- rep(c("a", NA), length.out = nrow(x))
  p <- plot_recovery(dplyr::filter(x, term == "kappa"),
    facet_by = "condition", annotate = TRUE
  )
  panels <- unique(ggplot2::layer_data(p, 3L)$PANEL)
  labelled <- unique(ggplot2::layer_data(p, text_layers(p))$PANEL)
  expect_setequal(labelled, panels)
})

test_that("a value that rounds to zero is labelled without a sign", {
  expect_identical(format_metric(-0.001), ".00")
  expect_identical(format_metric(-0.25), "-.25")
  expect_identical(format_metric(1), "1.00")
})

test_that("a panel spanning several summary rows cannot be annotated", {
  skip_if_not_installed("ggplot2")
  expect_error(
    plot_recovery(recovery_mixture2p, annotate = TRUE),
    "filter"
  )
})

test_that("annotate must be TRUE or FALSE", {
  skip_if_not_installed("ggplot2")
  message <- "annotate. must be"
  expect_error(plot_recovery(plot_example(), annotate = "yes"), message)
  expect_error(plot_recovery(plot_example(), annotate = NA), message)
})

test_that("plot_recovery rejects anything that is not a recovery object", {
  skip_if_not_installed("ggplot2")
  expect_error(plot_recovery(tibble::tibble(a = 1)), "bmmtools_recovery")
  expect_error(plot_recovery(1:10), "bmmtools_recovery")
})

# the cross-check method -------------------------------------------------

cross_check_plot_example <- function(level = "subject",
                                     n_subjects = 6L,
                                     intervals = TRUE) {
  ids <- as.character(seq_len(n_subjects))
  fit <- structure(
    list(parameters = c("kappa", "thetat"), ids = ids), class = "mockfit"
  )
  reference <- if (identical(level, "subject")) {
    tibble::tibble(
      term = rep(c("kappa", "thetat"), each = n_subjects),
      estimate = rep(seq_len(n_subjects) / 10, times = 2L),
      id = rep(ids, times = 2L)
    )
  } else {
    tibble::tibble(term = c("kappa", "thetat"), estimate = c(0.5, 1.5))
  }
  if (intervals) {
    reference$ci_low <- reference$estimate - 0.05
    reference$ci_high <- reference$estimate + 0.05
  }
  cross_check(fit, reference, scale = "link", level = level)
}

test_that("plot_recovery on a cross-check returns a ggplot", {
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot_recovery(cross_check_plot_example()), "ggplot")
})

test_that("the reference is on x and the estimate on y", {
  skip_if_not_installed("ggplot2")
  x <- cross_check_plot_example()
  built <- ggplot2::ggplot_build(plot_recovery(x))
  point_layer <- which(geom_classes(plot_recovery(x)) == "GeomPoint")

  expect_setequal(built$data[[point_layer]]$x, x$reference)
  expect_setequal(built$data[[point_layer]]$y, x$estimate)
})

test_that("both interval layers are drawn when the reference has intervals", {
  skip_if_not_installed("ggplot2")
  p <- plot_recovery(cross_check_plot_example())
  expect_equal(sum(geom_classes(p) == "GeomLinerange"), 2L)
  expect_true("GeomAbline" %in% geom_classes(p))
})

test_that("the horizontal layer is left out without reference intervals", {
  skip_if_not_installed("ggplot2")
  p <- plot_recovery(cross_check_plot_example(intervals = FALSE))
  expect_equal(sum(geom_classes(p) == "GeomLinerange"), 1L)
})

test_that("intervals = FALSE removes both interval layers", {
  skip_if_not_installed("ggplot2")
  p <- plot_recovery(cross_check_plot_example(), intervals = FALSE)
  expect_false("GeomLinerange" %in% geom_classes(p))
})

test_that("identity_line = FALSE removes the identity layer", {
  skip_if_not_installed("ggplot2")
  p <- plot_recovery(cross_check_plot_example(), identity_line = FALSE)
  expect_false("GeomAbline" %in% geom_classes(p))
})

test_that("subject level panels by term, population level does not", {
  skip_if_not_installed("ggplot2")
  subject <- plot_recovery(cross_check_plot_example("subject"))
  expect_s3_class(subject$facet, "FacetWrap")

  population <- plot_recovery(cross_check_plot_example("population"))
  expect_s3_class(population$facet, "FacetNull")

  # asking for panels overrides the level-dependent default
  asked <- plot_recovery(
    cross_check_plot_example("population"),
    facet_by = "term"
  )
  expect_s3_class(asked$facet, "FacetWrap")
})

test_that("the axes name the scale the comparison was made on", {
  skip_if_not_installed("ggplot2")
  p <- plot_recovery(cross_check_plot_example())
  expect_match(p$labels$x, "Reference")
  expect_match(p$labels$x, "link scale")
  expect_match(p$labels$y, "link scale")
})

test_that("a column that is not there is an error", {
  skip_if_not_installed("ggplot2")
  x <- cross_check_plot_example()
  expect_error(plot_recovery(x, facet_by = "nope"), "nope")
  expect_error(plot_recovery(x, color_by = "nope"), "nope")
})

test_that("colouring a cross-check by a column works", {
  skip_if_not_installed("ggplot2")
  p <- plot_recovery(cross_check_plot_example(), color_by = "source")
  expect_s3_class(p, "ggplot")
})

test_that("a cross-check with a broken contract is refused by the plot", {
  skip_if_not_installed("ggplot2")
  x <- cross_check_plot_example()
  broken <- structure(
    tibble::as_tibble(x)[setdiff(names(x), "reference")],
    class = class(x)
  )
  expect_error(plot_recovery(broken), "reference")
})
