# The recovery plot (spec section 7).
#
# ggplot2 is in Suggests, so the entry point checks for it rather than
# importing it: a scoring pipeline that never plots installs without it.

#' A readable axis label for the scale an object was scored on
#' @noRd
scale_label <- function(scale, what) {
  paste0(what, " (", scale, " scale)")
}

#' Check that a column a user asked to map exists
#' @noRd
check_plot_column <- function(x, column, arg, call = rlang::caller_env()) {
  if (is.null(column)) {
    return(invisible(NULL))
  }
  if (!is.character(column) || length(column) != 1L) {
    cli::cli_abort(
      "{.arg {arg}} must be a single column name or {.code NULL}.",
      call = call
    )
  }
  if (!column %in% names(x)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must name a column of {.arg x}, not {.val {column}}.",
        i = "Available: {.val {names(x)}}."
      ),
      call = call
    )
  }
  invisible(NULL)
}

#' Format a metric for a plot label: two decimals, no leading zero
#' @noRd
format_metric <- function(x) {
  if (is.na(x)) {
    return("NA")
  }
  # a negative value that rounds to zero would print as "-.00"
  if (round(x, 2L) == 0) x <- 0
  sub("^(-?)0\\.", "\\1.", sprintf("%.2f", x))
}

#' One label per panel, from summary() of that panel's rows
#'
#' Each panel must correspond to exactly one summary row; otherwise the
#' label would show a number for a group the panel does not isolate.
#'
#' @noRd
recovery_panel_labels <- function(x, facet_by, call = rlang::caller_env()) {
  # facet_wrap() draws NA as a panel of its own, so NA is kept as a group
  groups <- if (is.null(facet_by)) {
    list(seq_len(nrow(x)))
  } else {
    # split(drop = TRUE) would also drop a level that is entirely NA
    by_value <- split(
      seq_len(nrow(x)),
      addNA(factor(x[[facet_by]]), ifany = TRUE)
    )
    by_value[lengths(by_value) > 0L]
  }
  labels <- lapply(groups, function(rows) {
    piece <- x[rows, ]
    summarised <- summary(piece)
    if (nrow(summarised) != 1L) {
      problem <- if (is.null(facet_by)) {
        "The plot mixes terms, levels or conditions."
      } else {
        cli::format_inline(
          "Panel {.val {piece[[facet_by]][[1L]]}} mixes levels or conditions."
        )
      }
      cli::cli_abort(
        c(
          "Cannot annotate a panel that holds {nrow(summarised)} summary rows.",
          x = problem,
          i = "Use {.fn dplyr::filter} on {.field level} or \\
               {.field condition} first, or facet by another column."
        ),
        call = call
      )
    }
    out <- tibble::tibble(
      label = paste0(
        # the leading spaces inset both lines by the same amount; a
        # fractional hjust would shift each line by its own width
        " r = ", format_metric(summarised$r), "\n",
        " CCC = ", format_metric(summarised$ccc)
      )
    )
    if (!is.null(facet_by)) out[[facet_by]] <- piece[[facet_by]][[1L]]
    out
  })
  dplyr::bind_rows(labels)
}

#' Plot recovered estimates against their generating values
#'
#' The picture a recovery table is read through: the generating value on
#' the x axis, the posterior median on the y axis, the credible interval
#' as a bar and an identity line to compare against. A model that
#' recovers its parameters puts the points on the line; a model that is
#' biased puts them on a line beside it, and a model that cannot
#' distinguish them puts them on a cloud.
#'
#' For a correlation recovery from [recover_correlations()] the x axis is
#' the in-sample correlation of the simulated subjects (`truth =
#' "sample"`) or the generating correlation (`truth = "true"`), and the
#' points are coloured by estimator.
#'
#' For a cross-check from [cross_check()] the x axis is the reference,
#' the fit's interval runs vertically and the reference's own interval,
#' where it has one, runs horizontally. Panels default to one per term at
#' the subject level and to a single panel at the population level, where
#' each term contributes one point; `facet_by` overrides either.
#'
#' @param x A `bmmtools_recovery` object from [recover()] or
#'   [recover_subjects()], a `bmmtools_cor_recovery` object from
#'   [recover_correlations()], or a `bmmtools_cross_check` object from
#'   [cross_check()].
#' @param facet_by A column name to make panels from, or `NULL` for a
#'   single panel. Defaults to `"term"`: parameters usually live on
#'   scales too different to share an axis.
#' @param color_by A column name to colour points by, or `NULL`. `NULL`
#'   by default for a recovery object and `"estimator"` for a correlation
#'   recovery.
#' @param intervals Draw the intervals as bars.
#' @param identity_line Draw the line where the estimate equals the
#'   generating value.
#' @param scales Passed to [ggplot2::facet_wrap()]. `"free"` by default,
#'   for the same reason `facet_by` is.
#' @param annotate Label each panel with the Pearson correlation and
#'   Lin's concordance from [summary()][summary.bmmtools_recovery()] of
#'   the rows in that panel. Each panel must then hold a single term,
#'   level and condition; filter the object first if it does not.
#' @param ... Not used.
#'
#' @return A `ggplot` object.
#'
#' @examplesIf rlang::is_installed("ggplot2")
#' # subject-level recovery in one cell of the example grid
#' cell <- dplyr::filter(
#'   recovery_mixture2p,
#'   level == "subject", condition == "row-4"
#' )
#' plot_recovery(cell)
#' plot_recovery(cell, annotate = TRUE)
#'
#' # population-level estimates, coloured by design cell
#' population <- dplyr::filter(recovery_mixture2p, level == "population")
#' plot_recovery(population, color_by = "condition")
#'
#' @export
plot_recovery <- function(x, ...) {
  UseMethod("plot_recovery")
}

#' @rdname plot_recovery
#' @export
plot_recovery.default <- function(x, ...) {
  cli::cli_abort(
    "{.arg x} must be a {.cls bmmtools_recovery} object, \\
     not {.obj_type_friendly {x}}."
  )
}

#' @rdname plot_recovery
#' @export
plot_recovery.bmmtools_recovery <- function(x,
                                            facet_by = "term",
                                            color_by = NULL,
                                            intervals = TRUE,
                                            identity_line = TRUE,
                                            scales = "free",
                                            annotate = FALSE,
                                            ...) {
  rlang::check_dots_empty()
  rlang::check_installed("ggplot2", "to plot a recovery object.")
  check_recovery_contract(x)
  check_plot_column(x, facet_by, "facet_by")
  check_plot_column(x, color_by, "color_by")
  if (!rlang::is_bool(annotate)) {
    cli::cli_abort(
      "{.arg annotate} must be {.code TRUE} or {.code FALSE}, \\
       not {.obj_type_friendly {annotate}}."
    )
  }

  scale <- attr(x, "scale")
  if (is.null(scale)) scale <- unique(x$scale)

  mapping <- if (is.null(color_by)) {
    ggplot2::aes(x = .data$true_value, y = .data$estimate)
  } else {
    ggplot2::aes(
      x = .data$true_value, y = .data$estimate,
      colour = .data[[color_by]]
    )
  }

  p <- ggplot2::ggplot(tibble::as_tibble(x), mapping)

  # the identity line goes first so the data sits on top of it
  if (isTRUE(identity_line)) {
    p <- p + ggplot2::geom_abline(
      slope = 1, intercept = 0, linetype = "dashed", colour = "grey40"
    )
  }
  if (isTRUE(intervals)) {
    p <- p + ggplot2::geom_linerange(
      ggplot2::aes(ymin = .data$ci_low, ymax = .data$ci_high),
      alpha = 0.5
    )
  }
  p <- p + ggplot2::geom_point()

  if (annotate) {
    p <- p + ggplot2::geom_text(
      data = recovery_panel_labels(x, facet_by),
      mapping = ggplot2::aes(label = .data$label),
      x = -Inf, y = Inf, hjust = 0, vjust = 1.2,
      size = 3.2, lineheight = 0.9, inherit.aes = FALSE
    )
  }

  if (!is.null(facet_by)) {
    p <- p + ggplot2::facet_wrap(facet_by, scales = scales)
  }

  p +
    ggplot2::labs(
      x = scale_label(scale, "Generating value"),
      y = scale_label(scale, "Posterior median"),
      colour = color_by
    ) +
    ggplot2::theme_bw()
}

#' @rdname plot_recovery
#' @param truth For a correlation recovery, the value on the x axis:
#'   `"sample"`, the in-sample correlation, or `"true"`, the generating
#'   one. Rows without that value (a nonzero `true_value` on the natural
#'   scale is `NA`) are left out with a message.
#' @export
plot_recovery.bmmtools_cor_recovery <- function(x,
                                                truth = c("sample", "true"),
                                                facet_by = "term",
                                                color_by = "estimator",
                                                intervals = TRUE,
                                                identity_line = TRUE,
                                                ...) {
  rlang::check_dots_empty()
  rlang::check_installed("ggplot2", "to plot a recovery object.")
  check_cor_recovery_contract(x)
  truth <- rlang::arg_match(truth)
  check_plot_column(x, facet_by, "facet_by")
  check_plot_column(x, color_by, "color_by")

  column <- if (identical(truth, "sample")) "sample_value" else "true_value"
  scale <- attr(x, "scale") %||% unique(x$scale)
  data <- tibble::as_tibble(x)
  absent <- is.na(data[[column]])
  if (all(absent)) {
    cli::cli_abort(
      c(
        "No row has a {.field {column}} to plot.",
        i = "On the natural scale a nonzero generating correlation has no \\
             {.field true_value}; use {.code truth = \"sample\"}."
      )
    )
  }
  if (any(absent)) {
    cli::cli_inform(
      "{sum(absent)} row{?s} without a {.field {column}} {?is/are} not drawn."
    )
    data <- data[!absent, ]
  }

  mapping <- if (is.null(color_by)) {
    ggplot2::aes(x = .data[[column]], y = .data$estimate)
  } else {
    ggplot2::aes(
      x = .data[[column]], y = .data$estimate, colour = .data[[color_by]]
    )
  }
  p <- ggplot2::ggplot(data, mapping)
  if (isTRUE(identity_line)) {
    p <- p + ggplot2::geom_abline(
      slope = 1, intercept = 0, linetype = "dashed", colour = "grey40"
    )
  }
  if (isTRUE(intervals)) {
    p <- p + ggplot2::geom_linerange(
      ggplot2::aes(ymin = .data$ci_low, ymax = .data$ci_high),
      alpha = 0.5
    )
  }
  p <- p + ggplot2::geom_point()
  if (!is.null(facet_by)) {
    p <- p + ggplot2::facet_wrap(facet_by)
  }
  x_label <- if (identical(truth, "sample")) {
    "In-sample correlation"
  } else {
    "Generating correlation"
  }
  p +
    ggplot2::labs(
      x = scale_label(scale, x_label),
      y = scale_label(scale, "Estimated correlation"),
      colour = color_by
    ) +
    ggplot2::theme_bw()
}

#' @rdname plot_recovery
#' @export
plot_recovery.bmmtools_cross_check <- function(x,
                                               facet_by = "term",
                                               color_by = NULL,
                                               intervals = TRUE,
                                               identity_line = TRUE,
                                               scales = "free",
                                               ...) {
  rlang::check_dots_empty()
  rlang::check_installed("ggplot2", "to plot a cross-check.")
  check_cross_check_contract(x)
  # At the population level each term is a single point, and a panel per
  # point says less than one panel holding the whole comparison. Asking
  # for `facet_by` overrides this; `missing()` is what tells the two
  # apart, as in extract_estimates().
  if (missing(facet_by) && !any(x$level == "subject")) facet_by <- NULL
  check_plot_column(x, facet_by, "facet_by")
  check_plot_column(x, color_by, "color_by")

  scale <- attr(x, "scale") %||% unique(x$scale)
  data <- tibble::as_tibble(x)

  mapping <- if (is.null(color_by)) {
    ggplot2::aes(x = .data$reference, y = .data$estimate)
  } else {
    ggplot2::aes(
      x = .data$reference, y = .data$estimate, colour = .data[[color_by]]
    )
  }

  p <- ggplot2::ggplot(data, mapping)
  if (isTRUE(identity_line)) {
    p <- p + ggplot2::geom_abline(
      slope = 1, intercept = 0, linetype = "dashed", colour = "grey40"
    )
  }
  if (isTRUE(intervals)) {
    p <- p + ggplot2::geom_linerange(
      ggplot2::aes(ymin = .data$ci_low, ymax = .data$ci_high),
      alpha = 0.5
    )
    # the reference's own interval runs along x, and only where it exists:
    # a layer over all-NA bounds would draw nothing and warn about it
    if (any(!is.na(data$ref_low) & !is.na(data$ref_high))) {
      p <- p + ggplot2::geom_linerange(
        data = data[!is.na(data$ref_low) & !is.na(data$ref_high), ],
        mapping = ggplot2::aes(xmin = .data$ref_low, xmax = .data$ref_high),
        alpha = 0.5
      )
    }
  }
  p <- p + ggplot2::geom_point()

  if (!is.null(facet_by)) {
    p <- p + ggplot2::facet_wrap(facet_by, scales = scales)
  }

  p +
    ggplot2::labs(
      x = scale_label(scale, "Reference"),
      y = scale_label(scale, "Posterior median"),
      colour = color_by
    ) +
    ggplot2::theme_bw()
}

# the prior check --------------------------------------------------------

#' The prior-predictive draws of every set, one row per draw and column
#'
#' @return A data frame `prior`, `response`, `draw`, `value`.
#' @noRd
prior_draws_long <- function(x, draws) {
  model <- attr(x, "model")
  sets <- attr(x, "sets") %||% unique(x$prior)
  pieces <- lapply(sets, function(set) {
    shaped <- yrep_list(attr(x, "draws")[[set]], model)
    per_response <- lapply(names(shaped), function(name) {
      yrep <- shaped[[name]]
      keep <- seq_len(min(draws, nrow(yrep)))
      tibble::tibble(
        prior = set,
        response = name,
        draw = paste0(set, "-", rep(keep, each = ncol(yrep))),
        value = as.double(t(yrep[keep, , drop = FALSE]))
      )
    })
    dplyr::bind_rows(per_response)
  })
  dplyr::bind_rows(pieces)
}

#' The observed response, when the model names a column that `data` has
#'
#' A prior check run over a design skeleton has no observed data to
#' compare against, and a placeholder response column would draw a line
#' that means nothing, so the overlay is dropped rather than faked.
#'
#' @noRd
prior_observed <- function(x) {
  data <- attr(x, "data")
  model <- attr(x, "model")
  columns <- unlist(model$resp_vars, use.names = FALSE)
  columns <- intersect(columns, names(data))
  if (length(columns) == 0L) {
    return(NULL)
  }
  pieces <- lapply(columns, function(name) {
    tibble::tibble(
      response = name, value = as.double(data[[name]])
    )
  })
  out <- dplyr::bind_rows(pieces)
  out[is.finite(out$value), , drop = FALSE]
}

#' Plot what the priors say the data should look like
#'
#' The picture a prior check is read through. `"density"` draws each
#' prior-predictive draw as its own thin line, so the spread *between*
#' draws stays visible where a single pooled density would hide it;
#' `"histogram"` is for a discrete response, where a density is
#' misleading; `"statistic"` plots the summary table instead of the
#' distribution.
#'
#' @param x A `bmmtools_prior_check` object from [prior_check()].
#' @param type `"density"`, `"histogram"` or `"statistic"`.
#' @param observed Overlay the observed response from the data the check
#'   was run over. Dropped without comment when the model's response
#'   column is not in that data, since a prior check over a design
#'   skeleton has nothing observed to show. Ignored by `"statistic"`.
#' @param draws Prior-predictive draws to show per set, capped at what
#'   the object holds.
#' @param facet_by A column of `x` to make panels from, or `NULL`.
#'   `"response"` by default, so a model with several observables does
#'   not share an axis.
#' @param ... Not used.
#'
#' @return A `ggplot` object.
#'
#' @examplesIf rlang::is_installed("ggplot2")
#' plot_prior_check(prior_check_sdt_yn, type = "histogram")
#' plot_prior_check(prior_check_sdt_yn, type = "statistic")
#'
#' @export
plot_prior_check <- function(x,
                             type = c("density", "histogram", "statistic"),
                             observed = TRUE,
                             draws = 50,
                             facet_by = "response",
                             ...) {
  rlang::check_dots_empty()
  rlang::check_installed("ggplot2", "to plot a prior check.")

  if (!inherits(x, "bmmtools_prior_check")) {
    cli::cli_abort(
      "{.arg x} must be a {.cls bmmtools_prior_check} object, \\
       not {.obj_type_friendly {x}}."
    )
  }
  check_prior_check_contract(x)
  type <- rlang::arg_match(type)
  check_plot_column(x, facet_by, "facet_by")
  draws <- check_count(draws, "draws")

  p <- if (identical(type, "statistic")) {
    ggplot2::ggplot(
      tibble::as_tibble(x),
      ggplot2::aes(
        x = .data$value, y = .data$statistic, colour = .data$prior
      )
    ) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(x = "Value", y = NULL, colour = "Prior")
  } else {
    plot_prior_distribution(x, type, observed, draws)
  }

  if (!is.null(facet_by)) {
    p <- p + ggplot2::facet_wrap(facet_by, scales = "free")
  }
  p + ggplot2::theme_bw()
}

#' The distribution half of [plot_prior_check()]
#' @noRd
plot_prior_distribution <- function(x, type, observed, draws) {
  long <- prior_draws_long(x, draws)
  p <- ggplot2::ggplot(
    long, ggplot2::aes(x = .data$value, colour = .data$prior)
  )

  p <- if (identical(type, "histogram")) {
    p + ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(.data$density)),
      bins = 30, fill = NA, position = "identity"
    )
  } else {
    p + ggplot2::geom_density(
      ggplot2::aes(group = .data$draw),
      alpha = 0.3, linewidth = 0.2
    )
  }

  # the observed data go on top, so the prior-predictive layer stays
  # first and a colour scale reads off it
  if (isTRUE(observed)) {
    seen <- prior_observed(x)
    if (!is.null(seen) && nrow(seen) > 0L) {
      p <- p + ggplot2::geom_density(
        data = seen,
        mapping = ggplot2::aes(x = .data$value),
        inherit.aes = FALSE, linewidth = 1
      )
    }
  }

  p + ggplot2::labs(
    x = "Prior-predicted observable", y = "Density", colour = "Prior"
  )
}
