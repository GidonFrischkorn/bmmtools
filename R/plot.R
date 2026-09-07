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

#' Plot recovered estimates against their generating values
#'
#' The picture a recovery table is read through: the generating value on
#' the x axis, the posterior median on the y axis, the credible interval
#' as a bar and an identity line to compare against. A model that
#' recovers its parameters puts the points on the line; a model that is
#' biased puts them on a line beside it, and a model that cannot
#' distinguish them puts them on a cloud.
#'
#' @param x A `bmmtools_recovery` object from [recover()] or
#'   [recover_subjects()].
#' @param facet_by A column name to make panels from, or `NULL` for a
#'   single panel. Defaults to `"term"`: parameters usually live on
#'   scales too different to share an axis.
#' @param color_by A column name to colour points by, or `NULL`.
#' @param intervals Draw the credible intervals as bars.
#' @param identity_line Draw the line where the estimate equals the
#'   generating value.
#' @param scales Passed to [ggplot2::facet_wrap()]. `"free"` by default,
#'   for the same reason `facet_by` is.
#' @param ... Not used.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' recovery <- recover(fits, truth)
#' plot_recovery(recovery)
#' plot_recovery(recovery, facet_by = NULL, color_by = "term")
#' }
#'
#' @export
plot_recovery <- function(x,
                          facet_by = "term",
                          color_by = NULL,
                          intervals = TRUE,
                          identity_line = TRUE,
                          scales = "free",
                          ...) {
  rlang::check_dots_empty()
  rlang::check_installed("ggplot2", "to plot a recovery object.")

  if (!inherits(x, "bmmtools_recovery")) {
    cli::cli_abort(
      "{.arg x} must be a {.cls bmmtools_recovery} object, \\
       not {.obj_type_friendly {x}}."
    )
  }
  check_recovery_contract(x)
  check_plot_column(x, facet_by, "facet_by")
  check_plot_column(x, color_by, "color_by")

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
