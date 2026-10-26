# The design plot (Milestone 33, D102-D103): recovery across the cells of
# a grid, where plot_recovery() shows recovery within one set of rows.
#
# Every error bar has one stated source (see the roxygen of
# plot_design()). The bars that need the replications are computed from
# the rows here rather than added to summary(), so the summary contract
# and every stored summary stay as they are.

#' Okabe-Ito, in bmm-explorer's order, black last
#' @noRd
okabe_ito <- function() {
  c(
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#E69F00",
    "#F0E442", "#000000"
  )
}

#' Metrics read from the grid's cells rather than from its rows
#' @noRd
design_cell_metrics <- function() c("pass_rate", "fit_seconds")

#' Rates: drawn with a Wilson interval
#' @noRd
design_rate_metrics <- function() {
  c(
    "coverage", "coverage_50", "detected", "sign_recovery", "pass_rate",
    "rejection_rate", "false_positive_rate", "power", "coverage_sample"
  )
}

#' The metrics drawn when `metric` is `NULL`
#' @noRd
design_default_metrics <- function(level) {
  switch(level,
    subject = c("r", "reliability_model", "coverage"),
    correlation = c("bias", "coverage", "ci_width", "rejection_rate"),
    c("bias", "rmse", "coverage", "ci_width")
  )
}

#' Strip labels for the metrics bmmtools names; any other keeps its name
#' @noRd
design_metric_labels <- function() {
  c(
    bias = "Bias", rmse = "RMSE", mae = "MAE", coverage = "Coverage",
    coverage_50 = "Coverage (50 %)", ci_width = "CrI width",
    ci_width_50 = "CrI width (50 %)", contraction = "Contraction",
    r = "r", reliability_model = "Reliability (model)",
    reliability_true = "Reliability (truth)", z_mean = "Mean z",
    z_sd = "SD of z", pass_rate = "Share of fits passing",
    fit_seconds = "Seconds per fit", rejection_rate = "Share excluding 0",
    mean_estimate = "Mean estimate"
  )
}

#' The two ends of a Wilson score interval
#'
#' `wilson_lower_bound()` is the lower end; at a rate of 1 the upper end
#' is 1 and the lower n / (n + z^2).
#'
#' @noRd
wilson_interval <- function(p, n, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  spread <- z / (1 + z^2 / n) * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  list(lower = centre - spread, upper = centre + spread)
}

#' The design columns of a grid: its own columns, `condition` excluded
#' @noRd
design_columns_of <- function(grid) setdiff(names(grid), "condition")

#' Check `x_axis` and `color_by` against the design columns
#' @noRd
check_design_axes <- function(x_axis, color_by, design, call) {
  for (arg in c("x_axis", "color_by")) {
    value <- if (arg == "x_axis") x_axis else color_by
    if (!rlang::is_string(value) || !value %in% design) {
      cli::cli_abort(
        c(
          "{.arg {arg}} must name a design column of the grid.",
          i = "The design columns are {.val {design}}."
        ),
        call = call
      )
    }
  }
  if (identical(x_axis, color_by)) {
    cli::cli_abort(
      "{.arg color_by} must be another design column than {.arg x_axis}.",
      call = call
    )
  }
  invisible(NULL)
}

#' Refuse a design column named like a column of the plot's own data
#'
#' The plot data copy every design column next to their own columns, so
#' a design column of the same name would be overwritten and drawn as
#' something else, without an error.
#'
#' @noRd
check_design_names <- function(design, type, call) {
  own <- if (identical(type, "convergence")) {
    c("pass", "label", "dark")
  } else {
    c("metric", "value", "lower", "upper", "colour_group")
  }
  clash <- intersect(design, own)
  if (length(clash) > 0L) {
    cli::cli_abort(
      c(
        "The design column{?s} {.field {clash}} {?has/have} the name of a \\
         column the figure is drawn from.",
        i = "Rename {cli::qty(clash)}{?it/them} in the grid given to \\
             {.fn recovery_grid}."
      ),
      call = call
    )
  }
  invisible(NULL)
}

#' The replication-level part of the bars: what a summary has dropped
#'
#' One row per condition and term with the MCSE of the RMSE (the MSE's
#' MCSE from Morris et al., 2019, Table 6, through the delta method) and
#' of the mean CrI widths, contraction and z (SD over the replications'
#' means / sqrt(R)), and the count of rows behind each rate.
#'
#' @noRd
design_row_mcse <- function(rows) {
  # the rows a rate is a mean of: those where its indicator is known, and
  # for sign recovery those whose truth is not zero, as summary() counts
  indicators <- c(
    coverage = "covered", coverage_50 = "covered_50",
    coverage_sample = "covered_sample", rejection_rate = "excludes_zero"
  )
  rate_counts <- function(d) {
    known <- function(column) {
      if (column %in% names(d)) sum(!is.na(d[[column]])) else NA_integer_
    }
    counts <- lapply(indicators, known)
    counts$detected <- sum(!is.na(d$ci_low) & !is.na(d$ci_high))
    counts$sign_recovery <- sum(
      !is.na(d$true_value) & !is.na(d$estimate) & !d$true_value %in% 0
    )
    stats::setNames(counts, paste0("n_", names(counts)))
  }
  keys <- paste(rows$condition, rows$term, sep = "\r")
  pieces <- lapply(split(seq_len(nrow(rows)), keys), function(i) {
    d <- rows[i, ]
    e2 <- (d$estimate - d$true_value)^2
    e2 <- e2[!is.na(e2)]
    mse <- mean(e2)
    n <- length(e2)
    mse_mcse <- sqrt(sum((e2 - mse)^2) / (n * (n - 1)))
    per_rep <- function(values) {
      means <- tapply(values, d$replication, mean, na.rm = TRUE)
      means <- means[!is.na(means)]
      if (length(means) < 2L) {
        return(NA_real_)
      }
      stats::sd(means) / sqrt(length(means))
    }
    tibble::tibble(
      condition = d$condition[[1L]], term = d$term[[1L]],
      !!!rate_counts(d),
      rmse = if (n > 1L && mse > 0) mse_mcse / (2 * sqrt(mse)) else NA_real_,
      ci_width = per_rep(d$ci_high - d$ci_low),
      ci_width_50 = if (all(c("ci_low_50", "ci_high_50") %in% names(d))) {
        per_rep(d$ci_high_50 - d$ci_low_50)
      } else {
        NA_real_
      },
      contraction = if ("contraction" %in% names(d)) {
        per_rep(d$contraction)
      } else {
        NA_real_
      },
      z_mean = if ("z" %in% names(d)) per_rep(d$z) else NA_real_
    )
  })
  dplyr::bind_rows(pieces)
}

#' Lower and upper ends of one metric's bar, per summary row
#'
#' @param s Summary rows of one level.
#' @param level The level, which decides the count behind a rate.
#' @param row_mcse `design_row_mcse()` of the rows, or `NULL` for a
#'   summary.
#' @noRd
design_bounds <- function(metric, s, level, row_mcse) {
  value <- s[[metric]]
  none <- rep(NA_real_, nrow(s))
  at <- if (!is.null(row_mcse)) {
    match(
      paste(s$condition, s$term, sep = "\r"),
      paste(row_mcse$condition, row_mcse$term, sep = "\r")
    )
  }
  replications <- if (level %in% c("population", "sd")) {
    s$n
  } else {
    s$n_replications
  }
  if (metric %in% design_rate_metrics()) {
    # the count of independent intervals behind a rate: the rows that
    # entered it (one per replication), or at the subject level the count
    # behind the between-replication MCSE (D97's effective n). A summary
    # has dropped the rows; its complete pairs stand in.
    counted <- paste0("n_", metric)
    n <- if (identical(level, "subject")) {
      mcse <- s[[paste0(metric, "_mcse")]] %||% none
      effective_n(value, mcse, s)
    } else if (!is.null(row_mcse) && counted %in% names(row_mcse)) {
      row_mcse[[counted]][at]
    } else {
      replications
    }
    bounds <- wilson_interval(value, n)
    return(bounds)
  }
  plus_minus <- function(se) list(lower = value - se, upper = value + se)
  if (metric == "bias") {
    if (identical(level, "subject")) {
      return(list(lower = none, upper = none))
    }
    return(plus_minus(bias_mcse(s$bias, s$rmse, replications)))
  }
  if (metric == "r") {
    return(list(lower = s$r_low, upper = s$r_high))
  }
  if (metric %in% c("z_mean", "z_sd") && identical(level, "subject")) {
    # summary()'s z MCSEs count every subject as independent, which the
    # subjects of one replication, sharing its population posterior, are
    # not. The mean z takes the SD of the replications' means instead.
    # The SD of z gets no bar: it pools the subjects of every replication,
    # so it also holds the spread between the replications' means, and the
    # SD of the replications' own SDs is the MCSE of a smaller number
    if (metric == "z_mean" && !is.null(row_mcse)) {
      return(plus_minus(row_mcse$z_mean[at]))
    }
    return(list(lower = none, upper = none))
  }
  if (paste0(metric, "_mcse") %in% names(s)) {
    return(plus_minus(s[[paste0(metric, "_mcse")]]))
  }
  from_rows <- c("rmse", "ci_width", "ci_width_50", "contraction")
  # a subject RMSE pools the people of a replication, whose errors are
  # not independent, so it gets no bar
  subject_rmse <- metric == "rmse" && identical(level, "subject")
  if (!is.null(row_mcse) && metric %in% from_rows && !subject_rmse) {
    return(plus_minus(row_mcse[[metric]][at]))
  }
  list(lower = none, upper = none)
}

#' The long data of a profile: one row per design, term and metric
#' @noRd
design_profile_data <- function(s, metric, level, design, row_mcse,
                                term_order, call) {
  once <- is.character(metric) && length(metric) > 0L && !anyNA(metric) &&
    anyDuplicated(metric) == 0L
  if (!once) {
    cli::cli_abort(
      "{.arg metric} must name each metric once, or be {.code NULL}.",
      call = call
    )
  }
  numeric_columns <- names(s)[vapply(s, is.numeric, logical(1))]
  usable <- setdiff(numeric_columns, c(design, "n", "n_replications"))
  unknown <- setdiff(metric, usable)
  if (length(unknown) > 0L) {
    # the bounds, MCSEs and counts are columns too, but not metrics to read
    not_metrics <- "_(mcse|low|high)$|^(n_converged|ci_level)$"
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    offered <- usable[!grepl(not_metrics, usable)]
    cli::cli_abort(
      c(
        "{.arg metric} has {.val {unknown}}, which is not a metric of the \\
         {level} summary.",
        i = "Use one of {.val {offered}}, or {.val {design_cell_metrics()}}."
      ),
      call = call
    )
  }
  pieces <- lapply(metric, function(m) {
    bounds <- design_bounds(m, s, level, row_mcse)
    out <- tibble::as_tibble(s[c(design, "term")])
    out$metric <- m
    out$value <- s[[m]]
    out$lower <- bounds$lower
    out$upper <- bounds$upper
    out
  })
  out <- dplyr::bind_rows(pieces)
  out$term <- factor(out$term, levels = intersect(term_order, out$term))
  out$metric <- factor(out$metric, levels = metric)
  out
}

#' The cells' share of fits passing and mean seconds, per design
#' @noRd
design_cell_data <- function(cells, metric, design, call) {
  if ("fit_seconds" %in% metric && all(is.na(cells$fit_seconds))) {
    cli::cli_abort(
      "The grid's cells record no {.field fit_seconds} to draw.",
      call = call
    )
  }
  keys <- split(seq_len(nrow(cells)), cells$condition)
  rows <- lapply(keys, function(i) {
    d <- cells[i, ]
    passed <- d$converged %in% TRUE
    pass <- wilson_interval(mean(passed), length(passed))
    secs <- d$fit_seconds[!is.na(d$fit_seconds)]
    secs_se <- if (length(secs) > 1L) {
      stats::sd(secs) / sqrt(length(secs))
    } else {
      NA_real_
    }
    values <- list(
      pass_rate = c(mean(passed), pass$lower, pass$upper),
      fit_seconds = c(
        if (length(secs)) mean(secs) else NA_real_,
        if (length(secs)) mean(secs) + c(-1, 1) * secs_se else c(NA, NA)
      )
    )
    dplyr::bind_rows(lapply(metric, function(m) {
      out <- tibble::as_tibble(d[1L, design, drop = FALSE])
      out$metric <- m
      out$value <- values[[m]][[1L]]
      out$lower <- values[[m]][[2L]]
      out$upper <- values[[m]][[3L]]
      out
    }))
  })
  out <- dplyr::bind_rows(rows)
  out$metric <- factor(out$metric, levels = metric)
  out
}

#' A numeric design column with positive values gets a log axis
#' @noRd
design_x_scale <- function(values) {
  if (is.numeric(values) && all(values > 0, na.rm = TRUE)) {
    return(list(
      scale = ggplot2::scale_x_log10(breaks = sort(unique(values))),
      log = TRUE
    ))
  }
  list(scale = NULL, log = FALSE)
}

#' Draw a profile from its long data
#' @noRd
design_profile_plot <- function(data, x_axis, color_by, nominal, design,
                                call) {
  check_design_unique(
    data, intersect(c(x_axis, color_by, "term", "metric"), names(data)),
    design, call
  )
  groups <- sort(unique(data[[color_by]]))
  if (length(groups) > length(okabe_ito())) {
    cli::cli_abort(
      c(
        "{.arg color_by} has {length(groups)} values; a profile draws at \\
         most {length(okabe_ito())} lines.",
        i = "Filter the object to fewer values of {.field {color_by}}."
      ),
      call = call
    )
  }
  data$colour_group <- factor(data[[color_by]], levels = groups)
  x <- design_x_scale(data[[x_axis]])
  dodge <- ggplot2::position_dodge(width = if (x$log) 0.06 else 0.3)

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data[[x_axis]], y = .data$value, colour = .data$colour_group,
      group = .data$colour_group
    )
  )
  # the lines a metric is read against: no bias, the nominal coverage
  reference <- c(bias = 0, coverage = nominal, coverage_50 = 0.5)
  reference <- reference[names(reference) %in% levels(data$metric)]
  if (length(reference) > 0L) {
    p <- p + ggplot2::geom_hline(
      data = data.frame(
        metric = factor(names(reference), levels = levels(data$metric)),
        yintercept = unname(reference)
      ),
      mapping = ggplot2::aes(yintercept = .data$yintercept),
      colour = "grey40", linetype = "dashed", linewidth = 0.4,
      inherit.aes = FALSE
    )
  }
  p <- p +
    ggplot2::geom_line(position = dodge, linewidth = 0.7, na.rm = TRUE) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      position = dodge, linewidth = 0.6, na.rm = TRUE
    ) +
    ggplot2::geom_point(position = dodge, size = 2, na.rm = TRUE)

  labels <- design_metric_labels()
  metric_labeller <- function(values) {
    ifelse(values %in% names(labels), labels[values], values)
  }
  facets <- if ("term" %in% names(data)) {
    ggplot2::facet_grid(
      ggplot2::vars(.data$metric), ggplot2::vars(.data$term),
      scales = "free_y", switch = "y",
      labeller = ggplot2::labeller(metric = metric_labeller)
    )
  } else {
    ggplot2::facet_grid(
      ggplot2::vars(.data$metric),
      scales = "free_y", switch = "y",
      labeller = ggplot2::labeller(metric = metric_labeller)
    )
  }
  p <- p + facets + x$scale +
    ggplot2::scale_colour_manual(
      values = unname(okabe_ito()[seq_along(groups)]), name = color_by
    ) +
    ggplot2::labs(
      x = if (x$log) paste0(x_axis, " (log scale)") else x_axis, y = NULL
    )
  design_theme(p) +
    ggplot2::theme(strip.placement = "outside", legend.position = "top")
}

#' The shared theme of the design plots
#' @noRd
design_theme <- function(p) {
  p + ggplot2::theme_bw() + ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey95", colour = NA)
  )
}

#' Plot recovery across the designs of a grid
#'
#' The picture a grid summary is read through: how a metric changes with
#' the design. `type = "profile"` puts a metric on the y axis against
#' trials per subject (a log axis), one line per number of subjects and a
#' panel per term and metric. `type = "convergence"` draws one tile per
#' design with the share of fits that passed the convergence gate and the
#' checks the others failed. `type = "coverage"` pools the intervals of
#' every design per term and level, the coverage check that does not
#' depend on ten replications per design. [plot_recovery()] shows the
#' estimates of one design against their truth.
#'
#' A filter on `x` reaches the profiles of its rows and the model
#' correlations. The convergence tiles and the cell metrics draw the
#' designs `x` holds and count every replication of them: an errored fit
#' has no rows, so a filter on replications or on `converged` cannot be
#' read back from `x`, and a message says so. A design without rows
#' (every fit errored, or every fit in the replications `x` holds) is
#' still drawn, unless `x` was filtered by condition to designs whose
#' values it does not share; a message names it then. A filter that
#' removes every row of a design that ran in those replications removes
#' its tile.
#'
#' A figure refuses a design column named like a column it is drawn
#' from: `value`, `lower`, `upper`, `metric` or `colour_group`
#' in a profile, `pass`, `label` or `dark` in the tiles, `status`,
#' `converged`, `failed` or `fit_seconds` in the tiles and the cell
#' metrics, and a column of the correlation summary at the correlation
#' level. Rename it in the grid.
#'
#' @section Error bars:
#' Each metric's bar has one source, stated here so that a figure can be
#' read without the code:
#'
#' * `bias`: ±1 Monte Carlo standard error (MCSE), the SD of the errors
#'   over the replications / sqrt(R) (Morris et al., 2019, Table 6).
#'   None at the subject level, where the errors of one replication are
#'   not independent.
#' * `rmse`: ±1 MCSE of the MSE (Morris et al., 2019, Table 6) carried to
#'   the RMSE by the delta method, MCSE(MSE) / (2 RMSE). Not at the
#'   subject level.
#' * `ci_width`, `ci_width_50`, `contraction`: ±1 MCSE, the SD of the
#'   replications' means / sqrt(R).
#' * Rates (`coverage`, `coverage_50`, `detected`, `sign_recovery`,
#'   `rejection_rate`, `pass_rate`, ...): a 95 % Wilson score interval,
#'   which stays inside 0 to 1 and is not empty at a rate of 1. The count
#'   is the number of replications that entered the rate (for sign
#'   recovery, those whose truth is not zero); from a summary, which has
#'   dropped them, its complete pairs. At the subject level it is the
#'   count behind the between-replication MCSE that
#'   [summary()][summary.bmmtools_recovery()] reports (the effective n of
#'   [smallest_design()]), or the replications at a rate of 0 or 1. A
#'   subject rate of exactly 0 or 1 has a between-replication MCSE of 0,
#'   which carries no count, so the interval at the replications is
#'   conservative: it can be much wider than at a rate just below 1.
#' * `r`: the pooled interval `r_low` to `r_high` of
#'   [summary()][summary.bmmtools_recovery()].
#' * `reliability_model`, `reliability_true`, `reliability_gap`,
#'   `z_mean`, `z_sd`: ±1 of the MCSE column
#'   [summary()][summary.bmmtools_recovery()] reports. At the subject
#'   level that MCSE counts every subject as independent and understates,
#'   so `z_mean` takes ±1 MCSE, the SD of the replications' mean z /
#'   sqrt(R), and `z_sd` has no bar: it pools the subjects of every
#'   replication, so it also holds the spread between their mean z, which
#'   a bar from the replications' own SDs would leave out.
#' * `fit_seconds`: ±1 MCSE, the SD over fits / sqrt(fits).
#'
#' Any other metric is drawn without a bar. From a summary, which has
#' dropped the replications, `rmse`, the widths, contraction and the
#' subject-level `z_mean` have no bar.
#'
#' @param x A `bmmtools_recovery` object from [recovery_grid()], or its
#'   [summary()][summary.bmmtools_recovery()]. A summary draws profiles
#'   only.
#' @param type `"profile"`, `"convergence"` or `"coverage"`.
#' @param level For a profile, the level to draw: `"population"`,
#'   `"sd"`, `"subject"`, or `"correlation"` for the model correlations a
#'   grid run with `correlations` carries.
#' @param metric For a profile, the summary columns to draw, one panel
#'   row each, in this order. `NULL` draws bias, RMSE, coverage and CrI
#'   width for population means and SDs; r, the model-based reliability
#'   and coverage for subjects; bias, coverage, CrI width and the share of
#'   intervals excluding zero for correlations. `"pass_rate"` and
#'   `"fit_seconds"` come from the grid's cells and are drawn without
#'   `level` and without the other metrics.
#' @param x_axis,color_by Design columns of the grid: the one on the x
#'   axis and the one that makes the lines (the rows and columns of the
#'   convergence tiles).
#' @param ... Not used.
#'
#' @return A `ggplot` object. Its data holds the numbers drawn, with the
#'   bars in `lower` and `upper`.
#'
#' @references
#' Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
#' studies to evaluate statistical methods. *Statistics in Medicine,
#' 38*(11), 2074--2102. \doi{10.1002/sim.8086}
#'
#' @seealso [plot_recovery()] for the estimates of one design,
#'   [summary.bmmtools_recovery()] for the numbers.
#'
#' @examplesIf rlang::is_installed("ggplot2")
#' # population means across the four designs of the example grid
#' plot_design(recovery_mixture2p)
#'
#' # subject parameters, and r alone
#' plot_design(recovery_mixture2p, level = "subject", metric = "r")
#'
#' # how many fits passed the convergence gate, per design
#' plot_design(recovery_mixture2p, type = "convergence")
#'
#' @export
plot_design <- function(x, ...) {
  UseMethod("plot_design")
}

#' @rdname plot_design
#' @export
plot_design.default <- function(x, ...) {
  cli::cli_abort(
    "{.arg x} must be a {.cls bmmtools_recovery} object or its summary, \\
     not {.obj_type_friendly {x}}."
  )
}

#' @rdname plot_design
#' @export
plot_design.bmmtools_recovery <- function(x,
                                          type = c(
                                            "profile", "convergence",
                                            "coverage"
                                          ),
                                          level = "population",
                                          metric = NULL,
                                          x_axis = "n_trials",
                                          color_by = "n_subjects",
                                          ...) {
  rlang::check_dots_empty()
  rlang::check_installed("ggplot2", "to plot a recovery grid.")
  call <- rlang::current_env()
  check_recovery_contract(x)
  type <- rlang::arg_match(type)
  grid <- grid_design_columns(x)
  if (is.null(grid)) {
    cli::cli_abort(
      c(
        "{.arg x} carries no design grid.",
        i = "{.fn plot_design} draws the result of {.fn recovery_grid}; \\
             use {.fn plot_recovery} for one set of rows."
      ),
      call = call
    )
  }
  # pooled over every design, so neither axis applies
  if (type == "coverage") {
    return(design_coverage_plot(x, call))
  }
  design <- design_columns_of(grid)
  check_design_axes(x_axis, color_by, design, call)
  check_design_names(design, type, call)
  if (type == "convergence") {
    return(design_convergence_plot(x, grid, x_axis, color_by, call))
  }

  if (!is.null(metric) && any(metric %in% design_cell_metrics())) {
    check_cell_metrics(metric, call)
    cells <- design_cells(x, grid, call)
    data <- design_cell_data(cells, metric, design, call)
    return(design_profile_plot(data, x_axis, color_by, 0.95, design, call))
  }

  if (identical(level, "correlation")) {
    cors <- attr(x, "correlations")
    if (is.null(cors) || nrow(cors) == 0L) {
      cli::cli_abort(
        c(
          "{.arg x} carries no model correlations.",
          i = "Run {.fn recovery_grid} with {.arg correlations} to draw \\
               them."
        ),
        call = call
      )
    }
    rows <- tibble::as_tibble(cors)
    # the correlations sit in an attribute, which a filter on `x` does not
    # reach: draw the model's own estimate when there are several
    if (length(unique(rows$estimator)) > 1L) {
      if (!"model" %in% rows$estimator) {
        cli::cli_abort(
          "The model correlations hold the estimators \\
           {.val {unique(rows$estimator)}} and no {.val model}.",
          call = call
        )
      }
      cli::cli_inform(
        "Drawing the {.val model} estimate of the correlations; \\
         {.val {setdiff(unique(rows$estimator), 'model')}} {?is/are} left \\
         out."
      )
      rows <- rows[rows$estimator == "model", ]
    }
    # rows bound from two grids can hold both scales, which the summary
    # keeps apart and the figure cannot: draw the attribute's own scale
    scale <- attr(cors, "scale")
    scales <- unique(rows$scale)
    if (length(scales) > 1L) {
      if (!rlang::is_string(scale) || !scale %in% scales) {
        cli::cli_abort(
          c(
            "The model correlations hold the scales {.val {scales}}.",
            i = "Keep the correlations of one scale."
          ),
          call = call
        )
      }
      cli::cli_inform(
        "Drawing the model correlations on the {.val {scale}} scale; \\
         {.val {setdiff(scales, scale)}} {?is/are} left out."
      )
      rows <- rows[rows$scale == scale, ]
    }
    # the fits `x` still holds, so that a filter on replications reaches
    # the correlations too
    if (!is.null(rows[["replication"]]) && !is.null(x[["replication"]])) {
      held <- paste(rows$condition, rows$replication, sep = "\r") %in%
        paste(x$condition, x$replication, sep = "\r")
    } else {
      held <- rows$condition %in% x$condition
    }
    rows <- rows[held, ]
    s <- tibble::as_tibble(unclass(summary(new_bmmtools_cor_recovery(
      rows,
      scale = attr(cors, "scale"), ci_level = attr(cors, "ci_level")
    ))))
    # grid_design_columns() renames clashes with the recovery summary, not
    # with the correlation summary's own columns
    clash <- intersect(design, names(s))
    if (length(clash) > 0L) {
      cli::cli_abort(
        c(
          "The design column{?s} {.field {clash}} {?has/have} the name of \\
           a column of the correlation summary.",
          i = "Rename {cli::qty(clash)}{?it/them} in the grid given to \\
               {.fn recovery_grid}."
        ),
        call = call
      )
    }
    at <- match(s$condition, grid$condition)
    for (column in design) s[[column]] <- grid[[column]][at]
  } else {
    levels <- unique(x$level)
    if (!rlang::is_string(level) || !level %in% levels) {
      cli::cli_abort(
        c(
          "{.arg level} is {.val {level}}, which {.arg x} does not have.",
          i = "Its levels are {.val {levels}}, and {.val correlation} \\
               where the grid estimated correlations."
        ),
        call = call
      )
    }
    s <- tibble::as_tibble(unclass(summary(x)))
    s <- s[s$level == level, ]
    rows <- tibble::as_tibble(x)[x$level == level, ]
  }
  check_one_estimator(s, call)
  metric <- metric %||% design_default_metrics(level)
  data <- design_profile_data(
    s, metric, level, design, design_row_mcse(rows), unique(rows$term), call
  )
  design_profile_plot(
    data, x_axis, color_by, design_nominal(rows), design, call
  )
}

#' @rdname plot_design
#' @export
plot_design.bmmtools_recovery_summary <- function(x,
                                                  type = "profile",
                                                  level = "population",
                                                  metric = NULL,
                                                  x_axis = "n_trials",
                                                  color_by = "n_subjects",
                                                  ...) {
  rlang::check_dots_empty()
  rlang::check_installed("ggplot2", "to plot a recovery grid.")
  call <- rlang::current_env()
  if (!identical(type, "profile")) {
    cli::cli_abort(
      c(
        "A summary draws profiles only.",
        i = "Give the recovery object itself for \\
             {.code type = \"{type}\"}."
      ),
      call = call
    )
  }
  # a grid summary's design columns are those outside its contract that
  # hold one value per condition: by name, not position, so a relocated
  # column is still one, and a column added per term is a metric; so is
  # one asked for by `metric`
  design <- character()
  if ("condition" %in% names(x)) {
    extra <- setdiff(names(x), c(
      "condition", recovery_summary_columns(), precision_summary_columns()
    ))
    constant <- vapply(extra, function(column) {
      all(tapply(x[[column]], x$condition, function(v) {
        length(unique(v)) == 1L
      }))
    }, logical(1))
    design <- setdiff(extra[constant], metric)
  }
  if (length(design) == 0L) {
    cli::cli_abort(
      c(
        "{.arg x} has no design columns.",
        i = "Summarise the result of {.fn recovery_grid}, which carries \\
             its grid."
      ),
      call = call
    )
  }
  check_design_axes(x_axis, color_by, design, call)
  check_design_names(design, "profile", call)
  if (!is.null(metric) && any(metric %in% design_cell_metrics())) {
    cli::cli_abort(
      "{.val {design_cell_metrics()}} come from the grid's cells, which a \\
       summary does not have.",
      call = call
    )
  }
  s <- tibble::as_tibble(unclass(x))
  levels <- unique(s$level)
  if (!rlang::is_string(level) || !level %in% levels) {
    cli::cli_abort(
      "{.arg level} is {.val {level}}, which {.arg x} does not have; its \\
       levels are {.val {levels}}.",
      call = call
    )
  }
  s <- s[s$level == level, ]
  check_one_estimator(s, call)
  metric <- metric %||% design_default_metrics(level)
  data <- design_profile_data(
    s, metric, level, design, NULL, unique(s$term), call
  )
  design_profile_plot(data, x_axis, color_by, design_nominal(s), design, call)
}

#' The nominal mass of a summary's intervals, 0.95 when it has none
#' @noRd
design_nominal <- function(s) {
  # `[[`: a correlation summary has no ci_level column; its rows do
  nominal <- unique(s[["ci_level"]][!is.na(s[["ci_level"]])])
  if (length(nominal) == 1L) nominal else 0.95
}

#' One estimator per design and term, or the lines would overlap
#' @noRd
check_one_estimator <- function(s, call) {
  estimators <- unique(s$estimator)
  if (length(estimators) > 1L) {
    cli::cli_abort(
      c(
        "The rows hold {length(estimators)} estimators \\
         ({.val {estimators}}).",
        i = "Filter {.field estimator} to one first."
      ),
      call = call
    )
  }
  invisible(NULL)
}

#' The cell metrics are drawn alone
#' @noRd
check_cell_metrics <- function(metric, call) {
  if (!all(metric %in% design_cell_metrics())) {
    cli::cli_abort(
      c(
        "{.val {intersect(metric, design_cell_metrics())}} come{?s} from \\
         the grid's cells and cannot share a figure with \\
         {.val {setdiff(metric, design_cell_metrics())}}.",
        i = "Draw them in a figure of their own."
      ),
      call = call
    )
  }
  invisible(NULL)
}

#' The grid's cells joined to the design, or an error naming what is
#' missing
#'
#' Cells are kept by condition, every replication of each: an errored
#' fit has no rows, so which replications a filter meant to keep cannot
#' be read from `x` (a message says so). A design with no rows is kept
#' when its fits ran only in replications `x` does not hold (a filter on
#' replications removed them) or when no fit of it ran; `x` counts as
#' filtered by condition when a design ran in its replications and has
#' no rows, and an errored design is then kept only where each of its
#' design values is one a kept design has, with a message naming those
#' left out.
#'
#' Only the cell columns read here are kept. A design column named like
#' one of them is refused, because the figures read both under that
#' name. A component grid's cells have no `failed` or `fit_seconds`;
#' they are `NA`.
#'
#' @noRd
design_cells <- function(x, grid, call) {
  cells <- attr(x, "cells")
  if (is.null(cells) || nrow(cells) == 0L) {
    cli::cli_abort(
      "{.arg x} carries no cells, so there is nothing about its fits to \\
       draw.",
      call = call
    )
  }
  read <- c("status", "converged", "failed", "fit_seconds")
  clash <- intersect(design_columns_of(grid), read)
  if (length(clash) > 0L) {
    cli::cli_abort(
      c(
        "The grid has a design column named {.field {clash}}, which the \\
         cells also use.",
        i = "The cells cannot be drawn by design; the profiles of the \\
             summary metrics still can."
      ),
      call = call
    )
  }
  cells <- tibble::as_tibble(cells)
  ok <- (cells[["status"]] %||% rep("ok", nrow(cells))) %in% "ok"
  # by condition: an errored fit has no rows, so which replications a
  # filter on `x` meant to keep cannot be read from `x`, and every
  # replication of a design drawn is counted
  shown <- cells$condition %in% x$condition
  in_reps <- rep(TRUE, nrow(cells))
  if (!is.null(cells[["replication"]]) && !is.null(x[["replication"]])) {
    in_reps <- cells$replication %in% x$replication
  }
  # a design whose fits ran only in replications `x` does not hold lost
  # its rows to a filter on replications, not on conditions
  rowless <- !shown & cells$condition %in% cells$condition[ok] &
    !cells$condition %in% cells$condition[ok & in_reps]
  # a design with no fit that ran has no rows to filter
  errored <- !cells$condition %in% cells$condition[ok]
  if (any(errored) && any(!shown & ok & in_reps)) {
    # `x` was filtered by condition: keep an errored design only where
    # each of its design values is one a kept design has, and say which
    # it leaves out, since the filter could not name them
    at <- match(cells$condition, grid$condition)
    held <- grid[grid$condition %in% x$condition, , drop = FALSE]
    out_of_view <- errored
    for (column in design_columns_of(grid)) {
      errored <- errored & grid[[column]][at] %in% held[[column]]
    }
    left_out <- unique(cells$condition[out_of_view & !errored])
    if (length(left_out) > 0L) {
      cli::cli_inform(
        "Leaving out {.val {left_out}}: every fit there errored, and its \\
         design is not among those {.arg x} was filtered to."
      )
    }
  }
  kept <- shown | rowless | errored
  if (any(kept & ok & !in_reps)) {
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    drawn <- unique(cells$condition[rowless])
    cli::cli_inform(c(
      "The cells count every replication of a design; a filter on \\
       {.field replication} does not reach them.",
      i = if (length(drawn) > 0L) {
        "{.val {drawn}} {?has/have} no rows in {.arg x}, its fits in those \\
         replications having errored, and {?is/are} drawn."
      }
    ))
  }
  cells <- cells[kept, , drop = FALSE]
  out <- tibble::tibble(
    condition = cells$condition,
    status = cells[["status"]] %||% rep("ok", nrow(cells)),
    converged = cells[["converged"]] %||% rep(NA, nrow(cells)),
    failed = cells[["failed"]] %||% rep(NA_character_, nrow(cells)),
    fit_seconds = cells[["fit_seconds"]] %||% rep(NA_real_, nrow(cells))
  )
  # the grid, not the cells, carries the design
  at <- match(out$condition, grid$condition)
  for (column in design_columns_of(grid)) out[[column]] <- grid[[column]][at]
  out
}

#' Two designs that share the axes would be drawn as one point or tile
#' @noRd
check_design_unique <- function(data, keys, design, call) {
  if (anyDuplicated(data[keys]) == 0L) {
    return(invisible(NULL))
  }
  # nolint next: object_usage_linter. Used by cli's glue interpolation.
  others <- setdiff(design, keys)
  if (length(others) == 0L) {
    # the same design twice: the rows differ in something that is not a
    # design column, such as a scale or an estimator
    cli::cli_abort(
      c(
        "Several rows share the same {.field {keys}}.",
        x = "No design column tells them apart, so they would be drawn on \\
             top of each other.",
        i = "Keep one row of each: filter {.arg x}, for example to one \\
             {.field scale} or {.field estimator}."
      ),
      call = call
    )
  }
  cli::cli_abort(
    c(
      "Several designs share the same {.field {intersect(keys, design)}}.",
      x = "They differ in {.field {others}}, which the figure has no \\
           place for.",
      i = "Keep one value of each: filter {.arg x} by {.field condition} \\
           to the designs to draw."
    ),
    call = call
  )
}

#' The convergence tiles
#' @noRd
design_convergence_plot <- function(x, grid, x_axis, color_by, call) {
  cells <- design_cells(x, grid, call)
  design <- design_columns_of(grid)
  check_design_unique(
    cells[!duplicated(cells$condition), ], c(x_axis, color_by), design, call
  )
  checks <- c(
    ess_bulk = "bulk ESS", ess_tail = "tail ESS", rhat = "R-hat",
    divergent = "divergences"
  )
  tiles <- lapply(split(seq_len(nrow(cells)), cells$condition), function(i) {
    d <- cells[i, ]
    passed <- d$converged %in% TRUE
    failed <- d$failed[!passed & !is.na(d$failed) & nzchar(d$failed)]
    failed <- trimws(unlist(strsplit(failed, ",", fixed = TRUE)))
    counts <- table(factor(failed, levels = unique(failed)))
    why <- if (length(counts) > 0L) {
      named <- ifelse(names(counts) %in% names(checks),
        checks[names(counts)], names(counts)
      )
      paste0(named, ": ", as.integer(counts), collapse = "\n")
    }
    errors <- sum(!d$status %in% "ok")
    # a fit that ran but has no verdict: the gate assessed no variable
    unassessed <- sum(d$status %in% "ok" & is.na(d$converged))
    label <- paste(c(
      paste0(sum(passed), "/", length(passed), " pass"),
      why,
      if (errors > 0L) paste0("errors: ", errors),
      if (unassessed > 0L) paste0("not assessed: ", unassessed)
    ), collapse = "\n")
    out <- tibble::as_tibble(d[1L, design, drop = FALSE])
    out$pass <- mean(passed)
    out$label <- label
    out
  })
  data <- dplyr::bind_rows(tiles)
  data$dark <- data$pass >= 0.75
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = factor(.data[[x_axis]]), y = factor(.data[[color_by]]),
      fill = .data$pass
    )
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 2) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$label, colour = .data$dark),
      size = 3.2, lineheight = 0.95
    ) +
    ggplot2::scale_fill_gradient(
      low = "#DCEAF6", high = "#0B4F8A", limits = c(0, 1),
      name = "Share passing"
    ) +
    ggplot2::scale_colour_manual(
      values = c(`TRUE` = "white", `FALSE` = "grey10"), guide = "none"
    ) +
    ggplot2::labs(x = x_axis, y = color_by)
  design_theme(p) + ggplot2::theme(panel.grid = ggplot2::element_blank())
}

#' Coverage pooled over the designs, per level and term
#' @noRd
design_coverage_plot <- function(x, call) {
  rows <- tibble::as_tibble(x)[x$level %in% c("population", "sd"), ]
  if (nrow(rows) == 0L) {
    cli::cli_abort(
      "{.arg x} has no population or SD rows to pool.",
      call = call
    )
  }
  check_one_estimator(rows, call)
  nominal <- unique(rows$ci_level[!is.na(rows$ci_level)])
  if (length(nominal) != 1L) nominal <- attr(x, "ci_level") %||% 0.95
  keys <- paste(rows$level, rows$term, sep = "\r")
  pooled <- lapply(split(seq_len(nrow(rows)), keys), function(i) {
    d <- rows[i, ]
    one <- function(covered, mass) {
      covered <- covered[!is.na(covered)]
      if (length(covered) == 0L) {
        return(NULL)
      }
      bounds <- wilson_interval(mean(covered), length(covered))
      tibble::tibble(
        level = d$level[[1L]], term = d$term[[1L]],
        interval = paste0(format(100 * mass), " % intervals"),
        nominal = mass, value = mean(covered), lower = bounds$lower,
        upper = bounds$upper, n = length(covered)
      )
    }
    dplyr::bind_rows(one(d$covered, nominal), one(d$covered_50, 0.5))
  })
  data <- dplyr::bind_rows(pooled)
  data$term <- factor(data$term, levels = rev(unique(rows$term)))
  data$interval <- factor(
    data$interval,
    levels = unique(data$interval[order(data$nominal)])
  )
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$value, y = .data$term)) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = .data$nominal),
      colour = "grey40", linetype = "dashed"
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
      linewidth = 0.8, colour = okabe_ito()[[1L]]
    ) +
    ggplot2::geom_point(size = 2.5, colour = okabe_ito()[[1L]]) +
    ggplot2::facet_grid(
      ggplot2::vars(.data$level), ggplot2::vars(.data$interval),
      scales = "free_x"
    ) +
    ggplot2::labs(x = "Coverage pooled over every design", y = NULL)
  design_theme(p)
}
