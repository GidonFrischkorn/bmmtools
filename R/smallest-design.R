# The read-off: the smallest design of a grid that meets a target
# (Milestone 24, D97, with the 2026-10-24 amendments to Q3).

#' How each target column is bounded below
#'
#' Rates take a Wilson score bound, means an MCSE bound, and the two
#' columns that are already lower ends of an interval are read as they
#' are. Higher is better for all of them.
#'
#' @noRd
target_kinds <- function() {
  c(
    detected = "rate", sign_recovery = "rate", precision_met = "rate",
    reliability_true = "mean", reliability_model = "mean",
    r_low = "bound", ccc_low = "bound"
  )
}

#' Columns refused as targets, with the column to use instead
#' @noRd
target_siblings <- function() {
  c(r = "r_low", ccc = "ccc_low", ci_width = "precision_met")
}

#' The Wilson score lower bound of a rate
#'
#' Unlike `p - z * MCSE`, it is below 1 at a rate of 1: ten successes in
#' ten give 0.722 at 95 %.
#'
#' @noRd
wilson_lower_bound <- function(p, n, level) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  centre <- p + z^2 / (2 * n)
  spread <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  (centre - spread) / (1 + z^2 / n)
}

#' The count a rate's Wilson bound uses (D97, amended 2026-10-24)
#'
#' The effective n `p (1 - p) / MCSE^2`: the count of intervals behind a
#' binomial MCSE, and at the subject level a count corrected for the
#' replications' shared posteriors, capped at every interval there is.
#' At a rate of 0 or 1, or a zero MCSE, the MCSE carries no count, and
#' the replications stand in.
#'
#' @noRd
effective_n <- function(p, mcse, x) {
  intervals <- ifelse(x$level == "subject", x$n * x$n_replications,
    x$n_replications
  )
  usable <- !is.na(mcse) & mcse > 0 & p > 0 & p < 1
  out <- as.double(x$n_replications)
  out[usable] <- pmin(
    p[usable] * (1 - p[usable]) / mcse[usable]^2,
    intervals[usable]
  )
  out
}

#' The lower end of a target column's Monte Carlo interval
#' @noRd
target_bound <- function(x, column, level) {
  value <- x[[column]]
  switch(target_kinds()[[column]],
    bound = value,
    mean = value - stats::qnorm(1 - (1 - level) / 2) *
      x[[paste0(column, "_mcse")]],
    rate = wilson_lower_bound(
      value, effective_n(value, x[[paste0(column, "_mcse")]], x), level
    )
  )
}

#' The smallest value of `along` met there and at every larger one
#' @noRd
smallest_kept <- function(along, met) {
  ord <- order(along)
  along <- along[ord]
  met <- met[ord] %in% TRUE
  kept <- rev(cumprod(rev(met))) == 1
  if (!any(kept)) along[NA_integer_][1L] else along[which(kept)[[1L]]]
}

#' Find the smallest design that meets a target
#'
#' Reads a grid summary design by design and returns, for each group of
#' rows, the smallest value of `along` (by default the number of trials)
#' at which a target is met **and stays met at every larger value in the
#' grid**. A design meets a target when the lower end of the target's
#' Monte Carlo interval clears it, so the answer accounts for how many
#' replications the grid ran.
#'
#' @param x A grid summary that carries the design columns
#'   (`summary()` of a [recovery_grid()] result), or such a recovery
#'   object, which is summarised first. Filter it first to the rows the
#'   question is about, for example `level == "subject"`.
#' @param target A named numeric vector of thresholds on summary columns,
#'   such as `c(r_low = 0.8)` or `c(reliability_model = 0.8,
#'   precision_met = 0.8)`. Higher is better for every column allowed:
#'   `detected`, `sign_recovery`, `precision_met`, `reliability_true`,
#'   `reliability_model`, `r_low` and `ccc_low`.
#' @param along The design column to minimise.
#' @param by Design columns that define separate groups, each with its own
#'   answer; `NULL` for none. Every grid column that varies besides
#'   `along` belongs here: each group must have one row per value of
#'   `along` (and term, level and estimator), or the call is refused.
#' @param level The mass of the Monte Carlo interval whose lower end must
#'   clear the target.
#' @param coverage_check If `TRUE`, a design counts as meeting a target
#'   only when every population row at that design also has its coverage
#'   within one MCSE of nominal, `abs(coverage - ci_level) <=
#'   coverage_mcse`. This is a calibration condition, not a target. It
#'   reads the population rows left in `x`, so keep them (their groups
#'   drop out of the answer when the target column is `NA` there), and a
#'   design with none fails the check.
#'
#' @return A `bmmtools_smallest_design` tibble with one row per `by` group,
#'   term, level, estimator and target: the `by` columns, `term`, `level`,
#'   `estimator`, `target` (the column), `threshold`, `smallest` (the
#'   certified value of `along`, or `NA`), `smallest_point` (the same
#'   rule applied to the point estimate, which shows what Monte Carlo error
#'   costs), `status`, and `replications_needed`, the replications a
#'   perfect rate would need to certify the threshold when that is why a
#'   design is undecided (`NA` otherwise). Printing it states each answer
#'   as a sentence.
#'
#' @details
#' **The lower end.** For a rate (`detected`, `sign_recovery`,
#' `precision_met`) it is the Wilson score bound, because the plain MCSE
#' of a rate is 0 at 0 and 1 and would certify a target from three
#' successes in three. Its count is the effective sample size
#' `p (1 - p) / MCSE^2`: the number of intervals when the MCSE is
#' binomial, and at the subject level a number corrected for subjects of
#' one fit not being independent, capped at the number of intervals.
#' At a rate of 0 or 1, or an MCSE of 0, the number of replications is
#' used, so a perfect rate certifies 0.80 only from 16 replications on.
#' For a reliability it is `estimate - z * MCSE`. `r_low` and `ccc_low`
#' are already lower ends and are read as they are; `level` does not
#' change them. `r`, `ccc` and `ci_width` have no Monte Carlo interval
#' and are refused, naming the column to use. `coverage` is refused too:
#' "high enough" is the wrong question for a coverage, which
#' `coverage_check` asks properly.
#'
#' **Monotone in the grid.** Noise that lets 30 trials pass and 100 fail
#' does not produce 30. A consequence is that several targets, or several
#' terms, are met together from the largest of their single answers,
#' which is how a condition on every parameter is read off.
#'
#' **No interpolation, no extrapolation.** `status` is `"met"`, `"met at
#' the grid's smallest design"` (the true minimum may be lower),
#' `"not reached in this grid"`, or `"undecided: too few replications"`
#' (the point estimate meets the target, the bound does not). A group
#' whose target column is `NA` at every design, such as a reliability on
#' population rows, is left out.
#'
#' @examples
#' summary(recovery_mixture2p) |>
#'   dplyr::filter(level == "subject") |>
#'   smallest_design(c(r_low = 0.8))
#'
#' @export
smallest_design <- function(x,
                            target,
                            along = "n_trials",
                            by = "n_subjects",
                            level = 0.95,
                            coverage_check = FALSE) {
  if (inherits(x, "bmmtools_recovery")) x <- summary(x)
  check_design_target(target)
  check_level_arg(level)
  x <- check_design_input(x, target, along, by, isTRUE(coverage_check))
  design <- c(by, along)
  coverage_ok <- if (isTRUE(coverage_check)) design_coverage_ok(x, design)

  keys <- c(by, "term", "level", "estimator")
  group <- interaction(x[keys], drop = TRUE, lex.order = TRUE)
  rows <- list()
  for (column in names(target)) {
    bound <- target_bound(x, column, level)
    for (g in levels(group)) {
      at <- which(group == g)
      if (all(is.na(x[[column]][at]))) next
      calibrated <- if (is.null(coverage_ok)) {
        rep(TRUE, length(at))
      } else {
        coverage_ok[design_key(x[at, ], design)] %in% TRUE
      }
      rows[[length(rows) + 1L]] <- tibble::tibble(
        x[at[[1L]], keys],
        read_off_group(
          x[[along]][at], x[[column]][at], bound[at], calibrated,
          threshold = target[[column]], rate = is_rate_target(column)
        ),
        target = column
      )
    }
  }
  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0L) {
    cli::cli_abort("{.arg x} has no rows with a value for the target.")
  }
  out <- design_status(out, level)
  out <- out[order(match(out$target, names(target))), , drop = FALSE]
  structure(
    out[c(
      keys, "target", "threshold", "smallest", "smallest_point",
      "status", "replications_needed"
    )],
    class = c("bmmtools_smallest_design", class(tibble::tibble())),
    along = along,
    level = level
  )
}

#' @noRd
is_rate_target <- function(column) {
  identical(target_kinds()[[column]], "rate")
}

#' The read-off of one group of designs
#'
#' `perfect` says whether the answer is undecided only because a perfect
#' rate has too few replications: every design from the point answer on
#' that fails its bound is calibrated and has a rate of exactly 1.
#'
#' @noRd
read_off_group <- function(along, point, bound, calibrated, threshold,
                           rate) {
  met <- bound >= threshold & calibrated
  met_point <- point >= threshold & calibrated
  smallest_point <- smallest_kept(along, met_point)
  failing <- !is.na(smallest_point) & along >= smallest_point &
    !met %in% TRUE
  tibble::tibble(
    threshold = threshold,
    smallest = smallest_kept(along, met),
    smallest_point = smallest_point,
    least = min(along),
    perfect = rate && any(failing) &&
      all(point[failing] %in% 1 & calibrated[failing])
  )
}

#' Turn the read-off into one of the four statuses
#'
#' A perfect rate of R replications has the Wilson bound R / (R + z^2),
#' so the replications it needs to certify a threshold t are
#' `ceiling(z^2 t / (1 - t))`.
#'
#' @noRd
design_status <- function(out, level) {
  out$status <- ifelse(
    !is.na(out$smallest),
    ifelse(out$smallest == out$least, "met at the grid's smallest design",
      "met"
    ),
    ifelse(!is.na(out$smallest_point), "undecided: too few replications",
      "not reached in this grid"
    )
  )
  z <- stats::qnorm(1 - (1 - level) / 2)
  needs <- out$perfect & out$status == "undecided: too few replications"
  out$replications_needed <- NA_integer_
  out$replications_needed[needs] <- as.integer(ceiling(
    z^2 * out$threshold[needs] / (1 - out$threshold[needs])
  ))
  out
}

#' Check what `smallest_design()` reads before it reads it
#'
#' Each group of `by`, term, level and estimator must have one row per
#' value of `along`, and none missing: a second row at the same value
#' (a grid column left out of `by`) would make the answer depend on row
#' order.
#'
#' @return `x` as a tibble.
#' @noRd
check_design_input <- function(x, target, along, by, coverage_check,
                               call = rlang::caller_env()) {
  if (!is.data.frame(x)) {
    cli::cli_abort(
      "{.arg x} must be a grid summary, not {.obj_type_friendly {x}}.",
      call = call
    )
  }
  if (!is.character(along) || length(along) != 1L) {
    cli::cli_abort("{.arg along} must be one column name.", call = call)
  }
  if (!is.null(by) && !is.character(by)) {
    cli::cli_abort(
      "{.arg by} must be column names or {.code NULL}.",
      call = call
    )
  }
  missing <- setdiff(c(by, along), names(x))
  if (length(missing) > 0L) {
    cli::cli_abort(
      c(
        "{.arg x} has no design column{?s} {.val {missing}}.",
        i = "A grid summary carries the grid's columns after \\
             {.val condition}: summarise a {.fn recovery_grid} result."
      ),
      call = call
    )
  }
  kinds <- target_kinds()
  with_mcse <- intersect(names(target), names(kinds)[kinds != "bound"])
  needed <- c(
    names(target), if (length(with_mcse) > 0L) paste0(with_mcse, "_mcse"),
    "term", "level", "estimator", "n", "n_replications",
    if (coverage_check) c("coverage", "coverage_mcse", "ci_level")
  )
  absent <- setdiff(needed, names(x))
  if (length(absent) > 0L) {
    cli::cli_abort(
      c(
        "{.arg x} has no column{?s} {.val {absent}}.",
        i = if ("precision_met" %in% absent) {
          "{.val precision_met} comes from {.code summary(x, precision = )}."
        }
      ),
      call = call
    )
  }
  x <- tibble::as_tibble(x)
  if (!is.numeric(x[[along]]) || anyNA(x[[along]])) {
    cli::cli_abort(
      "{.arg along} column {.val {along}} must be numeric with no \\
       missing values.",
      call = call
    )
  }
  keys <- c(by, "term", "level", "estimator")
  if (anyNA(x[keys])) {
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    gaps <- keys[colSums(is.na(x[keys])) > 0]
    cli::cli_abort(
      "{.arg x} has missing values in {.val {gaps}}.",
      call = call
    )
  }
  twice <- duplicated(x[c(keys, along)])
  if (any(twice)) {
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    first <- x[which(twice)[[1L]], c(keys, along)]
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    where <- paste(names(first), unlist(first), sep = " = ", collapse = ", ")
    cli::cli_abort(
      c(
        "{.arg x} has more than one row for {where}.",
        i = "Each design needs one row per term, level and estimator; add \\
             the grid column that tells them apart to {.arg by}."
      ),
      call = call
    )
  }
  x
}

#' @noRd
check_design_target <- function(target, call = rlang::caller_env()) {
  # `c(r_low = NA)` is logical: a missing threshold, not a wrong type
  if (is.logical(target) && all(is.na(target))) storage.mode(target) <- "double"
  named <- !is.null(names(target)) && all(names(target) != "")
  if (!is.numeric(target) || !named) {
    cli::cli_abort(
      c(
        "{.arg target} must be a named numeric vector of thresholds.",
        i = "For example {.code target = c(r_low = 0.8)}."
      ),
      call = call
    )
  }
  siblings <- target_siblings()
  refused <- intersect(names(target), names(siblings))
  if (length(refused) > 0L) {
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    use <- siblings[refused]
    cli::cli_abort(
      c(
        "{.val {refused}} {?has/have} no Monte Carlo interval to certify \\
         a design with.",
        i = "Use {.val {use}}."
      ),
      call = call
    )
  }
  if ("coverage" %in% names(target)) {
    cli::cli_abort(
      c(
        "{.val coverage} is not a target: a coverage should be near its \\
         nominal mass, not above a threshold.",
        i = "Use {.code coverage_check = TRUE}."
      ),
      call = call
    )
  }
  unknown <- setdiff(names(target), names(target_kinds()))
  if (length(unknown) > 0L) {
    cli::cli_abort(
      c(
        "{.arg target} names {.val {unknown}}, which {?is/are} not a \\
         target column.",
        i = "Target columns: {.val {names(target_kinds())}}."
      ),
      call = call
    )
  }
  if (anyNA(target) || any(!is.finite(target))) {
    cli::cli_abort("Every {.arg target} must be a finite number.", call = call)
  }
  if (anyDuplicated(names(target)) > 0L) {
    cli::cli_abort(
      "{.arg target} names {.val {names(target)[duplicated(names(target))]}} \\
       more than once.",
      call = call
    )
  }
  rates <- names(target)[vapply(names(target), is_rate_target, TRUE)]
  outside <- rates[target[rates] <= 0 | target[rates] >= 1]
  if (length(outside) > 0L) {
    cli::cli_abort(
      c(
        "A rate target must lie strictly between 0 and 1; {.val {outside}} \\
         does not.",
        i = "No number of replications certifies a rate of 1."
      ),
      call = call
    )
  }
  invisible(target)
}

#' @noRd
check_level_arg <- function(level, call = rlang::caller_env()) {
  usable <- is.numeric(level) && length(level) == 1L && !is.na(level)
  if (!usable || !(level > 0 && level < 1)) {
    cli::cli_abort(
      "{.arg level} must be a single number between 0 and 1.",
      call = call
    )
  }
  invisible(level)
}

#' @noRd
design_key <- function(x, design) {
  do.call(paste, c(as.list(x[design]), sep = "\r"))
}

#' Whether every population row of each design has calibrated coverage
#'
#' D83's second leg: `abs(coverage - ci_level) <= coverage_mcse` for every
#' population row at the design. A missing value fails.
#'
#' @return A logical vector named by design.
#' @noRd
design_coverage_ok <- function(x, design, call = rlang::caller_env()) {
  population <- x[x$level == "population", , drop = FALSE]
  if (nrow(population) == 0L) {
    cli::cli_abort(
      c(
        "{.code coverage_check = TRUE} needs population rows, and {.arg x} \\
         has none.",
        i = "Keep the population rows in {.arg x}; the check reads their \\
             coverage."
      ),
      call = call
    )
  }
  ok <- abs(population$coverage - population$ci_level) <=
    population$coverage_mcse
  vapply(split(ok, design_key(population, design)), function(v) {
    all(v %in% TRUE)
  }, logical(1))
}

#' Print the read-off as sentences
#'
#' @param x A `bmmtools_smallest_design` object.
#' @param ... Not used.
#' @return `format()` returns a character vector; `print()` returns `x`
#'   invisibly.
#' @export
format.bmmtools_smallest_design <- function(x, ...) {
  along <- attr(x, "along") %||% "the design"
  level <- attr(x, "level") %||% 0.95
  keys <- setdiff(
    names(x),
    c(
      "term", "level", "estimator", "target", "threshold", "smallest",
      "smallest_point", "status", "replications_needed"
    )
  )
  several <- length(unique(x$estimator)) > 1L
  lines <- vapply(seq_len(nrow(x)), function(i) {
    where <- paste(
      vapply(keys, function(k) paste0(k, " = ", x[[k]][[i]]), ""),
      collapse = ", "
    )
    head <- paste0(
      if (nzchar(where)) paste0(where, ", "),
      x$term[[i]], if (several) paste0(" [", x$estimator[[i]], "]"),
      " (", x$level[[i]], "), ", x$target[[i]], " >= ",
      format(x$threshold[[i]]), ": "
    )
    body <- switch(x$status[[i]],
      "met" = paste0("met from ", along, " = ", x$smallest[[i]], "."),
      "met at the grid's smallest design" = paste0(
        "met from ", along, " = ", x$smallest[[i]],
        ", the grid's smallest; the true minimum may be lower."
      ),
      "not reached in this grid" = "not reached in this grid.",
      paste0(
        "undecided: the point estimate meets it from ", along, " = ",
        x$smallest_point[[i]], ", the ", format(100 * level),
        " % lower bound does not",
        if (!is.na(x$replications_needed[[i]])) {
          paste0(
            "; a perfect rate needs at least ", x$replications_needed[[i]],
            " replications to certify it"
          )
        },
        "."
      )
    )
    paste0(head, body)
  }, character(1))
  c("<bmmtools_smallest_design>", lines)
}

#' @rdname format.bmmtools_smallest_design
#' @export
print.bmmtools_smallest_design <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}
