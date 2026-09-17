# The recovery classes.
#
# `bmmtools_recovery` is a tibble subclass: one row per fit x parameter,
# carrying the apabayes `parameters` columns first and bmmtools's
# additions after. Subclassing a tibble
# means dplyr verbs keep working; a verb that drops a contract column
# drops the class, which is what dplyr_reconstruct() below implements.

#' The recovery contract
#'
#' Column names and their `typeof()`. `replication` is listed separately
#' because its type is the caller's choice --- an index, a name or a
#' grid identifier --- and only its presence is required.
#'
#' @noRd
recovery_contract <- function() {
  c(
    term = "character",
    estimate = "double",
    ci_low = "double",
    ci_high = "double",
    ci_method = "character",
    ci_level = "double",
    rhat = "double",
    ess_bulk = "double",
    ess_tail = "double",
    true_value = "double",
    bias = "double",
    covered = "logical",
    scale = "character",
    level = "character",
    id = "character",
    converged = "logical",
    condition = "character",
    estimator = "character"
  )
}

#' @noRd
recovery_contract_columns <- function() {
  c(names(recovery_contract()), "replication")
}

#' Contract columns a caller may leave out
#'
#' `converged` is unknown for a hand-built estimates tibble and
#' `condition` exists only for a grid; both are filled with `NA` rather
#' than demanded.
#'
#' @noRd
fill_optional_columns <- function(x) {
  if (!"converged" %in% names(x)) x$converged <- rep(NA, nrow(x))
  if (!"condition" %in% names(x)) x$condition <- rep(NA_character_, nrow(x))
  x
}

#' The columns `summary()` of a recovery object returns
#'
#' The apabayes request A1 columns, plus `n_replications`. A single
#' shape serves both [recover()] and [recover_subjects()] so that
#' downstream code reads one contract. At population level `n` and
#' `n_replications` differ only when a replication contributed no
#' complete pair, since one replication gives one pair per term.
#'
#' @noRd
recovery_summary_columns <- function() {
  c(
    "term", "estimator", "level", "scale", "n", "n_replications",
    "n_converged",
    "bias", "rmse", "coverage", "ci_width",
    "r", "r_low", "r_high", "rank_r",
    "ccc", "ccc_low", "ccc_high", "ccc_accuracy", "ccc_scale_shift",
    "ccc_location_shift", "calibration_slope", "truth_sd"
  )
}

#' Default a missing or `NA` estimator to `"bayes"`
#'
#' Deliberately *not* part of [fill_optional_columns()], which is shared
#' with `new_bmmtools_cor_recovery()`: the correlation contract requires
#' an `estimator` column of its own with a different vocabulary
#' (`"model"`, `"draws"`, `"point"`), and filling it in the shared helper
#' would stop a malformed correlation object from erroring and silently
#' stamp it `"bayes"` (milestone 8, trap 2).
#'
#' `NA` is filled as well as a missing column, because binding a labelled
#' tibble to an unlabelled one leaves `NA`, which would otherwise become
#' its own summary group.
#'
#' @noRd
fill_estimator <- function(x) {
  if (!"estimator" %in% names(x)) {
    return(rep("bayes", nrow(x)))
  }
  out <- as.character(x[["estimator"]])
  out[is.na(out)] <- "bayes"
  out
}

#' Construct a recovery object
#'
#' Validates the contract columns and their types, then attaches the
#' class and the attributes a print method and a plot need.
#'
#' @noRd
new_bmmtools_recovery <- function(x,
                                  scale,
                                  ci_level,
                                  call = NULL,
                                  error_call = rlang::caller_env()) {
  if (!is.data.frame(x)) {
    cli::cli_abort(
      "{.arg x} must be a data frame, not {.obj_type_friendly {x}}.",
      call = error_call
    )
  }

  x <- fill_optional_columns(x)
  # Not in fill_optional_columns(): that helper is shared with
  # new_bmmtools_cor_recovery(), whose contract *requires* an `estimator`
  # column with a different vocabulary ("model", "draws", "point").
  # Filling it there would stop a malformed correlation object from
  # erroring and stamp it "bayes" instead (milestone 8, trap 2).
  x$estimator <- fill_estimator(x)
  missing <- setdiff(recovery_contract_columns(), names(x))
  if (length(missing) > 0L) {
    cli::cli_abort(
      c(
        "A recovery object is missing the column{?s} {.val {missing}}.",
        i = "The contract is {.val {recovery_contract_columns()}}."
      ),
      call = error_call
    )
  }

  contract <- recovery_contract()
  for (column in names(contract)) {
    if (!identical(typeof(x[[column]]), contract[[column]])) {
      cli::cli_abort(
        "Column {.val {column}} must be {.cls {contract[[column]]}}, \\
         not {.cls {typeof(x[[column]])}}.",
        call = error_call
      )
    }
  }

  x <- tibble::as_tibble(x)[recovery_contract_columns()]
  structure(
    x,
    class = c("bmmtools_recovery", class(x)),
    scale = scale,
    ci_level = ci_level,
    call = call
  )
}

#' @noRd
new_bmmtools_recovery_summary <- function(x) {
  columns <- recovery_summary_columns()
  if ("condition" %in% names(x)) columns <- c("condition", columns)
  x <- tibble::as_tibble(x)[columns]
  structure(x, class = c("bmmtools_recovery_summary", class(x)))
}

# dplyr ------------------------------------------------------------------

#' Demote to a plain tibble when a contract column has gone
#'
#' The one place that decides whether a result still satisfies its
#' contract. A print, summary or plot method that assumed a missing
#' column would fail later and further from the cause.
#'
#' @param x The object a verb returned.
#' @param columns The contract it has to satisfy to keep its class.
#' @noRd
demote_if_incomplete <- function(x, columns) {
  if (!all(columns %in% names(x))) {
    return(tibble::as_tibble(x))
  }
  x
}

#' Keep the class only while the contract holds
#'
#' `dplyr::filter()` and friends return a recovery object; a `select()`
#' that drops a contract column returns a plain tibble.
#'
#' @param data,template See [dplyr::dplyr_reconstruct()].
#' @return `data`, as a recovery object when it still satisfies the
#'   contract and as a plain tibble otherwise.
#' @importFrom dplyr dplyr_reconstruct
#' @exportS3Method dplyr::dplyr_reconstruct
dplyr_reconstruct.bmmtools_recovery <- function(data, template) {
  if (!all(recovery_contract_columns() %in% names(data))) {
    return(tibble::as_tibble(data))
  }
  NextMethod()
}

#' Subset a recovery object
#'
#' `dplyr_reconstruct()` is not enough on its own. dplyr calls that
#' generic only when the result's class differs from the input's, and
#' `dplyr::select()` on a **tibble subclass** subsets through `[`, which
#' keeps the class; the generic is then never reached and a
#' column-dropping verb returns something still labelled a recovery
#' object with its contract broken. Measured on dplyr 1.2.1: without this
#' method `dplyr::select(x, "term", "estimate")` stays a
#' `bmmtools_recovery` and the next `print()` fails inside a metric with
#' a length mismatch. Demoting here is what makes the documented
#' behaviour --- a verb that drops a contract column drops the class ---
#' actually true.
#'
#' @param x A `bmmtools_recovery` object.
#' @param ... Passed to the tibble method.
#' @return A recovery object while the contract holds, a plain tibble
#'   once it does not.
#' @export
`[.bmmtools_recovery` <- function(x, ...) {
  demote_if_incomplete(NextMethod(), recovery_contract_columns())
}

#' Refuse to work on an object whose contract has been broken by hand
#'
#' Reached only when the class was attached to something that does not
#' satisfy the contract, since `[` and `dplyr_reconstruct()` demote. A
#' named error beats the length mismatch a metric would raise three
#' frames deeper.
#'
#' @noRd
check_recovery_contract <- function(x, call = rlang::caller_env()) {
  missing <- setdiff(recovery_contract_columns(), names(x))
  if (length(missing) > 0L) {
    cli::cli_abort(
      c(
        "This {.cls bmmtools_recovery} is missing the column{?s} \\
         {.val {missing}}.",
        i = "It was built by hand or altered in place; {.fn recover} and \\
             the {.pkg dplyr} verbs keep the contract or drop the class."
      ),
      call = call
    )
  }
  invisible(x)
}

# summary ----------------------------------------------------------------

#' Back-transform a Fisher z that may be undefined
#'
#' A perfect correlation puts z at infinity, and `tanh()` maps that back
#' to 1 or -1, which is right. Two of them with opposite signs average to
#' `NaN` instead: there is no combined correlation to report, and the
#' answer is `NA_real_` --- what the metric guard returns everywhere else
#' --- rather than `NaN`, which is a different missing value and prints
#' differently.
#'
#' @noRd
tanh_or_na <- function(z) {
  out <- tanh(z)
  if (length(out) == 1L && is.nan(out)) NA_real_ else out
}

#' Combine correlations across replications on Fisher's z scale
#'
#' Pooling subject-level pairs across replications mixes between-subject
#' with between-replication variance and inflates the correlation, so
#' each replication is correlated on its own and the results are combined
#' with inverse-variance weights (`n - 3`, the variance of Fisher's z).
#' With one replication this reduces exactly to `stats::cor.test()`'s
#' interval, and with equal subject counts to the plain average of the
#' transformed correlations.
#'
#' @noRd
fisher_z_combine <- function(r, n, ci_level = 0.95) {
  out <- list(r = NA_real_, r_low = NA_real_, r_high = NA_real_)
  keep <- !is.na(r) & !is.na(n)
  if (!any(keep)) {
    return(out)
  }
  r <- r[keep]
  n <- n[keep]

  z <- atanh(r)
  w <- pmax(n - 3, 0)
  if (sum(w) == 0) {
    out$r <- tanh_or_na(mean(z))
    return(out)
  }

  z_bar <- sum(w * z) / sum(w)
  out$r <- tanh_or_na(z_bar)
  # a perfect correlation puts z at infinity: the point estimate is 1 and
  # there is no interval, the same answer metric_r() gives
  if (!is.finite(z_bar)) {
    return(out)
  }
  se <- 1 / sqrt(sum(w))
  crit <- stats::qnorm(1 - (1 - ci_level) / 2)
  out$r_low <- tanh(z_bar - crit * se)
  out$r_high <- tanh(z_bar + crit * se)
  out
}

#' Combine concordance coefficients across replications on Lin's Z scale
#'
#' Each replication's `atanh(ccc)` is weighted by the inverse of Lin's
#' asymptotic variance (Lin, 1989, p. 259), which depends on `r` and the
#' location shift as well as on `n`, so `fisher_z_combine()`'s `n - 3`
#' weights would be wrong for it.
#'
#' All or nothing: a replication with a coefficient but no variance
#' (three subjects, `ccc` of exactly 0 or +/-1) is not given weight 0,
#' because those are exactly the extreme values and dropping them biases
#' the pooled value. The point estimate is then the unweighted Z mean
#' and there is no interval.
#'
#' @noRd
ccc_z_combine <- function(ccc, var_z, ci_level = 0.95) {
  out <- list(ccc = NA_real_, ccc_low = NA_real_, ccc_high = NA_real_)
  keep <- !is.na(ccc)
  if (!any(keep)) {
    return(out)
  }
  z <- atanh(ccc[keep])
  var_z <- var_z[keep]

  if (anyNA(var_z)) {
    out$ccc <- tanh_or_na(mean(z))
    return(out)
  }

  # every kept replication has a variance, so none has ccc = +/-1 and
  # z_bar is finite
  w <- 1 / var_z
  z_bar <- sum(w * z) / sum(w)
  out$ccc <- tanh(z_bar)
  crit <- stats::qnorm(1 - (1 - ci_level) / 2)
  se <- 1 / sqrt(sum(w))
  out$ccc_low <- tanh(z_bar - crit * se)
  out$ccc_high <- tanh(z_bar + crit * se)
  out
}

#' Replications whose fit passed the convergence gate
#'
#' `NA` when no row carries a verdict, never `0`: an unknown is not a
#' failure.
#'
#' @noRd
count_converged <- function(rows) {
  if (all(is.na(rows$converged))) {
    return(NA_integer_)
  }
  length(unique(rows$replication[rows$converged %in% TRUE]))
}

#' Mean of the values that are not missing
#' @noRd
mean_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}

#' Geometric mean of the values that are not missing
#'
#' For ratios, where `v` and `1/v` are equal departures from 1. `NA` when
#' nothing remains or any value is not positive: a negative calibration
#' slope has no geometric mean, and a plain mean would hide the sign.
#'
#' @noRd
geomean_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L || any(x <= 0)) NA_real_ else exp(mean(log(x)))
}

#' Metrics that do not depend on how rows are grouped
#' @noRd
summarise_errors <- function(rows) {
  list(
    bias = metric_bias(rows$estimate, rows$true_value),
    rmse = metric_rmse(rows$estimate, rows$true_value),
    coverage = metric_coverage(rows$true_value, rows$ci_low, rows$ci_high),
    ci_width = metric_ci_width(rows$ci_low, rows$ci_high)
  )
}

#' Population-level summary: one correlation across replications
#' @noRd
summarise_population <- function(rows) {
  r <- metric_r(rows$estimate, rows$true_value)
  ccc <- metric_ccc(rows$estimate, rows$true_value)
  c(
    list(
      n = as.double(r$n),
      n_replications = length(unique(rows$replication))
    ),
    summarise_errors(rows),
    list(
      r = r$r,
      r_low = r$r_low,
      r_high = r$r_high,
      rank_r = metric_rank_r(rows$estimate, rows$true_value),
      ccc = ccc$ccc,
      ccc_low = ccc$ccc_low,
      ccc_high = ccc$ccc_high,
      ccc_accuracy = ccc$accuracy,
      ccc_scale_shift = ccc$scale_shift,
      ccc_location_shift = ccc$location_shift,
      calibration_slope = ccc$calibration_slope,
      truth_sd = ccc$truth_sd
    )
  )
}

#' Subject-level summary: correlate within replication, then combine
#' @noRd
summarise_subject <- function(rows) {
  per <- lapply(split(rows, rows$replication), function(sub) {
    r <- metric_r(sub$estimate, sub$true_value)
    ccc <- metric_ccc(sub$estimate, sub$true_value)
    list(
      r = r$r,
      n = r$n,
      rank_r = metric_rank_r(sub$estimate, sub$true_value),
      ccc = ccc$ccc,
      var_z = ccc$var_z,
      ccc_accuracy = ccc$accuracy,
      ccc_scale_shift = ccc$scale_shift,
      ccc_location_shift = ccc$location_shift,
      calibration_slope = ccc$calibration_slope,
      truth_sd = ccc$truth_sd
    )
  })
  pull <- function(name) vapply(per, function(p) as.double(p[[name]]), 0)

  n <- pull("n")
  combined <- fisher_z_combine(pull("r"), n)
  ranks <- fisher_z_combine(pull("rank_r"), n)
  concordance <- ccc_z_combine(pull("ccc"), pull("var_z"))

  c(
    list(
      # subjects per replication; exact when the replications are
      # balanced, which is the case a simulation grid produces
      n = mean(n),
      n_replications = length(per)
    ),
    summarise_errors(rows),
    list(
      r = combined$r,
      r_low = combined$r_low,
      r_high = combined$r_high,
      rank_r = ranks$r,
      # Lin's Z, weighted by the inverse of his asymptotic variance;
      # the ratios combine on the log scale, where v and 1/v are
      # equally far from 1
      ccc = concordance$ccc,
      ccc_low = concordance$ccc_low,
      ccc_high = concordance$ccc_high,
      ccc_accuracy = mean_or_na(pull("ccc_accuracy")),
      ccc_scale_shift = geomean_or_na(pull("ccc_scale_shift")),
      ccc_location_shift = mean_or_na(pull("ccc_location_shift")),
      calibration_slope = geomean_or_na(pull("calibration_slope")),
      truth_sd = sqrt(mean_or_na(pull("truth_sd")^2))
    )
  )
}

#' Summarise a recovery object into per-parameter metrics
#'
#' One row per parameter and level, with these recovery metrics: bias,
#' RMSE, coverage, mean interval width, the Pearson correlation with a
#' Fisher-z interval, the Spearman correlation, and Lin's concordance
#' with its interval, its decomposition and a calibration slope (see
#' [recovery_ccc()] for how to read them).
#'
#' Correlation metrics are `NA`, never `0`, when fewer than three
#' complete pairs are available or when either side has no spread: `0`
#' would read as "no recovery" where the honest answer is "not estimable
#' from this design". The correlation interval is a 95% confidence
#' interval and does not follow `ci_level`, which is the mass of the
#' posterior interval and a different quantity.
#'
#' Population and SD rows are summarised across replications, one pair
#' per replication. Subject-level objects are summarised **within
#' replication and then combined**. Pooling subjects across replications
#' would mix between-subject with between-replication variance. `r` and `rank_r`
#' are combined on Fisher's z scale with weights `n - 3`; `ccc` on Lin's
#' Z scale with weights from his asymptotic variance, and without an
#' interval if any replication has no variance (three subjects, or a
#' coefficient of exactly 0 or 1); `ccc_scale_shift` and
#' `calibration_slope` as geometric means; `ccc_accuracy` and
#' `ccc_location_shift` as means; `truth_sd` as the root mean variance.
#' At this level `ccc = r * ccc_accuracy` holds only approximately.
#'
#' @param object A `bmmtools_recovery` object from [recover()] or
#'   [recover_subjects()].
#' @param ... Not used.
#'
#' @return A `bmmtools_recovery_summary` tibble with the columns `term`,
#'   `estimator`, `level`, `scale`, `n`, `n_replications`, `n_converged`,
#'   `bias`,
#'   `rmse`, `coverage`, `ci_width`, `r`, `r_low`, `r_high`, `rank_r`,
#'   `ccc`, `ccc_low`, `ccc_high`, `ccc_accuracy`, `ccc_scale_shift`,
#'   `ccc_location_shift`, `calibration_slope` and `truth_sd`, the
#'   standard deviation of the generating values. `r` and `ccc` both grow
#'   with `truth_sd` at a fixed measurement error, so compare them only
#'   at a similar spread.
#'
#' @details
#' Rows are grouped by `estimator` as well as by term and level, so two
#' estimators of the same parameter scored against one truth give two
#' rows rather than one pooled bias and RMSE. Their `n` may differ when
#' one of them failed on a subject the other estimated; [recover()] warns
#' when it does.
#'
#' `n_converged` is the number of replications whose fit passed
#' [check_convergence()], read from the `converged` column that
#' [extract_estimates()] fills. It is `NA` when no fit carried a
#' verdict, as with a hand-built estimates tibble: reporting it as equal
#' to `n` there would assert something that was never measured.
#'
#' @examples
#' summary(recovery_mixture2p)
#'
#' # one cell of the example grid, subject level only
#' recovery_mixture2p |>
#'   dplyr::filter(condition == "row-4", level == "subject") |>
#'   summary()
#'
#' @export
summary.bmmtools_recovery <- function(object, ...) {
  check_recovery_contract(object)
  if (nrow(object) == 0L) {
    return(new_bmmtools_recovery_summary(
      empty_recovery_summary()
    ))
  }

  # a grid carries its row label in `condition`; a plain recovery does
  # not, and the summary then has no such column
  by_condition <- !all(is.na(object$condition))
  keys <- paste(object$estimator, object$level, object$term, sep = "\r")
  if (by_condition) keys <- paste(object$condition, keys, sep = "\r")
  # factor(levels = unique(keys)) so the row order follows the order the
  # rows arrived in, rather than shifting to sort order the moment a
  # second estimator appears. The cor summary already does this.
  pieces <- lapply(split(
    seq_len(nrow(object)), factor(keys, levels = unique(keys))
  ), function(i) {
    rows <- object[i, ]
    body <- if (identical(rows$level[[1L]], "subject")) {
      summarise_subject(rows)
    } else {
      summarise_population(rows)
    }
    tibble::as_tibble(c(
      if (by_condition) list(condition = rows$condition[[1L]]),
      list(
        term = rows$term[[1L]],
        estimator = rows$estimator[[1L]],
        level = rows$level[[1L]],
        scale = rows$scale[[1L]]
      ),
      body,
      list(n_converged = count_converged(rows))
    ))
  })

  new_bmmtools_recovery_summary(dplyr::bind_rows(pieces))
}

#' @noRd
empty_recovery_summary <- function() {
  types <- c(
    term = "character", estimator = "character", level = "character",
    scale = "character",
    n = "double", n_replications = "integer", n_converged = "integer",
    bias = "double", rmse = "double", coverage = "double",
    ci_width = "double", r = "double", r_low = "double",
    r_high = "double", rank_r = "double", ccc = "double",
    ccc_low = "double", ccc_high = "double", ccc_accuracy = "double",
    ccc_scale_shift = "double", ccc_location_shift = "double",
    calibration_slope = "double", truth_sd = "double"
  )
  tibble::as_tibble(lapply(types, function(type) vector(type, 0L)))
}

#' @rdname summary.bmmtools_recovery
#' @export
summary.bmmtools_recovery_summary <- function(object, ...) {
  object
}

# printing ---------------------------------------------------------------

#' Format and print a recovery object
#'
#' The printed object leads with the scale it was scored on, because
#' bias, RMSE and the correlations are not invariant to the link
#' transform and a table of numbers without that label cannot be read.
#' Where the correlation metrics are `NA` the reason is printed as a
#' sentence rather than left as blank cells.
#'
#' @param x A `bmmtools_recovery` object.
#' @param ... Not used.
#'
#' @return `format()` returns a character vector; `print()` returns `x`
#'   invisibly.
#'
#' @export
format.bmmtools_recovery <- function(x, ...) {
  check_recovery_contract(x)
  scale <- attr(x, "scale")
  if (is.null(scale)) scale <- unique(x$scale)

  if (nrow(x) == 0L) {
    return(c("<bmmtools_recovery>", "No parameters scored."))
  }

  n_fits <- length(unique(x$replication))
  terms <- unique(x$term)
  levels <- unique(x$level)
  estimators <- unique(x$estimator)
  summarised <- summary(x)

  header <- c(
    "<bmmtools_recovery>",
    paste0("Scored on the ", scale, " scale."),
    paste0(
      n_fits, " fit", if (n_fits != 1L) "s", ", ",
      length(terms), " parameter", if (length(terms) != 1L) "s",
      ": ", paste(terms, collapse = ", "), "."
    ),
    paste0(
      "Level", if (length(levels) != 1L) "s", ": ",
      paste(levels, collapse = ", "), "."
    ),
    if ("sd" %in% levels) "SD rows are on the link scale.",
    # named only when there is a comparison to make; one estimator is the
    # ordinary case and saying "Estimators: bayes" is noise
    if (length(estimators) > 1L) {
      paste0("Estimators: ", paste(estimators, collapse = ", "), ".")
    },
    ""
  )

  note <- character(0)
  if (anyNA(summarised$r)) {
    note <- c(
      "",
      paste0(
        "r, rank_r and ccc are NA where they are not estimable: they ",
        "need at least 3 complete pairs and spread on both sides."
      )
    )
  }

  c(header, utils::capture.output(print(summarised)), note)
}

#' @rdname format.bmmtools_recovery
#' @export
print.bmmtools_recovery <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname format.bmmtools_recovery
#' @export
print.bmmtools_recovery_summary <- function(x, ...) {
  print(tibble::as_tibble(x), n = Inf, width = Inf)
  invisible(x)
}

# correlation recovery ---------------------------------------------------

#' The correlation-recovery contract
#'
#' As for `recovery_contract()`, `replication` is required but its type is
#' the caller's choice.
#'
#' @noRd
cor_recovery_contract <- function() {
  c(
    term = "character",
    var1 = "character",
    var2 = "character",
    estimator = "character",
    estimate = "double",
    ci_low = "double",
    ci_high = "double",
    ci_method = "character",
    ci_level = "double",
    rhat = "double",
    ess_bulk = "double",
    ess_tail = "double",
    true_value = "double",
    sample_value = "double",
    bias = "double",
    bias_sample = "double",
    covered = "logical",
    covered_sample = "logical",
    excludes_zero = "logical",
    scale = "character",
    n = "integer",
    converged = "logical",
    condition = "character"
  )
}

#' @noRd
cor_recovery_contract_columns <- function() {
  c(names(cor_recovery_contract()), "replication")
}

#' The columns `summary()` of a correlation recovery returns
#' @noRd
cor_recovery_summary_columns <- function() {
  c(
    "term", "estimator", "scale", "n_replications", "n_converged",
    "true_value", "sample_sd", "mean_estimate", "bias", "rmse",
    "bias_sample", "rmse_sample", "coverage", "coverage_sample", "ci_width",
    "rejection_rate", "false_positive_rate", "power",
    "r", "r_low", "r_high", "ccc", "ccc_low", "ccc_high"
  )
}

#' Construct a correlation-recovery object
#' @noRd
new_bmmtools_cor_recovery <- function(x,
                                      scale,
                                      ci_level,
                                      call = NULL,
                                      error_call = rlang::caller_env()) {
  if (!is.data.frame(x)) {
    cli::cli_abort(
      "{.arg x} must be a data frame, not {.obj_type_friendly {x}}.",
      call = error_call
    )
  }
  x <- fill_optional_columns(x)
  missing <- setdiff(cor_recovery_contract_columns(), names(x))
  if (length(missing) > 0L) {
    cli::cli_abort(
      c(
        "A correlation recovery object is missing the column{?s} \\
         {.val {missing}}.",
        i = "The contract is {.val {cor_recovery_contract_columns()}}."
      ),
      call = error_call
    )
  }
  contract <- cor_recovery_contract()
  for (column in names(contract)) {
    if (!identical(typeof(x[[column]]), contract[[column]])) {
      cli::cli_abort(
        "Column {.val {column}} must be {.cls {contract[[column]]}}, \\
         not {.cls {typeof(x[[column]])}}.",
        call = error_call
      )
    }
  }
  x <- tibble::as_tibble(x)[cor_recovery_contract_columns()]
  structure(
    x,
    class = c("bmmtools_cor_recovery", class(x)),
    scale = scale,
    ci_level = ci_level,
    call = call
  )
}

#' @noRd
check_cor_recovery_contract <- function(x, call = rlang::caller_env()) {
  missing <- setdiff(cor_recovery_contract_columns(), names(x))
  if (length(missing) > 0L) {
    cli::cli_abort(
      c(
        "This {.cls bmmtools_cor_recovery} is missing the column{?s} \\
         {.val {missing}}.",
        i = "It was built by hand or altered in place; \\
             {.fn recover_correlations} and the {.pkg dplyr} verbs keep \\
             the contract or drop the class."
      ),
      call = call
    )
  }
  invisible(x)
}

#' @noRd
#' @importFrom dplyr dplyr_reconstruct
#' @exportS3Method dplyr::dplyr_reconstruct
dplyr_reconstruct.bmmtools_cor_recovery <- function(data, template) {
  if (!all(cor_recovery_contract_columns() %in% names(data))) {
    return(tibble::as_tibble(data))
  }
  NextMethod()
}

#' Subset a correlation-recovery object
#'
#' As for [`[.bmmtools_recovery`]: `dplyr::select()` on a tibble subclass
#' subsets through `[` and never reaches `dplyr_reconstruct()`, so this
#' method is what drops the class once a contract column is gone.
#'
#' @param x A `bmmtools_cor_recovery` object.
#' @param ... Passed to the tibble method.
#' @return A correlation-recovery object while the contract holds, a
#'   plain tibble once it does not.
#' @export
`[.bmmtools_cor_recovery` <- function(x, ...) {
  demote_if_incomplete(NextMethod(), cor_recovery_contract_columns())
}

#' Standard deviation of the values that are not missing
#' @noRd
sd_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2L) NA_real_ else stats::sd(x)
}

#' The summary metrics of one term x estimator x scale group
#' @noRd
summarise_cor_rows <- function(rows) {
  truth <- unique(rows$true_value)
  rejection <- mean_or_na(as.double(rows$excludes_zero))
  r <- metric_r(rows$estimate, rows$sample_value)
  ccc <- metric_ccc(rows$estimate, rows$sample_value)
  list(
    n_replications = length(unique(rows$replication)),
    n_converged = count_converged(rows),
    true_value = if (length(truth) == 1L) as.double(truth) else NA_real_,
    sample_sd = sd_or_na(rows$sample_value),
    mean_estimate = mean_or_na(rows$estimate),
    bias = metric_bias(rows$estimate, rows$true_value),
    rmse = metric_rmse(rows$estimate, rows$true_value),
    bias_sample = metric_bias(rows$estimate, rows$sample_value),
    rmse_sample = metric_rmse(rows$estimate, rows$sample_value),
    coverage = metric_coverage(rows$true_value, rows$ci_low, rows$ci_high),
    coverage_sample = metric_coverage(
      rows$sample_value, rows$ci_low, rows$ci_high
    ),
    ci_width = metric_ci_width(rows$ci_low, rows$ci_high),
    rejection_rate = rejection,
    # an NA true_value is a nonzero correlation on the natural scale, so
    # it counts as "not 0"
    false_positive_rate = if (all(rows$true_value %in% 0)) {
      rejection
    } else {
      NA_real_
    },
    power = if (!any(rows$true_value %in% 0)) rejection else NA_real_,
    r = r$r,
    r_low = r$r_low,
    r_high = r$r_high,
    ccc = ccc$ccc,
    ccc_low = ccc$ccc_low,
    ccc_high = ccc$ccc_high
  )
}

#' @noRd
empty_cor_recovery_summary <- function() {
  columns <- cor_recovery_summary_columns()
  types <- stats::setNames(rep("double", length(columns)), columns)
  types[c("term", "estimator", "scale")] <- "character"
  types[c("n_replications", "n_converged")] <- "integer"
  tibble::as_tibble(lapply(types, function(type) vector(type, 0L)))
}

#' @noRd
new_cor_recovery_summary <- function(x) {
  columns <- cor_recovery_summary_columns()
  if ("condition" %in% names(x)) columns <- c("condition", columns)
  x <- tibble::as_tibble(x)[columns]
  structure(x, class = c("bmmtools_cor_recovery_summary", class(x)))
}

#' Summarise a correlation recovery into per-pair metrics
#'
#' One row per condition, pair, estimator and scale, summarised across
#' replications. Errors, coverage and the correlations are reported
#' against both truths: the generating correlation (`bias`, `rmse`,
#' `coverage`) and the correlation the simulated subjects had
#' (`bias_sample`, `rmse_sample`, `coverage_sample`, `r`, `ccc`).
#'
#' @param object A `bmmtools_cor_recovery` object from
#'   [recover_correlations()].
#' @param ... Not used.
#'
#' @return A `bmmtools_cor_recovery_summary` tibble with the columns
#'   `term`, `estimator`, `scale`, `n_replications`, `n_converged`,
#'   `true_value`, `sample_sd`, `mean_estimate`, `bias`, `rmse`,
#'   `bias_sample`, `rmse_sample`, `coverage`, `coverage_sample`,
#'   `ci_width`, `rejection_rate`, `false_positive_rate`, `power`, `r`,
#'   `r_low`, `r_high`, `ccc`, `ccc_low` and `ccc_high`, preceded by
#'   `condition` when the object came from a grid.
#'
#' @details
#' `true_value` is the generating correlation when it is the same in
#' every replication and `NA` otherwise. `sample_sd` is the standard
#' deviation of the in-sample correlations across replications.
#'
#' `rejection_rate` is the share of intervals that exclude zero. It is
#' reported as `false_positive_rate` when every `true_value` is 0 and as
#' `power` when none is. A `true_value` that is `NA`, which is what a
#' nonzero correlation becomes on the natural scale, counts as not 0, so
#' `power` is defined there. With a mix of zero and nonzero values both
#' are `NA`.
#'
#' `r` (with a Fisher-z interval) and `ccc` (Lin's concordance, see
#' [recovery_ccc()]) compare the estimates with `sample_value` across
#' replications, one pair per replication. They are `NA` with fewer than
#' three replications or when either side has no spread.
#'
#' @export
summary.bmmtools_cor_recovery <- function(object, ...) {
  check_cor_recovery_contract(object)
  if (nrow(object) == 0L) {
    return(new_cor_recovery_summary(empty_cor_recovery_summary()))
  }
  by_condition <- !all(is.na(object$condition))
  keys <- paste(object$term, object$estimator, object$scale, sep = "\r")
  if (by_condition) keys <- paste(object$condition, keys, sep = "\r")
  groups <- split(seq_len(nrow(object)), factor(keys, levels = unique(keys)))
  pieces <- lapply(groups, function(i) {
    rows <- object[i, ]
    tibble::as_tibble(c(
      if (by_condition) list(condition = rows$condition[[1L]]),
      list(
        term = rows$term[[1L]],
        estimator = rows$estimator[[1L]],
        scale = rows$scale[[1L]]
      ),
      summarise_cor_rows(rows)
    ))
  })
  new_cor_recovery_summary(dplyr::bind_rows(pieces))
}

#' @rdname summary.bmmtools_cor_recovery
#' @export
summary.bmmtools_cor_recovery_summary <- function(object, ...) {
  object
}

#' Format and print a correlation recovery
#'
#' The printed object leads with the scale and the estimators, because
#' the three estimators answer different questions and a correlation
#' table without that label cannot be read.
#'
#' @param x A `bmmtools_cor_recovery` object.
#' @param ... Not used.
#'
#' @return `format()` returns a character vector; `print()` returns `x`
#'   invisibly.
#'
#' @export
format.bmmtools_cor_recovery <- function(x, ...) {
  check_cor_recovery_contract(x)
  if (nrow(x) == 0L) {
    return(c("<bmmtools_cor_recovery>", "No correlations scored."))
  }
  scale <- attr(x, "scale") %||% unique(x$scale)
  estimators <- unique(x$estimator)
  terms <- unique(x$term)
  n_fits <- length(unique(x$replication))
  summarised <- summary(x)

  header <- c(
    "<bmmtools_cor_recovery>",
    paste0(
      "Scored on the ", paste(scale, collapse = ", "), " scale; estimator",
      if (length(estimators) != 1L) "s", ": ",
      paste(estimators, collapse = ", "), "."
    ),
    paste0(
      n_fits, " fit", if (n_fits != 1L) "s", ", ",
      length(terms), " correlation", if (length(terms) != 1L) "s",
      ": ", paste(terms, collapse = ", "), "."
    ),
    if (identical(scale, "natural")) {
      paste0(
        "On the natural scale true_value is 0 where the generating ",
        "correlation is 0 and NA otherwise; compare with sample_value."
      )
    },
    ""
  )
  note <- character(0)
  if (anyNA(summarised$r)) {
    note <- c(
      "",
      paste0(
        "r and ccc are NA where they are not estimable: they need at ",
        "least 3 replications and spread on both sides."
      )
    )
  }
  c(header, utils::capture.output(print(summarised)), note)
}

#' @rdname format.bmmtools_cor_recovery
#' @export
print.bmmtools_cor_recovery <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname format.bmmtools_cor_recovery
#' @export
print.bmmtools_cor_recovery_summary <- function(x, ...) {
  print(tibble::as_tibble(x), n = Inf, width = Inf)
  invisible(x)
}

# the cross-check --------------------------------------------------------

#' The cross-check contract
#'
#' The apabayes `parameters` columns first, then the reference and the
#' comparison. There is no
#' `replication` column: `cross_check()` takes one fit, and a list of
#' them is `recover()`'s shape.
#'
#' @noRd
cross_check_contract <- function() {
  c(
    term = "character",
    estimate = "double",
    ci_low = "double",
    ci_high = "double",
    ci_method = "character",
    ci_level = "double",
    rhat = "double",
    ess_bulk = "double",
    ess_tail = "double",
    reference = "double",
    ref_low = "double",
    ref_high = "double",
    source = "character",
    bias = "double",
    covered = "logical",
    overlap = "logical",
    scale = "character",
    level = "character",
    id = "character",
    converged = "logical"
  )
}

#' The columns `summary()` of a cross-check returns
#'
#' `bias`, `rmse` and `coverage` are `summary.bmmtools_recovery()`'s
#' words for the same quantities, so the three scorers read with one
#' vocabulary. `share_overlap` has no analogue there and keeps
#' its own name.
#'
#' @noRd
cross_check_summary_columns <- function() {
  c(
    "term", "level", "scale", "n", "n_converged",
    "bias", "rmse", "coverage", "share_overlap",
    "r", "r_low", "r_high", "ccc", "ccc_low", "ccc_high"
  )
}

#' Construct a cross-check object
#' @noRd
new_bmmtools_cross_check <- function(x,
                                     scale,
                                     ci_level,
                                     call = NULL,
                                     error_call = rlang::caller_env()) {
  if (!is.data.frame(x)) {
    cli::cli_abort(
      "{.arg x} must be a data frame, not {.obj_type_friendly {x}}.",
      call = error_call
    )
  }
  # a hand-built estimates tibble carries no verdict; an unknown is not a
  # failure, as for the recovery classes
  if (!"converged" %in% names(x)) x$converged <- rep(NA, nrow(x))

  contract <- cross_check_contract()
  missing <- setdiff(names(contract), names(x))
  if (length(missing) > 0L) {
    cli::cli_abort(
      c(
        "A cross-check object is missing the column{?s} {.val {missing}}.",
        i = "The contract is {.val {names(contract)}}."
      ),
      call = error_call
    )
  }
  for (column in names(contract)) {
    if (!identical(typeof(x[[column]]), contract[[column]])) {
      cli::cli_abort(
        "Column {.val {column}} must be {.cls {contract[[column]]}}, \\
         not {.cls {typeof(x[[column]])}}.",
        call = error_call
      )
    }
  }

  x <- tibble::as_tibble(x)[names(contract)]
  structure(
    x,
    class = c("bmmtools_cross_check", class(x)),
    scale = scale,
    ci_level = ci_level,
    call = call
  )
}

#' @noRd
check_cross_check_contract <- function(x, call = rlang::caller_env()) {
  missing <- setdiff(names(cross_check_contract()), names(x))
  if (length(missing) > 0L) {
    cli::cli_abort(
      c(
        "This {.cls bmmtools_cross_check} is missing the column{?s} \\
         {.val {missing}}.",
        i = "It was built by hand or altered in place; {.fn cross_check} \\
             and the {.pkg dplyr} verbs keep the contract or drop the \\
             class."
      ),
      call = call
    )
  }
  invisible(x)
}

#' @noRd
#' @importFrom dplyr dplyr_reconstruct
#' @exportS3Method dplyr::dplyr_reconstruct
dplyr_reconstruct.bmmtools_cross_check <- function(data, template) {
  if (!all(names(cross_check_contract()) %in% names(data))) {
    return(tibble::as_tibble(data))
  }
  NextMethod()
}

#' Subset a cross-check object
#'
#' As for [`[.bmmtools_recovery`]: `dplyr::select()` on a tibble subclass
#' subsets through `[` and never reaches `dplyr_reconstruct()`, so this
#' method is what drops the class once a contract column is gone.
#'
#' @param x A `bmmtools_cross_check` object.
#' @param ... Passed to the tibble method.
#' @return A cross-check object while the contract holds, a plain tibble
#'   once it does not.
#' @export
`[.bmmtools_cross_check` <- function(x, ...) {
  demote_if_incomplete(NextMethod(), names(cross_check_contract()))
}

#' Rows whose fit passed the convergence gate
#'
#' `count_converged()` counts replications, which a cross-check does not
#' have: it compares one fit. Counting rows keeps `n_converged` on the
#' same footing as `n`. `NA` when no row carries a verdict, never `0`.
#'
#' @noRd
count_converged_rows <- function(converged) {
  if (all(is.na(converged))) {
    return(NA_integer_)
  }
  sum(converged %in% TRUE)
}

#' @noRd
empty_cross_check_summary <- function() {
  columns <- cross_check_summary_columns()
  types <- stats::setNames(rep("double", length(columns)), columns)
  types[c("term", "level", "scale")] <- "character"
  types[["n_converged"]] <- "integer"
  tibble::as_tibble(lapply(types, function(type) vector(type, 0L)))
}

#' @noRd
new_cross_check_summary <- function(x) {
  x <- tibble::as_tibble(x)[cross_check_summary_columns()]
  structure(x, class = c("bmmtools_cross_check_summary", class(x)))
}

#' Summarise a cross-check into per-parameter metrics
#'
#' One row per term and level: the mean and root mean square of
#' `estimate - reference`, the share of fit intervals covering the
#' reference, the share of rows whose intervals overlap, and the
#' correlation and Lin's concordance between the two sides.
#'
#' `bias` here is a difference from a comparison value, not from a known
#' truth; see [cross_check()].
#'
#' @param object A `bmmtools_cross_check` object from [cross_check()].
#' @param ... Not used.
#'
#' @return A `bmmtools_cross_check_summary` tibble with the columns
#'   `term`, `level`, `scale`, `n`, `n_converged`, `bias`, `rmse`,
#'   `coverage`, `share_overlap`, `r`, `r_low`, `r_high`, `ccc`,
#'   `ccc_low` and `ccc_high`.
#'
#' @details
#' `r` and `ccc` come from the pairs of `estimate` and `reference` in the
#' group and are `NA`, never `0`, below three complete pairs or without
#' spread on either side. At the population level with one fit each term
#' contributes a single pair, so both are `NA` and the table is the
#' per-term difference and coverage; at the subject level they are the
#' comparison across subjects. `share_overlap` is `NA` for a reference
#' without intervals. `n_converged` counts the rows whose fit passed
#' [check_convergence()] and is `NA` when no row carries a verdict.
#'
#' @export
summary.bmmtools_cross_check <- function(object, ...) {
  check_cross_check_contract(object)
  if (nrow(object) == 0L) {
    return(new_cross_check_summary(empty_cross_check_summary()))
  }

  keys <- paste(object$level, object$term, sep = "\r")
  groups <- split(seq_len(nrow(object)), factor(keys, levels = unique(keys)))
  pieces <- lapply(groups, function(i) {
    rows <- object[i, ]
    r <- metric_r(rows$estimate, rows$reference)
    ccc <- metric_ccc(rows$estimate, rows$reference)
    tibble::as_tibble(list(
      term = rows$term[[1L]],
      level = rows$level[[1L]],
      scale = rows$scale[[1L]],
      n = as.double(r$n),
      n_converged = count_converged_rows(rows$converged),
      bias = metric_bias(rows$estimate, rows$reference),
      rmse = metric_rmse(rows$estimate, rows$reference),
      coverage = metric_coverage(rows$reference, rows$ci_low, rows$ci_high),
      share_overlap = mean_or_na(as.double(rows$overlap)),
      r = r$r,
      r_low = r$r_low,
      r_high = r$r_high,
      ccc = ccc$ccc,
      ccc_low = ccc$ccc_low,
      ccc_high = ccc$ccc_high
    ))
  })

  new_cross_check_summary(dplyr::bind_rows(pieces))
}

#' @rdname summary.bmmtools_cross_check
#' @export
summary.bmmtools_cross_check_summary <- function(object, ...) {
  object
}

#' Format and print a cross-check
#'
#' The header names the scale and where the reference came from, and says
#' that the reference is a comparison rather than a truth, because a
#' column called `bias` invites the other reading.
#'
#' @param x A `bmmtools_cross_check` object.
#' @param ... Not used.
#'
#' @return `format()` returns a character vector; `print()` returns `x`
#'   invisibly.
#'
#' @export
format.bmmtools_cross_check <- function(x, ...) {
  check_cross_check_contract(x)
  if (nrow(x) == 0L) {
    return(c("<bmmtools_cross_check>", "No parameters compared."))
  }
  scale <- attr(x, "scale") %||% unique(x$scale)
  terms <- unique(x$term)
  levels <- unique(x$level)
  sources <- unique(x$source)
  summarised <- summary(x)

  header <- c(
    "<bmmtools_cross_check>",
    paste0("Compared on the ", paste(scale, collapse = ", "), " scale."),
    paste0(
      length(terms), " parameter", if (length(terms) != 1L) "s",
      ": ", paste(terms, collapse = ", "), "."
    ),
    paste0(
      "Level", if (length(levels) != 1L) "s", ": ",
      paste(levels, collapse = ", "), "."
    ),
    paste0(
      "Reference source", if (length(sources) != 1L) "s", ": ",
      paste(sources, collapse = ", "), "."
    ),
    paste0(
      "The reference is a comparison, not a truth: bias is the signed ",
      "difference from it."
    ),
    ""
  )

  note <- character(0)
  if (anyNA(summarised$r)) {
    note <- c(
      "",
      paste0(
        "r and ccc are NA where they are not estimable: they need at ",
        "least 3 complete pairs and spread on both sides."
      )
    )
  }

  c(header, utils::capture.output(print(summarised)), note)
}

#' @rdname format.bmmtools_cross_check
#' @export
print.bmmtools_cross_check <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname format.bmmtools_cross_check
#' @export
format.bmmtools_cross_check_summary <- function(x, ...) {
  utils::capture.output(print(tibble::as_tibble(x), n = Inf, width = Inf))
}

#' @rdname format.bmmtools_cross_check
#' @export
print.bmmtools_cross_check_summary <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

# the prior check --------------------------------------------------------

#' The prior-check contract
#' @noRd
prior_check_contract <- function() {
  c(
    prior = "character",
    response = "character",
    statistic = "character",
    value = "double"
  )
}

#' Construct a prior-check object
#'
#' The summary rows are the object; the draws that produced them travel
#' as an attribute, because a table of statistics cannot be plotted as a
#' distribution and the distribution is the picture the check exists for.
#'
#' @noRd
new_bmmtools_prior_check <- function(x, ..., error_call = rlang::caller_env()) {
  contract <- prior_check_contract()
  missing <- setdiff(names(contract), names(x))
  if (length(missing) > 0L) {
    cli::cli_abort(
      c(
        "A prior check is missing the column{?s} {.val {missing}}.",
        i = "The contract is {.val {names(contract)}}."
      ),
      call = error_call
    )
  }
  x$value <- as.double(x$value)
  for (column in c("prior", "response", "statistic")) {
    x[[column]] <- as.character(x[[column]])
  }
  x <- tibble::as_tibble(x)[names(contract)]
  structure(x, class = c("bmmtools_prior_check", class(x)), ...)
}

#' @noRd
check_prior_check_contract <- function(x, call = rlang::caller_env()) {
  missing <- setdiff(names(prior_check_contract()), names(x))
  if (length(missing) > 0L) {
    cli::cli_abort(
      "This {.cls bmmtools_prior_check} is missing the column{?s} \\
       {.val {missing}}.",
      call = call
    )
  }
  invisible(x)
}

#' @noRd
#' @importFrom dplyr dplyr_reconstruct
#' @exportS3Method dplyr::dplyr_reconstruct
dplyr_reconstruct.bmmtools_prior_check <- function(data, template) {
  if (!all(names(prior_check_contract()) %in% names(data))) {
    return(tibble::as_tibble(data))
  }
  NextMethod()
}

#' Subset a prior-check object
#'
#' As for [`[.bmmtools_recovery`]: `dplyr::select()` on a tibble subclass
#' subsets through `[` and never reaches `dplyr_reconstruct()`, so
#' without this method a verb that drops a contract column would return a
#' broken object still wearing the class.
#'
#' @param x A `bmmtools_prior_check` object.
#' @param ... Passed to the tibble method.
#' @return A prior-check object while the contract holds, a plain tibble
#'   once it does not.
#' @export
`[.bmmtools_prior_check` <- function(x, ...) {
  demote_if_incomplete(NextMethod(), names(prior_check_contract()))
}

#' Compare prior sets on the observable scale
#'
#' One row per response and statistic, one column per prior set. With
#' exactly two sets a `difference` column is appended, the second minus
#' the first. Read it as an equivalence table: two priors are equivalent
#' on the observable scale where the differences are negligible. With three or
#' more sets there is no unambiguous contrast, so none is added.
#'
#' @param object A `bmmtools_prior_check` from [prior_check()].
#' @param ... Not used.
#'
#' @return A `bmmtools_prior_check_summary` tibble.
#'
#' @examples
#' summary(prior_check_sdt_yn)
#'
#' @export
summary.bmmtools_prior_check <- function(object, ...) {
  check_prior_check_contract(object)
  sets <- attr(object, "sets") %||% unique(object$prior)
  key <- paste(object$response, object$statistic, sep = "\r")
  keys <- unique(key)
  first <- match(keys, key)

  out <- tibble::tibble(
    response = object$response[first],
    statistic = object$statistic[first]
  )
  for (set in sets) {
    in_set <- ifelse(object$prior == set, key, NA_character_)
    out[[set]] <- object$value[match(keys, in_set)]
  }
  if (length(sets) == 2L) {
    out$difference <- out[[sets[[2L]]]] - out[[sets[[1L]]]]
  }

  structure(out, class = c("bmmtools_prior_check_summary", class(out)))
}

#' @rdname summary.bmmtools_prior_check
#' @export
summary.bmmtools_prior_check_summary <- function(object, ...) {
  object
}

#' Format and print a prior check
#'
#' The header names the model and the prior sets compared, because a
#' table of floor rates without them cannot be read. Where a rate is `NA`
#' the reason is printed as a sentence rather than left as a blank cell.
#'
#' @param x A `bmmtools_prior_check` object.
#' @param ... Not used.
#'
#' @return `format()` returns a character vector; `print()` returns `x`
#'   invisibly.
#'
#' @export
format.bmmtools_prior_check <- function(x, ...) {
  check_prior_check_contract(x)
  if (nrow(x) == 0L) {
    return(c("<bmmtools_prior_check>", "Nothing was summarised."))
  }
  model <- attr(x, "model")
  sets <- attr(x, "sets") %||% unique(x$prior)
  n_draws <- attr(x, "n_draws")

  header <- c(
    "<bmmtools_prior_check>",
    paste0(
      "Model: ", model$name %||% class(model)[[2L]] %||% "unknown",
      "; ", length(sets), " prior set", if (length(sets) != 1L) "s", ": ",
      paste(sets, collapse = ", "), "."
    ),
    paste0(
      "Prior-predictive draws per set: ",
      if (is.null(n_draws)) "unknown" else n_draws, "."
    ),
    ""
  )

  rates <- x$value[x$statistic %in% c("floor_rate", "ceiling_rate")]
  note <- character(0)
  if (length(rates) > 0L && anyNA(rates)) {
    note <- c(
      "",
      paste0(
        "The floor and ceiling rates are NA: no response range is known ",
        "for this model. Give `range`, or a `summary` of your own."
      )
    )
  }

  c(header, utils::capture.output(print(summary(x))), note)
}

#' @rdname format.bmmtools_prior_check
#' @export
print.bmmtools_prior_check <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname format.bmmtools_prior_check
#' @export
print.bmmtools_prior_check_summary <- function(x, ...) {
  print(tibble::as_tibble(x), n = Inf, width = Inf)
  invisible(x)
}
