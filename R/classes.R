# The recovery classes (spec section 6).
#
# `bmmtools_recovery` is a tibble subclass: one row per fit x parameter,
# carrying the apabayes `parameters` columns first and bmmtools's
# additions after (ARCHITECTURE.md decision 14). Subclassing a tibble
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
    converged = "logical"
  )
}

#' @noRd
recovery_contract_columns <- function() {
  c(names(recovery_contract()), "replication")
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
    "term", "level", "scale", "n", "n_replications", "n_converged",
    "bias", "rmse", "coverage", "ci_width",
    "r", "r_low", "r_high", "rank_r",
    "ccc", "ccc_scale_bias", "ccc_location_bias"
  )
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
  x <- tibble::as_tibble(x)[recovery_summary_columns()]
  structure(x, class = c("bmmtools_recovery_summary", class(x)))
}

# dplyr ------------------------------------------------------------------

#' Demote to a plain tibble when a contract column has gone
#'
#' The one place that decides whether a result is still a recovery
#' object. A print, summary or plot method that assumed a missing column
#' would fail later and further from the cause.
#'
#' @noRd
demote_if_incomplete <- function(x) {
  if (!all(recovery_contract_columns() %in% names(x))) {
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
  demote_if_incomplete(NextMethod())
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
      ccc_scale_bias = ccc$scale_bias,
      ccc_location_bias = ccc$location_bias
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
      ccc_scale_bias = ccc$scale_bias,
      ccc_location_bias = ccc$location_bias
    )
  })
  pull <- function(name) vapply(per, function(p) as.double(p[[name]]), 0)

  n <- pull("n")
  combined <- fisher_z_combine(pull("r"), n)
  ranks <- fisher_z_combine(pull("rank_r"), n)

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
      # concordance is not a correlation on Fisher's z scale, so the
      # within-replication values are averaged as they are
      ccc = mean_or_na(pull("ccc")),
      ccc_scale_bias = mean_or_na(pull("ccc_scale_bias")),
      ccc_location_bias = mean_or_na(pull("ccc_location_bias"))
    )
  )
}

#' Summarise a recovery object into per-parameter metrics
#'
#' One row per parameter and level, with the metrics of
#' `ARCHITECTURE.md` decision 18: bias, RMSE, coverage, mean interval
#' width, the Pearson correlation with a Fisher-z interval, the Spearman
#' correlation and Lin's concordance with its two components.
#'
#' Correlation metrics are `NA`, never `0`, when fewer than three
#' complete pairs are available or when either side has no spread: `0`
#' would read as "no recovery" where the honest answer is "not estimable
#' from this design". The correlation interval is a 95% confidence
#' interval and does not follow `ci_level`, which is the mass of the
#' posterior interval and a different quantity.
#'
#' Subject-level objects are summarised **within replication and then
#' combined** on Fisher's z scale. Pooling subjects across replications
#' would mix between-subject with between-replication variance.
#'
#' @param object A `bmmtools_recovery` object from [recover()] or
#'   [recover_subjects()].
#' @param ... Not used.
#'
#' @return A `bmmtools_recovery_summary` tibble with the columns `term`,
#'   `level`, `scale`, `n`, `n_replications`, `n_converged`, `bias`,
#'   `rmse`, `coverage`, `ci_width`, `r`, `r_low`, `r_high`, `rank_r`,
#'   `ccc`, `ccc_scale_bias` and `ccc_location_bias`.
#'
#' @details
#' `n_converged` is the number of replications whose fit passed
#' [check_convergence()], read from the `converged` column that
#' [extract_estimates()] fills. It is `NA` when no fit carried a
#' verdict, as with a hand-built estimates tibble: reporting it as equal
#' to `n` there would assert something that was never measured.
#'
#' @export
summary.bmmtools_recovery <- function(object, ...) {
  check_recovery_contract(object)
  if (nrow(object) == 0L) {
    return(new_bmmtools_recovery_summary(
      empty_recovery_summary()
    ))
  }

  keys <- paste(object$level, object$term, sep = "\r")
  pieces <- lapply(split(seq_len(nrow(object)), keys), function(i) {
    rows <- object[i, ]
    body <- if (identical(rows$level[[1L]], "subject")) {
      summarise_subject(rows)
    } else {
      summarise_population(rows)
    }
    tibble::as_tibble(c(
      list(
        term = rows$term[[1L]],
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
    term = "character", level = "character", scale = "character",
    n = "double", n_replications = "integer", n_converged = "integer",
    bias = "double", rmse = "double", coverage = "double",
    ci_width = "double", r = "double", r_low = "double",
    r_high = "double", rank_r = "double", ccc = "double",
    ccc_scale_bias = "double", ccc_location_bias = "double"
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
  summarised <- summary(x)

  header <- c(
    "<bmmtools_recovery>",
    paste0("Scored on the ", scale, " scale."),
    paste0(
      n_fits, " fit", if (n_fits != 1L) "s", ", ",
      length(terms), " parameter", if (length(terms) != 1L) "s",
      ": ", paste(terms, collapse = ", "), "."
    ),
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
