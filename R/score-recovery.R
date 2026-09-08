# Parameter recovery: score fits against the values that generated the
# data (spec sections 4 and 5).
#
# Both entry points work on the estimates tibble, so a fit is one
# accepted input among several rather than a requirement. That is what
# makes the score layer testable and installable with no fitting package
# present (ARCHITECTURE.md decisions 1 and 6).

#' Turn whatever was passed as `fits` into an estimates tibble
#'
#' A fit, a list of fits (one per replication), or an estimates tibble
#' that has been extracted already.
#'
#' @noRd
as_estimates_input <- function(fits, level, group, ci_level, drop_constants,
                               call = rlang::caller_env()) {
  if (inherits(fits, "brmsfit")) {
    out <- extract_estimates(
      fits,
      level = level, group = group,
      ci_level = ci_level, drop_constants = drop_constants
    )
    out$replication <- 1L
    return(out)
  }

  is_fit_list <- is.list(fits) && !is.data.frame(fits) &&
    length(fits) > 0L &&
    all(vapply(fits, inherits, logical(1), what = "brmsfit"))
  if (is_fit_list) {
    labels <- names(fits)
    if (is.null(labels)) labels <- seq_along(fits)
    pieces <- lapply(seq_along(fits), function(i) {
      out <- extract_estimates(
        fits[[i]],
        level = level, group = group,
        ci_level = ci_level, drop_constants = drop_constants
      )
      out$replication <- labels[[i]]
      out
    })
    return(dplyr::bind_rows(pieces))
  }

  if (is.data.frame(fits)) {
    missing <- setdiff(estimates_required_columns(), names(fits))
    if (length(missing) > 0L) {
      cli::cli_abort(
        c(
          "{.arg fits} is missing the estimates column{?s} {.val {missing}}.",
          i = "An estimates tibble comes from {.fn extract_estimates}."
        ),
        call = call
      )
    }
    out <- tibble::as_tibble(fits)
    if (!"replication" %in% names(out)) out$replication <- 1L
    # a hand-built tibble carries no verdict; an unknown is not a failure
    if (!"converged" %in% names(out)) out$converged <- NA
    return(out)
  }

  cli::cli_abort(
    c(
      "{.arg fits} must be a {.cls brmsfit}, a list of them, or an \\
       estimates tibble, not {.obj_type_friendly {fits}}.",
      i = "An estimates tibble comes from {.fn extract_estimates}."
    ),
    call = call
  )
}

#' Check that the truth tibble carries what the join needs
#' @noRd
check_truth <- function(truth, keys, call = rlang::caller_env()) {
  if (!is.data.frame(truth)) {
    cli::cli_abort(
      "{.arg truth} must be a data frame, not {.obj_type_friendly {truth}}.",
      call = call
    )
  }
  required <- c(keys, "true_value")
  missing <- setdiff(required, names(truth))
  if (length(missing) > 0L) {
    cli::cli_abort(
      "{.arg truth} is missing the column{?s} {.val {missing}}.",
      call = call
    )
  }
  invisible(truth)
}

#' Decide which link belongs to each term, and on which scale to score
#'
#' A `bmmfit` carries its own link table, so scoring on the natural scale
#' needs no argument. A bare `brmsfit` or a bare estimates tibble carries
#' none, and rather than guess an identity link and report numbers on an
#' unstated scale, scoring falls back to the link scale and says so
#' (ARCHITECTURE.md decision 2).
#'
#' @noRd
resolve_links <- function(fits, links, scale, call = rlang::caller_env()) {
  if (identical(scale, "link")) {
    return(list(links = NULL, scale = "link"))
  }

  if (!is.null(links)) {
    if (is.list(links)) links <- unlist(links)
    named_character <- is.character(links) &&
      !is.null(names(links)) &&
      !anyNA(names(links)) &&
      all(nzchar(names_or_empty(links)))
    if (!named_character) {
      cli::cli_abort(
        "{.arg links} must be a named character vector mapping a term \\
         to a link name.",
        call = call
      )
    }
    return(list(links = links, scale = "natural"))
  }

  model_links <- model_links_of(fits)
  if (!is.null(model_links)) {
    return(list(links = model_links, scale = "natural"))
  }

  cli::cli_inform(c(
    "No link information available, so scoring on the link scale.",
    i = "Supply {.arg links}, or pass a {.cls bmmfit}, to score on the \\
         natural scale."
  ))
  list(links = NULL, scale = "link")
}

#' @noRd
names_or_empty <- function(x) {
  nms <- names(x)
  if (is.null(nms)) character(0) else nms
}

#' The link table a bmm fit carries, if there is one
#' @noRd
model_links_of <- function(fits) {
  fit <- fits
  is_fit_list <- is.list(fits) && !is.data.frame(fits) &&
    length(fits) > 0L && inherits(fits[[1L]], "brmsfit")
  if (is_fit_list) {
    fit <- fits[[1L]]
  }
  if (!inherits(fit, "brmsfit")) {
    return(NULL)
  }
  links <- fit$bmm$model$links
  if (is.null(links) || length(links) == 0L) {
    return(NULL)
  }
  unlist(links)
}

#' Join estimates to truth, warning about what did not match
#'
#' A generating value with no estimate is a warning that lists what was
#' available, ported from `check_recovery()`; an estimate with no
#' generating value is dropped in silence, because fits routinely
#' estimate more than a grid varies. Every generating value unmatched is
#' an error, because it almost always means a naming mismatch and a
#' zero-row result would hide it until a plot came out blank.
#'
#' @noRd
join_truth <- function(estimates, truth, keys, call = rlang::caller_env()) {
  truth <- tibble::as_tibble(truth)[unique(c(keys, "true_value"))]

  unmatched <- dplyr::anti_join(
    dplyr::distinct(truth[keys]), estimates[keys],
    by = keys
  )
  joined <- dplyr::inner_join(
    estimates, truth,
    by = keys, relationship = "many-to-one"
  )

  if (nrow(joined) == 0L) {
    cli::cli_abort(
      c(
        "No term in {.arg truth} matches an estimated parameter.",
        i = "Estimated: {.val {unique(estimates$term)}}.",
        i = "In {.arg truth}: {.val {unique(truth$term)}}.",
        i = "Generating values use bmm's own parameter names."
      ),
      call = call
    )
  }

  if (nrow(unmatched) > 0L) {
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    labels <- if ("id" %in% keys) {
      paste0(unmatched$term, " (id ", unmatched$id, ")")
    } else {
      unmatched$term
    }
    cli::cli_warn(c(
      "{nrow(unmatched)} row{?s} of {.arg truth} had no matching \\
       estimate and {?was/were} dropped.",
      x = "Not found: {.val {labels}}.",
      i = "Estimated: {.val {unique(estimates$term)}}."
    ))
  }

  joined
}

#' Put estimates, both interval bounds and truth on the natural scale
#'
#' All four quantities go through the same inverse link, so that both
#' sides of the comparison are on one scale and `bias` and `covered` are
#' computed after the transform rather than before it.
#'
#' The bounds are reordered afterwards. Two links in the vocabulary are
#' decreasing (`inverse` and `loglog`), so the transform of the lower
#' bound is the upper bound of the transformed interval; leaving them
#' swapped would make every interval empty and coverage zero.
#'
#' Two links are not monotone across zero, and an interval that spans
#' zero on the link scale needs more than a reorder. Under `sqrt` the
#' natural value is `x^2`, whose minimum, 0, lies inside such an
#' interval, so the image is `[0, max(low^2, high^2)]`. Under `inverse`
#' the natural value `1 / x` diverges inside it, so the image is not an
#' interval at all: the bounds become `NA`, coverage becomes `NA`, and a
#' warning names the terms, because a plausible-looking wrong interval
#' would silently deflate coverage.
#'
#' @noRd
to_natural_scale <- function(x, links, call = rlang::caller_env()) {
  unbounded <- character()
  for (term in unique(x$term)) {
    link <- "identity"
    if (!is.null(links) && term %in% names(links)) link <- links[[term]]
    rows <- x$term == term

    spans_zero <- x$ci_low[rows] < 0 & x$ci_high[rows] > 0
    low <- inverse_link(x$ci_low[rows], link)
    high <- inverse_link(x$ci_high[rows], link)
    if (identical(link, "sqrt") && any(spans_zero)) {
      high[spans_zero] <- pmax(low, high)[spans_zero]
      low[spans_zero] <- 0
    }
    if (identical(link, "inverse") && any(spans_zero)) {
      low[spans_zero] <- NA_real_
      high[spans_zero] <- NA_real_
      unbounded <- c(unbounded, term)
    }
    x$estimate[rows] <- inverse_link(x$estimate[rows], link)
    x$true_value[rows] <- inverse_link(x$true_value[rows], link)
    x$ci_low[rows] <- pmin(low, high)
    x$ci_high[rows] <- pmax(low, high)
  }
  if (length(unbounded) > 0L) {
    cli::cli_warn(
      c(
        "Interval{?s} for {.val {unbounded}} span{?s/} zero on the link \\
         scale under the {.val inverse} link, so the natural-scale image \\
         is not an interval.",
        i = "Their bounds and {.field covered} are {.code NA}; score these \\
             terms with {.code scale = \"link\"}."
      ),
      call = call
    )
  }
  x
}

#' The shared body of recover() and recover_subjects()
#' @noRd
score_recovery <- function(fits, truth, level, group, scale, links,
                           ci_level, drop_constants, call, error_call) {
  estimates <- as_estimates_input(
    fits, level, group, ci_level, drop_constants,
    call = error_call
  )
  estimates <- estimates[estimates$level == level, , drop = FALSE]
  if (nrow(estimates) == 0L) {
    cli::cli_abort(
      "No estimates at level {.val {level}} to score.",
      call = error_call
    )
  }

  keys <- if (identical(level, "subject")) c("term", "id") else "term"
  check_truth(truth, keys, call = error_call)
  if ("replication" %in% names(truth)) keys <- c(keys, "replication")

  resolved <- resolve_links(fits, links, scale, call = error_call)
  joined <- join_truth(estimates, truth, keys, call = error_call)

  if (identical(resolved$scale, "natural")) {
    joined <- to_natural_scale(joined, resolved$links)
  }

  joined$bias <- joined$estimate - joined$true_value
  joined$covered <- joined$true_value >= joined$ci_low &
    joined$true_value <= joined$ci_high
  joined$scale <- resolved$scale

  new_bmmtools_recovery(
    joined,
    scale = resolved$scale,
    ci_level = joined$ci_level[[1L]],
    call = call,
    error_call = error_call
  )
}

#' Score parameter recovery against known generating values
#'
#' `recover()` answers the question a model developer has to answer
#' before anything else: when the data are generated from known
#' parameters, does fitting the model give those parameters back? It
#' takes fits and truth and returns one row per fit and parameter, with
#' the estimate, its interval, the generating value, the signed error and
#' whether the interval covered.
#'
#' `recover_subjects()` is the person-level variant. It asks whether the
#' model orders *people* correctly, which is the question a reliability
#' analysis needs and which population-level recovery cannot answer.
#'
#' @param fits A `brmsfit` (so also a `bmmfit`), a list of them with one
#'   element per replication, or an estimates tibble from
#'   [extract_estimates()]. The tibble is accepted so that a scoring
#'   pipeline runs with no fitting package installed.
#' @param truth A data frame of generating values on the **link scale**,
#'   under bmm's own parameter names. `recover()` needs `term` and
#'   `true_value`; `recover_subjects()` also needs `id`. A `replication`
#'   column is used as a join key when present, and when it is absent the
#'   same generating values are scored against every replication.
#' @param group The grouping factor subject-level estimates come from.
#' @param scale `"natural"` scores on the scale a reader interprets,
#'   `"link"` on the scale the model was estimated on. The choice
#'   matters: coverage is invariant to a monotone link but bias, RMSE and
#'   the correlations are not.
#' @param links A named character vector mapping a term to one of the
#'   link names [inverse_link()] understands. `NULL` reads the link table
#'   from a `bmmfit`; with no table available, scoring falls back to the
#'   link scale and says so. A term with no entry is treated as
#'   `"identity"`.
#' @param ci_level The interval mass, passed to [extract_estimates()]
#'   when `fits` is a fit.
#' @param drop_constants Passed to [extract_estimates()]. Parameters the
#'   model fixed are dropped by default: scoring them would report a
#'   parameter as perfectly recovered that was never estimated.
#' @param ... Not used.
#'
#' @return A `bmmtools_recovery` object: a tibble subclass with the
#'   columns `term`, `estimate`, `ci_low`, `ci_high`, `ci_method`,
#'   `ci_level`, `rhat`, `ess_bulk`, `ess_tail`, `true_value`, `bias`,
#'   `covered`, `scale`, `level`, `id`, `converged` and `replication`.
#'   `converged` is the verdict of [check_convergence()] with its
#'   default thresholds when `fits` are fit objects; to gate with other
#'   thresholds, call [extract_estimates()] with `converged =` first and
#'   pass the tibble. Call [summary()] on it for the per-parameter
#'   metrics.
#'
#' @details
#' A term in `truth` that no fit estimated produces a warning listing
#' what was available, and is dropped. A term a fit estimated that
#' `truth` does not mention is dropped silently, because fits routinely
#' estimate more than a simulation grid varies. Every term unmatched is
#' an error rather than an empty result.
#'
#' @examples
#' estimates <- tibble::tibble(
#'   term = c("kappa", "thetat"),
#'   estimate = c(2.1, 1.0), ci_low = c(1.6, 0.4), ci_high = c(2.6, 1.6),
#'   ci_method = "eti", ci_level = 0.95,
#'   rhat = 1, ess_bulk = 900, ess_tail = 900,
#'   level = "population", id = NA_character_
#' )
#' truth <- tibble::tibble(term = c("kappa", "thetat"), true_value = c(2, 1.1))
#'
#' recovery <- recover(
#'   estimates, truth,
#'   links = c(kappa = "log", thetat = "logit")
#' )
#' summary(recovery)
#'
#' @export
recover <- function(fits,
                    truth,
                    scale = c("natural", "link"),
                    links = NULL,
                    ci_level = 0.95,
                    drop_constants = TRUE,
                    ...) {
  rlang::check_dots_empty()
  scale <- rlang::arg_match(scale)
  score_recovery(
    fits, truth,
    level = "population", group = NULL, scale = scale, links = links,
    ci_level = ci_level, drop_constants = drop_constants,
    call = match.call(), error_call = rlang::current_env()
  )
}

#' @rdname recover
#' @export
recover_subjects <- function(fits,
                             truth,
                             group = "id",
                             scale = c("natural", "link"),
                             links = NULL,
                             ci_level = 0.95,
                             drop_constants = TRUE,
                             ...) {
  rlang::check_dots_empty()
  scale <- rlang::arg_match(scale)
  score_recovery(
    fits, truth,
    level = "subject", group = group, scale = scale, links = links,
    ci_level = ci_level, drop_constants = drop_constants,
    call = match.call(), error_call = rlang::current_env()
  )
}
