# Parameter recovery: score fits against the values that generated the
# data.
#
# Both entry points work on the estimates tibble, so a fit is one
# accepted input among several rather than a requirement. That is what
# makes the score layer testable and installable with no fitting package
# present.

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
    # a tibble that does not say which estimator produced it is a
    # posterior (milestone 8, decision 35). NA is filled too, not only a
    # missing column: binding a labelled tibble to an unlabelled one
    # leaves NA, and an NA estimator would become its own summary group.
    out$estimator <- fill_estimator(out)
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
  check_truth_frame(truth, call = call)
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
#' unstated scale, scoring falls back to the link scale and says so.
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
#' `cross_check()` joins its reference the same way, so `values`, `arg`
#' and `hint` carry the other side's column names and vocabulary. The
#' defaults are `recover()`'s, so its behaviour is unchanged.
#'
#' @param values The columns of `truth` the join carries over.
#' @param arg The name this side goes by in the messages.
#' @param hint The closing line of the nothing-matched error.
#' @noRd
join_truth <- function(estimates, truth, keys,
                       values = "true_value",
                       arg = "truth",
                       hint = "Generating values use bmm's parameter names.",
                       call = rlang::caller_env()) {
  truth <- tibble::as_tibble(truth)[unique(c(keys, values))]

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
        "No term in {.arg {arg}} matches an estimated parameter.",
        i = "Estimated: {.val {unique(estimates$term)}}.",
        i = "In {.arg {arg}}: {.val {unique(truth$term)}}.",
        i = hint
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
      "{nrow(unmatched)} row{?s} of {.arg {arg}} had no matching \\
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
#' `values` names the point columns to transform; the interval bounds are
#' always transformed. `cross_check()` passes `"estimate"` alone, because
#' its reference arrives on the comparison scale already and must not be
#' transformed a second time.
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
to_natural_scale <- function(x, links, values = c("estimate", "true_value"),
                             call = rlang::caller_env()) {
  unbounded <- character()
  for (term in unique(x$term)) {
    link <- link_of(term, links)
    rows <- x$term == term

    spans_zero <- x$ci_low[rows] < 0 & x$ci_high[rows] > 0
    # a row with no interval --- a failed ML fit keeps its row with
    # estimate and bounds NA (decision 40) --- neither spans zero nor
    # does not. Left NA it would both index NA rows here and make
    # `any()` return NA, which `if ()` refuses.
    spans_zero[is.na(spans_zero)] <- FALSE
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
    for (value in values) {
      x[[value]][rows] <- inverse_link(x[[value]][rows], link)
    }
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

#' Score the estimates of one level against its truth
#'
#' `sd` and `effect` rows are kept on the link scale whatever
#' `resolved$scale` says, because neither the SD of a link-scale random
#' effect nor a contrast between link-scale values has a natural-scale
#' counterpart that one inverse link would give: a difference of two
#' values is not the difference of their inverse links.
#'
#' @noRd
score_level <- function(estimates, truth, level, resolved, error_call) {
  estimates <- estimates[estimates$level == level, , drop = FALSE]
  keys <- level_keys(level)
  if ("replication" %in% names(truth)) keys <- c(keys, "replication")
  if ("condition" %in% names(truth)) keys <- c(keys, "condition")

  joined <- join_truth(estimates, truth, keys, call = error_call)

  scale <- if (level %in% c("sd", "effect")) "link" else resolved$scale
  if (identical(scale, "natural")) {
    joined <- to_natural_scale(joined, resolved$links)
  }

  joined$bias <- joined$estimate - joined$true_value
  joined$covered <- joined$true_value >= joined$ci_low &
    joined$true_value <= joined$ci_high
  joined$scale <- rep(scale, nrow(joined))
  joined
}

#' The columns a level's truth joins on, besides replication and condition
#' @noRd
level_keys <- function(level) {
  if (identical(level, "subject")) c("term", "id") else "term"
}

#' Match `truth` to the levels requested
#'
#' One level takes a data frame. Several take a named list with an
#' element per level, which is the shape of a `bmmtools_simulation`'s
#' `truth`, so it can be passed as it is.
#'
#' @return A named list of truth data frames, one per level.
#' @noRd
truth_by_level <- function(truth, level, call = rlang::caller_env()) {
  if (length(level) == 1L) {
    check_truth_frame(truth, call = call)
    return(stats::setNames(list(truth), level))
  }
  if (!is.list(truth) || is.data.frame(truth) || is.null(names(truth))) {
    cli::cli_abort(
      c(
        "With several levels, {.arg truth} must be a named list with an \\
         element for each level, not {.obj_type_friendly {truth}}.",
        i = "The {.field truth} of a {.cls bmmtools_simulation} has that \\
             shape."
      ),
      call = call
    )
  }
  missing <- setdiff(level, names(truth))
  if (length(missing) > 0L) {
    cli::cli_abort(
      "{.arg truth} has no element for the level{?s} {.val {missing}}.",
      call = call
    )
  }
  truth[level]
}

#' @noRd
check_truth_frame <- function(truth, call = rlang::caller_env()) {
  if (!is.data.frame(truth)) {
    cli::cli_abort(
      "{.arg truth} must be a data frame, not {.obj_type_friendly {truth}}.",
      call = call
    )
  }
  invisible(truth)
}

#' The shared body of recover() and recover_subjects()
#' @noRd
score_recovery <- function(fits, truth, level, group, scale, links,
                           ci_level, drop_constants, call, error_call) {
  estimates <- as_estimates_input(
    fits, level, group, ci_level, drop_constants,
    call = error_call
  )
  truths <- truth_by_level(truth, level, call = error_call)
  # every input is checked before resolve_links() can print its message
  for (lv in level) {
    if (!any(estimates$level == lv)) {
      cli::cli_abort(
        "No estimates at level {.val {lv}} to score.",
        call = error_call
      )
    }
    check_truth(truths[[lv]], level_keys(lv), call = error_call)
  }
  resolved <- resolve_links(fits, links, scale, call = error_call)

  pieces <- lapply(level, function(lv) {
    score_level(estimates, truths[[lv]], lv, resolved, error_call)
  })
  link_only <- intersect(c("sd", "effect"), level)
  if (length(link_only) > 0L && identical(resolved$scale, "natural")) {
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    what <- unname(c(
      sd = "Standard deviations", effect = "Effects"
    )[link_only])
    cli::cli_inform(c(
      "{.or {what}} are scored on the link scale.",
      i = "Their {.field scale} column reads {.val link}."
    ))
  }
  joined <- dplyr::bind_rows(pieces)
  check_estimator_balance(joined)

  new_bmmtools_recovery(
    joined,
    scale = resolved$scale,
    ci_level = joined$ci_level[[1L]],
    call = call,
    error_call = error_call
  )
}

#' Warn when two estimators were not scored on the same subjects
#'
#' Comparing two estimators is only meaningful when both were scored on
#' the same people: `metric_bias()` and `metric_rmse()` drop incomplete
#' pairs in silence, so an estimator that failed on the hard subjects
#' would otherwise report a flattering bias computed on the easy ones and
#' nothing would say so (milestone 8, decision 40).
#'
#' What counts is the pair a metric can use, not the row. Decision 40
#' keeps a failed ML subject *as a row* with `estimate = NA`, so comparing
#' which ids are present finds nothing --- the id is always there and the
#' estimate is what went missing.
#'
#' The comparison is also keyed by the cell and the parameter, not pooled.
#' Subject ids repeat in every cell of a grid, so a whole cell's failed ML
#' fit would intersect away against the same ids in another cell.
#'
#' @noRd
check_estimator_balance <- function(x) {
  subjects <- x[x$level == "subject" & !is.na(x$id), , drop = FALSE]
  estimators <- unique(subjects$estimator)
  if (length(estimators) < 2L) {
    return(invisible(x))
  }
  usable <- subjects[
    !is.na(subjects$estimate) & !is.na(subjects$true_value), ,
    drop = FALSE
  ]
  # `condition` is filled by fill_optional_columns() inside the
  # constructor, which has not run yet, so the cell keys are read
  # defensively and always contribute a field
  part <- function(nm) {
    if (nm %in% names(usable)) {
      as.character(usable[[nm]])
    } else {
      rep("", nrow(usable))
    }
  }
  key <- paste(
    part("condition"), part("replication"), usable$term, usable$id,
    sep = "\r"
  )
  keys <- lapply(estimators, function(e) unique(key[usable$estimator == e]))
  shared <- Reduce(intersect, keys)
  unbalanced <- setdiff(unique(unlist(keys)), shared)
  if (length(unbalanced) > 0L) {
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    ids <- unique(vapply(strsplit(unbalanced, "\r", fixed = TRUE), function(p) {
      p[[4L]]
    }, character(1)))
    cells <- length(unique(vapply(
      strsplit(unbalanced, "\r", fixed = TRUE),
      function(p) paste(p[1:2], collapse = "\r"), character(1)
    )))
    cli::cli_warn(c(
      "The estimators were not scored on the same subjects.",
      i = "Missing from at least one estimator: {.val {ids}}.",
      if (cells > 1L) {
        c(i = "In {cells} cells of the grid.")
      },
      i = "Metrics drop incomplete pairs, so the rows are not comparable \\
           until you filter to the subjects every estimator has."
    ))
  }
  invisible(x)
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
#'   same generating values are scored against every replication. When
#'   `recover()` scores several levels, a named list with a data frame for
#'   each, such as the `truth` of a `bmmtools_simulation`.
#' @param level The levels `recover()` scores: `"population"` (the
#'   default), `"sd"` (the between-subject standard deviations),
#'   `"effect"` (the contrasts of a `coding = "contrast"` design), or
#'   several of them. `"sd"` and `"effect"` are always scored on the link
#'   scale, whatever `scale` says.
#' @param group The grouping factor subject-level estimates come from.
#' @param scale `"natural"` scores on the scale a reader interprets,
#'   `"link"` on the scale the model was estimated on. The choice
#'   matters: coverage is invariant to a monotone link but bias, RMSE and
#'   the correlations are not.
#' @param links A named character vector mapping a term to one of the
#'   link names [inverse_link()] understands. `NULL` reads the link table
#'   from a `bmmfit`; with no table available, scoring falls back to the
#'   link scale and says so. A term with no entry takes the link of the
#'   longest entry it starts with followed by `_` (`kappa_task1` takes
#'   the link of `kappa`), and is treated as `"identity"` when there is
#'   none.
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
#'   `covered`, `scale`, `level`, `id`, `converged`, `condition`,
#'   `estimator` and `replication`.
#'   `converged` is the verdict of [check_convergence()] with its
#'   default thresholds when `fits` are fit objects; to gate with other
#'   thresholds, call [extract_estimates()] with `converged =` first and
#'   pass the tibble. Call [summary()] on it for the per-parameter
#'   metrics.
#'
#'   `estimator` names how each row was produced and defaults to
#'   `"bayes"`. [summary()] groups by it, so two estimators of the same
#'   parameter --- a hierarchical posterior and a subject-wise
#'   maximum-likelihood fit, say --- can be bound together and scored
#'   against one truth without being pooled into a single bias and RMSE.
#'   Set it with `extract_estimates(estimator = )` or as a column on a
#'   hand-built tibble, then `dplyr::bind_rows()` the two and pass the
#'   result here. Note that this `estimator` is unrelated to the one in
#'   [extract_correlations()], which names how a *correlation* was read.
#'
#' @details
#' Standard deviations are always scored on the link scale, the scale
#' the model estimates them on. Under `scale = "natural"` a message says
#' so, their `scale` column reads `"link"`, and the other rows and the
#' object's `scale` attribute stay on the natural scale.
#'
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
                    level = "population",
                    scale = c("natural", "link"),
                    links = NULL,
                    ci_level = 0.95,
                    drop_constants = TRUE,
                    ...) {
  rlang::check_dots_empty()
  level <- rlang::arg_match(
    level, c("population", "effect", "sd"),
    multiple = TRUE
  )
  scale <- rlang::arg_match(scale)
  score_recovery(
    fits, truth,
    level = level, group = NULL, scale = scale, links = links,
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
