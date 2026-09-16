# The reference cross-check.
#
# recover() asks whether a model gives back the values that generated the
# data. cross_check() asks the other question a new model has to answer:
# does it agree with what is already known? The right-hand side is a
# closed form, another implementation or a published value --- a
# comparison, not a truth, which is why `?cross_check` says so where the
# column is called `bias`.
#
# Everything the comparison needs already exists: extract_estimates() for
# both sides, resolve_links() and to_natural_scale() for the scale, and
# join_truth() for the warn-and-list behaviour when the two sides name
# different terms.

#' Refuse a list of fits, and anything without an extraction method
#'
#' A list against one fixed reference is `recover()`'s shape with the
#' reference as `true_value`, and the case that would justify a branch
#' here --- a reference that varies by replication --- a list does not
#' express.
#'
#' @noRd
check_cross_check_fit <- function(fit, call = rlang::caller_env()) {
  is_fit_list <- is.list(fit) && !is.object(fit) && length(fit) > 0L &&
    all(vapply(fit, is_one_fit, logical(1)))
  if (is_fit_list) {
    cli::cli_abort(
      c(
        "{.arg fit} must be one fit, not a list of {length(fit)}.",
        i = "{.fn recover} scores a list of fits against one truth and \\
             separates them with a {.field replication} column; \\
             {.fn cross_check} compares one fit with one reference."
      ),
      call = call
    )
  }
  if (!is_one_fit(fit)) {
    cli::cli_abort(
      "{.arg fit} must be a fit with an {.fn extract_estimates} method, \\
       not {.obj_type_friendly {fit}}.",
      call = call
    )
  }
  invisible(fit)
}

#' Which level the comparison is made at
#'
#' `extract_estimates()` has had `"sd"` and `"cor"` since Milestone 5.2,
#' but a correlation reference would need `var1`, `var2` and a
#' pair-ordering rule of its own, which is a second reference contract
#' for a rare case. Naming them explicitly beats the
#' bare "must be one of" that `arg_match()` would give.
#'
#' @noRd
check_cross_check_level <- function(level, call = rlang::caller_env()) {
  if (length(level) > 1L) {
    if (!identical(level, c("population", "subject"))) {
      cli::cli_abort(
        "{.arg level} must be a single level, {.val population} or \\
         {.val subject}.",
        call = call
      )
    }
    return("population")
  }
  deferred <- is.character(level) && length(level) == 1L &&
    !is.na(level) && level %in% c("sd", "cor")
  if (deferred) {
    cli::cli_abort(
      c(
        "{.arg level} must be {.val population} or {.val subject}, \\
         not {.val {level}}.",
        i = "{.fn extract_estimates} also returns {.val sd} and \\
             {.val cor}, but a reference for them would need \\
             {.field var1}, {.field var2} and a pair-ordering rule of \\
             its own."
      ),
      call = call
    )
  }
  rlang::arg_match0(
    level, c("population", "subject"), "level",
    error_call = call
  )
}

#' The estimates of one side, at one level, on the comparison scale
#'
#' The comparison **scale** is settled once, from `fit`, so the two
#' columns are always on one scale. The **links** are each side's own:
#' see `reference_links()`.
#'
#' @noRd
cross_check_side <- function(x, arg, level, group, ci_level, resolved, dots,
                             call = rlang::caller_env()) {
  estimates <- rlang::exec(
    extract_estimates, x,
    level = level, group = group, ci_level = ci_level, !!!dots
  )
  if (is.data.frame(estimates) && "level" %in% names(estimates)) {
    estimates <- estimates[estimates$level == level, , drop = FALSE]
  } else {
    estimates <- estimates[integer(0), , drop = FALSE]
  }
  if (nrow(estimates) == 0L) {
    cli::cli_abort(
      "{.arg {arg}} has no estimates at level {.val {level}} to compare.",
      call = call
    )
  }
  if (identical(resolved$scale, "natural")) {
    estimates <- to_natural_scale(
      estimates, resolved$links,
      values = "estimate", call = call
    )
  }
  estimates
}

#' Which links invert a fit reference's own draws
#'
#' A fit reference is another implementation, and another implementation
#' may hold the same parameter on a different link: `kappa` on a log
#' scale in one model and on the identity in another. Each side is
#' therefore inverted with **its own** table, and only then are the two
#' natural-scale numbers comparable. Using the left fit's table on the
#' reference's draws would put a silently mis-scaled number in the
#' `reference` column, which is the failure this function exists to
#' prevent.
#'
#' `links` given by hand still overrides both sides, as in `recover()`,
#' and a reference carrying no table of its own (a bare `brmsfit`) falls
#' back to the fit's. The comparison **scale** is not renegotiated here:
#' it stays whatever `fit` settled, so the two columns never end up on
#' different scales.
#'
#' @noRd
reference_links <- function(reference, resolved, links) {
  if (!identical(resolved$scale, "natural") || !is.null(links)) {
    return(resolved)
  }
  own <- model_links_of(reference)
  if (is.null(own)) {
    return(resolved)
  }
  resolved$links <- own
  resolved
}

#' A fit as the right-hand side of the comparison
#'
#' Extracted at the same level, put on the comparison scale through its
#' own links, then renamed to the reference columns. Its `source` is
#' `"fit"`, so a table that mixes a closed form with another
#' implementation says which row came from where.
#'
#' @noRd
reference_from_fit <- function(reference, level, group, ci_level, resolved,
                               links, dots, call = rlang::caller_env()) {
  estimates <- cross_check_side(
    reference, "reference", level, group, ci_level,
    reference_links(reference, resolved, links), dots,
    call = call
  )
  tibble::tibble(
    term = as.character(estimates$term),
    reference = as.double(estimates$estimate),
    ref_low = as.double(estimates$ci_low),
    ref_high = as.double(estimates$ci_high),
    source = "fit",
    id = as.character(estimates$id)
  )
}

#' The reference's interval bounds: both or neither, and both numeric
#'
#' A bound that is not numeric would be coerced to `NA` with nothing but
#' base R's "NAs introduced by coercion" to show for it, and an `NA`
#' bound is indistinguishable from the documented, valid state of a
#' reference with no interval at all. A mistyped column has to say so.
#'
#' @return `TRUE` when the reference carries an interval.
#' @noRd
check_reference_bounds <- function(reference, call = rlang::caller_env()) {
  bounds <- c("ci_low", "ci_high") %in% names(reference)
  if (xor(bounds[[1L]], bounds[[2L]])) {
    cli::cli_abort(
      c(
        "{.arg reference} has one interval bound without the other.",
        i = "Give both {.field ci_low} and {.field ci_high}, or neither."
      ),
      call = call
    )
  }
  if (!all(bounds)) {
    return(FALSE)
  }
  for (bound in c("ci_low", "ci_high")) {
    if (!is.numeric(reference[[bound]])) {
      cli::cli_abort(
        "{.field {bound}} of {.arg reference} must be numeric, \\
         not {.obj_type_friendly {reference[[bound]]}}.",
        call = call
      )
    }
  }
  TRUE
}

#' The reference's subject ids, required and character at that level
#' @noRd
check_reference_id <- function(reference, level, call = rlang::caller_env()) {
  id <- reference[["id"]]
  if (!identical(level, "subject")) {
    return(id)
  }
  if (is.null(id)) {
    cli::cli_abort(
      c(
        "{.arg reference} needs an {.field id} column at the \\
         {.val subject} level.",
        i = "Subject estimates are compared per subject, so both sides \\
             name the subject."
      ),
      call = call
    )
  }
  if (!is.character(id)) {
    cli::cli_abort(
      "{.field id} of {.arg reference} must be a character vector, \\
       not {.obj_type_friendly {id}}.",
      call = call
    )
  }
  id
}

#' A data frame as the right-hand side of the comparison
#'
#' `term` and `estimate` are required; the two interval bounds are both
#' present or both absent, because one of them alone cannot say whether
#' the intervals overlap. Values are read on the comparison `scale`: a
#' closed form and a published number are natural-scale quantities, and a
#' user with a link-scale reference says `scale = "link"`.
#'
#' @noRd
reference_from_tibble <- function(reference, level,
                                  call = rlang::caller_env()) {
  missing <- setdiff(c("term", "estimate"), names(reference))
  if (length(missing) > 0L) {
    cli::cli_abort(
      c(
        "{.arg reference} is missing the column{?s} {.val {missing}}.",
        i = "A reference names a {.field term} and gives an \\
             {.field estimate}, with {.field ci_low} and {.field ci_high} \\
             when it has an interval."
      ),
      call = call
    )
  }
  if (!is.character(reference$term)) {
    cli::cli_abort(
      "{.field term} of {.arg reference} must be a character vector, \\
       not {.obj_type_friendly {reference$term}}.",
      call = call
    )
  }
  if (!is.numeric(reference$estimate)) {
    cli::cli_abort(
      "{.field estimate} of {.arg reference} must be numeric, \\
       not {.obj_type_friendly {reference$estimate}}.",
      call = call
    )
  }

  bounds <- check_reference_bounds(reference, call = call)
  id <- check_reference_id(reference, level, call = call)

  # every column is built before tibble() is called: its data mask would
  # otherwise let the new `reference` column shadow the argument of the
  # same name, and `reference$ci_low` would look inside a double vector
  n <- nrow(reference)
  na_double <- rep(NA_real_, n)
  out <- list(
    term = reference$term,
    reference = as.double(reference$estimate),
    ref_low = if (bounds) as.double(reference$ci_low) else na_double,
    ref_high = if (bounds) as.double(reference$ci_high) else na_double,
    source = rep_len(
      as.character(reference[["source"]] %||% "reference"), n
    ),
    id = if (identical(level, "subject")) {
      as.character(id)
    } else {
      rep(NA_character_, n)
    }
  )
  tibble::as_tibble(out)
}

#' Normalise whatever was passed as `reference` into the reference columns
#' @noRd
check_reference <- function(reference, level, group, ci_level, resolved,
                            links, dots, call = rlang::caller_env()) {
  if (is.data.frame(reference)) {
    return(reference_from_tibble(reference, level, call = call))
  }
  if (is_one_fit(reference)) {
    return(reference_from_fit(
      reference, level, group, ci_level, resolved, links, dots,
      call = call
    ))
  }
  cli::cli_abort(
    c(
      "{.arg reference} must be a data frame or a fit with an \\
       {.fn extract_estimates} method, not \\
       {.obj_type_friendly {reference}}.",
      i = "A data frame needs {.field term} and {.field estimate}."
    ),
    call = call
  )
}

#' Join the two sides and compute the comparison
#'
#' `bias` and `covered` are `recover()`'s words for the same two
#' quantities, and `covered` points the same way in all three scorers:
#' the fit's interval covers the value it is compared with.
#' `overlap` is the extra question a reference with an interval can
#' answer, and is `NA` --- never `FALSE` --- for a row without one, so a
#' missing interval never reads as a disagreement.
#'
#' @noRd
cross_check_rows <- function(estimates, reference, level, scale,
                             call = rlang::caller_env()) {
  keys <- level_keys(level)
  joined <- join_truth(
    estimates, reference, keys,
    values = c("reference", "ref_low", "ref_high", "source"),
    arg = "reference",
    hint = "Both sides name parameters the way bmm does.",
    call = call
  )

  joined$bias <- joined$estimate - joined$reference
  joined$covered <- joined$reference >= joined$ci_low &
    joined$reference <= joined$ci_high

  has_interval <- !is.na(joined$ref_low) & !is.na(joined$ref_high)
  overlap <- joined$ref_low <= joined$ci_high &
    joined$ci_low <= joined$ref_high
  overlap[!has_interval] <- NA
  joined$overlap <- overlap

  joined$scale <- rep(scale, nrow(joined))
  joined
}

#' Compare a fit's estimates with a reference
#'
#' The question a model answers last: does it agree with what is already
#' known? `cross_check()` puts a fit's estimates beside a **reference**
#' --- a closed form such as `bmm::sdt_d()`, another implementation, or
#' published values --- and returns the comparison in the tibble shape
#' every other bmmtools scorer returns.
#'
#' **The reference is a comparison, not a truth.** `bias` is the signed
#' difference from it and `covered` says the fit's interval contains it;
#' neither asserts that the reference is correct. Where the two disagree,
#' the reference is as much a candidate for the error as the fit is, and
#' a reference derived under assumptions the fit does not make (an
#' equal-variance d′ against an unequal-variance fit, say) will disagree
#' for reasons that are nobody's bug.
#'
#' @param fit A fit with an [extract_estimates()] method, so a `brmsfit`
#'   and therefore a `bmmfit`. One fit: a list of them is an error
#'   pointing at [recover()], which scores replications against one truth
#'   and separates them with a `replication` column.
#' @param reference The right-hand side. A data frame with `term` and
#'   `estimate`, optionally `ci_low` and `ci_high` (both or neither),
#'   `source` and, at the subject level, `id`; its values are read on
#'   `scale`. Or another fit, which is extracted at the same level,
#'   transformed the same way, and marked `source = "fit"`.
#' @param scale `"natural"` compares on the scale a reader interprets,
#'   `"link"` on the scale the model was estimated on. The fit's
#'   estimates and interval bounds go through [inverse_link()] exactly as
#'   [recover()] transforms them, so a natural-scale `cross_check()` and a
#'   natural-scale `recover()` on the same fit agree to the last digit. A
#'   data-frame `reference` is taken to be on `scale` already and is not
#'   transformed.
#' @param links A named character vector mapping a term to one of the
#'   link names [inverse_link()] understands. `NULL` reads the link table
#'   from a `bmmfit`; with no table available the comparison falls back to
#'   the link scale and says so.
#' @param level `"population"` or `"subject"`. `"subject"` needs an `id`
#'   on both sides. `"sd"` and `"cor"` are not compared in this version:
#'   a correlation reference would need `var1`, `var2` and a
#'   pair-ordering rule of its own.
#' @param group The grouping variable subject-level estimates come from,
#'   passed to [extract_estimates()].
#' @param ci_level The interval mass, passed to [extract_estimates()].
#' @param ... Passed to [extract_estimates()] for both sides: `converged`,
#'   `ci_method`.
#'
#' @return A `bmmtools_cross_check` object: a tibble subclass with the
#'   columns `term`, `estimate`, `ci_low`, `ci_high`, `ci_method`,
#'   `ci_level`, `rhat`, `ess_bulk`, `ess_tail`, `reference`, `ref_low`,
#'   `ref_high`, `source`, `bias`, `covered`, `overlap`, `scale`,
#'   `level`, `id` and `converged`. `ref_low`, `ref_high` and `overlap`
#'   are `NA` for a reference without intervals. Call [summary()] on it
#'   for the per-parameter metrics and [plot_recovery()] for the picture.
#'
#' @details
#' A reference term the fit did not estimate is a warning listing what
#' the fit has, and is dropped; a fit term the reference does not mention
#' is dropped in silence, because a fit routinely estimates more than a
#' reference covers. No term matching at all is an error rather than an
#' empty result. At the subject level the same holds for ids.
#'
#' **d′ and d_a.** `bmm::sdt_d()` returns the equal-variance d′. bmm's
#' `sdt_yn` model estimates d_a = √2·δ/√(1 + r²), which equals d′ only
#' when `sdratio` is fixed at 0. Against a fit with a free `sdratio` the
#' reference is converted by the user, `d_a = d' * sqrt(2 / (1 +
#' exp(sdratio)^2))`; `cross_check()` compares what it is given and does
#' not guess a conversion, because the right one depends on which
#' `sdratio` the fit used.
#'
#' @examples
#' \dontrun{
#' # a signal-detection fit against the closed form, equal-variance case
#' fit <- bmm::bmm(
#'   bmm::bmf(d ~ 1, criterion ~ 1),
#'   data, bmm::sdt_yn(response = "resp", stimulus = "stim", n_trials = "n")
#' )
#' reference <- tibble::tibble(
#'   term = c("d", "criterion"),
#'   estimate = c(
#'     bmm::sdt_d(hit_rate = 0.8, fa_rate = 0.2),
#'     bmm::sdt_criterion(hit_rate = 0.8, fa_rate = 0.2)
#'   ),
#'   source = "closed form"
#' )
#' cross_check(fit, reference)
#'
#' # with sdratio free the fit estimates d_a, so the reference converts
#' sdratio <- 0.2
#' reference$estimate[1] <- reference$estimate[1] *
#'   sqrt(2 / (1 + exp(sdratio)^2))
#' cross_check(fit, reference)
#' }
#'
#' @export
cross_check <- function(fit,
                        reference,
                        scale = c("natural", "link"),
                        links = NULL,
                        level = c("population", "subject"),
                        group = "id",
                        ci_level = 0.95,
                        ...) {
  error_call <- rlang::current_env()
  scale <- rlang::arg_match(scale)
  level <- check_cross_check_level(level, call = error_call)
  check_ci_level(ci_level, call = error_call)
  check_cross_check_fit(fit, call = error_call)
  dots <- rlang::list2(...)

  # settled from the fit alone, so both sides land on one scale
  resolved <- resolve_links(fit, links, scale, call = error_call)
  estimates <- cross_check_side(
    fit, "fit", level, group, ci_level, resolved, dots,
    call = error_call
  )
  reference <- check_reference(
    reference, level, group, ci_level, resolved, links, dots,
    call = error_call
  )

  rows <- cross_check_rows(
    estimates, reference, level, resolved$scale,
    call = error_call
  )
  new_bmmtools_cross_check(
    rows,
    scale = resolved$scale,
    ci_level = ci_level,
    call = match.call(),
    error_call = error_call
  )
}
