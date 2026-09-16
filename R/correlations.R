# Correlation recovery.
#
# Three estimators of a between-subject correlation, all read from one
# fit: `model` is the group-level correlation the model estimates (the
# cor rows of extract_estimates()), `draws` correlates the subject values
# within each posterior draw, and `point` correlates the posterior means.
# A conjugate simulation of the three estimators is why all
# three are kept: under a posterior trade-off between two parameters,
# `point` is biased at rho = 0 while `draws` is not.
#
# The core works on a subject-draws array, a cor-estimates tibble and a
# covariate table, so that tests need no fit (as estimates_from_draws()
# does for extract_estimates()).

# subject draws -------------------------------------------------------------

#' Extract subject-level posterior draws as an array
#'
#' The per-draw subject values of every parameter with a group-level
#' effect: the population intercept plus the subject's deviation, summed
#' within each draw, on the link scale. This is what the `draws` and
#' `point` estimators of [extract_correlations()] read.
#'
#' @param fit A `brmsfit`, and so also a `bmmfit`. Other classes can
#'   provide a method that returns an array of the same shape.
#' @param group The grouping factor to read. `NULL` uses the fit's only
#'   grouping factor and errors if there is more than one.
#' @param ... Not used.
#'
#' @return A numeric array with the dimensions `iteration`, `chain`, `id`
#'   and `term`, named in its `dimnames`. `term` is the parameter name, or
#'   `<parameter>_<coefficient>` (`kappa_task1`) for a parameter with
#'   several group-level coefficients and no intercept, as in
#'   [extract_estimates()].
#'   Parameters the model fixed to a constant (zero posterior variance for
#'   every subject) are dropped. The grouping factor is stored in the
#'   attribute `group`.
#'
#' @examples
#' \dontrun{
#' draws <- extract_subject_draws(fit)
#' dim(draws)
#' dimnames(draws)$term
#' }
#'
#' @export
extract_subject_draws <- function(fit, group = NULL, ...) {
  UseMethod("extract_subject_draws")
}

#' @rdname extract_subject_draws
#' @export
extract_subject_draws.default <- function(fit, group = NULL, ...) {
  cli::cli_abort(
    c(
      "{.arg fit} must be a {.cls brmsfit}, not {.obj_type_friendly {fit}}.",
      i = "Its class is {.cls {class(fit)}}; another class needs an \\
           {.fn extract_subject_draws} method."
    )
  )
}

#' @rdname extract_subject_draws
#' @export
extract_subject_draws.brmsfit <- function(fit, group = NULL, ...) {
  rlang::check_installed("brms", "to extract subject draws from a fit.")
  rlang::check_dots_empty()
  subject_draws_from_draws(
    posterior::as_draws_array(fit), fit_groups(fit), group,
    call = rlang::current_env()
  )
}

#' The array behind extract_subject_draws(), from a draws object
#'
#' Built with `group_coefficients()` and `sum_group_draws()`, the same
#' route the subject level of `extract_estimates()` takes, so the two
#' cannot disagree about which draw belongs to which subject.
#'
#' @noRd
subject_draws_from_draws <- function(draws, groups, group = NULL,
                                     call = rlang::caller_env()) {
  group <- resolve_group(group, groups, call = call)
  coefficients <- group_coefficients(
    posterior::variables(draws), group,
    call = call
  )
  if (is.null(coefficients)) {
    cli::cli_abort(
      "The fit has no group-level coefficients for {.val {group}}.",
      call = call
    )
  }
  check_unique_terms(coefficients, c("term", "id"), call = call)

  summed <- unclass(sum_group_draws(draws, coefficients, call = call))
  ids <- unique(coefficients$id)
  terms <- unique(coefficients$term)
  n_iter <- dim(summed)[[1L]]
  n_chain <- dim(summed)[[2L]]
  out <- array(
    NA_real_,
    dim = c(n_iter, n_chain, length(ids), length(terms)),
    dimnames = list(
      iteration = as.character(seq_len(n_iter)),
      chain = as.character(seq_len(n_chain)),
      id = ids,
      term = terms
    )
  )
  for (k in seq_len(nrow(coefficients))) {
    out[, , coefficients$id[[k]], coefficients$term[[k]]] <- summed[, , k]
  }
  structure(drop_constant_terms(out), group = group)
}

#' Drop the terms a model fixed to a constant
#'
#' As in `drop_constant_rows()`: exactly zero posterior variance, and for
#' every subject, since a term is a column of the array. `NA` variance is
#' not zero and is kept.
#'
#' @noRd
drop_constant_terms <- function(x) {
  constant <- vapply(dimnames(x)[[4L]], function(term) {
    variances <- apply(x[, , , term, drop = FALSE], 3L, stats::var)
    all(vapply(variances, function(v) isTRUE(v == 0), logical(1)))
  }, logical(1))
  x[, , , !constant, drop = FALSE]
}

#' Validate a subject-draws array from any method
#' @noRd
check_subject_draws <- function(x, call = rlang::caller_env()) {
  if (!is.array(x) || !is.numeric(x) || length(dim(x)) != 4L) {
    cli::cli_abort(
      "Subject draws must be a 4-d numeric array, iteration by chain by \\
       id by term, not {.obj_type_friendly {x}}.",
      call = call
    )
  }
  ids <- dimnames(x)[[3L]]
  terms <- dimnames(x)[[4L]]
  bad <- is.null(ids) || is.null(terms) || anyNA(c(ids, terms)) ||
    anyDuplicated(ids) > 0L || anyDuplicated(terms) > 0L
  if (bad) {
    cli::cli_abort(
      "Subject draws must carry unique id and term names in their \\
       {.code dimnames}.",
      call = call
    )
  }
  x
}

# extract_correlations() ----------------------------------------------------

#' Extract between-subject correlations from a fit
#'
#' Three estimates of how two subject-level parameters, or a parameter
#' and an observed covariate, correlate across subjects:
#'
#' * `"model"`: the group-level correlation the model estimates, as in
#'   `extract_estimates(level = "cor")`. It exists only for parameters
#'   whose random effects share a correlation block (`(1 | p | id)`).
#' * `"draws"`: in each posterior draw, the Pearson correlation of the
#'   subject values across subjects; summarised like any other posterior
#'   (median, equal-tailed interval, rhat, ESS).
#' * `"point"`: the Pearson correlation of the subjects' posterior means,
#'   with a Fisher-z interval over the subjects.
#'
#' They answer different questions and do not agree. When the model
#' cannot separate two parameters within a subject, their posterior means
#' inherit that trade-off and `point` is biased even where the true
#' correlation is 0; `draws` is not. `point` is also attenuated by
#' shrinkage. See [recover_correlations()] for scoring them against the
#' truth.
#'
#' @param fit A `brmsfit` (so also a `bmmfit`), or any object with an
#'   [extract_subject_draws()] method. Or a set of separate fits, one per
#'   component: the result of [fit_components()], or a list of fits named
#'   by component; see "Separate fits" below.
#' @param estimator One or more of `"model"`, `"draws"` and `"point"`. All
#'   three by default.
#' @param pairs `NULL` for every pair of subject terms and covariates, or
#'   a character vector of pair terms such as `"kappa__thetat"`, in either
#'   order.
#' @param covariates `NULL`; a character vector of columns of `fit$data`
#'   that are constant within subject; or a data frame with the grouping
#'   column and one numeric column per covariate, one or several rows per
#'   subject.
#' @param group The grouping factor. `NULL` uses the fit's only one.
#' @param scale `"link"` or `"natural"`. On the natural scale each term is
#'   back-transformed with [inverse_link()] before correlating.
#' @param links A named character vector mapping a term to a link name.
#'   `NULL` reads the link table of a `bmmfit`; without one, the natural
#'   scale falls back to the link scale with a message.
#' @param ci_level The interval mass, a single number strictly between 0
#'   and 1.
#' @param converged Whether the fit passed the convergence gate. `NULL`
#'   computes it as [extract_estimates()] does for a `brmsfit` and leaves
#'   it `NA` for other classes; a logical scalar is used as given.
#' @param ... Not used.
#'
#' @return A tibble with the columns `term`, `var1`, `var2`, `estimator`,
#'   `estimate`, `ci_low`, `ci_high`, `ci_method`, `ci_level`, `rhat`,
#'   `ess_bulk`, `ess_tail`, `scale`, `n` (subjects) and `converged`.
#'   `term` is `var1__var2`, the two names sorted in the C locale.
#'
#' @details
#' **The model estimator** is on the link scale only. Under
#' `scale = "natural"` its rows are left out and a message says so. A
#' pair involving a covariate never has a model row, and neither does a
#' pair the model does not correlate.
#'
#' **The natural scale** for `point` takes the inverse link of each
#' subject's posterior mean on the link scale, not the mean of the inverse
#' links. Covariates are never transformed.
#'
#' **Intervals.** `ci_method` is `"eti"` for `model` and `draws` and
#' `"fisher_z"` for `point`, whose interval treats the posterior means as
#' if they were observed values; its `rhat` and ESS are `NA`.
#'
#' **Missing values.** A correlation is `NA`, never 0, with fewer than
#' three subjects or when one side has no spread across subjects. For
#' `draws` that is decided per draw, and the summary uses the draws that
#' remain.
#'
#' @section Separate fits:
#' With a set of fits, the subject draws of every fit are bound under
#' prefixed terms (`m3_c_task1`), subject by subject. Draws are paired by
#' iteration and chain, cut to the smallest counts among the fits with a
#' message. `model` rows exist only within a fit; `converged` of a pair is
#' whether every fit it involves converged, and may be given as a logical
#' named by component. `covariates` must be a data frame, because the
#' components' data carry none.
#'
#' No fit knows about the correlation between its parameters and another
#' model's, so correlations across fits are attenuated by the reliabilities
#' `rel_a` and `rel_b` of the two subject estimates: `point` is roughly
#' `rho * sqrt(rel_a * rel_b)` and `draws` roughly `rho * rel_a * rel_b`.
#' Neither recovers `rho`. A structural equation model of true against
#' estimated values, from [subject_table()], does; so would a joint
#' multivariate fit, which bmmtools does not support yet.
#'
#' @examples
#' \dontrun{
#' fit <- bmm::bmm(recovery_formula(model, re_cor = "all"), data, model)
#' extract_correlations(fit)
#' extract_correlations(fit, estimator = "draws", covariates = "G")
#' }
#'
#' @export
extract_correlations <- function(fit,
                                 estimator = c("model", "draws", "point"),
                                 pairs = NULL,
                                 covariates = NULL,
                                 group = NULL,
                                 scale = c("link", "natural"),
                                 links = NULL,
                                 ci_level = 0.95,
                                 converged = NULL,
                                 ...) {
  rlang::check_dots_empty()
  error_call <- rlang::current_env()
  estimator <- rlang::arg_match(
    estimator, c("model", "draws", "point"),
    multiple = TRUE
  )
  scale <- rlang::arg_match(scale)
  check_ci_level(ci_level, call = error_call)
  if (is_fit_set(fit, call = error_call)) {
    return(set_correlations(
      fit, estimator, pairs, covariates, group, scale, links, ci_level,
      converged,
      call = error_call
    ))
  }
  if (!is.null(converged)) check_converged(converged, call = error_call)
  resolved <- resolve_links(fit, links, scale, call = error_call)

  parts <- fit_correlation_parts(
    fit, estimator, resolved$scale, group, ci_level, converged,
    call = error_call
  )
  if (is.character(covariates)) {
    covariates <- covariates_from_data(
      covariates, fit$data, parts$group,
      call = error_call
    )
  }

  correlations_from_parts(
    parts$subject_draws, parts$cor_estimates, covariates,
    estimator = estimator, pairs = pairs,
    scale = resolved$scale, links = resolved$links,
    ci_level = ci_level, converged = parts$converged,
    n_subjects = parts$n_subjects, group = parts$group,
    call = error_call
  )
}

#' Read what the three estimators need from one fit
#'
#' A `brmsfit` is read once, as a draws array; any other class goes through
#' the `extract_subject_draws()` and `extract_estimates()` generics.
#'
#' @return A list with `subject_draws` and `cor_estimates` (each `NULL` when
#'   no requested estimator needs it), `converged`, `n_subjects` (counted in
#'   a brmsfit's data, else `NULL`) and `group`, the grouping column.
#' @noRd
fit_correlation_parts <- function(fit, estimator, scale, group, ci_level,
                                  converged, call = rlang::caller_env()) {
  is_brms <- inherits(fit, "brmsfit")
  draws <- NULL
  if (is_brms) {
    rlang::check_installed("brms", "to extract correlations from a fit.")
    draws <- posterior::as_draws_array(fit)
    group <- resolve_group(group, fit_groups(fit), call = call)
  }

  subject_draws <- NULL
  if (any(c("draws", "point") %in% estimator)) {
    subject_draws <- if (is_brms) {
      subject_draws_from_draws(draws, fit_groups(fit), group, call = call)
    } else {
      extract_subject_draws(fit, group = group)
    }
  }

  cor_estimates <- NULL
  if ("model" %in% estimator && identical(scale, "link")) {
    cor_estimates <- if (is_brms) {
      estimates_from_draws(
        draws, fit_groups(fit),
        level = "cor", group = group, ci_level = ci_level,
        converged = NA, ranef = fit$ranef, call = call
      )
    } else {
      extract_estimates(
        fit,
        level = "cor", group = group, ci_level = ci_level, converged = NA
      )
    }
  }

  if (is.null(converged)) {
    converged <- if (is_brms) fit_converged(fit, draws) else NA
  }
  group_column <- group %||% attr(subject_draws, "group") %||% "id"
  n_subjects <- NULL
  if (is_brms && is.data.frame(fit$data) && group_column %in% names(fit$data)) {
    n_subjects <- length(unique(fit$data[[group_column]]))
  }
  list(
    subject_draws = subject_draws, cor_estimates = cor_estimates,
    converged = converged, n_subjects = n_subjects, group = group_column
  )
}

# fit sets ----------------------------------------------------------------

#' Is `x` a set of separate fits, one per component?
#'
#' A `bmmtools_fit_set`, or a plain list of fits whose names are component
#' names. A plain list of fits with names that are not valid component names
#' is an error rather than a single fit, so that a typo does not end in a
#' message about the class of `fit`; an unnamed list is not a set.
#'
#' @noRd
is_fit_set <- function(x, call = rlang::caller_env()) {
  is_plain <- is.list(x) && !is.object(x) && length(x) > 0L &&
    all(vapply(x, is_one_fit, logical(1)))
  if (!inherits(x, "bmmtools_fit_set") && (!is_plain || is.null(names(x)))) {
    return(FALSE)
  }
  nms <- names(x)
  valid <- !is.null(nms) && anyDuplicated(nms) == 0L &&
    all(vapply(nms, valid_component_name, logical(1)))
  if (!valid) {
    cli::cli_abort(
      c(
        "A set of fits must be named by its component names, each unique, \\
         syntactic and without {.code _} or {.code .}.",
        i = "The names are {.val {nms}}."
      ),
      call = call
    )
  }
  TRUE
}

#' Every fit's bmm link table, under prefixed names, or `NULL`
#' @noRd
set_model_links <- function(fits) {
  tables <- lapply(names(fits), function(nm) {
    links <- model_links_of(fits[[nm]])
    if (is.null(links)) {
      return(NULL)
    }
    stats::setNames(as.character(links), prefix_terms(names(links), nm))
  })
  unlist(tables)
}

#' Prefix both names of pair terms (`kappa__thetat` to `a_kappa__a_thetat`)
#' @noRd
prefix_pair_terms <- function(terms, name) {
  vapply(strsplit(as.character(terms), "__", fixed = TRUE), function(p) {
    paste(prefix_terms(p, name), collapse = "__")
  }, character(1))
}

#' Validate `converged` for a set: a named logical per component, or `NULL`
#' @noRd
check_set_converged <- function(converged, components,
                                call = rlang::caller_env()) {
  if (is.null(converged)) {
    return(NULL)
  }
  nms <- names(converged)
  if (is.logical(converged) && length(converged) == 1L && is.null(nms)) {
    return(stats::setNames(rep(converged, length(components)), components))
  }
  by_component <- is.logical(converged) && !is.null(nms) &&
    length(converged) == length(components) && anyDuplicated(nms) == 0L &&
    setequal(nms, components)
  if (!by_component) {
    cli::cli_abort(
      c(
        "{.arg converged} must be {.code NULL}, a single logical, or a \\
         logical named by component.",
        i = "The components are {.val {components}}."
      ),
      call = call
    )
  }
  converged[components]
}

#' Bind the subject draws of separate fits along the term dimension
#'
#' Terms are prefixed with the component name. Subjects are put in the
#' first fit's order, and iterations and chains are cut to the smallest
#' counts, since draws of independent fits can be paired in any order.
#'
#' @param arrays A named list of subject-draws arrays, one per component.
#' @noRd
bind_subject_draws <- function(arrays, call = rlang::caller_env()) {
  nms <- names(arrays)
  arrays <- lapply(arrays, check_subject_draws, call = call)
  ids <- dimnames(arrays[[1L]])[[3L]]
  for (nm in nms[-1L]) {
    other <- dimnames(arrays[[nm]])[[3L]]
    if (length(other) != length(ids) || !setequal(other, ids)) {
      cli::cli_abort(
        c(
          "The fits of components {.val {nms[[1L]]}} and {.val {nm}} have \\
           different subjects.",
          i = "Separate fits are correlated subject by subject, so they \\
               must hold the same subjects."
        ),
        call = call
      )
    }
  }
  n_iter <- vapply(arrays, function(x) dim(x)[[1L]], integer(1))
  n_chain <- vapply(arrays, function(x) dim(x)[[2L]], integer(1))
  keep_iter <- min(n_iter)
  keep_chain <- min(n_chain)
  larger <- nms[n_iter > keep_iter | n_chain > keep_chain]
  if (length(larger) > 0L) {
    cli::cli_inform(c(
      "Kept {keep_iter} iterations and {keep_chain} chains per fit.",
      i = "Those are the smallest counts among the fits; {.val {larger}} \\
           had more. Draws of separate fits are paired by iteration and \\
           chain."
    ))
  }

  terms <- lapply(stats::setNames(nm = nms), function(nm) {
    prefix_terms(dimnames(arrays[[nm]])[[4L]], nm)
  })
  out <- array(
    NA_real_,
    dim = c(keep_iter, keep_chain, length(ids), length(unlist(terms))),
    dimnames = list(
      iteration = as.character(seq_len(keep_iter)),
      chain = as.character(seq_len(keep_chain)),
      id = ids,
      term = unname(unlist(terms))
    )
  )
  for (nm in nms) {
    if (length(terms[[nm]]) == 0L) next
    out[, , , terms[[nm]]] <- arrays[[nm]][
      seq_len(keep_iter), seq_len(keep_chain), ids, ,
      drop = FALSE
    ]
  }
  structure(out, group = attr(arrays[[1L]], "group"))
}

#' The subject draws of a set of fits, bound
#' @noRd
set_subject_draws <- function(fits, group = NULL, call = rlang::caller_env()) {
  arrays <- lapply(stats::setNames(nm = names(fits)), function(nm) {
    extract_subject_draws(fits[[nm]], group = group)
  })
  bind_subject_draws(arrays, call = call)
}

#' Whether each pair's fits converged
#'
#' A pair takes `all()` of the verdicts of the fits it involves; a pair of
#' two covariates involves none and takes every fit's.
#'
#' @noRd
pair_converged <- function(var1, var2, verdicts) {
  owner <- function(term) {
    prefix <- sub("_.*$", "", term)
    ifelse(
      grepl("_", term, fixed = TRUE) & prefix %in% names(verdicts),
      prefix, NA_character_
    )
  }
  first <- owner(var1)
  second <- owner(var2)
  vapply(seq_along(var1), function(k) {
    involved <- unique(c(first[[k]], second[[k]]))
    involved <- involved[!is.na(involved)]
    if (length(involved) == 0L) all(verdicts) else all(verdicts[involved])
  }, logical(1))
}

#' extract_correlations() for a set of separate fits
#' @noRd
set_correlations <- function(fits, estimator, pairs, covariates, group, scale,
                             links, ci_level, converged, call) {
  components <- names(fits)
  converged <- check_set_converged(converged, components, call = call)
  if (!is.null(covariates) && !is.data.frame(covariates)) {
    cli::cli_abort(
      c(
        "With a set of fits, {.arg covariates} must be a data frame with the \\
         grouping column and one numeric column per covariate.",
        i = "The components' data carry no covariates; a simulation set \\
             keeps them in {.field covariate_data}."
      ),
      call = call
    )
  }
  resolved <- resolve_links(
    NULL, links %||% set_model_links(fits), scale,
    call = call
  )

  parts <- lapply(stats::setNames(nm = components), function(nm) {
    fit_correlation_parts(
      fits[[nm]], estimator, resolved$scale, group, ci_level,
      converged[[nm]],
      call = call
    )
  })

  subject_draws <- NULL
  if (any(c("draws", "point") %in% estimator)) {
    subject_draws <- bind_subject_draws(
      lapply(parts, `[[`, "subject_draws"),
      call = call
    )
  }
  # model correlations exist within a fit only
  model_rows <- lapply(components, function(nm) {
    rows <- parts[[nm]]$cor_estimates
    if (is.null(rows) || nrow(rows) == 0L) {
      return(NULL)
    }
    rows <- tibble::as_tibble(rows)
    if ("level" %in% names(rows)) rows <- rows[rows$level %in% "cor", ]
    rows$term <- prefix_pair_terms(rows$term, nm)
    rows[setdiff(names(rows), c("var1", "var2", "level"))]
  })
  cor_estimates <- NULL
  if (!all(vapply(model_rows, is.null, logical(1)))) {
    cor_estimates <- dplyr::bind_rows(model_rows)
  }
  n_subjects <- NULL
  for (p in parts) n_subjects <- n_subjects %||% p$n_subjects

  out <- correlations_from_parts(
    subject_draws, cor_estimates, covariates,
    estimator = estimator, pairs = pairs,
    scale = resolved$scale, links = resolved$links,
    ci_level = ci_level, converged = NA,
    n_subjects = n_subjects, group = parts[[1L]]$group,
    call = call
  )
  verdicts <- vapply(parts, function(p) as.logical(p$converged), logical(1))
  out$converged <- pair_converged(out$var1, out$var2, verdicts)
  out
}

#' Validate `ci_level`: one number strictly between 0 and 1
#' @noRd
check_ci_level <- function(ci_level, call = rlang::caller_env()) {
  if (!is.numeric(ci_level) || length(ci_level) != 1L || is.na(ci_level)) {
    cli::cli_abort(
      "{.arg ci_level} must be a single number, \\
       not {.obj_type_friendly {ci_level}}.",
      call = call
    )
  }
  if (ci_level <= 0 || ci_level >= 1) {
    cli::cli_abort(
      "{.arg ci_level} must be between 0 and 1, not {ci_level}.",
      call = call
    )
  }
  invisible(ci_level)
}

#' Turn a character vector of covariate names into a covariate table
#' @noRd
covariates_from_data <- function(names, data, group,
                                 call = rlang::caller_env()) {
  if (!is.data.frame(data)) {
    cli::cli_abort(
      c(
        "The fit carries no data to read {.arg covariates} from.",
        i = "Pass a data frame with the {.val {group}} column instead."
      ),
      call = call
    )
  }
  missing <- setdiff(c(group, names), names(data))
  if (length(missing) > 0L) {
    cli::cli_abort(
      "The fit's data have no column{?s} {.val {missing}}.",
      call = call
    )
  }
  as.data.frame(data)[c(group, names)]
}

#' The correlation table from its parts
#'
#' @param subject_draws `NULL` or an array from `extract_subject_draws()`.
#' @param cor_estimates `NULL` or an estimates tibble whose `cor` rows are
#'   the model estimator.
#' @param covariates `NULL` or a data frame with the `group` column.
#' @param n_subjects Used for `n` of the model rows when there is no
#'   subject-draws array to count.
#' @noRd
correlations_from_parts <- function(subject_draws = NULL,
                                    cor_estimates = NULL,
                                    covariates = NULL,
                                    estimator = c("model", "draws", "point"),
                                    pairs = NULL,
                                    scale = "link",
                                    links = NULL,
                                    ci_level = 0.95,
                                    converged = NA,
                                    n_subjects = NULL,
                                    group = NULL,
                                    call = rlang::caller_env()) {
  if (!is.null(subject_draws)) {
    subject_draws <- check_subject_draws(subject_draws, call = call)
  }
  group <- group %||% attr(subject_draws, "group") %||% "id"
  model_rows <- model_cor_rows(cor_estimates, call = call)

  ids <- dimnames(subject_draws)[[3L]]
  terms <- if (is.null(subject_draws)) {
    unique(c(model_rows$var1, model_rows$var2))
  } else {
    dimnames(subject_draws)[[4L]]
  }
  cov <- check_cor_covariates(covariates, group, ids, terms, call = call)
  pair_table <- resolve_pairs(pairs, c(terms, colnames(cov)), call = call)
  n <- if (is.null(subject_draws)) {
    as.integer(n_subjects %||% NA_integer_)
  } else {
    length(ids)
  }

  pieces <- lapply(estimator, function(est) {
    switch(est,
      model = model_estimator_rows(model_rows, pair_table, scale, n),
      draws = draws_estimator_rows(
        subject_draws, cov, pair_table, scale, links, ci_level
      ),
      point = point_estimator_rows(
        subject_draws, cov, pair_table, scale, links, ci_level
      )
    )
  })
  if ("model" %in% estimator && identical(scale, "natural")) {
    cli::cli_inform(c(
      "Model correlations are estimated on the link scale only, so no \\
       {.val model} rows are returned on the natural scale.",
      i = "Use {.code scale = \"link\"} for them."
    ))
  }

  out <- dplyr::bind_rows(empty_correlation_rows(), pieces)
  out$ci_level <- rep(as.double(ci_level), nrow(out))
  out$scale <- rep(scale, nrow(out))
  out$converged <- rep(converged, nrow(out))
  out[correlation_columns()]
}

#' The columns of extract_correlations(), in order
#' @noRd
correlation_columns <- function() {
  c(
    "term", "var1", "var2", "estimator", "estimate", "ci_low", "ci_high",
    "ci_method", "ci_level", "rhat", "ess_bulk", "ess_tail", "scale", "n",
    "converged"
  )
}

#' @noRd
empty_correlation_rows <- function() {
  tibble::tibble(
    term = character(), var1 = character(), var2 = character(),
    estimator = character(), estimate = double(), ci_low = double(),
    ci_high = double(), ci_method = character(), rhat = double(),
    ess_bulk = double(), ess_tail = double(), n = integer()
  )
}

#' The cor rows of an estimates tibble, with their pair terms normalised
#'
#' A tibble without rows, as a fit class without a cor level returns, is
#' no model estimate at all rather than an error.
#'
#' @noRd
model_cor_rows <- function(cor_estimates, call = rlang::caller_env()) {
  empty <- tibble::tibble(
    term = character(), var1 = character(), var2 = character(),
    estimate = double(), ci_low = double(), ci_high = double(),
    rhat = double(), ess_bulk = double(), ess_tail = double()
  )
  if (is.null(cor_estimates) || nrow(cor_estimates) == 0L) {
    return(empty)
  }
  needed <- c(
    "term", "estimate", "ci_low", "ci_high", "rhat", "ess_bulk", "ess_tail"
  )
  missing <- setdiff(needed, names(cor_estimates))
  if (length(missing) > 0L) {
    cli::cli_abort(
      "The correlation estimates are missing the column{?s} \\
       {.val {missing}}.",
      call = call
    )
  }
  rows <- tibble::as_tibble(cor_estimates)
  if ("level" %in% names(rows)) rows <- rows[rows$level %in% "cor", ]
  parts <- strsplit(rows$term, "__", fixed = TRUE)
  rows <- rows[lengths(parts) == 2L, ]
  parts <- parts[lengths(parts) == 2L]
  named <- lapply(parts, function(p) pair_term(p[[1L]], p[[2L]]))
  tibble::tibble(
    term = vapply(named, `[[`, character(1), "term"),
    var1 = vapply(named, `[[`, character(1), "var1"),
    var2 = vapply(named, `[[`, character(1), "var2"),
    estimate = as.double(rows$estimate),
    ci_low = as.double(rows$ci_low),
    ci_high = as.double(rows$ci_high),
    rhat = as.double(rows$rhat),
    ess_bulk = as.double(rows$ess_bulk),
    ess_tail = as.double(rows$ess_tail)
  )
}

#' Validate the covariate table and reduce it to one row per subject
#'
#' @return `NULL`, or a numeric matrix subjects by covariates, in the
#'   order of `ids` when there are subject draws.
#' @noRd
check_cor_covariates <- function(covariates, group, ids, terms,
                                 call = rlang::caller_env()) {
  if (is.null(covariates)) {
    return(NULL)
  }
  if (!is.data.frame(covariates)) {
    cli::cli_abort(
      "{.arg covariates} must be {.code NULL}, a character vector of data \\
       columns, or a data frame, not {.obj_type_friendly {covariates}}.",
      call = call
    )
  }
  if (!group %in% names(covariates)) {
    cli::cli_abort(
      "{.arg covariates} must have the grouping column {.val {group}}.",
      call = call
    )
  }
  cov_names <- setdiff(names(covariates), group)
  if (length(cov_names) == 0L) {
    cli::cli_abort(
      "{.arg covariates} has no covariate column besides {.val {group}}.",
      call = call
    )
  }
  non_numeric <- cov_names[
    !vapply(cov_names, function(nm) is.numeric(covariates[[nm]]), logical(1))
  ]
  if (length(non_numeric) > 0L) {
    cli::cli_abort(
      "{.arg covariates} column{?s} {.val {non_numeric}} must be numeric.",
      call = call
    )
  }
  clash <- intersect(cov_names, terms)
  if (length(clash) > 0L) {
    cli::cli_abort(
      "{.arg covariates} column{?s} {.val {clash}} {?has/have} the name of \\
       a subject term.",
      call = call
    )
  }
  separator <- cov_names[grepl("__", cov_names, fixed = TRUE)]
  if (length(separator) > 0L) {
    cli::cli_abort(
      "{.arg covariates} name{?s} {.val {separator}} contain{?s/} {.code __}, \\
       the separator of a pair term.",
      call = call
    )
  }
  with_na <- cov_names[
    vapply(cov_names, function(nm) anyNA(covariates[[nm]]), logical(1))
  ]
  if (length(with_na) > 0L) {
    cli::cli_abort(
      "{.arg covariates} column{?s} {.val {with_na}} {?has/have} missing \\
       values.",
      call = call
    )
  }

  key <- as.character(covariates[[group]])
  varying <- cov_names[vapply(cov_names, function(nm) {
    per_subject <- tapply(covariates[[nm]], key, function(v) {
      length(unique(v)) > 1L
    })
    any(per_subject)
  }, logical(1))]
  if (length(varying) > 0L) {
    cli::cli_abort(
      c(
        "{.arg covariates} {.val {varying}} {?is/are} not constant within \\
         {.val {group}}.",
        i = "A covariate is a person variable: one value per subject."
      ),
      call = call
    )
  }

  subjects <- unique(key)
  first <- match(subjects, key)
  values <- matrix(
    vapply(cov_names, function(nm) {
      as.double(covariates[[nm]][first])
    }, double(length(first))),
    nrow = length(first),
    dimnames = list(subjects, cov_names)
  )
  if (is.null(ids)) {
    return(values)
  }
  absent <- setdiff(ids, subjects)
  if (length(absent) > 0L) {
    cli::cli_abort(
      "{.arg covariates} has no value for {.val {group}} {.val {absent}}.",
      call = call
    )
  }
  values[ids, , drop = FALSE]
}

#' Resolve `pairs` against the terms available
#'
#' @return A tibble `term`, `var1`, `var2`, one row per distinct pair.
#' @noRd
resolve_pairs <- function(pairs, universe, call = rlang::caller_env()) {
  empty <- tibble::tibble(
    term = character(), var1 = character(), var2 = character()
  )
  if (is.null(pairs)) {
    if (length(universe) < 2L) {
      return(empty)
    }
    combos <- utils::combn(length(universe), 2L)
    rows <- lapply(seq_len(ncol(combos)), function(k) {
      pair_term(universe[[combos[1L, k]]], universe[[combos[2L, k]]])
    })
    return(dplyr::bind_rows(empty, rows))
  }
  if (!is.character(pairs) || anyNA(pairs) || length(pairs) == 0L) {
    cli::cli_abort(
      "{.arg pairs} must be {.code NULL} or a character vector of pair \\
       terms such as {.val kappa__thetat}.",
      call = call
    )
  }
  parts <- strsplit(pairs, "__", fixed = TRUE)
  malformed <- pairs[vapply(parts, function(p) {
    length(p) != 2L || identical(p[[1L]], p[[2L]])
  }, logical(1))]
  if (length(malformed) > 0L) {
    cli::cli_abort(
      "{.arg pairs} entr{?y/ies} {.val {malformed}} must name two different \\
       terms joined by {.code __}.",
      call = call
    )
  }
  unknown <- setdiff(unlist(parts), universe)
  if (length(unknown) > 0L) {
    cli::cli_abort(
      c(
        "{.arg pairs} names {.val {unknown}}, which {?is/are} neither a \\
         subject term nor a covariate.",
        i = "Available: {.val {universe}}."
      ),
      call = call
    )
  }
  rows <- dplyr::bind_rows(
    empty,
    lapply(parts, function(p) pair_term(p[[1L]], p[[2L]]))
  )
  rows[!duplicated(rows$term), ]
}

#' @noRd
model_estimator_rows <- function(model_rows, pair_table, scale, n) {
  if (identical(scale, "natural")) {
    return(NULL)
  }
  keep <- intersect(pair_table$term, model_rows$term)
  rows <- model_rows[match(keep, model_rows$term), ]
  rows$estimator <- rep("model", nrow(rows))
  rows$ci_method <- rep("eti", nrow(rows))
  rows$n <- rep(n, nrow(rows))
  rows
}

#' The values of one term, draws by subjects, on the requested scale
#'
#' Rows are draws in the order iteration within chain, which is how a
#' vector of per-draw results folds back into an iteration by chain array.
#' A covariate is the same row in every draw.
#'
#' @noRd
term_draw_matrix <- function(subject_draws, cov, term, scale, links) {
  d <- dim(subject_draws)
  n_draws <- d[[1L]] * d[[2L]]
  if (!is.null(cov) && term %in% colnames(cov)) {
    return(matrix(cov[, term], nrow = n_draws, ncol = d[[3L]], byrow = TRUE))
  }
  values <- matrix(subject_draws[, , , term], nrow = n_draws, ncol = d[[3L]])
  if (identical(scale, "natural")) {
    values[] <- inverse_link(as.vector(values), link_of(term, links))
  }
  values
}

#' Pearson correlation of two matrices, row by row
#'
#' Centred, scaled and cross-multiplied, so every draw is done in one
#' pass. A row with no spread gives `NA`; "no spread" is judged relative
#' to the size of the values, because the mean of identical doubles need
#' not equal them exactly and would leave a rounding residue.
#'
#' @noRd
row_correlations <- function(a, b) {
  if (ncol(a) < 3L) {
    return(rep(NA_real_, nrow(a)))
  }
  ac <- a - rowMeans(a)
  bc <- b - rowMeans(b)
  sa <- sqrt(rowSums(ac^2))
  sb <- sqrt(rowSums(bc^2))
  tolerance <- sqrt(.Machine$double.eps)
  flat <- !(sa > tolerance * sqrt(rowSums(a^2))) |
    !(sb > tolerance * sqrt(rowSums(b^2)))
  r <- rowSums(ac * bc) / (sa * sb)
  r[flat] <- NA_real_
  r
}

#' @noRd
draws_estimator_rows <- function(subject_draws, cov, pair_table, scale,
                                 links, ci_level) {
  if (is.null(subject_draws) || nrow(pair_table) == 0L) {
    return(NULL)
  }
  d <- dim(subject_draws)
  per_draw <- vapply(seq_len(nrow(pair_table)), function(k) {
    row_correlations(
      term_draw_matrix(subject_draws, cov, pair_table$var1[[k]], scale, links),
      term_draw_matrix(subject_draws, cov, pair_table$var2[[k]], scale, links)
    )
  }, double(d[[1L]] * d[[2L]]))
  arr <- array(
    per_draw,
    dim = c(d[[1L]], d[[2L]], nrow(pair_table)),
    dimnames = list(
      iteration = NULL, chain = NULL, variable = pair_table$term
    )
  )
  summary <- summarise_selected(posterior::as_draws_array(arr), ci_level)
  tibble::tibble(
    term = pair_table$term,
    var1 = pair_table$var1,
    var2 = pair_table$var2,
    estimator = "draws",
    estimate = unname(summary$estimate),
    ci_low = unname(summary$ci_low),
    ci_high = unname(summary$ci_high),
    ci_method = "eti",
    rhat = unname(summary$rhat),
    ess_bulk = unname(summary$ess_bulk),
    ess_tail = unname(summary$ess_tail),
    n = d[[3L]]
  )
}

#' @noRd
point_estimator_rows <- function(subject_draws, cov, pair_table, scale,
                                 links, ci_level) {
  if (is.null(subject_draws) || nrow(pair_table) == 0L) {
    return(NULL)
  }
  means <- colMeans(subject_draws, dims = 2L)
  if (identical(scale, "natural")) {
    for (term in colnames(means)) {
      means[, term] <- inverse_link(means[, term], link_of(term, links))
    }
  }
  values <- cbind(means, cov)
  rows <- lapply(seq_len(nrow(pair_table)), function(k) {
    r <- metric_r(
      unname(values[, pair_table$var1[[k]]]),
      unname(values[, pair_table$var2[[k]]]),
      ci_level = ci_level
    )
    tibble::tibble(
      estimate = r$r, ci_low = r$r_low, ci_high = r$r_high,
      n = as.integer(r$n)
    )
  })
  stats <- dplyr::bind_rows(rows)
  tibble::tibble(
    term = pair_table$term,
    var1 = pair_table$var1,
    var2 = pair_table$var2,
    estimator = "point",
    estimate = stats$estimate,
    ci_low = stats$ci_low,
    ci_high = stats$ci_high,
    ci_method = "fisher_z",
    rhat = NA_real_,
    ess_bulk = NA_real_,
    ess_tail = NA_real_,
    n = stats$n
  )
}

# recover_correlations() ----------------------------------------------------

#' Score recovered correlations against the generating ones
#'
#' The correlation counterpart of [recover()]: one row per replication,
#' pair and estimator, with the estimate, its interval, the generating
#' correlation `true_value`, the correlation the simulated subjects
#' actually had `sample_value`, the errors against both and whether the
#' interval excludes zero.
#'
#' Two truths are carried because they answer different questions. A
#' finite sample of subjects does not have the generating correlation;
#' `sample_value` is what an estimator could at best recover from these
#' subjects, and `true_value` is what a reader wants to learn about.
#'
#' @param fits A fit (a `brmsfit`, or any object with an
#'   [extract_subject_draws()] method), a list of fits with one element
#'   per replication, or a tibble from [extract_correlations()] that may
#'   carry `replication` and `condition` columns. With a simulation set as
#'   `truth`, a set of fits from [fit_components()] (or a list of fits
#'   named by component), or an unnamed list of such sets, one per
#'   replication.
#' @param truth A `bmmtools_simulation`; a `bmmtools_simulation_set` from
#'   [simulate_components()]; a list of either, parallel to `fits`; or a
#'   list with the tibbles `cor` (`term`, `true_value`),
#'   `subjects` (`id`, `term`, `true_value`) and optionally `covariates`
#'   (`id`, `term`, `true_value`), all on the link scale and each of them
#'   optionally carrying `replication` or `condition`. Without a
#'   `replication` column, the same truth applies to every replication.
#' @param estimator One or more of `"model"`, `"draws"` and `"point"`.
#' @param scale `"link"` or `"natural"`. When `fits` is a tibble and
#'   `scale` is not given, the tibble's own scale is used.
#' @param links A named character vector mapping a term to a link name.
#'   `NULL` reads the link table of a `bmmfit`, then of the simulation's
#'   model; without one, the natural scale falls back to the link scale
#'   with a message.
#' @param pairs `NULL` for every pair in `truth`, or a character vector of
#'   pair terms; restricts both the estimates and the truth.
#' @param group The grouping factor, `"id"` by default.
#' @param ci_level The interval mass, passed to [extract_correlations()].
#' @param ... Passed to [extract_correlations()], for example
#'   `converged`.
#'
#' @return A `bmmtools_cor_recovery` object: a tibble subclass with the
#'   columns `term`, `var1`, `var2`, `estimator`, `estimate`, `ci_low`,
#'   `ci_high`, `ci_method`, `ci_level`, `rhat`, `ess_bulk`, `ess_tail`,
#'   `true_value`, `sample_value`, `bias`, `bias_sample`, `covered`,
#'   `covered_sample`, `excludes_zero`, `scale`, `n`, `converged`,
#'   `condition` and `replication`. Call
#'   [summary()][summary.bmmtools_cor_recovery()] on it for the
#'   per-pair metrics and [plot_recovery()] to draw it.
#'
#' @details
#' `bias` is `estimate - true_value` and `bias_sample` is `estimate -
#' sample_value`; `covered` and `covered_sample` use the closed interval;
#' `excludes_zero` is `ci_high < 0 | ci_low > 0`.
#'
#' **On the natural scale** `sample_value` is recomputed from the subject
#' values after [inverse_link()], covariates untransformed. A generating
#' correlation on the link scale has no natural-scale counterpart that
#' one transform gives, so `true_value` is `0` where the link-scale value
#' is 0 (independence survives a monotone transform) and `NA` otherwise.
#' The `model` estimator has no natural-scale rows; with fits as input
#' it is dropped with a message.
#'
#' A truth pair no estimate matches produces a warning and is dropped;
#' every pair unmatched is an error.
#'
#' **Separate fits.** A correlation between parameters of two separately
#' fitted models is attenuated by the reliabilities of both subject
#' estimates, so `bias` against `true_value` is expected even when every
#' fit is sound; see the section "Separate fits" of
#' [extract_correlations()].
#'
#' @examples
#' extracted <- tibble::tibble(
#'   term = "kappa__thetat", var1 = "kappa", var2 = "thetat",
#'   estimator = "draws", estimate = c(0.42, 0.18, 0.55),
#'   ci_low = c(0.1, -0.2, 0.2), ci_high = c(0.7, 0.5, 0.8),
#'   ci_method = "eti", ci_level = 0.95, rhat = 1, ess_bulk = 800,
#'   ess_tail = 800, scale = "link", n = 5L, replication = 1:3
#' )
#' subjects <- tibble::tibble(
#'   id = rep(as.character(1:5), 6),
#'   term = rep(rep(c("kappa", "thetat"), each = 5), 3),
#'   true_value = c(
#'     1.2, 0.4, 2.0, 1.1, 0.7, 0.3, -0.5, 0.9, 0.1, -0.2,
#'     0.8, 1.5, 1.1, 0.2, 1.9, 0.0, 0.6, 0.4, -0.3, 0.8,
#'     1.4, 0.9, 0.3, 1.8, 1.0, 0.5, 0.2, -0.4, 0.9, 0.1
#'   ),
#'   replication = rep(1:3, each = 10)
#' )
#' truth <- list(
#'   cor = tibble::tibble(term = "kappa__thetat", true_value = 0.5),
#'   subjects = subjects
#' )
#' recovery <- recover_correlations(extracted, truth)
#' summary(recovery)
#'
#' @export
recover_correlations <- function(fits,
                                 truth,
                                 estimator = c("model", "draws", "point"),
                                 scale = c("link", "natural"),
                                 links = NULL,
                                 pairs = NULL,
                                 group = "id",
                                 ci_level = 0.95,
                                 ...) {
  error_call <- rlang::current_env()
  estimator <- rlang::arg_match(
    estimator, c("model", "draws", "point"),
    multiple = TRUE
  )
  scale_given <- !missing(scale)
  scale <- rlang::arg_match(scale)
  check_ci_level(ci_level, call = error_call)

  kind <- if (is_simulation_set_truth(truth)) {
    cor_set_fits_kind(fits, call = error_call)
  } else {
    cor_fits_kind(fits, call = error_call)
  }
  truths <- cor_truth_input(truth, fits, kind, call = error_call)
  if (identical(kind, "tibble")) {
    fits <- check_extracted_correlations(fits, call = error_call)
    own_scale <- unique(fits$scale)
    if (!scale_given && length(own_scale) == 1L) scale <- own_scale
  }
  if (!is.null(pairs)) {
    pairs <- resolve_pairs(pairs, unique(unlist(
      strsplit(pairs, "__", fixed = TRUE)
    )), call = error_call)$term
  }

  model_links <- switch(kind,
    fit = ,
    list = model_links_of(fits),
    set = set_model_links(fits),
    sets = set_model_links(fits[[1L]]),
    NULL
  )
  if (is.null(links) && is.null(model_links) && length(truths$links) > 0L) {
    links <- unlist(truths$links)
  }
  # a set's tables are prefixed here, so they reach resolve_links() as
  # `links`; model_links_of() would read one fit's unprefixed table
  link_source <- NULL
  if (kind %in% c("fit", "list")) link_source <- fits
  if (kind %in% c("set", "sets")) links <- links %||% model_links
  resolved <- resolve_links(link_source, links, scale, call = error_call)

  estimates <- if (identical(kind, "tibble")) {
    tibble_cor_estimates(fits, estimator, resolved$scale, pairs, error_call)
  } else {
    if (identical(resolved$scale, "natural") && "model" %in% estimator) {
      estimator <- setdiff(estimator, "model")
      if (length(estimator) == 0L) {
        cli::cli_abort(
          c(
            "The {.val model} estimator is on the link scale only, so with \\
             {.code scale = \"natural\"} nothing is left to score.",
            i = "Use {.code scale = \"link\"}, or add {.val draws} or \\
                 {.val point}."
          ),
          call = error_call
        )
      }
      cli::cli_inform(c(
        "Model correlations are estimated on the link scale only, so no \\
         {.val model} rows are scored on the natural scale."
      ))
    }
    fit_cor_estimates(
      fits, kind, truths, estimator, resolved, pairs, group, ci_level,
      list(...), error_call
    )
  }

  truth_cor <- truths$cor
  if (!is.null(pairs)) truth_cor <- truth_cor[truth_cor$term %in% pairs, ]
  keys <- c("term", intersect(c("replication", "condition"), names(truth_cor)))
  # a truth split by condition cannot be matched to estimates that do not
  # say which condition they belong to; broadcasting would pair every fit
  # with every condition's truth
  unkeyed <- setdiff(keys, names(estimates))
  if (length(unkeyed) > 0L) {
    cli::cli_abort(
      c(
        "{.arg truth} carries the column{?s} {.val {unkeyed}}, but the \\
         correlation estimates do not.",
        i = "Add {.val {unkeyed}} to an extracted tibble passed as \\
             {.arg fits}, or drop {?it/them} from {.arg truth} when it \\
             describes a single condition."
      ),
      call = error_call
    )
  }
  joined <- join_truth(estimates, truth_cor, keys, call = error_call)

  samples <- sample_correlations(
    truth_cor, truths$subjects, truths$covariates,
    resolved$scale, resolved$links
  )
  sample_keys <- c(
    "term", intersect(c("replication", "condition"), names(samples))
  )
  joined <- dplyr::left_join(
    joined, samples,
    by = sample_keys, relationship = "many-to-one"
  )

  if (identical(resolved$scale, "natural")) {
    joined$true_value <- ifelse(joined$true_value == 0, 0, NA_real_)
  }
  joined$bias <- joined$estimate - joined$true_value
  joined$bias_sample <- joined$estimate - joined$sample_value
  joined$covered <- joined$true_value >= joined$ci_low &
    joined$true_value <= joined$ci_high
  joined$covered_sample <- joined$sample_value >= joined$ci_low &
    joined$sample_value <= joined$ci_high
  joined$excludes_zero <- joined$ci_high < 0 | joined$ci_low > 0
  joined$n <- as.integer(joined$n)
  if ("condition" %in% names(joined)) {
    joined$condition <- as.character(joined$condition)
  }

  new_bmmtools_cor_recovery(
    joined,
    scale = resolved$scale,
    ci_level = ci_level,
    call = match.call(),
    error_call = error_call
  )
}

#' Is this one fit, a list of fits or an extracted tibble?
#'
#' A fit is any classed object that is not a data frame, so that a class
#' with its own `extract_subject_draws()` method (the suite's mock) is
#' accepted like a `brmsfit`.
#'
#' @noRd
cor_fits_kind <- function(fits, call = rlang::caller_env()) {
  if (is.data.frame(fits)) {
    return("tibble")
  }
  if (inherits(fits, "bmmtools_fit_set")) {
    cli::cli_abort(
      c(
        "A set of fits needs a {.cls bmmtools_simulation_set} as \\
         {.arg truth}.",
        i = "Pass the result of {.fn simulate_components} as {.arg truth}."
      ),
      call = call
    )
  }
  if (is_one_fit(fits)) {
    return("fit")
  }
  is_list <- is.list(fits) && !is.object(fits) && length(fits) > 0L &&
    all(vapply(fits, is_one_fit, logical(1)))
  if (is_list) {
    return("list")
  }
  cli::cli_abort(
    c(
      "{.arg fits} must be a fit, a list of fits, or a tibble from \\
       {.fn extract_correlations}, not {.obj_type_friendly {fits}}."
    ),
    call = call
  )
}

#' Is `truth` a simulation set, or a list of them?
#' @noRd
is_simulation_set_truth <- function(truth) {
  if (inherits(truth, "bmmtools_simulation_set")) {
    return(TRUE)
  }
  is.list(truth) && !is.object(truth) && length(truth) > 0L &&
    all(vapply(truth, inherits, logical(1), what = "bmmtools_simulation_set"))
}

#' Is this one set of fits, a list of sets or an extracted tibble?
#'
#' Used when `truth` is a simulation set, so that a plain named list of fits
#' is read as one set rather than as replications named by the list.
#'
#' @noRd
cor_set_fits_kind <- function(fits, call = rlang::caller_env()) {
  if (is.data.frame(fits)) {
    return("tibble")
  }
  if (is_fit_set(fits, call = call)) {
    return("set")
  }
  is_set_list <- is.list(fits) && !is.object(fits) && length(fits) > 0L &&
    is.null(names(fits)) &&
    all(vapply(fits, function(f) {
      inherits(f, "bmmtools_fit_set") || (is.list(f) && !is.object(f))
    }, logical(1)))
  if (is_set_list && all(vapply(fits, is_fit_set, logical(1), call = call))) {
    return("sets")
  }
  cli::cli_abort(
    c(
      "With a simulation set as {.arg truth}, {.arg fits} must be a set of \\
       fits, a list of sets, or a tibble from {.fn extract_correlations}, \\
       not {.obj_type_friendly {fits}}.",
      i = "A set of fits comes from {.fn fit_components}, or is a list of \\
           fits named by component."
    ),
    call = call
  )
}

#' Check that each set of fits has the simulation set's components
#' @noRd
check_set_components <- function(sim_set, fit_sets, call) {
  expected <- names(sim_set$components)
  for (fits in fit_sets) {
    if (!setequal(names(fits), expected)) {
      cli::cli_abort(
        c(
          "The fits must be named by the simulation set's components.",
          i = "Fits: {.val {names(fits)}}; components: {.val {expected}}."
        ),
        call = call
      )
    }
  }
  invisible(NULL)
}

#' @noRd
is_one_fit <- function(x) {
  set_classes <- c(
    "bmmtools_simulation", "bmmtools_simulation_set", "bmmtools_fit_set",
    "bmmtools_component"
  )
  is.object(x) && !is.data.frame(x) && !inherits(x, set_classes)
}

#' Validate an extracted correlation tibble
#' @noRd
check_extracted_correlations <- function(fits, call = rlang::caller_env()) {
  required <- setdiff(correlation_columns(), "converged")
  missing <- setdiff(required, names(fits))
  if (length(missing) > 0L) {
    cli::cli_abort(
      c(
        "{.arg fits} is missing the column{?s} {.val {missing}}.",
        i = "A correlation tibble comes from {.fn extract_correlations}."
      ),
      call = call
    )
  }
  out <- tibble::as_tibble(fits)
  if (!"replication" %in% names(out)) out$replication <- rep(1L, nrow(out))
  if (!"converged" %in% names(out)) out$converged <- rep(NA, nrow(out))
  out
}

#' @noRd
tibble_cor_estimates <- function(fits, estimator, scale, pairs, call) {
  out <- fits[fits$estimator %in% estimator, ]
  if (!is.null(pairs)) out <- out[out$term %in% pairs, ]
  if (nrow(out) == 0L) {
    cli::cli_abort(
      "No correlation estimates for the estimator{?s} {.val {estimator}} \\
       to score.",
      call = call
    )
  }
  # nolint next: object_usage_linter. Used by cli's glue interpolation.
  other <- setdiff(unique(out$scale), scale)
  if (length(other) > 0L) {
    cli::cli_abort(
      c(
        "{.arg fits} was extracted on the {.val {other}} scale, not the \\
         {.val {scale}} scale.",
        i = "Extract again with {.code scale = \"{scale}\"}, or score with \\
             {.code scale = \"{other[[1L]]}\"}."
      ),
      call = call
    )
  }
  out
}

#' Extract correlations from one fit or a list of fits
#' @noRd
fit_cor_estimates <- function(fits, kind, truths, estimator, resolved,
                              pairs, group, ci_level, dots, call) {
  one <- function(fit, label) {
    covariates <- covariates_for(truths$covariates, label, group)
    out <- rlang::exec(
      extract_correlations, fit,
      estimator = estimator, pairs = pairs, covariates = covariates,
      group = group, scale = resolved$scale, links = resolved$links,
      ci_level = ci_level, !!!dots
    )
    out$replication <- rep(label, nrow(out))
    out
  }
  if (kind %in% c("fit", "set")) {
    return(one(fits, 1L))
  }
  labels <- names(fits) %||% seq_along(fits)
  dplyr::bind_rows(lapply(seq_along(fits), function(i) {
    one(fits[[i]], labels[[i]])
  }))
}

#' The covariate table of one replication, wide, from the long truth
#' @noRd
covariates_for <- function(covariates, label, group) {
  if (is.null(covariates) || nrow(covariates) == 0L) {
    return(NULL)
  }
  if ("replication" %in% names(covariates)) {
    covariates <- covariates[covariates$replication == label, ]
  }
  if (nrow(covariates) == 0L) {
    return(NULL)
  }
  ids <- unique(as.character(covariates$id))
  out <- data.frame(ids, stringsAsFactors = FALSE)
  names(out) <- group
  for (term in unique(covariates$term)) {
    rows <- covariates[covariates$term == term, ]
    out[[term]] <- as.double(rows$true_value[match(ids, rows$id)])
  }
  out
}

#' Normalise `truth` into correlation, subject and covariate tables
#'
#' @return A list with `cor`, `subjects`, `covariates` and `links` (the
#'   simulation model's link table, `NULL` when there is none).
#' @noRd
cor_truth_input <- function(truth, fits, kind, call = rlang::caller_env()) {
  empty <- tibble::tibble(
    id = character(), term = character(), true_value = double()
  )
  sim_parts <- function(sim) {
    tables <- sim$truth
    if (is.null(tables$cor) || is.null(tables$subjects)) {
      cli::cli_abort(
        c(
          "The simulation in {.arg truth} has no correlation truth.",
          i = "It was made before correlated truths existed; simulate it \\
               again."
        ),
        call = call
      )
    }
    list(
      cor = tables$cor, subjects = tables$subjects,
      covariates = tables$covariates %||% empty
    )
  }

  # a simulation set carries prefixed truth tables and a prefixed link
  # table; otherwise it is read as a simulation is
  links_of <- function(sim) {
    if (inherits(sim, "bmmtools_simulation_set")) sim$links else sim$model$links
  }
  fit_sets <- switch(kind,
    set = list(fits),
    sets = fits,
    list()
  )

  if (inherits(truth, c("bmmtools_simulation", "bmmtools_simulation_set"))) {
    if (inherits(truth, "bmmtools_simulation_set")) {
      check_set_components(truth, fit_sets, call)
    }
    parts <- sim_parts(truth)
    parts$links <- links_of(truth)
    return(normalise_cor_truth(parts, call))
  }

  sim_class <- "bmmtools_simulation"
  if (is_simulation_set_truth(truth)) sim_class <- "bmmtools_simulation_set"
  is_sim_list <- is.list(truth) && !is.object(truth) && length(truth) > 0L &&
    all(vapply(truth, inherits, logical(1), what = sim_class))
  if (is_sim_list) {
    n_fits <- switch(kind,
      fit = ,
      set = 1L,
      list = ,
      sets = length(fits),
      NA_integer_
    )
    if (!is.na(n_fits) && length(truth) != n_fits) {
      cli::cli_abort(
        "{.arg truth} has {length(truth)} simulation{?s} for \\
         {n_fits} fit{?s}; give one per fit.",
        call = call
      )
    }
    labels <- if (kind %in% c("list", "sets")) {
      names(fits) %||% seq_along(fits)
    } else {
      names(truth) %||% seq_along(truth)
    }
    if (identical(sim_class, "bmmtools_simulation_set")) {
      for (i in seq_along(truth)) {
        check_set_components(
          truth[[i]],
          if (identical(kind, "sets")) fit_sets[i] else fit_sets,
          call
        )
      }
    }
    pieces <- lapply(seq_along(truth), function(i) {
      parts <- sim_parts(truth[[i]])
      lapply(parts, function(table) {
        table$replication <- rep(labels[[i]], nrow(table))
        table
      })
    })
    bound <- lapply(c("cor", "subjects", "covariates"), function(name) {
      dplyr::bind_rows(lapply(pieces, `[[`, name))
    })
    parts <- stats::setNames(bound, c("cor", "subjects", "covariates"))
    parts$links <- links_of(truth[[1L]])
    return(normalise_cor_truth(parts, call))
  }

  is_table_list <- is.list(truth) && !is.data.frame(truth) &&
    !is.object(truth) && all(c("cor", "subjects") %in% names(truth))
  if (is_table_list) {
    parts <- list(
      cor = truth$cor, subjects = truth$subjects,
      covariates = truth$covariates %||% empty, links = NULL
    )
    return(normalise_cor_truth(parts, call))
  }

  cli::cli_abort(
    c(
      "{.arg truth} must be a {.cls bmmtools_simulation}, a list of them, \\
       or a list with {.field cor} and {.field subjects} tibbles, not \\
       {.obj_type_friendly {truth}}."
    ),
    call = call
  )
}

#' Check the truth tables and name every correlation pair under D27
#' @noRd
normalise_cor_truth <- function(parts, call) {
  needs <- list(
    cor = c("term", "true_value"),
    subjects = c("id", "term", "true_value"),
    covariates = c("id", "term", "true_value")
  )
  for (name in names(needs)) {
    table <- parts[[name]]
    if (!is.data.frame(table)) {
      cli::cli_abort(
        "{.field {name}} in {.arg truth} must be a data frame, \\
         not {.obj_type_friendly {table}}.",
        call = call
      )
    }
    missing <- setdiff(needs[[name]], names(table))
    if (length(missing) > 0L) {
      cli::cli_abort(
        "{.field {name}} in {.arg truth} is missing the column{?s} \\
         {.val {missing}}.",
        call = call
      )
    }
  }

  cor <- tibble::as_tibble(parts$cor)
  split_terms <- strsplit(as.character(cor$term), "__", fixed = TRUE)
  first <- vapply(split_terms, function(p) p[[1L]], character(1))
  second <- vapply(split_terms, function(p) {
    if (length(p) >= 2L) p[[2L]] else NA_character_
  }, character(1))
  if (!"var1" %in% names(cor)) cor$var1 <- first
  if (!"var2" %in% names(cor)) cor$var2 <- second
  named <- lapply(seq_len(nrow(cor)), function(k) {
    pair_term(cor$var1[[k]], cor$var2[[k]])
  })
  cor$term <- vapply(named, `[[`, character(1), "term")
  cor$var1 <- vapply(named, `[[`, character(1), "var1")
  cor$var2 <- vapply(named, `[[`, character(1), "var2")
  cor$true_value <- as.double(cor$true_value)

  subjects <- tibble::as_tibble(parts$subjects)
  subjects$id <- as.character(subjects$id)
  covariates <- tibble::as_tibble(parts$covariates)
  covariates$id <- as.character(covariates$id)
  list(
    cor = cor, subjects = subjects, covariates = covariates,
    links = parts$links
  )
}

#' The in-sample correlation of the true subject values, per replication
#'
#' Grouped by the `replication` and `condition` columns the subject table
#' carries; a covariate table without them is shared by every group.
#'
#' @return A tibble with those key columns, `term` and `sample_value`.
#' @noRd
sample_correlations <- function(cor, subjects, covariates, scale, links) {
  pairs <- cor[!duplicated(cor$term), c("term", "var1", "var2")]
  keys <- intersect(c("replication", "condition"), names(subjects))
  groups <- if (length(keys) == 0L) {
    list(seq_len(nrow(subjects)))
  } else {
    split(
      seq_len(nrow(subjects)),
      do.call(paste, c(unname(as.list(subjects[keys])), sep = "\r")),
      drop = TRUE
    )
  }
  cov_keys <- intersect(keys, names(covariates))

  pieces <- lapply(groups, function(rows) {
    sub <- subjects[rows, ]
    cov <- covariates
    for (key in cov_keys) {
      cov <- cov[cov[[key]] == sub[[key]][[1L]], ]
    }
    ids <- unique(c(sub$id, cov$id))
    value_of <- function(term) {
      if (term %in% cov$term) {
        at <- cov[cov$term == term, ]
        return(as.double(at$true_value[match(ids, at$id)]))
      }
      at <- sub[sub$term == term, ]
      values <- as.double(at$true_value[match(ids, at$id)])
      if (identical(scale, "natural")) {
        values <- inverse_link(values, link_of(term, links))
      }
      values
    }
    sample_value <- vapply(seq_len(nrow(pairs)), function(k) {
      metric_r(value_of(pairs$var1[[k]]), value_of(pairs$var2[[k]]))$r
    }, double(1))
    out <- tibble::tibble(term = pairs$term, sample_value = sample_value)
    for (key in rev(keys)) {
      out <- tibble::add_column(
        out, !!key := rep(sub[[key]][[1L]], nrow(out)),
        .before = 1L
      )
    }
    out
  })
  dplyr::bind_rows(pieces)
}
