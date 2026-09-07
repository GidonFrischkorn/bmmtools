# The estimates tibble: the single contract between the run layer and the
# score layer (ARCHITECTURE.md decision 1, spec section 2).
#
# extract_estimates() is the only score-layer function that touches a fit
# object. Everything downstream of it --- recover(), the metrics, the
# summary and the plot --- works on tibbles, which is what makes the
# layer testable and installable without brms.

#' The estimates tibble contract
#'
#' Column names and their `typeof()`, in contract order. The first nine
#' are the apabayes `parameters` columns; `level` and `id` are the
#' bmmtools additions (ARCHITECTURE.md decision 14).
#'
#' @return A named character vector: column name to storage type.
#' @noRd
estimates_contract <- function() {
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
    level = "character",
    id = "character"
  )
}

#' An empty estimates tibble
#'
#' Used when every term was dropped, so that a caller gets the contract
#' columns rather than a zero-column tibble it cannot index.
#'
#' @noRd
empty_estimates <- function() {
  contract <- estimates_contract()
  out <- lapply(contract, function(type) vector(type, 0L))
  tibble::as_tibble(out)
}

#' Reduce a coefficient name to its bare parameter name
#'
#' `b_kappa_Intercept` becomes `kappa`; so does `b_kappa_setsize2`.
#' Milestone 1 scores intercept-only fits (the `~ 1 + (1 | id)` default
#' of decision 5), so a design suffix is dropped rather than recorded ---
#' and because two design coefficients then collapse onto one name, the
#' caller's duplicate check turns such a fit into an error instead of a
#' silent join fan-out.
#'
#' @noRd
strip_design_suffix <- function(x) {
  x <- sub("^b_", "", x)
  sub("_[^_]+$", "", x)
}

#' The grouping factors a fit has
#'
#' Read from the fit's own `ranef` table (measured on brms 2.23.0: a data
#' frame with a `group` column, one row per group-level coefficient).
#' A fit without group-level effects gives `character(0)`, which
#' `resolve_group()` turns into a message a user can act on.
#'
#' @noRd
fit_groups <- function(fit) {
  re <- fit$ranef
  if (is.data.frame(re) && nrow(re) > 0L && !is.null(re$group)) {
    return(unique(as.character(re$group)))
  }
  character(0)
}

#' Decide which grouping factor subject-level estimates come from
#' @noRd
resolve_group <- function(group, groups, call = rlang::caller_env()) {
  if (is.null(group)) {
    if (length(groups) == 0L) {
      cli::cli_abort(
        c(
          "The fit has no group-level effects, so there are no \\
           subject-level estimates to extract.",
          i = "Use {.code level = \"population\"}."
        ),
        call = call
      )
    }
    if (length(groups) > 1L) {
      cli::cli_abort(
        c(
          "The fit has more than one grouping factor: {.val {groups}}.",
          i = "Supply {.arg group} to say which one to use."
        ),
        call = call
      )
    }
    return(groups)
  }
  if (!is.character(group) || length(group) != 1L || is.na(group)) {
    cli::cli_abort(
      "{.arg group} must be a single string, \\
       not {.obj_type_friendly {group}}.",
      call = call
    )
  }
  if (!group %in% groups) {
    cli::cli_abort(
      "{.arg group} must be a grouping factor of the fit, \\
       one of {.val {groups}}, not {.val {group}}.",
      call = call
    )
  }
  group
}

#' Summarise a draws array into the numeric half of the contract
#'
#' The point estimate is the posterior median (spec, API decision 1): it
#' is invariant under the link transform, so scoring on the link scale
#' and scoring on the natural scale use the same number.
#'
#' `post_var` is carried alongside because a parameter bmm fixed to a
#' constant is identified by zero posterior variance, never by a missing
#' rhat --- a chain that broke also has a missing rhat and must stay
#' visible. Draws that are `NA` are dropped from the point summary and
#' left to show up as a missing rhat.
#'
#' @noRd
summarise_selected <- function(draws, ci_level) {
  probs <- c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2)
  dm <- posterior::as_draws_matrix(draws)
  conv <- posterior::summarise_draws(
    draws, posterior::default_convergence_measures()
  )
  conv <- conv[match(colnames(dm), conv$variable), ]

  tibble::tibble(
    variable = colnames(dm),
    estimate = apply(dm, 2L, stats::median, na.rm = TRUE),
    ci_low = apply(
      dm, 2L, stats::quantile,
      probs = probs[1], names = FALSE, na.rm = TRUE
    ),
    ci_high = apply(
      dm, 2L, stats::quantile,
      probs = probs[2], names = FALSE, na.rm = TRUE
    ),
    rhat = as.double(conv$rhat),
    ess_bulk = as.double(conv$ess_bulk),
    ess_tail = as.double(conv$ess_tail),
    post_var = apply(dm, 2L, stats::var, na.rm = TRUE)
  )
}

#' Turn a summarised table into contract rows
#' @noRd
as_estimates <- function(x, level, ci_level, ci_method) {
  tibble::tibble(
    term = x$term,
    estimate = as.double(x$estimate),
    ci_low = as.double(x$ci_low),
    ci_high = as.double(x$ci_high),
    ci_method = ci_method,
    ci_level = as.double(ci_level),
    rhat = x$rhat,
    ess_bulk = x$ess_bulk,
    ess_tail = x$ess_tail,
    level = level,
    id = as.character(x$id)
  )
}

#' Drop the parameters a model fixed to constants
#'
#' A constant has exactly zero posterior variance. `NA` variance --- a
#' chain with missing draws --- is not zero and is kept.
#'
#' @noRd
drop_constant_rows <- function(x, drop_constants) {
  if (!drop_constants) {
    return(x)
  }
  x[!vapply(x$post_var, function(v) isTRUE(v == 0), logical(1)), ]
}

#' Error when terms collapse onto the same name
#' @noRd
check_unique_terms <- function(x, keys, call = rlang::caller_env()) {
  id <- do.call(paste, c(unname(as.list(x[keys])), sep = "\r"))
  is_duplicate <- duplicated(id) | duplicated(id, fromLast = TRUE)
  duplicated_terms <- unique(x$term[is_duplicate])
  # nolint next: object_usage_linter. Used by cli's glue interpolation.
  offending <- x$variable[is_duplicate]
  if (length(duplicated_terms) > 0L) {
    cli::cli_abort(
      c(
        "Parameter name{?s} {.val {duplicated_terms}} \\
         {?is/are} not unique in the fit.",
        i = "Milestone 1 scores intercept-only fits; a term with design \\
             structure reduces to the same bare name as its intercept.",
        i = "The offending coefficient{?s}: {.val {offending}}."
      ),
      call = call
    )
  }
  invisible(x)
}

#' Population-level rows
#' @noRd
population_estimates <- function(draws, ci_level, ci_method, drop_constants) {
  variables <- grep("^b_", posterior::variables(draws), value = TRUE)
  if (length(variables) == 0L) {
    return(empty_estimates())
  }

  out <- summarise_selected(
    posterior::subset_draws(draws, variable = variables), ci_level
  )
  out$term <- strip_design_suffix(out$variable)
  out$id <- NA_character_
  out <- drop_constant_rows(out, drop_constants)
  check_unique_terms(out, "term")

  as_estimates(out, "population", ci_level, ci_method)
}

#' Parse the group-level coefficient names of one grouping factor
#'
#' Measured on the mixture2p fixture (bmm 1.4.1.9000, brms 2.23.0):
#' subject draws are named `r_id__kappa[1,Intercept]`, so the parameter
#' sits in the distributional-parameter slot and the level in the
#' bracket. A model without distributional parameters names them
#' `r_id[1,Intercept]`, where the coefficient is the parameter.
#'
#' @return A tibble with `variable`, `term`, `id` and `population`, the
#'   name of the population coefficient the deviation belongs to.
#' @noRd
group_coefficients <- function(variables, group) {
  prefix <- paste0("r_", group)
  keep <- startsWith(variables, paste0(prefix, "__")) |
    startsWith(variables, paste0(prefix, "["))
  variables <- variables[keep]
  if (length(variables) == 0L) {
    return(NULL)
  }

  rest <- substring(variables, nchar(prefix) + 1L)
  dpar <- sub("^__", "", sub("\\[.*$", "", rest))
  inside <- sub("^.*\\[(.*)\\]$", "\\1", rest)
  parts <- strsplit(inside, ",", fixed = TRUE)
  id <- vapply(parts, function(p) p[[1L]], character(1))
  coefficient <- vapply(
    parts, function(p) paste(p[-1L], collapse = ","), character(1)
  )

  has_dpar <- nzchar(dpar)
  tibble::tibble(
    variable = variables,
    term = ifelse(has_dpar, dpar, coefficient),
    id = id,
    population = ifelse(
      has_dpar,
      paste0("b_", dpar, "_", coefficient),
      paste0("b_", coefficient)
    )
  )
}

#' Add the population intercept to each group-level deviation, per draw
#'
#' The order matters (spec, API decision 3): summing the draws and then
#' summarising gives the right interval, while summarising each side and
#' adding the two summaries gives one that is too narrow.
#'
#' @noRd
sum_group_draws <- function(draws, coefficients, call = rlang::caller_env()) {
  available <- posterior::variables(draws)
  missing <- setdiff(unique(coefficients$population), available)
  if (length(missing) > 0L) {
    cli::cli_abort(
      c(
        "No population coefficient{?s} {.val {missing}} to add the \\
         group-level deviation{?s} to.",
        i = "Subject-level estimates are the per-draw sum of the \\
             intercept and the deviation."
      ),
      call = call
    )
  }

  plain <- unclass(draws)
  summed <- array(
    NA_real_,
    dim = c(dim(plain)[1L], dim(plain)[2L], nrow(coefficients)),
    dimnames = list(
      iteration = NULL, chain = NULL, variable = coefficients$variable
    )
  )
  for (k in seq_len(nrow(coefficients))) {
    summed[, , k] <- plain[, , coefficients$population[[k]]] +
      plain[, , coefficients$variable[[k]]]
  }
  posterior::as_draws_array(summed)
}

#' Subject-level rows
#' @noRd
subject_estimates <- function(draws, groups, group, ci_level, ci_method,
                              drop_constants) {
  group <- resolve_group(group, groups)
  coefficients <- group_coefficients(posterior::variables(draws), group)
  if (is.null(coefficients)) {
    cli::cli_abort(
      "The fit has no group-level coefficients for {.val {group}}."
    )
  }

  out <- summarise_selected(sum_group_draws(draws, coefficients), ci_level)
  out$term <- coefficients$term[match(out$variable, coefficients$variable)]
  out$id <- coefficients$id[match(out$variable, coefficients$variable)]
  out <- drop_constant_rows(out, drop_constants)
  check_unique_terms(out, c("term", "id"))

  as_estimates(out, "subject", ci_level, ci_method)
}

#' Turn a draws object into the estimates tibble
#'
#' The workhorse behind [extract_estimates()]. It takes a draws object
#' and the fit's grouping factors rather than the fit itself, so that
#' structures a saved fixture does not contain --- two grouping factors,
#' a chain with missing draws --- can be tested without running Stan.
#'
#' @noRd
estimates_from_draws <- function(draws,
                                 groups,
                                 level = "population",
                                 group = NULL,
                                 ci_level = 0.95,
                                 ci_method = "eti",
                                 drop_constants = TRUE) {
  level <- rlang::arg_match(
    level, c("population", "subject"),
    multiple = TRUE
  )
  if (!is.numeric(ci_level) || length(ci_level) != 1L || is.na(ci_level)) {
    cli::cli_abort(
      "{.arg ci_level} must be a single number, \\
       not {.obj_type_friendly {ci_level}}."
    )
  }
  if (ci_level <= 0 || ci_level >= 1) {
    cli::cli_abort("{.arg ci_level} must be between 0 and 1, not {ci_level}.")
  }
  if (!identical(ci_method, "eti")) {
    cli::cli_abort(c(
      "{.arg ci_method} must be {.val eti} in this version of bmmtools.",
      i = "The column is carried so that the apabayes contract stays \\
           satisfied when other interval types arrive."
    ))
  }
  if (!rlang::is_bool(drop_constants)) {
    cli::cli_abort(
      "{.arg drop_constants} must be {.code TRUE} or {.code FALSE}."
    )
  }

  pieces <- list()
  if ("population" %in% level) {
    pieces$population <- population_estimates(
      draws, ci_level, ci_method, drop_constants
    )
  }
  if ("subject" %in% level) {
    pieces$subject <- subject_estimates(
      draws, groups, group, ci_level, ci_method, drop_constants
    )
  }
  out <- dplyr::bind_rows(pieces)

  had_terms <- any(grepl("^b_", posterior::variables(draws)))
  if (nrow(out) == 0L && had_terms && drop_constants) {
    cli::cli_warn(c(
      "Every parameter was dropped as a constant.",
      i = "Use {.code drop_constants = FALSE} to see them with their \\
           missing diagnostics."
    ))
    return(empty_estimates())
  }
  out
}

#' Extract a tidy table of parameter estimates from a fit
#'
#' The bridge between a fitted model and everything bmmtools scores. It
#' returns the *estimates tibble*: one row per parameter, with the
#' posterior median, a credible interval and the convergence diagnostics,
#' under bmm's own parameter names. [recover()] and
#' [recover_subjects()] take either a fit or one of these tibbles, so a
#' scoring pipeline can be built, tested and run with no fitting package
#' installed.
#'
#' The first nine columns are the apabayes `parameters` contract, in its
#' order, so `apabayes::apabayes_tidy(x, type = "parameters")` accepts
#' the result without renaming.
#'
#' @param fit A `brmsfit`, and so also a `bmmfit`.
#' @param level Which estimates to return: `"population"`, `"subject"`,
#'   or both, in which case they are stacked and the `level` column
#'   separates them.
#' @param group The grouping factor subject-level estimates come from.
#'   `NULL` uses the fit's only grouping factor and errors if there is
#'   more than one.
#' @param ci_level The interval mass, a single number strictly between 0
#'   and 1.
#' @param ci_method The interval type. Only `"eti"`, the equal-tailed
#'   interval, is available in this version; the column is carried so
#'   that adding others later is not a breaking change.
#' @param drop_constants Drop the parameters the model fixed to
#'   constants. A constant is identified by zero posterior variance, not
#'   by a missing rhat, so that a chain that broke is never silently
#'   dropped as though it had been fixed on purpose.
#' @param ... Not used. Present so the generic can gain arguments later;
#'   anything passed is an error.
#'
#' @return A tibble with the columns `term`, `estimate`, `ci_low`,
#'   `ci_high`, `ci_method`, `ci_level`, `rhat`, `ess_bulk`, `ess_tail`,
#'   `level` and `id`, in that order. `id` is `NA` for population rows.
#'
#' @details
#' Subject-level estimates are the **per-draw sum** of the population
#' intercept and the group-level deviation, summarised afterwards. That
#' puts them on the same scale as the population rows, so a per-subject
#' truth can be compared against them directly; taking the sum after
#' summarising instead would give intervals that are too narrow.
#'
#' Group-level standard deviations are not returned in this version.
#' `level = "sd"` is reserved for them, so that adding them later does
#' not change what the existing levels mean.
#'
#' @examples
#' \dontrun{
#' fit <- bmm::bmm(bmm::bmf(kappa ~ 1, thetat ~ 1), data, model)
#' extract_estimates(fit)
#' extract_estimates(fit, level = "subject", ci_level = 0.89)
#' }
#'
#' @export
extract_estimates <- function(fit, ...) {
  UseMethod("extract_estimates")
}

#' @rdname extract_estimates
#' @export
extract_estimates.default <- function(fit, ...) {
  cli::cli_abort(
    "{.arg fit} must be a {.cls brmsfit}, not {.obj_type_friendly {fit}}."
  )
}

#' @rdname extract_estimates
#' @export
extract_estimates.brmsfit <- function(fit,
                                      level = c("population", "subject"),
                                      group = NULL,
                                      ci_level = 0.95,
                                      ci_method = "eti",
                                      drop_constants = TRUE,
                                      ...) {
  rlang::check_installed("brms", "to extract estimates from a fit.")
  rlang::check_dots_empty()
  if (missing(level)) level <- "population"

  estimates_from_draws(
    draws = posterior::as_draws_array(fit),
    groups = fit_groups(fit),
    level = level,
    group = group,
    ci_level = ci_level,
    ci_method = ci_method,
    drop_constants = drop_constants
  )
}
