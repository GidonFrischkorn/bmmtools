# The estimates tibble: the single contract between the run layer and the
# score layer.
#
# extract_estimates() is the only score-layer function that touches a fit
# object. Everything downstream of it --- recover(), the metrics, the
# summary and the plot --- works on tibbles, which is what makes the
# layer testable and installable without brms.

#' The estimates tibble contract
#'
#' Column names and their `typeof()`, in contract order. The first nine
#' are the apabayes `parameters` columns; `level` and `id` are the
#' bmmtools additions.
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
    id = "character",
    converged = "logical",
    estimator = "character"
  )
}

#' The columns an estimates tibble must bring
#'
#' `converged` is filled by [extract_estimates()] and optional on input,
#' so a hand-built tibble with the other eleven columns still scores.
#' `estimator` is optional for the same reason: a tibble that does not say
#' which estimator produced it is a posterior, which is what `"bayes"`
#' means (milestone 8, decision 35).
#'
#' @noRd
estimates_required_columns <- function() {
  setdiff(names(estimates_contract()), c("converged", "estimator"))
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

#' Split a coefficient name into its parameter and its coefficient
#'
#' `kappa_Intercept` is parameter `kappa`, coefficient `Intercept`;
#' `kappa_task1` is `kappa` and `task1`. A name without `_` has no
#' parameter and is the coefficient itself (`Intercept` in a model without
#' distributional parameters). The split is at the last `_`, which is
#' unambiguous for task terms because neither `task_col` nor a task level
#' may contain one. A `b_` prefix is removed by the caller.
#'
#' @return A list with `par` (`""` where there is none) and `coef`.
#' @noRd
split_coefficient <- function(x) {
  has_par <- grepl("_", x, fixed = TRUE)
  list(
    par = ifelse(has_par, sub("_[^_]+$", "", x), ""),
    coef = ifelse(has_par, sub("^.*_", "", x), x)
  )
}

#' The terms of parameter and coefficient pairs
#'
#' A parameter with one coefficient is its bare name, whatever the
#' coefficient; one with several and no `Intercept`, as cell-means coding
#' (`0 + task`) gives, is `<par>_<coef>` per coefficient. An `Intercept`
#' together with other coefficients is a contrast design (`1 + task`,
#' `coding = "contrast"`): the intercept is the bare name and each other
#' coefficient is `<par>_<coef>`, which is what
#' [contrast_truth()] names the truths those coefficients estimate. A
#' coefficient without a parameter keeps its own name.
#'
#' The two designs cannot collide: cell-means terms carry the task
#' *level* (`kappa_task1` for a task called 1) and contrast terms the
#' contrast *column* (`kappa_task1` for the first contrast), and one fit
#' has only one of them.
#'
#' @param par,coef Parallel character vectors, one entry per coefficient
#'   (repeated entries, one per subject, are allowed).
#' @noRd
coefficient_terms <- function(par, coef, call = rlang::caller_env()) {
  term <- ifelse(nzchar(par), par, coef)
  for (p in unique(par[nzchar(par)])) {
    at <- par == p
    coefs <- unique(coef[at])
    if (length(coefs) < 2L) next
    if ("Intercept" %in% coefs) {
      contrast <- at & coef != "Intercept"
      term[contrast] <- paste0(p, "_", coef[contrast])
      next
    }
    term[at] <- paste0(p, "_", coef[at])
  }
  term
}

#' What each coefficient is: a population value or an effect
#'
#' A coefficient of a parameter that also has an `Intercept` is an effect
#' --- a contrast between tasks, not a value of the parameter --- and is
#' scored as one (`level = "effect"`, always on the link scale). Everything
#' else, the intercept of such a parameter included, is a population value.
#'
#' Read from the coefficient names alone, so an extraction needs nothing
#' from the simulation that produced the fit.
#'
#' @inheritParams coefficient_terms
#' @noRd
coefficient_kinds <- function(par, coef) {
  kind <- rep("population", length(coef))
  for (p in unique(par[nzchar(par)])) {
    at <- par == p
    if (!("Intercept" %in% coef[at])) next
    kind[at & coef != "Intercept"] <- "effect"
  }
  kind
}

#' Terms for brms coefficient names, as a lookup keyed by the name
#' @noRd
name_terms <- function(x, call = rlang::caller_env()) {
  x <- unique(x)
  parts <- split_coefficient(x)
  stats::setNames(coefficient_terms(parts$par, parts$coef, call = call), x)
}

#' Kinds for brms coefficient names, as a lookup keyed by the name
#' @noRd
name_kinds <- function(x) {
  x <- unique(x)
  parts <- split_coefficient(x)
  stats::setNames(coefficient_kinds(parts$par, parts$coef), x)
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
#' The point estimate is the posterior median: it
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
as_estimates <- function(x, level, ci_level, ci_method,
                         estimator = "bayes") {
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
    id = as.character(x$id),
    converged = NA,
    estimator = as.character(estimator)
  )
}

#' Validate the `estimator` argument: one non-missing string
#'
#' The vocabulary is open --- `"bayes"` and `"ml"` are what the package
#' produces, but a user comparing three of their own estimators names them
#' whatever they like. What is refused is a missing or non-scalar label,
#' which would silently split or merge summary rows.
#'
#' @noRd
check_estimator <- function(estimator, call = rlang::caller_env()) {
  bad <- !is.character(estimator) || length(estimator) != 1L ||
    is.na(estimator) || !nzchar(estimator)
  if (bad) {
    cli::cli_abort(
      "{.arg estimator} must be a single non-empty string, \\
       not {.obj_type_friendly {estimator}}.",
      call = call
    )
  }
  estimator
}

#' Validate the `converged` argument: one logical, `NA` allowed
#' @noRd
check_converged <- function(converged, call = rlang::caller_env()) {
  if (!is.logical(converged) || length(converged) != 1L) {
    cli::cli_abort(
      "{.arg converged} must be {.code TRUE}, {.code FALSE} or {.code NA}, \\
       not {.obj_type_friendly {converged}}.",
      call = call
    )
  }
  converged
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
        i = "Two coefficients reduce to the same term, so a join with the \\
             truth would fan out.",
        i = "The offending coefficient{?s}: {.val {offending}}."
      ),
      call = call
    )
  }
  invisible(x)
}

#' Population-level rows: the values of a parameter, its effects, or both
#'
#' Both kinds come from the same `b_` draws (see [coefficient_kinds()]),
#' so one pass produces them and `keep` says which are wanted.
#'
#' @noRd
population_estimates <- function(draws, ci_level, ci_method, drop_constants,
                                 keep = c("population", "effect"),
                                 call = rlang::caller_env()) {
  variables <- grep("^b_", posterior::variables(draws), value = TRUE)
  if (length(variables) == 0L) {
    return(structure(empty_estimates(), n_found = 0L))
  }

  # named before summarising, so a refused design costs no summary
  terms <- name_terms(sub("^b_", "", variables), call = call)
  kinds <- name_kinds(sub("^b_", "", variables))
  # selected before summarising as well, so that asking for one kind
  # neither summarises the other nor counts it as found and then dropped
  variables <- variables[unname(kinds[sub("^b_", "", variables)]) %in% keep]
  if (length(variables) == 0L) {
    return(structure(empty_estimates(), n_found = 0L))
  }
  out <- summarise_selected(
    posterior::subset_draws(draws, variable = variables), ci_level
  )
  out$term <- unname(terms[sub("^b_", "", out$variable)])
  out$id <- NA_character_
  n_found <- nrow(out)
  out <- drop_constant_rows(out, drop_constants)
  check_unique_terms(out, "term", call = call)

  structure(
    as_estimates(
      out, unname(kinds[sub("^b_", "", out$variable)]), ci_level, ci_method
    ),
    n_found = n_found
  )
}

#' Parse the group-level coefficient names of one grouping factor
#'
#' Measured on the mixture2p fixture (bmm 1.4.1.9000, brms 2.23.0):
#' subject draws are named `r_id__kappa[1,Intercept]`, so the parameter
#' sits in the distributional-parameter slot and the level in the
#' bracket. A model without distributional parameters names them
#' `r_id[1,Intercept]`, where the coefficient is the parameter. Cell-means
#' coding gives `r_id__kappa[1,task1]` and the term `kappa_task1` (see
#' `coefficient_terms()`).
#'
#' @return A tibble with `variable`, `term`, `id` and `population`, the
#'   name of the population coefficient the deviation belongs to.
#' @noRd
group_coefficients <- function(variables, group, call = rlang::caller_env()) {
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
    term = coefficient_terms(dpar, coefficient, call = call),
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
#' The order matters: summing the draws and then
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
                              drop_constants, call = rlang::caller_env()) {
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

  out <- summarise_selected(sum_group_draws(draws, coefficients), ci_level)
  out$term <- coefficients$term[match(out$variable, coefficients$variable)]
  out$id <- coefficients$id[match(out$variable, coefficients$variable)]
  n_found <- nrow(out)
  out <- drop_constant_rows(out, drop_constants)
  check_unique_terms(out, c("term", "id"), call = call)

  structure(
    as_estimates(out, "subject", ci_level, ci_method),
    n_found = n_found
  )
}

#' Validate a `ranef` table: `NULL` or a data frame with the columns read
#'
#' `resp`, `dpar` and `nlpar` are optional, because brms fills them with
#' empty strings where they do not apply and a non-brms source need not
#' carry them at all.
#'
#' @noRd
check_ranef <- function(ranef, call = rlang::caller_env()) {
  if (is.null(ranef)) {
    return(NULL)
  }
  if (!is.data.frame(ranef)) {
    cli::cli_abort(
      "{.arg ranef} must be a data frame, such as {.code fit$ranef}, \\
       or {.code NULL}, not {.obj_type_friendly {ranef}}.",
      call = call
    )
  }
  missing <- setdiff(c("id", "group", "coef", "cor"), names(ranef))
  if (length(missing) > 0L) {
    cli::cli_abort(
      "{.arg ranef} is missing the column{?s} {.val {missing}}.",
      call = call
    )
  }
  ranef
}

#' A column of a `ranef` table as character, empty strings when absent
#' @noRd
ranef_column <- function(ranef, column) {
  if (!column %in% names(ranef)) {
    return(rep("", nrow(ranef)))
  }
  out <- as.character(ranef[[column]])
  out[is.na(out)] <- ""
  out
}

#' The group-level coefficients of one grouping factor, from `ranef`
#'
#' brms names a coefficient inside `sd_` and `cor_` by its prefix ---
#' `resp`, then `dpar` (unless it is `mu`), then `nlpar`, each followed by
#' `_` when present --- and the coefficient: `kappa_Intercept` for a bmm
#' parameter. Building the names from the table rather than parsing them
#' out of the draws is what keeps a parameter name that itself contains
#' `_` intact.
#'
#' The term is the parameter (`nlpar`, else `dpar`), or the coefficient
#' when there is neither, as in `group_coefficients()`; a parameter with
#' several coefficients gets one term per coefficient through
#' `coefficient_terms()`. `resp` is not part of the term until Milestone 8
#' prepends it.
#'
#' Columns are read as vectors, never by subsetting the table: brms gives
#' it a class of its own whose `[` method needs brms loaded.
#'
#' @return A tibble with `variable` (the brms name without its `sd_`
#'   prefix), `term`, `cor` and `block` (the correlation block id).
#' @noRd
ranef_coefficients <- function(ranef, group, call = rlang::caller_env()) {
  rows <- as.character(ranef$group) == group
  resp <- ranef_column(ranef, "resp")[rows]
  dpar <- ranef_column(ranef, "dpar")[rows]
  dpar[dpar == "mu"] <- ""
  nlpar <- ranef_column(ranef, "nlpar")[rows]
  coef <- as.character(ranef$coef)[rows]

  prefix <- paste0(
    ifelse(nzchar(resp), paste0(resp, "_"), ""),
    ifelse(nzchar(dpar), paste0(dpar, "_"), ""),
    ifelse(nzchar(nlpar), paste0(nlpar, "_"), "")
  )
  par <- ifelse(nzchar(nlpar), nlpar, dpar)
  tibble::tibble(
    variable = paste0(prefix, coef),
    term = coefficient_terms(par, coef, call = call),
    cor = as.logical(ranef$cor)[rows] %in% TRUE,
    block = as.character(ranef$id)[rows]
  )
}

#' Terms for the coefficient names inside `sd_` and `cor_`, parsed
#'
#' Only for draws without a `ranef` table. The names of both kinds of
#' variable are pooled, so that a coefficient is named the same way at the
#' `sd` and the `cor` level, whichever draws a hand-built object carries.
#'
#' @return A named character vector: coefficient name to term.
#' @noRd
parsed_group_terms <- function(variables, group, call = rlang::caller_env()) {
  sd_prefix <- paste0("sd_", group, "__")
  cor_prefix <- paste0("cor_", group, "__")
  sds <- variables[startsWith(variables, sd_prefix)]
  cors <- variables[startsWith(variables, cor_prefix)]
  coefficients <- c(
    substring(sds, nchar(sd_prefix) + 1L),
    unlist(strsplit(substring(cors, nchar(cor_prefix) + 1L), "__",
      fixed = TRUE
    ))
  )
  name_terms(coefficients, call = call)
}

#' The `sd_` variables of one grouping factor and their terms
#'
#' With `ranef`, the names are built and looked up; without, every
#' `sd_<group>__` variable is parsed. brms reserves `__` in variable
#' names, so the prefix cannot also match another group's.
#'
#' @return A tibble with `variable` and `term`.
#' @noRd
sd_variables <- function(variables, group, ranef, call = rlang::caller_env()) {
  prefix <- paste0("sd_", group, "__")
  if (is.null(ranef)) {
    found <- variables[startsWith(variables, prefix)]
    rest <- substring(found, nchar(prefix) + 1L)
    terms <- parsed_group_terms(variables, group, call = call)
    return(tibble::tibble(variable = found, term = unname(terms[rest])))
  }
  coefficients <- ranef_coefficients(ranef, group, call = call)
  out <- tibble::tibble(
    variable = paste0(prefix, coefficients$variable),
    term = coefficients$term
  )
  out[out$variable %in% variables, ]
}

#' The `cor_` variables of one grouping factor and their pair terms
#'
#' With `ranef`, every pair of coefficients that brms correlates (`cor`
#' is `TRUE` and the two share a block id) is looked up under both orders
#' of its two names, because brms's order is not guaranteed. A pair
#' neither order finds is not estimated and gets no row.
#'
#' @return A tibble with `variable`, `term`, `var1` and `var2`.
#' @noRd
cor_variables <- function(variables, group, ranef,
                          call = rlang::caller_env()) {
  prefix <- paste0("cor_", group, "__")
  empty <- tibble::tibble(
    variable = character(), term = character(),
    var1 = character(), var2 = character()
  )

  if (is.null(ranef)) {
    found <- variables[startsWith(variables, prefix)]
    parts <- strsplit(substring(found, nchar(prefix) + 1L), "__", fixed = TRUE)
    terms <- parsed_group_terms(variables, group, call = call)
    pairs <- lapply(seq_along(found)[lengths(parts) == 2L], function(k) {
      names <- unname(terms[parts[[k]]])
      c(variable = found[[k]], pair_term(names[[1L]], names[[2L]]))
    })
    return(dplyr::bind_rows(empty, pairs))
  }

  coefficients <- ranef_coefficients(ranef, group, call = call)
  correlated <- coefficients[coefficients$cor, ]
  if (nrow(correlated) < 2L) {
    return(empty)
  }
  combos <- utils::combn(nrow(correlated), 2L)
  pairs <- lapply(seq_len(ncol(combos)), function(k) {
    a <- correlated[combos[1L, k], ]
    b <- correlated[combos[2L, k], ]
    if (!identical(a$block, b$block)) {
      return(NULL)
    }
    candidates <- paste0(
      prefix,
      c(
        paste0(a$variable, "__", b$variable),
        paste0(b$variable, "__", a$variable)
      )
    )
    hit <- candidates[candidates %in% variables]
    if (length(hit) == 0L) {
      return(NULL)
    }
    c(variable = hit[[1L]], pair_term(a$term, b$term))
  })
  dplyr::bind_rows(empty, pairs)
}

#' Standard-deviation and correlation rows of one grouping factor
#'
#' One body for both levels, since they differ only in which variables
#' they select. The terms are checked on all of the group's coefficients
#' first, so that an intercept with contrasts, or two coefficients reduced
#' to one term, is refused at the cor level too.
#'
#' @noRd
group_level_estimates <- function(draws, groups, group, ranef, level,
                                  ci_level, ci_method, drop_constants,
                                  call = rlang::caller_env()) {
  group <- resolve_group(group, groups, call = call)
  variables <- posterior::variables(draws)

  if (is.null(ranef)) {
    check_unique_terms(sd_variables(variables, group, NULL, call = call),
      "term",
      call = call
    )
  } else {
    check_unique_terms(ranef_coefficients(ranef, group, call = call), "term",
      call = call
    )
  }

  selected <- if (identical(level, "sd")) {
    sd_variables(variables, group, ranef, call = call)
  } else {
    cor_variables(variables, group, ranef, call = call)
  }
  if (nrow(selected) == 0L) {
    if (identical(level, "sd")) {
      cli::cli_abort(
        "The fit has no group-level standard deviations for {.val {group}}.",
        call = call
      )
    }
    return(structure(empty_estimates(), n_found = 0L))
  }

  out <- summarise_selected(
    posterior::subset_draws(draws, variable = selected$variable), ci_level
  )
  out$term <- selected$term[match(out$variable, selected$variable)]
  out$id <- NA_character_
  n_found <- nrow(out)
  out <- drop_constant_rows(out, drop_constants)

  structure(
    as_estimates(out, level, ci_level, ci_method),
    n_found = n_found
  )
}

#' Turn a draws object into the estimates tibble
#'
#' The workhorse behind [extract_estimates()]. It takes a draws object
#' and the fit's grouping factors rather than the fit itself, so that
#' structures a saved fixture does not contain --- two grouping factors,
#' a chain with missing draws --- can be tested without running Stan.
#'
#' `ranef` is the fit's `ranef` table. The `"sd"` and `"cor"` levels use
#' it to build brms's variable names; without it they parse the names,
#' which is enough for hand-built draws.
#'
#' @noRd
estimates_from_draws <- function(draws,
                                 groups,
                                 level = "population",
                                 group = NULL,
                                 ci_level = 0.95,
                                 ci_method = "eti",
                                 drop_constants = TRUE,
                                 converged = NA,
                                 ranef = NULL,
                                 estimator = "bayes",
                                 call = rlang::caller_env()) {
  converged <- check_converged(converged, call = call)
  estimator <- check_estimator(estimator, call = call)
  level <- rlang::arg_match(
    level, c("population", "subject", "effect", "sd", "cor"),
    multiple = TRUE,
    error_call = call
  )
  ranef <- check_ranef(ranef, call = call)
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
  if (!identical(ci_method, "eti")) {
    cli::cli_abort(
      c(
        "{.arg ci_method} must be {.val eti} in this version of bmmtools.",
        i = "The column is carried so that the apabayes contract stays \\
             satisfied when other interval types arrive."
      ),
      call = call
    )
  }
  if (!rlang::is_bool(drop_constants)) {
    cli::cli_abort(
      "{.arg drop_constants} must be {.code TRUE} or {.code FALSE}.",
      call = call
    )
  }

  pieces <- list()
  coefficients <- intersect(c("population", "effect"), level)
  if (length(coefficients) > 0L) {
    pieces$population <- population_estimates(
      draws, ci_level, ci_method, drop_constants,
      keep = coefficients, call = call
    )
  }
  if ("subject" %in% level) {
    pieces$subject <- subject_estimates(
      draws, groups, group, ci_level, ci_method, drop_constants,
      call = call
    )
  }
  for (grouped in intersect(c("sd", "cor"), level)) {
    pieces[[grouped]] <- group_level_estimates(
      draws, groups, group, ranef, grouped,
      ci_level, ci_method, drop_constants,
      call = call
    )
  }
  # Counted before constants are dropped. Testing the draws for any
  # `^(b|sd|cor)_` name instead would report a level that found nothing
  # to extract (no cor_ draws in an uncorrelated fit) as one whose every
  # parameter was dropped.
  n_found <- sum(vapply(pieces, function(p) attr(p, "n_found"), integer(1)))
  out <- dplyr::bind_rows(pieces)
  attr(out, "n_found") <- NULL

  if (nrow(out) == 0L && n_found > 0L && drop_constants) {
    cli::cli_warn(c(
      "Every parameter was dropped as a constant.",
      i = "Use {.code drop_constants = FALSE} to see them with their \\
           missing diagnostics."
    ))
    return(empty_estimates())
  }
  out$converged <- rep(converged, nrow(out))
  out$estimator <- rep(estimator, nrow(out))
  out
}

#' The convergence verdict of a fit under the default thresholds
#'
#' What `converged = NULL` means in [extract_estimates()] and
#' [extract_correlations()], in one place.
#'
#' @noRd
fit_converged <- function(fit, draws) {
  convergence_from_fit(
    fit, draws,
    thresholds = check_thresholds(),
    treedepth_max = fit_treedepth_max(fit)
  )$pass
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
#'   `"effect"` (the contrasts of a `1 + task` design, see
#'   [simulate_recovery()]'s `coding`), `"sd"` (the group-level standard
#'   deviations), `"cor"` (the group-level correlations), or several, in
#'   which case they are stacked and the `level` column separates them.
#'
#'   `"population"` and `"effect"` split the same `b_` coefficients: a
#'   parameter whose coefficients include an `Intercept` has that
#'   intercept as its population value and every other coefficient as an
#'   effect. A fit with no such parameter has no `"effect"` rows.
#' @param group The grouping factor subject-level, SD and correlation
#'   estimates come from. `NULL` uses the fit's only grouping factor and
#'   errors if there is more than one.
#' @param ci_level The interval mass, a single number strictly between 0
#'   and 1.
#' @param ci_method The interval type. Only `"eti"`, the equal-tailed
#'   interval, is available in this version; the column is carried so
#'   that adding others later is not a breaking change.
#' @param drop_constants Drop the parameters the model fixed to
#'   constants. A constant is identified by zero posterior variance, not
#'   by a missing rhat, so that a chain that broke is never silently
#'   dropped as though it had been fixed on purpose.
#' @param converged Whether the fit passed the convergence gate. `NULL`,
#'   the default, computes it as `check_convergence(fit)$pass` with the
#'   default thresholds; a logical scalar is used as given, so a verdict
#'   from [check_convergence()] with other thresholds can be passed in;
#'   `NA` marks it unknown.
#' @param estimator A label for how the estimates were produced, carried
#'   into the `estimator` column and used by [recover()] to keep two
#'   estimators of the same parameter apart in `summary()` and in
#'   [plot_recovery()]. `"bayes"`, the default, is the posterior of a
#'   sampled fit. Pass another label --- `"ml"` for a maximum-likelihood
#'   fit, or any name of your own --- when scoring several estimators
#'   against one truth. Without it, rows from two estimators would be
#'   pooled into a single bias and RMSE.
#' @param ... Not used. Present so the generic can gain arguments later;
#'   anything passed is an error.
#'
#' @return A tibble with the columns `term`, `estimate`, `ci_low`,
#'   `ci_high`, `ci_method`, `ci_level`, `rhat`, `ess_bulk`, `ess_tail`,
#'   `level`, `id`, `converged` and `estimator`, in that order. `id` is
#'   `NA` except on subject rows; `converged` and `estimator` are the same
#'   value on every row of a fit.
#'
#' @details
#' Subject-level estimates are the **per-draw sum** of the population
#' intercept and the group-level deviation, summarised afterwards. That
#' puts them on the same scale as the population rows, so a per-subject
#' truth can be compared against them directly; taking the sum after
#' summarising instead would give intervals that are too narrow.
#'
#' Standard deviations and correlations are returned on the link scale,
#' with `id` set to `NA`. An SD row carries the parameter's name
#' (`kappa`); a correlation row carries the two names joined by `__` and
#' sorted in the C locale (`kappa__thetat`), whichever order brms used.
#' A correlation the model does not estimate has no row.
#'
#' **Terms.** A parameter with one coefficient carries its own name. One
#' with several coefficients and no intercept, as cell-means coding
#' (`kappa ~ 0 + task`) gives, has a row per coefficient, named
#' `<parameter>_<coefficient>` (`kappa_task1`) at every level. An intercept
#' together with other coefficients (population effects or contrasts) is an
#' error in this version.
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
                                      level = c(
                                        "population", "subject", "effect",
                                        "sd", "cor"
                                      ),
                                      group = NULL,
                                      ci_level = 0.95,
                                      ci_method = "eti",
                                      drop_constants = TRUE,
                                      converged = NULL,
                                      estimator = "bayes",
                                      ...) {
  rlang::check_installed("brms", "to extract estimates from a fit.")
  rlang::check_dots_empty()
  if (missing(level)) level <- "population"
  draws <- posterior::as_draws_array(fit)
  if (is.null(converged)) {
    converged <- fit_converged(fit, draws)
  }

  estimates_from_draws(
    draws = draws,
    groups = fit_groups(fit),
    level = level,
    group = group,
    ci_level = ci_level,
    ci_method = ci_method,
    drop_constants = drop_constants,
    converged = converged,
    estimator = estimator,
    ranef = fit$ranef,
    # so a bad argument is reported against extract_estimates(), not
    # against the internal helper that happened to inspect it
    call = rlang::current_env()
  )
}
