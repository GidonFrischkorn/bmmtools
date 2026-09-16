# Simulation-based calibration: the naming layer (spec 6, section 1;
# local/ARCHITECTURE.md decision 15).
#
# `sbc()` builds two objects and hands them to `SBC::compute_SBC()`. The
# risky half is not the pipeline, it is the names. SBC matches the
# generator's `variables` against the fit's draws matrix; a name that
# matches nothing is not an error there, it is a variable that never
# gets ranked, and the run comes back looking like a short one rather
# than a wrong one. Measured 2026-09-16:
# `posterior::subset_draws(variable = "^zzz$", regex = TRUE)` returns a
# zero-column draws object and says nothing.
#
# So the names are built here, from the model and the formula, and
# asserted with `check_variable_names()`, which errors. Decision 15 said
# the generator's `variables` carry bmm's parameter names; that is wrong
# and was corrected on 2026-09-16 --- they carry the fit's brms draw
# names (`b_kappa_Intercept`, not `kappa`), because those are what SBC
# matches against.

#' Quote the regex metacharacters a parameter or group name may contain
#'
#' bmm's parameter names are syntactic, so they may hold a `.`, which
#' matches any character unquoted: a pattern for `k.a` would also match
#' `kXa` and rank the wrong variable.
#'
#' @noRd
escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

#' The brms draw-name patterns one level selects
#'
#' Over the model's **free** parameters only: a parameter bmm fixes has a
#' draw (`b_mu1_Intercept`) but no prior to rank it against, so the
#' generator has no truth for it and it is never ranked.
#'
#' The four shapes, measured on the committed fixtures:
#' `b_<par>_Intercept`, `sd_<g>__<par>_Intercept`,
#' `cor_<g>__<par>_Intercept__<par>_Intercept` and
#' `r_<g>__<par>[<i>,Intercept]`. A correlation gets a pattern for
#' **both** orders of the pair, because which one brms writes is not
#' guaranteed, and the two share one name because exactly one of them
#' will match.
#'
#' @param group The grouping factor, or `NULL` for a formula with no
#'   group term, which then has no `sd`, `cor` or `subject` patterns.
#' @return A character vector of regexes named `"<level>:<term>"`. The
#'   name says which variable each pattern is there to find, which is
#'   what lets `prior_parameter_draws()` report a pattern that found
#'   nothing instead of quietly returning a narrower matrix.
#' @noRd
sbc_variables <- function(model, group, level, call = rlang::caller_env()) {
  check_model(model, call = call)
  level <- rlang::arg_match(
    level, c("population", "sd", "cor", "subject"),
    multiple = TRUE, error_call = call
  )
  free <- model_parameters(model)$free
  esc <- escape_regex(free)
  out <- character(0)
  add <- function(out, patterns, lvl, terms) {
    c(out, stats::setNames(patterns, paste0(lvl, ":", terms)))
  }

  if ("population" %in% level) {
    out <- add(out, paste0("^b_", esc, "_Intercept$"), "population", free)
  }
  if (is.null(group)) {
    return(out)
  }
  g <- escape_regex(group)

  if ("sd" %in% level) {
    out <- add(out, paste0("^sd_", g, "__", esc, "_Intercept$"), "sd", free)
  }
  if ("cor" %in% level && length(free) > 1L) {
    pairs <- utils::combn(seq_along(free), 2L)
    a <- pairs[1L, ]
    b <- pairs[2L, ]
    terms <- vapply(
      seq_along(a),
      function(k) pair_term(free[[a[[k]]]], free[[b[[k]]]])$term,
      character(1)
    )
    forward <- paste0(
      "^cor_", g, "__", esc[a], "_Intercept__", esc[b], "_Intercept$"
    )
    backward <- paste0(
      "^cor_", g, "__", esc[b], "_Intercept__", esc[a], "_Intercept$"
    )
    out <- add(out, forward, "cor", terms)
    out <- add(out, backward, "cor", terms)
  }
  if ("subject" %in% level) {
    out <- add(
      out, paste0("^r_", g, "__", esc, "\\[.+,Intercept\\]$"),
      "subject", free
    )
  }
  out
}

#' The level each pattern of `sbc_variables()` belongs to
#' @noRd
sbc_variable_levels <- function(variables) {
  sub(":.*$", "", names(variables) %||% character(0))
}

#' Assert that the generated truth names are the names that will be ranked
#'
#' The one check that makes an SBC run mean anything. SBC matches
#' `variables` against the draws matrix by name and ranks what it finds;
#' a mismatch produces a shorter `results$stats`, not a complaint. This
#' is therefore an error, never a warning: a warning in a run of a
#' hundred fits is a line nobody reads, and the ranks it would leave
#' behind look like a result.
#'
#' @param names The names the generator emitted.
#' @param expected The names `sbc()` resolved and will ask SBC to rank.
#' @noRd
check_variable_names <- function(names, expected,
                                 call = rlang::caller_env()) {
  if (setequal(names, expected)) {
    return(invisible(names))
  }
  # nolint start: object_usage_linter. Used by cli's glue interpolation.
  only_generated <- setdiff(names, expected)
  only_expected <- setdiff(expected, names)
  # nolint end
  cli::cli_abort(
    c(
      "The generated truth names are not the names that would be ranked.",
      x = if (length(only_generated) > 0L) {
        "Generated but not ranked: {.val {only_generated}}."
      },
      x = if (length(only_expected) > 0L) {
        "Ranked but not generated: {.val {only_expected}}."
      },
      i = "Generated: {.val {names}}.",
      i = "Expected: {.val {expected}}.",
      i = "{.pkg SBC} matches these by name and says nothing when one \\
           does not match, so this is an error rather than an empty \\
           rank histogram. Please report it at \\
           {.url https://github.com/GidonFrischkorn/bmmtools/issues}."
    ),
    call = call
  )
}

#' Split one parameter formula into its fixed part and its group terms
#'
#' On the parse tree, not on the deparsed text. A regex for the
#' outermost parentheses cannot be written, and one for the innermost
#' reads `(1 | gr(id, cor = FALSE))` as the group term `(id, cor =
#' FALSE)` --- a grouping factor called `id, cor = FALSE`, for which
#' `sbc_variables()` would then build patterns that match nothing on the
#' real fit, in silence. That is the failure this whole layer exists to
#' prevent, so the parsing is structural.
#'
#' A top-level term is a group term when it is a parenthesis holding a
#' `|` or `||` call; everything else at the top level is a fixed term.
#'
#' @return A list with `fixed`, the deparsed fixed terms, and `groups`,
#'   one entry per group term as `sbc_group_term()` reads it.
#' @noRd
sbc_formula_parts <- function(f, parameter, call = rlang::caller_env()) {
  terms <- flatten_sum(f[[3L]])
  is_group <- vapply(terms, is_group_term, logical(1))
  list(
    fixed = vapply(terms[!is_group], deparse1, character(1)),
    groups = lapply(
      terms[is_group], sbc_group_term,
      parameter = parameter, call = call
    )
  )
}

#' The summands of a `+` tree, in order
#' @noRd
flatten_sum <- function(e) {
  if (is.call(e) && identical(e[[1L]], quote(`+`)) && length(e) == 3L) {
    return(c(flatten_sum(e[[2L]]), flatten_sum(e[[3L]])))
  }
  list(e)
}

#' Is this top-level term a parenthesis holding a bar?
#' @noRd
is_group_term <- function(e) {
  if (!is.call(e) || !identical(e[[1L]], quote(`(`))) {
    return(FALSE)
  }
  inner <- e[[2L]]
  is.call(inner) && length(inner) == 3L &&
    as.character(inner[[1L]]) %in% c("|", "||")
}

#' Read one group term: its effects, its grouping factor, its correlation
#'
#' `(1 | id)`, `(1 | p | id)` and `(1 || id)`. R parses `1 | p | id`
#' left-associatively as `(1 | p) | id`, so a `|` on the left of a `|`
#' is brms's correlation id, which is what puts `cor_` draws in the fit.
#'
#' The grouping factor must be a bare name. `gr(id, cor = FALSE)`,
#' `mm(id1, id2)` and `id:session` are all formulas brms accepts and
#' none of them is a name `sbc_variables()` can build `sd_`, `cor_` and
#' `r_` patterns from, so each is refused here rather than turned into a
#' pattern that matches nothing.
#'
#' @noRd
sbc_group_term <- function(e, parameter, call = rlang::caller_env()) {
  inner <- e[[2L]]
  lhs <- inner[[2L]]
  group <- inner[[3L]]
  correlated <- identical(as.character(inner[[1L]]), "|") &&
    is.call(lhs) && identical(as.character(lhs[[1L]]), "|")
  effects <- if (correlated) lhs[[2L]] else lhs

  if (!is.symbol(group)) {
    cli::cli_abort(
      c(
        "The grouping factor of {.val {parameter}} is \\
         {.code {deparse1(group)}}, not a name.",
        i = "{.fn sbc} builds the {.field sd_}, {.field cor_} and \\
             {.field r_} draw names from the grouping factor, so it has \\
             to be a bare column name: {.code (1 | id)}, not \\
             {.code (1 | gr(id))} or {.code (1 | id:session)}."
      ),
      call = call
    )
  }
  list(
    effects = deparse1(effects),
    group = as.character(group),
    correlated = correlated
  )
}

#' The grouping factor a formula varies over, and whether it correlates
#'
#' Read without checking the fixed part, which is what `sbc(generator =)`
#' needs: with a user generator bmmtools never interprets the draw row,
#' so any fixed part is fine, but `level` still needs the group name to
#' build the `sd_`, `cor_` and `r_` patterns.
#'
#' @return A list with `group` (`NULL` when there is no group term) and
#'   `correlated`.
#' @noRd
sbc_formula_groups <- function(formula, call = rlang::caller_env()) {
  if (!inherits(formula, "bmmformula")) {
    cli::cli_abort(
      "{.arg formula} must be a {.cls bmmformula}, \\
       not {.obj_type_friendly {formula}}.",
      call = call
    )
  }
  parameters <- names(formula) %||% seq_along(formula)
  found <- list()
  for (i in seq_along(formula)) {
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    parameter <- as.character(parameters[[i]])
    parts <- sbc_formula_parts(formula[[i]], parameter, call = call)
    if (length(parts$groups) == 0L) next
    if (length(parts$groups) > 1L) {
      # nolint next: object_usage_linter. Used by cli's glue interpolation.
      terms <- vapply(parts$groups, function(g) g$group, character(1))
      cli::cli_abort(
        c(
          "{.val {parameter}} varies over {length(parts$groups)} grouping \\
           factors: {.val {terms}}.",
          i = "{.fn sbc} simulates one set of subjects, so each parameter \\
               has at most one group term."
        ),
        call = call
      )
    }
    found[[parameter]] <- parts$groups[[1L]]
  }

  if (length(found) == 0L) {
    return(list(group = NULL, correlated = FALSE))
  }
  groups <- vapply(found, function(x) x$group, character(1))
  if (length(unique(groups)) > 1L) {
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    labels <- paste0(names(groups), " (", groups, ")")
    cli::cli_abort(
      c(
        "The parameter formulas vary over more than one grouping \\
         factor: {.val {unique(groups)}}.",
        x = "Found: {.val {labels}}.",
        i = "{.fn sbc} simulates one set of subjects, so every parameter \\
             varies over the same one."
      ),
      call = call
    )
  }
  list(
    group = unname(groups[[1L]]),
    correlated = any(vapply(found, function(x) x$correlated, logical(1)))
  )
}

#' Check that every parameter formula is one the default generator can read
#'
#' The default generator maps a prior draw row onto
#' `simulate_recovery()`'s `pars`, `sds` and `cors`. That map is exact
#' for the shapes `recovery_formula()` writes --- an intercept and at
#' most one group term --- and is guesswork for a covariate or a task
#' factor, so those are refused and the message points at `generator`,
#' which lifts the restriction by taking the draw row whole (spec,
#' decision (d)).
#'
#' @return As `sbc_formula_groups()`.
#' @noRd
check_sbc_formula <- function(formula, call = rlang::caller_env()) {
  out <- sbc_formula_groups(formula, call = call)
  # nolint next: object_usage_linter. Used by cli's glue interpolation.
  parameters <- names(formula) %||% seq_along(formula)

  for (i in seq_along(formula)) {
    parts <- sbc_formula_parts(
      formula[[i]], as.character(parameters[[i]]),
      call = call
    )
    effects <- c(
      parts$fixed,
      vapply(parts$groups, function(g) g$effects, character(1))
    )
    offending <- effects[effects != "1"]
    if (length(offending) > 0L) {
      cli::cli_abort(
        c(
          "{.val {parameters[[i]]}} is not intercept-only: \\
           {.val {offending}}.",
          i = "The default generator maps a prior draw onto \\
               {.fn simulate_recovery}, which is exact for the formulas \\
               {.fn recovery_formula} writes and guesswork for a \\
               covariate or a task factor.",
          i = "Supply {.arg generator} to fit any formula: it takes the \\
               draw row whole, so bmmtools never has to interpret it."
        ),
        call = call
      )
    }
  }
  out
}

#' The design `data` describes: how many subjects, how many trials each
#'
#' Only the layout is read; the response values are ignored, because the
#' generator replaces them. The design must be balanced, because
#' `simulate_recovery()` is.
#'
#' @return A list with `n_subjects`, `n_trials` and `group`.
#' @noRd
sbc_layout <- function(data, group, call = rlang::caller_env()) {
  if (!is.data.frame(data)) {
    cli::cli_abort(
      "{.arg data} must be a data frame, not {.obj_type_friendly {data}}.",
      call = call
    )
  }
  if (!group %in% names(data)) {
    cli::cli_abort(
      c(
        "{.arg data} has no column {.val {group}}.",
        i = "Available: {.val {names(data)}}."
      ),
      call = call
    )
  }
  counts <- table(as.character(data[[group]]))
  if (length(counts) == 0L) {
    cli::cli_abort(
      "{.arg data} has no rows, so there is no design to simulate over.",
      call = call
    )
  }

  trials <- as.integer(counts)
  modal <- trials[[which.max(tabulate(match(trials, unique(trials))))]]
  if (length(unique(trials)) > 1L) {
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    odd <- paste0(
      names(counts)[trials != modal], " (", trials[trials != modal], ")"
    )
    cli::cli_abort(
      c(
        "{.arg data} is unbalanced: every subject must have the same \\
         number of rows.",
        x = "Most have {modal}; these do not: {.val {odd}}.",
        i = "{.fn simulate_recovery} generates a balanced design, so the \\
             layout it is given has to be one."
      ),
      call = call
    )
  }

  list(n_subjects = length(counts), n_trials = modal, group = group)
}

#' The `...` arguments `sbc()` refuses, and where each one belongs
#'
#' `cores` collides inside `SBC_fit.SBC_backend_function()`, which sets
#' it from `cores_arg` and then fails in `do.call()` with "formal
#' argument matched by multiple actual arguments" (measured 2026-09-15).
#' The rest are set by `sbc()` itself or belong to `fit_cached()`, which
#' the per-dataset fits do not go through.
#'
#' @noRd
sbc_refused_dots <- function() {
  c(
    cores = "use {.arg cores_per_fit}, which SBC gives to each fit",
    sample_prior = "{.fn sbc} sets it for the prior fit",
    file = "the file arguments belong to {.fn fit_cached}; the dataset \\
            fits are cached by SBC through {.arg cache_mode}",
    file_refit = "the file arguments belong to {.fn fit_cached}; use \\
                  {.arg refit}, which reaches the prior fit",
    file_compress = "the file arguments belong to {.fn fit_cached}; the \\
                     dataset fits do not go through it"
  )
}

#' @noRd
check_sbc_dots <- function(dots, call = rlang::caller_env()) {
  if (length(dots) == 0L) {
    return(invisible(dots))
  }
  given <- names(dots) %||% rep("", length(dots))
  if (any(!nzchar(given))) {
    cli::cli_abort(
      c(
        "Every argument in {.arg ...} must be named.",
        i = "They are passed on to the fitter by name."
      ),
      call = call
    )
  }

  refused <- sbc_refused_dots()
  offending <- intersect(given, names(refused))
  if (length(offending) == 0L) {
    return(invisible(dots))
  }
  bullets <- stats::setNames(
    paste0("{.arg ", offending, "}: ", refused[offending]),
    rep("x", length(offending))
  )
  cli::cli_abort(
    c(
      "{.fn sbc} refuses {cli::qty(length(offending))}{?this/these} \\
       argument{?s} in {.arg ...}: {.val {offending}}.",
      bullets
    ),
    call = call
  )
}

#' Resolve the patterns to draw names, refusing any that found nothing
#'
#' `posterior::subset_draws(regex = TRUE)` returns a zero-column object
#' rather than erroring when a pattern matches nothing (measured
#' 2026-09-16), and it does so per pattern: a selection where one of four
#' patterns is wrong comes back three columns wide, with no error and no
#' message. Checking only that *something* matched would reproduce that
#' silence one level up, so the coverage is checked per **variable**.
#'
#' Per variable, not per pattern, because a correlation contributes two
#' patterns --- one per order of the pair --- of which exactly one can
#' match. `sbc_variables()` gives both the same `"<level>:<term>"` name,
#' and it is the names that have to be covered.
#'
#' @return The matched draw names, in the order the patterns were given.
#' @noRd
select_variables <- function(variables, available,
                             call = rlang::caller_env()) {
  matched <- lapply(variables, function(p) grep(p, available, value = TRUE))
  wanted <- names(variables) %||% variables
  found <- vapply(
    unique(wanted),
    function(w) sum(lengths(matched[wanted == w])),
    integer(1)
  )

  if (any(found == 0L)) {
    # nolint start: object_usage_linter. Used by cli's glue interpolation.
    missing <- names(found)[found == 0L]
    patterns <- unname(variables[wanted %in% missing])
    # nolint end
    cli::cli_abort(
      c(
        "No draw of the prior fit matches {cli::qty(length(missing))}\\
         {?this variable/these variables}: {.val {missing}}.",
        x = "Looked for: {.val {patterns}}.",
        i = "The fit has: {.val {available}}.",
        i = "{.pkg SBC} would rank what it found and say nothing about \\
             the rest, so this is an error rather than a shorter run."
      ),
      call = call
    )
  }
  unique(unlist(matched, use.names = FALSE))
}

#' Prior draws of the variables that will be ranked
#'
#' An internal generic for the same reason `prior_predict_draws()` is
#' one: bmm's mock backend returns a `bmmfit` whose `$fit` is a bare
#' number, so nothing that reads draws runs on it (measured 2026-09-08),
#' and the suite may not compile Stan. A test registers a `mockfit`
#' method; the default method is the real path.
#'
#' The spec sketches this as `prior_parameter_draws(fit, n_sims, seed)`
#' while saying in the same paragraph that the patterns come from
#' `sbc_variables()` and are used here, so `variables` is the third
#' formal: a generic that subsets has to be told what to subset to.
#'
#' @noRd
prior_parameter_draws <- function(fit, n_sims, variables, seed = NULL, ...) {
  UseMethod("prior_parameter_draws")
}

#' @noRd
#' @export
prior_parameter_draws.default <- function(fit, n_sims, variables,
                                          seed = NULL, ...) {
  rlang::check_installed("brms", "to read prior draws from a fit.")
  n_sims <- check_count(n_sims, "n_sims")

  draws <- posterior::as_draws_matrix(fit$fit)
  selected <- posterior::subset_draws(
    draws,
    variable = select_variables(variables, posterior::variables(draws))
  )

  available <- posterior::ndraws(selected)
  if (available < n_sims) {
    cli::cli_abort(c(
      "The prior fit has {available} draw{?s}, fewer than the \\
       {n_sims} {.arg n_sims} asks for.",
      i = "Raise {.arg iter} or {.arg chains}, or lower {.arg n_sims}."
    ))
  }

  # the sample is across chains, so chains stop meaning anything; saying
  # so here is what keeps posterior from saying it on every call
  selected <- posterior::merge_chains(selected)
  rows <- with_seed_if(seed, sample.int(available, n_sims))
  posterior::as_draws_matrix(posterior::subset_draws(selected, draw = rows))
}
