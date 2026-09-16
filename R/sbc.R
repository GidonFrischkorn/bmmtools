# Simulation-based calibration: the naming layer.
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
# asserted with `check_variable_names()`, which errors. The generator's
# `variables` carry the fit's brms draw names (`b_kappa_Intercept`, not
# `kappa`), because those are what SBC matches against (measured
# 2026-09-16).

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
#' @param broad `FALSE` builds the exact, intercept-only shapes above,
#'   for the default generator, whose formulas have no other shape.
#'   `TRUE` builds the prefixes a user `generator` needs (spec (i),
#'   decided A): every `b_<par>_*`, every `sd_<g>__<par>_*`, every
#'   `cor_<g>__*` (one pattern, named after the group) and every
#'   `r_<g>__<par>[*`, so that a cell-means `b_kappa_task1` or a random
#'   slope's `sd_id__kappa_cond` is a draw of its level. The trailing
#'   `_` after the parameter is what keeps `^b_kappa_` off
#'   `b_kappa2_Intercept`.
#' @return A character vector of regexes named `"<level>:<term>"`. The
#'   name says which variable each pattern is there to find, which is
#'   what lets `prior_parameter_draws()` report a pattern that found
#'   nothing instead of quietly returning a narrower matrix.
#' @noRd
sbc_variables <- function(model, group, level, broad = FALSE,
                          call = rlang::caller_env()) {
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
    out <- add(
      out,
      if (broad) paste0("^b_", esc, "_") else paste0("^b_", esc, "_Intercept$"),
      "population", free
    )
  }
  if (is.null(group)) {
    return(out)
  }
  g <- escape_regex(group)

  if ("sd" %in% level) {
    out <- add(
      out,
      if (broad) {
        paste0("^sd_", g, "__", esc, "_")
      } else {
        paste0("^sd_", g, "__", esc, "_Intercept$")
      },
      "sd", free
    )
  }
  if ("cor" %in% level && broad) {
    out <- add(out, paste0("^cor_", g, "__"), "cor", group)
  } else if ("cor" %in% level && length(free) > 1L) {
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
      out,
      if (broad) {
        paste0("^r_", g, "__", esc, "\\[")
      } else {
        paste0("^r_", g, "__", esc, "\\[.+,Intercept\\]$")
      },
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
#' @return A list with `group` (`NULL` when there is no group term),
#'   `correlated`, and `varying`, the parameters whose own formula
#'   carries the group term. `varying` is what keeps `sbc()` from asking
#'   for an `sd_` draw of a parameter that does not vary (spec N5): a
#'   formula where `kappa` varies and `thetat` does not is one the
#'   default generator simulates perfectly well, and building the
#'   patterns over every free parameter would refuse it.
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
    return(list(group = NULL, correlated = FALSE, varying = character(0)))
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
    correlated = any(vapply(found, function(x) x$correlated, logical(1))),
    varying = names(found)
  )
}

#' Check that every parameter formula is one the default generator can read
#'
#' The default generator maps a prior draw row onto
#' `simulate_recovery()`'s `pars`, `sds` and `cors`. That map is exact
#' for the shapes `recovery_formula()` writes --- an intercept and at
#' most one group term --- and is guesswork for a covariate or a task
#' factor, so those are refused and the message points at `generator`,
#' which lifts the restriction by taking the draw row whole.
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
#' `ids` are the labels of the grouping column in the order `data` gives
#' them. They matter for the subject level: the prior fit is fitted to
#' `data` and names its `r_` draws after these labels, while
#' `simulate_recovery()` numbers its subjects 1..n, so the generator
#' relabels the simulated column to these before the two are matched by
#' name.
#'
#' @param balanced `TRUE` refuses a design whose subjects have different
#'   numbers of rows, because `simulate_recovery()` cannot produce one.
#'   `FALSE`, for a user `generator` that owns the data, records the
#'   modal count as `n_trials` and lets the design be.
#' @return A list with `n_subjects`, `n_trials`, `group` and `ids`.
#' @noRd
sbc_layout <- function(data, group, balanced = TRUE,
                       call = rlang::caller_env()) {
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
  if (balanced && length(unique(trials)) > 1L) {
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

  list(
    n_subjects = length(counts), n_trials = modal, group = group,
    ids = unique(as.character(data[[group]]))
  )
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
  subsample_prior_draws(
    posterior::as_draws_matrix(fit$fit),
    check_count(n_sims, "n_sims"), variables, seed
  )
}

#' Subset a draws object to the variables and sample `n_sims` of its rows
#'
#' The body of `prior_parameter_draws.default()`, kept apart so that the
#' `mockfit` method of the test suite runs the same coverage check and
#' the same draw-count error rather than a second implementation of
#' them: a mock that is laxer than the real path would hide exactly the
#' mismatch this layer exists to catch.
#'
#' @noRd
subsample_prior_draws <- function(draws, n_sims, variables, seed = NULL,
                                  call = rlang::caller_env()) {
  selected <- posterior::subset_draws(
    draws,
    variable = select_variables(
      variables, posterior::variables(draws),
      call = call
    )
  )

  available <- posterior::ndraws(selected)
  if (available < n_sims) {
    cli::cli_abort(
      c(
        "The prior fit has {available} draw{?s}, fewer than the \\
         {n_sims} {.arg n_sims} asks for.",
        i = "Raise {.arg iter} or {.arg chains}, or lower {.arg n_sims}."
      ),
      call = call
    )
  }

  # the sample is across chains, so chains stop meaning anything; saying
  # so here is what keeps posterior from saying it on every call
  selected <- posterior::merge_chains(selected)
  rows <- with_seed_if(seed, sample.int(available, n_sims))
  posterior::as_draws_matrix(posterior::subset_draws(selected, draw = rows))
}

# Simulation-based calibration: the run.
#
# The half above builds and checks names; this half builds the two
# objects SBC asks for and hands them over. Nothing here computes a
# rank: `sbc()` returns SBC's own `SBC_results` so that its plots and
# its tests apply unchanged.

#' The variable names a set of varying parameters can carry
#'
#' `sbc_variables()` builds its `sd_`, `cor_` and `r_` patterns over
#' every free parameter of the model, which is right when the formula
#' gives every parameter a group term and wrong when it does not: a
#' formula where `kappa` varies and `thetat` does not would otherwise ask
#' for `sd_id__thetat_Intercept`, which no fit of it has, and
#' `select_variables()` would refuse a run the default generator can do
#' perfectly well (spec N5).
#'
#' @noRd
sbc_varying_names <- function(varying) {
  pairs <- character(0)
  if (length(varying) > 1L) {
    combos <- utils::combn(varying, 2L)
    pairs <- vapply(
      seq_len(ncol(combos)),
      function(k) pair_term(combos[[1L, k]], combos[[2L, k]])$term,
      character(1)
    )
  }
  c(
    paste0("sd:", varying),
    paste0("cor:", pairs),
    paste0("subject:", varying)
  )
}

#' Why a level the user asked for cannot be ranked, or `NA` when it can
#' @noRd
sbc_level_obstacle <- function(level, groups) {
  varying <- groups$varying
  if (level == "population" || length(varying) > 0L) {
    if (level != "cor") {
      return(NA_character_)
    }
    if (!groups$correlated) {
      return("its group terms are not correlated, so the fit has no \\
              {.field cor_} draws; write {.code (1 | p | id)} to \\
              correlate them")
    }
    if (length(varying) < 2L) {
      return("only {.val {varying}} varies between subjects, and a \\
              correlation needs two parameters that do")
    }
    return(NA_character_)
  }
  "no parameter of the formula has a group term, so the fit has no \\
   between-subject draws"
}

#' Which levels are drawn from the prior and which of them are ranked
#'
#' Two pattern sets, not one. `level` says what is **ranked**; what is
#' **drawn** is everything the formula implies, because SBC's uniformity
#' rests on the data coming from the joint prior. Holding the
#' between-subject SDs at zero while the fitted model has a prior on them
#' would break that joint and take the population ranks down with it ---
#' silently, as a run that looks short rather than wrong (spec N5).
#'
#' The subject level is a third case: ranked when asked, never drawn.
#' The subject values are produced by `simulate_recovery()` from the
#' drawn SDs and correlations, which is the prior conditional, so taking
#' them off the prior row as well would simulate from one draw and rank
#' against another. Their *names* are still taken from the prior fit,
#' because that is what makes the generator's constructed
#' `r_<g>__<par>[<label>,Intercept]` names checkable against real ones
#' before a single data set is fitted.
#'
#' @return A list with `draw` (the patterns whose values feed the
#'   simulation), `rank` (the subset SBC is asked to rank), `take` (every
#'   pattern to subset from the prior fit, the union of the two) and
#'   `level` (the levels that survived).
#' @noRd
sbc_pattern_sets <- function(model, groups, group, level,
                             call = rlang::caller_env()) {
  dropped <- character(0)
  for (lvl in level) {
    obstacle <- sbc_level_obstacle(lvl, groups)
    if (is.na(obstacle)) next
    dropped <- c(dropped, lvl)
    cli::cli_inform(c(
      "Dropping the {.val {lvl}} level: {obstacle}.",
      i = "The other levels are ranked as asked."
    ))
  }
  rank <- setdiff(level, dropped)

  # what is drawn is decided by the same obstacles, over all three
  # drawable levels --- not by `dropped`, which only ever holds levels
  # the caller asked to rank. Reading it off `dropped` would leave `cor`
  # in the drawn set whenever `level` omitted it, and the default
  # `recovery_formula()` writes uncorrelated group terms, so the
  # commonest call of all would ask a real fit for a `cor_` draw it does
  # not have --- after the prior fit had already run.
  drawable <- c("population", "sd", "cor")
  draw <- drawable[
    is.na(vapply(drawable, sbc_level_obstacle, character(1), groups = groups))
  ]
  has_group <- length(groups$varying) > 0L
  patterns <- sbc_variables(
    model, if (has_group) group else NULL, union(draw, rank),
    call = call
  )
  levels <- sbc_variable_levels(patterns)
  patterns <- patterns[
    levels == "population" | names(patterns) %in%
      sbc_varying_names(groups$varying)
  ]
  levels <- sbc_variable_levels(patterns)

  list(
    draw = patterns[levels %in% draw],
    rank = patterns[levels %in% rank],
    take = patterns,
    level = rank
  )
}

#' The subject truths as brms names them: deviations, under the labels
#'
#' brms's `r_<g>__<par>[<label>,Intercept]` draws are deviations from
#' the intercept --- measured on the fixture 2026-09-16, the eight
#' `r_id__kappa[i,Intercept]` posterior means average -0.0278 against a
#' `b_kappa_Intercept` of 1.8692 --- while `simulate_recovery()`'s
#' `truth$subjects` are absolute. Emitting the absolute values would give
#' ranks that look like a badly miscalibrated model and are in fact a
#' units error, which is why this subtraction has its own regression
#' test.
#'
#' Matched on `id` and `term`, never on position. The simulation numbers
#' its subjects 1..n and the prior fit names its draws after the labels
#' of the user's grouping column, so the n-th simulated subject is
#' emitted under the n-th label.
#'
#' @param truth A simulation's `truth`, with `population` and `subjects`.
#' @param group The grouping factor, as it appears in the draw names.
#' @param ids The labels, in the order the simulated subjects map onto.
#' @return A named numeric vector, one deviation per subject and varying
#'   parameter, on the link scale.
#' @noRd
sbc_subject_truths <- function(truth, group, ids) {
  subjects <- truth$subjects
  if (nrow(subjects) == 0L) {
    # paste0() would recycle the empty term to one name, not to none
    return(stats::setNames(numeric(0), character(0)))
  }
  population <- stats::setNames(
    truth$population$true_value, truth$population$term
  )
  deviations <- subjects$true_value - unname(population[subjects$term])
  labels <- ids[as.integer(subjects$id)]
  stats::setNames(
    as.double(deviations),
    paste0("r_", group, "__", subjects$term, "[", labels, ",Intercept]")
  )
}

#' Give the simulated group column the formula's name and the layout's labels
#'
#' `simulate_recovery()` always calls the column `id` and numbers the
#' subjects 1..n. The formula may call it something else, in which case
#' brms would be the one to notice, and the prior fit names its `r_`
#' draws after the labels the user's `data` had, which the data set fits
#' have to reproduce for SBC to match a subject truth to a subject draw.
#'
#' @noRd
relabel_subjects <- function(data, group, ids) {
  names(data)[names(data) == "id"] <- group
  # the simulated column is a factor with levels 1..n in order, so its
  # integer codes are the subject numbers
  data[[group]] <- factor(ids[as.integer(data[[group]])], levels = ids)
  data
}

#' The `cor_` draw name of each pair of varying parameters, in either order
#'
#' Which order brms writes is not guaranteed, so both are looked for and
#' the one the fit has is kept. A pair the fit has no draw for is left
#' out and its correlation stays 0.
#'
#' @noRd
sbc_cor_names <- function(varying, group, available) {
  if (length(varying) < 2L) {
    return(list())
  }
  combos <- utils::combn(varying, 2L)
  out <- list()
  for (k in seq_len(ncol(combos))) {
    a <- combos[[1L, k]]
    b <- combos[[2L, k]]
    found <- intersect(
      c(
        paste0("cor_", group, "__", a, "_Intercept__", b, "_Intercept"),
        paste0("cor_", group, "__", b, "_Intercept__", a, "_Intercept")
      ),
      available
    )
    if (length(found) > 0L) {
      out[[length(out) + 1L]] <- list(a = a, b = b, name = found[[1L]])
    }
  }
  out
}

#' The correlation matrix one prior draw implies, or `NULL` for none
#' @noRd
sbc_cor_matrix <- function(row, pairs, varying) {
  if (length(pairs) == 0L) {
    return(NULL)
  }
  out <- diag(length(varying))
  dimnames(out) <- list(varying, varying)
  for (pair in pairs) {
    out[pair$a, pair$b] <- row[[pair$name]]
    out[pair$b, pair$a] <- row[[pair$name]]
  }
  out
}

#' Simulate one data set from one prior draw, saying which one if it fails
#'
#' A prior draw is not a plausible parameter value, it is whatever the
#' prior allows, and a wide prior on a link-scale parameter allows a lot.
#' Measured 2026-09-16 on `bmm::mixture2p()` under bmm's own defaults:
#' the group-level SD of `kappa` has a half-`student_t(3, 0, 2.5)` prior,
#' whose draws reach 9.5 on the log scale, so a subject's concentration
#' reaches exp(2 + 2 x 9.5), and `rmixture2p()` dies inside its own
#' sampler with `node stack overflow` -- a message that names neither the
#' simulation nor a parameter.
#'
#' The failure is real and belongs to the prior, not to bmmtools, so this
#' does not catch it. It says which draw it was and what that draw held,
#' because the alternative is an hour of Stan ending in three words.
#'
#' @noRd
simulate_from_draw <- function(row, row_number, model, layout, population,
                               free, sds, varying, pairs,
                               call = rlang::caller_env()) {
  tryCatch(
    simulate_recovery(
      model,
      pars = stats::setNames(unname(row[population]), free),
      n_subjects = layout$n_subjects,
      n_trials = layout$n_trials,
      sds = if (length(sds) > 0L) {
        stats::setNames(unname(row[sds]), varying)
      },
      cors = sbc_cor_matrix(row, pairs, varying),
      seed = NULL
    ),
    error = function(e) {
      # nolint next: object_usage_linter. Used by cli's glue interpolation.
      shown <- paste0(names(row), " = ", format(row, digits = 4))
      cli::cli_abort(
        c(
          "Simulating data set {row_number} from its prior draw failed.",
          x = "The draw was: {.val {shown}}.",
          i = "These are values the {.emph prior} allows, on the link \
               scale, not values anyone would fit. A prior wide enough \
               to put mass where the model cannot generate is itself \
               the finding: tighten it and run {.fn sbc} again, or check \
               it first with {.fn prior_check}."
        ),
        parent = e, call = call
      )
    }
  )
}

#' The SBC generator: one prior draw, one simulated data set
#'
#' `SBC::generate_datasets()` calls the function with no arguments, once
#' per simulation, so the row it is on lives in this closure.
#' `future.chunk.size = Inf` pins SBC's sequential branch (`replicate()`,
#' measured 2026-09-16): under the futures branch each worker would get
#' its own copy of the counter and emit the same rows twice.
#'
#' @param draws The prior draws, `n_sims` rows under the fit's own draw
#'   names --- **every** variable the formula implies, not only the ones
#'   that will be ranked (spec N5).
#' @param rank The draw names to emit as `variables`, which is what SBC
#'   ranks. Any `r_<group>__` name among them is a subject truth and is
#'   emitted as the simulation's own deviation, never as the prior row's
#'   value (see `sbc_subject_truths()`).
#' @param layout `sbc_layout()`'s list. Without `ids` the labels are the
#'   subject numbers themselves.
#' @noRd
sbc_generator <- function(draws, rank, model, layout, correlated, group,
                          call = rlang::caller_env()) {
  # a plain matrix, not `as.matrix()`: that keeps the `draws_matrix`
  # class, whose `[` drops the variable names one row at a time and
  # would then index every value as NA
  draws <- posterior::as_draws_matrix(draws)
  available <- posterior::variables(draws)
  values <- matrix(
    as.numeric(draws),
    nrow = posterior::ndraws(draws), dimnames = list(NULL, available)
  )
  free <- model_parameters(model)$free
  population <- paste0("b_", free, "_Intercept")

  sds <- stats::setNames(
    paste0("sd_", group, "__", free, "_Intercept"), free
  )
  sds <- sds[sds %in% available]
  varying <- names(sds)
  pairs <- if (correlated) {
    sbc_cor_names(varying, group, available)
  } else {
    list()
  }
  ids <- layout$ids %||% as.character(seq_len(layout$n_subjects))
  subject_names <- grep(
    paste0("^r_", escape_regex(group), "__"), rank,
    value = TRUE
  )

  # the row counter lives in its own environment rather than behind
  # `<<-`: SBC calls the generator with no arguments, so which row it is
  # on has to be kept somewhere, and an explicit environment says where
  state <- rlang::env(row = 0L) # nolint: object_usage_linter. Used in `f`.
  f <- function() {
    state$row <- state$row + 1L
    row_number <- state$row
    if (row_number > nrow(values)) {
      cli::cli_abort(
        c(
          "The generator has only {nrow(values)} prior draw{?s} and was \\
           asked for a {row_number}{.strong th} data set.",
          i = "{.fn sbc} draws exactly {.arg n_sims} rows, so \\
               {.fn SBC::generate_datasets} has to be called with the \\
               same number."
        ),
        call = call
      )
    }
    row <- values[row_number, ]

    simulation <- simulate_from_draw(
      row, row_number, model, layout, population, free, sds, varying,
      pairs, call
    )
    data <- relabel_subjects(simulation$data, group, ids)

    truths <- row[rank]
    if (length(subject_names) > 0L) {
      deviations <- sbc_subject_truths(simulation$truth, group, ids)
      check_variable_names(names(deviations), subject_names, call = call)
      truths[subject_names] <- deviations[subject_names]
    }
    check_variable_names(names(truths), rank, call = call)
    list(variables = as.list(truths), generated = data)
  }

  SBC::SBC_generator_function(f, future.chunk.size = Inf)
}

#' Refuse a `generator` that is not a function of `(draws, data)`
#' @noRd
check_sbc_generator <- function(generator, call = rlang::caller_env()) {
  if (is.null(generator)) {
    return(invisible(NULL))
  }
  if (!is.function(generator)) {
    cli::cli_abort(
      c(
        "{.arg generator} must be a function or {.code NULL}, \\
         not {.obj_type_friendly {generator}}.",
        i = "{.fn sbc} calls it as {.code generator(draws, data)}."
      ),
      call = call
    )
  }
  arguments <- names(formals(generator))
  if (!"..." %in% arguments && length(arguments) < 2L) {
    cli::cli_abort(
      c(
        "{.arg generator} takes {length(arguments)} argument{?s}, and \\
         {.fn sbc} calls it as {.code generator(draws, data)}.",
        i = "{.arg draws} is one prior draw as a named numeric vector \\
             under the fit's own draw names; {.arg data} is the layout \\
             frame."
      ),
      call = call
    )
  }
  invisible(NULL)
}

#' What a user generator receives and what is ranked, from the fit's names
#'
#' Generator mode (spec (i), decided A) does not build names from the
#' formula, because the formula may have any shape: it reads the prior
#' fit's own variable names and takes every `b_`, `sd_`, `cor_` and
#' `r_` draw of the model's free parameters. The generator receives all
#' of them (`take`), because a generator handed less than the fit puts a
#' prior on would simulate from a point mass where the fit has a prior,
#' and every rank would be off in silence. `level` then picks which
#' classes are ranked, and a class the fit has no draw of is dropped
#' with a message --- the formula can say why when it has no group term,
#' and otherwise the fit is the reason.
#'
#' @param available The prior fit's variable names.
#' @return As `sbc_pattern_sets()`, but `take` and `rank` are resolved
#'   draw names rather than patterns.
#' @noRd
sbc_generator_sets <- function(model, groups, group, level, available,
                               call = rlang::caller_env()) {
  has_group <- !is.null(groups$group)
  classes <- c("population", "sd", "cor", "subject")
  patterns <- sbc_variables(
    model, if (has_group) group else NULL, classes,
    broad = TRUE, call = call
  )
  pattern_levels <- sbc_variable_levels(patterns)
  matched <- lapply(patterns, function(p) grep(p, available, value = TRUE))
  names_of <- function(lvls) {
    unique(unlist(matched[pattern_levels %in% lvls], use.names = FALSE))
  }
  found <- vapply(
    classes, function(lvl) length(names_of(lvl)) > 0L, logical(1)
  )

  if (!found[["population"]]) {
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    looked <- unname(patterns[pattern_levels == "population"])
    cli::cli_abort(
      c(
        "No draw of the prior fit is a population-level draw of a free \\
         parameter.",
        x = "Looked for: {.val {looked}}.",
        i = "The fit has: {.val {available}}.",
        i = "Please report it at \\
             {.url https://github.com/GidonFrischkorn/bmmtools/issues}."
      ),
      call = call
    )
  }

  # nolint start: object_usage_linter. Used by cli's glue interpolation.
  prefixes <- c(
    sd = paste0("sd_", group, "__"),
    cor = paste0("cor_", group, "__"),
    subject = paste0("r_", group, "__")
  )
  # nolint end
  dropped <- character(0)
  for (lvl in level) {
    if (found[[lvl]]) next
    dropped <- c(dropped, lvl)
    reason <- if (!has_group) {
      "no parameter of the formula has a group term, so the fit has no \\
       between-subject draws"
    } else {
      "the prior fit has no {.field {prefixes[[lvl]]}} draws"
    }
    cli::cli_inform(c(
      paste0("Dropping the {.val {lvl}} level: ", reason, "."),
      i = "The other levels are ranked as asked."
    ))
  }
  rank <- setdiff(level, dropped)

  list(
    take = names_of(classes[found]),
    rank = names_of(rank),
    level = rank
  )
}

#' Check the data frame a user generator returned, naming the data set
#'
#' The model's response columns and, when the formula has a group term,
#' the grouping column have to be there, or brms would be the one to
#' say so after the prior fit and a compile. When the subject level is
#' ranked the grouping column also has to carry the layout's labels:
#' the truths are named after the prior fit's labels, and a data set
#' fitted under other labels would give SBC nothing to rank them
#' against, silently.
#'
#' @noRd
check_generated_data <- function(generated, row_number, model, group,
                                 needs_group, ids, subject_ranked,
                                 call = rlang::caller_env()) {
  if (!is.data.frame(generated)) {
    cli::cli_abort(
      c(
        "{.arg generator} returned {.obj_type_friendly {generated}} for \\
         data set {row_number}, not a data frame.",
        i = "It is called as {.code generator(draws, data)} and returns \\
             the data frame to fit."
      ),
      call = call
    )
  }
  needed <- c(
    unlist(model$resp_vars, use.names = FALSE),
    if (needs_group) group
  )
  missing <- setdiff(needed, names(generated))
  if (length(missing) > 0L) {
    # nolint start: object_usage_linter. Used by cli's glue interpolation.
    responses <- unlist(model$resp_vars, use.names = FALSE)
    returned <- names(generated)
    # nolint end
    cli::cli_abort(
      c(
        "The data {.arg generator} returned for data set {row_number} \\
         lack the column{?s} {.val {missing}}.",
        i = "It returned {.val {returned}}.",
        i = "The model's response column{?s}: {.val {responses}}.",
        i = if (needs_group) "The grouping column: {.val {group}}."
      ),
      call = call
    )
  }
  if (subject_ranked) {
    labels <- unique(as.character(generated[[group]]))
    if (!setequal(labels, ids)) {
      cli::cli_abort(
        c(
          "The data {.arg generator} returned for data set {row_number} \\
           carry other subject labels than {.arg data}.",
          x = "Returned: {.val {labels}}.",
          x = "The prior fit's labels: {.val {ids}}.",
          i = "The subject level names its truths after the prior fit's \\
               labels, so the data set fits have to use the same ones, \\
               or {.pkg SBC} would have nothing to rank them against."
        ),
        call = call
      )
    }
  }
  invisible(generated)
}

#' The SBC generator around a user function: one row in, one data set out
#'
#' bmmtools never interprets the row. It hands every `b_`, `sd_`, `cor_`
#' and `r_` draw of the prior fit to `generator(draws, data)` as a named
#' numeric vector, checks what comes back, and emits the ranked entries
#' of the same row as `variables` --- so the truth names equal the fit's
#' draw names by construction, for any formula.
#'
#' @param draws The prior draws, `n_sims` rows, already subset to what
#'   the generator receives.
#' @param rank The draw names SBC ranks, a subset of the columns.
#' @noRd
sbc_user_generator <- function(draws, rank, generator, data, model, layout,
                               group, needs_group,
                               call = rlang::caller_env()) {
  draws <- posterior::as_draws_matrix(draws)
  values <- matrix(
    as.numeric(draws),
    nrow = posterior::ndraws(draws),
    dimnames = list(NULL, posterior::variables(draws))
  )
  subject_ranked <- any(grepl(
    paste0("^r_", escape_regex(group), "__"), rank
  ))

  state <- rlang::env(row = 0L) # nolint: object_usage_linter. Used in `f`.
  f <- function() {
    state$row <- state$row + 1L
    row_number <- state$row
    if (row_number > nrow(values)) {
      cli::cli_abort(
        c(
          "The generator has only {nrow(values)} prior draw{?s} and was \\
           asked for a {row_number}{.strong th} data set.",
          i = "{.fn sbc} draws exactly {.arg n_sims} rows, so \\
               {.fn SBC::generate_datasets} has to be called with the \\
               same number."
        ),
        call = call
      )
    }
    row <- values[row_number, ]

    generated <- tryCatch(
      generator(row, data),
      error = function(e) {
        # nolint next: object_usage_linter. Used by cli's glue interpolation.
        shown <- paste0(names(row), " = ", format(row, digits = 4))
        cli::cli_abort(
          c(
            "{.arg generator} failed on data set {row_number}.",
            x = "The draw was: {.val {shown}}.",
            i = "These are values the {.emph prior} allows, on the link \\
                 scale, not values anyone would fit."
          ),
          parent = e, call = call
        )
      }
    )
    check_generated_data(
      generated, row_number, model, group, needs_group, layout$ids,
      subject_ranked, call = call
    )

    variables <- as.list(row[rank])
    check_variable_names(names(variables), rank, call = call)
    list(variables = variables, generated = generated)
  }

  SBC::SBC_generator_function(f, future.chunk.size = Inf)
}

#' The SBC backend: one data set, one fit
#'
#' The lambda takes `cores` because `cores_arg = "cores"` makes
#' `SBC_fit.SBC_backend_function()` put it into the call, and a function
#' without the formal cannot receive it. The returned `bmmfit` is read
#' by SBC's own
#' `brmsfit` methods, so there is no backend class to write.
#'
#' @noRd
sbc_backend <- function(formula, model, prior, dots, fitter) {
  fit_one <- function(generated, cores) {
    rlang::exec(
      fitter,
      formula = formula, data = generated, model = model, prior = prior,
      cores = cores, !!!dots
    )
  }
  SBC::SBC_backend_function(
    fit_one,
    generated_arg = "generated", cores_arg = "cores"
  )
}

#' Count the fits past each convergence bar, and warn once
#'
#' Two bars apply to the same fits: bmmtools' Rhat 1.05 and SBC's own
#' 1.01. Shipping both silently is worse than either, so this
#' states bmmtools' and names SBC's. `results$default_diagnostics`
#' survives `keep_fits = FALSE` (measured 2026-09-16), so nothing has to
#' be kept to read it.
#'
#' `check_convergence()`'s `ess_bulk_min = 400` is deliberately not
#' applied here:
#' SBC-length fits are short by design and it would fire on nearly every
#' run, which is how a warning stops being read.
#'
#' @noRd
sbc_diagnostics <- function(results, rhat_max = 1.05,
                            ess_to_rank_min = 0.5) {
  out <- list(
    n_fits = NA_integer_,
    rhat_max = rhat_max,
    ess_to_rank_min = ess_to_rank_min,
    n_high_rhat = NA_integer_,
    n_low_ess_to_rank = NA_integer_,
    n_missing = NA_integer_
  )
  diagnostics <- results$default_diagnostics
  needed <- c("max_rhat", "min_ess_to_rank")
  if (!is.data.frame(diagnostics) || !all(needed %in% names(diagnostics))) {
    return(out)
  }

  out$n_fits <- nrow(diagnostics)
  out$n_high_rhat <- sum(diagnostics$max_rhat > rhat_max, na.rm = TRUE)
  out$n_low_ess_to_rank <- sum(
    diagnostics$min_ess_to_rank < ess_to_rank_min,
    na.rm = TRUE
  )
  out$n_missing <- sum(
    is.na(diagnostics$max_rhat) | is.na(diagnostics$min_ess_to_rank)
  )
  if (out$n_high_rhat + out$n_low_ess_to_rank + out$n_missing == 0L) {
    return(out)
  }

  cli::cli_warn(
    c(
      "{out$n_fits} fit{?s} of the calibration, and some did not converge.",
      x = "{out$n_high_rhat} {?has/have} an Rhat above {rhat_max}.",
      x = "{out$n_low_ess_to_rank} {?has/have} a tail ESS below \\
           {ess_to_rank_min} of the maximum rank, which skews the ranks \\
           themselves.",
      x = if (out$n_missing > 0L) {
        "{out$n_missing} {?has/have} no Rhat or no ESS at all."
      },
      i = "{.pkg SBC}'s own {.fn summary} uses the stricter Rhat 1.01; \\
           these counts are against {.pkg bmmtools}' {rhat_max}.",
      i = "Fit with {.code keep_fits = TRUE} and run \\
           {.fn check_convergence} per fit to see which."
    ),
    class = "bmmtools_sbc_diagnostics"
  )
  out
}

#' Refuse a grouping column with missing values
#'
#' `sbc_layout()` counts with `table()`, which drops `NA` without saying
#' so: a design with a missing id would be read off the rows that happen
#' to have one, and the simulation would be of a smaller study than the
#' user described (6.2's review, deferred to here).
#'
#' @noRd
check_group_ids <- function(data, group, call = rlang::caller_env()) {
  if (!group %in% names(data)) {
    return(invisible(data))
  }
  # nolint next: object_usage_linter. Used by cli's glue interpolation.
  missing <- sum(is.na(data[[group]]))
  if (missing > 0L) {
    cli::cli_abort(
      c(
        "{.arg data} has {missing} row{?s} whose {.val {group}} is \\
         {.code NA}.",
        i = "The layout is counted per subject, and a missing id would \\
             be dropped from that count without changing the design \\
             {.fn sbc} then simulates."
      ),
      call = call
    )
  }
  invisible(data)
}

#' Validate the arguments `sbc()` does not hand straight on
#' @noRd
check_sbc_args <- function(prior, seed, fitter, cache_mode, cache_location,
                           call = rlang::caller_env()) {
  if (!is.null(prior) && !inherits(prior, "brmsprior")) {
    cli::cli_abort(
      c(
        "{.arg prior} must be a {.cls brmsprior} or {.code NULL}, \\
         not {.obj_type_friendly {prior}}.",
        i = "{.fn sbc} calibrates one prior. {.fn prior_check} is what \\
             compares prior sets."
      ),
      call = call
    )
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L)) {
    cli::cli_abort(
      "{.arg seed} must be a single number or {.code NULL}.",
      call = call
    )
  }
  if (!is.null(fitter) && !is.function(fitter)) {
    cli::cli_abort(
      "{.arg .fitter} must be a function, \\
       not {.obj_type_friendly {fitter}}.",
      call = call
    )
  }
  if (identical(cache_mode, "results") && is.null(cache_location)) {
    cli::cli_abort(
      c(
        "{.arg cache_location} is needed when \\
         {.code cache_mode = \"results\"}.",
        i = "It is the directory {.pkg SBC} writes each fit's result to."
      ),
      call = call
    )
  }
  invisible(NULL)
}

#' Check that a model's implementation is calibrated
#'
#' Simulation-based calibration (Talts et al. 2018; Modrák et al. 2023)
#' is the check that a model's likelihood, its Stan code and its
#' post-processing agree with each other. If they do, the rank of a
#' parameter drawn from the prior, among the posterior draws of a fit to
#' data simulated from that draw, is uniform over the simulations. A
#' rank histogram that is not flat says the implementation is wrong; a
#' flat one says nothing about whether the model is a good one.
#'
#' `sbc()` computes no ranks. It fits the prior once, turns its draws
#' into data sets through [simulate_recovery()], builds the backend that
#' fits each of them, and hands both to `SBC::compute_SBC()`. What comes
#' back is SBC's own `SBC_results`, so `SBC::plot_rank_hist()`,
#' `SBC::plot_ecdf_diff()`, `SBC::plot_coverage()` and `results$stats`
#' all work as SBC documents them.
#'
#' @param model A `bmmodel`, built with the column names `data` uses.
#' @param formula The `bmmformula` under check. Without a `generator`,
#'   every parameter formula must be intercept-only, with or without one
#'   `(1 | id)` or `(1 | p | id)` term --- the shapes
#'   [recovery_formula()] writes. Anything else is an error naming the
#'   term: the default generator maps a prior draw onto
#'   [simulate_recovery()]'s arguments, and that map is exact for those
#'   shapes and guesswork for a covariate or a task factor. With a
#'   `generator` any formula works, as long as every parameter that has
#'   a group term has the same one, written as a bare column name.
#' @param data The design to simulate over, and only that: the subjects
#'   are the unique values of the grouping column and the trials are the
#'   rows each of them has. The response values are ignored, because the
#'   generator replaces them. Without a `generator`, every subject must
#'   have the same number of rows. The simulated data sets carry the
#'   same subject labels as `data`, which is what lets a subject truth
#'   meet its own `r_` draw. With a `generator`, `data` is handed to it
#'   as it is, columns and all.
#' @param prior A `brmsprior`, or `NULL` for bmm's defaults. This is the
#'   prior that is calibrated, so `NULL` calibrates bmm's own. A list is
#'   an error: [prior_check()] is what compares prior sets.
#' @param n_sims How many data sets to simulate and fit. Twenty is enough
#'   to see a run through; a verdict needs enough ranks to read a
#'   histogram over about twenty bins, which is what the default is for.
#'   **It is also `n_sims` model fits**, so a default run is hours of
#'   Stan and belongs in a script rather than at a prompt.
#' @param level Which draws are ranked: `"population"` the
#'   `b_<par>_Intercept` draws, `"sd"` the `sd_<group>__<par>_Intercept`
#'   draws, `"cor"` the `cor_<group>__…` draws and `"subject"` the
#'   `r_<group>__<par>[<id>,Intercept]` draws, one per subject and
#'   varying parameter. `"population"` cannot be dropped. A level the
#'   formula cannot produce --- `"sd"` or `"subject"` without a group
#'   term, `"cor"` without a correlated one --- is dropped with a
#'   message.
#'
#'   `"subject"` is off by default because it adds `n_subjects` times the
#'   number of varying parameters to the variables SBC ranks and plots:
#'   twenty subjects and two parameters are forty more rank histograms.
#'   It is the level that checks the partial pooling, and the one where a
#'   units error would show --- see Details.
#'
#'   `level` selects what is **ranked**, not what is **drawn**:
#'   everything the formula implies is always drawn from the prior and
#'   simulated from, because the ranks are only uniform when the data
#'   come from the joint prior. Holding the between-subject SDs at zero
#'   while the fitted model has a prior on them would take the
#'   population ranks down with it.
#'
#'   With a `generator`, each level means **every** draw of that class
#'   for the model's free parameters, read off the prior fit rather than
#'   built from the formula: `"population"` every `b_<par>_*` (so a
#'   cell-means formula's `b_kappa_task1` and a covariate's
#'   `b_kappa_cond`), `"sd"` every `sd_<group>__<par>_*`, `"cor"` every
#'   `cor_<group>__*` and `"subject"` every `r_<group>__<par>[*]`. A
#'   level the prior fit has no draw of is dropped with a message.
#' @param generator `NULL` simulates each data set with
#'   [simulate_recovery()], which is exact for the intercept-only
#'   formulas above and requires them. A function lifts that
#'   restriction: it is called once per data set as
#'   `generator(draws, data)`, where `draws` is that data set's prior
#'   draw as a **named numeric vector under the fit's own draw names**
#'   --- every `b_`, `sd_`, `cor_` and `r_` draw of the free parameters,
#'   whether or not it is ranked, because a generator handed less than
#'   the fit puts a prior on would simulate from a point mass where the
#'   fit has a prior --- and `data` is `data` as given. It returns the
#'   data frame to fit, which must carry the model's response columns
#'   and, when the formula has a group term, the grouping column with
#'   the labels of `data` if the subject level is ranked. bmmtools never
#'   interprets the row, so the truth names equal the draw names by
#'   construction. Note that brms's `r_` draws are deviations from the
#'   intercept: a generator that builds a subject's value adds them to
#'   the population value itself.
#' @param ... Passed to the fitter, for the prior fit and for every data
#'   set fit: `chains`, `iter`, `backend`, `init`, `control`. `cores`,
#'   `sample_prior`, `file`, `file_refit` and `file_compress` are
#'   refused, each with a message saying where it belongs.
#' @param seed Applied around the prior-draw subsample and around
#'   `SBC::generate_datasets()`, so the same seed gives the same data
#'   sets, and passed to the fitter of the prior fit, where it enters
#'   [fit_cached()]'s key. `NULL` leaves the random number generator
#'   alone and is recorded as `NA`. The data set fits take no seed:
#'   SBC's `future.seed` handles them, and one seed across fits would
#'   correlate them.
#' @param file Where to cache the **prior fit**, as in [prior_check()].
#'   `NULL` uses a temporary file. The data set fits are cached by SBC
#'   through `cache_mode`, not by [fit_cached()].
#' @param refit Passed to [fit_cached()] for the prior fit.
#' @param cores_per_fit Cores for each data set fit. `NULL` leaves
#'   `SBC::compute_SBC()`'s own default.
#' @param thin_ranks Thinning before ranking. `NULL` leaves SBC's default
#'   for a function backend.
#' @param keep_fits `TRUE` keeps every fit in the result, which is
#'   hundreds of MB for a default run. `FALSE` still keeps the ranks and
#'   the convergence diagnostics.
#' @param cache_mode,cache_location Passed to `SBC::compute_SBC()`.
#'   `"results"` needs a `cache_location`.
#' @param .fitter The fitting function, `bmm::bmm()` by default. Used for
#'   the prior fit and inside the backend; tests inject a stand-in so
#'   that nothing is compiled.
#'
#' @return `SBC::compute_SBC()`'s `SBC_results`, unchanged apart from one
#'   attribute, `bmmtools_sbc`, holding `model`, `n_sims`, `level`,
#'   `variables` (the draw names that were ranked), `seed` (`NA` when
#'   none was given), `prior` (`brms::prior_summary()` of the prior fit,
#'   which is what the draws came from), `layout`, `diagnostics` and
#'   `generator` (`"simulate_recovery"` or `"user"`).
#'
#' @details
#' Two convergence bars apply to the same fits. bmmtools warns once on
#' its own Rhat 1.05 (see [check_convergence()]) together with SBC's
#' rank-ESS 0.5, records both counts in the attribute, and names SBC's
#' stricter Rhat 1.01 in the message. An ESS below half the maximum rank
#' skews the ranks themselves, which is why it is a bar here and not
#' only a diagnostic.
#'
#' The subject truths are **deviations**. brms's `r_` draws are
#' deviations from the intercept, while [simulate_recovery()]'s
#' `truth$subjects` are absolute values, so the generator emits each
#' subject's value minus that data set's population value, on the link
#' scale, under the label the subject has in `data`. The subject values
#' themselves are not read off the prior draw: they are drawn by
#' [simulate_recovery()] from the drawn SDs and correlations, which is
#' the same conditional the prior puts on them.
#'
#' @references
#' Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A.
#' (2018). Validating Bayesian inference algorithms with simulation-based
#' calibration. \doi{10.48550/arXiv.1804.06788}
#'
#' Modrák, M., Moon, A. H., Kim, S., Bürkner, P.-C., Huurre, N.,
#' Faltejsková, K., Gelman, A., & Vehtari, A. (2023). Simulation-based
#' calibration checking for Bayesian computation: The choice of test
#' quantities shapes sensitivity. *Bayesian Analysis*.
#' \doi{10.1214/23-BA1404}
#'
#' @examples
#' \dontrun{
#' model <- bmm::mixture2p(resp_error = "y")
#' results <- sbc(
#'   model,
#'   recovery_formula(model),
#'   data.frame(id = rep(1:20, each = 50), y = 0),
#'   n_sims = 20,
#'   chains = 2, iter = 500, backend = "cmdstanr",
#'   seed = 1
#' )
#' SBC::plot_rank_hist(results)
#' attr(results, "bmmtools_sbc")$variables
#'
#' # any formula, with a generator that reads the draw row itself: here a
#' # condition effect on kappa, so `draws` carries `b_kappa_cond` too
#' design <- data.frame(
#'   id = rep(1:20, each = 50), cond = rep(c(-0.5, 0.5), 500), y = 0
#' )
#' generate <- function(draws, data) {
#'   # one call per subject and condition: rmixture2p() takes scalars
#'   cells <- split(seq_len(nrow(data)), list(data$id, data$cond))
#'   for (rows in cells) {
#'     id <- data$id[[rows[[1]]]]
#'     cond <- data$cond[[rows[[1]]]]
#'     kappa <- exp(
#'       draws[["b_kappa_Intercept"]] + draws[["b_kappa_cond"]] * cond +
#'         draws[[paste0("r_id__kappa[", id, ",Intercept]")]]
#'     )
#'     p_mem <- plogis(
#'       draws[["b_thetat_Intercept"]] +
#'         draws[[paste0("r_id__thetat[", id, ",Intercept]")]]
#'     )
#'     data$y[rows] <- bmm::rmixture2p(
#'       length(rows), kappa = kappa, p_mem = p_mem
#'     )
#'   }
#'   data
#' }
#' results <- sbc(
#'   model,
#'   bmm::bmf(kappa ~ 1 + cond + (1 | id), thetat ~ 1 + (1 | id)),
#'   design,
#'   generator = generate,
#'   n_sims = 20,
#'   chains = 2, iter = 500, backend = "cmdstanr"
#' )
#' }
#'
#' @export
sbc <- function(model,
                formula,
                data,
                prior = NULL,
                n_sims = 100,
                level = c("population", "sd"),
                generator = NULL,
                ...,
                seed = NULL,
                file = NULL,
                refit = c("on_change", "never", "always"),
                cores_per_fit = NULL,
                thin_ranks = NULL,
                keep_fits = FALSE,
                cache_mode = c("none", "results"),
                cache_location = NULL,
                .fitter = NULL) {
  rlang::check_installed(
    c("SBC", "posterior"), "to run simulation-based calibration."
  )
  check_model(model)
  level <- rlang::arg_match(
    level, c("population", "sd", "cor", "subject"),
    multiple = TRUE
  )
  refit <- rlang::arg_match(refit)
  cache_mode <- rlang::arg_match(cache_mode)
  n_sims <- check_count(n_sims, "n_sims")
  if (!"population" %in% level) {
    cli::cli_abort(c(
      "{.arg level} must include {.val population}.",
      i = "The population parameters are the ones every model has, and a \\
           calibration that ranks none of them checks nothing."
    ))
  }
  if (!is.data.frame(data)) {
    cli::cli_abort(
      "{.arg data} must be a data frame, not {.obj_type_friendly {data}}."
    )
  }
  check_sbc_args(prior, seed, .fitter, cache_mode, cache_location)
  check_sbc_generator(generator)
  dots <- rlang::list2(...)
  check_sbc_dots(dots)

  # with a user generator the formula may have any shape, so only its
  # group term is read; the default generator needs the whole check
  user <- !is.null(generator)
  groups <- if (user) {
    sbc_formula_groups(formula)
  } else {
    check_sbc_formula(formula)
  }
  group <- groups$group %||% "id"
  check_group_ids(data, group)
  layout <- sbc_layout(data, group, balanced = !user)
  sets <- if (!user) sbc_pattern_sets(model, groups, group, level)

  check_improper_priors(formula, data, model, list(sbc = prior))
  base <- file %||% tempfile(pattern = "bmmtools-sbc-")
  fit <- prior_fit(
    "sbc", prior, formula, data, model, base, refit, seed, dots, .fitter
  )
  if (user) {
    # every draw first, because which of them are the four classes is
    # read off the fit's own names, not built from the formula
    draws <- prior_parameter_draws(fit, n_sims, c(all = "."), seed = seed)
    sets <- sbc_generator_sets(
      model, groups, group, level, posterior::variables(draws)
    )
    draws <- posterior::subset_draws(draws, variable = sets$take)
    ranked <- sets$rank
  } else {
    draws <- prior_parameter_draws(fit, n_sims, sets$take, seed = seed)
    ranked <- select_variables(sets$rank, posterior::variables(draws))
  }

  fitter <- .fitter
  if (is.null(fitter)) {
    rlang::check_installed("bmm", "to fit the simulated data sets.")
    fitter <- bmm::bmm
  }
  sbc_generator_object <- if (user) {
    sbc_user_generator(
      draws, ranked, generator, data, model, layout, group,
      needs_group = !is.null(groups$group)
    )
  } else {
    sbc_generator(draws, ranked, model, layout, groups$correlated, group)
  }
  datasets <- with_seed_if(
    seed, SBC::generate_datasets(sbc_generator_object, n_sims)
  )

  options <- list(
    datasets = datasets,
    backend = sbc_backend(formula, model, prior, dots, fitter),
    keep_fits = keep_fits,
    cache_mode = cache_mode
  )
  # each NULL option is left out so that SBC's own default applies
  options$cores_per_fit <- cores_per_fit
  options$thin_ranks <- thin_ranks
  options$cache_location <- cache_location
  results <- rlang::exec(SBC::compute_SBC, !!!options)

  attr(results, "bmmtools_sbc") <- list(
    model = class(model)[[length(class(model))]],
    n_sims = n_sims,
    level = sets$level,
    variables = ranked,
    seed = if (is.null(seed)) NA_real_ else as.double(seed),
    prior = tryCatch(brms::prior_summary(fit), error = function(e) NULL),
    layout = layout,
    diagnostics = sbc_diagnostics(results),
    generator = if (user) "user" else "simulate_recovery"
  )
  results
}
