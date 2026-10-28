# validation_plan(): a declarative plan of one model's validation
# (Milestone 26, D79-D85).
#
# The plan fits nothing. It names the truth, the grid, the budget, the
# SBC design point and its prior, the prior-check layout, the reference,
# the criteria, the seed and the owner, and it lists the stages with the
# status each has before anything runs. It can be printed, written next
# to a card and diffed against an earlier version; run_validation()
# executes it.

#' Plan the validation of a bmm model
#'
#' `r lifecycle::badge("experimental")`
#'
#' Builds a `bmmtools_validation_plan`: everything [run_validation()]
#' needs to validate one model, resolved before anything is fitted. Every
#' argument left `NULL` is filled from the model's [validation_defaults()]
#' entry, or for a model without one from the generic grid, in which case
#' `pars` is required.
#'
#' @param model A `bmmodel`, or the name of a model with a
#'   [validation_defaults()] entry, which is then built with the entry's
#'   variable names (bmm is needed for that).
#' @param pars Named population values on the link scale, as in
#'   [simulate_recovery()].
#' @param sds Named between-subject SDs on the link scale; `NULL` takes the
#'   entry's.
#' @param cors A correlation matrix between subject values, as in
#'   [simulate_recovery()]; `NULL` takes the entry's, which may be none.
#' @param grid A data frame with the columns `n_subjects` and `n_trials`,
#'   one row per design.
#' @param reps Replications per design.
#' @param sampler A named list passed to the fitter of every recovery cell
#'   (`chains`, `iter`, `warmup`, `cores`, `backend`, ...). A list given
#'   here changes the entries it names and keeps the others. It may not
#'   set `prior` or an abbreviation of it (`pri`), nor may `sbc$sampler`
#'   or `prior_check$sampler`: the calibration's prior is `sbc$prior`, and
#'   the recovery stage fits under bmm's default priors.
#' @param sbc A named list for the calibration: `n_subjects`, `n_trials`,
#'   `n_sims`, `prior` (a `brmsprior`, a data frame with the columns of
#'   `brms::set_prior()`, or `NULL` for bmm's defaults), `prior_reason`
#'   (why it is not bmm's default, printed on the card), `level`,
#'   `sampler`, `generator`, `cores_per_fit` and `thin_ranks` (see
#'   [sbc()]). Named entries replace the entry's, the rest are kept,
#'   except that a `prior` other than the entry's, given without a
#'   `prior_reason`, drops the entry's reason. The prior is compared as
#'   given, column order and types included: a `brmsprior` with the same
#'   rows as an entry's data frame, the entry's columns in another order,
#'   or a factor column counts as another prior. A `prior_reason` without
#'   a `prior` is an error.
#' @param prior_check A named list for the prior predictive check:
#'   `n_subjects`, `n_trials`, `n_draws`, `range`, `range_reason` (why the
#'   range is what it is, printed on the card), `summary` and `sampler`
#'   (see [prior_check()]). Merged as `sbc` is: a `range` other than the
#'   entry's, given without a `range_reason`, drops the entry's reason,
#'   and a `range_reason` without a `range` is an error. The range is
#'   compared as given, names included, so `c(lo = 0.15, hi = 5)` counts
#'   as another range than the entry's `c(0.15, 5)`.
#' @param reference `NULL`, or a reference implementation's estimates for
#'   the agreement stage: a data frame with `term` and `estimate` (and
#'   optionally `ci_low`, `ci_high`) on the natural scale, as
#'   [cross_check()] takes it, plus the provenance columns `package`,
#'   `version`, `call`, `data_hash`, `date` and `source`.
#' @param convergence `NULL` for [check_convergence()]'s defaults, or a
#'   named list of its arguments, as in [recovery_grid()].
#' @param criteria A named list replacing some of the default criteria
#'   (see "Criteria").
#' @param stages `NULL` for every stage that exists, or the names of the
#'   stages to run (see "Stages"). The run record always runs; a stage
#'   that reads the recovery grid needs `"recovery"` too.
#' @param formula,trial_design,generator As in [recovery_grid()]. `NULL`
#'   `formula` means [recovery_formula()] with correlated subject effects
#'   when there are at least two subject-varying parameters, so the
#'   individual-differences stage has correlations to score.
#' @param seed The master seed, required. Each stage derives its own from
#'   it (see "Seeds").
#' @param owner,version The entry's owner and plan version; `NULL` takes
#'   the entry's.
#'
#' @return A `bmmtools_validation_plan`: a named list with every field
#'   above resolved, plus `model_name`, `defaults` (`"bmmtools <version>
#'   entry"` or `"none"`), `re_cor`, `correlations`, `seeds`, `stages` (a
#'   tibble: `stage`, `order`, `title`, `required`, `status`, `issue`,
#'   `reason`), `hash` and `bmmtools_version`.
#'
#' @section Stages:
#' In card order: `run_record`, `prior_predictive`, `sbc`, `recovery`,
#' `computation`, `identifiability`, `individual_differences`,
#' `data_requirements`, `agreement`, `model_recovery` and
#' `posterior_checks`. A stage is `"planned"`; `"not built"` with the
#' issue that will build it (model recovery, #6; posterior checks, #20);
#' `"not applicable"` with the reason (individual differences with fewer
#' than two subject-varying parameters; data requirements without a
#' subject-varying parameter, or without two numbers of trials at
#' `read_off_n_subjects` subjects in the grid; agreement without a
#' reference); or `"excluded"` when
#' `stages` leaves it out. All but agreement, model recovery and posterior
#' checks are required for a card to be called complete.
#'
#' @section Criteria:
#' `prior_bound_rate_max` (0.10): the largest floor or ceiling rate of the
#' prior predictive draws under the calibration's prior; under bmm's
#' defaults the rates are reported, not judged, unless the calibration
#' uses them. `sbc_ecdf_prob` (0.95): the simultaneous band of SBC's
#' ECDF-difference check no ranked variable may leave.
#' `convergence_rate_min` (0.95): the share of recovery cells that pass
#' [check_convergence()]. `coverage_level` (0.95) and `coverage_alpha`
#' (0.05): the coverage of every population parameter, its intervals
#' pooled over the designs, must not fall below `coverage_level` by a
#' one-sided binomial test at `coverage_alpha` divided by the number of
#' parameters; the generating correlations are held to the same test.
#' `contraction_min` (0) and `z_mean_max` (2): the identifiability bars at
#' the largest design.
#' `read_off_r_low` (0.80) and `read_off_n_subjects` (50): the data
#' requirements read-off, the smallest `n_trials` at which every
#' subject-varying parameter's person-level `r_low` reaches the bar at
#' that many subjects; 0.80 is a convention, and the read-off at 0.70 and
#' 0.90 is reported beside it.
#' `tolerance`: by metric, how far a number without a Monte Carlo
#' standard error may move between two runs before
#' [compare_validations()] flags it.
#'
#' @section Seeds:
#' The recovery grid runs under `seed`, the prior check under `seed + 1`
#' and the calibration under `seed + 2`. The seeds, the owner and the
#' reasons (`sbc$prior_reason`, `prior_check$range_reason`) are left out
#' of the plan's hash and of the stages' keys, so the hash changes only
#' when what is planned changes, and rewording a reason reruns nothing.
#'
#' @seealso [run_validation()], [validation_defaults()].
#'
#' @examples
#' plan <- validation_plan(
#'   structure(
#'     list(links = list(kappa = "log", thetat = "logit")),
#'     class = c("bmmodel", "circular", "mixture2p")
#'   ),
#'   seed = 1
#' )
#' plan
#' \dontrun{
#' plan <- validation_plan("ezdm", seed = 20261024)
#' }
#'
#' @export
validation_plan <- function(model,
                            pars = NULL,
                            sds = NULL,
                            cors = NULL,
                            grid = NULL,
                            reps = NULL,
                            sampler = NULL,
                            sbc = NULL,
                            prior_check = NULL,
                            reference = NULL,
                            convergence = NULL,
                            criteria = NULL,
                            stages = NULL,
                            formula = NULL,
                            trial_design = NULL,
                            generator = NULL,
                            seed,
                            owner = NULL,
                            version = NULL) {
  if (missing(seed)) {
    cli::cli_abort(c(
      "{.arg seed} is required.",
      i = "A plan is a record of what was run; every stage derives its \\
           seed from this one."
    ))
  }
  seed <- check_plan_seed(seed)
  name <- validation_model_name(model)
  entry <- validation_defaults(name)
  model <- resolve_plan_model(model, name, entry)
  defaults <- if (is.null(entry)) {
    "none"
  } else {
    paste("bmmtools", package_version_string("bmmtools"), "entry")
  }
  base <- entry %||% generic_validation_defaults()

  if (is.null(pars)) {
    if (is.null(entry)) {
      cli::cli_abort(c(
        "{.arg pars} is required: {.val {name}} has no validation defaults.",
        i = "{.fn validation_defaults} lists the models that have them."
      ))
    }
    pars <- entry$pars
  }
  pars <- check_plan_numbers(pars, "pars")
  sds <- if (is.null(sds)) base$sds else check_plan_numbers(sds, "sds")
  cors <- cors %||% base$cors
  grid <- check_plan_grid(grid %||% base$grid)
  reps <- check_count(reps %||% base$reps, "reps")
  sampler <- merge_plan_list(base$sampler, sampler, "sampler")
  given_sbc <- names(sbc)
  sbc <- merge_plan_list(base$sbc, sbc, "sbc")
  given_pc <- names(prior_check)
  prior_check <- merge_plan_list(base$prior_check, prior_check, "prior_check")
  # the card prints the reasons, so a prior or range other than the
  # entry's, given without a reason, does not inherit the entry's, which
  # explains a prior that did not run or a range that is not judged; the
  # entry's own given again keeps it
  sbc <- drop_inherited_reason(sbc, base$sbc, given_sbc, "prior")
  prior_check <- drop_inherited_reason(
    prior_check, base$prior_check, given_pc, "range"
  )
  sbc["prior_reason"] <- list(
    check_plan_reason(sbc[["prior_reason"]], "sbc$prior_reason")
  )
  prior_check["range_reason"] <- list(check_plan_reason(
    prior_check[["range_reason"]], "prior_check$range_reason"
  ))
  if (is.null(sbc[["prior"]]) && !is.null(sbc[["prior_reason"]])) {
    cli::cli_abort(c(
      "{.arg sbc} gives a {.field prior_reason} but no {.field prior}.",
      i = "The reason says why the calibration's prior is not bmm's \\
           default, and the card prints it; without a prior the \\
           calibration runs under bmm's defaults. Give the prior, or \\
           leave out the reason."
    ))
  }
  no_range <- is.null(prior_check[["range"]])
  if (no_range && !is.null(prior_check[["range_reason"]])) {
    cli::cli_abort(c(
      "{.arg prior_check} gives a {.field range_reason} but no \\
       {.field range}.",
      i = "The reason says why the range is what it is, and the card \\
           prints it under the prior check. Give the range, or leave out \\
           the reason."
    ))
  }
  sbc$sampler <- merge_plan_list(
    base$sbc$sampler, sbc$sampler, "sbc$sampler"
  )
  check_plan_sampler_prior(sampler, "sampler")
  check_plan_sampler_prior(sbc$sampler, "sbc$sampler")
  check_plan_sampler_prior(prior_check$sampler, "prior_check$sampler")
  check_plan_design(sbc, "sbc")
  check_plan_design(prior_check, "prior_check")
  sbc$n_sims <- check_count(sbc$n_sims, "sbc$n_sims")
  prior_check$n_draws <- check_count(prior_check$n_draws, "prior_check$n_draws")
  check_plan_prior(sbc[["prior"]])
  reference <- check_plan_reference(reference %||% base$reference)
  convergence <- check_convergence_arg(convergence)
  criteria <- merge_plan_criteria(criteria)
  owner <- check_plan_string(owner %||% base$owner, "owner")
  version <- check_plan_string(version %||% base$version, "version")

  varying <- names(sds)[!is.na(sds) & sds > 0]
  individual <- length(varying) >= 2L
  requirements <- data_requirements_reason(varying, grid, criteria)
  re_cor <- if (individual) "all" else "none"
  correlations <- if (individual) "model" else NULL

  plan <- list(
    model = model,
    model_name = name,
    defaults = defaults,
    pars = pars,
    sds = sds,
    cors = cors,
    grid = grid,
    reps = reps,
    sampler = sampler,
    sbc = sbc,
    prior_check = prior_check,
    reference = reference,
    convergence = convergence,
    criteria = criteria,
    formula = formula,
    trial_design = trial_design,
    generator = generator,
    re_cor = re_cor,
    correlations = correlations,
    sections = base$sections %||% list(),
    seed = seed,
    seeds = c(recovery = seed, prior_predictive = seed + 1, sbc = seed + 2),
    owner = owner,
    version = version,
    stages = plan_stages(stages, individual, reference, requirements),
    bmmtools_version = package_version_string("bmmtools")
  )
  plan$hash <- plan_hash(plan)
  structure(plan, class = "bmmtools_validation_plan")
}

#' The bmm call that builds a model of this name, for a hint
#'
#' `bmm::<name>(...)` when bmm exports `name`; for a version class such as
#' `ezdm_4par`, the constructor with its `version` argument, read from the
#' constructor's own choices; `NULL` when bmm is not installed or neither
#' applies.
#' @noRd
plan_constructor_call <- function(name) {
  if (!requireNamespace("bmm", quietly = TRUE)) {
    return(NULL)
  }
  exported <- function(x) {
    fn <- tryCatch(getExportedValue("bmm", x), error = function(e) NULL)
    if (is.function(fn)) fn else NULL
  }
  if (!is.null(exported(name))) {
    return(paste0("bmm::", name, "(...)"))
  }
  parts <- regmatches(name, regexec("^(.+)_([^_]+)$", name))[[1L]]
  if (length(parts) != 3L) {
    return(NULL)
  }
  constructor <- exported(parts[[2L]])
  if (is.null(constructor) || !"version" %in% names(formals(constructor))) {
    return(NULL)
  }
  versions <- tryCatch(
    eval(formals(constructor)$version, environment(constructor)),
    error = function(e) NULL
  )
  if (!is.character(versions) || !parts[[3L]] %in% versions) {
    return(NULL)
  }
  paste0("bmm::", parts[[2L]], "(..., version = \"", parts[[3L]], "\")")
}

#' Build a named model from its entry
#' @noRd
resolve_plan_model <- function(model, name, entry,
                               call = rlang::caller_env()) {
  if (inherits(model, "bmmodel")) {
    return(model)
  }
  if (is.null(entry)) {
    constructor_call <- plan_constructor_call(name)
    hint <- if (is.null(constructor_call)) {
      "Pass the model object with {.arg pars}."
    } else {
      "Pass the model object, such as {.code {constructor_call}}, with \\
       {.arg pars}."
    }
    cli::cli_abort(
      c(
        "{.val {name}} has no validation defaults, so it cannot be built \\
         from its name.",
        i = hint
      ),
      call = call
    )
  }
  rlang::check_installed("bmm", "to build a model from its name.")
  constructor <- tryCatch(
    getExportedValue("bmm", name),
    error = function(e) NULL
  )
  if (!is.function(constructor)) {
    cli::cli_abort(
      "bmm {package_version_string('bmm')} has no model {.val {name}}.",
      call = call
    )
  }
  # a call to bmm::<name>(), not exec() on the function: bmm's check for
  # missing arguments reads the function's name off sys.call(), which a
  # function object does not have (measured on bmm 1.4.1, `ezdm`)
  eval(rlang::call2(name, !!!entry$model_args, .ns = "bmm"))
}

#' @noRd
check_plan_seed <- function(seed, call = rlang::caller_env()) {
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
    cli::cli_abort(
      "{.arg seed} must be a single finite number, \\
       not {.obj_type_friendly {seed}}.",
      call = call
    )
  }
  as.double(seed)
}

#' @noRd
check_plan_numbers <- function(x, arg, call = rlang::caller_env()) {
  bad <- !is.numeric(x) || length(x) == 0L || !rlang::is_named(x) ||
    anyNA(x) || anyDuplicated(names(x)) > 0L
  if (bad) {
    cli::cli_abort(
      "{.arg {arg}} must be a named numeric vector without missing values \\
       or repeated names.",
      call = call
    )
  }
  x
}

#' @noRd
check_plan_grid <- function(grid, call = rlang::caller_env()) {
  needed <- c("n_subjects", "n_trials")
  bad <- !is.data.frame(grid) || !all(needed %in% names(grid)) ||
    nrow(grid) == 0L
  if (bad) {
    cli::cli_abort(
      c(
        "{.arg grid} must be a data frame with the columns \\
         {.field n_subjects} and {.field n_trials} and at least one row.",
        i = "One row per design, such as \\
             {.code expand.grid(n_subjects = c(20, 50), n_trials = 100)}."
      ),
      call = call
    )
  }
  for (column in needed) {
    values <- grid[[column]]
    bad <- !is.numeric(values) || anyNA(values) || any(values < 1) ||
      any(values != round(values))
    if (bad) {
      cli::cli_abort(
        "{.field {column}} in {.arg grid} must hold positive whole numbers.",
        call = call
      )
    }
  }
  data.frame(
    n_subjects = as.integer(grid$n_subjects),
    n_trials = as.integer(grid$n_trials)
  )
}

#' Replace the entries of a defaults list that a caller names
#'
#' `[<-` with `list()`, so that an explicit `NULL` (bmm's default prior,
#' say) replaces the entry rather than deleting it.
#'
#' @noRd
merge_plan_list <- function(base, given, arg, call = rlang::caller_env()) {
  if (is.null(given)) {
    return(base)
  }
  bad <- !is.list(given) || is.data.frame(given) ||
    (length(given) > 0L && !rlang::is_named(given))
  if (bad) {
    cli::cli_abort(
      "{.arg {arg}} must be a named list or {.code NULL}, \\
       not {.obj_type_friendly {given}}.",
      call = call
    )
  }
  out <- base %||% list()
  for (name in names(given)) out[name] <- list(given[[name]])
  out
}

#' Drop the entry's reason for a field the call changed without one
#'
#' `field` is `"prior"` or `"range"`, its reason `<field>_reason`; the
#' field counts as changed when it is given and differs from the entry's
#' as `plan_key()` stores it.
#' @noRd
drop_inherited_reason <- function(merged, base, given, field) {
  reason <- paste0(field, "_reason")
  changed <- field %in% given &&
    !identical(plan_key(merged[[field]]), plan_key(base[[field]]))
  if (changed && !reason %in% given) merged[reason] <- list(NULL)
  merged
}

#' A sampler list may not carry a prior
#'
#' In `sampler` it would reach every recovery fit, which the card says fit
#' under bmm's defaults; in the calibration's or the prior check's sampler
#' it collides with the prior the stage passes itself. An abbreviation
#' (`pri`) reaches `bmm::bmm()`'s `prior` by partial matching.
#' @noRd
check_plan_sampler_prior <- function(x, arg, call = rlang::caller_env()) {
  offending <- prior_abbreviations(names(x))
  if (length(offending) == 0L) {
    return(invisible(x))
  }
  cli::cli_abort(
    c(
      "{.arg {arg}} may not set {.arg prior} or an abbreviation of it: \\
       {.val {offending}}.",
      x = validation_refused_dots()[["prior"]]
    ),
    call = call
  )
}

#' A reason printed on the card: `NULL`, or one string with something in it
#' @noRd
check_plan_reason <- function(x, arg, call = rlang::caller_env()) {
  if (is.null(x)) {
    return(NULL)
  }
  string <- is.character(x) && length(x) == 1L && !is.na(x)
  if (!string || !nzchar(trimws(x))) {
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    what <- if (string) {
      "a blank string"
    } else {
      cli::format_inline("{.obj_type_friendly {x}}")
    }
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a single non-empty string or {.code NULL}, \\
         not {what}.",
        i = "The card prints it as a sentence; {.code NULL} prints none."
      ),
      call = call
    )
  }
  x
}

#' @noRd
check_plan_design <- function(x, arg, call = rlang::caller_env()) {
  for (field in c("n_subjects", "n_trials")) {
    check_count(x[[field]], paste0(arg, "$", field), call = call)
  }
  invisible(x)
}

#' @noRd
check_plan_prior <- function(prior, call = rlang::caller_env()) {
  if (is.null(prior) || inherits(prior, "brmsprior")) {
    return(invisible(prior))
  }
  if (!is.data.frame(prior) || !all(c("prior", "class") %in% names(prior))) {
    cli::cli_abort(
      c(
        "{.arg sbc$prior} must be a {.cls brmsprior}, a data frame with \\
         the columns of {.fn brms::set_prior}, or {.code NULL}.",
        i = "A data frame needs at least {.field prior} and \\
             {.field class}."
      ),
      call = call
    )
  }
  invisible(prior)
}

#' The provenance columns a reference table carries (D85)
#' @noRd
reference_provenance_columns <- function() {
  c("package", "version", "call", "data_hash", "date", "source")
}

#' @noRd
check_plan_reference <- function(reference, call = rlang::caller_env()) {
  if (is.null(reference)) {
    return(NULL)
  }
  needed <- c("term", "estimate", reference_provenance_columns())
  if (!is.data.frame(reference)) {
    cli::cli_abort(
      "{.arg reference} must be a data frame or {.code NULL}, \\
       not {.obj_type_friendly {reference}}.",
      call = call
    )
  }
  missing <- setdiff(needed, names(reference))
  if (length(missing) > 0L) {
    cli::cli_abort(
      c(
        "{.arg reference} lacks the column{?s} {.val {missing}}.",
        i = "A reference names where its numbers come from: \\
             {.val {reference_provenance_columns()}}."
      ),
      call = call
    )
  }
  tibble::as_tibble(reference)
}

#' @noRd
merge_plan_criteria <- function(criteria, call = rlang::caller_env()) {
  out <- default_validation_criteria()
  if (is.null(criteria)) {
    return(out)
  }
  if (!is.list(criteria) || !rlang::is_named(criteria)) {
    cli::cli_abort(
      "{.arg criteria} must be a named list or {.code NULL}.",
      call = call
    )
  }
  unknown <- setdiff(names(criteria), names(out))
  if (length(unknown) > 0L) {
    cli::cli_abort(
      c(
        "{.arg criteria} has no criterion {.val {unknown}}.",
        i = "The criteria are {.val {names(out)}}."
      ),
      call = call
    )
  }
  if ("tolerance" %in% names(criteria)) {
    tolerance <- criteria$tolerance
    if (!is.numeric(tolerance) || !rlang::is_named(tolerance)) {
      cli::cli_abort(
        "{.code criteria$tolerance} must be a named numeric vector.",
        call = call
      )
    }
    out$tolerance[names(tolerance)] <- tolerance
    criteria$tolerance <- NULL
  }
  for (name in names(criteria)) {
    value <- criteria[[name]]
    if (!is.numeric(value) || length(value) != 1L || is.na(value)) {
      cli::cli_abort(
        "{.code criteria${name}} must be a single number.",
        call = call
      )
    }
    out[[name]] <- value
  }
  out
}

#' @noRd
check_plan_string <- function(x, arg, call = rlang::caller_env()) {
  if (is.null(x)) {
    return(NA_character_)
  }
  if (!is.character(x) || length(x) != 1L) {
    cli::cli_abort(
      "{.arg {arg}} must be a single string, not {.obj_type_friendly {x}}.",
      call = call
    )
  }
  x
}

#' Why the data requirements stage cannot run on a plan, or `NA`
#'
#' D83 reads the smallest `n_trials` off the person-level correlation at
#' `criteria$read_off_n_subjects`: it needs a subject-varying parameter,
#' that many subjects in the grid, and two numbers of trials there.
#' @noRd
data_requirements_reason <- function(varying, grid, criteria) {
  n <- criteria$read_off_n_subjects
  trials <- unique(grid$n_trials[grid$n_subjects == n])
  if (length(varying) == 0L) {
    "no subject-varying parameter"
  } else if (length(trials) == 0L) {
    sprintf("the grid has no design with %s subjects", format(n))
  } else if (length(trials) == 1L) {
    sprintf("the grid has one number of trials at %s subjects", format(n))
  } else {
    NA_character_
  }
}

#' The stage table of a plan, before anything runs
#' @noRd
plan_stages <- function(stages, individual, reference,
                        requirements = NA_character_,
                        call = rlang::caller_env()) {
  table <- validation_stage_table()
  known <- table$stage
  if (!is.null(stages)) {
    if (!is.character(stages) || anyNA(stages)) {
      cli::cli_abort(
        "{.arg stages} must be a character vector of stage names.",
        call = call
      )
    }
    unknown <- setdiff(stages, known)
    if (length(unknown) > 0L) {
      cli::cli_abort(
        c(
          "{.arg stages} names no stage {.val {unknown}}.",
          i = "The stages are {.val {known}}."
        ),
        call = call
      )
    }
    stages <- union("run_record", stages)
    readers <- table$stage[!is.na(table$reads) & table$built]
    orphans <- intersect(stages, readers)
    if (length(orphans) > 0L && !"recovery" %in% stages) {
      cli::cli_abort(
        c(
          "{.val {orphans}} read{?s/} the recovery grid, which {.arg stages} \\
           leaves out.",
          i = "Add {.val recovery} to {.arg stages}."
        ),
        call = call
      )
    }
  }

  status <- ifelse(table$built, "planned", "not built")
  reason <- ifelse(
    table$built, NA_character_, sprintf("not run (#%d)", table$issue)
  )
  if (!individual) {
    at <- table$stage == "individual_differences"
    status[at] <- "not applicable"
    reason[at] <- "fewer than two subject-varying parameters"
  }
  if (!is.na(requirements)) {
    at <- table$stage == "data_requirements"
    status[at] <- "not applicable"
    reason[at] <- requirements
  }
  if (is.null(reference)) {
    at <- table$stage == "agreement"
    status[at] <- "not applicable"
    reason[at] <- "no reference implementation"
  }
  if (!is.null(stages)) {
    at <- !table$stage %in% stages & status == "planned"
    status[at] <- "excluded"
    reason[at] <- "not in `stages`"
  }
  tibble::tibble(
    stage = table$stage,
    order = table$order,
    title = table$title,
    required = table$required,
    status = status,
    issue = table$issue,
    reason = reason
  )
}

#' What a plan's hash is computed over
#'
#' Formulas enter deparsed and functions as their text, as in the cache
#' key, so that an environment never changes the hash. The seeds, the
#' owner, the version string, the bmmtools version and the reasons are
#' left out (D84): the hash answers whether what is planned changed.
#'
#' @noRd
plan_hash_input <- function(plan) {
  excluded <- c(
    "seed", "seeds", "owner", "version", "bmmtools_version", "hash",
    "defaults"
  )
  fields <- unclass(plan)[setdiff(names(plan), excluded)]
  plan_key(without_plan_reasons(fields))
}

#' Plan fields without their reasons
#'
#' `sbc$prior_reason` and `prior_check$range_reason` are prose printed on
#' the card; they say why, not what is run, so neither the plan's hash nor
#' a stage's key reads them.
#'
#' @noRd
without_plan_reasons <- function(fields) {
  if (is.list(fields$sbc)) fields$sbc["prior_reason"] <- NULL
  if (is.list(fields$prior_check)) {
    fields$prior_check["range_reason"] <- NULL
  }
  fields
}

#' Reduce a plan field to what can be hashed reproducibly
#' @noRd
plan_key <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (inherits(x, "formula") || inherits(x, "bmmformula")) {
    return(formula_key(x))
  }
  if (is.function(x)) {
    return(function_key(x))
  }
  if (is.data.frame(x)) {
    return(lapply(as.list(x), plan_key))
  }
  if (is.list(x)) {
    out <- lapply(unclass(x), plan_key)
    names(out) <- names(x)
    return(out)
  }
  x
}

#' @noRd
plan_hash <- function(plan) {
  rlang::hash(plan_hash_input(plan))
}

#' Format and print a validation plan
#'
#' @param x A `bmmtools_validation_plan`.
#' @param ... Not used.
#' @return `format()` returns the lines; `print()` writes them and returns
#'   `x` invisibly.
#' @export
format.bmmtools_validation_plan <- function(x, ...) {
  numbers <- function(v) {
    if (is.null(v) || length(v) == 0L) {
      return("none")
    }
    paste(names(v), "=", signif(v, 3), collapse = ", ")
  }
  sampler <- function(s) {
    if (length(s) == 0L) {
      return("the fitter's defaults")
    }
    paste(names(s), "=", vapply(s, format_plan_value, ""), collapse = ", ")
  }
  stages <- x$stages
  stage_lines <- sprintf(
    "  %2d %-24s %-15s%s",
    stages$order, stages$stage, stages$status,
    ifelse(is.na(stages$reason), "", paste0(" ", stages$reason))
  )
  sbc_prior <- if (is.null(x$sbc[["prior"]])) {
    "bmm's default priors"
  } else {
    "the plan's own prior"
  }
  c(
    sprintf(
      "<bmmtools validation plan> %s, version %s, hash %s",
      x$model_name, format_plan_version(x$version), substr(x$hash, 1L, 8L)
    ),
    sprintf("Owner: %s; defaults: %s", na_or(x$owner, "none"), x$defaults),
    sprintf("Truth (link scale): %s", numbers(x$pars)),
    sprintf("Subject SDs (link scale): %s", numbers(x$sds)),
    sprintf(
      "Grid: %d design%s x %d replication%s (n_subjects %s; n_trials %s)",
      nrow(x$grid), if (nrow(x$grid) == 1L) "" else "s",
      x$reps, if (x$reps == 1L) "" else "s",
      paste(unique(x$grid$n_subjects), collapse = ", "),
      paste(unique(x$grid$n_trials), collapse = ", ")
    ),
    sprintf("Sampler: %s", sampler(x$sampler)),
    sprintf(
      "SBC: %d simulations at %d subjects x %d trials under %s",
      x$sbc$n_sims, x$sbc$n_subjects, x$sbc$n_trials, sbc_prior
    ),
    sprintf(
      "Prior check: %d draws at %d subjects x %d trials",
      x$prior_check$n_draws, x$prior_check$n_subjects,
      x$prior_check$n_trials
    ),
    sprintf(
      "Reference: %s",
      if (is.null(x$reference)) {
        "none"
      } else {
        paste(unique(x$reference$package), collapse = ", ")
      }
    ),
    sprintf("Seed: %s", format(x$seed)),
    "Stages:",
    stage_lines
  )
}

#' @rdname format.bmmtools_validation_plan
#' @export
print.bmmtools_validation_plan <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @noRd
format_plan_value <- function(v) {
  if (is.null(v)) {
    return("NULL")
  }
  if (is.character(v)) {
    return(paste0("\"", paste(v, collapse = ", "), "\""))
  }
  paste(format(v), collapse = ", ")
}

#' @noRd
format_plan_version <- function(version) {
  if (is.na(version)) "unversioned" else version
}

#' @noRd
na_or <- function(x, y) {
  if (is.null(x) || (length(x) == 1L && is.na(x))) y else x
}

#' Write a validation plan to a file and read it back
#'
#' A plan is a list, so it round-trips through [saveRDS()]. Writing it
#' next to the card it produced keeps a record of what was planned that
#' can be compared with a later version by its `hash`.
#'
#' @param plan A `bmmtools_validation_plan`.
#' @param file A path ending in `.rds`.
#' @return `write_plan()` returns `file` invisibly; `read_plan()` the plan.
#' @examples
#' plan <- validation_plan(
#'   structure(
#'     list(links = list(kappa = "log", thetat = "logit")),
#'     class = c("bmmodel", "circular", "mixture2p")
#'   ),
#'   seed = 1
#' )
#' file <- tempfile(fileext = ".rds")
#' write_plan(plan, file)
#' identical(read_plan(file)$hash, plan$hash)
#' @export
write_plan <- function(plan, file) {
  check_validation_plan(plan)
  check_plan_file(file)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write_atomic(file, function(tmp) saveRDS(plan, tmp))
  invisible(file)
}

#' @rdname write_plan
#' @export
read_plan <- function(file) {
  check_plan_file(file)
  if (!file.exists(file)) {
    cli::cli_abort("{.file {file}} does not exist.")
  }
  plan <- readRDS(file)
  if (!inherits(plan, "bmmtools_validation_plan")) {
    cli::cli_abort(
      "{.file {file}} holds {.obj_type_friendly {plan}}, not a validation \\
       plan."
    )
  }
  plan
}

#' @noRd
check_plan_file <- function(file, call = rlang::caller_env()) {
  if (!is.character(file) || length(file) != 1L || is.na(file)) {
    cli::cli_abort(
      "{.arg file} must be a single path, not {.obj_type_friendly {file}}.",
      call = call
    )
  }
  invisible(file)
}

#' @noRd
check_validation_plan <- function(plan, arg = "plan",
                                  call = rlang::caller_env()) {
  if (!inherits(plan, "bmmtools_validation_plan")) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a validation plan, \\
         not {.obj_type_friendly {plan}}.",
        i = "{.fn validation_plan} makes one."
      ),
      call = call
    )
  }
  invisible(plan)
}
