# Components: several models, one set of simulated people (spec 5, section
# 5.5; local/ARCHITECTURE.md decisions 27 and 33).
#
# A component is one model with its own trials, tasks, generator and
# formula. The subject values of every component, and the covariates, are
# one multivariate normal draw, so parameters of different models can be
# correlated across people. Each component is then simulated and fitted on
# its own. Terms are prefixed with the component name (`m3_c_task1`), which
# is why a name may hold no `_`: the name is everything before the first
# `_` of a prefixed term.

#' Is this a valid component name?
#'
#' The rule is bmm's sanitising of response names in its multivariate
#' branches (decision 33), so a component name survives a later joint fit
#' unchanged.
#'
#' @noRd
valid_component_name <- function(name) {
  is.character(name) && length(name) == 1L && !is.na(name) && nzchar(name) &&
    identical(gsub("\\.|_", "", make.names(name)), name)
}

#' @noRd
check_component_name <- function(name, call = rlang::caller_env()) {
  if (valid_component_name(name)) {
    return(name)
  }
  suggestion <- NULL
  if (is.character(name) && length(name) == 1L && !is.na(name) &&
        nzchar(name)) {
    suggestion <- gsub("\\.|_", "", make.names(name))
  }
  cli::cli_abort(
    c(
      "{.arg name} must be a single syntactic name without {.code _} or \\
       {.code .}.",
      i = "It prefixes the component's terms, as in {.val m3_kappa}.",
      i = if (!is.null(suggestion)) "Use {.val {suggestion}}, for example."
    ),
    call = call
  )
}

#' Validate `pars` and `sds` of a component when they are numbers
#'
#' A function is checked for having no arguments only; its value is checked
#' when `simulate_components()` evaluates it.
#'
#' @noRd
check_component_values <- function(pars, sds, model, tasks, task_col,
                                   call = rlang::caller_env()) {
  for (arg in c("pars", "sds")) {
    value <- if (identical(arg, "pars")) pars else sds
    if (is.function(value) && length(formals(value)) > 0L) {
      cli::cli_abort(
        "{.arg {arg}} is a function, so it must take no arguments.",
        call = call
      )
    }
  }
  if (!is.function(pars)) {
    pars <- expand_task_values(pars, model, tasks, task_col, "pars",
      call = call
    )
    pars <- check_pars(pars, model, tasks, task_col, call = call)
  }
  if (!is.null(sds) && !is.function(sds)) {
    # with a function for `pars`, the SDs are checked against every term
    # the model could have
    reference <- pars
    if (is.function(pars)) {
      info <- model_parameters(model)
      terms <- c(task_terms(info$free, tasks, task_col), info$fixed)
      reference <- stats::setNames(rep(0, length(terms)), terms)
    }
    sds <- expand_task_values(sds, model, tasks, task_col, "sds", call = call)
    check_sds(sds, reference, model, call = call)
  }
  invisible(NULL)
}

#' Describe one model of a multi-model simulation
#'
#' A component is one model with its own population values, trials, tasks,
#' generator and formula. [simulate_components()] simulates several of them
#' for the same people, so that parameters of different models can be
#' correlated across subjects, and [fit_components()] fits each on its own.
#'
#' @param model A `bmmodel`, built with the column names the fit will use.
#' @param pars,n_trials,sds,tasks,task_col,generator As in
#'   [simulate_recovery()], for this component. `pars` and `sds` may be
#'   functions with no arguments, evaluated by [simulate_components()] under
#'   its seed.
#' @param formula `NULL`, or the `bmmformula` [fit_components()] fits this
#'   component with. `NULL` uses [recovery_formula()] with the task column
#'   when there are tasks.
#' @param name The component's name, which prefixes its terms: `name = "m3"`
#'   turns `c_task1` into `m3_c_task1`. A single syntactic name without `_`
#'   or `.`, unique within a set.
#'
#' @return A list of class `bmmtools_component` with the arguments as
#'   validated; `tasks` and `task_col` are `NULL` without tasks.
#'
#' @details
#' What can be checked without drawing is checked here: the model, the
#' number of trials, the tasks, that a generator exists, the formula, and
#' numeric `pars` and `sds` against the model. Correlations and covariates
#' belong to the set and are given to [simulate_components()].
#'
#' @examples
#' \dontrun{
#' recovery_component(
#'   bmm::sdm(resp_error = "y"),
#'   pars = c(c = log(3), kappa = log(5)),
#'   n_trials = 60, sds = c(c = 0.4), name = "sdm"
#' )
#' }
#'
#' @export
recovery_component <- function(model,
                               pars,
                               n_trials,
                               sds = NULL,
                               tasks = NULL,
                               task_col = "task",
                               generator = NULL,
                               formula = NULL,
                               name) {
  check_model(model)
  if (missing(name)) {
    cli::cli_abort(c(
      "{.arg name} is missing.",
      i = "A component needs a name to prefix its terms."
    ))
  }
  name <- check_component_name(name)
  n_trials <- check_count(n_trials, "n_trials")
  design <- check_tasks(tasks, task_col, model)

  if (is.null(generator)) {
    if (is.null(generator_for(model))) {
      cli::cli_abort(c(
        "No generator is known for a model of class {.cls {class(model)}}.",
        i = "Supply {.arg generator}, a function \\
             {.code (pars, n_trials, model)} returning one subject's rows."
      ))
    }
  } else if (!is.function(generator)) {
    cli::cli_abort(
      "{.arg generator} must be a function, \\
       not {.obj_type_friendly {generator}}."
    )
  }
  if (!is.null(formula) && !inherits(formula, "bmmformula")) {
    cli::cli_abort(
      "{.arg formula} must be {.code NULL} or a {.cls bmmformula}, \\
       not {.obj_type_friendly {formula}}."
    )
  }
  check_component_values(pars, sds, model, design$tasks, design$task_col)

  structure(
    list(
      name = name,
      model = model,
      pars = pars,
      sds = sds,
      n_trials = n_trials,
      tasks = design$tasks,
      task_col = design$task_col,
      generator = generator,
      formula = formula
    ),
    class = "bmmtools_component"
  )
}

#' The short label of a model: its most specific class
#'
#' Not `model$name`, which is a full sentence ("Two-parameter mixture model
#' by Zhang and Luck (2008).") and would not fit several components on one
#' line.
#'
#' @noRd
model_label <- function(model) {
  cls <- class(model)
  cls[[length(cls)]]
}

#' @export
print.bmmtools_component <- function(x, ...) {
  extra <- character()
  if (length(x$tasks) > 0L) {
    extra <- paste0("; tasks: ", paste(x$tasks, collapse = ", "))
  }
  cat(
    "<bmmtools_component> ", x$name, ": ", model_label(x$model), ", ",
    x$n_trials, " trials", extra,
    "; generator: ", if (is.null(x$generator)) "adapter" else "user",
    "; formula: ", if (is.null(x$formula)) "default" else "given",
    "\n",
    sep = ""
  )
  invisible(x)
}

# simulate_components() -----------------------------------------------------

#' Validate the list of components and name it by the components
#' @noRd
check_components <- function(components, call = rlang::caller_env()) {
  is_list <- is.list(components) && !is.object(components) &&
    length(components) > 0L
  if (!is_list) {
    cli::cli_abort(
      c(
        "{.arg components} must be a non-empty list of \\
         {.cls bmmtools_component}s, not {.obj_type_friendly {components}}.",
        i = "Wrap a single component in {.code list()}."
      ),
      call = call
    )
  }
  bad <- which(!vapply(
    components, inherits, logical(1),
    what = "bmmtools_component"
  ))
  if (length(bad) > 0L) {
    cli::cli_abort(
      "{.arg components} element{?s} {bad} {?is/are} not a \\
       {.cls bmmtools_component}; build one with {.fn recovery_component}.",
      call = call
    )
  }
  nms <- vapply(components, `[[`, character(1), "name")
  duplicated_names <- unique(nms[duplicated(nms)])
  if (length(duplicated_names) > 0L) {
    cli::cli_abort(
      "Component name{?s} {.val {duplicated_names}} {?is/are} used more \\
       than once.",
      call = call
    )
  }
  given <- names(components)
  if (!is.null(given) && !identical(unname(given), unname(nms))) {
    cli::cli_abort(
      c(
        "The names of {.arg components} must be the component names.",
        i = "List names: {.val {given}}; component names: {.val {nms}}."
      ),
      call = call
    )
  }
  stats::setNames(components, nms)
}

#' Prefix names with a component name
#' @noRd
prefix_terms <- function(x, name) {
  if (length(x) == 0L) {
    return(character())
  }
  paste0(name, "_", x)
}

#' Refuse a bare parameter of a tasked component in the set's `cors`
#' @noRd
check_component_cors <- function(cors, specs, call = rlang::caller_env()) {
  if (!is.matrix(cors)) {
    return(invisible(cors))
  }
  given <- unique(c(rownames(cors), colnames(cors)))
  for (spec in specs) {
    if (is.null(spec$tasks)) next
    free <- model_parameters(spec$model)$free
    bare <- intersect(given, prefix_terms(free, spec$name))
    if (length(bare) > 0L) {
      par <- sub(paste0("^", spec$name, "_"), "", bare[[1L]])
      # nolint next: object_usage_linter. Used by cli's glue interpolation.
      full <- prefix_terms(
        task_terms(par, spec$tasks, spec$task_col),
        spec$name
      )
      cli::cli_abort(
        c(
          "{.arg cors} names {.val {bare}}, but component {.val {spec$name}} \\
           has tasks, so a correlation is between its task terms.",
          i = "Use the full terms, such as {.val {full}}."
        ),
        call = call
      )
    }
  }
  invisible(cors)
}

#' Run code for one component, naming the component in any error
#' @noRd
in_component <- function(name, expr, call) {
  withCallingHandlers(
    expr,
    error = function(e) {
      cli::cli_abort(
        "In component {.val {name}}:",
        parent = e,
        call = call
      )
    }
  )
}

#' Simulate several models for the same simulated people
#'
#' Draws the subject values of every component, and any covariates, from
#' one multivariate normal distribution, then simulates each component's
#' data with its own model and generator. Parameters of different models
#' can therefore be correlated across subjects, which is what a study of
#' individual differences across tasks or models needs.
#'
#' @param components A list of [recovery_component()]s. If the list has
#'   names, they must be the component names.
#' @param n_subjects The number of subjects, shared by every component.
#' @param cors A correlation matrix over prefixed terms (`<name>_<term>`,
#'   such as `m3_c_task1`) and covariates; terms it does not name are
#'   uncorrelated. `NULL` means none are correlated. May be a function with
#'   no arguments, as in [simulate_recovery()].
#' @param covariates Observed, error-free person variables drawn with the
#'   parameters, as in [simulate_recovery()]. A covariate name must differ
#'   from every component name and from every model's parameters and
#'   columns.
#' @param seed A seed applied with `withr::with_seed()` around everything
#'   that is drawn; `NULL` leaves the random number generator alone.
#'
#' @return A list of class `bmmtools_simulation_set` with
#'   * `components`: a named list of `bmmtools_simulation`s, one per
#'     component, with its own unprefixed terms, no covariates and `seed`
#'     `NA`;
#'   * `truth`: the tables of [simulate_recovery()] (`population`,
#'     `subjects`, `sd`, `cor`, `covariates`) over every component, with
#'     prefixed terms; `cor` holds every pair of varying terms and
#'     covariates;
#'   * `pars`, `sds`: the realised values, prefixed and expanded over tasks;
#'   * `cors`: the realised correlation matrix over the varying terms, then
#'     the covariates (`NULL` with fewer than two);
#'   * `covariates` and `covariate_data`, a tibble with `id` and one column
#'     per covariate;
#'   * `links`: every model's link table under prefixed names;
#'   * `specs`: the components as given;
#'   * `n_subjects` and `seed` (`NA` when none).
#'
#' @details
#' **Order of the draws.** Inside the seed, each component's `pars` and then
#' `sds` are evaluated, in component order, then `cors`. The subject values
#' of all components are one draw, `Z %*% chol(cors)`, over the varying
#' terms, component by component. The generators then run in component
#' order, and the covariates are drawn last, so adding a covariate never
#' changes the responses. A set with one component and no covariates gives
#' exactly the simulation [simulate_recovery()] gives with the same
#' arguments and seed. Adding a component after another leaves the first
#' one's subject values unchanged, but not its responses, because the draw
#' for the added terms comes before the generators.
#'
#' **Covariates** are not added to the components' data; they are in
#' `covariate_data` and in `truth$covariates`.
#'
#' **What separate fits can recover.** Each component is fitted on its own
#' by [fit_components()], so no fit knows about the correlation between its
#' parameters and another model's. Correlating the estimates across fits
#' attenuates the correlation by the reliabilities `rel_a` and `rel_b` of
#' the two subject estimates: the `point` estimator of
#' [extract_correlations()] is roughly `rho * sqrt(rel_a * rel_b)`, and the
#' `draws` estimator roughly `rho * rel_a * rel_b`, so neither recovers
#' `rho`. A structural equation model of true against estimated values, from
#' [subject_table()], does; so would a joint multivariate fit of all
#' components, which bmmtools does not support yet.
#'
#' @examples
#' \dontrun{
#' components <- list(
#'   recovery_component(
#'     bmm::mixture2p(resp_error = "y"),
#'     pars = c(kappa = log(8), thetat = qlogis(0.75)),
#'     n_trials = 60, sds = c(kappa = 0.3, thetat = 0.5), name = "m2p"
#'   ),
#'   recovery_component(
#'     bmm::sdm(resp_error = "y"),
#'     pars = c(c = log(3), kappa = log(5)),
#'     n_trials = 60, sds = c(c = 0.4), name = "sdm"
#'   )
#' )
#' cors <- diag(2)
#' dimnames(cors) <- rep(list(c("m2p_thetat", "sdm_c")), 2)
#' cors[1, 2] <- cors[2, 1] <- 0.6
#' set <- simulate_components(components, n_subjects = 100, cors = cors,
#'   seed = 1
#' )
#' fits <- fit_components(set, dir = "fits")
#' recover_correlations(fits, set, estimator = c("draws", "point"))
#' }
#'
#' @export
simulate_components <- function(components,
                                n_subjects,
                                cors = NULL,
                                covariates = NULL,
                                seed = NULL) {
  error_call <- rlang::current_env()
  specs <- check_components(components, call = error_call)
  n_subjects <- check_count(n_subjects, "n_subjects")
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L)) {
    cli::cli_abort("{.arg seed} must be a single number or {.code NULL}.")
  }
  checked <- NULL
  for (spec in specs) {
    checked <- check_covariates(covariates, spec$model, call = error_call)
  }
  covariates <- checked
  clash <- intersect(names(covariates), names(specs))
  if (length(clash) > 0L) {
    cli::cli_abort(
      "{.arg covariates} name{?s} {.val {clash}} {?is/are} also {?a/} \\
       component name{?s}.",
      call = error_call
    )
  }
  for (spec in specs) {
    check_tasks(spec$tasks, spec$task_col, spec$model, covariates,
      call = error_call
    )
  }

  out <- with_seed_if(seed, {
    realised <- lapply(specs, function(spec) {
      in_component(spec$name, {
        pars <- resolve_truth_arg(spec$pars, "pars")
        pars <- expand_task_values(
          pars, spec$model, spec$tasks, spec$task_col, "pars"
        )
        pars <- check_pars(pars, spec$model, spec$tasks, spec$task_col)
        sds <- resolve_truth_arg(spec$sds, "sds")
        sds <- expand_task_values(
          sds, spec$model, spec$tasks, spec$task_col, "sds"
        )
        sds <- check_sds(sds, pars, spec$model)
        list(pars = pars, sds = sds)
      }, call = error_call)
    })
    all_pars <- unlist(unname(lapply(specs, function(spec) {
      values <- realised[[spec$name]]$pars
      stats::setNames(values, prefix_terms(names(values), spec$name))
    })))
    all_sds <- unlist(unname(lapply(specs, function(spec) {
      values <- realised[[spec$name]]$sds
      stats::setNames(values, prefix_terms(names(values), spec$name))
    })))
    cors <- resolve_truth_arg(cors, "cors", want = "matrix", call = error_call)
    check_component_cors(cors, specs, call = error_call)
    full <- check_cors(cors, all_pars, all_sds, covariates, call = error_call)

    parts <- draw_parameter_values(all_pars, all_sds, full, n_subjects)
    sims <- lapply(specs, function(spec) {
      simulate_component(spec, realised[[spec$name]], parts$values, full,
        n_subjects,
        call = error_call
      )
    })
    # covariates are drawn after every generator, so their random numbers
    # never shift the ones that produced the responses
    values <- draw_covariate_values(parts, full, covariates)
    list(
      sims = sims, values = values, pars = all_pars, sds = all_sds,
      cors = full
    )
  })

  ids <- factor(seq_len(n_subjects), levels = seq_len(n_subjects))
  covariate_data <- tibble::tibble(id = ids)
  for (g in names(covariates)) {
    covariate_data[[g]] <- unname(out$values[, g])
  }
  links <- component_links(specs)

  structure(
    list(
      components = out$sims,
      truth = truth_tables(
        out$pars, out$sds, out$values, out$cors, covariates
      ),
      pars = out$pars,
      sds = out$sds,
      cors = if (nrow(out$cors) < 2L) NULL else out$cors,
      covariates = covariates,
      covariate_data = covariate_data,
      links = links,
      specs = specs,
      n_subjects = n_subjects,
      seed = if (is.null(seed)) NA_real_ else as.double(seed)
    ),
    class = "bmmtools_simulation_set"
  )
}

#' Every component model's link table under prefixed names
#'
#' Needs no draw, so `recovery_grid()` knows a set's links before a cell
#' runs.
#'
#' @noRd
component_links <- function(specs) {
  unlist(unname(lapply(specs, function(spec) {
    table <- spec$model$links
    stats::setNames(
      as.character(unlist(table, use.names = FALSE)),
      prefix_terms(names(table), spec$name)
    )
  })))
}

#' Simulate one component from its slice of the joint draw
#'
#' The slice goes in as `subject_pars`, so `simulate_recovery()` draws
#' nothing and only runs the generator; its principal sub-matrix of the
#' correlations labels the component's own `truth$cor`.
#'
#' @noRd
simulate_component <- function(spec, realised, values, full, n_subjects,
                               call) {
  varying <- names(realised$sds)[realised$sds > 0]
  prefixed <- prefix_terms(varying, spec$name)
  ids <- rownames(values)
  subject_pars <- tibble::tibble(
    id = rep(ids, times = length(varying)),
    term = rep(varying, each = length(ids)),
    true_value = as.double(values[, prefixed, drop = TRUE])
  )
  own_cors <- NULL
  if (length(varying) >= 2L) {
    own_cors <- full[prefixed, prefixed, drop = FALSE]
    dimnames(own_cors) <- list(varying, varying)
  }
  in_component(spec$name, {
    simulate_recovery(
      spec$model, realised$pars,
      n_subjects = n_subjects, n_trials = spec$n_trials,
      sds = realised$sds, cors = own_cors,
      tasks = spec$tasks, task_col = spec$task_col %||% "task",
      subject_pars = subject_pars,
      generator = spec$generator,
      seed = NULL
    )
  }, call = call)
}

#' @export
print.bmmtools_simulation_set <- function(x, ...) {
  described <- vapply(names(x$components), function(nm) {
    sim <- x$components[[nm]]
    paste0(nm, " (", model_label(sim$model), ", ", sim$n_trials, " trials)")
  }, character(1))
  varying <- names(x$sds)[x$sds > 0]
  extra <- character()
  if (length(x$covariates) > 0L) {
    extra <- c(extra, paste0(
      "; covariates: ", paste(names(x$covariates), collapse = ", ")
    ))
  }
  nonzero <- x$truth$cor[x$truth$cor$true_value != 0, , drop = FALSE]
  if (nrow(nonzero) > 0L) {
    extra <- c(extra, paste0(
      "; nonzero correlations: ",
      paste0(
        nonzero$term, " = ", signif(nonzero$true_value, 3),
        collapse = ", "
      )
    ))
  }
  cat(
    "<bmmtools_simulation_set> ", length(described), " component",
    if (length(described) == 1L) "" else "s", ": ",
    paste(described, collapse = ", "), "; ",
    x$n_subjects, " subjects; varying terms: ", length(varying),
    extra,
    "\n",
    sep = ""
  )
  invisible(x)
}

# fit_components() ----------------------------------------------------------

#' Validate `prior` of fit_components(): `NULL` or a list by component
#' @noRd
check_component_prior <- function(prior, components,
                                  call = rlang::caller_env()) {
  if (is.null(prior)) {
    return(NULL)
  }
  nms <- names(prior)
  is_named_list <- is.list(prior) && !is.object(prior) && !is.null(nms) &&
    !anyNA(nms) && all(nzchar(nms))
  if (!is_named_list) {
    cli::cli_abort(
      c(
        "{.arg prior} must be {.code NULL} or a named list with one prior \\
         per component, not {.obj_type_friendly {prior}}.",
        i = "For a prior on component {.val {components[[1L]]}}, use \\
             {.code list({components[[1L]]} = <prior>)}."
      ),
      call = call
    )
  }
  unknown <- setdiff(nms, components)
  if (length(unknown) > 0L) {
    cli::cli_abort(
      c(
        "{.arg prior} names {.val {unknown}}, which {?is/are} not a \\
         component.",
        i = "Components: {.val {components}}."
      ),
      call = call
    )
  }
  prior
}

#' The formula a component is fitted with: its own, else the default
#' @noRd
component_formula <- function(spec) {
  spec$formula %||% recovery_formula(spec$model, task_col = spec$task_col)
}

#' Fit every component of a simulation set, each on its own
#'
#' Fits each component's data with its own model through [fit_cached()], so
#' a second call reuses the saved fits while nothing that determines them
#' has changed.
#'
#' @param sim_set A `bmmtools_simulation_set` from [simulate_components()].
#' @param dir The directory the fits are saved in, one file per component,
#'   named after it (`<dir>/<name>.rds`, with its cache key next to it).
#' @param prior `NULL`, or a named list with a prior for some of the
#'   components, such as `list(sdm = <brmsprior>)`. A component not named is
#'   fitted with bmm's default priors.
#' @param ... Passed to [fit_cached()] and on to the fitter, for every
#'   component: `refit`, and sampler arguments such as `chains`, `iter` or
#'   `backend`. When the set was simulated with a seed and `...` has no
#'   `seed`, the set's seed is passed as `seed`, the same for every component.
#' @param .fitter The fitting function, `bmm::bmm()` by default; see
#'   [fit_cached()].
#'
#' @return A named list of fits, one per component, of class
#'   `bmmtools_fit_set`, which [extract_correlations()],
#'   [recover_correlations()] and [subject_table()] read as one set.
#'
#' @details
#' A component fitted without a `formula` uses [recovery_formula()], with
#' the task column when the component has tasks.
#'
#' @examples
#' \dontrun{
#' fits <- fit_components(set, dir = "fits", chains = 4, iter = 2000)
#' fits
#' }
#'
#' @export
fit_components <- function(sim_set, dir, prior = NULL, ..., .fitter = NULL) {
  error_call <- rlang::current_env()
  if (!inherits(sim_set, "bmmtools_simulation_set")) {
    cli::cli_abort(
      "{.arg sim_set} must be a {.cls bmmtools_simulation_set} from \\
       {.fn simulate_components}, not {.obj_type_friendly {sim_set}}."
    )
  }
  bad_dir <- !is.character(dir) || length(dir) != 1L || is.na(dir) ||
    !nzchar(dir)
  if (bad_dir) {
    cli::cli_abort(
      "{.arg dir} must be a single path, not {.obj_type_friendly {dir}}."
    )
  }
  components <- names(sim_set$components)
  prior <- check_component_prior(prior, components, call = error_call)
  dots <- rlang::list2(...)
  if ("file" %in% names(dots)) {
    cli::cli_abort(
      c(
        "{.arg file} cannot be passed to {.fn fit_components}.",
        i = "Each fit is saved as {.file <dir>/<component>.rds}; choose \\
             {.arg dir} instead."
      )
    )
  }
  # as run_cell() does for a grid cell: the simulation's seed, unless the
  # caller gave one
  if (!is.na(sim_set$seed) && !"seed" %in% names(dots)) {
    dots <- c(list(seed = sim_set$seed), dots)
  }

  fits <- lapply(stats::setNames(nm = components), function(nm) {
    sim <- sim_set$components[[nm]]
    formula <- component_formula(sim_set$specs[[nm]])
    rlang::exec(
      fit_cached,
      formula = formula, data = sim$data, model = sim$model,
      file = file.path(dir, nm), prior = prior[[nm]], !!!dots,
      .fitter = .fitter
    )
  })
  structure(fits, class = c("bmmtools_fit_set", "list"))
}

#' @export
print.bmmtools_fit_set <- function(x, ...) {
  cat(
    "<bmmtools_fit_set> ", length(x), " fit",
    if (length(x) == 1L) "" else "s", "\n",
    sep = ""
  )
  for (nm in names(x)) {
    cache <- attr(x[[nm]], "bmmtools_cache")
    how <- ""
    if (!is.null(cache)) {
      how <- paste0(
        if (isTRUE(cache$reused)) ", reused " else ", fitted ",
        cache$file
      )
    }
    cat("  ", nm, ": <", class(x[[nm]])[[1L]], ">", how, "\n", sep = "")
  }
  invisible(x)
}
