# The generate layer (spec 3, sections 1 and 2; local/ARCHITECTURE.md
# decisions 3, 5, 13, 17).
#
# A bmm model plus population values on the link scale become data and
# a truth table. Subject values are drawn on the link scale, converted
# through inverse_link() and handed, one subject at a time, to the
# model's own generator or to one the user supplies.

#' The parameters a model estimates and the ones it fixes
#'
#' Derived from the model object's slots, which is what
#' `bmm::parameters()` reads too (measured 2026-09-08 on six models: the
#' free set is the names of `model$links` minus the names of
#' `model$fixed_parameters`). Reading the slots directly keeps this
#' layer usable on a model object built without bmm loaded.
#'
#' @noRd
model_parameters <- function(model) {
  links <- model$links
  fixed <- model$fixed_parameters
  list(
    free = setdiff(names(links), names(fixed)),
    fixed = names(fixed),
    links = links,
    fixed_values = fixed
  )
}

#' @noRd
check_model <- function(model, call = rlang::caller_env()) {
  if (!inherits(model, "bmmodel")) {
    cli::cli_abort(
      "{.arg model} must be a {.cls bmmodel}, \\
       not {.obj_type_friendly {model}}.",
      call = call
    )
  }
  invisible(model)
}

#' @noRd
check_count <- function(x, name, call = rlang::caller_env()) {
  bad <- !is.numeric(x) || length(x) != 1L || is.na(x) || x < 1 ||
    x != round(x)
  if (bad) {
    cli::cli_abort(
      "{.arg {name}} must be a single whole number of at least 1, \\
       not {.obj_type_friendly {x}}.",
      call = call
    )
  }
  as.integer(x)
}

#' Validate the population values against the model
#'
#' With tasks, `pars` has already been expanded to full terms by
#' `expand_task_values()`, and every free parameter needs a value per task.
#'
#' @noRd
check_pars <- function(pars, model, tasks = NULL, task_col = NULL,
                       call = rlang::caller_env()) {
  info <- model_parameters(model)
  info$free <- task_terms(info$free, tasks, task_col)
  bad <- !is.numeric(pars) || is.null(names(pars)) ||
    anyNA(names(pars)) || !all(nzchar(names(pars)))
  if (bad) {
    cli::cli_abort(
      "{.arg pars} must be a named numeric vector of population values \\
       on the link scale.",
      call = call
    )
  }
  known <- c(info$free, info$fixed)
  unknown <- setdiff(names(pars), known)
  if (length(unknown) > 0L) {
    cli::cli_abort(
      c(
        "{.arg pars} names {.val {unknown}}, which the model does not have.",
        i = "Its parameters are {.val {known}}."
      ),
      call = call
    )
  }
  missing <- setdiff(info$free, names(pars))
  if (length(missing) > 0L) {
    cli::cli_abort(
      c(
        "{.arg pars} is missing {.val {missing}}.",
        i = "Every estimated parameter needs a value: {.val {info$free}}."
      ),
      call = call
    )
  }
  pars
}

#' Validate the between-subject SDs and fill the ones not given with 0
#' @noRd
check_sds <- function(sds, pars, model, call = rlang::caller_env()) {
  info <- model_parameters(model)
  out <- stats::setNames(rep(0, length(pars)), names(pars))
  if (is.null(sds)) {
    return(out)
  }
  if (!is.numeric(sds) || is.null(names(sds)) || anyNA(names(sds))) {
    cli::cli_abort(
      "{.arg sds} must be a named numeric vector of between-subject SDs \\
       on the link scale, or {.code NULL}.",
      call = call
    )
  }
  fixed <- intersect(names(sds), info$fixed)
  if (length(fixed) > 0L) {
    cli::cli_abort(
      "{.arg sds} names {.val {fixed}}, which the model fixes; a fixed \\
       parameter cannot vary between subjects.",
      call = call
    )
  }
  unknown <- setdiff(names(sds), names(pars))
  if (length(unknown) > 0L) {
    cli::cli_abort(
      "{.arg sds} names {.val {unknown}}, which {.arg pars} does not give.",
      call = call
    )
  }
  negative <- names(sds)[is.na(sds) | sds < 0]
  if (length(negative) > 0L) {
    cli::cli_abort(
      "{.arg sds} must be non-negative; {.val {negative}} {?is/are} not.",
      call = call
    )
  }
  out[names(sds)] <- sds
  out
}

#' Validate subject values supplied instead of drawn
#'
#' Covariates, when there are any, must be supplied too: drawing them
#' afresh beside given subject values would break the correlation between
#' the two.
#'
#' @return A matrix, subjects by parameters then covariates, on the link
#'   scale.
#' @noRd
check_subject_pars <- function(subject_pars, pars, sds, n_subjects,
                               covariates = NULL,
                               call = rlang::caller_env()) {
  needed <- c("id", "term", "true_value")
  bad <- !is.data.frame(subject_pars) ||
    !all(needed %in% names(subject_pars))
  if (bad) {
    cli::cli_abort(
      "{.arg subject_pars} must be a data frame with the columns \\
       {.val {needed}}.",
      call = call
    )
  }
  ids <- as.character(seq_len(n_subjects))
  varying <- c(names(sds)[sds > 0], names(covariates))
  bad_id <- setdiff(unique(as.character(subject_pars$id)), ids)
  # a term outside `varying` would drive the data and be absent from the
  # subject truth, so it is refused rather than silently dropped
  bad_term <- setdiff(unique(subject_pars$term), varying)
  if (length(bad_id) > 0L || length(bad_term) > 0L) {
    cli::cli_abort(
      c(
        "{.arg subject_pars} does not match the design.",
        i = if (length(bad_id) > 0L) "Unknown id{?s}: {.val {bad_id}}.",
        i = if (length(bad_term) > 0L) {
          "Term{?s} {.val {bad_term}} {?does/do} not vary; only parameters \\
           with a positive {.arg sds} entry can take subject values."
        }
      ),
      call = call
    )
  }
  values <- matrix(
    c(pars, rep(NA_real_, length(covariates))),
    nrow = n_subjects, ncol = length(pars) + length(covariates), byrow = TRUE,
    dimnames = list(ids, c(names(pars), names(covariates)))
  )
  for (term in unique(subject_pars$term)) {
    rows <- subject_pars[subject_pars$term == term, , drop = FALSE]
    have <- as.character(rows$id)
    if (!setequal(have, ids) || anyDuplicated(have) > 0L) {
      cli::cli_abort(
        "{.arg subject_pars} must give {.val {term}} once for every \\
         subject 1 to {n_subjects}.",
        call = call
      )
    }
    values[have, term] <- rows$true_value
  }
  missing <- setdiff(varying, unique(subject_pars$term))
  if (length(missing) > 0L) {
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    what <- if (all(missing %in% names(covariates))) {
      "covariate"
    } else {
      "varying term"
    }
    cli::cli_abort(
      "{.arg subject_pars} must give the {what}{?s} {.val {missing}}.",
      call = call
    )
  }
  values
}

#' One subject's parameters on the natural scale, fixed ones included
#' @noRd
natural_pars <- function(link_values, model) {
  info <- model_parameters(model)
  all_values <- as.list(link_values)
  for (p in setdiff(info$fixed, names(all_values))) {
    all_values[[p]] <- info$fixed_values[[p]]
  }
  lapply(stats::setNames(names(all_values), names(all_values)), function(p) {
    inverse_link(all_values[[p]], info$links[[p]] %||% "identity")
  })
}

#' Check that the generator produced the model's response columns
#'
#' Only `resp_vars` are checked: `other_vars` mixes column names with
#' constants (`sdt_yn` carries `dist = "normal"`, measured), so it
#' cannot be read as a list of columns. bmm's own data checks, reached
#' through the mock backend in the tests, cover the rest.
#'
#' @noRd
check_generated <- function(data, model, covariate_names = NULL,
                            task_col = NULL, call = rlang::caller_env()) {
  if (!is.data.frame(data)) {
    cli::cli_abort(
      "The generator must return a data frame, \\
       not {.obj_type_friendly {data}}.",
      call = call
    )
  }
  clash <- intersect(covariate_names, names(data))
  if (length(clash) > 0L) {
    cli::cli_abort(
      "The generator returned column{?s} {.val {clash}}, which {?is/are} \\
       also {?a/} covariate name{?s}.",
      call = call
    )
  }
  if (!is.null(task_col) && task_col %in% names(data)) {
    cli::cli_abort(
      "The generator returned a column {.val {task_col}}, which is the \\
       task column; choose another {.arg task_col}.",
      call = call
    )
  }
  needed <- unlist(model$resp_vars, use.names = FALSE)
  missing <- setdiff(needed, names(data))
  if (length(missing) > 0L) {
    cli::cli_abort(
      c(
        "The generated data lack the model's column{?s} {.val {missing}}.",
        i = "The generator returned {.val {names(data)}}."
      ),
      call = call
    )
  }
  invisible(data)
}

# tasks (spec 5, section 5.4; local/ARCHITECTURE.md decisions 27 and 30) ----

#' Validate `tasks` and, when there are tasks, `task_col`
#'
#' A task level becomes part of a brms coefficient name (`task1`) and so of
#' a term, which is why it may hold letters and digits only. A single level
#' is refused: brms would estimate one coefficient per parameter, which the
#' extraction reads back under the bare parameter name, so the truth terms
#' (`kappa_task1`) would never match.
#'
#' @return A list with `tasks` and `task_col`, both `NULL` without tasks.
#' @noRd
check_tasks <- function(tasks, task_col, model, covariates = NULL,
                        call = rlang::caller_env()) {
  if (is.null(tasks)) {
    return(list(tasks = NULL, task_col = NULL))
  }
  if (!is.character(tasks) || anyNA(tasks) || length(tasks) < 2L) {
    cli::cli_abort(
      "{.arg tasks} must be {.code NULL} or a character vector of at least \\
       two task levels, not {.obj_type_friendly {tasks}}.",
      call = call
    )
  }
  bad_level <- tasks[!grepl("^[A-Za-z0-9]+$", tasks)]
  if (length(bad_level) > 0L) {
    cli::cli_abort(
      "{.arg tasks} level{?s} {.val {bad_level}} must consist of letters \\
       and digits only.",
      call = call
    )
  }
  if (anyDuplicated(tasks) > 0L) {
    cli::cli_abort(
      "{.arg tasks} names {.val {unique(tasks[duplicated(tasks)])}} more \\
       than once.",
      call = call
    )
  }
  bad_col <- !is.character(task_col) || length(task_col) != 1L ||
    is.na(task_col) || !grepl("^[A-Za-z][A-Za-z0-9]*$", task_col)
  if (bad_col) {
    cli::cli_abort(
      "{.arg task_col} must be a single name of letters and digits that \\
       starts with a letter.",
      call = call
    )
  }
  info <- model_parameters(model)
  # other_vars also holds constants (sdt_yn's dist = "normal"), so a name
  # equal to one is refused too; that is deliberate, as in
  # check_covariates(): it catches a clash with a real column such as
  # n_trials before any generator runs, at the price of a few odd names
  columns <- Filter(is.character, c(model$resp_vars, model$other_vars))
  taken <- c(
    "id", info$free, info$fixed, unlist(columns, use.names = FALSE),
    names(covariates)
  )
  if (task_col %in% taken) {
    cli::cli_abort(
      c(
        "{.arg task_col} {.val {task_col}} cannot be used.",
        i = "It must differ from {.val id}, the model's parameters and data \\
             columns, and the covariates."
      ),
      call = call
    )
  }
  list(tasks = tasks, task_col = task_col)
}

#' The full terms of parameters over tasks, parameter-major
#'
#' `c("kappa", "thetat")` over tasks 1 and 2 of `task` gives `kappa_task1`,
#' `kappa_task2`, `thetat_task1`, `thetat_task2`: brms's coefficient name
#' for `0 + task` after the parameter (D27). Without tasks, the parameters.
#'
#' @noRd
task_terms <- function(pars, tasks = NULL, task_col = NULL) {
  if (is.null(tasks) || length(pars) == 0L) {
    return(pars)
  }
  as.vector(t(outer(pars, paste0("_", task_col, tasks), paste0)))
}

#' Expand a `pars` or `sds` vector to full task terms
#'
#' A name is a free parameter, shared by every task, or a full term, which
#' overrides the shared value for its task. A fixed parameter stays bare
#' (in `sds` it is left for `check_sds()` to refuse). Free parameters come
#' out in the order they first appear, then any not named at all (so that
#' `check_pars()` can report them missing), each over the tasks in order;
#' a term with no value is left out.
#'
#' @return `x` unchanged without tasks or when it is not a named numeric
#'   vector (the checks after it report that); otherwise the expanded
#'   vector.
#' @noRd
expand_task_values <- function(x, model, tasks, task_col, arg,
                               call = rlang::caller_env()) {
  shaped <- is.numeric(x) && !is.null(names(x)) && !anyNA(names(x))
  if (is.null(tasks) || !shaped) {
    return(x)
  }
  info <- model_parameters(model)
  given <- names(x)
  per_task <- lapply(stats::setNames(nm = info$free), function(p) {
    task_terms(p, tasks, task_col)
  })
  owner <- rep(NA_character_, length(given))
  for (p in info$free) {
    owner[given == p | given %in% per_task[[p]]] <- p
  }
  unknown <- given[is.na(owner) & !given %in% info$fixed]
  if (length(unknown) > 0L) {
    example <- task_terms(info$free, tasks, task_col)[1L]
    cli::cli_abort(
      c(
        "{.arg {arg}} names {.val {unknown}}, which {?is/are} neither a \\
         parameter of the model nor one of its task terms.",
        i = paste0(
          "Name a free parameter for every task, or one task of it as ",
          "{.val ", example, "}; fixed parameters are not per task. ",
          "Tasks: {.val {tasks}}."
        )
      ),
      call = call
    )
  }

  out <- numeric()
  for (p in unique(c(owner[!is.na(owner)], info$free))) {
    for (term in per_task[[p]]) {
      if (term %in% given) {
        out[[term]] <- x[[term]]
      } else if (p %in% given) {
        out[[term]] <- x[[p]]
      }
    }
  }
  c(out, x[given %in% info$fixed])
}

#' Refuse bare parameter names in `cors` when there are tasks
#' @noRd
check_task_cors <- function(cors, model, tasks, task_col,
                            call = rlang::caller_env()) {
  if (is.null(tasks) || !is.matrix(cors)) {
    return(invisible(cors))
  }
  free <- model_parameters(model)$free
  bare <- intersect(unique(c(rownames(cors), colnames(cors))), free)
  if (length(bare) > 0L) {
    cli::cli_abort(
      c(
        "{.arg cors} names {.val {bare}}, but with tasks a correlation is \\
         between task terms.",
        i = "Use the full terms, such as \\
             {.val {task_terms(bare[[1L]], tasks, task_col)}}."
      ),
      call = call
    )
  }
  invisible(cors)
}

#' Which columns of the subject values each generator call reads
#'
#' One element per task (a single one without tasks): a character vector
#' of value columns named by the bare parameter the generator sees.
#'
#' @noRd
task_layout <- function(pars, model, tasks, task_col) {
  terms <- names(pars)
  if (is.null(tasks)) {
    return(list(stats::setNames(terms, terms)))
  }
  info <- model_parameters(model)
  fixed <- terms[terms %in% info$fixed]
  lapply(tasks, function(task) {
    suffix <- paste0("_", task_col, task)
    own <- terms[terms %in% paste0(info$free, suffix)]
    bare <- substr(own, 1L, nchar(own) - nchar(suffix))
    stats::setNames(c(own, fixed), c(bare, fixed))
  })
}

#' Run the generator for every subject, and every task within a subject
#' @noRd
generate_data <- function(values, pars, model, generator, n_subjects,
                          n_trials, tasks, task_col, cov_names) {
  layout <- task_layout(pars, model, tasks, task_col)
  pieces <- lapply(seq_len(n_subjects), function(i) {
    lapply(seq_along(layout), function(k) {
      # keep the names when the matrix has a single column
      link_values <- stats::setNames(values[i, layout[[k]]], names(layout[[k]]))
      rows <- generator(natural_pars(link_values, model), n_trials, model)
      check_generated(rows, model, cov_names, task_col)
      rows <- tibble::as_tibble(rows)
      if (!is.null(tasks)) {
        rows <- tibble::add_column(
          rows,
          !!task_col := factor(rep(tasks[[k]], nrow(rows)), levels = tasks),
          .before = 1L
        )
      }
      tibble::add_column(
        rows,
        id = factor(rep(i, nrow(rows)), levels = seq_len(n_subjects)),
        .before = 1L
      )
    })
  })
  dplyr::bind_rows(do.call(c, pieces))
}

#' Run an expression under a seed, or as is when there is none
#' @noRd
with_seed_if <- function(seed, expr) {
  if (is.null(seed)) {
    return(expr)
  }
  withr::with_seed(seed, expr)
}

#' The truth tables of a simulation, on the link scale
#' @noRd
truth_tables <- function(pars, sds, values, cors = NULL, covariates = NULL) {
  varying <- names(sds)[sds > 0]
  cov_names <- names(covariates)
  population <- tibble::tibble(
    term = names(pars),
    true_value = as.double(unname(pars))
  )
  subjects <- tibble::tibble(
    id = rep(rownames(values), times = length(varying)),
    term = rep(varying, each = nrow(values)),
    true_value = as.double(values[, varying, drop = TRUE])
  )
  empty <- tibble::tibble(
    id = character(), term = character(), true_value = double()
  )
  if (length(varying) == 0L) subjects <- empty
  covariate_rows <- empty
  if (length(cov_names) > 0L) {
    covariate_rows <- tibble::tibble(
      id = rep(rownames(values), times = length(cov_names)),
      term = rep(cov_names, each = nrow(values)),
      true_value = as.double(values[, cov_names, drop = TRUE])
    )
  }
  list(
    population = population,
    subjects = subjects,
    sd = sd_table(sds),
    cor = cor_table(cors),
    covariates = covariate_rows
  )
}

#' Simulate data from a bmm model with known parameters
#'
#' The generate half of the engine: a model plus population values on
#' the link scale become a data set and the truth that produced it, so
#' that [recover()] and [recover_subjects()] can score a fit of that
#' data. Subject values are drawn on the link scale around the
#' population values, converted to the natural scale through
#' [inverse_link()], and handed to the model's own `r<model>()` generator
#' (decision 13) or to a function you supply.
#'
#' @param model A `bmmodel`, built with the column names the fit will
#'   use, so the generated columns match.
#' @param pars Named numeric: population values **on the link scale
#'   under bmm's parameter names**, one per estimated parameter. A
#'   fixed parameter may be given too, in which case the generator uses
#'   that value. May also be a function with no arguments returning such
#'   a vector, evaluated under `seed`, for random hyperparameters.
#' @param n_subjects,n_trials Subjects, and trials per subject and per
#'   row of the generator's layout (for `sdt_yn`, per stimulus class).
#' @param sds Named numeric: between-subject standard deviations on the
#'   link scale. A parameter not named does not vary. `NULL` means no
#'   parameter varies. May be a function, as `pars`.
#' @param cors A correlation matrix between subject values on the link
#'   scale, with the parameter and covariate names as dimnames; terms it
#'   does not name are uncorrelated. `NULL` means none are correlated. May
#'   be a function, as `pars`. [cors_from_factors()] builds one from
#'   factor loadings.
#' @param covariates Observed, error-free person variables drawn jointly
#'   with the parameters: a named list of `c(mean = , sd = )`. Each
#'   becomes a subject-constant data column after `id`. A name must be
#'   syntactic, contain no `_`, and differ from the model's parameters and
#'   columns.
#' @param tasks `NULL`, or a character vector of at least two task levels
#'   (letters and digits) for a design in which every subject does every
#'   task. Each free parameter then has one value per task, under the term
#'   `<parameter>_<task_col><level>` (`kappa_task1`), brms's coefficient
#'   name for `0 + task`. See the details.
#' @param task_col The name of the task column in `data`, `"task"` by
#'   default. Used only with `tasks`.
#' @param subject_pars A data frame `id`, `term`, `true_value` (link
#'   scale) of subject values to use instead of drawing them, so that
#'   replications can share the same simulated people. With `covariates`
#'   it must give their values too.
#' @param generator A function `(pars, n_trials, model)` returning one
#'   subject's rows as a data frame with the model's column names;
#'   `pars` is a named list on the natural scale, fixed parameters
#'   included. `NULL` uses the adapter bmmtools ships for the model.
#' @param seed A seed applied with `withr::with_seed()` around the draws
#'   and the generator; `NULL` leaves the random number generator alone.
#'
#' @return A list of class `bmmtools_simulation` with `data` (a tibble,
#'   `id` first, then any covariates, then the task column), `truth` (a
#'   list of tibbles on the link scale: `population` with `term` and
#'   `true_value`; `subjects` with `id`, `term` and `true_value`; `sd` with
#'   `term` and `true_value`; `cor` with `term`, `var1`, `var2` and
#'   `true_value`;
#'   `covariates` with `id`, `term` and `true_value`), the realised
#'   `pars`, `sds`, `cors` (the full matrix over varying parameters then
#'   covariates, `NULL` with fewer than two) and `covariates`,
#'   `n_subjects`, `n_trials`, `seed` (`NA` when none), `model`,
#'   `generator`, `tasks` and `task_col` (both `NULL` without tasks).
#'
#' @details
#' Adapters exist for `sdt_yn`, `sdt_mafc`, `ezdm` (three parameters),
#' `ddm`, `mixture2p` and `sdm`. Every other model takes a `generator`.
#' The truth for the subjects and for the SDs lists only the parameters
#' that vary, because a parameter that does not vary has nothing
#' person-level to recover.
#'
#' **Correlated draws.** Subject values and covariates are one
#' multivariate normal draw, `Z %*% chol(cors)`, over the varying
#' parameters in `pars` order and then the covariates, scaled by the SDs
#' afterwards. Without correlations this gives exactly the values an
#' uncorrelated draw gave in earlier versions, so seeded simulations
#' stay reproducible; adding a covariate never changes the parameter
#' values. A correlation pair is named `<a>__<b>`, the two names sorted
#' in the C locale.
#'
#' **Tasks.** With `tasks`, `pars` and `sds` may name a parameter, which
#' gives every task that value, or a full term such as `kappa_task2`, which
#' overrides it for that task; the realised `pars` and `sds` are stored
#' with full terms. `cors`, `subject_pars` and the truth tables use full
#' terms only, so a correlation between tasks is, for example,
#' `kappa_task1__kappa_task2`. The generator is called once per subject and
#' task with that task's values under the bare parameter names, and
#' `n_trials` is per task. The task values are cell means: fit them with
#' `recovery_formula(model, task_col = "task")`. Fixed parameters are the
#' same in every task. Adding tasks changes the random numbers drawn
#' compared with a simulation without them; `tasks = NULL` gives exactly
#' the simulation of earlier versions.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_recovery(
#'   bmm::mixture2p(resp_error = "y"),
#'   pars = c(kappa = log(8), thetat = qlogis(0.75)),
#'   n_subjects = 30, n_trials = 60,
#'   sds = c(kappa = 0.3, thetat = 0.5), seed = 1
#' )
#' fit <- bmm::bmm(recovery_formula(sim$model), sim$data, sim$model)
#' recover(fit, sim$truth$population)
#' recover_subjects(fit, sim$truth$subjects)
#'
#' # two tasks, kappa lower in the second, correlated .6 across tasks
#' cors <- diag(2)
#' dimnames(cors) <- rep(list(c("kappa_task1", "kappa_task2")), 2)
#' cors[1, 2] <- cors[2, 1] <- 0.6
#' two_tasks <- simulate_recovery(
#'   bmm::mixture2p(resp_error = "y"),
#'   pars = c(kappa = log(8), kappa_task2 = log(5), thetat = qlogis(0.75)),
#'   n_subjects = 30, n_trials = 60,
#'   sds = c(kappa_task1 = 0.3, kappa_task2 = 0.3), cors = cors,
#'   tasks = c("1", "2"), seed = 1
#' )
#' formula <- recovery_formula(
#'   two_tasks$model,
#'   re_cor = "within", task_col = "task"
#' )
#' }
#'
#' @export
simulate_recovery <- function(model,
                              pars,
                              n_subjects,
                              n_trials,
                              sds = NULL,
                              cors = NULL,
                              covariates = NULL,
                              tasks = NULL,
                              task_col = "task",
                              subject_pars = NULL,
                              generator = NULL,
                              seed = NULL) {
  check_model(model)
  n_subjects <- check_count(n_subjects, "n_subjects")
  n_trials <- check_count(n_trials, "n_trials")
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L)) {
    cli::cli_abort("{.arg seed} must be a single number or {.code NULL}.")
  }

  generator_name <- "user"
  if (is.null(generator)) {
    generator <- generator_for(model)
    if (is.null(generator)) {
      cli::cli_abort(c(
        "No generator is known for a model of class \\
         {.cls {class(model)}}.",
        i = "Supply {.arg generator}, a function \\
             {.code (pars, n_trials, model)} returning one subject's rows."
      ))
    }
    generator_name <- paste0("adapter:", adapter_name(model))
  }
  if (!is.function(generator)) {
    cli::cli_abort(
      "{.arg generator} must be a function, \\
       not {.obj_type_friendly {generator}}."
    )
  }

  covariates <- check_covariates(covariates, model)
  cov_names <- names(covariates)
  design <- check_tasks(tasks, task_col, model, covariates)
  tasks <- design$tasks
  task_col <- design$task_col

  # functions first, under the seed and in a fixed order, so that random
  # hyperparameters continue the stream the subject draws then take
  out <- with_seed_if(seed, {
    pars <- resolve_truth_arg(pars, "pars")
    pars <- expand_task_values(pars, model, tasks, task_col, "pars")
    pars <- check_pars(pars, model, tasks, task_col)
    sds <- resolve_truth_arg(sds, "sds")
    sds <- expand_task_values(sds, model, tasks, task_col, "sds")
    sds <- check_sds(sds, pars, model)
    cors <- resolve_truth_arg(cors, "cors", want = "matrix")
    check_task_cors(cors, model, tasks, task_col)
    cors <- check_cors(cors, pars, sds, covariates)
    parts <- NULL
    values <- if (is.null(subject_pars)) {
      parts <- draw_parameter_values(pars, sds, cors, n_subjects)
      parts$values
    } else {
      check_subject_pars(subject_pars, pars, sds, n_subjects, covariates)
    }
    data <- generate_data(
      values, pars, model, generator, n_subjects, n_trials,
      tasks, task_col, cov_names
    )
    # covariates are drawn after the generator, so their random numbers
    # never shift the ones that produced the responses
    if (!is.null(parts)) {
      values <- draw_covariate_values(parts, cors, covariates)
    }
    for (g in rev(cov_names)) {
      data <- tibble::add_column(
        data,
        !!g := unname(values[as.integer(data$id), g]),
        .after = "id"
      )
    }
    list(
      data = data, values = values,
      pars = pars, sds = sds, cors = cors
    )
  })

  structure(
    list(
      data = out$data,
      truth = truth_tables(
        out$pars, out$sds, out$values, out$cors, covariates
      ),
      pars = out$pars,
      sds = out$sds,
      cors = if (nrow(out$cors) < 2L) NULL else out$cors,
      covariates = covariates,
      n_subjects = n_subjects,
      n_trials = n_trials,
      seed = if (is.null(seed)) NA_real_ else as.double(seed),
      model = model,
      generator = generator_name,
      tasks = tasks,
      task_col = task_col
    ),
    class = "bmmtools_simulation"
  )
}

#' @export
print.bmmtools_simulation <- function(x, ...) {
  varying <- names(x$sds)[x$sds > 0]
  extra <- character()
  if (length(x$tasks) > 0L) {
    extra <- c(extra, paste0("; tasks: ", paste(x$tasks, collapse = ", ")))
  }
  if (length(x$covariates) > 0L) {
    extra <- c(extra, paste0(
      "; covariates: ", paste(names(x$covariates), collapse = ", ")
    ))
  }
  nonzero <- x$truth$cor[x$truth$cor$true_value != 0, , drop = FALSE]
  if (!is.null(nonzero) && nrow(nonzero) > 0L) {
    extra <- c(extra, paste0(
      "; nonzero correlations: ",
      paste0(
        nonzero$term, " = ", signif(nonzero$true_value, 3),
        collapse = ", "
      )
    ))
  }
  cat(
    "<bmmtools_simulation> ", x$model$name %||% class(x$model)[[2L]], ": ",
    x$n_subjects, " subjects, ", x$n_trials, " trials, ",
    nrow(x$data), " rows; varying: ",
    if (length(varying) == 0L) "none" else paste(varying, collapse = ", "),
    extra,
    "\n",
    sep = ""
  )
  invisible(x)
}

#' The default recovery formula: every free parameter gets a random intercept
#'
#' The formula the validation scripts in bmm converged on (decision 5):
#' `<parameter> ~ 1 + (1 | id)` for every parameter the model estimates,
#' and nothing for the ones it fixes. [recovery_grid()] uses it when no
#' formula is given; a user who wants to change one term can start
#' from it.
#'
#' @param model A `bmmodel`.
#' @param group The grouping variable, `"id"` by default.
#' @param re_cor Which random effects are correlated. `"none"` gives
#'   independent random effects: `(1 | id)`, or `(0 + task || id)` with
#'   tasks. `"within"` correlates the tasks of each parameter,
#'   `(0 + task | id)`, and needs `task_col`. `"all"` gives one correlation
#'   matrix across every parameter (and task), `(1 | p | id)` or
#'   `(0 + task | p | id)`, which is what the `model` estimator of
#'   correlation recovery reads.
#' @param task_col `NULL`, or the task column of a simulation with `tasks`
#'   (see [simulate_recovery()]). Each free parameter then gets one
#'   population value per task, `<parameter> ~ 0 + task`, and one
#'   subject-level effect per task: cell means, whose terms match the truth
#'   of the simulation.
#'
#' @return A `bmmformula`.
#'
#' @examples
#' \dontrun{
#' recovery_formula(bmm::mixture2p(resp_error = "y"))
#' recovery_formula(bmm::mixture2p(resp_error = "y"), re_cor = "all")
#' recovery_formula(
#'   bmm::mixture2p(resp_error = "y"),
#'   re_cor = "within", task_col = "task"
#' )
#' }
#'
#' @export
recovery_formula <- function(model,
                             group = "id",
                             re_cor = c("none", "within", "all"),
                             task_col = NULL) {
  check_model(model)
  re_cor <- rlang::arg_match(re_cor)
  valid_col <- is.character(task_col) && length(task_col) == 1L &&
    !is.na(task_col) && make.names(task_col) == task_col &&
    !grepl("[_.]", task_col)
  if (!is.null(task_col) && !valid_col) {
    cli::cli_abort(
      "{.arg task_col} must be {.code NULL} or a single syntactic name \\
       without {.code _} or {.code .}."
    )
  }
  rlang::check_installed("bmm", "to build a bmm formula.")
  free <- model_parameters(model)$free
  if (identical(re_cor, "within") && is.null(task_col)) {
    cli::cli_inform(
      "Without {.arg task_col} each parameter has a single random intercept, \\
       so there is nothing to correlate within a parameter; using \\
       {.code re_cor = \"none\"}."
    )
    re_cor <- "none"
  }
  if (identical(re_cor, "all") && length(free) < 2L) {
    if (is.null(task_col)) {
      cli::cli_inform(
        "The model estimates only one parameter, so there is nothing to \\
         correlate; using {.code re_cor = \"none\"}."
      )
      re_cor <- "none"
    } else {
      cli::cli_inform(
        "The model estimates only one parameter, so correlating across \\
         parameters is correlating its tasks; using \\
         {.code re_cor = \"within\"}."
      )
      re_cor <- "within"
    }
  }
  effects <- if (is.null(task_col)) "1" else paste0("0 + ", task_col)
  bar <- switch(re_cor,
    none = if (is.null(task_col)) " | " else " || ",
    within = " | ",
    all = " | p | "
  )
  formulas <- lapply(free, function(p) {
    stats::as.formula(
      paste0(p, " ~ ", effects, " + (", effects, bar, group, ")"),
      env = globalenv()
    )
  })
  do.call(bmm::bmf, formulas)
}
