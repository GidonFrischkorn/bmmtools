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
#' @noRd
check_pars <- function(pars, model, call = rlang::caller_env()) {
  info <- model_parameters(model)
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
#' @return A matrix, subjects by parameters, on the link scale.
#' @noRd
check_subject_pars <- function(subject_pars, pars, sds, n_subjects,
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
  varying <- names(sds)[sds > 0]
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
    pars,
    nrow = n_subjects, ncol = length(pars), byrow = TRUE,
    dimnames = list(ids, names(pars))
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
    cli::cli_abort(
      "{.arg subject_pars} must give the varying parameter{?s} \\
       {.val {missing}}.",
      call = call
    )
  }
  values
}

#' Draw subject values on the link scale
#'
#' @return A matrix, subjects by parameters.
#' @noRd
draw_subject_pars <- function(pars, sds, n_subjects) {
  values <- vapply(names(pars), function(p) {
    pars[[p]] + stats::rnorm(n_subjects, 0, sds[[p]])
  }, numeric(n_subjects))
  values <- matrix(
    values,
    nrow = n_subjects,
    dimnames = list(as.character(seq_len(n_subjects)), names(pars))
  )
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
check_generated <- function(data, model, call = rlang::caller_env()) {
  if (!is.data.frame(data)) {
    cli::cli_abort(
      "The generator must return a data frame, \\
       not {.obj_type_friendly {data}}.",
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
truth_tables <- function(pars, sds, values) {
  varying <- names(sds)[sds > 0]
  population <- tibble::tibble(
    term = names(pars),
    true_value = as.double(unname(pars))
  )
  subjects <- tibble::tibble(
    id = rep(rownames(values), times = length(varying)),
    term = rep(varying, each = nrow(values)),
    true_value = as.double(values[, varying, drop = TRUE])
  )
  if (length(varying) == 0L) {
    subjects <- tibble::tibble(
      id = character(), term = character(), true_value = double()
    )
  }
  list(population = population, subjects = subjects)
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
#'   that value.
#' @param n_subjects,n_trials Subjects, and trials per subject and per
#'   row of the generator's layout (for `sdt_yn`, per stimulus class).
#' @param sds Named numeric: between-subject standard deviations on the
#'   link scale. A parameter not named does not vary. `NULL` means no
#'   parameter varies.
#' @param subject_pars A data frame `id`, `term`, `true_value` (link
#'   scale) of subject values to use instead of drawing them, so that
#'   replications can share the same simulated people.
#' @param generator A function `(pars, n_trials, model)` returning one
#'   subject's rows as a data frame with the model's column names;
#'   `pars` is a named list on the natural scale, fixed parameters
#'   included. `NULL` uses the adapter bmmtools ships for the model.
#' @param seed A seed applied with `withr::with_seed()` around the draws
#'   and the generator; `NULL` leaves the random number generator alone.
#'
#' @return A list of class `bmmtools_simulation` with `data` (a tibble,
#'   `id` first), `truth` (a list of two tibbles, `population` with
#'   `term` and `true_value`, `subjects` with `id`, `term` and
#'   `true_value`, both on the link scale), the validated `pars` and
#'   `sds`, `n_subjects`, `n_trials`, `seed` (`NA` when none), `model`
#'   and `generator`.
#'
#' @details
#' Adapters exist for `sdt_yn`, `sdt_mafc`, `ezdm` (three parameters),
#' `ddm`, `mixture2p` and `sdm`. Every other model takes a `generator`.
#' The truth for the subjects lists only the parameters that vary,
#' because a parameter that does not vary has nothing person-level to
#' recover.
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
#' }
#'
#' @export
simulate_recovery <- function(model,
                              pars,
                              n_subjects,
                              n_trials,
                              sds = NULL,
                              subject_pars = NULL,
                              generator = NULL,
                              seed = NULL) {
  check_model(model)
  pars <- check_pars(pars, model)
  n_subjects <- check_count(n_subjects, "n_subjects")
  n_trials <- check_count(n_trials, "n_trials")
  sds <- check_sds(sds, pars, model)
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

  supplied <- NULL
  if (!is.null(subject_pars)) {
    supplied <- check_subject_pars(subject_pars, pars, sds, n_subjects)
  }

  out <- with_seed_if(seed, {
    values <- supplied
    if (is.null(values)) values <- draw_subject_pars(pars, sds, n_subjects)
    pieces <- lapply(seq_len(n_subjects), function(i) {
      # keep the names when the matrix has a single column
      link_values <- stats::setNames(values[i, ], colnames(values))
      rows <- generator(natural_pars(link_values, model), n_trials, model)
      check_generated(rows, model)
      rows <- tibble::as_tibble(rows)
      tibble::add_column(
        rows,
        id = factor(rep(i, nrow(rows)), levels = seq_len(n_subjects)),
        .before = 1L
      )
    })
    list(data = dplyr::bind_rows(pieces), values = values)
  })

  structure(
    list(
      data = out$data,
      truth = truth_tables(pars, sds, out$values),
      pars = pars,
      sds = sds,
      n_subjects = n_subjects,
      n_trials = n_trials,
      seed = if (is.null(seed)) NA_real_ else as.double(seed),
      model = model,
      generator = generator_name
    ),
    class = "bmmtools_simulation"
  )
}

#' @export
print.bmmtools_simulation <- function(x, ...) {
  varying <- names(x$sds)[x$sds > 0]
  cat(
    "<bmmtools_simulation> ", x$model$name %||% class(x$model)[[2L]], ": ",
    x$n_subjects, " subjects, ", x$n_trials, " trials, ",
    nrow(x$data), " rows; varying: ",
    if (length(varying) == 0L) "none" else paste(varying, collapse = ", "),
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
#'
#' @return A `bmmformula`.
#'
#' @examples
#' \dontrun{
#' recovery_formula(bmm::mixture2p(resp_error = "y"))
#' }
#'
#' @export
recovery_formula <- function(model, group = "id") {
  check_model(model)
  rlang::check_installed("bmm", "to build a bmm formula.")
  free <- model_parameters(model)$free
  formulas <- lapply(free, function(p) {
    stats::as.formula(paste0(p, " ~ 1 + (1 | ", group, ")"), env = globalenv())
  })
  do.call(bmm::bmf, formulas)
}
