# True and estimated subject values side by side (spec 5, section 5.3;
# local/ARCHITECTURE.md decision 31).
#
# The table a structural equation model reads: one row per subject, the
# true value and the posterior point estimate of every parameter, and the
# covariates. The grid stores the per-subject means and medians in its
# sidecars, so a table over a whole grid needs no fit.

#' Posterior means and medians per subject, long, on the link scale
#'
#' @param subject_draws An array from `extract_subject_draws()`.
#' @param covariates `NULL` or a long covariate truth (`id`, `term`,
#'   `true_value`); each value becomes its own mean and median.
#' @return A tibble `id`, `term`, `mean`, `median`.
#' @noRd
subject_means_from_draws <- function(subject_draws, covariates = NULL,
                                     call = rlang::caller_env()) {
  subject_draws <- check_subject_draws(subject_draws, call = call)
  ids <- dimnames(subject_draws)[[3L]]
  terms <- dimnames(subject_draws)[[4L]]
  means <- colMeans(subject_draws, dims = 2L)
  medians <- apply(subject_draws, c(3L, 4L), stats::median)
  out <- tibble::tibble(
    id = rep(ids, times = length(terms)),
    term = rep(terms, each = length(ids)),
    mean = as.double(means),
    median = as.double(medians)
  )
  if (is.null(covariates) || nrow(covariates) == 0L) {
    return(out)
  }
  values <- as.double(covariates$true_value)
  dplyr::bind_rows(out, tibble::tibble(
    id = as.character(covariates$id),
    term = as.character(covariates$term),
    mean = values,
    median = values
  ))
}

#' @noRd
empty_subject_means <- function() {
  tibble::tibble(
    id = character(), term = character(), mean = double(), median = double()
  )
}

#' Label subject means with their cell and attach the true values
#'
#' A term's true value is the subject's value where the truth lists one,
#' the covariate value for a covariate, and otherwise the population value,
#' which is every subject's value of a parameter that does not vary.
#'
#' @return A tibble `condition`, `replication`, `id`, `term`, `covariate`,
#'   `mean`, `median`, `true_value`.
#' @noRd
label_subject_means <- function(means, truth, covariate_names, condition,
                                replication) {
  means <- dplyr::bind_rows(empty_subject_means(), means)
  key <- paste(means$id, means$term, sep = "\r")
  lookup <- function(table) {
    if (is.null(table) || nrow(table) == 0L) {
      return(rep(NA_real_, length(key)))
    }
    at <- match(key, paste(table$id, table$term, sep = "\r"))
    as.double(table$true_value)[at]
  }
  population <- truth$population
  from_population <- rep(NA_real_, nrow(means))
  if (!is.null(population) && nrow(population) > 0L) {
    from_population <- as.double(population$true_value)[
      match(means$term, population$term)
    ]
  }
  n <- nrow(means)
  tibble::tibble(
    condition = rep(as.character(condition), n),
    replication = rep(replication, n),
    id = means$id,
    term = means$term,
    covariate = means$term %in% covariate_names,
    mean = means$mean,
    median = means$median,
    true_value = dplyr::coalesce(
      lookup(truth$subjects), lookup(truth$covariates), from_population
    )
  )
}

#' True and estimated subject values, one row per subject
#'
#' Puts every subject's generating value next to its posterior point
#' estimate, one column pair per parameter, with the covariates alongside.
#' This is the input for a structural equation model of true against
#' estimated values, for example in lavaan, which recovers a correlation
#' that the posterior means attenuate by shrinkage (see
#' [extract_correlations()]). bmmtools does not fit that model itself and
#' does not depend on lavaan.
#'
#' @param x A `bmmtools_simulation`, together with `fit`; a
#'   `bmmtools_simulation_set` from [simulate_components()], together with
#'   its fits; or the result of [recovery_grid()], whose attribute
#'   `subject_means` carries the posterior means, medians and true values of
#'   every cell.
#' @param fit The fit of `x` when `x` is a simulation: a `brmsfit`, or any
#'   object with an [extract_subject_draws()] method. For a simulation set,
#'   the fits from [fit_components()], one per component, whose terms appear
#'   prefixed with the component name. Not used with a grid result.
#' @param point The point estimate: the posterior `"mean"` or `"median"`
#'   of each subject's value.
#' @param scale `"link"` or `"natural"`. On the natural scale the true
#'   values and the point estimates of each parameter go through
#'   [inverse_link()]; covariates are left as they are.
#' @param links A named character vector mapping a term to a link name.
#'   `NULL` reads the simulation model's link table, or the one the grid
#'   stored; without one, the natural scale falls back to the link scale
#'   with a message.
#'
#' @return A tibble with one row per condition, replication and subject:
#'   `condition` (`NA` for a single simulation), `replication` (`1` for a
#'   single simulation), `id`, then `true_<term>` and `est_<term>` for each
#'   parameter with subject-level draws, then one column per covariate,
#'   under its own name.
#'
#' @details
#' A parameter that does not vary in the simulation has no subject truth;
#' its `true_<term>` column is the population value for every subject.
#'
#' On the natural scale the estimate is the inverse link of the point
#' estimate on the link scale, not the mean of the inverse links, as for
#' the `point` estimator of [extract_correlations()]. For the median the
#' two are the same.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_recovery(
#'   bmm::mixture2p(resp_error = "y"),
#'   pars = c(kappa = log(8), thetat = qlogis(0.75)),
#'   n_subjects = 100, n_trials = 50,
#'   sds = c(kappa = 0.3, thetat = 0.5),
#'   covariates = list(G = c(mean = 0, sd = 1)), seed = 1
#' )
#' fit <- bmm::bmm(recovery_formula(sim$model), sim$data, sim$model)
#' subject_table(sim, fit)
#'
#' # one table over a grid, from its sidecars
#' out <- recovery_grid(model, grid, pars, dir = "fits", sds = sds)
#' subject_table(out, point = "median")
#' }
#'
#' @export
subject_table <- function(x,
                          fit = NULL,
                          point = c("mean", "median"),
                          scale = c("link", "natural"),
                          links = NULL) {
  error_call <- rlang::current_env()
  point <- rlang::arg_match(point)
  scale <- rlang::arg_match(scale)

  if (inherits(x, "bmmtools_simulation")) {
    if (is.null(fit)) {
      cli::cli_abort(
        "A simulation needs its {.arg fit} to read the estimates from.",
        call = error_call
      )
    }
    means <- subject_means_from_draws(
      extract_subject_draws(fit), x$truth$covariates,
      call = error_call
    )
    long <- label_subject_means(
      means, x$truth, names(x$covariates),
      condition = NA_character_, replication = 1L
    )
    model_links <- x$model$links
  } else if (inherits(x, "bmmtools_simulation_set")) {
    # separate fits, bound subject by subject under prefixed terms
    is_set <- !is.null(fit) && is_fit_set(fit, call = error_call)
    if (!is_set) {
      cli::cli_abort(
        c(
          "A simulation set needs its fits in {.arg fit}, one per component.",
          i = "They come from {.fn fit_components}."
        ),
        call = error_call
      )
    }
    check_set_components(x, list(fit), call = error_call)
    means <- subject_means_from_draws(
      set_subject_draws(fit, call = error_call), x$truth$covariates,
      call = error_call
    )
    long <- label_subject_means(
      means, x$truth, names(x$covariates),
      condition = NA_character_, replication = 1L
    )
    model_links <- x$links
  } else if (!is.null(attr(x, "subject_means"))) {
    if (!is.null(fit)) {
      cli::cli_abort(
        c(
          "{.arg fit} is only used with a simulation.",
          i = "A grid result carries its subject means already."
        ),
        call = error_call
      )
    }
    long <- attr(x, "subject_means")
    model_links <- attr(long, "links")
  } else if (inherits(x, "bmmtools_recovery")) {
    cli::cli_abort(
      c(
        "{.arg x} carries no subject means.",
        i = "They come from {.fn recovery_grid}; for a single fit pass its \\
             simulation and {.arg fit}."
      ),
      call = error_call
    )
  } else {
    cli::cli_abort(
      "{.arg x} must be a {.cls bmmtools_simulation} or the result of \\
       {.fn recovery_grid}, not {.obj_type_friendly {x}}.",
      call = error_call
    )
  }
  if (nrow(long) == 0L) {
    cli::cli_abort(
      "{.arg x} has no subject means: no cell had a parameter that varies.",
      call = error_call
    )
  }

  resolved <- resolve_links(NULL, links %||% model_links, scale,
    call = error_call
  )
  subject_table_wide(long, point, resolved)
}

#' The wide table from the long subject means
#' @noRd
subject_table_wide <- function(long, point, resolved) {
  estimate <- long[[point]]
  true_value <- long$true_value
  if (identical(resolved$scale, "natural")) {
    for (term in unique(long$term[!long$covariate])) {
      at <- long$term == term & !long$covariate
      link <- link_of(term, resolved$links)
      estimate[at] <- inverse_link(estimate[at], link)
      true_value[at] <- inverse_link(true_value[at], link)
    }
  }

  key <- paste(long$condition, long$replication, long$id, sep = "\r")
  first <- !duplicated(key)
  out <- tibble::tibble(
    condition = as.character(long$condition[first]),
    replication = long$replication[first],
    id = as.character(long$id[first])
  )
  rows <- key[first]
  column <- function(values, term) {
    at <- long$term == term
    values[at][match(rows, key[at])]
  }
  for (term in unique(long$term[!long$covariate])) {
    out[[paste0("true_", term)]] <- column(true_value, term)
    out[[paste0("est_", term)]] <- column(estimate, term)
  }
  for (term in unique(long$term[long$covariate])) {
    out[[term]] <- column(true_value, term)
  }
  out
}
