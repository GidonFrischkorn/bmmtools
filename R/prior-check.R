# Prior predictive checks on the observable scale.
#
# A prior on kappa is uninterpretable; the same prior read as "89% of
# prior draws imply chance-level responding" is not. The draws come from
# bmm(sample_prior = "only") through fit_cached(), not from an R-side
# parser of prior strings, so the prior being checked is exactly the
# prior bmm fits with, for every prior brms accepts.

#' Normalise the `prior` argument into a named list of sets
#'
#' A `brmsprior` is a data frame and a data frame is a list, so the class
#' is tested before the list.
#'
#' @noRd
prior_sets <- function(prior, call = rlang::caller_env()) {
  if (is.null(prior)) {
    return(list(default = NULL))
  }
  if (inherits(prior, "brmsprior")) {
    return(list(user = prior))
  }
  if (is.data.frame(prior) || !is.list(prior)) {
    cli::cli_abort(
      "{.arg prior} must be a {.cls brmsprior}, a named list of them, \\
       or {.code NULL}, not {.obj_type_friendly {prior}}.",
      call = call
    )
  }
  names(prior) <- names(prior) %||% rep("", length(prior))
  if (length(prior) == 1L && !nzchar(names(prior))) {
    names(prior) <- "user"
  }
  if (any(!nzchar(names(prior)))) {
    cli::cli_abort(
      c(
        "Every element of {.arg prior} needs a name.",
        i = "The name labels every row the set contributes."
      ),
      call = call
    )
  }
  reserved <- intersect(names(prior), c("response", "statistic", "difference"))
  if (length(reserved) > 0L) {
    cli::cli_abort(
      "A prior set cannot be called {.val {reserved}}: {.fn summary} makes \\
       one column per set and those names are taken.",
      call = call
    )
  }
  prior
}

#' Prior-predictive draws from a fit
#'
#' An internal generic rather than a direct call to
#' `brms::posterior_predict()`. bmm's mock backend returns a `bmmfit`
#' whose `$fit` is a bare number, so `posterior_predict()` cannot run on
#' it (measured 2026-09-08); without an injection point the prediction
#' step would either go untested or force the suite to compile Stan.
#'
#' @noRd
prior_predict_draws <- function(fit, ndraws, seed = NULL, ...) {
  UseMethod("prior_predict_draws")
}

#' @noRd
#' @export
prior_predict_draws.default <- function(fit, ndraws, seed = NULL, ...) {
  rlang::check_installed("brms", "to draw from the prior predictive.")
  with_seed_if(seed, brms::posterior_predict(fit, ndraws = ndraws))
}

#' Shape what the prediction step returned into one matrix per response
#'
#' The array branch is written from brms's documented multivariate shape,
#' not from a measurement: no model in the adapter table was observed
#' returning one.
#'
#' @noRd
yrep_list <- function(yrep, model, call = rlang::caller_env()) {
  if (is.matrix(yrep)) {
    name <- unlist(model$resp_vars, use.names = FALSE)[[1L]] %||% "response"
    return(stats::setNames(list(yrep), name))
  }
  if (is.array(yrep) && length(dim(yrep)) == 3L) {
    k <- dim(yrep)[[3L]]
    names_k <- dimnames(yrep)[[3L]] %||% paste0("response_", seq_len(k))
    out <- lapply(seq_len(k), function(i) yrep[, , i, drop = TRUE])
    return(stats::setNames(out, names_k))
  }
  cli::cli_abort(
    c(
      "The prior-predictive draws must be a matrix or a 3-dimensional \\
       array, not {.obj_type_friendly {yrep}}.",
      i = "{.fn prior_predict_draws} returns what \\
           {.fn brms::posterior_predict} gives for the model."
    ),
    call = call
  )
}

#' Floor and ceiling of a model's response, or `NA` when it has none
#'
#' A table keyed on the model class, like [generator_for()]. `ceiling`
#' may be a vector along the rows of `data`, which is what a count model
#' needs: the ceiling is that row's number of trials.
#'
#' @noRd
response_range <- function(model, data, call = rlang::caller_env()) {
  none <- function(why) {
    cli::cli_inform(c(
      why,
      i = "Give {.arg range}, or a {.arg summary} of your own, for floor \\
           and ceiling rates."
    ))
    list(floor = NA_real_, ceiling = NA_real_)
  }
  key <- adapter_name(model)
  if (is.na(key) || key %in% c("ddm", "ezdm")) {
    return(none(
      "No response range is known for a model of class {.cls {class(model)}}."
    ))
  }
  if (key %in% c("mixture2p", "sdm")) {
    return(list(floor = -pi, ceiling = pi))
  }
  # sdt_yn and sdt_mafc: a count out of the trials in that row
  column <- model$other_vars$n_trials
  if (is.null(column) || !is.character(column) || !column %in% names(data)) {
    return(none(
      "The column {.val {column %||% 'n_trials'}} is not in {.arg data}, \\
       so the ceiling of the response count is unknown."
    ))
  }
  list(floor = 0, ceiling = as.double(data[[column]]))
}

#' Accept a range as a pair or as a list, and give back a list
#' @noRd
as_range <- function(range, call = rlang::caller_env()) {
  if (is.list(range)) {
    return(list(floor = range$floor, ceiling = range$ceiling))
  }
  # `c(NA, NA)` is logical, and "no bounds are known" is a range a caller
  # may reasonably state
  numeric_enough <- is.numeric(range) || all(is.na(range))
  if (!numeric_enough || length(range) != 2L) {
    cli::cli_abort(
      "{.arg range} must be a numeric vector of two values, the floor and \\
       the ceiling of the response, or {.code NULL}.",
      call = call
    )
  }
  list(floor = range[[1L]], ceiling = range[[2L]])
}

#' The share of draws at or beyond a bound that may vary by column
#' @noRd
bound_rate <- function(yrep, bound, side) {
  if (is.null(bound) || all(is.na(bound))) {
    return(NA_real_)
  }
  limit <- matrix(bound, nrow = nrow(yrep), ncol = ncol(yrep), byrow = TRUE)
  hit <- if (identical(side, "floor")) yrep <= limit else yrep >= limit
  mean(hit, na.rm = TRUE)
}

#' The default prior-predictive summary
#'
#' Returns the `function(yrep, data)` the user contract names; the range
#' and the model reach it through the closure rather than through extra
#' arguments.
#'
#' For a discrete response the two rates are the share of
#' prior-predictive draws implying floor or ceiling performance. For a
#' continuous response exact equality to a boundary has probability zero,
#' so they report the mass outside the support, which is `0` unless the
#' model is wrong, and the quantiles are what carries the information.
#'
#' @noRd
default_prior_summary <- function(range, model) {
  range <- as_range(range)
  function(yrep, data) {
    pieces <- lapply(names(yrep_list(yrep, model)), function(name) {
      values <- yrep_list(yrep, model)[[name]]
      quantiles <- stats::quantile(
        values, c(0.5, 0.9, 0.95),
        na.rm = TRUE, names = FALSE
      )
      tibble::tibble(
        response = name,
        statistic = c("floor_rate", "ceiling_rate", "q50", "q90", "q95"),
        value = c(
          bound_rate(values, range$floor, "floor"),
          bound_rate(values, range$ceiling, "ceiling"),
          quantiles
        )
      )
    })
    dplyr::bind_rows(pieces)
  }
}

#' Check that a summary function returned what the contract asks for
#' @noRd
check_summary_result <- function(out, set, call = rlang::caller_env()) {
  if (!is.data.frame(out)) {
    cli::cli_abort(
      "{.arg summary} must return a data frame, \\
       not {.obj_type_friendly {out}}.",
      call = call
    )
  }
  missing <- setdiff(c("statistic", "value"), names(out))
  if (length(missing) > 0L) {
    cli::cli_abort(
      c(
        "{.arg summary} must return the column{?s} {.val {missing}}.",
        i = "It returned {.val {names(out)}} for the prior set \\
             {.val {set}}."
      ),
      call = call
    )
  }
  invisible(out)
}

#' The class a prior row belongs to, as one string
#'
#' A coefficient row and the class-level row that covers it agree on
#' class, group, parameter and distributional parameter, and differ only
#' in `coef`. Pasting those four gives a key the two share. Columns brms
#' did not return count as empty, so a table without `nlpar` still keys.
#'
#' @noRd
prior_class_key <- function(priors) {
  cols <- c("class", "group", "nlpar", "dpar")
  parts <- lapply(cols, function(col) {
    if (col %in% names(priors)) as.character(priors[[col]]) else ""
  })
  do.call(paste, c(parts, list(sep = "\r")))
}

#' Warn about population-level slopes with no proper prior
#'
#' `sample_prior = "only"` samples every parameter, so a flat slope has
#' nothing to sample from. Whether brms refuses outright was not measured
#' --- establishing it needs a compile, which no test here may do --- so
#' this warns and lets the fitter have the last word.
#'
#' brms leaves a coefficient row's `prior` empty when the class-level row
#' of the same class carries the prior, so an empty string alone does not
#' mean flat. A coefficient is flat only when nothing at its own level and
#' nothing at its class level is proper.
#'
#' @noRd
check_improper_priors <- function(formula, data, model, sets) {
  if (!rlang::is_installed("bmm")) {
    return(invisible(NULL))
  }
  for (set in names(sets)) {
    priors <- tryCatch(
      as.data.frame(
        bmm::default_prior(formula, data, model, prior = sets[[set]])
      ),
      error = function(e) NULL
    )
    needed <- c("class", "coef", "prior")
    if (is.null(priors) || !all(needed %in% names(priors))) {
      next
    }
    proper <- nzchar(priors$prior) & priors$prior != "(flat)"
    key <- prior_class_key(priors)
    covered <- key %in% key[proper & !nzchar(priors$coef)]
    flat <- priors$class == "b" & nzchar(priors$coef) &
      priors$coef != "Intercept" & !proper & !covered
    parameter <- priors$nlpar[flat] %||% rep("", sum(flat))
    slopes <- ifelse(
      nzchar(parameter),
      paste0(priors$coef[flat], " (", parameter, ")"),
      priors$coef[flat]
    )
    if (length(slopes) > 0L) {
      cli::cli_warn(c(
        "The prior set {.val {set}} leaves {.val {slopes}} without a \\
         proper prior.",
        i = "{.code sample_prior = \"only\"} samples every parameter; a \\
             flat prior has nothing to sample from."
      ))
    }
  }
  invisible(NULL)
}

#' A file name that survives being a prior set's name
#' @noRd
set_path <- function(base, set) {
  paste0(base, "-", gsub("[^A-Za-z0-9._-]+", "_", set))
}

#' Fit one prior set from the prior alone
#' @noRd
prior_fit <- function(set, prior, formula, data, model, base, refit, seed,
                      dots, fitter, call = rlang::caller_env()) {
  args <- c(
    list(
      formula = formula, data = data, model = model,
      file = set_path(base, set), prior = prior, refit = refit,
      sample_prior = "only", .fitter = fitter
    ),
    if (!is.null(seed)) list(seed = seed),
    dots
  )
  withCallingHandlers(
    tryCatch(
      rlang::exec(fit_cached, !!!args),
      error = function(e) {
        cli::cli_abort(
          c(
            "Fitting the prior set {.val {set}} failed.",
            i = "{.code sample_prior = \"only\"} needs a proper prior on \\
                 every parameter it samples."
          ),
          parent = e, call = call
        )
      }
    ),
    bmmtools_cache_message = function(m) NULL
  )
}

#' Check what the priors say the data should look like
#'
#' The question to ask before fitting anything: a prior on a link-scale
#' parameter is not interpretable, and the same prior expressed as "this
#' share of prior draws implies chance-level responding" is. Draws come
#' from `bmm(sample_prior = "only")` through [fit_cached()], so the prior
#' being checked is exactly the prior bmm fits with --- exact for every
#' prior brms accepts, and it gives predictions on the observable scale
#' for free.
#'
#' @param model A `bmmodel`, built with the column names `data` uses.
#' @param formula The `bmmformula` you intend to fit. The prior a model
#'   needs depends on it, so there is no default.
#' @param data The design the prior is checked over. The response column
#'   must exist because brms builds its Stan data from it; its values do
#'   not enter the prior draws, and are used only as the observed
#'   reference in [plot_prior_check()].
#' @param prior `NULL` checks bmm's own defaults. A `brmsprior` checks
#'   that one. A **named** list runs one fit per element and compares
#'   them, a `NULL` element meaning bmm's defaults, so
#'   `list(default = NULL, tight = my_prior)` is the usual comparison.
#' @param summary A function `(yrep, data)` returning a data frame with
#'   at least `statistic` and `value`, and optionally `response`. `NULL`
#'   uses the default: the share of draws at the floor and at the ceiling
#'   of the response, and the 50th, 90th and 95th percentiles of the
#'   prior-predicted observable.
#' @param n_draws Prior-predictive draws kept per prior set. More than
#'   the fit holds is reported and capped, not refused.
#' @param range The floor and the ceiling of the response, as two
#'   numbers. `NULL` derives them from the model where they are known:
#'   `0` and the trial count for `sdt_yn` and `sdt_mafc`, `-pi` and `pi`
#'   for `mixture2p` and `sdm`. Where they are not known the two rates
#'   are `NA` rather than invented.
#' @param file Where to cache the fits. `NULL` uses a temporary file, so
#'   the compile is reused within the session and nothing is left behind;
#'   a path caches durably and each set gets `<file>-<set>`.
#' @param refit Passed to [fit_cached()].
#' @param seed Passed to the fitter, where it enters the cache key, and
#'   used around the draw subsample so that the same seed keeps the same
#'   `n_draws`. `NULL` leaves the random number generator alone and is
#'   recorded as `NA`.
#' @param ... Passed to [fit_cached()] and on to the fitter: `chains`,
#'   `iter`, `backend`, `cores`. `sample_prior` is set here and giving it
#'   is an error.
#' @param .fitter The fitting function, `bmm::bmm()` by default. Tests
#'   inject a stand-in so that no model is compiled.
#'
#' @return A `bmmtools_prior_check` tibble, one row per prior set,
#'   response and statistic, with the columns `prior`, `response`,
#'   `statistic` and `value`. It carries the draws themselves as an
#'   attribute, because a summary cannot be plotted as a distribution,
#'   along with `brms::prior_summary()` of each fit --- the prior that
#'   actually produced the draws, rather than the one that was asked for.
#'
#' @details
#' A population-level slope left without a proper prior is warned about
#' before fitting: `sample_prior = "only"` samples every parameter and a
#' flat prior has nothing to sample from. It is a warning rather than an
#' error because bmm's defaults cover the intercepts and the group-level
#' SDs, and only a slope you added is likely to be missing.
#'
#' A floor or ceiling rate is meaningful for a **discrete** response --- a
#' count out of `n_trials` --- where it is the share of prior draws
#' implying perfect or floor performance. For a **continuous** response
#' exact equality to a boundary has probability zero and the quantiles
#' are what to read. A circular response error is best read through a
#' summary of your own on the absolute error, because the quantiles of a
#' signed error are symmetric about zero and say nothing:
#'
#' ```r
#' prior_check(model, formula, data, summary = function(yrep, data) {
#'   tibble::tibble(
#'     statistic = c("q50", "q90"),
#'     value = unname(stats::quantile(abs(yrep), c(0.5, 0.9)))
#'   )
#' })
#' ```
#'
#' @examples
#' \dontrun{
#' model <- bmm::mixture2p(resp_error = "y")
#' checked <- prior_check(
#'   model, recovery_formula(model), my_data,
#'   n_draws = 500, seed = 1
#' )
#' }
#'
#' # a precomputed check comparing two prior sets
#' prior_check_sdt_yn
#' summary(prior_check_sdt_yn)
#'
#' @export
prior_check <- function(model,
                        formula,
                        data,
                        prior = NULL,
                        summary = NULL,
                        n_draws = 1000,
                        range = NULL,
                        file = NULL,
                        refit = c("on_change", "never", "always"),
                        seed = NULL,
                        ...,
                        .fitter = NULL) {
  check_model(model)
  refit <- rlang::arg_match(refit)
  if (!is.data.frame(data)) {
    cli::cli_abort(
      "{.arg data} must be a data frame, not {.obj_type_friendly {data}}."
    )
  }
  if (nrow(data) == 0L) {
    cli::cli_abort("{.arg data} must have at least one row.")
  }
  n_draws <- check_count(n_draws, "n_draws")
  if (!is.null(summary) && !is.function(summary)) {
    cli::cli_abort(
      "{.arg summary} must be a function {.code (yrep, data)} or \\
       {.code NULL}, not {.obj_type_friendly {summary}}."
    )
  }
  if (!is.null(summary) && length(formals(summary)) < 2L) {
    cli::cli_abort(
      "{.arg summary} must take two arguments, the draws and the data."
    )
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L)) {
    cli::cli_abort("{.arg seed} must be a single number or {.code NULL}.")
  }
  dots <- rlang::list2(...)
  if ("sample_prior" %in% names(dots)) {
    cli::cli_abort(
      c(
        "{.arg sample_prior} cannot be passed to {.fn prior_check}.",
        i = "It is set to {.val only}: that is what makes this a prior
             check."
      )
    )
  }

  if (is.null(range)) range <- response_range(model, data)
  range <- as_range(range)
  summarise <- summary %||% default_prior_summary(range, model)
  sets <- prior_sets(prior)
  check_improper_priors(formula, data, model, sets)

  base <- file %||% tempfile(pattern = "bmmtools-prior-check-")
  draws <- list()
  priors <- list()
  pieces <- list()
  for (set in names(sets)) {
    fit <- prior_fit(
      set, sets[[set]], formula, data, model, base, refit, seed, dots,
      .fitter
    )
    yrep <- prior_predict_draws(fit, ndraws = n_draws, seed = seed)
    shaped <- yrep_list(yrep, model)
    kept <- nrow(shaped[[1L]])
    if (kept < n_draws) {
      cli::cli_inform(
        "The prior set {.val {set}} has {kept} draw{?s}, fewer than the \\
         {n_draws} asked for."
      )
    }
    rows <- check_summary_result(summarise(yrep, data), set)
    if (!"response" %in% names(rows)) {
      rows$response <- if (length(shaped) == 1L) {
        names(shaped)[[1L]]
      } else {
        NA_character_
      }
    }
    rows$prior <- set
    draws[[set]] <- yrep
    priors[[set]] <- tryCatch(
      brms::prior_summary(fit),
      error = function(e) NULL
    )
    pieces[[set]] <- rows
  }

  new_bmmtools_prior_check(
    dplyr::bind_rows(pieces),
    sets = names(sets),
    draws = draws,
    priors = priors,
    data = data,
    model = model,
    formula = formula,
    n_draws = n_draws,
    seed = if (is.null(seed)) NA_real_ else as.double(seed),
    range = range,
    files = vapply(names(sets), function(s) set_path(base, s), character(1)),
    call = match.call()
  )
}
