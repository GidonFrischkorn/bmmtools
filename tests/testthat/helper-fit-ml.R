# A mock fitter for fit_ml().
#
# fit_ml() runs ONE no-pooling fit whose per-subject estimates arrive as
# population coefficients named `<par>_<by><level>`. The mock returns an
# object shaped like that, so the formula reduction, the rename back to
# (term, id) and the convergence rule can be tested without a compiler.
# It registers extract_estimates.mlmockfit the way helper-generate.R
# registers extract_estimates.mockfit.

#' Build a mock fitter for `fit_ml()`
#'
#' @param estimates Optional named numeric vector keyed by the coefficient
#'   name (`"kappa_id1"`). Anything not named is 0.
#' @param sd The half-width used for the interval of every row.
#' @return A list with `fitter` and `calls`, as [mock_fitter()].
#' @noRd
ml_mock_fitter <- function(estimates = NULL, sd = 0.5) {
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  calls$last <- NULL
  fitter <- function(formula, data, model, prior = NULL, ...) {
    calls$n <- calls$n + 1L
    calls$last <- list(
      formula = formula, data = data, model = model, prior = prior,
      dots = list(...)
    )
    structure(
      list(
        coefficients = ml_mock_coefficients(formula, data),
        estimates = estimates,
        sd = sd
      ),
      class = c("mlmockfit", "mockfit")
    )
  }
  list(fitter = fitter, calls = calls)
}

#' The coefficient names a no-pooling fit of `formula` would carry
#'
#' `<par>_<by><level>` per parameter and level, which is what brms names
#' the columns of `model.matrix(~ 0 + id)`. With a single level there is
#' one coefficient per parameter and `coefficient_terms()` leaves the bare
#' parameter name, which is the edge case the rename has to handle.
#'
#' @noRd
ml_mock_coefficients <- function(formula, data) {
  pars <- names(formula)
  by <- all.vars(formula[[1L]][[3L]])[1L]
  levels <- levels(as.factor(data[[by]]))
  out <- expand.grid(
    level = levels, par = pars,
    stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE
  )
  out$term <- if (length(levels) == 1L) {
    out$par
  } else {
    paste0(out$par, "_", by, out$level)
  }
  out[, c("par", "level", "term")]
}

#' Population estimates for a mock no-pooling fit
#' @noRd
extract_estimates_mlmockfit <- function(fit,
                                        level = "population",
                                        group = NULL,
                                        ci_level = 0.95,
                                        ci_method = "eti",
                                        drop_constants = TRUE,
                                        converged = NULL,
                                        estimator = "bayes",
                                        ...) {
  if (!identical(level, "population")) {
    cli::cli_abort(
      "The fit has no group-level effects, so there are no subject-level \\
       estimates to extract."
    )
  }
  co <- fit$coefficients
  value <- rep(0, nrow(co))
  if (!is.null(fit$estimates)) {
    hit <- match(co$term, names(fit$estimates))
    value[!is.na(hit)] <- unname(fit$estimates[hit[!is.na(hit)]])
  }
  tibble::tibble(
    term = co$term,
    estimate = value,
    ci_low = value - fit$sd,
    ci_high = value + fit$sd,
    ci_method = "eti",
    ci_level = ci_level,
    rhat = 1,
    ess_bulk = 1000,
    ess_tail = 1000,
    level = "population",
    id = NA_character_,
    converged = TRUE,
    estimator = estimator
  )
}

registerS3method(
  "extract_estimates", "mlmockfit", extract_estimates_mlmockfit,
  envir = asNamespace("bmmtools")
)

#' A two-parameter stand-in model with the shape `model_parameters()` reads
#' @noRd
fake_ml_model <- function(free = c("kappa", "thetat"),
                          fixed = list(mu1 = 0)) {
  links <- as.list(stats::setNames(
    rep("log", length(free) + length(fixed)), c(free, names(fixed))
  ))
  structure(
    list(name = "toy", links = links, fixed_parameters = fixed),
    class = "bmmodel"
  )
}

#' A data set with `n` subjects and `n_trials` rows each
#' @noRd
fake_ml_data <- function(n = 4L, n_trials = 5L, by = "id", levels = NULL) {
  levels <- levels %||% as.character(seq_len(n))
  out <- data.frame(
    x = factor(rep(levels, each = n_trials), levels = levels),
    y = seq_len(length(levels) * n_trials) / 10
  )
  names(out)[1L] <- by
  out
}

# Helpers for the optim route.

#' An objective whose optimum and curvature are known in closed form
#'
#' `sum((theta - y)^2) / 2` on the link scale, so the minimum is `mean(y)`
#' and the Hessian is the number of rows --- giving a standard error of
#' `1 / sqrt(n)`. The toy model's link is `log`, so `log(pars$kappa)` is
#' `theta` and the objective is written through the natural-scale value it
#' is handed, exactly as a real density would be.
#'
#' @noRd
ml_square_nll <- function(par = "kappa") {
  function(pars, data, model) {
    sum((log(pars[[par]]) - data$y)^2) / 2
  }
}

#' One case per model that has both a generator and a density adapter
#'
#' Each case carries a `raw` closure that calls the same bmm density with
#' `log = FALSE`, so the adapter can be checked against it. That is the
#' check for the trap recorded in R/densities.R: dmixture2p(), dsdm(),
#' dsdt_yn() and dsdt_mafc() default to `log = FALSE` while dddm() and
#' dezdm() default to TRUE, so an adapter that leaves the argument out is
#' wrong for four of the six and right for the other two.
#'
#' The signal-detection cases appear only when the installed bmm has them;
#' CRAN bmm 1.3.2 does not, and the four that remain are what a runner
#' exercises.
#'
#' @noRd
ml_density_cases <- function() {
  cases <- list(
    list(
      label = "mixture2p", model = bmm::mixture2p(resp_error = "y"),
      link = c(kappa = log(8), thetat = stats::qlogis(0.75)),
      raw = function(pars, data, model) {
        bmm::dmixture2p(
          data$y,
          mu = pars$mu1, kappa = pars$kappa, p_mem = pars$thetat
        )
      }
    ),
    list(
      label = "sdm", model = bmm::sdm(resp_error = "y"),
      link = c(c = log(4), kappa = log(3.5)),
      raw = function(pars, data, model) {
        bmm::dsdm(data$y, mu = pars$mu, c = pars$c, kappa = pars$kappa)
      }
    ),
    list(
      label = "ddm", model = bmm::ddm(rt = "rt", response = "response"),
      link = c(
        drift = log(1.5), bound = log(1.2), ndt = log(0.3),
        zr = stats::qlogis(0.5)
      ),
      raw = function(pars, data, model) {
        exp(bmm::dddm(
          data$rt, data$response,
          drift = pars$drift, bound = pars$bound,
          ndt = pars$ndt, zr = pars$zr, log = TRUE
        ))
      }
    ),
    list(
      label = "ezdm",
      model = bmm::ezdm(
        mean_rt = "mean_rt", var_rt = "var_rt", n_upper = "n_upper",
        n_trials = "n_trials"
      ),
      link = c(drift = log(1.5), bound = log(1.2), ndt = log(0.3)),
      raw = function(pars, data, model) {
        exp(bmm::dezdm(
          data$mean_rt, data$var_rt, data$n_upper, data$n_trials,
          drift = pars$drift, bound = pars$bound, ndt = pars$ndt,
          s = pars$s, log = TRUE
        ))
      }
    )
  )
  if (!exists("sdt_yn", envir = asNamespace("bmm"), inherits = FALSE)) {
    return(cases)
  }
  c(cases, list(
    list(
      label = "sdt_yn",
      model = bmm::sdt_yn(response = "hits", stimulus = "stim", n_trials = "n"),
      link = c(d = 1.5, criterion = 0.2),
      raw = function(pars, data, model) {
        bmm::dsdt_yn(
          data$hits, data$n, data$stim,
          d = pars$d,
          criterion = pars$criterion, sdratio = pars$sdratio,
          dist = model$other_vars$dist
        )
      }
    ),
    list(
      label = "sdt_mafc",
      model = bmm::sdt_mafc(response = "k", n_trials = "n", m = 4L),
      link = c(d = 1.5),
      raw = function(pars, data, model) {
        bmm::dsdt_mafc(
          data$k, data$n,
          m = model$other_vars$m, d = pars$d,
          dist = model$other_vars$dist
        )
      }
    )
  ))
}

#' Population rows as `extract_estimates(level = "population")` returns them
#'
#' The input `ml_subject_rows()` renames. Only the columns it reads matter.
#'
#' @noRd
fake_ml_population <- function(terms) {
  tibble::tibble(
    term = terms,
    estimate = seq_along(terms) + 0,
    ci_low = 0, ci_high = 1, ci_method = "eti", ci_level = 0.95,
    rhat = 1, ess_bulk = 1000, ess_tail = 1000,
    level = "population", id = NA_character_, converged = TRUE,
    estimator = "bayes"
  )
}
