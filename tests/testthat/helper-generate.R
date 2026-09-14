# Helpers for the generate-layer and grid tests.
#
# The grid is tested end to end without a sampler: the mock fitter from
# helper-cache.R stands in for bmm::bmm(), and the method below turns
# what it returns into an estimates tibble, so extract_estimates(),
# recover() and summary() all run on it.

#' A mock fitter whose fits remember the model and the subject ids
#' @noRd
grid_mock_fitter <- function(fail_on = NULL) {
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  calls$log <- list()
  fitter <- function(formula, data, model, prior = NULL, ...) {
    calls$n <- calls$n + 1L
    dots <- list(...)
    calls$log[[calls$n]] <- list(
      n_rows = nrow(data), ids = levels(data$id), seed = dots$seed,
      chains = dots$chains, iter = dots$iter
    )
    if (!is.null(fail_on) && fail_on(calls$n, data, dots)) {
      stop("mock sampler failure")
    }
    structure(
      list(
        parameters = {
          p <- bmm::parameters(model)
          p$parameter[!p$fixed]
        },
        ids = levels(data$id)
      ),
      class = "mockfit"
    )
  }
  list(fitter = fitter, calls = calls)
}

#' Estimates for a mock fit: 0 with a wide interval, converged
#' @noRd
extract_estimates_mockfit <- function(fit,
                                      level = c("population", "subject"),
                                      group = NULL,
                                      ci_level = 0.95,
                                      ci_method = "eti",
                                      drop_constants = TRUE,
                                      converged = NULL,
                                      ...) {
  if (missing(level)) level <- "population"
  one <- function(term, id, lvl) {
    tibble::tibble(
      term = term, estimate = 0, ci_low = -10, ci_high = 10,
      ci_method = "eti", ci_level = ci_level, rhat = 1, ess_bulk = 1000,
      ess_tail = 1000, level = lvl, id = id, converged = TRUE
    )
  }
  pieces <- list()
  if ("population" %in% level) {
    pieces$population <- one(fit$parameters, NA_character_, "population")
  }
  if ("subject" %in% level) {
    pieces$subject <- dplyr::bind_rows(lapply(fit$ids, function(i) {
      one(fit$parameters, i, "subject")
    }))
  }
  dplyr::bind_rows(pieces)
}

registerS3method(
  "extract_estimates", "mockfit", extract_estimates_mockfit,
  envir = asNamespace("bmmtools")
)

#' A long subject truth (`id`, `term`, `true_value`) as one column per term
#' @noRd
subjects_wide <- function(subjects) {
  terms <- unique(subjects$term)
  out <- lapply(terms, function(t) subjects$true_value[subjects$term == t])
  stats::setNames(as.data.frame(out), terms)
}

#' A two-row grid over subjects and trials
#' @noRd
small_grid <- function() {
  data.frame(n_subjects = c(3L, 4L), n_trials = c(10L, 20L))
}
