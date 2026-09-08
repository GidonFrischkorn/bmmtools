# A mock fitter for the cache tests.
#
# fit_cached() takes the fitting function through `.fitter` so that its
# reuse logic can be tested without bmm, brms or a compiler. The mock
# records what it was called with and how often, and returns an object
# that saveRDS() round-trips.

#' Build a mock fitter that counts its calls
#'
#' @return A list with `fitter`, the function to inject, and `calls`, an
#'   environment whose `n` is the number of times it was called and whose
#'   `last` holds the arguments of the last call.
#' @noRd
mock_fitter <- function() {
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  calls$last <- NULL
  fitter <- function(formula, data, model, prior = NULL, ...) {
    calls$n <- calls$n + 1L
    calls$last <- list(
      formula = formula, data = data, model = model, prior = prior,
      dots = list(...)
    )
    # only names go into the saved object: a closure among the dots
    # would be serialised with its environment
    structure(
      list(call_number = calls$n, dot_names = names(list(...))),
      class = "mockfit"
    )
  }
  list(fitter = fitter, calls = calls)
}

#' A formula-like object with the shape of a bmmformula
#'
#' A named list of formulas carrying the class name, so the key logic
#' sees what it would see from bmm without bmm being installed.
#'
#' @noRd
fake_bmmformula <- function(...) {
  structure(list(...), class = "bmmformula")
}

#' A small data set and a stand-in model object
#' @noRd
fake_data <- function(n = 6L) {
  data.frame(id = factor(rep(1:2, each = n / 2)), y = seq_len(n) / 10)
}

fake_model <- function(name = "toy") {
  structure(list(name = name, links = list(a = "log")), class = "bmmodel")
}
