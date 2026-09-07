# Hand-built draws objects for the score-layer tests.
#
# estimates_from_draws() is the workhorse behind extract_estimates(), and
# it takes a draws object rather than a fit precisely so that structures
# the saved fixture does not contain --- two grouping factors, a chain
# with a missing draw, a model with no group-level effects --- can be
# tested without compiling any Stan.

#' Build a draws_array from a named list of variables
#'
#' Each element of `values` is either a single number (a constant, which
#' is what a parameter bmm fixes looks like in the draws) or a numeric
#' vector of length `n_iter * n_chain`.
#'
#' @noRd
fake_draws <- function(values, n_iter = 40L, n_chain = 2L) {
  n <- n_iter * n_chain
  cols <- lapply(values, function(v) {
    if (length(v) == 1L) {
      return(rep(as.double(v), n))
    }
    stopifnot(length(v) == n)
    as.double(v)
  })
  arr <- array(
    unlist(cols, use.names = FALSE),
    dim = c(n_iter, n_chain, length(cols)),
    dimnames = list(iteration = NULL, chain = NULL, variable = names(cols))
  )
  posterior::as_draws_array(arr)
}

#' The saved mixture2p fixture
#'
#' Built once by `fixtures/make-fixtures.R`; the suite reads it and
#' compiles nothing.
#'
#' @noRd
mixture2p_fit <- function() {
  readRDS(test_path("fixtures", "mixture2p-fit.rds"))
}

#' The generating values that go with the fixture
#' @noRd
mixture2p_truth <- function() {
  readRDS(test_path("fixtures", "mixture2p-truth.rds"))
}
