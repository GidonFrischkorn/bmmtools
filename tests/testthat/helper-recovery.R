# Hand-built estimates and truth tibbles for the recovery tests.
#
# recover() takes an estimates tibble as well as a fit (ARCHITECTURE.md
# decision 1), so the whole scoring path can be exercised without brms
# and without a saved fit. These helpers build the two inputs.

#' A hand-built estimates tibble
#'
#' Carries the eleven contract columns of `extract_estimates()` plus the
#' `replication` column `recover()` reads when it is given a list.
#'
#' @noRd
fake_estimates <- function(term,
                           estimate,
                           ci_low = estimate - 1,
                           ci_high = estimate + 1,
                           level = "population",
                           id = NA_character_,
                           replication = 1L,
                           ci_level = 0.95,
                           rhat = 1) {
  tibble::tibble(
    term = as.character(term),
    estimate = as.double(estimate),
    ci_low = as.double(ci_low),
    ci_high = as.double(ci_high),
    ci_method = "eti",
    ci_level = as.double(ci_level),
    rhat = as.double(rhat),
    ess_bulk = 1000,
    ess_tail = 1000,
    level = as.character(level),
    id = as.character(id),
    replication = replication
  )
}

#' A hand-built truth tibble
#' @noRd
fake_truth <- function(term, true_value, ...) {
  tibble::tibble(
    term = as.character(term),
    true_value = as.double(true_value),
    ...
  )
}
