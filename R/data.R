#' Example recovery object: a mixture2p grid
#'
#' A precomputed [recovery_grid()] result, shipped so that `summary()`,
#' [plot_recovery()] and dplyr verbs can be tried without compiling a
#' model. The grid crosses 20 and 50 subjects with 30 and 100 trials per
#' subject, with 5 replications per cell, for [bmm::mixture2p()] with
#' population values `kappa = log(8)` and `thetat = qlogis(0.75)` on the
#' link scale and between-subject SDs of 0.3 and 0.5. Every fit used 4
#' chains of 1000 iterations with cmdstanr, and the master seed was 2026.
#' Scores are on the natural scale.
#'
#' @format A `bmmtools_recovery` tibble at both the population and the
#'   subject level, with the columns documented in [recover()] and a
#'   `condition` column naming the grid row. The attribute `cells` holds
#'   one row per cell with its seed, status, runtime and convergence
#'   verdict; `grid` holds the design.
#'
#' @source `data-raw/example-objects.R` in the package's GitHub
#'   repository.
#'
#' @seealso [recovery_grid()], [summary.bmmtools_recovery()]
#'
#' @examples
#' summary(recovery_mixture2p)
#' attr(recovery_mixture2p, "cells")
"recovery_mixture2p"

#' Example prior check: two priors for sdt_yn
#'
#' A precomputed [prior_check()] result, shipped so that `summary()` and
#' [plot_prior_check()] can be tried without compiling a model. Two
#' prior sets are compared for [bmm::sdt_yn()] with a random intercept on
#' `d` and `criterion`: bmm's defaults (`default`) and the defaults with
#' `normal(0, 1)` priors on both between-subject SDs (`narrow_sd`). The
#' data are 20 simulated subjects with 50 trials per stimulus class; 100
#' prior-predictive draws are kept per set, from 2 chains of 1000
#' iterations. Seed 2026.
#'
#' @format A `bmmtools_prior_check` tibble with the columns `prior`,
#'   `response`, `statistic` and `value`, carrying the draws, the prior
#'   summaries of both fits, the data and the model as attributes.
#'
#' @source `data-raw/example-objects.R` in the package's GitHub
#'   repository.
#'
#' @seealso [prior_check()], [summary.bmmtools_prior_check()]
#'
#' @examples
#' summary(prior_check_sdt_yn)
"prior_check_sdt_yn"
