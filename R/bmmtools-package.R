#' bmmtools: validate cognitive measurement models fitted with bmm
#'
#' bmmtools checks a bmm model before it is used on real data: that it
#' can be estimated from data of a given size, and that its priors imply
#' plausible data. Every check is the same loop, simulate, fit, score.
#'
#' @section Simulate:
#' [simulate_recovery()] draws a data set and the truth that produced it
#' from a bmm model and population values on the link scale.
#' [recovery_formula()] gives every free parameter a random intercept.
#'
#' @section Fit:
#' [fit_cached()] fits once and reuses the saved fit while nothing that
#' determines it has changed. [check_convergence()] is the convergence
#' gate. [recovery_grid()] runs simulate, fit and score over a design
#' grid with replications, one file per cell.
#'
#' @section Score:
#' [recover()] and [recover_subjects()] score population-level and
#' person-level estimates against the truth; `summary()` of the result
#' gives bias, RMSE, coverage, interval width and the correlations.
#' [prior_check()] summarises prior-predictive draws on the observable
#' scale. [plot_recovery()] and [plot_prior_check()] draw both.
#'
#' The articles on the package website walk through each step:
#' <https://www.gfrischkorn.org/bmmtools/>.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom rlang .data %||%
## usethis namespace: end
NULL
