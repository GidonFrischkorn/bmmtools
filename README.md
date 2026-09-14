
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bmmtools <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

bmmtools holds the validation half of the model-development workflow of
the [bmm](https://github.com/venpopov/bmm) package.
`bmm::use_model_template()` scaffolds a new cognitive measurement model;
bmmtools checks it. Every check is the same loop, simulate, fit, score,
and only the scorer changes:

| Scorer | truth comes from | scored against | answers |
|----|----|----|----|
| `recover()` | a fixed grid of generating values | posterior vs. truth: correlation with CI, RMSE, bias, CrI coverage | can this model be estimated from data of this size |
| `prior_check()` | the prior | the observable scale: floor and ceiling rates, quantile profile | are these priors sane on the scale a reader understands |
| `sbc()` (*planned*) | the prior | posterior rank of the truth (a thin adapter over the SBC package) | is the implementation correct |
| `cross_check()` (*planned*) | a fixed grid, or real data | a closed-form estimator, another implementation, or published values | does the new model agree with what is already known |

`recover_subjects()` scores person-level posterior parameters, which no
other package does. Scoring dispatches on `brmsfit`, so a plain brms fit
can be scored too; the model-aware generate layer needs bmm.

The package is under construction: `recover()`, `recover_subjects()`,
`prior_check()` and the simulation grid are built, `sbc()` and
`cross_check()` are not yet.

## Who it is for

bmmtools is written for people who develop cognitive measurement models
or study how well they can be estimated: the author of a new bmm model,
a methodologist running a parameter recovery study, a thesis that asks
how many trials a model needs. If you fit bmm models to your own data,
you do not need it; bmm’s own documentation covers that.

## A recovery check in brief

``` r
library(bmmtools)
model <- bmm::mixture2p(resp_error = "y")

sim <- simulate_recovery(
  model,
  pars = c(kappa = log(8), thetat = qlogis(0.75)),
  sds = c(kappa = 0.3, thetat = 0.5),
  n_subjects = 30, n_trials = 60, seed = 1
)
fit <- fit_cached(
  recovery_formula(model), sim$data, model,
  file = "fits/mixture2p", seed = 1
)
check_convergence(fit)

recovery <- recover(fit, sim$truth$population)
summary(recovery)
recover_subjects(fit, sim$truth$subjects) |> plot_recovery()
```

`recovery_grid()` runs the same loop over a grid of subjects, trials and
replications, with one file per cell and a resume, and `prior_check()`
shows what the priors imply for the data before any model is fitted to
them. The articles on the [package
website](https://gidonfrischkorn.github.io/bmmtools/) walk through each
step.

## What bmmtools is not

- Not a simulation-based calibration package. `sbc()` builds a generator
  from a bmm model specification and hands it to
  [SBC](https://hyunjimoon.github.io/SBC/), which owns the ranks, the
  ECDF diagnostics and the plots.
- Not a prior *sensitivity* package.
  [priorsense](https://CRAN.R-project.org/package=priorsense) answers
  how much the posterior moves when the prior is power-scaled.
  `prior_check()` answers a different question: what the prior implies
  on the observable scale before any data are seen.
- Not a plotting layer. [bayesplot](https://mc-stan.org/bayesplot/) owns
  posterior plots. bmmtools adds the two plot families bayesplot has no
  concept of: the recovery scatter (truth against posterior estimate
  with credible intervals, one panel per parameter) and the
  prior-predictive panels on the observable scale.
- Not [gp3bayes](https://CRAN.R-project.org/package=gp3bayes), which
  also offers parameter recovery and SBC vignettes but is written for
  pupillometry and eye-tracking data. bmmtools is written for the
  measurement models bmm fits: mixture, signal detection, evidence
  accumulation, multinomial processing tree and memory measurement
  models.
- Not a reporting package. Every scorer returns a tidy tibble whose
  columns follow the
  [apabayes](https://github.com/GidonFrischkorn/apabayes) contract, so a
  recovery table drops into a manuscript through `apabayes` with no glue
  code. Neither package imports the other.
- No model fitting of its own (that is `bmm()`), and no Bayesian power
  or design analysis.

## Installation

bmmtools is not on CRAN yet.

``` r
# install.packages("pak")
pak::pak("GidonFrischkorn/bmmtools")
```
