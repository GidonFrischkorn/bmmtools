# Package index

## About bmmtools

- [`bmmtools`](https://www.gfrischkorn.org/bmmtools/reference/bmmtools-package.md)
  [`bmmtools-package`](https://www.gfrischkorn.org/bmmtools/reference/bmmtools-package.md)
  : bmmtools: validate cognitive measurement models fitted with bmm

## Simulate

Data sets and the truth that produced them.

- [`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md)
  : Simulate data from a bmm model with known parameters
- [`cors_from_factors()`](https://www.gfrischkorn.org/bmmtools/reference/cors_from_factors.md)
  : Build a correlation matrix from factor loadings
- [`recovery_formula()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_formula.md)
  : The default recovery formula: every free parameter gets a random
  intercept
- [`recovery_component()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_component.md)
  : Describe one model of a multi-model simulation
- [`simulate_components()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_components.md)
  : Simulate several models for the same simulated people

## Fit and run a study

Cached fitting, the convergence gate and the design grid.

- [`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md)
  : Fit a bmm model once and reuse it while nothing has changed
- [`fit_components()`](https://www.gfrischkorn.org/bmmtools/reference/fit_components.md)
  : Fit every component of a simulation set, each on its own
- [`fit_ml()`](https://www.gfrischkorn.org/bmmtools/reference/fit_ml.md)
  : Subject-wise maximum likelihood estimates of a bmm model
- [`check_convergence()`](https://www.gfrischkorn.org/bmmtools/reference/check_convergence.md)
  : Check whether a fit passes the convergence gate
- [`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md)
  : Run a parameter-recovery grid

## Score recovery

Estimates against generating values, per fit and per parameter.

- [`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
  [`recover_subjects()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
  : Score parameter recovery against known generating values
- [`summary(`*`<bmmtools_recovery>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/summary.bmmtools_recovery.md)
  [`summary(`*`<bmmtools_recovery_summary>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/summary.bmmtools_recovery.md)
  : Summarise a recovery object into per-parameter metrics
- [`recovery_ccc()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_ccc.md)
  : Concordance between recovered and generating values
- [`extract_estimates()`](https://www.gfrischkorn.org/bmmtools/reference/extract_estimates.md)
  : Extract a tidy table of parameter estimates from a fit
- [`recover_correlations()`](https://www.gfrischkorn.org/bmmtools/reference/recover_correlations.md)
  : Score recovered correlations against the generating ones
- [`summary(`*`<bmmtools_cor_recovery>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/summary.bmmtools_cor_recovery.md)
  [`summary(`*`<bmmtools_cor_recovery_summary>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/summary.bmmtools_cor_recovery.md)
  : Summarise a correlation recovery into per-pair metrics
- [`extract_correlations()`](https://www.gfrischkorn.org/bmmtools/reference/extract_correlations.md)
  : Extract between-subject correlations from a fit
- [`extract_subject_draws()`](https://www.gfrischkorn.org/bmmtools/reference/extract_subject_draws.md)
  : Extract subject-level posterior draws as an array
- [`subject_table()`](https://www.gfrischkorn.org/bmmtools/reference/subject_table.md)
  : True and estimated subject values, one row per subject
- [`inverse_link()`](https://www.gfrischkorn.org/bmmtools/reference/inverse_link.md)
  : Transform values from the link scale to the natural scale

## Prior predictive checks

What the priors imply on the observable scale.

- [`prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/prior_check.md)
  : Check what the priors say the data should look like
- [`summary(`*`<bmmtools_prior_check>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/summary.bmmtools_prior_check.md)
  [`summary(`*`<bmmtools_prior_check_summary>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/summary.bmmtools_prior_check.md)
  : Compare prior sets on the observable scale

## Correctness

Estimates against what is already known.

- [`cross_check()`](https://www.gfrischkorn.org/bmmtools/reference/cross_check.md)
  : Compare a fit's estimates with a reference
- [`summary(`*`<bmmtools_cross_check>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/summary.bmmtools_cross_check.md)
  [`summary(`*`<bmmtools_cross_check_summary>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/summary.bmmtools_cross_check.md)
  : Summarise a cross-check into per-parameter metrics
- [`sbc()`](https://www.gfrischkorn.org/bmmtools/reference/sbc.md) :
  Check that a model's implementation is calibrated

## Plot

- [`plot_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/plot_recovery.md)
  : Plot recovered estimates against their generating values
- [`plot_prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/plot_prior_check.md)
  : Plot what the priors say the data should look like

## Example objects

Precomputed results, so examples run without Stan.

- [`recovery_mixture2p`](https://www.gfrischkorn.org/bmmtools/reference/recovery_mixture2p.md)
  : Example recovery object: a mixture2p grid
- [`prior_check_sdt_yn`](https://www.gfrischkorn.org/bmmtools/reference/prior_check_sdt_yn.md)
  : Example prior check: two priors for sdt_yn

## Object methods

Printing, subsetting and dplyr support for the result classes.

- [`format(`*`<bmmtools_cor_recovery>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/format.bmmtools_cor_recovery.md)
  [`print(`*`<bmmtools_cor_recovery>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/format.bmmtools_cor_recovery.md)
  [`print(`*`<bmmtools_cor_recovery_summary>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/format.bmmtools_cor_recovery.md)
  : Format and print a correlation recovery
- [`format(`*`<bmmtools_cross_check>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/format.bmmtools_cross_check.md)
  [`print(`*`<bmmtools_cross_check>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/format.bmmtools_cross_check.md)
  [`format(`*`<bmmtools_cross_check_summary>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/format.bmmtools_cross_check.md)
  [`print(`*`<bmmtools_cross_check_summary>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/format.bmmtools_cross_check.md)
  : Format and print a cross-check
- [`format(`*`<bmmtools_ml>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/format.bmmtools_ml.md)
  [`print(`*`<bmmtools_ml>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/format.bmmtools_ml.md)
  : Format and print subject-wise ML estimates
- [`format(`*`<bmmtools_prior_check>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/format.bmmtools_prior_check.md)
  [`print(`*`<bmmtools_prior_check>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/format.bmmtools_prior_check.md)
  [`print(`*`<bmmtools_prior_check_summary>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/format.bmmtools_prior_check.md)
  : Format and print a prior check
- [`format(`*`<bmmtools_recovery>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/format.bmmtools_recovery.md)
  [`print(`*`<bmmtools_recovery>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/format.bmmtools_recovery.md)
  [`print(`*`<bmmtools_recovery_summary>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/format.bmmtools_recovery.md)
  : Format and print a recovery object
- [`` `[`( ``*`<bmmtools_recovery>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/sub-.bmmtools_recovery.md)
  : Subset a recovery object
- [`` `[`( ``*`<bmmtools_cor_recovery>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/sub-.bmmtools_cor_recovery.md)
  : Subset a correlation-recovery object
- [`` `[`( ``*`<bmmtools_cross_check>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/sub-.bmmtools_cross_check.md)
  : Subset a cross-check object
- [`` `[`( ``*`<bmmtools_prior_check>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/sub-.bmmtools_prior_check.md)
  : Subset a prior-check object
- [`` `[`( ``*`<bmmtools_ml>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/sub-.bmmtools_ml.md)
  : Subset an ML estimates object
- [`dplyr_reconstruct(`*`<bmmtools_recovery>`*`)`](https://www.gfrischkorn.org/bmmtools/reference/dplyr_reconstruct.bmmtools_recovery.md)
  : Keep the class only while the contract holds
