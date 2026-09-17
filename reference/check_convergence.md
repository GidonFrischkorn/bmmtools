# Check whether a fit passes the convergence gate

One row of worst-case diagnostics and a verdict. The thresholds are
arguments with the defaults the validation scripts in bmm used: rhat at
most 1.05, bulk ESS at least 400, at most ten divergent transitions. A
fit that fails is meant to be **scored but flagged**, never dropped:
[`extract_estimates()`](https://www.gfrischkorn.org/bmmtools/reference/extract_estimates.md)
carries the verdict as its `converged` column and
[`summary()`](https://rdrr.io/r/base/summary.html) of a recovery object
counts it as `n_converged`.

## Usage

``` r
check_convergence(fit, ...)

# Default S3 method
check_convergence(fit, ...)

# S3 method for class 'brmsfit'
check_convergence(
  fit,
  rhat_max = 1.05,
  ess_bulk_min = 400,
  ess_tail_min = NULL,
  divergent_max = 10,
  treedepth_max = NULL,
  variables = NULL,
  ...
)
```

## Arguments

- fit:

  A `brmsfit`, and so also a `bmmfit`.

- ...:

  Passed to methods.

- rhat_max:

  Largest acceptable rhat over the variables considered.

- ess_bulk_min:

  Smallest acceptable bulk effective sample size.

- ess_tail_min:

  Smallest acceptable tail effective sample size. `NULL`, the default,
  reports it without gating on it.

- divergent_max:

  Largest acceptable number of divergent transitions after warmup,
  summed over chains.

- treedepth_max:

  The sampler's maximum tree depth, used to count the post-warmup
  iterations that hit it. `NULL` reads it from the fit and falls back to
  Stan's default of 10. Tree-depth hits are reported and never gated on:
  hitting the maximum is an efficiency warning, not a validity failure.

- variables:

  Names of the draws variables to consider. `NULL` means every variable
  [`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html)
  returns. Variables with a missing rhat, which is what a parameter
  fixed to a constant has, are dropped rather than failed on.

## Value

A one-row tibble with `max_rhat`, `min_ess_bulk`, `min_ess_tail`,
`n_divergent`, `n_max_treedepth`, `n_variables`, `pass` and `failed`,
the last a comma-separated list of the criteria that failed (`""` when
none did). `pass` is `NA`, not `TRUE`, when no variable had a finite
rhat. The thresholds used are stored in the attribute `thresholds`.

## Details

The numbers come from
[`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html)
and
[`brms::nuts_params()`](https://mc-stan.org/bayesplot/reference/bayesplot-extractors.html).
A fit without sampler diagnostics (a variational fit, or one whose
sampler parameters were not saved) leaves `n_divergent` and
`n_max_treedepth` as `NA` and does not fail the divergence criterion, so
a reader can see it was not checked.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- bmm::bmm(bmm::bmf(kappa ~ 1, thetat ~ 1), data, model)
check_convergence(fit)
check_convergence(fit, rhat_max = 1.01, ess_tail_min = 400)
} # }
```
