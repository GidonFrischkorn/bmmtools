# Measure one fit: compile, per-chain time, ESS per second, gradient cost

What a fit cost, in units that name their phase and their denominator.
Sampling time comes from the sampler's own clock per chain
([`rstan::get_elapsed_time()`](https://mc-stan.org/rstan/reference/stanfit-class.html),
which works under the cmdstanr backend too), and the work from
`n_leapfrog__`
([`brms::nuts_params()`](https://mc-stan.org/bayesplot/reference/bayesplot-extractors.html)).
Every rate divides by seconds **summed over chains**, never by wall
time: per-chain times are chain-local work, so a wall-time denominator
would make the same fit look cheaper on more cores.

## Usage

``` r
benchmark_metrics(fit, ...)

# Default S3 method
benchmark_metrics(fit, ...)

# S3 method for class 'brmsfit'
benchmark_metrics(
  fit,
  wall_seconds = NA_real_,
  compile_seconds = NA_real_,
  compile_source = NA_character_,
  variables = NULL,
  inc_warmup = NULL,
  rhat_max = 1.05,
  ess_bulk_min = 400,
  ess_tail_min = NULL,
  divergent_max = 10,
  ...
)
```

## Arguments

- fit:

  A fitted model, `bmmfit` or `brmsfit`.

- ...:

  Not used.

- wall_seconds:

  The wall time of the fitting call, which the caller must measure: no
  fit records it.

- compile_seconds:

  Compile time from
  [`benchmark_compile()`](https://www.gfrischkorn.org/bmmtools/reference/benchmark_compile.md);
  no fit records that either.

- compile_source:

  How `compile_seconds` was obtained, for the reader of the row:
  `"forced"`, `"fit"` or `"none"`.

- variables:

  Restrict the rhat and ESS set to these variables.

- inc_warmup:

  Read the warmup phase's leapfrog steps as well. `NULL`, the default,
  uses them when the fit saved its warmup draws, which happens only with
  `save_warmup = TRUE`. A fit that did not save them reports `NA` for
  the warmup work whatever this is set to: asking brms for warmup draws
  it does not have returns the post-warmup ones, which would read as
  zero warmup work.

- rhat_max, ess_bulk_min, ess_tail_min, divergent_max:

  Thresholds for the `converged` column, as in
  [`check_convergence()`](https://www.gfrischkorn.org/bmmtools/reference/check_convergence.md).

## Value

A one-row tibble of class `bmmtools_benchmark`; see
[`benchmark()`](https://www.gfrischkorn.org/bmmtools/reference/benchmark.md)
for the columns. The per-chain table is in the attribute `chains`.

## Details

The minimum ESS is taken over the **estimated** parameters: a parameter
fixed by a constant prior has no posterior variance, no rhat and no ESS,
and `n_variables` records how many were assessed. This is
[`check_convergence()`](https://www.gfrischkorn.org/bmmtools/reference/check_convergence.md)'s
rule, and the convergence columns come from the same code.

## See also

[`benchmark()`](https://www.gfrischkorn.org/bmmtools/reference/benchmark.md),
[`benchmark_compile()`](https://www.gfrischkorn.org/bmmtools/reference/benchmark_compile.md),
[`check_convergence()`](https://www.gfrischkorn.org/bmmtools/reference/check_convergence.md)

## Examples

``` r
if (FALSE) { # \dontrun{
started <- Sys.time()
fit <- bmm::bmm(formula, data, model, backend = "cmdstanr")
row <- benchmark_metrics(
  fit,
  wall_seconds = as.double(difftime(Sys.time(), started, units = "secs"))
)
row$us_per_gradient
} # }
```
