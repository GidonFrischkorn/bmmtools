# Compare implementations on the same data, in interleaved runs

One row per run, with compile time and sampling time kept apart, ESS per
chain-second for bulk and tail, and the cost of one gradient evaluation.
Runs are interleaved — A, B, A, B — because running every A before every
B would confound the comparison with drift in machine load, and the
first run of each implementation is discarded: it carries the fresh
compile and a cold machine. The discarded run is still returned, marked
`warmup_run`, so that nothing is dropped silently.

## Usage

``` r
benchmark(
  implementations,
  data,
  ...,
  runs = 5L,
  warmup_runs = 1L,
  compile = TRUE,
  on_error = c("record", "stop"),
  .metrics = NULL,
  .compiler = NULL
)
```

## Arguments

- implementations:

  A named list, one element per implementation. Each is either a
  function `(data, ...) -> fit` or a list with a `fit` function and,
  optionally, a `compile` function `(data, ...) -> Stan code` whose
  compile is timed once.

- data:

  A data frame used for every pair, or a list of data frames, one per
  pair. A list is not recycled: it must hold `runs + warmup_runs` sets,
  because a reused data set would make the paired ratios dependent.

- ...:

  Passed to every implementation's fitter (`chains`, `iter`, `cores`,
  `seed`, `backend`, ...). Pass the same arguments to all of them;
  [`summary()`](https://rdrr.io/r/base/summary.html) warns when the
  recorded settings differ.

- runs:

  Recorded pairs per implementation.

- warmup_runs:

  Leading runs per implementation that are recorded and then excluded
  from [`summary()`](https://rdrr.io/r/base/summary.html).

- compile:

  Measure a forced compile once per implementation, for those that
  supply a `compile` hook.

- on_error:

  `"record"` keeps a row whose measurements are `NA`, puts the message
  in `error` and goes on; `"stop"` aborts.

- .metrics:

  The measuring function,
  [`benchmark_metrics()`](https://www.gfrischkorn.org/bmmtools/reference/benchmark_metrics.md)
  by default. Tests inject a stand-in so that no model is fitted.

- .compiler:

  Passed to
  [`benchmark_compile()`](https://www.gfrischkorn.org/bmmtools/reference/benchmark_compile.md).

## Value

A tibble of class `bmmtools_benchmark`, one row per run:

- identity — `implementation`, `pair` (0 for a warmup pair), `run`,
  `warmup_run`, `data_id`, `data_rows`;

- settings as run — `backend`, `chains`, `iter`, `warmup`, `thin`,
  `fit_cores`, `fit_threads`, `seed`;

- time in seconds — `wall_seconds` (the fitting call), `compile_seconds`
  with `compile_source`, `warmup_chain_seconds`,
  `sampling_chain_seconds`, `total_chain_seconds` (each summed over
  chains), `max_chain_seconds` and `overhead_seconds`
  (`wall_seconds - max_chain_seconds`, which contains a per-fit compile
  under a backend that recompiles, and R-side overhead under cmdstanr
  with a cached executable);

- work — `leapfrog_sampling`, `leapfrog_warmup` (`NA` unless warmup
  draws were saved), `us_per_gradient` (the median over chains of chain
  sampling time over chain leapfrog steps), with `us_per_gradient_min`,
  `us_per_gradient_max` and `us_per_gradient_warmup`;

- efficiency — `ess_bulk_min`, `ess_tail_min` and the four rates
  `ess_{bulk,tail}_per_chain_sec_{sampling,total}`, per chain-second,
  over post-warmup or over warmup plus sampling;

- what the speed must be read with — `max_rhat`, `n_divergent`,
  `n_max_treedepth`, `n_variables`, `converged`, `error`;

- every
  [`run_info()`](https://www.gfrischkorn.org/bmmtools/reference/run_info.md)
  column, per run.

The per-chain table of every run is in the attribute `chains`, the
compile table in `compile`, and the call and its arguments in
`settings`.

## Details

No fit is written. The module never calls
[`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md),
so a benchmark fit can never be served to a recovery grid, and no cache
key changes because of it.

## See also

[`benchmark_metrics()`](https://www.gfrischkorn.org/bmmtools/reference/benchmark_metrics.md),
[`benchmark_compile()`](https://www.gfrischkorn.org/bmmtools/reference/benchmark_compile.md),
[`run_info()`](https://www.gfrischkorn.org/bmmtools/reference/run_info.md)

## Examples

``` r
if (FALSE) { # \dontrun{
bench <- benchmark(
  implementations = list(
    current = function(data, ...) brms::brm(current_formula, data, ...),
    rebuilt = function(data, ...) {
      bmm::bmm(
        bmm::bmf(kappa ~ 1, thetat ~ 1), data,
        bmm::mixture2p(resp_error = "y"), ...
      )
    }
  ),
  data = data_sets, runs = 5,
  backend = "cmdstanr", chains = 4, iter = 2000, cores = 4
)
summary(bench)$ratios
} # }
```
