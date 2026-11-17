# Spread beside the centre, and the ratio with an interval

Spread beside the centre, and the ratio with an interval

## Usage

``` r
# S3 method for class 'bmmtools_benchmark'
summary(object, metrics = NULL, reference = NULL, ci_level = 0.95, ...)
```

## Arguments

- object:

  A
  [`benchmark()`](https://www.gfrischkorn.org/bmmtools/reference/benchmark.md)
  result.

- metrics:

  Columns to summarise; by default `us_per_gradient`, the two sampling
  ESS rates, `sampling_chain_seconds`, `wall_seconds` and
  `compile_seconds`.

- reference:

  The implementation every ratio is taken against, the first by default.

- ci_level:

  The interval's level.

- ...:

  Not used.

## Value

An object of class `bmmtools_benchmark_summary`: a list with
`per_implementation` (one row per implementation and metric, with
`n_runs`, `n_measured`, `median`, `mad`, `min`, `max`, `n_not_converged`
and `n_error`) and `ratios` (one row per implementation and metric
against `reference`, with `ratio`, `ci_low`, `ci_high`, `n_pairs` and
`ci_method`).

## Details

`mad` is [`stats::mad()`](https://rdrr.io/r/stats/mad.html), the median
absolute deviation scaled by 1.4826 so that it estimates the standard
deviation of a normal sample. `compile_seconds` is measured **once per
implementation**, so its row reports `n_measured = 1` and its ratio says
`"single measurement"` instead of an interval.

Every ratio is the other implementation **over** the reference, for all
metrics. For a time or a cost, below 1 means the other implementation is
faster; for an ESS rate, where more is better, above 1 does.

Warmup runs are excluded. Runs that did not converge are **kept** and
counted: dropping the slow non-converged runs of one implementation
would bias the ratio in its favour. The ratio is the geometric mean of
the per-pair ratios with a t interval on the log ratios, which uses the
pairing the interleaved design creates; with one pair the bounds are
`NA` and `ci_method` says so, so a bare quotient is never reported.

## See also

[`benchmark()`](https://www.gfrischkorn.org/bmmtools/reference/benchmark.md)

## Examples

``` r
if (FALSE) { # \dontrun{
s <- summary(bench)
s$per_implementation
s$ratios
} # }
```
