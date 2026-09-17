# Summarise a correlation recovery into per-pair metrics

One row per condition, pair, estimator and scale, summarised across
replications. Errors, coverage and the correlations are reported against
both truths: the generating correlation (`bias`, `rmse`, `coverage`) and
the correlation the simulated subjects had (`bias_sample`,
`rmse_sample`, `coverage_sample`, `r`, `ccc`).

## Usage

``` r
# S3 method for class 'bmmtools_cor_recovery'
summary(object, ...)

# S3 method for class 'bmmtools_cor_recovery_summary'
summary(object, ...)
```

## Arguments

- object:

  A `bmmtools_cor_recovery` object from
  [`recover_correlations()`](https://www.gfrischkorn.org/bmmtools/reference/recover_correlations.md).

- ...:

  Not used.

## Value

A `bmmtools_cor_recovery_summary` tibble with the columns `term`,
`estimator`, `scale`, `n_replications`, `n_converged`, `true_value`,
`sample_sd`, `mean_estimate`, `bias`, `rmse`, `bias_sample`,
`rmse_sample`, `coverage`, `coverage_sample`, `ci_width`,
`rejection_rate`, `false_positive_rate`, `power`, `r`, `r_low`,
`r_high`, `ccc`, `ccc_low` and `ccc_high`, preceded by `condition` when
the object came from a grid.

## Details

`true_value` is the generating correlation when it is the same in every
replication and `NA` otherwise. `sample_sd` is the standard deviation of
the in-sample correlations across replications.

`rejection_rate` is the share of intervals that exclude zero. It is
reported as `false_positive_rate` when every `true_value` is 0 and as
`power` when none is. A `true_value` that is `NA`, which is what a
nonzero correlation becomes on the natural scale, counts as not 0, so
`power` is defined there. With a mix of zero and nonzero values both are
`NA`.

`r` (with a Fisher-z interval) and `ccc` (Lin's concordance, see
[`recovery_ccc()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_ccc.md))
compare the estimates with `sample_value` across replications, one pair
per replication. They are `NA` with fewer than three replications or
when either side has no spread.
