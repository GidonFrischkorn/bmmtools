# The verdict of every stage against the plan's criteria

**\[experimental\]**

## Usage

``` r
verdict_table(result)
```

## Arguments

- result:

  A `bmmtools_validation`.

## Value

A tibble with `stage`, `title`, `required`, `status` (the stage's),
`criterion` (the bar in words), `value` (the number judged), `pass`
(`TRUE`, `FALSE`, or `NA` when the stage did not run or the criterion
cannot be evaluated) and `note`.

## Details

One row per criterion, with the bar written out, so a card can print the
criterion next to its verdict rather than a bare pass or fail.

The criteria are the plan's (see
[`validation_plan()`](https://www.gfrischkorn.org/bmmtools/reference/validation_plan.md),
"Criteria"):

- prior predictive: the largest floor or ceiling rate under the
  calibration's prior (or bmm's defaults when the calibration uses them)
  at most `prior_bound_rate_max`; not judged for a response without a
  known range.

- SBC: no ranked variable outside the `sbc_ecdf_prob` simultaneous ECDF
  band, and no failed check (a high Rhat, a low rank ESS or missing
  diagnostics) in any fit. The value is the number of failed checks, so
  a fit that fails two counts twice.

- computation: the share of converged cells at least
  `convergence_rate_min`.

- recovery: the coverage of every population parameter, pooled over the
  designs, not below `coverage_level` by a one-sided binomial test at
  `coverage_alpha` divided by the number of parameters. The coverage of
  each design is reported, not judged.

- identifiability: at the largest design, every population parameter's
  mean contraction above `contraction_min` and its mean z-score within
  `z_mean_max` of 0.

- individual differences: the coverage of every generating correlation
  held to the same pooled test, divided by the number of correlations.

- data requirements: every subject-varying parameter reaches
  `r_low >= read_off_r_low` within the grid at `read_off_n_subjects`
  subjects, the value being the largest read-off over the parameters; a
  parameter the grid gave no person-level correlation is not judged. Its
  coverage leg is the recovery criterion's pooled test, repeated as a
  second row: a per-design coverage bar changed or removed a calibrated
  `ezdm` model's read-off in 24 % to 56 % of runs at 10 replications.

- agreement: reported, not judged.

The ECDF band is simultaneous over the ranks of one variable; it is not
corrected for the number of variables, so with more ranked variables a
calibrated model leaves some band more often.
