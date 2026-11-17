# Find the smallest design that meets a target

Reads a grid summary design by design and returns, for each group of
rows, the smallest value of `along` (by default the number of trials) at
which a target is met **and stays met at every larger value in the
grid**. A design meets a target when the lower end of the target's Monte
Carlo interval clears it, so the answer accounts for how many
replications the grid ran.

## Usage

``` r
smallest_design(
  x,
  target,
  along = "n_trials",
  by = "n_subjects",
  level = 0.95,
  coverage_check = FALSE
)
```

## Arguments

- x:

  A grid summary that carries the design columns
  ([`summary()`](https://rdrr.io/r/base/summary.html) of a
  [`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md)
  result), or such a recovery object, which is summarised first. Filter
  it first to the rows the question is about, for example
  `level == "subject"`.

- target:

  A named numeric vector of thresholds on summary columns, such as
  `c(r_low = 0.8)` or `c(reliability_model = 0.8, precision_met = 0.8)`.
  Higher is better for every column allowed: `detected`,
  `sign_recovery`, `precision_met`, `reliability_true`,
  `reliability_model`, `r_low` and `ccc_low`.

- along:

  The design column to minimise.

- by:

  Design columns that define separate groups, each with its own answer;
  `NULL` for none. Every grid column that varies besides `along` belongs
  here: each group must have one row per value of `along` (and term,
  level and estimator), or the call is refused.

- level:

  The mass of the Monte Carlo interval whose lower end must clear the
  target.

- coverage_check:

  If `TRUE`, a design counts as meeting a target only when every
  population row at that design also has its coverage within one MCSE of
  nominal, `abs(coverage - ci_level) <= coverage_mcse`. This is a
  calibration condition, not a target. It reads the population rows left
  in `x`, so keep them (their groups drop out of the answer when the
  target column is `NA` there), and a design with none fails the check.

## Value

A `bmmtools_smallest_design` tibble with one row per `by` group, term,
level, estimator and target: the `by` columns, `term`, `level`,
`estimator`, `target` (the column), `threshold`, `smallest` (the
certified value of `along`, or `NA`), `smallest_point` (the same rule
applied to the point estimate, which shows what Monte Carlo error
costs), `status`, and `replications_needed`, the replications a perfect
rate would need to certify the threshold when that is why a design is
undecided (`NA` otherwise). Printing it states each answer as a
sentence.

## Details

**The lower end.** For a rate (`detected`, `sign_recovery`,
`precision_met`) it is the Wilson score bound, because the plain MCSE of
a rate is 0 at 0 and 1 and would certify a target from three successes
in three. Its count is the effective sample size `p (1 - p) / MCSE^2`:
the number of intervals when the MCSE is binomial, and at the subject
level a number corrected for subjects of one fit not being independent,
capped at the number of intervals. At a rate of 0 or 1, or an MCSE of 0,
the number of replications is used, so a perfect rate certifies 0.80
only from 16 replications on. For a reliability it is
`estimate - z * MCSE`. `r_low` and `ccc_low` are already lower ends and
are read as they are; `level` does not change them. `r`, `ccc` and
`ci_width` have no Monte Carlo interval and are refused, naming the
column to use. `coverage` is refused too: "high enough" is the wrong
question for a coverage, which `coverage_check` asks properly.

**Monotone in the grid.** Noise that lets 30 trials pass and 100 fail
does not produce 30. A consequence is that several targets, or several
terms, are met together from the largest of their single answers, which
is how a condition on every parameter is read off.

**No interpolation, no extrapolation.** `status` is `"met"`,
`"met at the grid's smallest design"` (the true minimum may be lower),
`"not reached in this grid"`, or `"undecided: too few replications"`
(the point estimate meets the target, the bound does not). A group whose
target column is `NA` at every design, such as a reliability on
population rows, is left out.

## Examples

``` r
summary(recovery_mixture2p) |>
  dplyr::filter(level == "subject") |>
  smallest_design(c(r_low = 0.8))
#> <bmmtools_smallest_design>
#> n_subjects = 20, kappa (subject), r_low >= 0.8: not reached in this grid.
#> n_subjects = 20, thetat (subject), r_low >= 0.8: met from n_trials = 100.
#> n_subjects = 50, kappa (subject), r_low >= 0.8: not reached in this grid.
#> n_subjects = 50, thetat (subject), r_low >= 0.8: met from n_trials = 100.
```
