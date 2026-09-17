# Extract a tidy table of parameter estimates from a fit

The bridge between a fitted model and everything bmmtools scores. It
returns the *estimates tibble*: one row per parameter, with the
posterior median, a credible interval and the convergence diagnostics,
under bmm's own parameter names.
[`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
and
[`recover_subjects()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
take either a fit or one of these tibbles, so a scoring pipeline can be
built, tested and run with no fitting package installed.

## Usage

``` r
extract_estimates(fit, ...)

# Default S3 method
extract_estimates(fit, ...)

# S3 method for class 'brmsfit'
extract_estimates(
  fit,
  level = c("population", "subject", "sd", "cor"),
  group = NULL,
  ci_level = 0.95,
  ci_method = "eti",
  drop_constants = TRUE,
  converged = NULL,
  estimator = "bayes",
  ...
)
```

## Arguments

- fit:

  A `brmsfit`, and so also a `bmmfit`.

- ...:

  Not used. Present so the generic can gain arguments later; anything
  passed is an error.

- level:

  Which estimates to return: `"population"`, `"subject"`, `"sd"` (the
  group-level standard deviations), `"cor"` (the group-level
  correlations), or several, in which case they are stacked and the
  `level` column separates them.

- group:

  The grouping factor subject-level, SD and correlation estimates come
  from. `NULL` uses the fit's only grouping factor and errors if there
  is more than one.

- ci_level:

  The interval mass, a single number strictly between 0 and 1.

- ci_method:

  The interval type. Only `"eti"`, the equal-tailed interval, is
  available in this version; the column is carried so that adding others
  later is not a breaking change.

- drop_constants:

  Drop the parameters the model fixed to constants. A constant is
  identified by zero posterior variance, not by a missing rhat, so that
  a chain that broke is never silently dropped as though it had been
  fixed on purpose.

- converged:

  Whether the fit passed the convergence gate. `NULL`, the default,
  computes it as `check_convergence(fit)$pass` with the default
  thresholds; a logical scalar is used as given, so a verdict from
  [`check_convergence()`](https://www.gfrischkorn.org/bmmtools/reference/check_convergence.md)
  with other thresholds can be passed in; `NA` marks it unknown.

- estimator:

  A label for how the estimates were produced, carried into the
  `estimator` column and used by
  [`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
  to keep two estimators of the same parameter apart in
  [`summary()`](https://rdrr.io/r/base/summary.html) and in
  [`plot_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/plot_recovery.md).
  `"bayes"`, the default, is the posterior of a sampled fit. Pass
  another label — `"ml"` for a maximum-likelihood fit, or any name of
  your own — when scoring several estimators against one truth. Without
  it, rows from two estimators would be pooled into a single bias and
  RMSE.

## Value

A tibble with the columns `term`, `estimate`, `ci_low`, `ci_high`,
`ci_method`, `ci_level`, `rhat`, `ess_bulk`, `ess_tail`, `level`, `id`,
`converged` and `estimator`, in that order. `id` is `NA` except on
subject rows; `converged` and `estimator` are the same value on every
row of a fit.

## Details

The first nine columns are the apabayes `parameters` contract, in its
order, so `apabayes::apabayes_tidy(x, type = "parameters")` accepts the
result without renaming.

Subject-level estimates are the **per-draw sum** of the population
intercept and the group-level deviation, summarised afterwards. That
puts them on the same scale as the population rows, so a per-subject
truth can be compared against them directly; taking the sum after
summarising instead would give intervals that are too narrow.

Standard deviations and correlations are returned on the link scale,
with `id` set to `NA`. An SD row carries the parameter's name (`kappa`);
a correlation row carries the two names joined by `__` and sorted in the
C locale (`kappa__thetat`), whichever order brms used. A correlation the
model does not estimate has no row.

**Terms.** A parameter with one coefficient carries its own name. One
with several coefficients and no intercept, as cell-means coding
(`kappa ~ 0 + task`) gives, has a row per coefficient, named
`<parameter>_<coefficient>` (`kappa_task1`) at every level. An intercept
together with other coefficients (population effects or contrasts) is an
error in this version.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- bmm::bmm(bmm::bmf(kappa ~ 1, thetat ~ 1), data, model)
extract_estimates(fit)
extract_estimates(fit, level = "subject", ci_level = 0.89)
} # }
```
