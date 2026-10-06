# Score parameter recovery against known generating values

`recover()` answers the question a model developer has to answer before
anything else: when the data are generated from known parameters, does
fitting the model give those parameters back? It takes fits and truth
and returns one row per fit and parameter, with the estimate, its
interval, the generating value, the signed error and whether the
interval covered.

## Usage

``` r
recover(
  fits,
  truth,
  level = "population",
  scale = c("natural", "link"),
  links = NULL,
  ci_level = 0.95,
  drop_constants = TRUE,
  prior_sd = NULL,
  ...
)

recover_subjects(
  fits,
  truth,
  group = "id",
  scale = c("natural", "link"),
  links = NULL,
  ci_level = 0.95,
  drop_constants = TRUE,
  ...
)
```

## Arguments

- fits:

  A `brmsfit` (so also a `bmmfit`), a list of them with one element per
  replication, or an estimates tibble from
  [`extract_estimates()`](https://www.gfrischkorn.org/bmmtools/reference/extract_estimates.md).
  The tibble is accepted so that a scoring pipeline runs with no fitting
  package installed.

- truth:

  A data frame of generating values on the **link scale**, under bmm's
  own parameter names. `recover()` needs `term` and `true_value`;
  `recover_subjects()` also needs `id`. A `replication` column is used
  as a join key when present, and when it is absent the same generating
  values are scored against every replication. When `recover()` scores
  several levels, a named list with a data frame for each, such as the
  `truth` of a `bmmtools_simulation`.

- level:

  The levels `recover()` scores: `"population"` (the default), `"sd"`
  (the between-subject standard deviations), `"effect"` (the contrasts
  of a `coding = "contrast"` design), or several of them. `"sd"` and
  `"effect"` are always scored on the link scale, whatever `scale` says.

- scale:

  `"natural"` scores on the scale a reader interprets, `"link"` on the
  scale the model was estimated on. The choice matters: coverage is
  invariant to a monotone link but bias, RMSE and the correlations are
  not.

- links:

  A named character vector mapping a term to one of the link names
  [`inverse_link()`](https://www.gfrischkorn.org/bmmtools/reference/inverse_link.md)
  understands. `NULL` reads the link table from a `bmmfit`; with no
  table available, scoring falls back to the link scale and says so. A
  term with no entry takes the link of the longest entry it starts with
  followed by `_` (`kappa_task1` takes the link of `kappa`), and is
  treated as `"identity"` when there is none.

- ci_level:

  The interval mass, passed to
  [`extract_estimates()`](https://www.gfrischkorn.org/bmmtools/reference/extract_estimates.md)
  when `fits` is a fit.

- drop_constants:

  Passed to
  [`extract_estimates()`](https://www.gfrischkorn.org/bmmtools/reference/extract_estimates.md).
  Parameters the model fixed are dropped by default: scoring them would
  report a parameter as perfectly recovered that was never estimated.

- prior_sd:

  `NULL` (the default), or the prior SD on the link scale that each
  posterior SD is compared with for `contraction`. One of: a table with
  the columns `term`, `level` and `prior_sd_link` (and optionally
  `condition`), such as
  [`prior_sd_table()`](https://www.gfrischkorn.org/bmmtools/reference/prior_sd_table.md)
  returns; a fit of the same model with `sample_prior = "only"`, whose
  posterior SDs are the prior SDs, read with
  [`extract_estimates()`](https://www.gfrischkorn.org/bmmtools/reference/extract_estimates.md)
  so its terms match; or a
  [`prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/prior_check.md)
  with one prior set, whose cached fit is read. A fit that sampled the
  data is refused. Only `recover()` takes it: subject rows have no
  contraction (see `contraction` below).

- ...:

  Not used.

- group:

  The grouping factor subject-level estimates come from.

## Value

A `bmmtools_recovery` object: a tibble subclass with the columns `term`,
`estimate`, `ci_low`, `ci_high`, `ci_method`, `ci_level`, `rhat`,
`ess_bulk`, `ess_tail`, `true_value`, `bias`, `covered`, `scale`,
`level`, `id`, `converged`, `condition`, `estimator`, `ci_low_50`,
`ci_high_50`, `covered_50`, `post_mean_link`, `post_sd_link`, `z`,
`contraction` and `replication`.

`post_mean_link` and `post_sd_link` are the mean and SD of the posterior
draws **on the link scale, whatever `scale` is**: neither is invariant
under a nonlinear link, so neither is transformed with the estimate and
its interval. `estimate` stays the posterior median. For a
maximum-likelihood row they are the point estimate and its standard
error; for a hand-built tibble without them, `NA`.

`z` is the posterior z-score, `(post_mean_link - truth) / post_sd_link`
with the truth on the link scale. Across replications whose generating
values are drawn from the prior, a calibrated posterior gives z mean 0
and SD 1; with generating values held fixed, as a recovery grid holds
them, z also carries the prior's pull toward its centre, so a calibrated
posterior can show a nonzero mean and an SD below 1 (see
[`summary.bmmtools_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/summary.bmmtools_recovery.md)).
`contraction` is `1 - post_sd_link^2 / prior_sd_link^2`, the share of
the prior variance the data removed: near 1 when the data decide the
estimate, near 0 when the prior does, and below 0 when the posterior is
wider than the prior. Both stay on the link scale. `z` is `NA` where the
posterior SD is 0 or missing. `contraction` is `NA` without `prior_sd`,
where it has no row for the term, for subject rows (a subject's
posterior narrows through pooling as well as through its own data, so
the ratio would not measure what the data did), and for
maximum-likelihood rows, which have no prior.

`converged` is the verdict of
[`check_convergence()`](https://www.gfrischkorn.org/bmmtools/reference/check_convergence.md)
with its default thresholds when `fits` are fit objects; to gate with
other thresholds, call
[`extract_estimates()`](https://www.gfrischkorn.org/bmmtools/reference/extract_estimates.md)
with `converged =` first and pass the tibble. Call
[`summary()`](https://rdrr.io/r/base/summary.html) on it for the
per-parameter metrics.

`estimator` names how each row was produced and defaults to `"bayes"`.
[`summary()`](https://rdrr.io/r/base/summary.html) groups by it, so two
estimators of the same parameter — a hierarchical posterior and a
subject-wise maximum-likelihood fit, say — can be bound together and
scored against one truth without being pooled into a single bias and
RMSE. Set it with `extract_estimates(estimator = )` or as a column on a
hand-built tibble, then
[`dplyr::bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html)
the two and pass the result here. Note that this `estimator` is
unrelated to the one in
[`extract_correlations()`](https://www.gfrischkorn.org/bmmtools/reference/extract_correlations.md),
which names how a *correlation* was read.

## Details

`recover_subjects()` is the person-level variant. It asks whether the
model orders *people* correctly, which is the question a reliability
analysis needs and which population-level recovery cannot answer.

Standard deviations are always scored on the link scale, the scale the
model estimates them on. Under `scale = "natural"` a message says so,
their `scale` column reads `"link"`, and the other rows and the object's
`scale` attribute stay on the natural scale.

So is a term whose link is `softmax`, such as `mixture3p`'s `thetat` and
`thetant`: they are log weights against a guessing weight of 0, and a
trial's probabilities are the softmax of the group, so no inverse of one
term gives its natural value. Under `scale = "natural"` a message names
them and their `scale` column reads `"link"`.

A term in `truth` that no fit estimated produces a warning listing what
was available, and is dropped. A term a fit estimated that `truth` does
not mention is dropped silently, because fits routinely estimate more
than a simulation grid varies. Every term unmatched is an error rather
than an empty result.

## Examples

``` r
estimates <- tibble::tibble(
  term = c("kappa", "thetat"),
  estimate = c(2.1, 1.0), ci_low = c(1.6, 0.4), ci_high = c(2.6, 1.6),
  ci_method = "eti", ci_level = 0.95,
  rhat = 1, ess_bulk = 900, ess_tail = 900,
  level = "population", id = NA_character_
)
truth <- tibble::tibble(term = c("kappa", "thetat"), true_value = c(2, 1.1))

recovery <- recover(
  estimates, truth,
  links = c(kappa = "log", thetat = "logit")
)
summary(recovery)
#> # A tibble: 2 × 15
#>   term   level      scale       n    bias   rmse coverage coverage_50 z_mean
#>   <chr>  <chr>      <chr>   <dbl>   <dbl>  <dbl>    <dbl>       <dbl>  <dbl>
#> 1 kappa  population natural     1  0.777  0.777         1          NA     NA
#> 2 thetat population natural     1 -0.0192 0.0192        1          NA     NA
#>    z_sd contraction     r   ccc calibration_slope calibration_intercept
#>   <dbl>       <dbl> <dbl> <dbl>             <dbl>                 <dbl>
#> 1    NA          NA    NA    NA                NA                    NA
#> 2    NA          NA    NA    NA                NA                    NA
#> # 19 more columns; `tibble::as_tibble()` prints them all.
```
