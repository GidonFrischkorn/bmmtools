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
  (the between-subject standard deviations), or `c("population", "sd")`
  for both.

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

- ...:

  Not used.

- group:

  The grouping factor subject-level estimates come from.

## Value

A `bmmtools_recovery` object: a tibble subclass with the columns `term`,
`estimate`, `ci_low`, `ci_high`, `ci_method`, `ci_level`, `rhat`,
`ess_bulk`, `ess_tail`, `true_value`, `bias`, `covered`, `scale`,
`level`, `id`, `converged`, `condition`, `estimator` and `replication`.
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
#> # A tibble: 2 × 23
#>   term   estimator level      scale       n n_replications n_converged    bias
#>   <chr>  <chr>     <chr>      <chr>   <dbl>          <int>       <int>   <dbl>
#> 1 kappa  bayes     population natural     1              1          NA  0.777 
#> 2 thetat bayes     population natural     1              1          NA -0.0192
#>     rmse coverage ci_width     r r_low r_high rank_r   ccc ccc_low ccc_high
#>    <dbl>    <dbl>    <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>   <dbl>    <dbl>
#> 1 0.777         1    8.51     NA    NA     NA     NA    NA      NA       NA
#> 2 0.0192        1    0.233    NA    NA     NA     NA    NA      NA       NA
#>   ccc_accuracy ccc_scale_shift ccc_location_shift calibration_slope truth_sd
#>          <dbl>           <dbl>              <dbl>             <dbl>    <dbl>
#> 1           NA              NA                 NA                NA        0
#> 2           NA              NA                 NA                NA        0
```
