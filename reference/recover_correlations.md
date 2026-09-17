# Score recovered correlations against the generating ones

The correlation counterpart of
[`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md):
one row per replication, pair and estimator, with the estimate, its
interval, the generating correlation `true_value`, the correlation the
simulated subjects actually had `sample_value`, the errors against both
and whether the interval excludes zero.

## Usage

``` r
recover_correlations(
  fits,
  truth,
  estimator = c("model", "draws", "point"),
  scale = c("link", "natural"),
  links = NULL,
  pairs = NULL,
  group = "id",
  ci_level = 0.95,
  ...
)
```

## Arguments

- fits:

  A fit (a `brmsfit`, or any object with an
  [`extract_subject_draws()`](https://www.gfrischkorn.org/bmmtools/reference/extract_subject_draws.md)
  method), a list of fits with one element per replication, or a tibble
  from
  [`extract_correlations()`](https://www.gfrischkorn.org/bmmtools/reference/extract_correlations.md)
  that may carry `replication` and `condition` columns. With a
  simulation set as `truth`, a set of fits from
  [`fit_components()`](https://www.gfrischkorn.org/bmmtools/reference/fit_components.md)
  (or a list of fits named by component), or an unnamed list of such
  sets, one per replication.

- truth:

  A `bmmtools_simulation`; a `bmmtools_simulation_set` from
  [`simulate_components()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_components.md);
  a list of either, parallel to `fits`; or a list with the tibbles `cor`
  (`term`, `true_value`), `subjects` (`id`, `term`, `true_value`) and
  optionally `covariates` (`id`, `term`, `true_value`), all on the link
  scale and each of them optionally carrying `replication` or
  `condition`. Without a `replication` column, the same truth applies to
  every replication.

- estimator:

  One or more of `"model"`, `"draws"` and `"point"`.

- scale:

  `"link"` or `"natural"`. When `fits` is a tibble and `scale` is not
  given, the tibble's own scale is used.

- links:

  A named character vector mapping a term to a link name. `NULL` reads
  the link table of a `bmmfit`, then of the simulation's model; without
  one, the natural scale falls back to the link scale with a message.

- pairs:

  `NULL` for every pair in `truth`, or a character vector of pair terms;
  restricts both the estimates and the truth.

- group:

  The grouping factor, `"id"` by default.

- ci_level:

  The interval mass, passed to
  [`extract_correlations()`](https://www.gfrischkorn.org/bmmtools/reference/extract_correlations.md).

- ...:

  Passed to
  [`extract_correlations()`](https://www.gfrischkorn.org/bmmtools/reference/extract_correlations.md),
  for example `converged`.

## Value

A `bmmtools_cor_recovery` object: a tibble subclass with the columns
`term`, `var1`, `var2`, `estimator`, `estimate`, `ci_low`, `ci_high`,
`ci_method`, `ci_level`, `rhat`, `ess_bulk`, `ess_tail`, `true_value`,
`sample_value`, `bias`, `bias_sample`, `covered`, `covered_sample`,
`excludes_zero`, `scale`, `n`, `converged`, `condition` and
`replication`. Call
[summary()](https://www.gfrischkorn.org/bmmtools/reference/summary.bmmtools_cor_recovery.md)
on it for the per-pair metrics and
[`plot_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/plot_recovery.md)
to draw it.

## Details

Two truths are carried because they answer different questions. A finite
sample of subjects does not have the generating correlation;
`sample_value` is what an estimator could at best recover from these
subjects, and `true_value` is what a reader wants to learn about.

`bias` is `estimate - true_value` and `bias_sample` is
`estimate - sample_value`; `covered` and `covered_sample` use the closed
interval; `excludes_zero` is `ci_high < 0 | ci_low > 0`.

**On the natural scale** `sample_value` is recomputed from the subject
values after
[`inverse_link()`](https://www.gfrischkorn.org/bmmtools/reference/inverse_link.md),
covariates untransformed. A generating correlation on the link scale has
no natural-scale counterpart that one transform gives, so `true_value`
is `0` where the link-scale value is 0 (independence survives a monotone
transform) and `NA` otherwise. The `model` estimator has no
natural-scale rows; with fits as input it is dropped with a message.

A truth pair no estimate matches produces a warning and is dropped;
every pair unmatched is an error.

**Separate fits.** A correlation between parameters of two separately
fitted models is attenuated by the reliabilities of both subject
estimates, so `bias` against `true_value` is expected even when every
fit is sound; see the section "Separate fits" of
[`extract_correlations()`](https://www.gfrischkorn.org/bmmtools/reference/extract_correlations.md).

## Examples

``` r
extracted <- tibble::tibble(
  term = "kappa__thetat", var1 = "kappa", var2 = "thetat",
  estimator = "draws", estimate = c(0.42, 0.18, 0.55),
  ci_low = c(0.1, -0.2, 0.2), ci_high = c(0.7, 0.5, 0.8),
  ci_method = "eti", ci_level = 0.95, rhat = 1, ess_bulk = 800,
  ess_tail = 800, scale = "link", n = 5L, replication = 1:3
)
subjects <- tibble::tibble(
  id = rep(as.character(1:5), 6),
  term = rep(rep(c("kappa", "thetat"), each = 5), 3),
  true_value = c(
    1.2, 0.4, 2.0, 1.1, 0.7, 0.3, -0.5, 0.9, 0.1, -0.2,
    0.8, 1.5, 1.1, 0.2, 1.9, 0.0, 0.6, 0.4, -0.3, 0.8,
    1.4, 0.9, 0.3, 1.8, 1.0, 0.5, 0.2, -0.4, 0.9, 0.1
  ),
  replication = rep(1:3, each = 10)
)
truth <- list(
  cor = tibble::tibble(term = "kappa__thetat", true_value = 0.5),
  subjects = subjects
)
recovery <- recover_correlations(extracted, truth)
summary(recovery)
#> # A tibble: 1 × 24
#>   term          estimator scale n_replications n_converged true_value sample_sd
#>   <chr>         <chr>     <chr>          <int>       <int>      <dbl>     <dbl>
#> 1 kappa__thetat draws     link               3          NA        0.5   0.00513
#>   mean_estimate   bias  rmse bias_sample rmse_sample coverage coverage_sample
#>           <dbl>  <dbl> <dbl>       <dbl>       <dbl>    <dbl>           <dbl>
#> 1         0.383 -0.117 0.193      -0.607       0.626        1               0
#>   ci_width rejection_rate false_positive_rate power     r r_low r_high     ccc
#>      <dbl>          <dbl>               <dbl> <dbl> <dbl> <dbl>  <dbl>   <dbl>
#> 1    0.633          0.667                  NA 0.667 0.511    NA     NA 0.00167
#>   ccc_low ccc_high
#>     <dbl>    <dbl>
#> 1      NA       NA
```
