# Prior standard deviations from a fit's priors

The standard deviation on the link scale that each parameter's prior
implies, for the priors that have one in closed form. Pass the result as
`prior_sd` to
[`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
or
[`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md)
to get the contraction of each posterior against its prior without
fitting the prior alone.

## Usage

``` r
prior_sd_table(x, group = "id", ...)
```

## Arguments

- x:

  A fit (a `brmsfit` or `bmmfit`, read through
  [`brms::prior_summary()`](https://mc-stan.org/rstantools/reference/prior_summary.html)),
  or a prior table with the columns of one: `prior`, `class`, `coef`,
  `group`, and `nlpar` or `dpar`, as
  `bmm::default_prior(formula, data, model)` returns before anything is
  fitted.

- group:

  The grouping factor whose standard deviations are wanted, as in
  [`extract_estimates()`](https://www.gfrischkorn.org/bmmtools/reference/extract_estimates.md).

- ...:

  Not used.

## Value

A tibble with the columns `term`, `level` (`"population"`, `"effect"` or
`"sd"`), `prior` (the prior as brms writes it) and `prior_sd_link`, one
row per estimated parameter, named as
[`extract_estimates()`](https://www.gfrischkorn.org/bmmtools/reference/extract_estimates.md)
names it.

## Details

A coefficient whose own prior is empty takes the prior of its class, as
brms does. The closed forms are, for coefficients (class `b`):
`normal(mu, sigma)` gives `sigma`, `logistic(mu, s)` gives
`s * pi / sqrt(3)`, and `student_t(nu, mu, sigma)` gives
`sigma * sqrt(nu / (nu - 2))` for `nu > 2`. For standard deviations
(class `sd`), which brms bounds at 0, the prior is the half of a
distribution centred on 0, except `exponential(lambda)`, bmm's default
there, which gives `1 / lambda`. Half-normal gives
`sigma * sqrt(1 - 2 / pi)` and half-Student-t, for `nu > 2`, the square
root of `sigma^2 * nu / (nu - 2)` minus its squared mean
`2 * sigma * sqrt(nu) * gamma((nu + 1) / 2) / (sqrt(pi) * gamma(nu / 2) * (nu - 1))`.

Any other prior is `NA` — a flat prior, a prior with bounds on a
coefficient, a Student t with `nu <= 2` (whose variance is infinite), a
half prior not centred on 0, another family — and one message names
those terms. A prior-only fit covers them: pass it, or a
[`prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/prior_check.md),
as `prior_sd` instead. Fixed parameters (`constant()`) have no row, as
they have no estimate.

A Student t prior on a standard deviation with `nu = 3`, brms's default,
has an infinite fourth moment, so the SD of a prior-only fit's draws
converges slowly; the closed form is exact.

## See also

[`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
for contraction, which this feeds.

## Examples

``` r
prior_sd_table(bmm::default_prior(
  bmm::bmf(kappa ~ 1 + (1 | id), thetat ~ 1 + (1 | id)),
  data = data.frame(
    y = c(0.1, -0.2, 0.3, 0.05), id = factor(c(1, 1, 2, 2))
  ),
  model = bmm::mixture2p(resp_error = "y")
))
#> # A tibble: 4 × 4
#>   term   level      prior                prior_sd_link
#>   <chr>  <chr>      <chr>                        <dbl>
#> 1 kappa  population normal(2, 1)                  1   
#> 2 thetat population logistic(0, 1)                1.81
#> 3 kappa  sd         student_t(3, 0, 2.5)          3.34
#> 4 thetat sd         student_t(3, 0, 2.5)          3.34
```
