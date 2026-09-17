# Example prior check: two priors for sdt_yn

A precomputed
[`prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/prior_check.md)
result, shipped so that
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`plot_prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/plot_prior_check.md)
can be tried without compiling a model. Two prior sets are compared for
`bmm::sdt_yn()` with a random intercept on `d` and `criterion`: bmm's
defaults (`default`) and the defaults with `normal(0, 1)` priors on both
between-subject SDs (`narrow_sd`). The data are 20 simulated subjects
with 50 trials per stimulus class; 100 prior-predictive draws are kept
per set, from 2 chains of 1000 iterations. Seed 2026.

## Usage

``` r
prior_check_sdt_yn
```

## Format

A `bmmtools_prior_check` tibble with the columns `prior`, `response`,
`statistic` and `value`, carrying the draws, the prior summaries of both
fits, the data and the model as attributes.

## Source

`data-raw/example-objects.R` in the package's GitHub repository.

## See also

[`prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/prior_check.md),
[`summary.bmmtools_prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/summary.bmmtools_prior_check.md)

## Examples

``` r
summary(prior_check_sdt_yn)
#> # A tibble: 5 × 5
#>   response statistic    default narrow_sd difference
#>   <chr>    <chr>          <dbl>     <dbl>      <dbl>
#> 1 hits     floor_rate     0.392     0.152     -0.241
#> 2 hits     ceiling_rate   0.420     0.146     -0.274
#> 3 hits     q50           30        23         -7    
#> 4 hits     q90           50        50          0    
#> 5 hits     q95           50        50          0    
```
