# Compare prior sets on the observable scale

One row per response and statistic, one column per prior set. With
exactly two sets a `difference` column is appended, the second minus the
first. Read it as an equivalence table: two priors are equivalent on the
observable scale where the differences are negligible. With three or
more sets there is no unambiguous contrast, so none is added.

## Usage

``` r
# S3 method for class 'bmmtools_prior_check'
summary(object, ...)

# S3 method for class 'bmmtools_prior_check_summary'
summary(object, ...)
```

## Arguments

- object:

  A `bmmtools_prior_check` from
  [`prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/prior_check.md).

- ...:

  Not used.

## Value

A `bmmtools_prior_check_summary` tibble.

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
