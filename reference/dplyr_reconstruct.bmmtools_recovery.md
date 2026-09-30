# Keep the class only while the contract holds

[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
and friends return a recovery object; a `select()` that drops a contract
column returns a plain tibble.

## Usage

``` r
# S3 method for class 'bmmtools_recovery'
dplyr_reconstruct(data, template)
```

## Arguments

- data, template:

  See
  [`dplyr::dplyr_reconstruct()`](https://dplyr.tidyverse.org/reference/dplyr_extending.html).

## Value

`data`, as a recovery object when it still satisfies the contract and as
a plain tibble otherwise.
