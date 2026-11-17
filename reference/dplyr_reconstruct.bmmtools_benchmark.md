# Keep the class only while the contract holds

Keep the class only while the contract holds

## Usage

``` r
# S3 method for class 'bmmtools_benchmark'
dplyr_reconstruct(data, template)
```

## Arguments

- data, template:

  See
  [`dplyr::dplyr_reconstruct()`](https://dplyr.tidyverse.org/reference/dplyr_extending.html).

## Value

`data`, as a benchmark while it satisfies the contract and as a plain
tibble otherwise.
