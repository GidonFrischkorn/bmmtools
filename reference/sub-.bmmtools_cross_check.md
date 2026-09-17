# Subset a cross-check object

As for \[`[.bmmtools_recovery`\]:
[`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html)
on a tibble subclass subsets through `[` and never reaches
`dplyr_reconstruct()`, so this method is what drops the class once a
contract column is gone.

## Usage

``` r
# S3 method for class 'bmmtools_cross_check'
x[...]
```

## Arguments

- x:

  A `bmmtools_cross_check` object.

- ...:

  Passed to the tibble method.

## Value

A cross-check object while the contract holds, a plain tibble once it
does not.
