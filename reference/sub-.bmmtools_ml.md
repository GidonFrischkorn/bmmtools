# Subset an ML estimates object

As for \[`[.bmmtools_recovery`\]:
[`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html)
on a tibble subclass subsets through `[` and never reaches
`dplyr_reconstruct()`, so this method is what drops the class once a
contract column is gone.

## Usage

``` r
# S3 method for class 'bmmtools_ml'
x[...]
```

## Arguments

- x:

  A `bmmtools_ml` object.

- ...:

  Passed to the tibble method.

## Value

An ML object while the contract holds, a plain tibble once it does not.
