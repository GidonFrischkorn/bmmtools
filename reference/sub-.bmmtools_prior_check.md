# Subset a prior-check object

As for \[`[.bmmtools_recovery`\]:
[`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html)
on a tibble subclass subsets through `[` and never reaches
`dplyr_reconstruct()`, so without this method a verb that drops a
contract column would return a broken object still wearing the class.

## Usage

``` r
# S3 method for class 'bmmtools_prior_check'
x[...]
```

## Arguments

- x:

  A `bmmtools_prior_check` object.

- ...:

  Passed to the tibble method.

## Value

A prior-check object while the contract holds, a plain tibble once it
does not.
