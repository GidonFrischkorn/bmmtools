# Subset a recovery object

`dplyr_reconstruct()` is not enough on its own. dplyr calls that generic
only when the result's class differs from the input's, and
[`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html)
on a **tibble subclass** subsets through `[`, which keeps the class; the
generic is then never reached and a column-dropping verb returns
something still labelled a recovery object with its contract broken.
Measured on dplyr 1.2.1: without this method
`dplyr::select(x, "term", "estimate")` stays a `bmmtools_recovery` and
the next [`print()`](https://rdrr.io/r/base/print.html) fails inside a
metric with a length mismatch. Demoting here is what makes the
documented behaviour — a verb that drops a contract column drops the
class — actually true.

## Usage

``` r
# S3 method for class 'bmmtools_recovery'
x[...]
```

## Arguments

- x:

  A `bmmtools_recovery` object.

- ...:

  Passed to the tibble method.

## Value

A recovery object while the contract holds, a plain tibble once it does
not.
