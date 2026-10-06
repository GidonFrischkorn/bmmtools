# Subset a benchmark

A verb that drops a contract column returns a plain tibble; the
per-chain and compile tables describe the object as returned, so they
are dropped by a subset rather than left to describe rows that are no
longer there.

## Usage

``` r
# S3 method for class 'bmmtools_benchmark'
x[...]
```

## Arguments

- x:

  A `bmmtools_benchmark` object.

- ...:

  Passed to the tibble method.

## Value

A benchmark while the contract holds, a plain tibble once it does not.
