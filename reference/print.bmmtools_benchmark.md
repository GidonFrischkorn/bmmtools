# Print a benchmark or its summary

A header line naming the implementations, the runs and how many were
recorded, then the tibble itself.

## Usage

``` r
# S3 method for class 'bmmtools_benchmark'
print(x, ...)

# S3 method for class 'bmmtools_benchmark_summary'
print(x, ...)
```

## Arguments

- x:

  A
  [`benchmark()`](https://www.gfrischkorn.org/bmmtools/reference/benchmark.md)
  result, or a
  [`summary.bmmtools_benchmark()`](https://www.gfrischkorn.org/bmmtools/reference/summary.bmmtools_benchmark.md)
  one.

- ...:

  Passed to the tibble method.

## Value

`x`, invisibly.
