# Format and print subject-wise ML estimates

The header says how many cells converged without being asked, because
the difference between noticing three failed subjects and not noticing
them is the difference between a comparison and a wrong number.

## Usage

``` r
# S3 method for class 'bmmtools_ml'
format(x, ...)

# S3 method for class 'bmmtools_ml'
print(x, ...)
```

## Arguments

- x:

  A `bmmtools_ml` object.

- ...:

  Not used.

## Value

[`format()`](https://rdrr.io/r/base/format.html) returns a character
vector; [`print()`](https://rdrr.io/r/base/print.html) returns `x`
invisibly.
