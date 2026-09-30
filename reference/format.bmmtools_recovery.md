# Format and print a recovery object

The printed object leads with the scale it was scored on, because bias,
RMSE and the correlations are not invariant to the link transform and a
table of numbers without that label cannot be read. Where the
correlation metrics are `NA` the reason is printed as a sentence rather
than left as blank cells.

## Usage

``` r
# S3 method for class 'bmmtools_recovery'
format(x, ...)

# S3 method for class 'bmmtools_recovery'
print(x, ...)

# S3 method for class 'bmmtools_recovery_summary'
print(x, ...)
```

## Arguments

- x:

  A `bmmtools_recovery` object.

- ...:

  Not used.

## Value

[`format()`](https://rdrr.io/r/base/format.html) returns a character
vector; [`print()`](https://rdrr.io/r/base/print.html) returns `x`
invisibly.
