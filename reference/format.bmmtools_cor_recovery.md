# Format and print a correlation recovery

The printed object leads with the scale and the estimators, because the
three estimators answer different questions and a correlation table
without that label cannot be read.

## Usage

``` r
# S3 method for class 'bmmtools_cor_recovery'
format(x, ...)

# S3 method for class 'bmmtools_cor_recovery'
print(x, ...)

# S3 method for class 'bmmtools_cor_recovery_summary'
print(x, ...)
```

## Arguments

- x:

  A `bmmtools_cor_recovery` object.

- ...:

  Not used.

## Value

[`format()`](https://rdrr.io/r/base/format.html) returns a character
vector; [`print()`](https://rdrr.io/r/base/print.html) returns `x`
invisibly.
