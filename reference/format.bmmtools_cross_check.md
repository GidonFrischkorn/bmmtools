# Format and print a cross-check

The header names the scale and where the reference came from, and says
that the reference is a comparison rather than a truth, because a column
called `bias` invites the other reading.

## Usage

``` r
# S3 method for class 'bmmtools_cross_check'
format(x, ...)

# S3 method for class 'bmmtools_cross_check'
print(x, ...)

# S3 method for class 'bmmtools_cross_check_summary'
format(x, ...)

# S3 method for class 'bmmtools_cross_check_summary'
print(x, ...)
```

## Arguments

- x:

  A `bmmtools_cross_check` object.

- ...:

  Not used.

## Value

[`format()`](https://rdrr.io/r/base/format.html) returns a character
vector; [`print()`](https://rdrr.io/r/base/print.html) returns `x`
invisibly.
