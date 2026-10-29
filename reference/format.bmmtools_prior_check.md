# Format and print a prior check

The header names the model and the prior sets compared, because a table
of floor rates without them cannot be read. Where a rate is `NA` the
reason is printed as a sentence rather than left as a blank cell.

## Usage

``` r
# S3 method for class 'bmmtools_prior_check'
format(x, ...)

# S3 method for class 'bmmtools_prior_check'
print(x, ...)

# S3 method for class 'bmmtools_prior_check_summary'
print(x, ...)
```

## Arguments

- x:

  A `bmmtools_prior_check` object.

- ...:

  Not used.

## Value

[`format()`](https://rdrr.io/r/base/format.html) returns a character
vector; [`print()`](https://rdrr.io/r/base/print.html) returns `x`
invisibly.
