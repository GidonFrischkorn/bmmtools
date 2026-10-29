# The numbers of a validation, one per row

**\[experimental\]**

## Usage

``` r
card_numbers(result)
```

## Arguments

- result:

  A `bmmtools_validation`.

## Value

A tibble with the columns `stage`, `condition` (the design as
`"<subjects> x <trials>"`, the prior set, `"all"`,
`"<subjects> subjects"` for the data requirements read-off, or `NA`),
`level` (the estimate level or the correlation estimator, or `NA`),
`term`, `metric`, `value` and `mcse`, the Monte Carlo standard error
where one follows from the stored results: binomial for a rate, at the
rate `(x + 1) / (n + 2)` so that a rate of 0 or 1 keeps an error above
0, `sqrt((rmse^2 - bias^2) / (n - 1))` for a bias, and the summary's own
for a mean z-score and its SD. `NA` elsewhere.

## Details

Every number a validation card prints, in one long tibble, so a card can
be checked against it and two runs compared with
[`compare_validations()`](https://www.gfrischkorn.org/bmmtools/reference/compare_validations.md).

## See also

[`run_validation()`](https://www.gfrischkorn.org/bmmtools/reference/run_validation.md),
which stores it as `numbers`.
