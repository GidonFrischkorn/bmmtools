# Summarise a cross-check into per-parameter metrics

One row per term and level: the mean and root mean square of
`estimate - reference`, the share of fit intervals covering the
reference, the share of rows whose intervals overlap, and the
correlation and Lin's concordance between the two sides.

## Usage

``` r
# S3 method for class 'bmmtools_cross_check'
summary(object, ...)

# S3 method for class 'bmmtools_cross_check_summary'
summary(object, ...)
```

## Arguments

- object:

  A `bmmtools_cross_check` object from
  [`cross_check()`](https://www.gfrischkorn.org/bmmtools/reference/cross_check.md).

- ...:

  Not used.

## Value

A `bmmtools_cross_check_summary` tibble with the columns `term`,
`level`, `scale`, `n`, `n_converged`, `bias`, `rmse`, `coverage`,
`share_overlap`, `r`, `r_low`, `r_high`, `ccc`, `ccc_low` and
`ccc_high`.

## Details

`bias` here is a difference from a comparison value, not from a known
truth; see
[`cross_check()`](https://www.gfrischkorn.org/bmmtools/reference/cross_check.md).

`r` and `ccc` come from the pairs of `estimate` and `reference` in the
group and are `NA`, never `0`, below three complete pairs or without
spread on either side. At the population level with one fit each term
contributes a single pair, so both are `NA` and the table is the
per-term difference and coverage; at the subject level they are the
comparison across subjects. `share_overlap` is `NA` for a reference
without intervals. `n_converged` counts the rows whose fit passed
[`check_convergence()`](https://www.gfrischkorn.org/bmmtools/reference/check_convergence.md)
and is `NA` when no row carries a verdict.
