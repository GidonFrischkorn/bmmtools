# Compare the numbers of two validations

**\[experimental\]**

## Usage

``` r
compare_validations(old, new)
```

## Arguments

- old, new:

  Each a `bmmtools_validation`, a path to its `result.rds` or its
  directory, or a
  [`card_numbers()`](https://www.gfrischkorn.org/bmmtools/reference/card_numbers.md)
  tibble.

## Value

A tibble with the key columns of
[`card_numbers()`](https://www.gfrischkorn.org/bmmtools/reference/card_numbers.md)
(`stage`, `condition`, `level`, `term`, `metric`), then `old`, `new`,
`diff` (`new - old`), `mcse_old`, `mcse_new` and `moved`: where both
runs have an MCSE, `TRUE` when `|diff|` exceeds
`2 * sqrt(mcse_old^2 + mcse_new^2)`; otherwise when it exceeds the new
plan's tolerance for the metric (its `criteria$tolerance`, 0 by default,
so any change of such a number is flagged; timings never are). A number
of both runs that is missing (`NA`) in one of them moved, such as a
read-off that reaches its bar in one run only; a number present in one
run only has `moved = NA`. The attributes `old_bmm` and `new_bmm` name
the bmm versions of the two runs, `n_compared` counts the numbers
present in both, and `mcse_multiplier` is the 2 of the rule.

## Details

Joins the
[`card_numbers()`](https://www.gfrischkorn.org/bmmtools/reference/card_numbers.md)
of an earlier and a later run and flags the numbers that moved by more
than their Monte Carlo error, so a card regenerated after a bmm release
can say what changed.
