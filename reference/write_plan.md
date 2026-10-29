# Write a validation plan to a file and read it back

A plan is a list, so it round-trips through
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html). Writing it next to
the card it produced keeps a record of what was planned that can be
compared with a later version by its `hash`.

## Usage

``` r
write_plan(plan, file)

read_plan(file)
```

## Arguments

- plan:

  A `bmmtools_validation_plan`.

- file:

  A path ending in `.rds`.

## Value

`write_plan()` returns `file` invisibly; `read_plan()` the plan.

## Examples

``` r
plan <- validation_plan(
  structure(
    list(links = list(kappa = "log", thetat = "logit")),
    class = c("bmmodel", "circular", "mixture2p")
  ),
  seed = 1
)
file <- tempfile(fileext = ".rds")
write_plan(plan, file)
identical(read_plan(file)$hash, plan$hash)
#> [1] TRUE
```
