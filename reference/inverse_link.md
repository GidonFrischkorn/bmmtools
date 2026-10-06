# Transform values from the link scale to the natural scale

Cognitive measurement models are estimated on a link scale that keeps
parameters unconstrained — a memory precision is sampled as
`log(kappa)`, a mixture weight as `logit(theta)` — while readers
interpret them on the natural scale. `inverse_link()` is the
back-transform, applied by
[`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
to estimates, interval bounds and generating values alike so that both
sides of a recovery comparison are on the same scale.

## Usage

``` r
inverse_link(x, link)
```

## Arguments

- x:

  A numeric vector of values on the link scale. `NA` is passed through
  at its own position.

- link:

  A single string naming the link: one of `"identity"`, `"log"`,
  `"softplus"`, `"log1p"`, `"logm1"`, `"inverse"`, `"sqrt"`, `"logit"`,
  `"probit"`, `"tan_half"`, `"loglog"` or `"cloglog"`. `NULL` is a
  synonym for `"identity"`.

## Value

A double vector the same length as `x`, on the natural scale.

## Details

The vocabulary and the closed forms are taken from bmm's internal
`link_transform()` (`bmm/R/helpers-parameters.R`), so the two cannot
disagree about what a link name means. bmmtools ships its own copy
because `link_transform()` is not exported.

Values outside a link's natural domain are transformed as R computes
them rather than being replaced by `NA`: `inverse_link(0, "inverse")` is
`Inf` and `inverse_link(-2, "sqrt")` is `4`. In a scoring pipeline a
visible `Inf` is safer than a silent `NA`, which would be absorbed by a
metric's missing-data guard and reported as a smaller sample rather than
as a problem.

## Examples

``` r
inverse_link(c(-1, 0, 1), "log")
#> [1] 0.3678794 1.0000000 2.7182818
inverse_link(0, "logit")
#> [1] 0.5

# a precision sampled on the log scale, back on its own scale
inverse_link(c(1.6, 2.3), "log")
#> [1] 4.953032 9.974182
```
