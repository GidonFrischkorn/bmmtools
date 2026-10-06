# Plot what the priors say the data should look like

The picture a prior check is read through. `"density"` draws each
prior-predictive draw as its own thin line, so the spread *between*
draws stays visible where a single pooled density would hide it;
`"histogram"` is for a discrete response, where a density is misleading;
`"statistic"` plots the summary table instead of the distribution.

## Usage

``` r
plot_prior_check(
  x,
  type = c("density", "histogram", "statistic"),
  observed = TRUE,
  draws = 50,
  facet_by = "response",
  ...
)
```

## Arguments

- x:

  A `bmmtools_prior_check` object from
  [`prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/prior_check.md).

- type:

  `"density"`, `"histogram"` or `"statistic"`.

- observed:

  Overlay the observed response from the data the check was run over.
  Dropped without comment when the model's response column is not in
  that data, since a prior check over a design skeleton has nothing
  observed to show. Ignored by `"statistic"`.

- draws:

  Prior-predictive draws to show per set, capped at what the object
  holds.

- facet_by:

  A column of `x` to make panels from, or `NULL`. `"response"` by
  default, so a model with several observables does not share an axis.

- ...:

  Not used.

## Value

A `ggplot` object.

## Examples

``` r
plot_prior_check(prior_check_sdt_yn, type = "histogram")

plot_prior_check(prior_check_sdt_yn, type = "statistic")
```
