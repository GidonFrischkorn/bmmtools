# Extract subject-level posterior draws as an array

The per-draw subject values of every parameter with a group-level
effect: the population intercept plus the subject's deviation, summed
within each draw, on the link scale. This is what the `draws` and
`point` estimators of
[`extract_correlations()`](https://www.gfrischkorn.org/bmmtools/reference/extract_correlations.md)
read.

## Usage

``` r
extract_subject_draws(fit, group = NULL, ...)

# Default S3 method
extract_subject_draws(fit, group = NULL, ...)

# S3 method for class 'brmsfit'
extract_subject_draws(fit, group = NULL, ...)
```

## Arguments

- fit:

  A `brmsfit`, and so also a `bmmfit`. Other classes can provide a
  method that returns an array of the same shape.

- group:

  The grouping factor to read. `NULL` uses the fit's only grouping
  factor and errors if there is more than one.

- ...:

  Not used.

## Value

A numeric array with the dimensions `iteration`, `chain`, `id` and
`term`, named in its `dimnames`. `term` is the parameter name, or
`<parameter>_<coefficient>` (`kappa_task1`) for a parameter with several
group-level coefficients and no intercept, as in
[`extract_estimates()`](https://www.gfrischkorn.org/bmmtools/reference/extract_estimates.md).
Parameters the model fixed to a constant (zero posterior variance for
every subject) are dropped. The grouping factor is stored in the
attribute `group`.

## Examples

``` r
if (FALSE) { # \dontrun{
draws <- extract_subject_draws(fit)
dim(draws)
dimnames(draws)$term
} # }
```
