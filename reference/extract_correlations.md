# Extract between-subject correlations from a fit

Three estimates of how two subject-level parameters, or a parameter and
an observed covariate, correlate across subjects:

## Usage

``` r
extract_correlations(
  fit,
  estimator = c("model", "draws", "point"),
  pairs = NULL,
  covariates = NULL,
  group = NULL,
  scale = c("link", "natural"),
  links = NULL,
  ci_level = 0.95,
  converged = NULL,
  ...
)
```

## Arguments

- fit:

  A `brmsfit` (so also a `bmmfit`), or any object with an
  [`extract_subject_draws()`](https://www.gfrischkorn.org/bmmtools/reference/extract_subject_draws.md)
  method. Or a set of separate fits, one per component: the result of
  [`fit_components()`](https://www.gfrischkorn.org/bmmtools/reference/fit_components.md),
  or a list of fits named by component; see "Separate fits" below.

- estimator:

  One or more of `"model"`, `"draws"` and `"point"`. All three by
  default.

- pairs:

  `NULL` for every pair of subject terms and covariates, or a character
  vector of pair terms such as `"kappa__thetat"`, in either order.

- covariates:

  `NULL`; a character vector of columns of `fit$data` that are constant
  within subject; or a data frame with the grouping column and one
  numeric column per covariate, one or several rows per subject.

- group:

  The grouping factor. `NULL` uses the fit's only one.

- scale:

  `"link"` or `"natural"`. On the natural scale each term is
  back-transformed with
  [`inverse_link()`](https://www.gfrischkorn.org/bmmtools/reference/inverse_link.md)
  before correlating.

- links:

  A named character vector mapping a term to a link name. `NULL` reads
  the link table of a `bmmfit`; without one, the natural scale falls
  back to the link scale with a message.

- ci_level:

  The interval mass, a single number strictly between 0 and 1.

- converged:

  Whether the fit passed the convergence gate. `NULL` computes it as
  [`extract_estimates()`](https://www.gfrischkorn.org/bmmtools/reference/extract_estimates.md)
  does for a `brmsfit` and leaves it `NA` for other classes; a logical
  scalar is used as given.

- ...:

  Not used.

## Value

A tibble with the columns `term`, `var1`, `var2`, `estimator`,
`estimate`, `ci_low`, `ci_high`, `ci_method`, `ci_level`, `rhat`,
`ess_bulk`, `ess_tail`, `scale`, `n` (subjects) and `converged`. `term`
is `var1__var2`, the two names sorted in the C locale.

## Details

- `"model"`: the group-level correlation the model estimates, as in
  `extract_estimates(level = "cor")`. It exists only for parameters
  whose random effects share a correlation block (`(1 | p | id)`).

- `"draws"`: in each posterior draw, the Pearson correlation of the
  subject values across subjects; summarised like any other posterior
  (median, equal-tailed interval, rhat, ESS).

- `"point"`: the Pearson correlation of the subjects' posterior means,
  with a Fisher-z interval over the subjects.

They answer different questions and do not agree. When the model cannot
separate two parameters within a subject, their posterior means inherit
that trade-off and `point` is biased even where the true correlation is
0; `draws` is not. `point` is also attenuated by shrinkage. See
[`recover_correlations()`](https://www.gfrischkorn.org/bmmtools/reference/recover_correlations.md)
for scoring them against the truth.

**The model estimator** is on the link scale only. Under
`scale = "natural"` its rows are left out and a message says so. A pair
involving a covariate never has a model row, and neither does a pair the
model does not correlate.

**The natural scale** for `point` takes the inverse link of each
subject's posterior mean on the link scale, not the mean of the inverse
links. Covariates are never transformed.

**Intervals.** `ci_method` is `"eti"` for `model` and `draws` and
`"fisher_z"` for `point`, whose interval treats the posterior means as
if they were observed values; its `rhat` and ESS are `NA`.

**Missing values.** A correlation is `NA`, never 0, with fewer than
three subjects or when one side has no spread across subjects. For
`draws` that is decided per draw, and the summary uses the draws that
remain.

## Separate fits

With a set of fits, the subject draws of every fit are bound under
prefixed terms (`m3_c_task1`), subject by subject. Draws are paired by
iteration and chain, cut to the smallest counts among the fits with a
message. `model` rows exist only within a fit; `converged` of a pair is
whether every fit it involves converged, and may be given as a logical
named by component. `covariates` must be a data frame, because the
components' data carry none.

No fit knows about the correlation between its parameters and another
model's, so correlations across fits are attenuated by the reliabilities
`rel_a` and `rel_b` of the two subject estimates: `point` is roughly
`rho * sqrt(rel_a * rel_b)` and `draws` roughly `rho * rel_a * rel_b`.
Neither recovers `rho`. A structural equation model of true against
estimated values, from
[`subject_table()`](https://www.gfrischkorn.org/bmmtools/reference/subject_table.md),
does; so would a joint multivariate fit, which bmmtools does not support
yet.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- bmm::bmm(recovery_formula(model, re_cor = "all"), data, model)
extract_correlations(fit)
extract_correlations(fit, estimator = "draws", covariates = "G")
} # }
```
