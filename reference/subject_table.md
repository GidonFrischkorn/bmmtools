# True and estimated subject values, one row per subject

Puts every subject's generating value next to its posterior point
estimate, one column pair per parameter, with the covariates alongside.
This is the input for a structural equation model of true against
estimated values, for example in lavaan, which recovers a correlation
that the posterior means attenuate by shrinkage (see
[`extract_correlations()`](https://www.gfrischkorn.org/bmmtools/reference/extract_correlations.md)).
bmmtools does not fit that model itself and does not depend on lavaan.

## Usage

``` r
subject_table(
  x,
  fit = NULL,
  point = c("mean", "median"),
  scale = c("link", "natural"),
  links = NULL
)
```

## Arguments

- x:

  A `bmmtools_simulation`, together with `fit`; a
  `bmmtools_simulation_set` from
  [`simulate_components()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_components.md),
  together with its fits; or the result of
  [`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md),
  whose attribute `subject_means` carries the posterior means, medians
  and true values of every cell.

- fit:

  The fit of `x` when `x` is a simulation: a `brmsfit`, or any object
  with an
  [`extract_subject_draws()`](https://www.gfrischkorn.org/bmmtools/reference/extract_subject_draws.md)
  method. For a simulation set, the fits from
  [`fit_components()`](https://www.gfrischkorn.org/bmmtools/reference/fit_components.md),
  one per component, whose terms appear prefixed with the component
  name. Not used with a grid result.

- point:

  The point estimate: the posterior `"mean"` or `"median"` of each
  subject's value.

- scale:

  `"link"` or `"natural"`. On the natural scale the true values and the
  point estimates of each parameter go through
  [`inverse_link()`](https://www.gfrischkorn.org/bmmtools/reference/inverse_link.md);
  covariates are left as they are.

- links:

  A named character vector mapping a term to a link name. `NULL` reads
  the simulation model's link table, or the one the grid stored; without
  one, the natural scale falls back to the link scale with a message.

## Value

A tibble with one row per condition, replication and subject:
`condition` (`NA` for a single simulation), `replication` (`1` for a
single simulation), `id`, then `true_<term>` and `est_<term>` for each
parameter with subject-level draws, then one column per covariate, under
its own name.

## Details

A parameter that does not vary in the simulation has no subject truth;
its `true_<term>` column is the population value for every subject.

On the natural scale the estimate is the inverse link of the point
estimate on the link scale, not the mean of the inverse links, as for
the `point` estimator of
[`extract_correlations()`](https://www.gfrischkorn.org/bmmtools/reference/extract_correlations.md).
For the median the two are the same.

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simulate_recovery(
  bmm::mixture2p(resp_error = "y"),
  pars = c(kappa = log(8), thetat = qlogis(0.75)),
  n_subjects = 100, n_trials = 50,
  sds = c(kappa = 0.3, thetat = 0.5),
  covariates = list(G = c(mean = 0, sd = 1)), seed = 1
)
fit <- bmm::bmm(recovery_formula(sim$model), sim$data, sim$model)
subject_table(sim, fit)

# one table over a grid, from its sidecars
out <- recovery_grid(model, grid, pars, dir = "fits", sds = sds)
subject_table(out, point = "median")
} # }
```
