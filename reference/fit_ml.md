# Subject-wise maximum likelihood estimates of a bmm model

Fits a model with no pooling — one independent set of parameters per
subject — and returns the estimates in the same tibble shape a
hierarchical fit gives, so that the two can be scored against one truth
and compared. It exists to measure what hierarchical estimation buys,
not to recommend maximum likelihood for inference.

## Usage

``` r
fit_ml(
  model,
  data,
  formula = NULL,
  by = "id",
  method = c("stan", "optim"),
  prior = c("flat", "default"),
  ci_level = 0.95,
  draws = 1000,
  max_abs_link = 20,
  start = NULL,
  nll = NULL,
  file = NULL,
  refit = c("on_change", "never", "always"),
  seed = NULL,
  ...,
  .fitter = NULL
)
```

## Arguments

- model:

  A `bmmodel`, as
  [`bmm::mixture2p()`](https://venpopov.com/bmm/reference/mixture2p.html)
  returns.

- data:

  The data to fit.

- formula:

  An already-reduced `bmmformula` with no group-level terms. `NULL`, the
  default, writes `p ~ 0 + <by>` for every free parameter. Supply one
  for a design `fit_ml()` does not write itself. `method = "optim"`
  builds no formula and does not accept one.

- by:

  The column defining one independent set of parameters. A single column
  in this version.

- method:

  `"stan"` optimises bmm's own generated likelihood through
  `algorithm = "laplace"`, re-implementing nothing, and needs CmdStan.
  `"optim"` maximises bmm's R density for the model directly with
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html) and needs no
  compiler. See the two-routes section.

- prior:

  `"flat"`, the default, puts an empty prior on every free parameter so
  the mode is the maximum-likelihood estimate. `"default"` keeps bmm's
  priors, which makes the mode a penalised estimate rather than an ML
  one. `method = "optim"` has no priors at all and accepts `"flat"`
  only.

- ci_level:

  The interval level. The interval is the Laplace approximation's on the
  Stan route (`ci_method = "laplace"`) and a Wald interval from the
  optimiser's Hessian on the other (`ci_method = "wald"`).

- draws:

  How many draws to take from the Laplace approximation.
  `method = "stan"` only.

- max_abs_link:

  The largest `abs(estimate)` on the link scale that still counts as
  converged.

- start:

  Named link-scale start values for `method = "optim"`, one per free
  parameter. `NULL`, the default, starts at zero, which was measured to
  reach the same optimum as a start at the Stan mode on every subject of
  the reference run.

- nll:

  A function of `(pars, data, model)` returning the negative
  log-likelihood of one subject's `data` as a single number, where
  `pars` is a named list on the natural scale with the model's fixed
  parameters included. Supply it for a model bmmtools carries no density
  for. `method = "optim"` only.

- file, refit:

  Passed to
  [`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md).
  `NULL` uses a temporary file, so nothing is cached between sessions.
  `method = "optim"` does not cache and ignores both: it produces no fit
  object, and at well under a second per call a cache would add an
  invalidation surface without saving anything.

- seed:

  Passed to the fitter. `method = "optim"` is deterministic and ignores
  it.

- ...:

  Passed to
  [`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md)
  and on to the fitter.

- .fitter:

  As in
  [`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md).
  `method = "stan"` only.

## Value

A `bmmtools_ml`: a tibble satisfying the estimates contract, one row per
subject and term, with `level = "subject"` and `estimator = "ml"`.
`attr(x, "ml_cells")` holds one row per subject with its term and
convergence counts.

## Two routes

`method = "stan"` optimises bmm's own generated likelihood with
`algorithm = "laplace"`. It re-implements nothing, so it works for any
bmm model, and it needs a working CmdStan.

`method = "optim"` maximises the likelihood directly with
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) on bmm's R
density for the model. It needs no compiler and is roughly twenty times
faster — 0.35 s against 8 s for 40 subjects of `mixture2p` — but it is
capped at the models bmmtools carries a density for, listed in the error
it raises otherwise, and `nll` is the escape hatch for anything else.

The two routes return the same columns, the same class and the same
`estimator`, and differ only in `ci_method`. Measured on 40 subjects of
`mixture2p` at 100 trials, they agree to Monte Carlo error on the
estimate and their Wald and Laplace standard errors agree to a mean
ratio of 1.000 (range 0.95 to 1.05), giving the same `coverage` and
`calibration_slope` to three decimals. The route is not a moderator of a
scored number.

## One fit, not one per subject (`method = "stan"`)

Under `p ~ 0 + <by>` every free parameter belongs to one subject and
bmm's remaining parameters are constants, so the log posterior is a sum
of per-subject terms and the joint mode is the vector of per-subject
modes. The Stan route therefore runs a single fit for the whole data
set. This was checked rather than assumed: over 40 subjects, the joint
Stan mode and independent optimisation of bmm's own density agreed to
3.0e-06 on the link scale.

The cost is that convergence is no longer separable by subject — one
subject at a boundary can stall the whole optimisation — so the per-cell
verdict in `ml_cells` is the range check described below. The `optim`
route does not have this property: it is one independent optimisation
per subject.

## What counts as converged

`converged` is `TRUE` where the estimate is finite and
`abs(estimate) <= max_abs_link` on the link scale, and, on the `optim`
route, where [`optim()`](https://rdrr.io/r/stats/optim.html) also
returned its success code. brms keeps no optimiser return code on a
Laplace fit, so on the Stan route the range check is the only verdict
available — and not a weak one: fitting a single subject, Stan reported
convergence at a point 0.064 log-likelihood short of the mode.

A row that fails is kept with `estimate = NA` rather than dropped,
because `metric_bias()` and `metric_rmse()` drop incomplete pairs
silently: dropping would compare the Bayesian estimator on every subject
against maximum likelihood on the easiest ones and report the difference
as the estimator's. Report `n` alongside any comparison.

On the `optim` route a Hessian that cannot be inverted gives a row with
no interval: `ci_low` and `ci_high` are `NA` while the estimate and
`converged` stand, so the subject counts in `bias`, `rmse` and `r` and
drops out of `coverage` only.

`rhat`, `ess_bulk` and `ess_tail` are always `NA`. Neither route
produces MCMC draws, so the values `posterior` would return describe the
approximation and not the estimate.

## See also

[`recover_subjects()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
to score the result against a known truth, and
[`extract_estimates()`](https://www.gfrischkorn.org/bmmtools/reference/extract_estimates.md)
for the hierarchical side of the comparison. `bind_rows()` of the two is
the join; there is no separate comparison function.
[`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md)
with `ml = TRUE` runs both fits on every cell of a design and does the
join itself.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- bmm::mixture2p(resp_error = "y")
sim <- simulate_recovery(
  model,
  pars = c(kappa = 2, thetat = 1), sds = c(kappa = 0.3, thetat = 0.5),
  n_subjects = 40, n_trials = 100, seed = 1
)
ml <- fit_ml(model, sim$data)
summary(recover_subjects(ml, sim$truth$subjects))

# the same thing without a compiler, and without Stan
ml2 <- fit_ml(model, sim$data, method = "optim")
} # }
```
