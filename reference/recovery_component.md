# Describe one model of a multi-model simulation

A component is one model with its own population values, trials, tasks,
generator and formula.
[`simulate_components()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_components.md)
simulates several of them for the same people, so that parameters of
different models can be correlated across subjects, and
[`fit_components()`](https://www.gfrischkorn.org/bmmtools/reference/fit_components.md)
fits each on its own.

## Usage

``` r
recovery_component(
  model,
  pars,
  n_trials,
  sds = NULL,
  tasks = NULL,
  task_col = "task",
  generator = NULL,
  formula = NULL,
  name
)
```

## Arguments

- model:

  A `bmmodel`, built with the column names the fit will use.

- pars, n_trials, sds, tasks, task_col, generator:

  As in
  [`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md),
  for this component. `pars` and `sds` may be functions with no
  arguments, evaluated by
  [`simulate_components()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_components.md)
  under its seed.

- formula:

  `NULL`, or the `bmmformula`
  [`fit_components()`](https://www.gfrischkorn.org/bmmtools/reference/fit_components.md)
  fits this component with. `NULL` uses
  [`recovery_formula()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_formula.md)
  with the task column when there are tasks.

- name:

  The component's name, which prefixes its terms: `name = "m3"` turns
  `c_task1` into `m3_c_task1`. A single syntactic name without `_` or
  `.`, unique within a set.

## Value

A list of class `bmmtools_component` with the arguments as validated;
`tasks` and `task_col` are `NULL` without tasks.

## Details

What can be checked without drawing is checked here: the model, the
number of trials, the tasks, that a generator exists, the formula, and
numeric `pars` and `sds` against the model. Correlations and covariates
belong to the set and are given to
[`simulate_components()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_components.md).

## Examples

``` r
if (FALSE) { # \dontrun{
recovery_component(
  bmm::sdm(resp_error = "y"),
  pars = c(c = log(3), kappa = log(5)),
  n_trials = 60, sds = c(c = 0.4), name = "sdm"
)
} # }
```
