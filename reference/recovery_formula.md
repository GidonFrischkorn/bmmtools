# The default recovery formula: every free parameter gets a random intercept

The formula the validation scripts in bmm converged on:
`<parameter> ~ 1 + (1 | id)` for every parameter the model estimates,
and nothing for the ones it fixes.
[`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md)
uses it when no formula is given; a user who wants to change one term
can start from it.

## Usage

``` r
recovery_formula(
  model,
  group = "id",
  re_cor = c("none", "within", "all"),
  task_col = NULL
)
```

## Arguments

- model:

  A `bmmodel`.

- group:

  The grouping variable, `"id"` by default.

- re_cor:

  Which random effects are correlated. `"none"` gives independent random
  effects: `(1 | id)`, or `(0 + task || id)` with tasks. `"within"`
  correlates the tasks of each parameter, `(0 + task | id)`, and needs
  `task_col`. `"all"` gives one correlation matrix across every
  parameter (and task), `(1 | p | id)` or `(0 + task | p | id)`, which
  is what the `model` estimator of correlation recovery reads.

- task_col:

  `NULL`, or the task column of a simulation with `tasks` (see
  [`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md)).
  Each free parameter then gets one population value per task,
  `<parameter> ~ 0 + task`, and one subject-level effect per task: cell
  means, whose terms match the truth of the simulation.

## Value

A `bmmformula`.

## Examples

``` r
if (FALSE) { # \dontrun{
recovery_formula(bmm::mixture2p(resp_error = "y"))
recovery_formula(bmm::mixture2p(resp_error = "y"), re_cor = "all")
recovery_formula(
  bmm::mixture2p(resp_error = "y"),
  re_cor = "within", task_col = "task"
)
} # }
```
