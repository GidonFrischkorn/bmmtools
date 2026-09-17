# Simulate data from a bmm model with known parameters

The generate half of the engine: a model plus population values on the
link scale become a data set and the truth that produced it, so that
[`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
and
[`recover_subjects()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
can score a fit of that data. Subject values are drawn on the link scale
around the population values, converted to the natural scale through
[`inverse_link()`](https://www.gfrischkorn.org/bmmtools/reference/inverse_link.md),
and handed to the model's own `r<model>()` generator, or to the function
you pass as `generator`, which always takes precedence over the built-in
one.

## Usage

``` r
simulate_recovery(
  model,
  pars,
  n_subjects,
  n_trials,
  sds = NULL,
  cors = NULL,
  covariates = NULL,
  tasks = NULL,
  task_col = "task",
  subject_pars = NULL,
  generator = NULL,
  seed = NULL
)
```

## Arguments

- model:

  A `bmmodel`, built with the column names the fit will use, so the
  generated columns match.

- pars:

  Named numeric: population values **on the link scale under bmm's
  parameter names**, one per estimated parameter. A fixed parameter may
  be given too, in which case the generator uses that value. May also be
  a function with no arguments returning such a vector, evaluated under
  `seed`, for random hyperparameters.

- n_subjects, n_trials:

  Subjects, and trials per subject and per row of the generator's layout
  (for `sdt_yn`, per stimulus class).

- sds:

  Named numeric: between-subject standard deviations on the link scale.
  A parameter not named does not vary. `NULL` means no parameter varies.
  May be a function, as `pars`.

- cors:

  A correlation matrix between subject values on the link scale, with
  the parameter and covariate names as dimnames; terms it does not name
  are uncorrelated. `NULL` means none are correlated. May be a function,
  as `pars`.
  [`cors_from_factors()`](https://www.gfrischkorn.org/bmmtools/reference/cors_from_factors.md)
  builds one from factor loadings.

- covariates:

  Observed, error-free person variables drawn jointly with the
  parameters: a named list of `c(mean = , sd = )`. Each becomes a
  subject-constant data column after `id`. A name must be syntactic,
  contain no `_`, and differ from the model's parameters and columns.

- tasks:

  `NULL`, or a character vector of at least two task levels (letters and
  digits) for a design in which every subject does every task. Each free
  parameter then has one value per task, under the term
  `<parameter>_<task_col><level>` (`kappa_task1`), brms's coefficient
  name for `0 + task`. See the details.

- task_col:

  The name of the task column in `data`, `"task"` by default. Used only
  with `tasks`.

- subject_pars:

  A data frame `id`, `term`, `true_value` (link scale) of subject values
  to use instead of drawing them, so that replications can share the
  same simulated people. With `covariates` it must give their values
  too.

- generator:

  A function `(pars, n_trials, model)` returning one subject's rows as a
  data frame with the model's column names; `pars` is a named list on
  the natural scale, fixed parameters included. `NULL` uses the adapter
  bmmtools ships for the model.

- seed:

  A seed applied with
  [`withr::with_seed()`](https://withr.r-lib.org/reference/with_seed.html)
  around the draws and the generator; `NULL` leaves the random number
  generator alone.

## Value

A list of class `bmmtools_simulation` with `data` (a tibble, `id` first,
then any covariates, then the task column), `truth` (a list of tibbles
on the link scale: `population` with `term` and `true_value`; `subjects`
with `id`, `term` and `true_value`; `sd` with `term` and `true_value`;
`cor` with `term`, `var1`, `var2` and `true_value`; `covariates` with
`id`, `term` and `true_value`), the realised `pars`, `sds`, `cors` (the
full matrix over varying parameters then covariates, `NULL` with fewer
than two) and `covariates`, `n_subjects`, `n_trials`, `seed` (`NA` when
none), `model`, `generator`, `tasks` and `task_col` (both `NULL` without
tasks).

## Details

Adapters exist for `sdt_yn`, `sdt_mafc`, `ezdm` (three parameters),
`ddm`, `mixture2p` and `sdm`. Every other model takes a `generator`. The
truth for the subjects and for the SDs lists only the parameters that
vary, because a parameter that does not vary has nothing person-level to
recover.

**Correlated draws.** Subject values and covariates are one multivariate
normal draw, `Z %*% chol(cors)`, over the varying parameters in `pars`
order and then the covariates, scaled by the SDs afterwards. Without
correlations this gives exactly the values an uncorrelated draw gave in
earlier versions, so seeded simulations stay reproducible; adding a
covariate never changes the parameter values. A correlation pair is
named `<a>__<b>`, the two names sorted in the C locale.

**Tasks.** With `tasks`, `pars` and `sds` may name a parameter, which
gives every task that value, or a full term such as `kappa_task2`, which
overrides it for that task; the realised `pars` and `sds` are stored
with full terms. `cors`, `subject_pars` and the truth tables use full
terms only, so a correlation between tasks is, for example,
`kappa_task1__kappa_task2`. The generator is called once per subject and
task with that task's values under the bare parameter names, and
`n_trials` is per task. The task values are cell means: fit them with
`recovery_formula(model, task_col = "task")`. Fixed parameters are the
same in every task. Adding tasks changes the random numbers drawn
compared with a simulation without them; `tasks = NULL` gives exactly
the simulation of earlier versions.

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simulate_recovery(
  bmm::mixture2p(resp_error = "y"),
  pars = c(kappa = log(8), thetat = qlogis(0.75)),
  n_subjects = 30, n_trials = 60,
  sds = c(kappa = 0.3, thetat = 0.5), seed = 1
)
fit <- bmm::bmm(recovery_formula(sim$model), sim$data, sim$model)
recover(fit, sim$truth$population)
recover_subjects(fit, sim$truth$subjects)

# two tasks, kappa lower in the second, correlated .6 across tasks
cors <- diag(2)
dimnames(cors) <- rep(list(c("kappa_task1", "kappa_task2")), 2)
cors[1, 2] <- cors[2, 1] <- 0.6
two_tasks <- simulate_recovery(
  bmm::mixture2p(resp_error = "y"),
  pars = c(kappa = log(8), kappa_task2 = log(5), thetat = qlogis(0.75)),
  n_subjects = 30, n_trials = 60,
  sds = c(kappa_task1 = 0.3, kappa_task2 = 0.3), cors = cors,
  tasks = c("1", "2"), seed = 1
)
formula <- recovery_formula(
  two_tasks$model,
  re_cor = "within", task_col = "task"
)
} # }
```
