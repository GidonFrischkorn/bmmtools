# Simulate data from a bmm model with known parameters

The generate half of the engine: a model plus population values on the
link scale become a data set and the truth that produced it, so that
[`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
and
[`recover_subjects()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
can score a fit of that data. Subject values are drawn on the link scale
around the population values, converted to the natural scale through
[`inverse_link()`](https://www.gfrischkorn.org/bmmtools/reference/inverse_link.md)
(a `softmax` term excepted, see the details), and handed to the model's
own `r<model>()` generator, or to the function you pass as `generator`,
which always takes precedence over the built-in one.

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
  coding = c("cell", "contrast"),
  contrasts = NULL,
  subject_pars = NULL,
  formula = NULL,
  trial_design = NULL,
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
  (for `sdt_yn`, per stimulus class). With a `trial_design`, `n_trials`
  is the number of design rows per subject, which for a count model are
  rows of data, not trials (see the section "Trial design").

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

- coding:

  What the task column means to a fit. `"cell"`, the default, leaves it
  as it was: a fit of `0 + task` estimates one value per task and the
  truth is those values. `"contrast"` attaches `contrasts` to the task
  column, so a fit of `1 + task` estimates an intercept and
  `length(tasks) - 1` contrasts, and `truth` holds those instead:
  `<parameter>` for the intercept and `<parameter>_<task_col><j>` for
  the contrasts, with the SD and correlation truths transformed with
  them.

  **What is drawn does not change.** Subject values are drawn around the
  cell means exactly as under `"cell"`, from the same RNG stream, so
  `sds`, `cors` and `subject_pars` keep their meaning and a simulation
  is comparable across codings. Only the design the fit sees, and the
  truth it is scored against, differ.

- contrasts:

  The contrast matrix for `coding = "contrast"`: `length(tasks)` rows
  and one column fewer, or a function `(n)` returning one.
  [stats::contr.treatment](https://rdrr.io/r/stats/contrast.html) by
  default, so the intercept is the first task and each contrast a
  difference from it. For a study in which bmm's default prior should
  mean the same for every task, use an orthonormal coding such as
  `bayestestR::contr.equalprior`, whose intercept is the grand mean and
  whose contrasts are uncorrelated with it. Ignored under
  `coding = "cell"`.

- subject_pars:

  A data frame `id`, `term`, `true_value` (link scale) of subject values
  to use instead of drawing them, so that replications can share the
  same simulated people. With `covariates` it must give their values
  too.

- formula:

  `NULL`, or the `bmmformula` the model is fitted with, for a built-in
  generator that reads it: only the custom version of `m3`, whose
  activation formulas exist nowhere else (see the section "The m3
  model"). A formula for any other model, or with a `generator`, is an
  error: a generator you write closes over its own.
  [`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md),
  [`recovery_component()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_component.md)
  and [`sbc()`](https://www.gfrischkorn.org/bmmtools/reference/sbc.md)
  pass on the formula they fit with.

- trial_design:

  `NULL`, or the per-trial variables of the design (set sizes,
  non-target locations, option counts): what a generator needs per trial
  and the model reads from the data but does not hold. For a count
  model, the number of trials (and stimulus) of each row instead; see
  the section "Trial design". One of

  - a data frame of `n_trials` rows, one trial list used for every
    subject (and every task);

  - a data frame with an `id` column holding the subjects `1` to
    `n_subjects`, `n_trials` rows each, one trial list per subject;

  - a function `(n_trials)` returning a data frame of `n_trials` rows,
    called once per subject and task under `seed`, for designs drawn at
    random, such as non-target locations.

  Each generator call receives its rows as `trial_design`, and the same
  rows are bound into `data` after `id`, the covariates and the task
  column. See the section "Trial design".

- generator:

  A function `(pars, n_trials, model)` returning one subject's rows as a
  data frame with the model's column names; `pars` is a named list on
  the natural scale, fixed parameters included, except a parameter whose
  link is `softmax`, which arrives on the link scale (see the details).
  With a `trial_design` it is called as
  `(pars, n_trials, model, trial_design)` and must have that fourth
  argument or `...`. `NULL` uses the adapter bmmtools ships for the
  model.

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
tasks), `trial_design` (`NULL`, the data frame, or a function's text),
`trial_design_columns` (the columns it contributed to `data`) and
`formula` (`NULL`, or the formula that reached the generator, each
element deparsed to text).

## Details

Adapters exist for `sdt_yn`, `sdt_mafc`, `ezdm` (three parameters),
`ddm`, `cswald` (both versions), `mixture2p`, `sdm`, `mixture3p`, `imm`
(`full`, `bsc` and `abc`) and `m3` (`ss`, `cs` and `custom`). Every
other model takes a `generator`. The `mixture3p` and `imm` adapters need
a `trial_design`; see the section "Trial design". The custom `m3` needs
`formula`; see the section "The m3 model". The truth for the subjects
and for the SDs lists only the parameters that vary, because a parameter
that does not vary has nothing person-level to recover.

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

## Trial design

A generator learns one subject's parameters and a number of trials. A
model whose responses depend on what was shown on each trial —
`mixture3p` and `imm` need a set size and non-target locations — gets
those through `trial_design`. The engine hands each generator call its
rows and binds the same rows into `data`, so the fit sees exactly the
trials the responses were drawn from; a generator returns the response
columns only, and returning a design column is an error. A design column
must differ from the model's response columns, the covariates and the
task column.

**`mixture3p` and `imm`.** Their adapters read, from the columns the
model names: the set size, when `set_size` is a column name (a number is
every trial's set size); the non-target locations `nt_features`,
relative to the target, in radians; and for `imm` `full` and `bsc` the
distances `nt_distances`. A trial of set size `k` has its lures in the
first `k - 1` non-target columns, which must not be `NA`; columns beyond
are ignored and are `NA` by bmm's convention. Each trial is one call of
[`bmm::rmixture3p()`](https://venpopov.com/bmm/reference/mixture3p_dist.html)
or [`bmm::rimm()`](https://venpopov.com/bmm/reference/IMMdist.html),
with the weights bmm's likelihood gives that trial.

**The `softmax` link.** `mixture3p`'s `thetat` and `thetant` are log
weights against a guessing weight of 0: on a trial with lures,
`(p_mem, p_nt, p_guess)` is the softmax of `(thetat, thetant, 0)`, the
lures sharing `p_nt` equally, and on a set-size-1 trial there is no lure
weight. Since no inverse of one term gives its natural value, a
parameter whose link is `softmax` reaches the generator on the link
scale, and
[`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
scores it there (its `scale` column says so).

**`m3` with its numbers of options as columns.** When `num_options`
names columns, those columns are a trial design and each of its rows is
one trial: see the section "The m3 model".

**Count models.** One row of `sdt_yn`, `sdt_mafc`, `ezdm` or `m3` with
numbers of options holds many trials. Without a design, every row holds
`n_trials` of them. A design is optional for these models; with one,
each design row is one row of `data` and carries its own number of
trials, so `n_trials` counts the design rows. The trial count is read
from the model's `n_trials` column, and for `m3`, which has none, from a
design column `n_trials`. `sdt_yn` also reads its stimulus column (0 or
1), so that 100 signal and 50 noise trials are
`data.frame(stimulus = c(1, 0), n_trials = c(100, 50))` under the
model's column names. A count must be a whole number of at least 1, and
of at least 3 for `ezdm`, whose generator refuses fewer; a design that
breaks this is refused before anything is drawn.

Without `trial_design` nothing changes: the generator is called with
three arguments, as in bmmtools 0.2.0, and a seeded simulation gives the
same data, and so the same
[`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md)
key. A data frame draws no random numbers; a function draws its own,
interleaved with the generator's, once per subject and task. What a
function returns can only be checked once it has been called, so an
error in it (wrong row count, a missing column) stops the simulation
after some random numbers have been drawn; under `seed` that changes
nothing, and without one the random number stream has moved.

## The m3 model

`m3`'s response is a count per category, in the columns `resp_cats`
names. With numbers of options on the model
(`num_options = c(1, 4, 5)`), a subject has **one row of `n_trials`
trials** per task, drawn with one
[`bmm::rm3()`](https://venpopov.com/bmm/reference/m3dist.html) call, as
`sdt_mafc` has one row of counts. With `num_options` naming columns, the
options differ between trials, so those columns come from `trial_design`
and **each row is one trial**: `n_trials` rows per subject, each
counting one response in its category. The probabilities are bmm's: each
category's activation, exponentiated under `choice_rule = "softmax"`,
times its number of options, normalised. The parameters' links are
elementwise, and which links the model has depends on the choice rule,
so take `pars` on the link scale of the model as built:
[`bmm::parameters()`](https://venpopov.com/bmm/reference/parameters.html)
lists them for `ss` and `cs`, and the model's `links` for `custom`.

The versions `ss` and `cs` take their activation formulas from bmm. The
`custom` version knows its activations only from the formula it is
fitted with, so it needs `formula`, holding one activation per category
(`corr ~ b + a + c`), and a model built with a link for each parameter
they use (`m3(..., links = list(c = "log", a = "log"))`); without links
bmm knows no parameter to give a value to until it fits. Every parameter
with a link must appear in an activation. Extra parameter formulas in
`formula` (`c ~ 1 + (1 | id)`) are ignored here, and a grid records only
the activations.

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

# a per-trial design: a condition the generator reads on every trial
shifted <- function(pars, n_trials, model, trial_design) {
  mu <- ifelse(trial_design$cond == "a", 0, 0.5)
  data.frame(y = vapply(mu, function(m) {
    bmm::rmixture2p(1, mu = m, kappa = pars$kappa, p_mem = pars$thetat)
  }, numeric(1)))
}
sim <- simulate_recovery(
  bmm::mixture2p(resp_error = "y"),
  pars = c(kappa = log(8), thetat = qlogis(0.75)),
  n_subjects = 10, n_trials = 40,
  trial_design = function(n_trials) {
    data.frame(cond = sample(c("a", "b"), n_trials, replace = TRUE))
  },
  generator = shifted, seed = 1
)

# mixture3p: set size 2 to 4, non-target locations drawn per trial
lures <- function(n_trials) {
  ss <- sample(2:4, n_trials, replace = TRUE)
  nt <- matrix(runif(3 * n_trials, -pi, pi), n_trials)
  nt[col(nt) >= ss] <- NA
  data.frame(ss = ss, nt1 = nt[, 1], nt2 = nt[, 2], nt3 = nt[, 3])
}
sim <- simulate_recovery(
  bmm::mixture3p(
    resp_error = "y", nt_features = paste0("nt", 1:3), set_size = "ss"
  ),
  pars = c(kappa = log(8), thetat = 1.5, thetant = 0),
  n_subjects = 10, n_trials = 60, trial_design = lures, seed = 1
)

# a custom m3: the activation formulas reach the adapter through formula
model <- bmm::m3(
  resp_cats = c("corr", "other", "npl"), num_options = c(1, 4, 5),
  choice_rule = "simple", version = "custom",
  links = list(c = "log", a = "log")
)
formula <- bmm::bmf(
  corr ~ b + a + c, other ~ b + a, npl ~ b,
  c ~ 1 + (1 | id), a ~ 1 + (1 | id)
)
sim <- simulate_recovery(
  model,
  pars = c(c = log(3), a = log(0.5)),
  n_subjects = 30, n_trials = 100, sds = c(c = 0.3),
  formula = formula, seed = 1
)
} # }
```
