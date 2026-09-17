# Run a parameter-recovery grid

The loop the validation scripts in bmm wrote by hand, with the
requirements their overnight runs taught: one durable file per cell
written as soon as the cell finishes, a resume that reads those files
and needs no temporary state, a smoke mode in its own directory, a
preflight fit that catches a compile or init error before any cell runs,
and cells ordered so that the first completed block spans the design
rather than exhausting one level.

## Usage

``` r
recovery_grid(
  model,
  grid,
  pars,
  dir,
  reps = 1,
  sds = NULL,
  cors = NULL,
  covariates = NULL,
  tasks = NULL,
  task_col = "task",
  formula = NULL,
  prior = NULL,
  generator = NULL,
  seed = NULL,
  subjects = c("redraw", "fixed"),
  re_cor = c("none", "within", "all"),
  scale = c("natural", "link"),
  levels = c("population", "subject"),
  correlations = NULL,
  cor_scale = "link",
  ml = FALSE,
  smoke = FALSE,
  preflight = TRUE,
  ...,
  .fitter = NULL
)
```

## Arguments

- model:

  A `bmmodel`, or a function of one grid row (a one-row data frame)
  returning one, for models whose constructor depends on the design. Or
  a list of
  [`recovery_component()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_component.md)s,
  or a function of the row returning such a list, to simulate several
  models for the same people and fit each on its own; see the section
  "Components".

- grid:

  A data frame with the columns `n_subjects` and `n_trials` (with
  components, `n_subjects` only; see "Components" for their columns). A
  column named after a parameter gives that cell's population value on
  the link scale, overriding `pars`; a column `sd_<parameter>` overrides
  `sds`; a column `cor_<a>__<b>` sets the correlation of two parameters
  or covariates (the names in either order). With `tasks`, these columns
  may also use full task terms (`kappa_task2`, `sd_kappa_task2`,
  `cor_kappa_task1__kappa_task2`); a bare parameter column sets every
  task, and a full-term column in the same row overrides it for its
  task. A SimDesign design is a data frame and works as is.

- pars, sds, cors:

  Defaults for every cell, as in
  [`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md).
  Each may also be a `function(row)` of the one-row grid data frame,
  evaluated under the cell's seed, which draws new hyperparameters for
  every data set; the grid columns above are applied to its result.

- dir:

  Directory for the per-cell files; created if missing.

- reps:

  Replications per cell.

- covariates:

  As in
  [`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md),
  the same for every cell.

- tasks, task_col:

  As in
  [`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md),
  the same for every cell. With `tasks`, the default formula is
  `recovery_formula(model, re_cor = re_cor, task_col = task_col)`.

- formula:

  A `bmmformula`; `NULL` means
  [`recovery_formula()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_formula.md)
  of the cell's model; or a `function(row)` of the one-row grid data
  frame returning a `bmmformula`, called once per row, for designs whose
  formula depends on the row.

- prior:

  Passed to
  [`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md);
  with components, `NULL` or a list with a prior per component.

- generator:

  As in
  [`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md).

- seed:

  Master seed. Each cell derives its own from it and the cell's row and
  replication, so a cell is reproducible on its own. `NULL` leaves
  everything unseeded and records `NA`.

- subjects:

  `"redraw"` draws new subject values in every replication; `"fixed"`
  draws them once per row and reuses them, so replications become a
  simulated retest of the same people. The realised population values,
  SDs, correlations and covariate values of replication 1 are reused
  too.

- re_cor:

  Passed to
  [`recovery_formula()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_formula.md)
  when `formula` is `NULL`.

- scale:

  The scale recovery is scored on, as in
  [`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md);
  natural by default.

- levels:

  The estimate levels extracted from every fit and scored: one or more
  of `"population"`, `"subject"` and `"sd"` (the between-subject
  standard deviations, always scored on the link scale).

- correlations:

  `NULL` (the default) for no correlations, or one or more of the
  estimators of
  [`extract_correlations()`](https://www.gfrischkorn.org/bmmtools/reference/extract_correlations.md),
  `"model"`, `"draws"` and `"point"`, to extract and score the
  between-subject correlations of every cell, pairs with covariates
  included.

- cor_scale:

  The scale the correlations are extracted and scored on, `"link"` or
  `"natural"`. The `"model"` estimator exists on the link scale only.

- ml:

  `FALSE`, the default, fits every cell hierarchically only. `TRUE` also
  fits every cell subject by subject with
  [`fit_ml()`](https://www.gfrischkorn.org/bmmtools/reference/fit_ml.md)
  on the same simulated data, and scores both estimators against one
  truth at the subject level; the ML rows carry `estimator = "ml"`. A
  named list gives arguments to
  [`fit_ml()`](https://www.gfrischkorn.org/bmmtools/reference/fit_ml.md)
  (`prior`, `draws`, `max_abs_link`, `ci_level`, `formula`, `.fitter`,
  and through its `...` the fitter's own, such as
  `refresh = 0, silent = 2` to quiet the optimiser); `model`, `data`,
  `file`, `refit`, `seed` and `by` are the grid's. Needs `"subject"` in
  `levels`. See the section "Subject-wise ML".

- smoke:

  `TRUE` runs the first two rows with two replications into
  `<dir>/smoke`, so a smoke run never overwrites a full one.

- preflight:

  Run the first cell once with one chain and 200 iterations into
  `<dir>/preflight` before the loop. Skipped when the first cell already
  has a cached fit or a sidecar that can be used.

- ...:

  Passed to
  [`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md)
  and so to the fitter: `chains`, `iter`, `warmup`, `backend`, `cores`,
  ...

- .fitter:

  As in
  [`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md).

## Value

A `bmmtools_recovery` object over every cell that produced a fit, at the
requested `levels`, with a `condition` column (`"row-<i>"`) naming the
grid row. [`summary()`](https://rdrr.io/r/base/summary.html) groups by
condition. Its attributes:

- `cells`: a tibble with one row per cell (`condition`, `replication`,
  `n_subjects`, `n_trials`, `seed`, `file`, `status`, `elapsed` in
  seconds, `converged`).

- `grid`: the grid as run.

- `correlations`: with `correlations` requested, a
  `bmmtools_cor_recovery` from
  [`recover_correlations()`](https://www.gfrischkorn.org/bmmtools/reference/recover_correlations.md)
  over every cell, with `condition` and `replication`; otherwise `NULL`.

- `subject_means`: a tibble with one row per cell, subject and term:
  `condition`, `replication`, `id`, `term`, `covariate` (whether the
  term is a covariate), the posterior `mean` and `median` on the link
  scale, and `true_value`, the generating value on the link scale (a
  covariate's value is its own mean, median and true value). Its
  attribute `links` is the model's link table.
  [`subject_table()`](https://www.gfrischkorn.org/bmmtools/reference/subject_table.md)
  turns it into one row per subject.

- `ml_cells`: with `ml`, a tibble with one row per cell (`condition`,
  `replication`, `status`, `elapsed`, `n_subjects`, `n_converged`,
  `converged`, `file`) for the ML fit; see "Subject-wise ML". `cells`
  stays one row per cell and describes the hierarchical fit.

## Details

A cell whose fit errors is recorded with `status = "error"` and the grid
goes on; a warning at the end names the failed cells. A cell whose data
cannot be generated stops the grid, because that is a design error
rather than a sampler accident. Cell files are
`cell-<row>-rep-<rep>-sim.rds` (data and truth) and
`cell-<row>-rep-<rep>.rds` with its `.key` (the fit, through
[`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md));
a resume reads an existing simulation file rather than regenerating, so
a cached fit and its truth always match.

**The sidecar.** While a fit is in memory, everything scored from it is
written to `cell-<row>-rep-<rep>-est.rds`: the estimates at `levels`,
the correlations, the subject means, and what they were extracted with
(the cell's cache key, the bmmtools version, `levels`, `correlations`,
`cor_scale`). A resume uses the sidecar without reading the fit when the
key and the version match and it holds every requested level and
estimator, so the fit files may be deleted once a grid has run.
Otherwise the cell goes through
[`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md)
again, which reuses a cached fit, and the sidecar is rewritten.

## Subject-wise ML

With `ml`, every cell is fitted twice on the same simulated data: the
hierarchical fit, and one
[`fit_ml()`](https://www.gfrischkorn.org/bmmtools/reference/fit_ml.md)
call with no pooling, cached as `cell-<row>-rep-<rep>-ml.rds` under the
cell's seed. The ML rows are bound to the hierarchical subject rows
before scoring, so the result has both estimators at the subject level,
`"bayes"` and `"ml"`, against one truth;
[`summary()`](https://rdrr.io/r/base/summary.html) reports them as
separate rows and `plot_recovery(color_by = "estimator")` colours them.
An ML fit has no population or `sd` level, so those levels stay
hierarchical only, and its subject means are not added to
`subject_means`. The `...` of the grid (`chains`, `iter`, `backend`,
`refresh`, ...) belong to the hierarchical fit and do not reach
[`fit_ml()`](https://www.gfrischkorn.org/bmmtools/reference/fit_ml.md);
`ml = list(draws = , refresh = 0)` and the other
[`fit_ml()`](https://www.gfrischkorn.org/bmmtools/reference/fit_ml.md)
arguments do.

The sidecar records the ML request (its prior, draws, interval level,
link-scale range and formula) with the ML rows, and a resume reuses them
only when the request is the same; a change to `max_abs_link` or
`ci_level` re-extracts from the cached ML fit without refitting, and a
change to `prior`, `draws` or `formula` refits through
[`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md).
A cell whose ML fit fails keeps its hierarchical rows, is recorded with
`status = "error"` in `ml_cells`, and is named in a warning at the end;
its sidecar is written without ML rows, so a resume tries the ML fit
again. There is no ML preflight: a compile error of the no-pooling model
shows up as every cell's ML status.

Within a cell, a subject whose ML estimate leaves the link-scale range
keeps its row with `estimate = NA` and `converged = FALSE`, as
[`fit_ml()`](https://www.gfrischkorn.org/bmmtools/reference/fit_ml.md)
documents; the recovery rows carry that verdict per subject and term,
and `ml_cells` counts the subjects per cell. Report `n` next to any
comparison of the two estimators.

## Components

With `model` a list of
[`recovery_component()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_component.md)s
(or a `function(row)` returning one, the same components by name in
every row), each cell simulates one set with
[`simulate_components()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_components.md)
under the cell's seed and fits every component on its own through
[`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md).
`pars`, `sds`, `generator`, `tasks`, `task_col`, `re_cor` and `formula`
belong to the components and are errors here; `cors`, `covariates`,
`seed`, `reps`, `scale`, `levels`, `correlations`, `cor_scale`, `smoke`,
`preflight` and `...` keep their meaning, and `prior` is `NULL` or a
list by component, as in
[`fit_components()`](https://www.gfrischkorn.org/bmmtools/reference/fit_components.md).
`subjects = "fixed"` is not supported for components yet. A component
cannot be named `sim`, `est`, `cor` or `sd`.

Grid columns: `n_subjects`; `n_trials_<comp>` for a component's number
of trials (a plain `n_trials` column is an error); `<comp>_<par>` or
`<comp>_<par>_<task_col><level>` for population values, a bare parameter
setting every task and a full term overriding it; `sd_<comp>_<par>` for
SDs; `cor_<a>__<b>` with prefixed terms or covariates, in either order.
They apply to numeric `pars` and `sds` of a component and to the result
of a function-valued one. A missing value (`NA`) leaves the component's
value, so a column may name a task that only some rows have. A column
whose component or term is unknown is an error; so is a column whose
prefix is no component but whose remainder is a term of one, such as a
misspelt component name. Other columns are left for `model` and `cors`
functions.

Files per cell: `cell-<row>-rep-<rep>-sim.rds` holds the set;
`cell-<row>-rep-<rep>-<comp>.rds` (with its `.key`) each fit; a sidecar
`cell-<row>-rep-<rep>-<comp>-est.rds` per component holds its estimates
and, with `"model"` requested, its model correlations; and
`cell-<row>-rep-<rep>-cor.rds` holds what needs every fit, the
correlations of
[`extract_correlations()`](https://www.gfrischkorn.org/bmmtools/reference/extract_correlations.md)
on the set of fits (covariates included) and the subject means, with the
components' keys. A resume in which every sidecar matches reads no fit;
otherwise the fits are obtained through
[`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md),
which reuses cached ones, and the stale sidecars are rewritten. The
preflight fits every component of the first cell into
`<dir>/preflight-<comp>`, and is skipped when each of them has a cached
fit or a usable sidecar.

The result is one `bmmtools_recovery` with prefixed terms (`a_kappa`,
`b_c_task1`); `converged` is that component fit's verdict, and
`scale = "natural"` uses every component's links. The `cells` attribute
has one row per cell and component, with a `component` column; the
`correlations` and `subject_means` attributes use prefixed terms, so
[`subject_table()`](https://www.gfrischkorn.org/bmmtools/reference/subject_table.md)
works on the result. A component whose fit fails is recorded as
`"error"` for that cell, the cell's other components are still scored,
its correlations and subject means are skipped, and the grid goes on.

Correlations between parameters of separately fitted components are
attenuated by the reliabilities of both estimates; see "Separate fits"
in
[`extract_correlations()`](https://www.gfrischkorn.org/bmmtools/reference/extract_correlations.md).
The `"model"` estimator works within one fit only and needs a component
`formula` with correlated terms.

## Examples

``` r
if (FALSE) { # \dontrun{
out <- recovery_grid(
  bmm::mixture2p(resp_error = "y"),
  grid = expand.grid(n_subjects = c(20, 50), n_trials = c(30, 100)),
  pars = c(kappa = log(8), thetat = qlogis(0.75)),
  sds = c(kappa = 0.3, thetat = 0.5),
  dir = "fits/mixture2p", reps = 10, seed = 1,
  chains = 4, iter = 2000, backend = "cmdstanr"
)
summary(out)
} # }
```
