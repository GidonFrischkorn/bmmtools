# bmmtools (development version)

* A benchmark module compares implementations of a model on the same
  data: `benchmark()` runs them interleaved (A, B, A, B) and returns one
  row per run, `benchmark_metrics()` measures a single fit, and
  `benchmark_compile()` times a forced compile on its own. Compile time
  and sampling time are kept apart, and every rate divides by seconds
  summed over chains rather than by wall time, which would make the same
  fit look cheaper on more cores: a run reports `us_per_gradient` (the
  median over chains of chain sampling time over that chain's leapfrog
  steps) and ESS per chain-second for bulk and tail, over the
  post-warmup phase and over warmup plus sampling. `summary()` reports
  the spread beside the centre per implementation, and the ratio between
  implementations as a geometric mean of per-pair ratios with a t
  interval on the logs -- never a bare quotient. No benchmark fit is
  cached, so none can be served to a recovery grid (#4).
* The README and DESCRIPTION name both groups bmmtools serves: model
  developers and methodologists, and researchers planning a study. A new
  README section, "Planning a study", shows how `recovery_grid()` and
  `summary()` already answer how many subjects and trials a design needs
  (#14).
* `simulate_recovery()`, `recovery_grid()` and `recovery_component()`
  gain `trial_design`, the per-trial variables a generator needs and the
  model does not hold, such as set sizes and non-target locations: a data
  frame shared by every subject, a data frame with an `id` column giving
  each subject its own trials, or a function of `n_trials` called once
  per subject under the seed. Each generator call receives its rows as a
  fourth argument, `trial_design`, and the same rows are bound into the
  data. Without it nothing changes: a generator written for 0.2.0 is
  called as before, and a seeded simulation gives the same data and the
  same cache key. `sbc()` takes the design from its `data` for a model
  whose built-in generator reads one (#15).
* `simulate_recovery()` has built-in generators for `mixture3p` and for
  `imm` in all three versions (`full`, `bsc`, `abc`). Both read the set
  size, the non-target locations and, for `imm`, the distances per trial
  from `trial_design`, and draw each trial with the weights bmm's
  likelihood gives it. `prior_check()` knows their response range.
  `mixture3p`'s `thetat` and `thetant` have bmm's `softmax` link, which
  has no inverse for one term: they reach a generator on the link scale,
  and `recover()`, `recover_subjects()`, `cross_check()`,
  `subject_table()` and the correlation scorers keep them there under
  `scale = "natural"` with a message, where they used to stop with an
  error; `print()` and `plot_recovery()` name those terms (#15).
* `simulate_recovery()` has a built-in generator for `m3` in all three
  versions (`ss`, `cs`, `custom`). With numbers of options on the model,
  a subject has one row of category counts over `n_trials` trials; with
  `num_options` naming columns, those columns come from `trial_design`
  and each row is one trial. The custom version's activation formulas
  exist only in the formula it is fitted with, so `simulate_recovery()`
  gains `formula`, read by that generator alone; `recovery_grid()`,
  `recovery_component()` and `sbc()` pass on the formula they fit with.
  A generator you write never receives it and is called as before.
  Categories land in the columns `resp_cats` names however the
  activations or the numbers of options are ordered. `prior_check()`
  takes each row's total as the ceiling of every category's count (#15).
* `sbc()` simulates a count model (`sdt_yn`, `sdt_mafc`, `ezdm`, and
  `m3` with numbers of options) with the trials each row of `data`
  holds. Before, its default generator gave every row as many trials as
  `data` had rows per subject: 2 for `sdt_yn` and 1 for `sdt_mafc` and
  `m3`, so **SBC results for these three models from 0.2.0 were computed
  on a smaller design than `data` described**, and `ezdm` could not run,
  with an error that blamed the prior. Rows may now hold different
  numbers of trials, such as 100 signal and 50 noise trials, and a count
  that cannot be simulated is an error before the prior is fitted. The
  per-trial models are unaffected: their simulated data sets are the
  same as in 0.2.0. For this, `simulate_recovery()`, `recovery_grid()`
  and `recovery_component()` take an optional `trial_design` for these
  four models, one design row per row of data with its own number of
  trials (and, for `sdt_yn`, its stimulus); without one, nothing changes
  (#26).
* `fit_cached()` keys a fit run without a `backend` on the Stan toolchain
  bmm actually uses: the `brms.backend` option if set, otherwise cmdstanr
  if it is installed, otherwise rstan. Before, it assumed rstan, so a fit
  bmm ran with CmdStan was keyed on the rstan version: a CmdStan upgrade
  did not invalidate it, and an rstan upgrade refitted it. **Every cache
  written without a `backend`, without the `brms.backend` option, and
  with cmdstanr installed is refitted once**, with a message naming the
  `toolchain` component; this includes `recovery_grid()` cells and the
  prior fit of `sbc()` and `prior_check()`. A grid cell's stored
  estimates carry its key too, so a grid whose fit files were deleted
  after it ran is refitted, not re-read. Calls that pass `backend`, or
  set the option, as every article and `data-raw/` script does, keep
  their keys (#24).
* Every estimate now carries a central 50 % interval, `ci_low_50` and
  `ci_high_50`, beside the `ci_level` interval: the 25th and 75th
  percentiles of the draws, a Wald interval at `qnorm(0.75)` standard
  errors on `fit_ml(method = "optim")`, and a Fisher-z interval for the
  `point` correlation estimator. Recovery rows gain `covered_50`, and
  `summary()` of a recovery or a correlation recovery gains
  `coverage_50` and `ci_width_50`: the 95 % interval checks the tails of
  the posterior, the 50 % interval its centre. The mass is fixed and
  carried in the name, so it does not follow `ci_level`. An estimates
  tibble without the inner bounds still scores, with `coverage_50` `NA`.
  `cross_check()` keeps one interval.
* `recovery_grid()` stores the 50 % interval in each cell's file. A cell
  file written before it is re-extracted from its cached fit on the next
  run, without refitting; `collect_grid()` reads such files with the
  inner columns `NA` and says how many there were.
* Printing a recovery summary shows a core set of columns (the row, `n`,
  `bias`, `rmse`, both coverages, `r`, `ccc` and the calibration
  columns) and names how many it hides; `tibble::as_tibble()` prints all
  of them. The summary itself still has every column.
  Tracked in [#1](https://github.com/GidonFrischkorn/bmmtools/issues/1).
* `summary()` of a recovery object and `recovery_ccc()` gain
  `calibration_intercept`, after `calibration_slope`: the intercept of
  the generating value regressed on the estimate, so that the two
  columns give the whole calibration line. It is read on the scale of
  the row, and at subject level the replications' intercepts are
  averaged. It is `NA` wherever the slope is
  ([#1](https://github.com/GidonFrischkorn/bmmtools/issues/1)).
* `coding = "contrast"` with the default contrasts
  (`stats::contr.treatment`), or any contrast matrix with column names,
  now recovers: the matrix's column names made brms call the contrast
  `task2` where the truth says `task1`, and scoring stopped with "No
  term in `truth` matches an estimated parameter". The column names are
  dropped, so the contrast is always `<task_col>1` … `<task_col>(k-1)`.
  Unnamed contrasts such as `bayestestR::contr.equalprior` were not
  affected ([#1](https://github.com/GidonFrischkorn/bmmtools/issues/1)).
* "Running a recovery study" gains a section on a combined design: one
  grid that recovers population intercepts, a contrast-coded task
  effect, between-subject SDs, correlated subject values and their
  correlations, reported as one table per estimand family, from 20
  sampled fits of `mixture2p`
  ([#1](https://github.com/GidonFrischkorn/bmmtools/issues/1)).
* The README states what bmmtools reads off a fitted bmm model beyond
  its exported API (`$links`, `$fixed_parameters`, `$resp_vars`,
  `$other_vars`, and the generator adapters), and why a bmm parameter
  rename can break an adapter. Tracked in
  [#10](https://github.com/GidonFrischkorn/bmmtools/issues/10).
* The R-CMD-check workflow now checks bmmtools against both released
  bmm and bmm's `develop` branch (weekly on a schedule, and on every
  push and pull request), and installs SBC so `test-sbc.R`'s 47
  previously-skipped tests run on the runner. Signal-detection tests
  guard per model with the new `skip_if_no_bmm_model()`, replacing
  `skip_if_no_bmm_sdt()`, which only ever checked for `sdt_yn`. Tracked
  in [#7](https://github.com/GidonFrischkorn/bmmtools/issues/7).

# bmmtools 0.2.0

## Simulate

* `simulate_recovery(coding = "contrast")` and
  `recovery_grid(coding = "contrast")` generate a task design as an
  intercept and a contrast rather than as cell means, with `contrasts`
  taking any contrast matrix --- a design's own, or one of
  `bayestestR`'s equal-prior codings. The truth is transformed with the
  design, so the intercept and the effect have true values, true SDs and
  a true correlation matrix of their own, and `recovery_formula(coding =
  "contrast")` writes the matching formula. `contrasts` is deliberately
  not an argument of `recovery_formula()`: a factor's contrast matrix
  reaches brms through the data, not through the formula.
* `"effect"` joins the `level` vocabulary of `extract_estimates()`,
  `recover()` and `recovery_grid()`, so a contrast coefficient can be
  scored on its own, away from the intercept it shares a parameter with.
  It is kept on the link scale, as `"sd"` is.
* A generator and a density adapter for bmm's censored-shifted Wald
  (`cswald`), both the `simple` and the `crisk` version, bringing the
  adapter table to seven models. `bmm::rcswald()` takes the boundary
  separation of the diffusion it draws from, while the `simple`
  version's own `bound` is the distance from an unbiased start to one
  boundary, so the adapter converts between them: generated data
  profiles back to the model's `bound`, not to twice it.

## Running a study

* `run_info()` records the machine and the toolchain a study ran on ---
  R version, platform, host, cores, CPU model, memory, the installed
  versions of bmm, brms, cmdstanr, posterior and bmmtools, and the
  CmdStan version --- as one row, so a runtime table can say what the
  times it reports were measured on. It never errors: a field whose
  source is missing on the platform is `NA`.

* `collect_grid()` rebuilds what `recovery_grid()` returned from the
  files it left behind, without fitting anything and without reading a
  fit. The per-cell sidecars hold everything that was scored, so a study
  can run on a server and be assembled on a laptop that has neither the
  fits nor the versions of bmm, brms and Stan that produced them. It is
  a reader, not a resume: it checks neither the cache key nor the
  version, and a resume is still rerunning the identical call in the
  same directory. `scale`, `levels` and `correlations` re-score the
  stored rows, so a grid run on the link scale can be read again on the
  natural one.
* `recovery_grid(convergence = )` gives `check_convergence()` its
  thresholds, checked before any cell runs. The gate is computed once
  per cell and its whole diagnostics row is stored in the cell's
  sidecar, so `attr(x, "cells")` now carries `max_rhat`,
  `min_ess_bulk`, `min_ess_tail`, `n_divergent`, `n_max_treedepth`,
  `n_variables` and `failed`, together with `fit_seconds`, `chains`,
  `iter` and `threads` --- all of it without reading a fit. `elapsed`
  keeps its meaning: how long the cell took in this run, read time
  included, while `fit_seconds` is how long the fit itself took.
* `recovery_grid()` writes `<dir>/grid.rds`, the record of what the grid
  was called with, which is what `collect_grid()` reads. A directory
  reused for a different design has its record rewritten, with a message
  naming what differs.
* `fit_cached()` records how long the fit took. The time is written to
  `<file>.meta.rds` beside the fit and comes back in the
  `bmmtools_cache` attribute as `seconds`, including when the call
  reused a cached fit, so a resumed study still reports the time each
  fit once took rather than the time it took to read it. The time is
  deliberately not part of the cache key: a key component that changed
  with every run would differ from every stored key and refit
  everything. A fit cached before this release has no meta file and
  reports `NA`.

## Comparing estimators

* `fit_ml()` estimates a model subject by subject with no pooling and
  returns the estimates in the same tibble shape a hierarchical fit
  gives, so `recover_subjects()` can score both against one truth. It
  optimises bmm's own generated likelihood through
  `algorithm = "laplace"` and re-implements no model. Flat priors are the
  default, so the mode is a maximum-likelihood estimate rather than a
  penalised one; `prior = "default"` keeps bmm's priors. It exists to
  measure what hierarchical estimation buys, not to recommend maximum
  likelihood for inference.
* `fit_ml(method = "optim")` is a second route to the same estimates,
  maximising bmm's R density for the model directly with `optim()` and
  taking its interval from the Hessian (`ci_method = "wald"`). It needs
  no compiler and no Stan, which makes it the route that runs anywhere,
  and it optimises each subject independently, so one subject at a
  boundary cannot stall the rest. It is capped at the models bmmtools
  carries a density for --- the same seven the simulation adapters cover
  --- and `nll` supplies one for anything else. Measured against the
  Stan route on 40 subjects of `mixture2p`: the same estimates to Monte
  Carlo error, standard errors agreeing to a mean ratio of 1.000, and
  the same `coverage` and `calibration_slope` to three decimals, in 0.35
  s against 8 s.
* With `method = "stan"`, the default, `fit_ml()` runs **one** fit for
  the whole data set rather than one per subject. Under `p ~ 0 + id` the
  log posterior is a sum of per-subject terms, so the joint mode is the
  vector of per-subject modes; this was checked against independent
  optimisation of bmm's own density and agreed to 3e-06 on the link
  scale.
* A subject whose estimate leaves a finite range on the link scale, or
  whose optimiser reports failure on the `optim` route, is reported with
  `estimate = NA` and `converged = FALSE`, never dropped:
  the recovery metrics drop incomplete pairs silently, so dropping would
  score maximum likelihood on the easiest subjects and the hierarchical
  fit on all of them. `attr(x, "ml_cells")` and `print()` carry the
  per-subject counts.
* `recovery_grid(ml = TRUE)` fits every cell subject by subject as well,
  on the same simulated data, and scores both estimators against one
  truth at the subject level, so the shrinkage comparison can be read
  over a whole design rather than one data set. A Stan-route ML fit is
  cached beside the cell's other files; either route's rows go into the
  sidecar with the request that made them, and a resume reads them back
  without refitting. A cell whose ML fit
  fails keeps its hierarchical rows and is named in a warning, with the
  per-cell record in `attr(x, "ml_cells")`. `ml = list(...)` passes
  arguments to `fit_ml()`, so `ml = list(method = "optim")` runs the
  whole design without a compiler. The sidecar records which route made
  its rows and a resume that asks for the other one refits.

* A new article, "Hierarchical estimation against subject-wise maximum
  likelihood", runs the comparison end to end on `mixture2p`: one data
  set fitted both ways and scored against one truth, then the same
  contrast over a grid of trial counts. It is where the estimator table
  of "Recovery summary columns" stops being a conjugate toy.
* Recovery objects carry an `estimator` column, and `summary()` groups by
  it. Two estimators of the same parameter — a hierarchical posterior and
  a subject-wise maximum-likelihood fit, say — can be bound together and
  scored against one truth without being pooled into a single bias and
  RMSE. `extract_estimates()` gains an `estimator` argument (default
  `"bayes"`), and a hand-built estimates tibble may carry the column
  itself.
* `plot_recovery(color_by = "estimator", annotate = TRUE)` labels each
  panel once per estimator instead of refusing to annotate it.
* `recover_subjects()` warns when two estimators were not scored on the
  same subjects, because the metrics drop incomplete pairs in silence and
  the rows would otherwise not be comparable.

### Breaking

* `summary()` on a recovery object now returns its rows in the order the
  estimates arrived rather than sorted by term, so that adding a second
  estimator does not reshuffle the rows of the first. A summary of one
  estimator holds the same numbers as before, in a different order.
* `estimator` is a required column of the recovery contract, so a
  `bmmtools_recovery` object saved by 0.1.0 no longer satisfies it.
  Restore one with `dplyr::mutate(old, estimator = "bayes")`. Objects
  built by `recover()`, `recover_subjects()` and `recovery_grid()` are
  filled automatically and are unaffected.

## Score recovery

* `summary()` reports `mae`, the mean absolute error, beside `rmse`. The
  two differ in how much one badly recovered subject moves them, which
  is worth seeing rather than inferring.
* Population and effect rows gain `detected`, the share of intervals
  that exclude zero, and `sign_recovery`, the share that exclude zero
  *and* fall on the side of it the truth is on. At a true value of zero
  `sign_recovery` is `NA`, because zero has no sign, and `detected` is
  then a false-positive rate rather than power --- a page reporting them
  has to say which of the two a row is. Subject rows carry both columns
  as `NA`.

# bmmtools 0.1.0

First release. bmmtools validates cognitive measurement models fitted
with bmm: one engine — simulate, fit, score — and a scorer per question.

## Simulate

* `simulate_recovery()` turns a bmm model and population values on the
  link scale into a data set and the truth that produced it. Subject
  values can be correlated (`cors`) and drawn together with observed
  covariates (`covariates`); `pars`, `sds` and `cors` may be functions
  evaluated under the seed. `tasks` and `task_col` give every subject
  every task, so each free parameter has one value per task under the
  term `kappa_task1`, brms's name for the cell mean. The truth lists the
  population values, the between-subject SDs, every correlation pair,
  the subject values and the covariates.
* `cors_from_factors()` builds a correlation matrix from factor loadings.
* `recovery_formula()` writes the formula a recovery study needs: a
  random intercept over `id` for every free parameter, correlated across
  parameters with `re_cor = "all"`, within each parameter's tasks with
  `re_cor = "within"`. With `task_col` it writes cell-means formulas such
  as `kappa ~ 0 + task + (0 + task || id)`.
* `recovery_component()` and `simulate_components()` simulate several
  models for the same simulated people, so parameters of different models
  can be correlated. The subject values of all components and the
  covariates are one multivariate normal draw, and terms are prefixed
  with the component name.
* Six generator adapters ship — `sdt_yn`, `sdt_mafc`, `ezdm`, `ddm`,
  `mixture2p` and `sdm` — and a model without one works the same way once
  you supply `generator`.

## Fit and run a study

* `fit_cached()` fits once and refits only when something that determines
  the fit changes. The key covers the formula, data, model, prior, seed,
  sampler settings, backend and the installed versions of bmm, brms and
  the Stan toolchain, and the message names what changed.
* `fit_components()` fits each component of a set on its own, with an
  optional prior per component.
* `check_convergence()` summarises rhat, bulk and tail ESS, divergent
  transitions and tree-depth hits into one row and a pass verdict. A fit
  that fails is still scored, and the verdict travels with the estimates
  as the `converged` column.
* `recovery_grid()` runs the loop over a design grid of subjects, trials
  and replications: one durable file per cell, a resume that needs no
  temporary state, a smoke mode, a preflight fit that catches a compile
  error before any cell runs, and cells ordered so the first completed
  block spans the design. Each cell writes a sidecar holding everything
  scored from the fit, so a resume can read the sidecars and the fits can
  be deleted. Grid columns set population values, SDs, task terms and
  correlations per row, and `model`, `formula`, `pars`, `sds` and `cors`
  may be functions of the row. A list of components as `model` simulates
  one set per cell and fits every component separately.

## Score recovery

* `extract_estimates()` returns the population- and subject-level
  estimates of a fit as a tibble; `level = "sd"` and `level = "cor"` add
  the between-subject standard deviations and correlations on the link
  scale, with correlations named `kappa__thetat` whichever order brms
  used.
* `extract_subject_draws()` returns the per-draw subject values as an
  array of iterations, chains, subjects and parameters.
* `recover()` and `recover_subjects()` score estimates against the
  generating values, on the natural scale by default.
  `recover(level = c("population", "sd"))` also scores the SDs, always on
  the link scale.
* `summary()` of a recovery reports bias, RMSE, credible-interval
  coverage, interval width and correlations, including Lin's concordance
  with a 95% interval, its accuracy factor, the scale and location
  shifts, a calibration slope and the spread of the generating values. At
  subject level the concordance is pooled across replications on Lin's Z
  scale; the geometric-mean pooling of the scale shift and slope has not
  been checked by simulation.
* `extract_correlations()` estimates between-subject correlations three
  ways — the model's own group-level correlation, the per-draw
  correlation across subjects, and the correlation of posterior means —
  and `recover_correlations()` scores them against both the generating
  correlation and the one the simulated subjects actually had. Its
  `summary()` adds the rate at which intervals exclude zero, as a
  false-positive rate or as power. Across separate fits these
  correlations are attenuated by the reliabilities of both estimates, and
  the documentation says which routes are not.
* `subject_table()` puts true and estimated subject values side by side,
  one row per subject, for a simulation and its fit or for a whole grid.
  It is the input for a structural equation model of true against
  estimated values.
* `recovery_ccc()` computes the same concordance columns for any pair of
  vectors.
* `inverse_link()` transforms values from the link scale for the twelve
  links bmm uses.

## Prior predictive checks

* `prior_check()` samples from the prior only and summarises the
  prior-predictive draws on the scale of the response: floor and ceiling
  rates, a quantile profile, or a statistic you write yourself. Several
  prior sets can be compared in one call.

## Correctness

* `sbc()` runs simulation-based calibration for a bmm model: it fits the
  prior once, turns its draws into data sets, fits each of them, and
  returns the SBC package's own `SBC_results`, so `SBC::plot_rank_hist()`
  and the rest apply unchanged. `level` picks which draws are ranked —
  the population parameters, the between-subject SDs, their correlations,
  and with `level = "subject"` the subject-level draws, one per subject
  and varying parameter. Whatever is ranked, everything the formula
  implies is always drawn from the prior and simulated from, because the
  ranks are uniform only when the data come from the joint prior.
* `sbc(generator =)` lifts the restriction on `formula`. The function
  receives every `b_`, `sd_`, `cor_` and `r_` draw of the prior fit as
  one named vector, under the fit's own draw names, and returns the data
  frame to fit, so any formula works and the truth names equal the draw
  names by construction.
* `sbc()` says what went wrong. When a prior draw is one the model cannot
  generate from, the error names the simulation and what that draw held.
  When any fit exceeds bmmtools' rhat 1.05 or falls below SBC's rank-ESS
  0.5, one warning reports both counts and names SBC's stricter rhat.
* `cross_check()` compares a fit's estimates with a reference — a closed
  form such as `bmm::sdt_d()`, another implementation, or published
  values — and returns `bias`, whether the fit's interval covers the
  reference and whether the two intervals overlap, in the same tibble
  shape the other scorers return. The reference is a comparison, not a
  truth, and the documentation says so.

## Plots

* `plot_recovery()` plots estimates against generating values, one panel
  per parameter, with credible intervals and a line at equality, and with
  `annotate = TRUE` labels each panel with r and the concordance. It is a
  generic and also plots a correlation recovery, coloured by estimator,
  and a cross-check, with the reference on x.
* `plot_prior_check()` plots prior-predictive draws against the observed
  data.

## Data, documentation and infrastructure

* `recovery_mixture2p` and `prior_check_sdt_yn` are example results, so
  the examples run without Stan.
* A pkgdown website with six articles:
  <https://www.gfrischkorn.org/bmmtools/>.
* A hex logo, `man/figures/logo.png`.
* DESCRIPTION carries `Remotes: hyunjimoon/SBC`, so `pak` and `remotes`
  resolve the SBC entry in Suggests; SBC is not on CRAN.
* Licensed under GPL (>= 2), compatible with bmm's GPL-2.
* The signal-detection adapters resolve bmm's `rsdt_yn()`, `rsdt_mafc()`,
  `dsdt_yn()` and `dsdt_mafc()` by name at call time. No released bmm
  exports them, so a literal `bmm::rsdt_yn()` made `R CMD check` report a
  missing object wherever a released bmm was installed. Behaviour is
  unchanged where the functions exist, and an absent one is now a named
  error rather than R's bare "not an exported object".

## Known limitations

* bmm's default prior on the group-level SD of a log-link parameter such
  as `kappa` is wide enough that `bmm::rmixture2p()` can fail on draws
  from it, so `sbc()` on `mixture2p` needs a prior given explicitly. The
  walkthrough article shows one.
* `cross_check()` takes one fit and compares at the population level;
  `"sd"` and `"cor"` references, and a list of fits, are not supported.
* `subjects = "fixed"` is an error for component sets.
