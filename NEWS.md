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
