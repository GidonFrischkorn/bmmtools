# bmmtools (development version)

* `plot_design()` plots recovery across the designs of a grid: a metric
  against trials per subject on a log axis, one line per number of
  subjects and a panel per term (`type = "profile"`), for population
  means, SDs, subject parameters or model correlations; the share of
  fits passing the convergence gate and the checks the others failed,
  per design (`"convergence"`); and coverage pooled over every design
  per term and level (`"coverage"`). Seconds per fit and the pass rate
  are profiles too. Each error bar has one documented source: ±1 Monte
  Carlo standard error for bias, RMSE, interval width and contraction,
  a 95 % Wilson interval for rates, and `r_low` to `r_high` for `r`.

* `render_card()` renders the validation card of a `run_validation()`
  result to one self-contained HTML file: the headline, the run record
  and the plan, each stage's verdict beside its criterion, a section per
  stage, and the numbers that moved since an earlier run. Every number
  on the card is read off the stored result; the template types none.
  It needs the quarto package and the Quarto CLI. The recovery and
  individual-differences criteria pool each parameter's intervals over
  the designs and fail it by a one-sided binomial test at
  `coverage_alpha` (0.05) over the number of parameters; held design by
  design, a calibrated model failed them most of the time at 10
  replications. The `ezdm` entry checks prior predictive mean response
  times against a plausibility range of 0.15 to 5 s (#21).

* `run_validation()` runs the data requirements stage: the smallest
  number of trials at which every subject-varying parameter's
  person-level correlation with its generating value has a lower bound
  (`r_low`) of at least 0.80 at 50 subjects, read off the recovery grid
  by `smallest_design()`, with the read-off at 0.70 and 0.90 beside it.
  Both the bar and the number of subjects are criteria of the plan. Its
  coverage condition is the recovery stage's pooled test: held design by
  design, it changed or removed the read-off of a calibrated `ezdm`
  model in 24 % to 56 % of runs at 10 replications. The stage is not
  applicable to a model without a subject-varying parameter or a grid
  without two numbers of trials at that many subjects (#21, #17).

# bmmtools 0.3.0

## Breaking changes

* **A recovery saved by 0.2.0 no longer satisfies the recovery
  contract.** A `bmmtools_recovery` lacks `ci_low_50`, `ci_high_50` and
  `covered_50` (#1) and `post_mean_link`, `post_sd_link`, `z` and
  `contraction` (#16); a `bmmtools_cor_recovery` lacks the same columns
  but `contraction`. Printing or summarising one is an error that names
  the missing columns. The link-scale moments need the draws, so there
  is no in-place fix: run `recovery_grid()` again on the grid's
  directory, which re-extracts every cell from its cached fit without
  refitting, or call `recover()`, `recover_subjects()` or
  `recover_correlations()` on the fits again.
* **A summary saved by 0.2.0 still prints but lacks the columns 0.3.0
  adds**: `coverage_50`, `ci_width_50` and `calibration_intercept` (#1),
  the z-score and contraction columns (#16), the reliabilities, the
  Monte Carlo standard errors of the rates, `ci_level` and a grid's
  design columns (#17). Summarise the recovery again to get them.
* **`criterion` of `sdt_yn` means something else under bmm 1.4.1.** A
  simulation or grid stored under 0.2.0 with a development build of bmm
  holds `criterion` truths measured from the midpoint between the noise
  and signal means; bmm 1.4.1 estimates `criterion` from the noise mean,
  so scoring such truths against a new fit reads as bias. Simulate again
  with `bias` in `pars` to keep the old quantity (#30).
* SBC results for `sdt_yn`, `sdt_mafc` and `m3` with numbers of options
  from 0.2.0 were computed on fewer trials than `data` described (#26),
  and a fit cached without a `backend` is refitted once (#24); both are
  described under their changes below.
* `summary()` of a recovery refuses an argument it does not know instead
  of ignoring it (#17).

## Planning a study

* `simulate_recovery()`, `recovery_grid()`, `sbc()`, `prior_check()` and
  `fit_ml(method = "optim")` now work without a `generator` for two more
  bmm models (#8). The four-parameter `ezdm` simulates its response-time
  summaries per boundary from `zr`; a boundary reached fewer than twice
  leaves them `NA`, and bmm (from 1.4.1) keeps such rows. `sdt_ranking`
  takes its set size `m` from the model, or per row from a column of
  `trial_design` when `m` names one, and each row's number of trials
  from `n_trials` there. Its `sdratio` is the log of the SD ratio, so
  the adapter hands bmm's generator the ratio itself. `sdt_rating` and
  `sdt_cdp` still need a `generator` of your own: bmm does not yet
  export the map from their parameters to thresholds. A model version
  bmmtools does not know now gets no adapter rather than another
  version's, so a later version of `ezdm` asks for a `generator`.

* `sbc()` now reads an `NA` count in `m3` data with numbers of options
  as 0, as bmm does, and takes the row's number of trials as its total
  over the other counts. It used to stop with an error that the row's
  total was not a whole number of trials.

* A new article, "Recovering the EZ diffusion model on the standard
  grid", reports a recovery study of the three-parameter `ezdm` under
  bmm 1.4.1: 20 or 50 subjects, 30, 100 or 300 trials, 10 replications
  each. It covers population values, between-subject SDs, subject
  parameters and the drift-boundary correlation. It also says how many
  trials a person-level correlation of 0.70, 0.80 or 0.90 needs. The
  Monte Carlo standard errors, reliabilities and that read-off are
  computed by hand in the article until `summary()` reports them (#17),
  and a validation card will replace the article (#21).
* A validation of a bmm model is now one plan and one run.
  `validation_plan()` resolves everything a validation needs before
  anything is fitted: the generating values, a grid of designs with its
  replications, the sampler budget, the design point and prior of the
  simulation-based calibration, the prior-check layout, an optional
  reference implementation, the criteria, the seed, the owner and a
  version, with a content hash that changes only when what is planned
  changes. `validation_defaults()` holds the values a plan fills from;
  entries ship for `mixture2p` and `ezdm`, and a model without one needs
  `pars`. `run_validation()` runs the stages in order (run record, prior
  predictive check, SBC, recovery grid, computation, identifiability,
  individual differences, agreement) into one directory, each cached
  under its own folder. A rerun reads every stage whose inputs did not
  change, a stage that fails is recorded with its message while the
  others run, and stages that do not exist yet (data requirements,
  model recovery, posterior checks) are listed as not run with their
  issue. The result keeps every stage's own object and adds
  `card_numbers()`, every number a card prints with its Monte Carlo
  standard error, `verdict_table()`, each criterion in words with its
  value and verdict, and a headline that says "validated" only with the
  plan version beside it. `compare_validations()` flags the numbers that
  moved between two runs by more than their Monte Carlo error (#21, #9).

* bmmtools works with bmm 1.4.1, which offers `criterion` and `bias` of
  `sdt_yn()` and `sdt_rating()` as alternative parameterizations, of
  which a fit estimates one. It read both as free parameters, so
  `simulate_recovery()` demanded a value for `bias` and
  `recovery_formula()` wrote a formula bmm rejects. Now one member of
  each alternative group is free: the one `pars` names, or bmm's
  default, `criterion`. Naming both is an error, and the other member
  appears in no truth table, formula or score. `recovery_formula()`
  gains `pars`, whose names choose the member
  (`recovery_formula(sim$model, pars = sim$pars)`); `recovery_grid()`,
  the components and `fit_ml()` in a grid follow each cell's `pars`,
  `sbc()` and `fit_ml()` follow their formula (or `fit_ml()`'s `start`),
  and a formula that estimates the other member is refused before
  anything is fitted. The `sdt_yn` adapter and density pass `bias` to
  bmm when it is the chosen member (#30).
* bmmtools now declares `bmm (>= 1.4.1)` in Suggests, the first CRAN
  release with the signal-detection models. Note what bmm 1.4.1 changed
  for `sdt_yn`: `criterion` is measured from the noise mean, and the
  criterion measured from the midpoint between the noise and signal
  means, which development builds of bmm called `criterion`, is now
  `bias`. A script written against such a build that simulates
  `c(d = , criterion = )` keeps running and simulates a different
  quantity; write `bias` to keep the old meaning (#30).

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
* Estimates, recovery rows and correlation estimates gain
  `post_mean_link` and `post_sd_link`, the posterior mean and SD of the
  draws on the link scale. They stay on the link scale whatever `scale`
  a recovery is scored on, because neither is invariant under a
  nonlinear link; `estimate` stays the posterior median. Maximum-
  likelihood rows carry the point estimate and its standard error, and
  `point` correlations `NA`. They are the groundwork for posterior
  z-scores and contraction in `summary()`. A grid cell file without them
  is re-extracted from its cached fit on the next `recovery_grid()` run,
  without refitting; `collect_grid()` reads it with `NA` and says so
  (#16).
* Recovery rows gain `z`, the posterior z-score `(post_mean_link -
  truth) / post_sd_link`, and `contraction`, `1 - post_sd_link^2 /
  prior_sd_link^2`, both on the link scale; `summary()` reports
  `z_mean`, `z_sd` and the mean `contraction`, with `z_mean_mcse` and
  `z_sd_mcse`, and its print shows the first three after
  `coverage_50`. Together they tell a posterior that covers the truth
  because the data located it from one that covers it because it never
  left a prior that already held it. The prior SD comes from the new
  `prior_sd` argument of `recover()` and `recovery_grid()`: a table, a
  fit with `sample_prior = "only"`, or a one-set `prior_check()`;
  `recovery_grid(prior_sd = "fit")` fits the prior once per distinct
  prior set and caches it, and `"analytic"` uses the new
  `prior_sd_table()`, which gives the prior SD in closed form for
  normal, logistic and Student t priors, their halves on standard
  deviations, and the exponential prior bmm puts on them.
  `collect_grid()` reuses the table the grid scored with and takes
  `prior_sd` to replace it. Without a source, `contraction` is `NA`;
  subject rows and maximum-likelihood rows never have one. Correlation
  recoveries gain `z` and the z summary columns. The article "Recovery
  summary columns" defines them and shows, on a grid with one to 100
  trials per subject, intervals that cover the truth while the
  posterior is still mostly the prior (#16).
* `summary()` of a recovery reports how reliably a design measures
  individual differences, on subject rows: `reliability_true`, the
  squared correlation of the posterior means with the generating
  values, and `reliability_model`, which needs no truth and so is the
  number real data can report -- `var(mean) / (var(mean) + mean(sd^2))`
  for posterior means and `1 - mean(se^2) / var(estimate)` for
  maximum-likelihood rows. Both are computed per replication on the
  link scale and averaged; `reliability_gap` is their paired
  difference, which shows whether the model's own reliability can be
  trusted for that design. Each comes with a Monte Carlo standard error,
  and so now do the rates: `coverage_mcse` and `coverage_50_mcse`
  (binomial at the nominal mass, so a coverage of 100 % still has one;
  between replications on subject rows), `detected_mcse` and
  `sign_recovery_mcse`. The print shows both reliabilities when there
  are subject rows (#17).
* `summary()` of a recovery takes `precision`, a target width for the
  `ci_level` interval, and reports `precision_target`, the share of
  intervals no wider than it (`precision_met`) and that share's Monte
  Carlo standard error. Targets are a table of `term`, `level` and
  `width` (optionally per `condition`), or a named vector for
  population rows (`precision = c(kappa = 1)`), and are read on each
  row's own scale. On subject rows the share is over every subject
  interval -- the chance that one participant's interval is that narrow
  -- with an MCSE between replications. `summary()` now refuses an
  argument it does not know instead of ignoring it (#17).
* `smallest_design()` reads the smallest design of a grid that meets a
  target, such as `c(r_low = 0.8)` or `c(reliability_model = 0.8,
  precision_met = 0.8)`: the fewest trials (or any other grid column)
  at which the target is met, per number of subjects, and met at every
  larger value in the grid too, so noise that lets 30 trials pass and
  100 fail does not answer 30. A design meets a target when the lower
  end of its 95 % Monte Carlo interval clears it -- the Wilson bound
  for rates, so a perfect rate certifies 0.80 only from 16 replications
  on -- and the answer says when the grid's smallest design already
  met it, when the target was not reached, and when the point estimate
  meets it but the replications do not, with the number of replications
  a perfect rate would need. `coverage_check = TRUE` also requires
  population coverage within one MCSE of nominal. To support it, a
  grid summary now carries the grid's design columns after `condition`
  (`n_subjects`, `n_trials`, ...), and shows them in its print, and
  every summary a `ci_level` column (#17).
* The README and DESCRIPTION name both groups bmmtools serves: model
  developers and methodologists, and researchers planning a study. A new
  README section, "Planning a study", shows how `recovery_grid()` and
  `summary()` already answer how many subjects and trials a design needs
  (#14).
* A new article, "Recovering the EZ diffusion model on the standard
  grid", reports a recovery study of the three-parameter `ezdm` under
  bmm 1.4.1: 20 or 50 subjects, 30, 100 or 300 trials, 10 replications
  each. It covers population values, between-subject SDs, subject
  parameters and the drift-boundary correlation. It also says how many
  trials a person-level correlation of 0.70, 0.80 or 0.90 needs. The
  article computes its Monte Carlo standard errors, reliabilities and
  that read-off by hand, as it was written before `summary()` and
  `smallest_design()` reported them (#17); a validation card will
  replace it (#21).

## Score recovery

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
* `summary()` of a recovery object and `recovery_ccc()` gain
  `calibration_intercept`, after `calibration_slope`: the intercept of
  the generating value regressed on the estimate, so that the two
  columns give the whole calibration line. It is read on the scale of
  the row, and at subject level the replications' intercepts are
  averaged. It is `NA` wherever the slope is
  ([#1](https://github.com/GidonFrischkorn/bmmtools/issues/1)).
* Printing a recovery summary shows a core set of columns (the row, `n`,
  `bias`, `rmse`, both coverages, `r`, `ccc` and the calibration
  columns) and names how many it hides; `tibble::as_tibble()` prints all
  of them. The summary itself still has every column.
  Tracked in [#1](https://github.com/GidonFrischkorn/bmmtools/issues/1).
* Estimates, recovery rows and correlation estimates gain
  `post_mean_link` and `post_sd_link`, the posterior mean and SD of the
  draws on the link scale. They stay on the link scale whatever `scale`
  a recovery is scored on, because neither is invariant under a
  nonlinear link; `estimate` stays the posterior median. Maximum-
  likelihood rows carry the point estimate and its standard error, and
  `point` correlations `NA`. They are the groundwork for posterior
  z-scores and contraction in `summary()`. A grid cell file without them
  is re-extracted from its cached fit on the next `recovery_grid()` run,
  without refitting; `collect_grid()` reads it with `NA` and says so
  (#16).
* Recovery rows gain `z`, the posterior z-score `(post_mean_link -
  truth) / post_sd_link`, and `contraction`, `1 - post_sd_link^2 /
  prior_sd_link^2`, both on the link scale; `summary()` reports
  `z_mean`, `z_sd` and the mean `contraction`, with `z_mean_mcse` and
  `z_sd_mcse`, and its print shows the first three after
  `coverage_50`. Together they tell a posterior that covers the truth
  because the data located it from one that covers it because it never
  left a prior that already held it. The prior SD comes from the new
  `prior_sd` argument of `recover()` and `recovery_grid()`: a table, a
  fit with `sample_prior = "only"`, or a one-set `prior_check()`;
  `recovery_grid(prior_sd = "fit")` fits the prior once per distinct
  prior set and caches it, and `"analytic"` uses the new
  `prior_sd_table()`, which gives the prior SD in closed form for
  normal, logistic and Student t priors, their halves on standard
  deviations, and the exponential prior bmm puts on them.
  `collect_grid()` reuses the table the grid scored with and takes
  `prior_sd` to replace it. Without a source, `contraction` is `NA`;
  subject rows and maximum-likelihood rows never have one. Correlation
  recoveries gain `z` and the z summary columns. The article "Recovery
  summary columns" defines them and shows, on a grid with one to 100
  trials per subject, intervals that cover the truth while the
  posterior is still mostly the prior (#16).

## Simulate

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
* `coding = "contrast"` with the default contrasts
  (`stats::contr.treatment`), or any contrast matrix with column names,
  now recovers: the matrix's column names made brms call the contrast
  `task2` where the truth says `task1`, and scoring stopped with "No
  term in `truth` matches an estimated parameter". The column names are
  dropped, so the contrast is always `<task_col>1` … `<task_col>(k-1)`.
  Unnamed contrasts such as `bayestestR::contr.equalprior` were not
  affected ([#1](https://github.com/GidonFrischkorn/bmmtools/issues/1)).

## bmm 1.4.1

* bmmtools now declares `bmm (>= 1.4.1)` in Suggests, the first CRAN
  release with the signal-detection models. Note what bmm 1.4.1 changed
  for `sdt_yn`: `criterion` is measured from the noise mean, and the
  criterion measured from the midpoint between the noise and signal
  means, which development builds of bmm called `criterion`, is now
  `bias`. A script written against such a build that simulates
  `c(d = , criterion = )` keeps running and simulates a different
  quantity; write `bias` to keep the old meaning (#30).
* bmmtools works with bmm 1.4.1, which offers `criterion` and `bias` of
  `sdt_yn()` and `sdt_rating()` as alternative parameterizations, of
  which a fit estimates one. It read both as free parameters, so
  `simulate_recovery()` demanded a value for `bias` and
  `recovery_formula()` wrote a formula bmm rejects. Now one member of
  each alternative group is free: the one `pars` names, or bmm's
  default, `criterion`. Naming both is an error, and the other member
  appears in no truth table, formula or score. `recovery_formula()`
  gains `pars`, whose names choose the member
  (`recovery_formula(sim$model, pars = sim$pars)`); `recovery_grid()`,
  the components and `fit_ml()` in a grid follow each cell's `pars`,
  `sbc()` and `fit_ml()` follow their formula (or `fit_ml()`'s `start`),
  and a formula that estimates the other member is refused before
  anything is fitted. The `sdt_yn` adapter and density pass `bias` to
  bmm when it is the chosen member (#30).

## Running a study

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
* The help page of `recovery_grid()` gains a section "Backend", `?sbc`
  points to it, and "Running a recovery study" a paragraph: use cmdstanr
  for grids and SBC. Under rstan every fit compiles its model again
  (19–23 s per fit, measured for a small `mixture2p` grid); cmdstanr
  compiles once and reuses the executable, and
  `options(cmdstanr_write_stan_file_dir = )` keeps it across R sessions
  (#18).
* "Running a recovery study" gains a section on a combined design: one
  grid that recovers population intercepts, a contrast-coded task
  effect, between-subject SDs, correlated subject values and their
  correlations, reported as one table per estimand family, from 20
  sampled fits of `mixture2p`
  ([#1](https://github.com/GidonFrischkorn/bmmtools/issues/1)).

## Comparing implementations

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

## Validating a model (experimental)

* A validation of a bmm model is now one plan and one run.
  `validation_plan()` resolves everything a validation needs before
  anything is fitted: the generating values, a grid of designs with its
  replications, the sampler budget, the design point and prior of the
  simulation-based calibration, the prior-check layout, an optional
  reference implementation, the criteria, the seed, the owner and a
  version, with a content hash that changes only when what is planned
  changes. `validation_defaults()` holds the values a plan fills from;
  entries ship for `mixture2p` and `ezdm`, and a model without one needs
  `pars`. `run_validation()` runs the stages in order (run record, prior
  predictive check, SBC, recovery grid, computation, identifiability,
  individual differences, agreement) into one directory, each cached
  under its own folder. A rerun reads every stage whose inputs did not
  change, a stage that fails is recorded with its message while the
  others run, and stages that do not exist yet (data requirements,
  model recovery, posterior checks) are listed as not run with their
  issue. The result keeps every stage's own object and adds
  `card_numbers()`, every number a card prints with its Monte Carlo
  standard error, `verdict_table()`, each criterion in words with its
  value and verdict, and a headline that says "validated" only with the
  plan version beside it. `compare_validations()` flags the numbers that
  moved between two runs by more than their Monte Carlo error (#21, #9).
* The six validation functions above are experimental, and their help
  pages say so: no validation card has been made with them yet, so
  their arguments, the stages a run holds and the shape of its result
  can still change.

## Infrastructure and documentation

* The R-CMD-check workflow now checks bmmtools against both released
  bmm and bmm's `develop` branch (weekly on a schedule, and on every
  push and pull request), and installs SBC so `test-sbc.R`'s 47
  previously-skipped tests run on the runner. Signal-detection tests
  guard per model with the new `skip_if_no_bmm_model()`, replacing
  `skip_if_no_bmm_sdt()`, which only ever checked for `sdt_yn`. Tracked
  in [#7](https://github.com/GidonFrischkorn/bmmtools/issues/7).
* The README states what bmmtools reads off a fitted bmm model beyond
  its exported API (`$links`, `$fixed_parameters`, `$resp_vars`,
  `$other_vars`, and the generator adapters), and why a bmm parameter
  rename can break an adapter. Tracked in
  [#10](https://github.com/GidonFrischkorn/bmmtools/issues/10).

## Known limitations

* A grid cell is not reproducible from its `seed` alone. Two runs of the
  same `ezdm` grid with the same `seed` simulated identical data and
  passed the same sampler seed, but gave different chains: estimates
  differed by up to 0.036, and 2 of 360 population interval verdicts
  changed. The likely cause, not yet confirmed, is that bmm draws the
  initial values from R's random number generator, which `fit_cached()`
  does not seed. A cached fit is read back as it was; a refit is a new
  draw.

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
