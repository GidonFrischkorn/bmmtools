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

## Known limitations

* bmm's default prior on the group-level SD of a log-link parameter such
  as `kappa` is wide enough that `bmm::rmixture2p()` can fail on draws
  from it, so `sbc()` on `mixture2p` needs a prior given explicitly. The
  walkthrough article shows one.
* `cross_check()` takes one fit and compares at the population level;
  `"sd"` and `"cor"` references, and a list of fits, are not supported.
* `subjects = "fixed"` is an error for component sets.
