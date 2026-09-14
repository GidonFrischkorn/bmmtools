# bmmtools 0.0.0.9000

* A package website, <https://gidonfrischkorn.github.io/bmmtools/>, with
  five articles: validating a new model end to end, running a recovery
  grid, prior predictive checks, validating a model without a built-in
  adapter, and the definitions of the recovery summary columns. The
  articles render from saved results, so the site builds without Stan.
* Two example objects ship with the package so that `summary()`,
  `plot_recovery()` and `plot_prior_check()` can be tried without
  compiling a model: `recovery_mixture2p`, a 2 × 2 grid with five
  replications, and `prior_check_sdt_yn`, comparing bmm's default priors
  with narrower between-subject SD priors. The reference examples now
  run on them.

* `prior_check()` shows what your priors say the data should look like.
  The draws come from `bmm(sample_prior = "only")` through
  `fit_cached()`, so the prior being checked is exactly the prior bmm
  fits with, and the result is on the observable scale: the share of
  prior-predictive draws at the floor and the ceiling of the response and
  the 50th, 90th and 95th percentiles, or any `summary(yrep, data)` you
  write. The floor and the ceiling are derived from the model where they
  are known (`0` and the trial count for the signal-detection models,
  plus or minus pi for the circular ones) and are `NA` where they are
  not, rather than invented. Pass a named list of priors to fit each one
  and compare them: `summary()` then gives one column per prior set and,
  with exactly two, their difference.
* `plot_prior_check()` draws the prior-predictive distribution, one thin
  line per draw with the observed response over it, as a histogram for a
  discrete response, or the summary statistics themselves.
* `simulate_recovery()` turns a bmm model and population values on the
  link scale into a data set and the truth that produced it: subject
  values are drawn on the link scale, converted with `inverse_link()`,
  and handed to the model's own `r<model>()` generator through a small
  adapter (shipped for `sdt_yn`, `sdt_mafc`, three-parameter `ezdm`,
  `ddm`, `mixture2p` and `sdm`) or to a generator you supply. Subject
  values can be passed in instead of drawn, so replications can share
  the same simulated people. `recovery_formula()` gives every free
  parameter a random intercept.
* `recovery_grid()` runs a design grid over subjects and trials with
  replications: one durable file per cell, a resume that reads those
  files, a smoke mode in its own directory, a preflight fit that catches
  a compile or init error before any cell, cells ordered so the first
  completed block spans the design, one seed per cell derived from the
  master seed, and scoring at both levels into one recovery object
  whose `condition` column names the grid row. `summary()` groups by
  condition when it is present.
* `check_convergence()` is the convergence gate: worst-case rhat and
  effective sample sizes, divergent transitions and tree-depth hits in
  one row, with a `pass` verdict under thresholds that are arguments
  (rhat at most 1.05, bulk ESS at least 400, at most ten divergences).
  Parameters with missing diagnostics are dropped rather than failed
  on, tree-depth hits are reported but not gated, and a fit with
  nothing assessable gets `NA`, not a pass.
* `extract_estimates()` carries the verdict as a `converged` column, and
  `summary()` of a recovery object reports `n_converged`, the number of
  replications whose fit passed. A hand-built estimates tibble still
  scores, with `NA` where no verdict exists.
* `fit_cached()` fits a model once and reuses the saved fit while
  nothing that determines it has changed. Its key covers the formula
  (deparsed, so the environment does not matter), data, model, prior,
  seed, chains, iterations, warmup, thinning, `control`, `init` (by its
  text when it is a function), backend, the installed bmm and brms
  versions and the Stan toolchain of the backend in use, which is more
  than brms's own `file` cache compares. A sidecar `.key` file records
  each component's hash, so a refit says which component changed; the
  fit and the key are written atomically.
* `recover()` and `recover_subjects()` score fits against known
  generating values and return a `bmmtools_recovery` object: one row per
  fit and parameter, with the estimate, its interval, the generating
  value, the signed error and whether the interval covered. Both take a
  `brmsfit`, a list of them, or an estimates tibble, so a scoring
  pipeline runs with no fitting package installed.
* `summary()` of a recovery object gives the per-parameter metrics: bias,
  RMSE, coverage, mean interval width, the Pearson correlation with a
  Fisher-z interval, the Spearman correlation, and Lin's concordance with
  its scale and location components. Where fewer than three complete
  pairs or no spread makes a correlation undefined it is `NA`, never `0`.
  Subject-level correlations are computed within replication and combined
  on Fisher's z scale.
* `extract_estimates()` turns a fit into the estimates tibble, the
  contract between the run layer and the score layer. Subject-level rows
  are the per-draw sum of the population intercept and the group-level
  deviation; parameters the model fixed are identified by zero posterior
  variance, not by a missing rhat, and dropped by default.
* Scoring on the natural scale handles the two links whose inverse is
  not monotone across zero: a `sqrt`-link interval that spans zero on
  the link scale becomes `[0, max]`, and an `inverse`-link interval that
  spans zero has no natural-scale image, so its bounds and `covered` are
  `NA` and a warning names the term.
* A recovery object survives the `dplyr` verbs that keep its contract
  columns and becomes a plain tibble under one that drops them, so a
  reduced object prints as a tibble instead of failing inside a metric.
* `plot_recovery()` plots estimates against generating values, one panel
  per parameter, with interval bars and an identity line.
* `inverse_link()` transforms a vector from the link scale to the natural
  scale for the twelve links bmm uses.
* Package skeleton, Milestone 0.
