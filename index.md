# bmmtools

Writing a new cognitive measurement model has two halves. bmm builds and
fits it:
[`bmm::use_model_template()`](https://venpopov.com/bmm/reference/use_model_template.html)
writes the skeleton of a new model, and
[`bmm::bmm()`](https://venpopov.com/bmm/reference/bmm.html) fits it.
bmmtools is the other half: it checks that the model works before anyone
trusts it on real data.

Every check is the same loop, simulate data from known values, fit the
model, score the fit, and only the scorer changes:

| Scorer | truth comes from | scored against | answers |
|----|----|----|----|
| [`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md) | a fixed grid of generating values | bias, RMSE, how often the credible interval covers the generating value, and the correlation between generating value and estimate with a confidence interval | can this model be estimated from data of this size |
| [`prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/prior_check.md) | the prior | the observable scale: floor and ceiling rates, quantile profile | are these priors sane on the scale a reader understands |
| [`recover_correlations()`](https://www.gfrischkorn.org/bmmtools/reference/recover_correlations.md) | correlated generating values, or a factor model | between-subject correlations, estimated three ways and scored against both the correlation used to generate and the one the simulated people actually had | can this model measure individual differences and their correlates |
| [`sbc()`](https://www.gfrischkorn.org/bmmtools/reference/sbc.md) | the prior | posterior rank of the truth, through the SBC package | is the implementation correct |
| [`cross_check()`](https://www.gfrischkorn.org/bmmtools/reference/cross_check.md) | a fixed grid, or real data | a closed-form estimator, another implementation, or published values | does the new model agree with what is already known |

[`sbc()`](https://www.gfrischkorn.org/bmmtools/reference/sbc.md) is
simulation-based calibration. Draw parameters from the prior, simulate
data from them, fit, and record where the drawn value falls among the
posterior draws. Averaged over the prior the posterior is the prior
again, so those ranks have to be uniform; a histogram with a peak or a
slope means the implementation and the model on paper are not the same
thing.

[`recover_subjects()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
scores person-level posterior parameters, which no other package does.
Scoring dispatches on `brmsfit`, so a plain brms fit can be scored too;
only the simulation step, which reads a bmm model object to know what to
generate, needs bmm itself.

Version 0.2.0 adds what running a study needs:
[`collect_grid()`](https://www.gfrischkorn.org/bmmtools/reference/collect_grid.md),
[`run_info()`](https://www.gfrischkorn.org/bmmtools/reference/run_info.md),
a convergence gate and a per-cell diagnostics record on
[`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md),
contrast-coded effect designs, and
[`fit_ml()`](https://www.gfrischkorn.org/bmmtools/reference/fit_ml.md)
for comparing hierarchical estimation with subject-wise maximum
likelihood. Every function this README names is built and documented:
the five scorers of the table,
[`recover_subjects()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md),
and underneath them
[`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md),
[`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md)
and
[`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md).
It is an early release even so. bmmtools is not on CRAN, and the
lifecycle badge says experimental because argument names and return
shapes can still change before 1.0.

## Who it is for

bmmtools serves two groups.

- **Model developers and methodologists**, who check that a model is
  implemented correctly and how well its parameters can be estimated:
  the author of a new bmm model, or a methodologist running a parameter
  recovery study.
- **Researchers planning a study**, who want to know how many subjects
  and trials their design needs to estimate a parameter or an effect
  with a given precision, before they collect any data.

Both questions are answered the same way: simulate data from known
values, fit the model with
[`bmm::bmm()`](https://venpopov.com/bmm/reference/bmm.html), and score
the fit against the values that generated it.

## A recovery check in brief

The model below is bmm’s two-parameter mixture model for continuous
reproduction, the task where someone reports a remembered colour or
orientation on a continuous scale. `kappa` is the precision of the
memory responses and `thetat` the probability that a response comes from
memory rather than a guess. Both are given on the link scale, `kappa`
logged and `thetat` on the log-odds scale, because that is the scale the
model is estimated on.

[`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md)
returns two things: `data`, ready to fit, and `truth`, which holds the
generating values. `truth$population` has the values above,
`truth$subjects` each simulated person’s own.

``` r

library(bmmtools)
model <- bmm::mixture2p(resp_error = "y")

sim <- simulate_recovery(
  model,
  pars = c(kappa = log(8), thetat = qlogis(0.75)),
  sds = c(kappa = 0.3, thetat = 0.5),
  n_subjects = 30, n_trials = 60, seed = 1
)
fit <- fit_cached(
  recovery_formula(model), sim$data, model,
  file = "fits/mixture2p", seed = 1
)
check_convergence(fit)

recovery <- recover(fit, sim$truth$population)
summary(recovery)
recover_subjects(fit, sim$truth$subjects) |> plot_recovery()
```

Line by line:
[`recovery_formula()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_formula.md)
writes the bmm formula a recovery study needs, a random intercept over
`id` for every free parameter.
[`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md)
fits with [`bmm::bmm()`](https://venpopov.com/bmm/reference/bmm.html)
and saves the result, so calling it again with everything unchanged
reads the fit back instead of sampling for an hour.
[`check_convergence()`](https://www.gfrischkorn.org/bmmtools/reference/check_convergence.md)
returns one row of diagnostics, the worst rhat, the smallest effective
sample sizes and the divergent transitions, together with a pass
verdict.
[`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
scores the population values and
[`recover_subjects()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
each simulated person’s, both returning a tibble that
[`summary()`](https://rdrr.io/r/base/summary.html) reduces to one row
per parameter and
[`plot_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/plot_recovery.md)
draws as estimate against generating value.

[`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md)
runs the same loop over a grid of subjects, trials and replications,
writing one file per cell so an interrupted run can be resumed, and
[`prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/prior_check.md)
shows what the priors imply for the data before any model is fitted to
them. The articles on the [package
website](https://www.gfrischkorn.org/bmmtools/) walk through each step.

## Planning a study

There is no dedicated design function yet, but
[`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md)
already answers how precise the estimates of a planned design will be.
Give it a grid over `n_subjects` and `n_trials`, the parameter values
and effects you expect, on the link scale, and a number of replications;
with `coding = "contrast"` it also scores a condition effect.
[`summary()`](https://rdrr.io/r/base/summary.html) then returns, per
design cell, the width of the credible intervals (`ci_width`), the share
of replications whose interval excludes zero (`detected`, which at a
true effect is the power to detect it), and the correlation between each
simulated person’s true and estimated values (`r`). [Running a recovery
study](https://www.gfrischkorn.org/bmmtools/articles/recovery-grid.html#the-result)
shows the subject-level `r` of `mixture2p` over 20 and 50 subjects
crossed with 30 and 100 trials, and [One grid, every
estimand](https://www.gfrischkorn.org/bmmtools/articles/recovery-grid.html#one-grid-every-estimand)
shows how the interval of an experimental effect narrows from 40 to 100
subjects.
[`?summary.bmmtools_recovery`](https://www.gfrischkorn.org/bmmtools/reference/summary.bmmtools_recovery.md)
defines every column.

Two planned additions make this easier: `design_check()`
([\#5](https://github.com/GidonFrischkorn/bmmtools/issues/5)) is
proposed to take a design in subjects, conditions and trials with
plausible parameter values and run the grid for you, and
[\#17](https://github.com/GidonFrischkorn/bmmtools/issues/17) adds
precision targets, the reliability of subject estimates and the smallest
design that meets a target.

## What bmmtools is not

- Not a simulation-based calibration package.
  [`sbc()`](https://www.gfrischkorn.org/bmmtools/reference/sbc.md)
  builds a generator from a bmm model specification and hands it to
  [SBC](https://hyunjimoon.github.io/SBC/), which owns the ranks, the
  ECDF diagnostics and the plots.
- Not a prior *sensitivity* package.
  [priorsense](https://CRAN.R-project.org/package=priorsense) answers
  how much the posterior moves when the prior is power-scaled.
  [`prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/prior_check.md)
  answers a different question: what the prior implies on the observable
  scale before any data are seen.
- Not a plotting layer. [bayesplot](https://mc-stan.org/bayesplot/) owns
  posterior plots. bmmtools adds the two plot families bayesplot has no
  concept of: the recovery scatter (truth against posterior estimate
  with credible intervals, one panel per parameter) and the
  prior-predictive panels on the observable scale.
- Not [gp3bayes](https://CRAN.R-project.org/package=gp3bayes), which
  also offers parameter recovery and SBC vignettes but is written for
  pupillometry and eye-tracking data. bmmtools is written for the
  measurement models bmm fits: mixture, signal detection, evidence
  accumulation, multinomial processing tree and memory measurement
  models.
- Not a reporting package. Every scorer returns a tidy tibble whose
  columns follow the apabayes contract, so a recovery table drops into a
  manuscript through `apabayes` with no glue code. Neither package
  imports the other.
- No model fitting of its own (that is `bmm()`).

## What bmmtools reads from bmm

bmmtools reads more of a fitted model than bmm’s exported API promises.
The scorers take `$links`, `$fixed_parameters`, `$resp_vars` and
`$other_vars` directly off the `bmmodel` object, and ten generator
adapters translate a model’s own parameters into the arguments of its
`r*()` function — `thetat` on the logit scale becomes `p_mem` for
`rmixture2p()`, for example. None of that is part of bmm’s documented
contract, so a bmm rename can break an adapter without a deprecation
warning on either side: when the rebuilt circular mixture models ship,
`mu1` becomes `mu`, and the adapter that reads `mu1` has to follow.
[Issue \#10](https://github.com/GidonFrischkorn/bmmtools/issues/10)
tracks what bmmtools reads and what it has asked bmm for in return.

## Installation

bmmtools is not on CRAN. Install the current release from GitHub:

``` r

# install.packages("pak")
pak::pak("GidonFrischkorn/bmmtools@v0.2.0")
```

Without the tag you get the development version, which is whatever is on
`main`:

``` r

pak::pak("GidonFrischkorn/bmmtools")
```

Fitting anything needs bmm, brms and a working
[cmdstanr](https://mc-stan.org/cmdstanr/) installation; they are
Suggests, so bmmtools installs and its score layer runs without them.

[`sbc()`](https://www.gfrischkorn.org/bmmtools/reference/sbc.md)
additionally needs the [SBC](https://hyunjimoon.github.io/SBC/) package,
which is not on CRAN either. `pak` resolves it from the `Remotes:` line
in DESCRIPTION, so the commands above already bring it in. To install it
on its own:

``` r

pak::pak("hyunjimoon/SBC")
```
