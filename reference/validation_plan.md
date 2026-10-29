# Plan the validation of a bmm model

**\[experimental\]**

## Usage

``` r
validation_plan(
  model,
  pars = NULL,
  sds = NULL,
  cors = NULL,
  grid = NULL,
  reps = NULL,
  sampler = NULL,
  sbc = NULL,
  prior_check = NULL,
  reference = NULL,
  convergence = NULL,
  criteria = NULL,
  stages = NULL,
  formula = NULL,
  trial_design = NULL,
  generator = NULL,
  seed,
  owner = NULL,
  version = NULL
)
```

## Arguments

- model:

  A `bmmodel`, or the name of a model with a
  [`validation_defaults()`](https://www.gfrischkorn.org/bmmtools/reference/validation_defaults.md)
  entry, which is then built with the entry's variable names (bmm is
  needed for that).

- pars:

  Named population values on the link scale, as in
  [`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md).

- sds:

  Named between-subject SDs on the link scale; `NULL` takes the entry's.

- cors:

  A correlation matrix between subject values, as in
  [`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md);
  `NULL` takes the entry's, which may be none.

- grid:

  A data frame with the columns `n_subjects` and `n_trials`, one row per
  design.

- reps:

  Replications per design.

- sampler:

  A named list passed to the fitter of every recovery cell (`chains`,
  `iter`, `warmup`, `cores`, `backend`, ...). A list given here changes
  the entries it names and keeps the others. It may not set `prior` or
  an abbreviation of it (`pri`), nor may `sbc$sampler` or
  `prior_check$sampler`: the calibration's prior is `sbc$prior`, and the
  recovery stage fits under bmm's default priors.

- sbc:

  A named list for the calibration: `n_subjects`, `n_trials`, `n_sims`,
  `prior` (a `brmsprior`, a data frame with the columns of
  [`brms::set_prior()`](https://paulbuerkner.com/brms/reference/set_prior.html),
  or `NULL` for bmm's defaults), `prior_reason` (why it is not bmm's
  default, printed on the card), `level`, `sampler`, `generator`,
  `cores_per_fit` and `thin_ranks` (see
  [`sbc()`](https://www.gfrischkorn.org/bmmtools/reference/sbc.md)).
  Named entries replace the entry's, the rest are kept, except that a
  `prior` other than the entry's, given without a `prior_reason`, drops
  the entry's reason. The prior is compared as given, column order and
  types included: a `brmsprior` with the same rows as an entry's data
  frame, the entry's columns in another order, or a factor column counts
  as another prior. A `prior_reason` without a `prior` is an error.

- prior_check:

  A named list for the prior predictive check: `n_subjects`, `n_trials`,
  `n_draws`, `range`, `range_reason` (why the range is what it is,
  printed on the card), `summary` and `sampler` (see
  [`prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/prior_check.md)).
  Merged as `sbc` is: a `range` other than the entry's, given without a
  `range_reason`, drops the entry's reason, and a `range_reason` without
  a `range` is an error. The range is compared as given, names included,
  so `c(lo = 0.15, hi = 5)` counts as another range than the entry's
  `c(0.15, 5)`.

- reference:

  `NULL`, or a reference implementation's estimates for the agreement
  stage: a data frame with `term` and `estimate` (and optionally
  `ci_low`, `ci_high`) on the natural scale, as
  [`cross_check()`](https://www.gfrischkorn.org/bmmtools/reference/cross_check.md)
  takes it, plus the provenance columns `package`, `version`, `call`,
  `data_hash`, `date` and `source`.

- convergence:

  `NULL` for
  [`check_convergence()`](https://www.gfrischkorn.org/bmmtools/reference/check_convergence.md)'s
  defaults, or a named list of its arguments, as in
  [`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md).

- criteria:

  A named list replacing some of the default criteria (see "Criteria").

- stages:

  `NULL` for every stage that exists, or the names of the stages to run
  (see "Stages"). The run record always runs; a stage that reads the
  recovery grid needs `"recovery"` too.

- formula, trial_design, generator:

  As in
  [`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md).
  `NULL` `formula` means
  [`recovery_formula()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_formula.md)
  with correlated subject effects when there are at least two
  subject-varying parameters, so the individual-differences stage has
  correlations to score.

- seed:

  The master seed, required. Each stage derives its own from it (see
  "Seeds").

- owner, version:

  The entry's owner and plan version; `NULL` takes the entry's.

## Value

A `bmmtools_validation_plan`: a named list with every field above
resolved, plus `model_name`, `defaults` (`"bmmtools <version> entry"` or
`"none"`), `re_cor`, `correlations`, `seeds`, `stages` (a tibble:
`stage`, `order`, `title`, `required`, `status`, `issue`, `reason`),
`hash` and `bmmtools_version`.

## Details

Builds a `bmmtools_validation_plan`: everything
[`run_validation()`](https://www.gfrischkorn.org/bmmtools/reference/run_validation.md)
needs to validate one model, resolved before anything is fitted. Every
argument left `NULL` is filled from the model's
[`validation_defaults()`](https://www.gfrischkorn.org/bmmtools/reference/validation_defaults.md)
entry, or for a model without one from the generic grid, in which case
`pars` is required.

## Stages

In card order: `run_record`, `prior_predictive`, `sbc`, `recovery`,
`computation`, `identifiability`, `individual_differences`,
`data_requirements`, `agreement`, `model_recovery` and
`posterior_checks`. A stage is `"planned"`; `"not built"` with the issue
that will build it (model recovery, \#6; posterior checks, \#20);
`"not applicable"` with the reason (individual differences with fewer
than two subject-varying parameters; data requirements without a
subject-varying parameter, or without two numbers of trials at
`read_off_n_subjects` subjects in the grid; agreement without a
reference); or `"excluded"` when `stages` leaves it out. All but
agreement, model recovery and posterior checks are required for a card
to be called complete.

## Criteria

`prior_bound_rate_max` (0.10): the largest floor or ceiling rate of the
prior predictive draws under the calibration's prior; under bmm's
defaults the rates are reported, not judged, unless the calibration uses
them. `sbc_ecdf_prob` (0.95): the simultaneous band of SBC's
ECDF-difference check no ranked variable may leave.
`convergence_rate_min` (0.95): the share of recovery cells that pass
[`check_convergence()`](https://www.gfrischkorn.org/bmmtools/reference/check_convergence.md).
`coverage_level` (0.95) and `coverage_alpha` (0.05): the coverage of
every population parameter, its intervals pooled over the designs, must
not fall below `coverage_level` by a one-sided binomial test at
`coverage_alpha` divided by the number of parameters; the generating
correlations are held to the same test. `contraction_min` (0) and
`z_mean_max` (2): the identifiability bars at the largest design.
`read_off_r_low` (0.80) and `read_off_n_subjects` (50): the data
requirements read-off, the smallest `n_trials` at which every
subject-varying parameter's person-level `r_low` reaches the bar at that
many subjects; 0.80 is a convention, and the read-off at 0.70 and 0.90
is reported beside it. `tolerance`: by metric, how far a number without
a Monte Carlo standard error may move between two runs before
[`compare_validations()`](https://www.gfrischkorn.org/bmmtools/reference/compare_validations.md)
flags it.

## Seeds

The recovery grid runs under `seed`, the prior check under `seed + 1`
and the calibration under `seed + 2`. The seeds, the owner and the
reasons (`sbc$prior_reason`, `prior_check$range_reason`) are left out of
the plan's hash and of the stages' keys, so the hash changes only when
what is planned changes, and rewording a reason reruns nothing.

## See also

[`run_validation()`](https://www.gfrischkorn.org/bmmtools/reference/run_validation.md),
[`validation_defaults()`](https://www.gfrischkorn.org/bmmtools/reference/validation_defaults.md).

## Examples

``` r
plan <- validation_plan(
  structure(
    list(links = list(kappa = "log", thetat = "logit")),
    class = c("bmmodel", "circular", "mixture2p")
  ),
  seed = 1
)
plan
#> <bmmtools validation plan> mixture2p, version 0.1.0, hash c6c12584
#> Owner: Gidon Frischkorn; defaults: bmmtools 0.3.0 entry
#> Truth (link scale): kappa = 2.08, thetat = 1.1
#> Subject SDs (link scale): kappa = 0.3, thetat = 0.5
#> Grid: 6 designs x 10 replications (n_subjects 20, 50; n_trials 30, 100, 300)
#> Sampler: chains = 4, iter = 2000, warmup = 1000, cores = 4, backend = "cmdstanr"
#> SBC: 100 simulations at 20 subjects x 50 trials under the plan's own prior
#> Prior check: 500 draws at 20 subjects x 50 trials
#> Reference: none
#> Seed: 1
#> Stages:
#>    1 run_record               planned        
#>    2 prior_predictive         planned        
#>    3 sbc                      planned        
#>    4 recovery                 planned        
#>    5 computation              planned        
#>    6 identifiability          planned        
#>    7 individual_differences   planned        
#>    8 data_requirements        planned        
#>    9 agreement                not applicable  no reference implementation
#>   10 model_recovery           not built       not run (#6)
#>   11 posterior_checks         not built       not run (#20)
if (FALSE) { # \dontrun{
plan <- validation_plan("ezdm", seed = 20261024)
} # }
```
