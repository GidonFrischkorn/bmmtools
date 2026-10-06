# Validating a model without a built-in adapter

bmmtools ships generator adapters for ten bmm models: `sdt_yn`,
`sdt_mafc`, the three-parameter `ezdm`, `ddm`, `cswald` (both the
`simple` and the `crisk` version), `mixture2p`, `sdm`, `mixture3p`,
`imm` (all three versions) and `m3` (`ss`, `cs` and `custom`). A model
you are developing, for example one started with
[`bmm::use_model_template()`](https://venpopov.com/bmm/reference/use_model_template.html),
has none, and
[`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md)
stops with an error asking for `generator`. This article shows how to
write one, including one that needs to know what was shown on each
trial, what the rest of the loop needs from it, and how to score fits
and estimates that did not come from the generate layer at all.

The example is bmm’s three-parameter mixture model, `mixture3p`. It has
a built-in adapter, so you would not write this generator for it; it is
here because its responses depend on the non-target locations of each
trial, which is the case a generator of your own most often has to
handle. The code was run once and its results saved, so this page was
built without Stan.

## What a generator is

A generator is a function `function(pars, n_trials, model)` that returns
the rows of **one subject** as a data frame, with the column names the
model expects.
[`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md)
calls it once per subject and adds the `id` column itself.

- `pars` is a named list of that subject’s parameter values on the
  **natural scale**, fixed parameters included.
- `n_trials` is the number of trials per subject, as given to
  [`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md)
  or in a grid row.
- `model` is the model object, so the generator can read column names or
  other settings from it.
- `trial_design`, a fourth argument, holds that subject’s rows of the
  per-trial design when
  [`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md)
  is given one (below). A generator for a model that needs no per-trial
  input leaves it out.

[`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md)
takes `pars` on the link scale, draws subject values there and converts
them with
[`inverse_link()`](https://www.gfrischkorn.org/bmmtools/reference/inverse_link.md)
using the model’s link table before calling the generator. Check that
table first:

``` r

library(bmmtools)

model <- bmm::mixture3p(
  resp_error = "y", nt_features = c("nt1", "nt2"), set_size = 3
)
model$links
```

    #> $mu1
    #> [1] "tan_half"
    #> 
    #> $kappa
    #> [1] "log"
    #> 
    #> $thetat
    #> [1] "softmax"
    #> 
    #> $thetant
    #> [1] "softmax"

`kappa` arrives in the generator as a precision. `thetat` and `thetant`
have bmm’s `softmax` link: they are log weights against a guessing
weight fixed at 0, and a trial’s probabilities are the softmax of all
three, so no inverse of one of them alone gives its natural value.
[`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md)
therefore passes them to the generator **on the link scale**, and the
generator turns them into probabilities itself.

## The trial design

Each response depends on where the non-targets were on that trial, and
the fit needs those locations as data columns. They are the trial
design: what was shown, as opposed to what the subject did. Give it to
[`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md)
as `trial_design`, and it reaches both the generator and the data:

``` r

# non-target locations relative to the target, one pair per trial
lures <- function(n_trials) {
  nt <- matrix(stats::runif(2 * n_trials, -pi, pi), ncol = 2)
  data.frame(nt1 = nt[, 1], nt2 = nt[, 2])
}
```

`trial_design` takes a data frame of `n_trials` rows (one trial list for
every subject), a data frame with an `id` column (one list per subject),
or, as here, a function of `n_trials`, called once per subject under the
seed so that every subject sees new locations.

## Writing the generator

``` r

my_mixture3p <- function(pars, n_trials, model, trial_design) {
  # thetat and thetant arrive on the link scale: log weights against a
  # guessing weight fixed at 0
  weights <- exp(c(pars$thetat, pars$thetant, 0))
  weights <- weights / sum(weights)
  y <- vapply(seq_len(n_trials), function(i) {
    bmm::rmixture3p(
      1,
      mu = c(pars$mu1, trial_design$nt1[[i]], trial_design$nt2[[i]]),
      kappa = pars$kappa, p_mem = weights[[1]], p_nt = weights[[2]]
    )
  }, numeric(1))
  data.frame(y = y)
}
```

The generator reads its subject’s trials from `trial_design` and returns
the response column only.
[`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md)
binds the design columns into the data itself, so the fit sees exactly
the trials the responses were drawn from; a generator that returned
`nt1` as well would be stopped with an error. Its output is checked
against the model’s response columns before anything is fitted.

## Simulate, fit and score

From here the loop is the same as for a model with an adapter:

``` r

sim <- simulate_recovery(
  model,
  pars = c(kappa = log(8), thetat = 1.5, thetant = 0),
  sds = c(kappa = 0.3, thetat = 0.5),
  n_subjects = 30, n_trials = 100,
  trial_design = lures,
  generator = my_mixture3p,
  seed = 1
)
sim$data
sim$truth$subjects
```

    #> # A tibble: 3,000 × 4
    #>    id        nt1     nt2        y
    #>    <fct>   <dbl>   <dbl>    <dbl>
    #>  1 1      3.09   -1.49    0.265  
    #>  2 1     -0.0277 -2.10    0.0348 
    #>  3 1     -0.0983 -1.12   -0.679  
    #>  4 1     -2.05    0.0636 -0.429  
    #>  5 1      1.60    2.66    0.272  
    #>  6 1     -0.290   0.0689  1.11   
    #>  7 1      0.0702 -1.52    0.345  
    #>  8 1     -1.84   -2.85   -0.553  
    #>  9 1     -1.70   -0.516  -0.653  
    #> 10 1      0.601   2.22   -0.00603
    #> # ℹ 2,990 more rows
    #> # A tibble: 60 × 3
    #>    id    term  true_value
    #>    <chr> <chr>      <dbl>
    #>  1 1     kappa       1.89
    #>  2 2     kappa       2.13
    #>  3 3     kappa       1.83
    #>  4 4     kappa       2.56
    #>  5 5     kappa       2.18
    #>  6 6     kappa       1.83
    #>  7 7     kappa       2.23
    #>  8 8     kappa       2.30
    #>  9 9     kappa       2.25
    #> 10 10    kappa       1.99
    #> # ℹ 50 more rows

`thetant` is in `pars` but not in `sds`, so it does not vary between
subjects and has no rows in `truth$subjects`: a parameter that does not
vary has nothing person-level to recover.

``` r

fit <- fit_cached(
  recovery_formula(model), sim$data, model,
  file = file.path(fits_dir, "mixture3p-30x100"),
  seed = 1, chains = 4, iter = 2000, cores = 4, backend = "cmdstanr"
)
check_convergence(fit)
```

    #> # A tibble: 1 × 8
    #>   max_rhat min_ess_bulk min_ess_tail n_divergent n_max_treedepth n_variables pass  failed
    #>      <dbl>        <dbl>        <dbl>       <int>           <int>       <int> <lgl> <chr> 
    #> 1     1.01         840.        1442.           0               0          98 TRUE  ""

``` r

recover(fit, sim$truth$population)
summary(recover_subjects(fit, sim$truth$subjects))
```

    #> <bmmtools_recovery>
    #> Scored on the natural scale.
    #> 1 fit, 3 parameters: kappa, thetat, thetant.
    #> Level: population.
    #> thetat, thetant are on the link scale.
    #> 
    #> # A tibble: 3 × 15
    #>   term    level      scale       n    bias   rmse coverage coverage_50 z_mean  z_sd
    #>   <chr>   <chr>      <chr>   <dbl>   <dbl>  <dbl>    <dbl>       <dbl>  <dbl> <dbl>
    #> 1 kappa   population natural     1 -0.448  0.448         1           0 -0.814    NA
    #> 2 thetat  population link        1  0.0712 0.0712        1           1  0.579    NA
    #> 3 thetant population link        1 -0.0331 0.0331        1           1 -0.228    NA
    #>   contraction     r   ccc calibration_slope calibration_intercept
    #>         <dbl> <dbl> <dbl>             <dbl>                 <dbl>
    #> 1          NA    NA    NA                NA                    NA
    #> 2          NA    NA    NA                NA                    NA
    #> 3          NA    NA    NA                NA                    NA
    #> # 19 more columns; `tibble::as_tibble()` prints them all.
    #> 
    #> r, rank_r and ccc are NA where they are not estimable: they need at least 3 complete pairs and spread on both sides.
    #> # A tibble: 2 × 15
    #>   term   level   scale       n     bias  rmse coverage coverage_50  z_mean  z_sd
    #>   <chr>  <chr>   <chr>   <dbl>    <dbl> <dbl>    <dbl>       <dbl>   <dbl> <dbl>
    #> 1 kappa  subject natural    30 -0.712   1.38     0.967       0.567 -0.460  0.772
    #> 2 thetat subject link       30  0.00280 0.200    1           0.6    0.0158 0.746
    #>   contraction     r   ccc calibration_slope calibration_intercept
    #>         <dbl> <dbl> <dbl>             <dbl>                 <dbl>
    #> 1          NA 0.831 0.764              1.01                0.643 
    #> 2          NA 0.860 0.844              1.04               -0.0719
    #> # 19 more columns; `tibble::as_tibble()` prints them all.

The same generator runs a whole design through
[`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md):

``` r

recovery_grid(
  model,
  grid = expand.grid(n_subjects = c(30, 60), n_trials = c(50, 100)),
  pars = c(kappa = log(8), thetat = 1.5, thetant = 0),
  sds = c(kappa = 0.3, thetat = 0.5),
  trial_design = lures,
  generator = my_mixture3p,
  dir = "fits/mixture3p", reps = 10, seed = 1
)
```

With the built-in adapter, the same grid needs no `generator`, and the
design names the set size too when it varies; see the section “Trial
design” of
[`?simulate_recovery`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md).

`thetat` and `thetant` stay on the link scale when they are scored:
[`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
says so in a message, and their `scale` column reads `"link"`.

## Scoring estimates you already have

[`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
and
[`recover_subjects()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
accept an **estimates tibble** in place of a fit. That is how to score
estimates from another fitting routine, from an older simulation, or
from any pipeline where no fit object is at hand. No fitting package
needs to be installed.

The tibble needs the columns `term`, `estimate`, `ci_low`, `ci_high`,
`ci_method`, `ci_level`, `rhat`, `ess_bulk`, `ess_tail`, `level` and
`id`, and optionally `converged` and `replication`. Estimates and truth
are on the link scale; `links` maps each term to its link so that
scoring can happen on the natural scale:

``` r

estimates <- tibble::tibble(
  term = rep(c("kappa", "thetat"), times = 3),
  estimate = c(2.10, 1.00, 1.70, 0.60, 2.45, 1.35),
  ci_low = c(1.60, 0.40, 1.20, 0.05, 1.95, 0.80),
  ci_high = c(2.60, 1.60, 2.20, 1.15, 2.95, 1.90),
  ci_method = "eti", ci_level = 0.95,
  rhat = 1, ess_bulk = 900, ess_tail = 900,
  level = "population", id = NA_character_,
  replication = rep(1:3, each = 2)
)
truth <- tibble::tibble(
  term = rep(c("kappa", "thetat"), times = 3),
  true_value = c(2.0, 1.1, 1.8, 0.5, 2.4, 1.4),
  replication = rep(1:3, each = 2)
)

recovered <- recover(estimates, truth, links = c(kappa = "log", thetat = "logit"))
summary(recovered)
#> # A tibble: 2 × 15
#>   term   level      scale       n     bias   rmse coverage coverage_50 z_mean  z_sd
#>   <chr>  <chr>      <chr>   <dbl>    <dbl>  <dbl>    <dbl>       <dbl>  <dbl> <dbl>
#> 1 kappa  population natural     3  0.256   0.647         1          NA     NA    NA
#> 2 thetat population natural     3 -0.00135 0.0180        1          NA     NA    NA
#>   contraction     r   ccc calibration_slope calibration_intercept
#>         <dbl> <dbl> <dbl>             <dbl>                 <dbl>
#> 1          NA 0.982 0.961             0.824                 1.22 
#> 2          NA 0.988 0.966             1.23                 -0.163
#> # 19 more columns; `tibble::as_tibble()` prints them all.
```

`n_converged` is `NA` because the tibble carries no convergence verdict.
When `truth` has a `replication` column it is used as a join key;
without one, the same truth is scored against every replication.

[`extract_estimates()`](https://www.gfrischkorn.org/bmmtools/reference/extract_estimates.md)
produces this tibble from a fit, so a scoring pipeline can be written
and tested on saved tibbles and pointed at fits later.

## Scoring a brms fit

Scoring dispatches on `brmsfit`, and a `bmmfit` is one, so a plain brms
fit can be scored as well. Without bmm’s link table, scoring falls back
to the link scale and says so, unless `links` is given:

``` r

recover(brms_fit, truth, links = c(kappa = "log", thetat = "logit"))
```

Term names are the coefficient names with the `b_` prefix and the
coefficient suffix removed, so `b_kappa_Intercept` becomes `kappa`.
`extract_estimates(brms_fit)$term` shows what `truth` has to match. A
term in `truth` that the fit did not estimate is dropped with a warning
that lists the available terms.

## Limits of this version

- Each parameter is scored through its intercept, or through one cell
  mean per task (`kappa ~ 0 + task`, the terms `kappa_task1`,
  `kappa_task2`, …), which is what
  [`recovery_formula()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_formula.md)
  writes with and without `task_col`. A formula with an intercept and
  other coefficients for one parameter, say `kappa ~ 1 + setsize`, stops
  with an error rather than join the estimates wrongly.
- Group-level standard deviations are scored on the link scale only
  (`recover(level = "sd")`). Group-level correlations are extracted
  (`extract_estimates(level = "cor")`) but not yet scored.
- [`prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/prior_check.md)
  derives the floor and ceiling of the response only for the models with
  an adapter; for your model, pass `range` or a `summary` function of
  your own.
- A parameter with a `softmax` link is scored on the link scale only.

## Contributing an adapter

If your model is part of bmm and its generator fits the pattern above,
an adapter in bmmtools saves every user from writing it. Adapters live
in `R/adapters.R`: a function per model that maps one subject’s
natural-scale parameters onto the model’s `r*()` function and names the
output columns from the model object, plus an entry in the table that
looks them up by class, and, when it reads a trial design, the columns
it reads in `trial_design_columns()`. A response range for
[`prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/prior_check.md)
goes in `response_range()` in `R/prior-check.R`. Open an issue or a pull
request on [GitHub](https://github.com/GidonFrischkorn/bmmtools).
