# Validating a model without a built-in adapter

bmmtools ships generator adapters for six bmm models: `sdt_yn`,
`sdt_mafc`, the three-parameter `ezdm`, `ddm`, `mixture2p` and `sdm`. A
model you are developing, for example one started with
[`bmm::use_model_template()`](https://venpopov.com/bmm/reference/use_model_template.html),
has none, and
[`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md)
stops with an error asking for `generator`. This article shows how to
write one, what the rest of the loop needs from it, and how to score
fits and estimates that did not come from the generate layer at all.

The example is bmm’s three-parameter mixture model, `mixture3p`, which
has no adapter. The code was run once and its results saved, so this
page was built without Stan.

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
    #> [1] "identity"
    #> 
    #> $thetant
    #> [1] "identity"

`kappa` arrives in the generator as a precision, but `thetat` and
`thetant` have an identity link: in `mixture3p` they are softmax weights
against a guessing weight fixed at 0. The generator has to turn them
into probabilities itself.

## Writing the generator

``` r

generate_mixture3p <- function(pars, n_trials, model) {
  # non-target locations relative to the target, one pair per trial
  nt <- matrix(stats::runif(2 * n_trials, -pi, pi), ncol = 2)
  # thetat and thetant are softmax weights against guessing, fixed at 0
  weights <- exp(c(pars$thetat, pars$thetant, 0))
  weights <- weights / sum(weights)
  y <- vapply(seq_len(n_trials), function(i) {
    bmm::rmixture3p(
      1,
      mu = c(pars$mu1, nt[i, ]), kappa = pars$kappa,
      p_mem = weights[[1]], p_nt = weights[[2]]
    )
  }, numeric(1))
  data.frame(y = y, nt1 = nt[, 1], nt2 = nt[, 2])
}
```

The non-target locations are drawn inside the generator, so they are
part of the simulated data.
[`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md)
applies its `seed` around the generator too, so these draws are
reproducible with the rest. The generator’s output is checked against
the model’s response columns before anything is fitted.

## Simulate, fit and score

From here the loop is the same as for a model with an adapter:

``` r

sim <- simulate_recovery(
  model,
  pars = c(kappa = log(8), thetat = 1.5, thetant = 0),
  sds = c(kappa = 0.3, thetat = 0.5),
  n_subjects = 30, n_trials = 100,
  generator = generate_mixture3p,
  seed = 1
)
sim$data
sim$truth$subjects
```

    #> # A tibble: 3,000 × 4
    #>    id           y     nt1     nt2
    #>    <fct>    <dbl>   <dbl>   <dbl>
    #>  1 1     -0.361    3.09   -1.49  
    #>  2 1     -0.193   -0.0277 -2.10  
    #>  3 1     -0.396   -0.0983 -1.12  
    #>  4 1     -0.410   -2.05    0.0636
    #>  5 1      0.00534  1.60    2.66  
    #>  6 1     -0.150   -0.290   0.0689
    #>  7 1     -1.95     0.0702 -1.52  
    #>  8 1      0.340   -1.84   -2.85  
    #>  9 1      0.922   -1.70   -0.516 
    #> 10 1     -0.127    0.601   2.22  
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
    #> 1     1.00         886.        1833.           0               0          98 TRUE  ""

``` r

recover(fit, sim$truth$population)
summary(recover_subjects(fit, sim$truth$subjects))
```

    #> <bmmtools_recovery>
    #> Scored on the natural scale.
    #> 1 fit, 3 parameters: kappa, thetat, thetant.
    #> Level: population.
    #> 
    #> # A tibble: 3 × 23
    #>   term    estimator level      scale       n n_replications n_converged    bias   rmse
    #>   <chr>   <chr>     <chr>      <chr>   <dbl>          <int>       <int>   <dbl>  <dbl>
    #> 1 kappa   bayes     population natural     1              1           1  0.120  0.120 
    #> 2 thetat  bayes     population natural     1              1           1  0.0867 0.0867
    #> 3 thetant bayes     population natural     1              1           1 -0.0157 0.0157
    #>   coverage ci_width     r r_low r_high rank_r   ccc ccc_low ccc_high ccc_accuracy
    #>      <dbl>    <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>   <dbl>    <dbl>        <dbl>
    #> 1        1    2.06     NA    NA     NA     NA    NA      NA       NA           NA
    #> 2        1    0.460    NA    NA     NA     NA    NA      NA       NA           NA
    #> 3        1    0.510    NA    NA     NA     NA    NA      NA       NA           NA
    #>   ccc_scale_shift ccc_location_shift calibration_slope truth_sd
    #>             <dbl>              <dbl>             <dbl>    <dbl>
    #> 1              NA                 NA                NA        0
    #> 2              NA                 NA                NA        0
    #> 3              NA                 NA                NA        0
    #> 
    #> r, rank_r and ccc are NA where they are not estimable: they need at least 3 complete pairs and spread on both sides.
    #> # A tibble: 2 × 23
    #>   term   estimator level   scale       n n_replications n_converged    bias  rmse coverage
    #>   <chr>  <chr>     <chr>   <chr>   <dbl>          <int>       <int>   <dbl> <dbl>    <dbl>
    #> 1 kappa  bayes     subject natural    30              1           1 -0.229  1.26     0.967
    #> 2 thetat bayes     subject natural    30              1           1  0.0188 0.250    0.967
    #>   ci_width     r r_low r_high rank_r   ccc ccc_low ccc_high ccc_accuracy ccc_scale_shift
    #>      <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>   <dbl>    <dbl>        <dbl>           <dbl>
    #> 1     5.57 0.818 0.649  0.910  0.788 0.769   0.604    0.871        0.940           0.717
    #> 2     1.01 0.771 0.569  0.886  0.689 0.742   0.548    0.860        0.962           0.759
    #>   ccc_location_shift calibration_slope truth_sd
    #>                <dbl>             <dbl>    <dbl>
    #> 1            -0.128               1.14    2.12 
    #> 2             0.0551              1.02    0.391

The same generator runs a whole design through
[`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md):

``` r

recovery_grid(
  model,
  grid = expand.grid(n_subjects = c(30, 60), n_trials = c(50, 100)),
  pars = c(kappa = log(8), thetat = 1.5, thetant = 0),
  sds = c(kappa = 0.3, thetat = 0.5),
  generator = generate_mixture3p,
  dir = "fits/mixture3p", reps = 10, seed = 1
)
```

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
#> # A tibble: 2 × 23
#>   term   estimator level      scale       n n_replications n_converged     bias   rmse
#>   <chr>  <chr>     <chr>      <chr>   <dbl>          <int>       <int>    <dbl>  <dbl>
#> 1 kappa  bayes     population natural     3              3          NA  0.256   0.647 
#> 2 thetat bayes     population natural     3              3          NA -0.00135 0.0180
#>   coverage ci_width     r r_low r_high rank_r   ccc ccc_low ccc_high ccc_accuracy
#>      <dbl>    <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>   <dbl>    <dbl>        <dbl>
#> 1        1    8.76  0.982    NA     NA      1 0.961      NA       NA        0.979
#> 2        1    0.220 0.988    NA     NA      1 0.966      NA       NA        0.977
#>   ccc_scale_shift ccc_location_shift calibration_slope truth_sd
#>             <dbl>              <dbl>             <dbl>    <dbl>
#> 1           1.19              0.111              0.824   2.10  
#> 2           0.806            -0.0200             1.23    0.0755
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

## Contributing an adapter

If your model is part of bmm and its generator fits the pattern above,
an adapter in bmmtools saves every user from writing it. Adapters live
in `R/adapters.R`: a function per model that maps one subject’s
natural-scale parameters onto the model’s `r*()` function and names the
output columns from the model object, plus an entry in the table that
looks them up by class. A response range for
[`prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/prior_check.md)
goes in `response_range()` in `R/prior-check.R`. Open an issue or a pull
request on [GitHub](https://github.com/GidonFrischkorn/bmmtools).
