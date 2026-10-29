# Validation defaults for a model

**\[experimental\]**

## Usage

``` r
validation_defaults(model = NULL)
```

## Arguments

- model:

  `NULL` for every entry; a model name (`"mixture2p"`) or a `bmmodel`
  for that model's.

## Value

With `model = NULL`, a named list of entries. Otherwise the entry, or
`NULL` with no error for a model without one. An entry is a list:
`model_args` (the arguments its constructor is called with when
[`validation_plan()`](https://www.gfrischkorn.org/bmmtools/reference/validation_plan.md)
is given the name), `pars`, `sds` and `cors` (the truth on the link
scale), `grid` (`n_subjects` by `n_trials`), `reps`, `sampler` (passed
to the fitter), `sbc` (the calibration's design point, number of
simulations, prior and the reason that prior is not bmm's default,
levels and sampler), `prior_check` (its design point, draws and
sampler), `reference`, `owner`, `version` and `sections` (model-specific
card sections).

## Details

The generating values, the standard grid and the budget a
[`validation_plan()`](https://www.gfrischkorn.org/bmmtools/reference/validation_plan.md)
fills from when an argument is left `NULL`. An entry exists for a model
once its validation card is made; a model without one, including a model
of your own, is planned with the generic grid and needs `pars`.

The standard grid is 20 and 50 subjects by 30, 100 and 300 trials per
subject, ten replications each, fitted with 4 chains of 2000 iterations
(1000 warmup) under cmdstanr; SBC runs 100 simulations at 20 subjects by
50 trials with 2 chains of 500. Ten replications give a coverage a Monte
Carlo standard error of about 0.07 to 0.1, which the card prints next to
every rate.

Two entries ship:

- `mixture2p`: `kappa = log(8)`, `thetat = qlogis(0.75)`, SDs 0.3 and
  0.5. Its calibration runs under a narrower prior than bmm's. The
  default SD prior on log `kappa` of a bmm development build drew values
  at which
  [`bmm::rmixture2p()`](https://venpopov.com/bmm/reference/mixture2p_dist.html)
  overflows the stack; the defaults of bmm 1.4.1 from CRAN rarely reach
  them, and the narrower prior keeps the calibration from moving with
  bmm's defaults. The entry records the reason and the card prints it.

- `ezdm` (three parameters): `drift = 1.5`, `bound = 0`,
  `ndt = log(0.3)`, SDs 0.3, 0.2 and 0.2, a correlation of 0.4 between
  drift and bound. Its calibration also runs under a narrower prior than
  bmm's, because bmm's `cauchy(0, 1)` drift intercept and exponential SD
  priors draw parameter values (a drift of -31, an SD of 2 on log bound)
  at which some fits fail the card's convergence criterion.

Each entry names its `owner`, who changes its values by pull request;
the card prints the owner and the plan version.

## See also

[`validation_plan()`](https://www.gfrischkorn.org/bmmtools/reference/validation_plan.md),
which fills a plan from an entry.

## Examples

``` r
names(validation_defaults())
#> [1] "mixture2p" "ezdm"     
validation_defaults("ezdm")$pars
#>     drift     bound       ndt 
#>  1.500000  0.000000 -1.203973 
validation_defaults("sdm")
#> NULL
```
