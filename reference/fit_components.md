# Fit every component of a simulation set, each on its own

Fits each component's data with its own model through
[`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md),
so a second call reuses the saved fits while nothing that determines
them has changed.

## Usage

``` r
fit_components(sim_set, dir, prior = NULL, ..., .fitter = NULL)
```

## Arguments

- sim_set:

  A `bmmtools_simulation_set` from
  [`simulate_components()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_components.md).

- dir:

  The directory the fits are saved in, one file per component, named
  after it (`<dir>/<name>.rds`, with its cache key next to it).

- prior:

  `NULL`, or a named list with a prior for some of the components, such
  as `list(sdm = <brmsprior>)`. A component not named is fitted with
  bmm's default priors.

- ...:

  Passed to
  [`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md)
  and on to the fitter, for every component: `refit`, and sampler
  arguments such as `chains`, `iter` or `backend`. When the set was
  simulated with a seed and `...` has no `seed`, the set's seed is
  passed as `seed`, the same for every component.

- .fitter:

  The fitting function,
  [`bmm::bmm()`](https://venpopov.com/bmm/reference/bmm.html) by
  default; see
  [`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md).

## Value

A named list of fits, one per component, of class `bmmtools_fit_set`,
which
[`extract_correlations()`](https://www.gfrischkorn.org/bmmtools/reference/extract_correlations.md),
[`recover_correlations()`](https://www.gfrischkorn.org/bmmtools/reference/recover_correlations.md)
and
[`subject_table()`](https://www.gfrischkorn.org/bmmtools/reference/subject_table.md)
read as one set.

## Details

A component fitted without a `formula` uses
[`recovery_formula()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_formula.md),
with the task column when the component has tasks.

## Examples

``` r
if (FALSE) { # \dontrun{
fits <- fit_components(set, dir = "fits", chains = 4, iter = 2000)
fits
} # }
```
