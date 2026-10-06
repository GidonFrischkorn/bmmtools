# Check what the priors say the data should look like

The question to ask before fitting anything: a prior on a link-scale
parameter is not interpretable, and the same prior expressed as "this
share of prior draws implies chance-level responding" is. Draws come
from `bmm(sample_prior = "only")` through
[`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md),
so the prior being checked is exactly the prior bmm fits with — exact
for every prior brms accepts, and it gives predictions on the observable
scale for free.

## Usage

``` r
prior_check(
  model,
  formula,
  data,
  prior = NULL,
  summary = NULL,
  n_draws = 1000,
  range = NULL,
  file = NULL,
  refit = c("on_change", "never", "always"),
  seed = NULL,
  ...,
  .fitter = NULL
)
```

## Arguments

- model:

  A `bmmodel`, built with the column names `data` uses.

- formula:

  The `bmmformula` you intend to fit. The prior a model needs depends on
  it, so there is no default.

- data:

  The design the prior is checked over. The response column must exist
  because brms builds its Stan data from it; its values do not enter the
  prior draws, and are used only as the observed reference in
  [`plot_prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/plot_prior_check.md).

- prior:

  `NULL` checks bmm's own defaults. A `brmsprior` checks that one. A
  **named** list runs one fit per element and compares them, a `NULL`
  element meaning bmm's defaults, so
  `list(default = NULL, tight = my_prior)` is the usual comparison.

- summary:

  A function `(yrep, data)` returning a data frame with at least
  `statistic` and `value`, and optionally `response`. `NULL` uses the
  default: the share of draws at the floor and at the ceiling of the
  response, and the 50th, 90th and 95th percentiles of the
  prior-predicted observable.

- n_draws:

  Prior-predictive draws kept per prior set. More than the fit holds is
  reported and capped, not refused.

- range:

  The floor and the ceiling of the response, as two numbers. `NULL`
  derives them from the model where they are known: `0` and the trial
  count for `sdt_yn` and `sdt_mafc`, `-pi` and `pi` for `mixture2p` and
  `sdm`. Where they are not known the two rates are `NA` rather than
  invented.

- file:

  Where to cache the fits. `NULL` uses a temporary file, so the compile
  is reused within the session and nothing is left behind; a path caches
  durably and each set gets `<file>-<set>`.

- refit:

  Passed to
  [`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md).

- seed:

  Passed to the fitter, where it enters the cache key, and used around
  the draw subsample so that the same seed keeps the same `n_draws`.
  `NULL` leaves the random number generator alone and is recorded as
  `NA`.

- ...:

  Passed to
  [`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md)
  and on to the fitter: `chains`, `iter`, `backend`, `cores`.
  `sample_prior` is set here and giving it is an error.

- .fitter:

  The fitting function,
  [`bmm::bmm()`](https://venpopov.com/bmm/reference/bmm.html) by
  default. Tests inject a stand-in so that no model is compiled.

## Value

A `bmmtools_prior_check` tibble, one row per prior set, response and
statistic, with the columns `prior`, `response`, `statistic` and
`value`. It carries the draws themselves as an attribute, because a
summary cannot be plotted as a distribution, along with
[`brms::prior_summary()`](https://mc-stan.org/rstantools/reference/prior_summary.html)
of each fit — the prior that actually produced the draws, rather than
the one that was asked for.

## Details

A population-level slope left without a proper prior is warned about
before fitting: `sample_prior = "only"` samples every parameter and a
flat prior has nothing to sample from. It is a warning rather than an
error because bmm's defaults cover the intercepts and the group-level
SDs, and only a slope you added is likely to be missing.

A floor or ceiling rate is meaningful for a **discrete** response — a
count out of `n_trials` — where it is the share of prior draws implying
perfect or floor performance. For a **continuous** response exact
equality to a boundary has probability zero and the quantiles are what
to read. A circular response error is best read through a summary of
your own on the absolute error, because the quantiles of a signed error
are symmetric about zero and say nothing:

    prior_check(model, formula, data, summary = function(yrep, data) {
      tibble::tibble(
        statistic = c("q50", "q90"),
        value = unname(stats::quantile(abs(yrep), c(0.5, 0.9)))
      )
    })

## Examples

``` r
if (FALSE) { # \dontrun{
model <- bmm::mixture2p(resp_error = "y")
checked <- prior_check(
  model, recovery_formula(model), my_data,
  n_draws = 500, seed = 1
)
} # }

# a precomputed check comparing two prior sets
prior_check_sdt_yn
#> <bmmtools_prior_check>
#> Model: Signal Detection Theory (Yes/No); 2 prior sets: default, narrow_sd.
#> Prior-predictive draws per set: 100.
#> 
#> # A tibble: 5 × 5
#>   response statistic    default narrow_sd difference
#>   <chr>    <chr>          <dbl>     <dbl>      <dbl>
#> 1 hits     floor_rate     0.392     0.152     -0.241
#> 2 hits     ceiling_rate   0.420     0.146     -0.274
#> 3 hits     q50           30        23         -7    
#> 4 hits     q90           50        50          0    
#> 5 hits     q95           50        50          0    
summary(prior_check_sdt_yn)
#> # A tibble: 5 × 5
#>   response statistic    default narrow_sd difference
#>   <chr>    <chr>          <dbl>     <dbl>      <dbl>
#> 1 hits     floor_rate     0.392     0.152     -0.241
#> 2 hits     ceiling_rate   0.420     0.146     -0.274
#> 3 hits     q50           30        23         -7    
#> 4 hits     q90           50        50          0    
#> 5 hits     q95           50        50          0    
```
