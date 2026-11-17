# Fit a bmm model once and reuse it while nothing has changed

A cache for the fits a recovery study produces. brms's own `file`
argument reuses a saved fit when the Stan code, the Stan data and the
algorithm are unchanged, and does not look at the seed, the number of
chains or iterations, the warmup, thinning, `control`, `init` or the
backend (measured on brms 2.23.0). In a simulation grid the cells that
share a design differ in exactly those, so `fit_cached()` keys on all of
them, and on the formula, data, model, prior, the installed versions of
bmm and brms and the Stan toolchain of the backend in use, and records
each component's hash next to the fit so that a mismatch is reported by
name. Without a `backend`, the backend in use is the one bmm picks: the
`brms.backend` option if set, otherwise cmdstanr if it is installed,
otherwise rstan.

## Usage

``` r
fit_cached(
  formula,
  data,
  model,
  file,
  prior = NULL,
  refit = c("on_change", "never", "always"),
  ...,
  .fitter = NULL
)
```

## Arguments

- formula, data, model, prior:

  Passed to the fitter as they are; the formula enters the key deparsed,
  so the environment it was built in does not matter.

- file:

  Path of the cache file. `.rds` is appended when the path has no
  extension; the key is written next to it as `<file>.key`. The
  directory is created before fitting, so an unwritable path fails
  before any sampling.

- refit:

  `"on_change"` (the default) fits when nothing is cached or when any
  key component changed; `"never"` returns whatever is cached and warns
  if it is stale; `"always"` fits and overwrites.

- ...:

  Passed to the fitter. `seed`, `chains`, `iter`, `warmup`, `thin`,
  `control`, `init`, `backend`, `sample_prior`, `algorithm` and
  `stanvars` enter the key, and so does `draws` (the number of draws of
  a Laplace or variational fit) when it is given; anything else
  (`cores`, `refresh`, `silent`, ...) changes how the fit is run, not
  what it is, and does not. An `init` given as a function enters the key
  by its text, so a closure built afresh in every iteration of a loop
  still matches. brms's `file_refit` and `file_compress` are refused:
  there is one cache, not two. The fitter runs under `seed`, see "Seed"
  below.

- .fitter:

  The fitting function,
  [`bmm::bmm()`](https://venpopov.com/bmm/reference/bmm.html) by
  default. Tests inject a stand-in so that no model is compiled.

## Value

The fit, with an attribute `bmmtools_cache`: a list of `file`, `key`,
`reused`, `seconds`, the wall time the fit took when it was run, which
is the first fit's time when this call reused one and `NA` for a fit
cached before bmmtools recorded it, and `rng_seeded`, whether the fit
ran under a seed (`NA` for a fit cached before bmmtools recorded it; see
"Seed").

## Details

The fit and the key are each written to a temporary name and renamed
into place, so an interrupted write never leaves a truncated file under
the final name; the key is written after the fit. A cached fit that
cannot be read is reported as such rather than as a parse error.

How long the fit took is written to `<file>.meta.rds`, beside the fit
and outside the key. It is deliberately not a key component: a component
that changes with every run would differ from every stored key and refit
everything a resumed study had already fitted.

## Seed

A `seed` reaches the sampler and also R's random number generator: the
fitter is called under `withr::with_seed(seed)`, because bmm draws the
initial values of every chain from R's generator inside the fit. The
same `seed` therefore gives the same initial values and the same draws,
and the caller's random number stream is left as it was. This was
measured with both backends, each of which calls bmm's initial-value
function in the R session before it starts the chains (with brms's
`future = TRUE` each chain calls it inside its own future instead, which
was not measured): with cmdstanr on an `ezdm` cell with 1 or 4 `cores`,
and with rstan 2.32.7 and cmdstanr 0.9.0 on a `mixture2p` fit (10
subjects of 30 trials, 2 chains of 600 iterations, `cores = 2`) in two
fresh R sessions (without the seeded fitter, the two rstan runs differed
in initial values and draws). Without a `seed` (or with `NA`) nothing is
seeded and the initial values are drawn from the caller's stream. A
`seed` larger than `.Machine$integer.max` in absolute value is refused
with an error when the fit is run, because R's generator takes no larger
seed; a fit already cached under such a seed is still read back. Three
things the seed does not fix: within-chain threading through
[`brms::threading()`](https://paulbuerkner.com/brms/reference/threading.html)
with its default `static = FALSE`, whose chains differ between runs once
the data are split into more than one partial sum (`static = TRUE`
repeats); the caller's
[`RNGkind()`](https://rdrr.io/r/base/Random.html), which
[`withr::with_seed()`](https://withr.r-lib.org/reference/with_seed.html)
keeps, so a fit inside a future with `future.seed = TRUE`
(L'Ecuyer-CMRG) draws other initial values than the same fit run
sequentially under the default kind; and the random number arithmetic of
another platform or R version.

Fits cached before bmmtools seeded the initial values (0.3.0 and
earlier) still match their key and are read back as they were; the
attribute's `rng_seeded` is `NA` for them, and `TRUE` or `FALSE` for a
fit run since. Only the cache's meta file and this attribute record it:
grid sidecars,
[`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md)
results and
[`run_validation()`](https://www.gfrischkorn.org/bmmtools/reference/run_validation.md)
records do not, so a grid resumed over fits cached by 0.3.0 or earlier
mixes fits drawn both ways and its result cannot say which.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_cached(
  bmm::bmf(kappa ~ 1 + (1 | id), thetat ~ 1 + (1 | id)),
  data, bmm::mixture2p(resp_error = "y"),
  file = "fits/cell-01", seed = 1, chains = 4, iter = 2000
)
attr(fit, "bmmtools_cache")$reused
} # }
```
