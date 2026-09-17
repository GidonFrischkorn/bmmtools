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
name.

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
  there is one cache, not two.

- .fitter:

  The fitting function,
  [`bmm::bmm()`](https://venpopov.com/bmm/reference/bmm.html) by
  default. Tests inject a stand-in so that no model is compiled.

## Value

The fit, with an attribute `bmmtools_cache`: a list of `file`, `key` and
`reused`.

## Details

The fit and the key are each written to a temporary name and renamed
into place, so an interrupted write never leaves a truncated file under
the final name; the key is written after the fit. A cached fit that
cannot be read is reported as such rather than as a parse error.

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
