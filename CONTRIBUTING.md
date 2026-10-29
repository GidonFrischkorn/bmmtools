# Contributing to bmmtools

Bug reports and pull requests are welcome at
<https://github.com/GidonFrischkorn/bmmtools/issues>. For a bug, a
minimal example that runs without fitting a model helps most: an
estimates tibble and a truth tibble reproduce almost every scoring
problem.

## Adding a generator adapter for a bmm model

[`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md)
and
[`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md)
need a generator for each model. bmmtools ships adapters for a few bmm
models in `R/adapters.R`; any other model takes a `generator =` argument
(see the article “Validating a model without a built-in adapter”). To
add an adapter:

1.  Write `generate_<model>(pars, n_trials, model)` in `R/adapters.R`.
    `pars` is one subject’s parameters on the natural scale, fixed
    parameters included. Return one subject’s rows as a data frame and
    take every column name from the model object (`model$resp_vars`,
    `model$other_vars`), never from literals.
2.  Register it in `adapter_classes()` and `generator_for()`.
3.  If the response has a known floor and ceiling, add it to
    `response_range()` in `R/prior-check.R`, so
    [`prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/prior_check.md)
    reports floor and ceiling rates for the model.
4.  Add tests to `tests/testthat/test-generate.R`: the generated columns
    match the model’s, and the data pass bmm’s checks through the mock
    backend (`bmm::bmm(backend = "mock", mock_fit = 1)`).

## Tests compile no Stan

The test suite never compiles or samples a model. Scoring is tested on
estimates tibbles and the saved fixture in `tests/testthat/fixtures/`;
fitting is tested by injecting a stand-in through `.fitter`; generated
data are validated through bmm’s mock backend. Tests that need bmm use
`skip_if_not_installed("bmm")`.

## Package code

- No [`set.seed()`](https://rdrr.io/r/base/Random.html) in package code.
  Randomness takes a `seed` argument and is applied locally with
  [`withr::with_seed()`](https://withr.r-lib.org/reference/with_seed.html).
- bmm and brms stay in Suggests.
- Run `devtools::document()`, `devtools::test()` and `devtools::check()`
  before opening a pull request, and add a bullet to `NEWS.md` for any
  user-facing change.

## Articles and example objects

The articles in `vignettes/articles/` render from saved results, so the
website builds without Stan. Code shown in an article is read from the
script that produced its results: `vignettes/articles/precompute/*.R`
for the articles and `data-raw/example-objects.R` for the objects in
`data/`. After changing either, rerun the script from the package root
and rebuild the site with
[`pkgdown::build_site()`](https://pkgdown.r-lib.org/reference/build_site.html).
