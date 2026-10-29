# Time the compile step on its own

No fit records how long it took to compile, so a benchmark times a
compile of its own. `force = TRUE` is what makes the number a compile
rather than a cache hit: under cmdstanr a second call to the same code
reuses the executable in 0.02 s against 4.7 s for a compile (measured
2026-10-05), while rstan recompiles every time (23.3 s).

## Usage

``` r
benchmark_compile(
  code,
  backend = "cmdstanr",
  force = TRUE,
  ...,
  .compiler = NULL
)
```

## Arguments

- code:

  Stan code, as
  [`bmm::stancode()`](https://paulbuerkner.com/brms/reference/stancode.html)
  or
  [`brms::make_stancode()`](https://paulbuerkner.com/brms/reference/stancode.html)
  returns it.

- backend:

  `"cmdstanr"` or `"rstan"`.

- force:

  Recompile even when an executable is cached. Leaving it `TRUE` is the
  point; `FALSE` measures a reuse.

- ...:

  Passed to the backend's compiler.

- .compiler:

  The compiling function, the backend's by default. Tests inject a
  stand-in so that nothing is compiled.

## Value

A one-row tibble with `compile_seconds`, `backend`, `forced`,
`code_hash` and `exe`.

## See also

[`benchmark()`](https://www.gfrischkorn.org/bmmtools/reference/benchmark.md),
[`benchmark_metrics()`](https://www.gfrischkorn.org/bmmtools/reference/benchmark_metrics.md)

## Examples

``` r
if (FALSE) { # \dontrun{
code <- bmm::stancode(formula, data, model, backend = "cmdstanr")
benchmark_compile(code, backend = "cmdstanr")
} # }
```
