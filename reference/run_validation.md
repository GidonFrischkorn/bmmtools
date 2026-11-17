# Run a validation plan

**\[experimental\]**

## Usage

``` r
run_validation(
  plan,
  dir,
  previous = NULL,
  refit = c("on_change", "never", "always"),
  workers = NULL,
  ...,
  .fitter = NULL
)
```

## Arguments

- plan:

  A `bmmtools_validation_plan`.

- dir:

  The directory every stage caches under, created if missing.

- previous:

  `NULL`, or an earlier result to compare with: a `bmmtools_validation`,
  its `result.rds`, or the directory it was run in. The comparison is
  the result's `changes`.

- refit:

  `"on_change"` (the default) reruns a stage when its inputs differ from
  those of its stored result, or when it has none; `"never"` reads what
  is stored and marks every stage without a stored result `"skipped"`;
  `"always"` reruns every stage, with `refit = "always"` for the prior
  fits of
  [`prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/prior_check.md)
  and [`sbc()`](https://www.gfrischkorn.org/bmmtools/reference/sbc.md).
  The recovery grid keeps its cells by their cache keys in every case
  (see
  [`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md));
  delete `<dir>/recovery/cells` to refit them.

- workers:

  Reserved for parallel recovery cells (#19). Only `NULL`, which runs
  the cells one after another, is accepted.

- ...:

  Passed to the fitter of every stage, after the plan's sampler lists,
  so an argument given here wins (`backend`, `iter`, `refresh`, ...).
  The arguments
  [`sbc()`](https://www.gfrischkorn.org/bmmtools/reference/sbc.md)
  refuses in its `...`, such as `cores`, do not reach the calibration. A
  `prior`, or an abbreviation of it such as `pri`, is refused, here and
  in the plan's sampler lists: the calibration's prior is the plan's
  `sbc$prior`, and the recovery stage fits under bmm's default priors.

- .fitter:

  The fitting function of every stage,
  [`bmm::bmm()`](https://venpopov.com/bmm/reference/bmm.html) by
  default. Tests inject a stand-in so that nothing is compiled.

## Value

A `bmmtools_validation`: a list with `plan`; `run_info`, the
[`run_info()`](https://www.gfrischkorn.org/bmmtools/reference/run_info.md)
rows written when the run started and ended, with `when` and the elapsed
`seconds`; `stages`, a tibble with one row per stage (`stage`, `order`,
`title`, `required`, `status`, `issue`, `reason`, `started`, `finished`,
`seconds`, `resumed`, `message`, `warnings`); one element per stage
holding its result object as the stage function returned it, or `NULL`;
`numbers`, the
[`card_numbers()`](https://www.gfrischkorn.org/bmmtools/reference/card_numbers.md)
tibble; `verdicts`, the
[`verdict_table()`](https://www.gfrischkorn.org/bmmtools/reference/verdict_table.md);
`headline`, a list with `state` (`"validated"`, `"not validated"`,
`"not judged"` or `"incomplete"`), `stages` (the stages the state names)
and `text`; and `changes`, the
[`compare_validations()`](https://www.gfrischkorn.org/bmmtools/reference/compare_validations.md)
of `previous` and this run, or `NULL`. The same object is written to
`<dir>/result.rds`.

## Details

Executes a
[`validation_plan()`](https://www.gfrischkorn.org/bmmtools/reference/validation_plan.md)
stage by stage into `dir`, where every stage keeps its result and its
status. A rerun reads the stages whose inputs have not changed and runs
the rest, so an interrupted run is resumed by calling it again.

## Stages

`status` is `"ok"`; `"failed"`, with the condition message in `message`;
`"skipped"`, when a stage it reads produced nothing or `refit = "never"`
found nothing stored; `"stale"`, when `refit = "never"` read a result
made from other inputs or software versions than this plan's (the result
is kept, but it is not this plan's); or the plan's own status for a
stage that does not run (`"not built"`, `"not applicable"`,
`"excluded"`). A failed stage is never an error of `run_validation()`.

The headline is `"incomplete"` while a required stage that exists is not
`"ok"`; `"not validated"` when a required criterion failed;
`"not judged"` when none failed but one could not be evaluated (a
response without a known range, say); and `"validated"`, printed with
the plan version and hash, only when every required criterion was
evaluated and met.

What each stage runs:

- `run_record`:
  [`run_info()`](https://www.gfrischkorn.org/bmmtools/reference/run_info.md),
  written to `<dir>/run-info.rds` before any other stage and again, with
  the end row, after the last.

- `prior_predictive`:
  [`prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/prior_check.md)
  on data simulated at the plan's prior-check design point, under bmm's
  default priors and, when the calibration has its own, under that one
  too.

- `sbc`:
  [`sbc()`](https://www.gfrischkorn.org/bmmtools/reference/sbc.md) at
  the plan's design point under the plan's prior.

- `recovery`:
  [`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md)
  over the plan's grid at the population, SD and subject levels, with
  the model correlations where there are at least two subject-varying
  parameters, every cell gated by
  [`check_convergence()`](https://www.gfrischkorn.org/bmmtools/reference/check_convergence.md)
  and contracted against a prior-only fit (`prior_sd = "fit"`).

- `computation`, `identifiability`, `individual_differences`: tables
  read off the recovery grid's cells and summary.

- `data_requirements`: the read-off,
  [`smallest_design()`](https://www.gfrischkorn.org/bmmtools/reference/smallest_design.md)
  of the subject-level `r_low` along `n_trials` at the plan's
  `read_off_n_subjects`, at `read_off_r_low` and at 0.70 and 0.90; one
  row per subject-varying parameter and threshold with `n_subjects`,
  `term`, `threshold`, `smallest` and `status`. It has no coverage
  condition of its own (see
  [`verdict_table()`](https://www.gfrischkorn.org/bmmtools/reference/verdict_table.md)).

- `agreement`:
  [`cross_check()`](https://www.gfrischkorn.org/bmmtools/reference/cross_check.md)
  of the fit of the largest design's first replication against the
  plan's reference.

## Resuming

A stage's inputs are the plan fields it reads, its seed, the `...` of
the call, the fitter, the versions of bmm, brms, bmmtools, the Stan
toolchain (and SBC for the calibration) and, for a stage read off the
grid, the grid's inputs. They are hashed into the stage's status file,
and a stage is read rather than rerun while its hash is unchanged and
its status is `"ok"`, so a new bmm release reruns every stage. A grid
with failed cells is rerun, which retries those cells and reads the
others from their cache, and a stage read off the grid is computed again
whenever the grid was not read. The whole plan is written to
`<dir>/plan.rds`; a plan whose hash differs from the stored one is
announced and replaces it.

## See also

[`validation_plan()`](https://www.gfrischkorn.org/bmmtools/reference/validation_plan.md),
[`card_numbers()`](https://www.gfrischkorn.org/bmmtools/reference/card_numbers.md),
[`verdict_table()`](https://www.gfrischkorn.org/bmmtools/reference/verdict_table.md),
[`compare_validations()`](https://www.gfrischkorn.org/bmmtools/reference/compare_validations.md).

## Examples

``` r
if (FALSE) { # \dontrun{
plan <- validation_plan("ezdm", seed = 20261024)
result <- run_validation(plan, "validation/ezdm")
result
verdict_table(result)
} # }
```
