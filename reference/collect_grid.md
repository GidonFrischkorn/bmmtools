# Read a finished recovery grid from its cell files

What
[`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md)
returns, rebuilt from the files it left in `dir`, without fitting
anything and without reading a fit. The per-cell sidecars hold
everything that was scored, so a study can run on one machine and be
assembled on another.

## Usage

``` r
collect_grid(dir, scale = NULL, levels = NULL, correlations = NULL)
```

## Arguments

- dir:

  The directory a
  [`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md)
  run wrote, holding `grid.rds` and the cell files. Written by bmmtools
  0.2.0 or later; an earlier run left no `grid.rds` and cannot be
  collected.

- scale:

  The scale to score on, as in
  [`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md).
  `NULL`, the default, uses the scale the grid was run on. Giving it
  re-scores the stored rows, so a grid run on the link scale can be read
  again on the natural one without refitting.

- levels, correlations:

  What to score, as in
  [`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md).
  `NULL` uses what the grid requested. Either may be narrowed to a
  subset of what the cells hold; asking for something they do not hold
  is an error, because a level is extracted from a fit and the fits are
  not read here.

## Value

A `bmmtools_recovery` with the attributes
[`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md)
gives it. `cells$elapsed` is `NA`: no cell was run, so none took any
time. `cells$fit_seconds` is the time each fit took when it was run,
which is what a runtime table wants and what survives the trip between
machines.

## Details

This is a reader, not a resume. It checks neither the cache key nor the
version of bmmtools that wrote the files, because on the machine that
collects the results neither can match: the key covers the installed
bmm, brms and Stan toolchain. A resume is what it has always been —
rerunning the identical
[`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md)
call in the same directory — and it is the resume, not this function,
that guarantees the files describe the call.

Cell files written before bmmtools stored the central 50 % interval are
read with `ci_low_50`, `ci_high_50` and `covered_50` set to `NA`, so
their `coverage_50` is `NA`, and a message says how many there were.
Rerunning the grid with the fits still in `dir` re-extracts them without
refitting.

## See also

[`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md),
which writes what this reads.

## Examples

``` r
if (FALSE) { # \dontrun{
# on the server
recovery_grid(model, grid, pars, dir = "runs/mixture2p", reps = 20)

# on a laptop, with only runs/mixture2p/*-sim.rds and *-est.rds copied
out <- collect_grid("runs/mixture2p")
summary(out)
collect_grid("runs/mixture2p", scale = "link")
} # }
```
