# Split a grid's cells into tasks for separate processes

The cells of a
[`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md)
run, numbered as the grid numbers them, cut into `n` runs of consecutive
cells for `recovery_grid(cells = )`. The grid orders its cells
replication by replication, so a task of at least `nrow(grid)` cells
spans every design row, and with `n` at most `reps` the tasks cost about
the same whatever the rows cost. With more tasks than replications, a
task holds only some rows and the tasks of the costly rows take longer.

## Usage

``` r
grid_tasks(grid, reps, n)
```

## Arguments

- grid:

  The grid that is passed to
  [`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md).

- reps:

  The replications that are passed to
  [`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md).

- n:

  The number of tasks, at most the number of cells.

## Value

A list of `n` integer vectors that hold every cell exactly once. Their
lengths differ by at most one, and the longer ones come first.

## See also

[`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md),
whose `cells` argument takes one element of the list.

## Examples

``` r
grid <- data.frame(n_subjects = c(20, 50), n_trials = c(30, 100))
grid_tasks(grid, reps = 10, n = 4)
#> [[1]]
#> [1] 1 2 3 4 5
#> 
#> [[2]]
#> [1]  6  7  8  9 10
#> 
#> [[3]]
#> [1] 11 12 13 14 15
#> 
#> [[4]]
#> [1] 16 17 18 19 20
#> 
```
