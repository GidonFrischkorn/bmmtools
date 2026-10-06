# The machine and package versions this session runs on

One row, meant for the runtime table of a study's parameter-recovery
page: what ran the fits, not what they estimated. Never errors — a
source it cannot read (no `/proc/cpuinfo`, no CmdStan installation,
`cmdstanr` not installed) degrades to `NA` rather than stopping the
call, because a runtime record is exactly the thing wanted after a run
that partly failed.

## Usage

``` r
run_info()
```

## Value

A one-row tibble: `r_version`, `platform`, `hostname`, `cores`
([`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html)),
`cpu_model`, `memory_gb`, the installed versions of bmm, brms, cmdstanr,
posterior and bmmtools, `cmdstan_version` (from
`cmdstanr::cmdstan_version()` when cmdstanr is installed and does not
error), and `timestamp`.

## Examples

``` r
run_info()
#> # A tibble: 1 × 13
#>   r_version platform hostname cores cpu_model memory_gb bmm_version brms_version
#>   <chr>     <chr>    <chr>    <int> <chr>         <dbl> <chr>       <chr>       
#> 1 4.6.1     x86_64-… runnerv…     4 AMD EPYC…      15.6 1.3.2       2.23.0      
#> # ℹ 5 more variables: cmdstanr_version <chr>, posterior_version <chr>,
#> #   bmmtools_version <chr>, cmdstan_version <chr>, timestamp <dttm>
```
