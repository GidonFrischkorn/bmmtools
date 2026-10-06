# Example recovery object: a mixture2p grid

A precomputed
[`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md)
result, shipped so that
[`summary()`](https://rdrr.io/r/base/summary.html),
[`plot_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/plot_recovery.md)
and dplyr verbs can be tried without compiling a model. The grid crosses
20 and 50 subjects with 30 and 100 trials per subject, with 5
replications per cell, for
[`bmm::mixture2p()`](https://venpopov.com/bmm/reference/mixture2p.html)
with population values `kappa = log(8)` and `thetat = qlogis(0.75)` on
the link scale and between-subject SDs of 0.3 and 0.5. Every fit used 4
chains of 1000 iterations with cmdstanr, and the master seed was 2026.
Scores are on the natural scale.

## Usage

``` r
recovery_mixture2p
```

## Format

A `bmmtools_recovery` tibble at both the population and the subject
level, with the columns documented in
[`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
and a `condition` column naming the grid row. The attribute `cells`
holds one row per cell with its seed, status, runtime and convergence
verdict; `grid` holds the design.

## Source

`data-raw/example-objects.R` in the package's GitHub repository.

## See also

[`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md),
[`summary.bmmtools_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/summary.bmmtools_recovery.md)

## Examples

``` r
summary(recovery_mixture2p)
#> # A tibble: 16 × 16
#>    condition term   level      scale       n      bias    rmse coverage
#>    <chr>     <chr>  <chr>      <chr>   <dbl>     <dbl>   <dbl>    <dbl>
#>  1 row-1     kappa  population natural     5  0.427    1.45       0.8  
#>  2 row-1     thetat population natural     5  0.0349   0.0678     0.8  
#>  3 row-2     kappa  population natural     5  0.656    0.671      1    
#>  4 row-2     thetat population natural     5 -0.000662 0.00505    1    
#>  5 row-3     kappa  population natural     5  0.125    0.512      1    
#>  6 row-3     thetat population natural     5 -0.00646  0.0240     1    
#>  7 row-4     kappa  population natural     5  0.121    0.225      1    
#>  8 row-4     thetat population natural     5 -0.00299  0.00909    1    
#>  9 row-1     kappa  subject    natural    20 -0.275    2.38       0.87 
#> 10 row-1     thetat subject    natural    20  0.0210   0.0869     0.92 
#> 11 row-2     kappa  subject    natural    50  0.0974   2.39       0.924
#> 12 row-2     thetat subject    natural    50  0.00457  0.0712     0.92 
#> 13 row-3     kappa  subject    natural    20 -0.144    1.58       0.96 
#> 14 row-3     thetat subject    natural    20  0.00453  0.0500     0.97 
#> 15 row-4     kappa  subject    natural    50 -0.259    1.69       0.932
#> 16 row-4     thetat subject    natural    50  0.00164  0.0502     0.948
#>    coverage_50   z_mean  z_sd contraction      r     ccc calibration_slope
#>          <dbl>    <dbl> <dbl>       <dbl>  <dbl>   <dbl>             <dbl>
#>  1       0.2    0.217   1.73           NA NA     NA                 NA    
#>  2       0      0.655   1.52           NA NA     NA                 NA    
#>  3       0.2    1.10    0.290          NA NA     NA                 NA    
#>  4       1     -0.0321  0.284          NA NA     NA                 NA    
#>  5       0.6    0.246   0.714          NA NA     NA                 NA    
#>  6       0.6   -0.213   1.01           NA NA     NA                 NA    
#>  7       0.8    0.283   0.509          NA NA     NA                 NA    
#>  8       0.8   -0.157   0.603          NA NA     NA                 NA    
#>  9       0.37  -0.168   1.31           NA  0.581  0.0508             1.77 
#> 10       0.45   0.246   1.00           NA  0.757  0.587              0.825
#> 11       0.396  0.214   1.19           NA  0.605  0.255              1.59 
#> 12       0.48   0.00563 1.06           NA  0.660  0.453              0.982
#> 13       0.52  -0.0107  1.01           NA  0.856  0.799              1.11 
#> 14       0.52   0.0416  1.00           NA  0.889  0.854              1.05 
#> 15       0.5   -0.0448  1.05           NA  0.751  0.676              1.11 
#> 16       0.476 -0.0250  1.07           NA  0.858  0.825              1.05 
#>    calibration_intercept
#>                    <dbl>
#>  1               NA     
#>  2               NA     
#>  3               NA     
#>  4               NA     
#>  5               NA     
#>  6               NA     
#>  7               NA     
#>  8               NA     
#>  9               -6.84  
#> 10                0.0675
#> 11               -8.26  
#> 12               -0.0393
#> 13               -1.08  
#> 14               -0.0501
#> 15               -0.699 
#> 16               -0.0443
#> # 19 more columns; `tibble::as_tibble()` prints them all.
attr(recovery_mixture2p, "cells")
#> # A tibble: 20 × 20
#>    condition replication n_subjects n_trials  seed file           status elapsed
#>    <chr>           <int>      <int>    <int> <dbl> <chr>          <chr>    <dbl>
#>  1 row-1               1         20       30  9946 data-raw/fits… ok     0.00263
#>  2 row-2               1         50       30 17865 data-raw/fits… ok     0.0158 
#>  3 row-3               1         20      100 25784 data-raw/fits… ok     0.00222
#>  4 row-4               1         50      100 33703 data-raw/fits… ok     0.00248
#>  5 row-1               2         20       30  9947 data-raw/fits… ok     0.00243
#>  6 row-2               2         50       30 17866 data-raw/fits… ok     0.00284
#>  7 row-3               2         20      100 25785 data-raw/fits… ok     0.00207
#>  8 row-4               2         50      100 33704 data-raw/fits… ok     0.00237
#>  9 row-1               3         20       30  9948 data-raw/fits… ok     0.00208
#> 10 row-2               3         50       30 17867 data-raw/fits… ok     0.00226
#> 11 row-3               3         20      100 25786 data-raw/fits… ok     0.00200
#> 12 row-4               3         50      100 33705 data-raw/fits… ok     0.00218
#> 13 row-1               4         20       30  9949 data-raw/fits… ok     0.00204
#> 14 row-2               4         50       30 17868 data-raw/fits… ok     0.00212
#> 15 row-3               4         20      100 25787 data-raw/fits… ok     0.00214
#> 16 row-4               4         50      100 33706 data-raw/fits… ok     0.00206
#> 17 row-1               5         20       30  9950 data-raw/fits… ok     0.00193
#> 18 row-2               5         50       30 17869 data-raw/fits… ok     0.00213
#> 19 row-3               5         20      100 25788 data-raw/fits… ok     0.00212
#> 20 row-4               5         50      100 33707 data-raw/fits… ok     0.00213
#> # ℹ 12 more variables: converged <lgl>, max_rhat <dbl>, min_ess_bulk <dbl>,
#> #   min_ess_tail <dbl>, n_divergent <int>, n_max_treedepth <int>,
#> #   n_variables <int>, failed <chr>, fit_seconds <dbl>, chains <int>,
#> #   iter <int>, threads <int>
```
