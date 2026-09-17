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
#> # A tibble: 16 × 24
#>    condition term   estimator level      scale       n n_replications
#>    <chr>     <chr>  <chr>     <chr>      <chr>   <dbl>          <int>
#>  1 row-1     kappa  bayes     population natural     5              5
#>  2 row-1     thetat bayes     population natural     5              5
#>  3 row-2     kappa  bayes     population natural     5              5
#>  4 row-2     thetat bayes     population natural     5              5
#>  5 row-3     kappa  bayes     population natural     5              5
#>  6 row-3     thetat bayes     population natural     5              5
#>  7 row-4     kappa  bayes     population natural     5              5
#>  8 row-4     thetat bayes     population natural     5              5
#>  9 row-1     kappa  bayes     subject    natural    20              5
#> 10 row-1     thetat bayes     subject    natural    20              5
#> 11 row-2     kappa  bayes     subject    natural    50              5
#> 12 row-2     thetat bayes     subject    natural    50              5
#> 13 row-3     kappa  bayes     subject    natural    20              5
#> 14 row-3     thetat bayes     subject    natural    20              5
#> 15 row-4     kappa  bayes     subject    natural    50              5
#> 16 row-4     thetat bayes     subject    natural    50              5
#>    n_converged      bias    rmse coverage ci_width      r  r_low r_high rank_r
#>          <int>     <dbl>   <dbl>    <dbl>    <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
#>  1           3  0.427    1.45       0.8     3.58   NA     NA     NA     NA    
#>  2           3  0.0349   0.0678     0.8     0.151  NA     NA     NA     NA    
#>  3           4  0.656    0.671      1       2.45   NA     NA     NA     NA    
#>  4           4 -0.000662 0.00505    1       0.0850 NA     NA     NA     NA    
#>  5           5  0.125    0.512      1       2.93   NA     NA     NA     NA    
#>  6           5 -0.00646  0.0240     1       0.112  NA     NA     NA     NA    
#>  7           3  0.121    0.225      1       1.63   NA     NA     NA     NA    
#>  8           3 -0.00299  0.00909    1       0.0645 NA     NA     NA     NA    
#>  9           3 -0.275    2.38       0.87    7.93    0.581  0.423  0.705  0.538
#> 10           3  0.0210   0.0869     0.92    0.292   0.757  0.650  0.834  0.671
#> 11           4  0.0974   2.39       0.924   8.49    0.605  0.517  0.679  0.615
#> 12           4  0.00457  0.0712     0.92    0.272   0.660  0.582  0.726  0.617
#> 13           5 -0.144    1.58       0.96    6.25    0.856  0.788  0.903  0.819
#> 14           5  0.00453  0.0500     0.97    0.197   0.889  0.835  0.926  0.830
#> 15           3 -0.259    1.69       0.932   5.76    0.751  0.690  0.802  0.776
#> 16           3  0.00164  0.0502     0.948   0.189   0.858  0.820  0.888  0.825
#>        ccc ccc_low ccc_high ccc_accuracy ccc_scale_shift ccc_location_shift
#>      <dbl>   <dbl>    <dbl>        <dbl>           <dbl>              <dbl>
#>  1 NA      NA       NA            NA              NA                NA     
#>  2 NA      NA       NA            NA              NA                NA     
#>  3 NA      NA       NA            NA              NA                NA     
#>  4 NA      NA       NA            NA              NA                NA     
#>  5 NA      NA       NA            NA              NA                NA     
#>  6 NA      NA       NA            NA              NA                NA     
#>  7 NA      NA       NA            NA              NA                NA     
#>  8 NA      NA       NA            NA              NA                NA     
#>  9  0.0508  0.0148   0.0867        0.474           0.261            -0.689 
#> 10  0.587   0.483    0.676         0.813           0.880             0.207 
#> 11  0.255   0.207    0.302         0.652           0.375             0.0871
#> 12  0.453   0.379    0.521         0.866           0.649             0.0546
#> 13  0.799   0.731    0.852         0.944           0.763            -0.0450
#> 14  0.854   0.796    0.897         0.962           0.840             0.0450
#> 15  0.676   0.615    0.730         0.904           0.669            -0.126 
#> 16  0.825   0.786    0.858         0.965           0.814             0.0206
#>    calibration_slope truth_sd
#>                <dbl>    <dbl>
#>  1            NA       0     
#>  2            NA       0     
#>  3            NA       0     
#>  4            NA       0     
#>  5            NA       0     
#>  6            NA       0     
#>  7            NA       0     
#>  8            NA       0     
#>  9             1.77    2.44  
#> 10             0.825   0.0942
#> 11             1.59    2.80  
#> 12             0.982   0.0917
#> 13             1.11    2.87  
#> 14             1.05    0.0990
#> 15             1.11    2.45  
#> 16             1.05    0.0956
attr(recovery_mixture2p, "cells")
#> # A tibble: 20 × 9
#>    condition replication n_subjects n_trials  seed file           status elapsed
#>    <chr>           <int>      <int>    <int> <dbl> <chr>          <chr>    <dbl>
#>  1 row-1               1         20       30  9946 data-raw/fits… ok       0.630
#>  2 row-2               1         50       30 17865 data-raw/fits… ok       0.656
#>  3 row-3               1         20      100 25784 data-raw/fits… ok       0.229
#>  4 row-4               1         50      100 33703 data-raw/fits… ok       0.595
#>  5 row-1               2         20       30  9947 data-raw/fits… ok       0.240
#>  6 row-2               2         50       30 17866 data-raw/fits… ok       0.682
#>  7 row-3               2         20      100 25785 data-raw/fits… ok       0.235
#>  8 row-4               2         50      100 33704 data-raw/fits… ok       0.502
#>  9 row-1               3         20       30  9948 data-raw/fits… ok       0.230
#> 10 row-2               3         50       30 17867 data-raw/fits… ok       0.501
#> 11 row-3               3         20      100 25786 data-raw/fits… ok       0.232
#> 12 row-4               3         50      100 33705 data-raw/fits… ok       0.679
#> 13 row-1               4         20       30  9949 data-raw/fits… ok       0.232
#> 14 row-2               4         50       30 17868 data-raw/fits… ok       0.499
#> 15 row-3               4         20      100 25787 data-raw/fits… ok       0.314
#> 16 row-4               4         50      100 33706 data-raw/fits… ok       0.494
#> 17 row-1               5         20       30  9950 data-raw/fits… ok       0.231
#> 18 row-2               5         50       30 17869 data-raw/fits… ok       0.601
#> 19 row-3               5         20      100 25788 data-raw/fits… ok       0.233
#> 20 row-4               5         50      100 33707 data-raw/fits… ok       0.600
#> # ℹ 1 more variable: converged <lgl>
```
