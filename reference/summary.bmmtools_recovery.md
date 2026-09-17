# Summarise a recovery object into per-parameter metrics

One row per parameter and level, with these recovery metrics: bias,
RMSE, coverage, mean interval width, the Pearson correlation with a
Fisher-z interval, the Spearman correlation, and Lin's concordance with
its interval, its decomposition and a calibration slope (see
[`recovery_ccc()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_ccc.md)
for how to read them).

## Usage

``` r
# S3 method for class 'bmmtools_recovery'
summary(object, ...)

# S3 method for class 'bmmtools_recovery_summary'
summary(object, ...)
```

## Arguments

- object:

  A `bmmtools_recovery` object from
  [`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
  or
  [`recover_subjects()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md).

- ...:

  Not used.

## Value

A `bmmtools_recovery_summary` tibble with the columns `term`,
`estimator`, `level`, `scale`, `n`, `n_replications`, `n_converged`,
`bias`, `rmse`, `coverage`, `ci_width`, `r`, `r_low`, `r_high`,
`rank_r`, `ccc`, `ccc_low`, `ccc_high`, `ccc_accuracy`,
`ccc_scale_shift`, `ccc_location_shift`, `calibration_slope` and
`truth_sd`, the standard deviation of the generating values. `r` and
`ccc` both grow with `truth_sd` at a fixed measurement error, so compare
them only at a similar spread.

## Details

Correlation metrics are `NA`, never `0`, when fewer than three complete
pairs are available or when either side has no spread: `0` would read as
"no recovery" where the honest answer is "not estimable from this
design". The correlation interval is a 95% confidence interval and does
not follow `ci_level`, which is the mass of the posterior interval and a
different quantity.

Population and SD rows are summarised across replications, one pair per
replication. Subject-level objects are summarised **within replication
and then combined**. Pooling subjects across replications would mix
between-subject with between-replication variance. `r` and `rank_r` are
combined on Fisher's z scale with weights `n - 3`; `ccc` on Lin's Z
scale with weights from his asymptotic variance, and without an interval
if any replication has no variance (three subjects, or a coefficient of
exactly 0 or 1); `ccc_scale_shift` and `calibration_slope` as geometric
means; `ccc_accuracy` and `ccc_location_shift` as means; `truth_sd` as
the root mean variance. At this level `ccc = r * ccc_accuracy` holds
only approximately.

Rows are grouped by `estimator` as well as by term and level, so two
estimators of the same parameter scored against one truth give two rows
rather than one pooled bias and RMSE. Their `n` may differ when one of
them failed on a subject the other estimated;
[`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
warns when it does.

`n_converged` is the number of replications whose fit passed
[`check_convergence()`](https://www.gfrischkorn.org/bmmtools/reference/check_convergence.md),
read from the `converged` column that
[`extract_estimates()`](https://www.gfrischkorn.org/bmmtools/reference/extract_estimates.md)
fills. It is `NA` when no fit carried a verdict, as with a hand-built
estimates tibble: reporting it as equal to `n` there would assert
something that was never measured.

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

# one cell of the example grid, subject level only
recovery_mixture2p |>
  dplyr::filter(condition == "row-4", level == "subject") |>
  summary()
#> # A tibble: 2 × 24
#>   condition term   estimator level   scale       n n_replications n_converged
#>   <chr>     <chr>  <chr>     <chr>   <chr>   <dbl>          <int>       <int>
#> 1 row-4     kappa  bayes     subject natural    50              5           3
#> 2 row-4     thetat bayes     subject natural    50              5           3
#>       bias   rmse coverage ci_width     r r_low r_high rank_r   ccc ccc_low
#>      <dbl>  <dbl>    <dbl>    <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>   <dbl>
#> 1 -0.259   1.69      0.932    5.76  0.751 0.690  0.802  0.776 0.676   0.615
#> 2  0.00164 0.0502    0.948    0.189 0.858 0.820  0.888  0.825 0.825   0.786
#>   ccc_high ccc_accuracy ccc_scale_shift ccc_location_shift calibration_slope
#>      <dbl>        <dbl>           <dbl>              <dbl>             <dbl>
#> 1    0.730        0.904           0.669            -0.126               1.11
#> 2    0.858        0.965           0.814             0.0206              1.05
#>   truth_sd
#>      <dbl>
#> 1   2.45  
#> 2   0.0956
```
