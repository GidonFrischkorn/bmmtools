# Summarise a recovery object into per-parameter metrics

One row per parameter and level, with these recovery metrics: bias,
RMSE, the coverage and mean width of the `ci_level` interval and of the
central 50 % interval, the Pearson correlation with a Fisher-z interval,
the Spearman correlation, and Lin's concordance with its interval, its
decomposition and a calibration line, slope and intercept (see
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
`bias`, `rmse`, `mae`, `coverage`, `coverage_50`, `ci_width`,
`ci_width_50`, `detected`, `sign_recovery`, `r`, `r_low`, `r_high`,
`rank_r`, `ccc`, `ccc_low`, `ccc_high`, `ccc_accuracy`,
`ccc_scale_shift`, `ccc_location_shift`, `calibration_slope`,
`calibration_intercept` and `truth_sd`, the standard deviation of the
generating values. `r` and `ccc` both grow with `truth_sd` at a fixed
measurement error, so compare them only at a similar spread.

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
means; `ccc_accuracy`, `ccc_location_shift` and `calibration_intercept`
as means (an intercept can be 0 or negative); `truth_sd` as the root
mean variance. `calibration_intercept` is read on the scale of the row:
the same fits scored on the link and on the natural scale give two
different intercepts, each the value of its own calibration line at 0.
At this level `ccc = r * ccc_accuracy` holds only approximately.

Rows are grouped by `estimator` as well as by term and level, so two
estimators of the same parameter scored against one truth give two rows
rather than one pooled bias and RMSE. Their `n` may differ when one of
them failed on a subject the other estimated;
[`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
warns when it does.

`detected` is the share of replications whose interval excludes zero and
`sign_recovery` the share whose estimate has the truth's sign. What
`detected` measures depends on the truth beside it: a false-positive
rate at a true value of zero, power at a true effect, and on a parameter
whose true effect is zero while another carries one, the rate at which
an effect leaks into a parameter that has none. `sign_recovery` is `NA`
when the truth is zero, which has no sign. Both are `NA` at the subject
level, where a row summarises a correlation across people rather than
one estimate against one truth.

**Two coverages.** `coverage` is the share of `ci_level` intervals (95 %
by default) that contain the generating value and `coverage_50` the
share of central 50 % intervals that do. A calibrated posterior has both
at their nominal mass; the 95 % interval checks the tails and the 50 %
interval the centre, so a posterior that is too narrow near its median
and too wide in its tails can be nominal on one and not on the other.
`coverage_50` is `NA` when the estimates carry no 50 % interval, as with
a hand-built tibble that has no `ci_low_50`/`ci_high_50`.

**The print.** Printing the summary shows a core set of columns — the
condition, the term, the estimator when there is more than one, the
level, the scale, `n`, `bias`, `rmse`, `coverage`, `coverage_50`, `r`,
`ccc` and the calibration columns — and says how many it hides. The
object itself keeps all of them;
[`tibble::as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
prints every one.

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
#> # A tibble: 16 × 13
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
#>    coverage_50      r     ccc calibration_slope calibration_intercept
#>          <dbl>  <dbl>   <dbl>             <dbl>                 <dbl>
#>  1       0.2   NA     NA                 NA                   NA     
#>  2       0     NA     NA                 NA                   NA     
#>  3       0.2   NA     NA                 NA                   NA     
#>  4       1     NA     NA                 NA                   NA     
#>  5       0.6   NA     NA                 NA                   NA     
#>  6       0.6   NA     NA                 NA                   NA     
#>  7       0.8   NA     NA                 NA                   NA     
#>  8       0.8   NA     NA                 NA                   NA     
#>  9       0.37   0.581  0.0508             1.77                -6.84  
#> 10       0.45   0.757  0.587              0.825                0.0675
#> 11       0.396  0.605  0.255              1.59                -8.26  
#> 12       0.48   0.660  0.453              0.982               -0.0393
#> 13       0.52   0.856  0.799              1.11                -1.08  
#> 14       0.52   0.889  0.854              1.05                -0.0501
#> 15       0.5    0.751  0.676              1.11                -0.699 
#> 16       0.476  0.858  0.825              1.05                -0.0443
#> # 17 more columns; `tibble::as_tibble()` prints them all.

# one cell of the example grid, subject level only
recovery_mixture2p |>
  dplyr::filter(condition == "row-4", level == "subject") |>
  summary()
#> # A tibble: 2 × 13
#>   condition term   level   scale       n     bias   rmse coverage coverage_50
#>   <chr>     <chr>  <chr>   <chr>   <dbl>    <dbl>  <dbl>    <dbl>       <dbl>
#> 1 row-4     kappa  subject natural    50 -0.259   1.69      0.932       0.5  
#> 2 row-4     thetat subject natural    50  0.00164 0.0502    0.948       0.476
#>       r   ccc calibration_slope calibration_intercept
#>   <dbl> <dbl>             <dbl>                 <dbl>
#> 1 0.751 0.676              1.11               -0.699 
#> 2 0.858 0.825              1.05               -0.0443
#> # 17 more columns; `tibble::as_tibble()` prints them all.
```
