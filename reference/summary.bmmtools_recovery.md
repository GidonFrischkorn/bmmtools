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
summary(object, precision = NULL, ...)

# S3 method for class 'bmmtools_recovery_summary'
summary(object, ...)
```

## Arguments

- object:

  A `bmmtools_recovery` object from
  [`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
  or
  [`recover_subjects()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md).

- precision:

  Optional targets for the width of the `ci_level` interval: a data
  frame with the columns `term`, `level` and `width`, and optionally
  `condition`, or a named numeric vector, which is shorthand for
  population-level targets (`c(kappa = 1)`). Each width is read on the
  row's own `scale`. Adds `precision_target`, `precision_met` and
  `precision_met_mcse`; see "Precision" below.

- ...:

  Not used; a misspelt argument is an error.

## Value

A `bmmtools_recovery_summary` tibble with the columns `term`,
`estimator`, `level`, `scale`, `n`, `n_replications`, `n_converged`,
`bias`, `rmse`, `mae`, `coverage`, `coverage_50`, `ci_width`,
`ci_width_50`, `detected`, `sign_recovery`, `r`, `r_low`, `r_high`,
`rank_r`, `ccc`, `ccc_low`, `ccc_high`, `ccc_accuracy`,
`ccc_scale_shift`, `ccc_location_shift`, `calibration_slope`,
`calibration_intercept`, `truth_sd`, `z_mean`, `z_sd`, `contraction`,
`z_mean_mcse`, `z_sd_mcse`, `reliability_true`, `reliability_true_mcse`,
`reliability_model`, `reliability_model_mcse`, `reliability_gap`,
`reliability_gap_mcse`, `coverage_mcse`, `coverage_50_mcse`,
`detected_mcse`, `sign_recovery_mcse` and `ci_level` (the mass the
coverages are nominal against, `NA` for a group that mixes masses),
followed, when `precision` is given, by `precision_target`,
`precision_met` and `precision_met_mcse`. A grid summary also has
`condition` first, followed by the grid's own columns (`n_subjects`,
`n_trials` or `n_trials_<component>`, and any other column of the grid),
so that each row says which design it describes; a grid column whose
name is already a summary column, such as a parameter called `bias`, is
renamed `grid_<name>`.
[`smallest_design()`](https://www.gfrischkorn.org/bmmtools/reference/smallest_design.md)
reads them. `truth_sd` is the standard deviation of the generating
values. `r` and `ccc` both grow with `truth_sd` at a fixed measurement
error, so compare them only at a similar spread.

**z and contraction.** `z_mean` and `z_sd` are the mean and SD of the
rows' posterior z-scores (see
[`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)),
pooled over every row of the group with no within-replication step,
since each z is already standardised. When the generating values are
drawn from the prior, a calibrated posterior gives 0 and 1, and a `z_sd`
above 1 says the posterior is too narrow for its error, below 1 too
wide. When they are held fixed, as in a
[`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md),
the prior's pull toward its centre enters both: a calibrated posterior
then has `z_mean` away from 0 for a value far from the prior's centre,
and `z_sd` below 1 in proportion to how much the prior still contributes
(for a conjugate normal, the posterior SD over the sampling error). The
reading is clean where contraction is near 1, and that is why the two
are reported together. `contraction` is the mean of the rows'
contractions. Read together they separate what coverage alone cannot: a
posterior that covers the truth because the data located it (contraction
near 1, `z_sd` near 1) from one that covers it because it never moved
from a prior that already contained it (contraction near 0).
`z_mean_mcse` is `z_sd / sqrt(n)` and `z_sd_mcse` is
`z_sd / sqrt(2 * (n - 1))`, with `n` the rows that have a z; both are
exact for independent normal z and understate when subject rows of one
replication share a population posterior or when z has heavy tails. With
fewer than two z, `z_sd` and both MCSEs are `NA`.

**Reliability of subject estimates.** Two estimates of how reliably a
design measures individual differences, on subject rows only (`NA`
elsewhere). Both are computed per replication on the link scale from the
posterior mean, whatever `scale` says, and then averaged over
replications. `reliability_true` is the squared correlation of the
posterior means with the generating values: with no change between
sessions, the test-retest correlation a second session would show
(Katahira et al., 2024). `reliability_model` needs no truth, so it is
the one real data can report: for posterior means
`var(mean) / (var(mean) + mean(sd^2))`, for an `"ml"` row the separation
form `1 - mean(se^2) / var(estimate)` (Adams, 2005), which is not
clamped at 0. Only the estimator `"ml"` gets the separation form; any
other label is read as a posterior mean. `reliability_gap` is the mean
of the paired difference, model minus truth, and says whether the
model's own number can be trusted for this design. Read it with
`reliability_gap_mcse`, not with the two marginal MCSEs: it is a paired
difference, so its MCSE is small where the two move together across
replications and can exceed theirs where they do not, and a gap smaller
than its MCSE is no evidence against the model's own number. A
replication with fewer than three complete subjects or no spread
contributes nothing, and each column uses the replications it can: one
whose generating values have no spread enters `reliability_model` but
neither `reliability_true` nor the gap, so the gap need not equal the
difference of the two columns. Likewise a natural-scale subject whose z
is `NA` (a zero posterior SD) has no link-scale truth and drops out of
`reliability_true` only. When subjects differ in how much data they
have, both weight the well-measured subjects more. Every `_mcse` here is
the SD across replications over the square root of their number, `NA`
with fewer than two.

**Monte Carlo error of the rates.** `coverage_mcse` and
`coverage_50_mcse` are binomial at the nominal mass,
`sqrt(p (1 - p) / R)` with `p` the `ci_level` or 0.5, so that a coverage
of 100 % still has an error to compare against. `detected_mcse` and
`sign_recovery_mcse` are binomial at the observed rate (Morris et al.,
2019), which is 0 at a rate of 0 or 1: ten successes in ten still leave
a 95 % interval reaching down to 0.72. At the subject level the coverage
MCSEs are the SD of the per-replication coverages over the square root
of their number, because the subjects of one fit share a population
posterior and are not independent intervals; that is the error of the
pooled coverage when every replication has the same number of subjects,
as a grid has. `detected_mcse` and `sign_recovery_mcse` are `NA` there,
with their rates. A group whose rows mix `ci_level`s has no nominal
mass, and `coverage_mcse` is `NA`.

**Precision.** With `precision`, `precision_met` is the share of
`ci_level` intervals no wider than `precision_target` (inclusive), and
`precision_met_mcse` its Monte Carlo standard error. The width is
`ci_high - ci_low` on the row's own `scale`, which is the one printed
beside it: natural for population and subject rows scored natural, link
for sd, effect and softmax rows. A natural-scale width is not translated
onto a link-scale row, because a width does not map through a nonlinear
link. On a population row the share is over replications and the MCSE
binomial at the observed rate. On a subject row it is over every subject
interval of every replication — the chance that one participant's
interval is this narrow — and the MCSE is the SD of the per-replication
shares over the square root of their number. Rows without a target have
`NA`. A term is named with its level, because one name can be a
population, an sd and a subject row at once; the named-vector shorthand
means the population row. An unnamed vector, a target for a term and
level the recovery does not have, a key given twice, and a width that is
not positive and finite are errors. A target applies to every estimator
of its term and level, and, without a `condition` column, to every
condition. A group whose rows mix scales or `ci_level`s has no share
(`NA`).

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
condition and, on a grid, the design columns after it, the term, the
estimator when there is more than one, the level, the scale, `n`,
`bias`, `rmse`, `coverage`, `coverage_50`, `z_mean`, `z_sd`,
`contraction`, `r`, `ccc`, `reliability_true` and `reliability_model`
when there are subject rows, the calibration columns, and
`precision_target` and `precision_met` when `precision` was given — and
says how many it hides. The object itself keeps all of them;
[`tibble::as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
prints every one.

`n_converged` is the number of replications whose fit passed
[`check_convergence()`](https://www.gfrischkorn.org/bmmtools/reference/check_convergence.md),
read from the `converged` column that
[`extract_estimates()`](https://www.gfrischkorn.org/bmmtools/reference/extract_estimates.md)
fills. It is `NA` when no fit carried a verdict, as with a hand-built
estimates tibble: reporting it as equal to `n` there would assert
something that was never measured.

## References

Adams, R. J. (2005). Reliability as a measurement design effect.
*Studies in Educational Evaluation, 31*(2–3), 162–172.
[doi:10.1016/j.stueduc.2005.05.008](https://doi.org/10.1016/j.stueduc.2005.05.008)

Katahira, K., Oba, T., & Toyama, A. (2024). Does the reliability of
computational models truly improve with hierarchical modeling? Some
recommendations and considerations for the assessment of model parameter
reliability. *Psychonomic Bulletin & Review, 31*(6), 2465–2486.
[doi:10.3758/s13423-024-02490-8](https://doi.org/10.3758/s13423-024-02490-8)

Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
studies to evaluate statistical methods. *Statistics in Medicine,
38*(11), 2074–2102.
[doi:10.1002/sim.8086](https://doi.org/10.1002/sim.8086)

## Examples

``` r
summary(recovery_mixture2p)
#> # A tibble: 16 × 20
#>    condition n_subjects n_trials term   level      scale       n      bias
#>    <chr>          <dbl>    <dbl> <chr>  <chr>      <chr>   <dbl>     <dbl>
#>  1 row-1             20       30 kappa  population natural     5  0.427   
#>  2 row-1             20       30 thetat population natural     5  0.0349  
#>  3 row-2             50       30 kappa  population natural     5  0.656   
#>  4 row-2             50       30 thetat population natural     5 -0.000662
#>  5 row-3             20      100 kappa  population natural     5  0.125   
#>  6 row-3             20      100 thetat population natural     5 -0.00646 
#>  7 row-4             50      100 kappa  population natural     5  0.121   
#>  8 row-4             50      100 thetat population natural     5 -0.00299 
#>  9 row-1             20       30 kappa  subject    natural    20 -0.275   
#> 10 row-1             20       30 thetat subject    natural    20  0.0210  
#> 11 row-2             50       30 kappa  subject    natural    50  0.0974  
#> 12 row-2             50       30 thetat subject    natural    50  0.00457 
#> 13 row-3             20      100 kappa  subject    natural    20 -0.144   
#> 14 row-3             20      100 thetat subject    natural    20  0.00453 
#> 15 row-4             50      100 kappa  subject    natural    50 -0.259   
#> 16 row-4             50      100 thetat subject    natural    50  0.00164 
#>       rmse coverage coverage_50   z_mean  z_sd contraction      r     ccc
#>      <dbl>    <dbl>       <dbl>    <dbl> <dbl>       <dbl>  <dbl>   <dbl>
#>  1 1.45       0.8         0.2    0.217   1.73           NA NA     NA     
#>  2 0.0678     0.8         0      0.655   1.52           NA NA     NA     
#>  3 0.671      1           0.2    1.10    0.290          NA NA     NA     
#>  4 0.00505    1           1     -0.0321  0.284          NA NA     NA     
#>  5 0.512      1           0.6    0.246   0.714          NA NA     NA     
#>  6 0.0240     1           0.6   -0.213   1.01           NA NA     NA     
#>  7 0.225      1           0.8    0.283   0.509          NA NA     NA     
#>  8 0.00909    1           0.8   -0.157   0.603          NA NA     NA     
#>  9 2.38       0.87        0.37  -0.168   1.31           NA  0.581  0.0508
#> 10 0.0869     0.92        0.45   0.246   1.00           NA  0.757  0.587 
#> 11 2.39       0.924       0.396  0.214   1.19           NA  0.605  0.255 
#> 12 0.0712     0.92        0.48   0.00563 1.06           NA  0.660  0.453 
#> 13 1.58       0.96        0.52  -0.0107  1.01           NA  0.856  0.799 
#> 14 0.0500     0.97        0.52   0.0416  1.00           NA  0.889  0.854 
#> 15 1.69       0.932       0.5   -0.0448  1.05           NA  0.751  0.676 
#> 16 0.0502     0.948       0.476 -0.0250  1.07           NA  0.858  0.825 
#>    reliability_true reliability_model calibration_slope calibration_intercept
#>               <dbl>             <dbl>             <dbl>                 <dbl>
#>  1           NA                NA                NA                   NA     
#>  2           NA                NA                NA                   NA     
#>  3           NA                NA                NA                   NA     
#>  4           NA                NA                NA                   NA     
#>  5           NA                NA                NA                   NA     
#>  6           NA                NA                NA                   NA     
#>  7           NA                NA                NA                   NA     
#>  8           NA                NA                NA                   NA     
#>  9            0.333             0.200             1.77                -6.84  
#> 10            0.522             0.511             0.825                0.0675
#> 11            0.384             0.282             1.59                -8.26  
#> 12            0.413             0.436             0.982               -0.0393
#> 13            0.710             0.628             1.11                -1.08  
#> 14            0.754             0.721             1.05                -0.0501
#> 15            0.601             0.567             1.11                -0.699 
#> 16            0.703             0.697             1.05                -0.0443
#> # 28 more columns; `tibble::as_tibble()` prints them all.

# the share of intervals on mean kappa at most 2 wide, per design
summary(recovery_mixture2p, precision = c(kappa = 2))
#> # A tibble: 16 × 22
#>    condition n_subjects n_trials term   level      scale       n      bias
#>    <chr>          <dbl>    <dbl> <chr>  <chr>      <chr>   <dbl>     <dbl>
#>  1 row-1             20       30 kappa  population natural     5  0.427   
#>  2 row-1             20       30 thetat population natural     5  0.0349  
#>  3 row-2             50       30 kappa  population natural     5  0.656   
#>  4 row-2             50       30 thetat population natural     5 -0.000662
#>  5 row-3             20      100 kappa  population natural     5  0.125   
#>  6 row-3             20      100 thetat population natural     5 -0.00646 
#>  7 row-4             50      100 kappa  population natural     5  0.121   
#>  8 row-4             50      100 thetat population natural     5 -0.00299 
#>  9 row-1             20       30 kappa  subject    natural    20 -0.275   
#> 10 row-1             20       30 thetat subject    natural    20  0.0210  
#> 11 row-2             50       30 kappa  subject    natural    50  0.0974  
#> 12 row-2             50       30 thetat subject    natural    50  0.00457 
#> 13 row-3             20      100 kappa  subject    natural    20 -0.144   
#> 14 row-3             20      100 thetat subject    natural    20  0.00453 
#> 15 row-4             50      100 kappa  subject    natural    50 -0.259   
#> 16 row-4             50      100 thetat subject    natural    50  0.00164 
#>       rmse coverage coverage_50   z_mean  z_sd contraction      r     ccc
#>      <dbl>    <dbl>       <dbl>    <dbl> <dbl>       <dbl>  <dbl>   <dbl>
#>  1 1.45       0.8         0.2    0.217   1.73           NA NA     NA     
#>  2 0.0678     0.8         0      0.655   1.52           NA NA     NA     
#>  3 0.671      1           0.2    1.10    0.290          NA NA     NA     
#>  4 0.00505    1           1     -0.0321  0.284          NA NA     NA     
#>  5 0.512      1           0.6    0.246   0.714          NA NA     NA     
#>  6 0.0240     1           0.6   -0.213   1.01           NA NA     NA     
#>  7 0.225      1           0.8    0.283   0.509          NA NA     NA     
#>  8 0.00909    1           0.8   -0.157   0.603          NA NA     NA     
#>  9 2.38       0.87        0.37  -0.168   1.31           NA  0.581  0.0508
#> 10 0.0869     0.92        0.45   0.246   1.00           NA  0.757  0.587 
#> 11 2.39       0.924       0.396  0.214   1.19           NA  0.605  0.255 
#> 12 0.0712     0.92        0.48   0.00563 1.06           NA  0.660  0.453 
#> 13 1.58       0.96        0.52  -0.0107  1.01           NA  0.856  0.799 
#> 14 0.0500     0.97        0.52   0.0416  1.00           NA  0.889  0.854 
#> 15 1.69       0.932       0.5   -0.0448  1.05           NA  0.751  0.676 
#> 16 0.0502     0.948       0.476 -0.0250  1.07           NA  0.858  0.825 
#>    reliability_true reliability_model calibration_slope calibration_intercept
#>               <dbl>             <dbl>             <dbl>                 <dbl>
#>  1           NA                NA                NA                   NA     
#>  2           NA                NA                NA                   NA     
#>  3           NA                NA                NA                   NA     
#>  4           NA                NA                NA                   NA     
#>  5           NA                NA                NA                   NA     
#>  6           NA                NA                NA                   NA     
#>  7           NA                NA                NA                   NA     
#>  8           NA                NA                NA                   NA     
#>  9            0.333             0.200             1.77                -6.84  
#> 10            0.522             0.511             0.825                0.0675
#> 11            0.384             0.282             1.59                -8.26  
#> 12            0.413             0.436             0.982               -0.0393
#> 13            0.710             0.628             1.11                -1.08  
#> 14            0.754             0.721             1.05                -0.0501
#> 15            0.601             0.567             1.11                -0.699 
#> 16            0.703             0.697             1.05                -0.0443
#>    precision_target precision_met
#>               <dbl>         <dbl>
#>  1                2             0
#>  2               NA            NA
#>  3                2             0
#>  4               NA            NA
#>  5                2             0
#>  6               NA            NA
#>  7                2             1
#>  8               NA            NA
#>  9               NA            NA
#> 10               NA            NA
#> 11               NA            NA
#> 12               NA            NA
#> 13               NA            NA
#> 14               NA            NA
#> 15               NA            NA
#> 16               NA            NA
#> # 29 more columns; `tibble::as_tibble()` prints them all.

# one cell of the example grid, subject level only
recovery_mixture2p |>
  dplyr::filter(condition == "row-4", level == "subject") |>
  summary()
#> # A tibble: 2 × 20
#>   condition n_subjects n_trials term   level   scale       n     bias   rmse
#>   <chr>          <dbl>    <dbl> <chr>  <chr>   <chr>   <dbl>    <dbl>  <dbl>
#> 1 row-4             50      100 kappa  subject natural    50 -0.259   1.69  
#> 2 row-4             50      100 thetat subject natural    50  0.00164 0.0502
#>   coverage coverage_50  z_mean  z_sd contraction     r   ccc reliability_true
#>      <dbl>       <dbl>   <dbl> <dbl>       <dbl> <dbl> <dbl>            <dbl>
#> 1    0.932       0.5   -0.0448  1.05          NA 0.751 0.676            0.601
#> 2    0.948       0.476 -0.0250  1.07          NA 0.858 0.825            0.703
#>   reliability_model calibration_slope calibration_intercept
#>               <dbl>             <dbl>                 <dbl>
#> 1             0.567              1.11               -0.699 
#> 2             0.697              1.05               -0.0443
#> # 28 more columns; `tibble::as_tibble()` prints them all.
```
