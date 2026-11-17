# Concordance between recovered and generating values

Lin's concordance correlation coefficient with its 95% interval and
decomposition, for any pair of vectors: subject estimates against their
generating values, or recovered correlations and effects against the
ones that were simulated.
[`summary()`](https://rdrr.io/r/base/summary.html) of a recovery object
reports the same columns per parameter.

## Usage

``` r
recovery_ccc(estimate, truth)
```

## Arguments

- estimate, truth:

  Numeric vectors of the same length. Incomplete pairs are dropped.

## Value

A one-row tibble with `ccc`, `ccc_low`, `ccc_high`, `ccc_accuracy`,
`ccc_scale_shift`, `ccc_location_shift`, `calibration_slope`,
`calibration_intercept` and `n`, the number of complete pairs. The
coefficient and its components are `NA` with fewer than three pairs or
when either vector has no spread.

## Details

The coefficient is
`1 - E[(estimate - truth)^2] / (sd_e^2 + sd_t^2 + (mean_e - mean_t)^2)`,
a rescaled mean squared deviation, and factors into Pearson's `r`
(precision) times `ccc_accuracy` (Lin, 1989, p. 258). Accuracy depends
on the scale shift `v = sd_e / sd_t` and the location shift
`u = (mean_e - mean_t) / sqrt(sd_e sd_t)`; all moments divide by `n`.

**Reading the scale shift.** Estimates from a hierarchical model are
shrunk towards the population mean, and a well-calibrated posterior mean
has `ccc_scale_shift` equal to `r`, not 1. `calibration_slope` (r
divided by `v`, the slope of truth regressed on the estimate) is 1 in
that case, above 1 when estimates are shrunk too much and below 1 when
they are shrunk too little. Because accuracy treats `v` and `1/v` alike,
`ccc` rewards too little shrinkage: read it together with the RMSE and
the calibration slope, and prefer `r` when only the rank order of
subjects matters.

**The calibration line.** `calibration_intercept` is the intercept of
the same regression, `mean_t - calibration_slope * mean_e`: the truth
the line predicts at an estimate of 0. It is read on the scale of the
input, so it depends on where 0 lies on that scale, and with the slope
it gives the whole line. It is not calibration-in-the-large
(`mean_t - mean_e`, the slope fixed at 1), which is minus the bias; the
two are related by
`calibration_intercept = (mean_t - mean_e) - (calibration_slope - 1) * mean_e`,
so at a slope of 1 they agree.

**The interval** uses Lin's Z transformation, `atanh(ccc)`, with the
asymptotic variance for random bivariate normal pairs, and is `NA` with
three or fewer pairs, at `ccc = +/-1` and when the Pearson correlation
is 0. When the generating values are fixed by design rather than drawn,
it is an approximation. It is always a 95% interval.

## References

Lin, L. I.-K. (1989). A concordance correlation coefficient to evaluate
reproducibility. *Biometrics, 45*(1), 255–268.
[doi:10.2307/2532051](https://doi.org/10.2307/2532051)

Lin, L. I.-K. (2000). A note on the concordance correlation coefficient.
*Biometrics, 56*(1), 324–325. <https://www.jstor.org/stable/2677159>

## Examples

``` r
truth <- c(0.2, 0.5, 0.9, 1.4, 1.8, 2.3)
# shrunk towards the mean and shifted up
estimate <- 0.7 * truth + 0.5
recovery_ccc(estimate, truth)
#> # A tibble: 1 × 9
#>     ccc ccc_low ccc_high ccc_accuracy ccc_scale_shift ccc_location_shift
#>   <dbl>   <dbl>    <dbl>        <dbl>           <dbl>              <dbl>
#> 1 0.915   0.794    0.966        0.915             0.7              0.238
#> # ℹ 3 more variables: calibration_slope <dbl>, calibration_intercept <dbl>,
#> #   n <int>
```
