# bmmtools: validate cognitive measurement models fitted with bmm

bmmtools checks a bmm model before it is used on real data: that it can
be estimated from data of a given size, and that its priors imply
plausible data. Every check is the same loop, simulate, fit, score.

## Simulate

[`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md)
draws a data set and the truth that produced it from a bmm model and
population values on the link scale.
[`recovery_formula()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_formula.md)
gives every free parameter a random intercept.

## Fit

[`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md)
fits once and reuses the saved fit while nothing that determines it has
changed.
[`check_convergence()`](https://www.gfrischkorn.org/bmmtools/reference/check_convergence.md)
is the convergence gate.
[`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md)
runs simulate, fit and score over a design grid with replications, one
file per cell.

## Score

[`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
and
[`recover_subjects()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
score population-level and person-level estimates against the truth;
[`summary()`](https://rdrr.io/r/base/summary.html) of the result gives
bias, RMSE, coverage, interval width and the correlations.
[`extract_correlations()`](https://www.gfrischkorn.org/bmmtools/reference/extract_correlations.md)
and
[`recover_correlations()`](https://www.gfrischkorn.org/bmmtools/reference/recover_correlations.md)
do the same for between-subject correlations.
[`prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/prior_check.md)
summarises prior-predictive draws on the observable scale.
[`plot_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/plot_recovery.md)
and
[`plot_prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/plot_prior_check.md)
draw both.

The articles on the package website walk through each step:
<https://www.gfrischkorn.org/bmmtools/>.

## See also

Useful links:

- <https://github.com/GidonFrischkorn/bmmtools>

- <https://www.gfrischkorn.org/bmmtools/>

- Report bugs at <https://github.com/GidonFrischkorn/bmmtools/issues>

## Author

**Maintainer**: Gidon T. Frischkorn <gfrischkorn@icloud.com>
([ORCID](https://orcid.org/0000-0002-5055-9764)) \[copyright holder\]

Authors:

- Gidon T. Frischkorn <gfrischkorn@icloud.com>
  ([ORCID](https://orcid.org/0000-0002-5055-9764)) \[copyright holder\]
