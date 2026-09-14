# bmmtools 0.0.0.9000

First release: tools to validate cognitive measurement models fitted with bmm.

* `simulate_recovery()` simulates data with known generating values from a bmm model. Subject values can be correlated (`cors`) and drawn together with observed covariates (`covariates`), and `pars`, `sds` and `cors` may be functions evaluated under the seed. The truth now also lists the between-subject SDs, every correlation pair and the covariate values. Without correlations or covariates, seeded simulations are unchanged. A covariate leaves the parameter values and the generated responses unchanged. Function-valued truths draw from the same random number stream and so change what follows.
* `cors_from_factors()` builds a correlation matrix from factor loadings.
* `recovery_formula()` builds a bmm formula with a random intercept for every free parameter, correlated across parameters with `re_cor = "all"`.
* `recovery_grid()` runs a resumable parameter recovery study over a design grid. Grid columns `cor_<a>__<b>` set a correlation per row, and `pars`, `sds` and `cors` may be functions of the row.
* `fit_cached()` fits a model once and refits only when something that determines the fit changes.
* `check_convergence()` summarises rhat, ESS, divergences and tree-depth hits into a pass verdict.
* `extract_estimates()` returns population- and subject-level estimates of a fit as a tibble. `level = "sd"` and `level = "cor"` add the between-subject standard deviations and correlations on the link scale, with correlations named `kappa__thetat` whichever order brms used.
* `recover()` and `recover_subjects()` score estimates against the generating values. `recover(level = c("population", "sd"))` also scores the standard deviations, always on the link scale, and takes the `truth` of a simulation as it is.
* `summary()` of a recovery object reports bias, RMSE, coverage, interval width and correlations, including Lin's concordance with a 95% interval, its accuracy factor, the scale and location shifts, a calibration slope and the spread of the generating values. At subject level the concordance is pooled across replications on Lin's Z scale; the geometric-mean pooling of the scale shift and slope has not been checked by simulation.
* `extract_subject_draws()` returns the per-draw subject values of a fit as an array of iterations, chains, subjects and parameters on the link scale.
* `extract_correlations()` estimates between-subject correlations in three ways: the model's own group-level correlation, the per-draw correlation across subjects and the correlation of the posterior means. Pairs can include observed covariates, and the last two estimators also work on the natural scale.
* `recover_correlations()` scores those correlations against the generating correlation and against the correlation the simulated subjects actually had. Its `summary()` reports bias, RMSE and coverage against both, the rate at which intervals exclude zero as a false-positive rate or power, and r and the concordance with the in-sample correlations.
* `recovery_grid()` writes a sidecar file per cell with everything it scores from the fit, so a resume reads the sidecars and the fit files can be deleted. `levels = "sd"` also scores the between-subject SDs, and `correlations` scores the between-subject correlations of every cell with `recover_correlations()`.
* `subject_table()` puts the true and estimated values of every subject side by side, one row per subject, for a simulation and its fit or for a whole grid. The table is the input for a structural equation model of true against estimated values.
* `plot_recovery()` is now a generic and also plots a correlation recovery, coloured by estimator.
* `recovery_ccc()` computes the same concordance columns for any pair of vectors, such as recovered correlations or effects.
* `plot_recovery()` plots estimates against generating values and, with `annotate = TRUE`, labels each panel with r and the concordance.
* `prior_check()` summarises the prior predictive distribution on the scale of the response.
* `plot_prior_check()` plots prior predictive draws against the observed data.
* `inverse_link()` transforms values from the link scale for the twelve links bmm uses.
* `recovery_mixture2p` and `prior_check_sdt_yn` are example results that need no Stan.
* A pkgdown website with five articles: <https://www.gfrischkorn.org/bmmtools/>.
* Licensed under GPL (>= 2), compatible with bmm's GPL-2.
