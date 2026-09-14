# bmmtools 0.0.0.9000

First release: tools to validate cognitive measurement models fitted with bmm.

* `simulate_recovery()` simulates data with known generating values from a bmm model.
* `recovery_formula()` builds a bmm formula with a random intercept for every free parameter.
* `recovery_grid()` runs a resumable parameter recovery study over a design grid.
* `fit_cached()` fits a model once and refits only when something that determines the fit changes.
* `check_convergence()` summarises rhat, ESS, divergences and tree-depth hits into a pass verdict.
* `extract_estimates()` returns population- and subject-level estimates of a fit as a tibble.
* `recover()` and `recover_subjects()` score estimates against the generating values.
* `summary()` of a recovery object reports bias, RMSE, coverage, interval width and correlations.
* `plot_recovery()` plots estimates against generating values.
* `prior_check()` summarises the prior predictive distribution on the scale of the response.
* `plot_prior_check()` plots prior predictive draws against the observed data.
* `inverse_link()` transforms values from the link scale for the twelve links bmm uses.
* `recovery_mixture2p` and `prior_check_sdt_yn` are example results that need no Stan.
* A pkgdown website with five articles: <https://www.gfrischkorn.org/bmmtools/>.
* Licensed under GPL (>= 2), compatible with bmm's GPL-2.
