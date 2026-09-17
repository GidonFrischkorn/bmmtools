# Check that a model's implementation is calibrated

Simulation-based calibration (Talts et al. 2018; Modrák et al. 2023) is
the check that a model's likelihood, its Stan code and its
post-processing agree with each other. If they do, the rank of a
parameter drawn from the prior, among the posterior draws of a fit to
data simulated from that draw, is uniform over the simulations. A rank
histogram that is not flat says the implementation is wrong; a flat one
says nothing about whether the model is a good one.

## Usage

``` r
sbc(
  model,
  formula,
  data,
  prior = NULL,
  n_sims = 100,
  level = c("population", "sd"),
  generator = NULL,
  ...,
  seed = NULL,
  file = NULL,
  refit = c("on_change", "never", "always"),
  cores_per_fit = NULL,
  thin_ranks = NULL,
  keep_fits = FALSE,
  cache_mode = c("none", "results"),
  cache_location = NULL,
  .fitter = NULL
)
```

## Arguments

- model:

  A `bmmodel`, built with the column names `data` uses.

- formula:

  The `bmmformula` under check. Without a `generator`, every parameter
  formula must be intercept-only, with or without one `(1 | id)` or
  `(1 | p | id)` term — the shapes
  [`recovery_formula()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_formula.md)
  writes. Anything else is an error naming the term: the default
  generator maps a prior draw onto
  [`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md)'s
  arguments, and that map is exact for those shapes and guesswork for a
  covariate or a task factor. With a `generator` any formula works, as
  long as every parameter that has a group term has the same one,
  written as a bare column name.

- data:

  The design to simulate over, and only that: the subjects are the
  unique values of the grouping column and the trials are the rows each
  of them has. The response values are ignored, because the generator
  replaces them. Without a `generator`, every subject must have the same
  number of rows. The simulated data sets carry the same subject labels
  as `data`, which is what lets a subject truth meet its own `r_` draw.
  With a `generator`, `data` is handed to it as it is, columns and all.

- prior:

  A `brmsprior`, or `NULL` for bmm's defaults. This is the prior that is
  calibrated, so `NULL` calibrates bmm's own. A list is an error:
  [`prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/prior_check.md)
  is what compares prior sets.

- n_sims:

  How many data sets to simulate and fit. Twenty is enough to see a run
  through; a verdict needs enough ranks to read a histogram over about
  twenty bins, which is what the default is for. **It is also `n_sims`
  model fits**, so a default run is hours of Stan and belongs in a
  script rather than at a prompt.

- level:

  Which draws are ranked: `"population"` the `b_<par>_Intercept` draws,
  `"sd"` the `sd_<group>__<par>_Intercept` draws, `"cor"` the
  `cor_<group>__…` draws and `"subject"` the
  `r_<group>__<par>[<id>,Intercept]` draws, one per subject and varying
  parameter. `"population"` cannot be dropped. A level the formula
  cannot produce — `"sd"` or `"subject"` without a group term, `"cor"`
  without a correlated one — is dropped with a message.

  `"subject"` is off by default because it adds `n_subjects` times the
  number of varying parameters to the variables SBC ranks and plots:
  twenty subjects and two parameters are forty more rank histograms. It
  is the level that checks the partial pooling, and the one where a
  units error would show — see Details.

  `level` selects what is **ranked**, not what is **drawn**: everything
  the formula implies is always drawn from the prior and simulated from,
  because the ranks are only uniform when the data come from the joint
  prior. Holding the between-subject SDs at zero while the fitted model
  has a prior on them would take the population ranks down with it.

  With a `generator`, each level means **every** draw of that class for
  the model's free parameters, read off the prior fit rather than built
  from the formula: `"population"` every `b_<par>_*` (so a cell-means
  formula's `b_kappa_task1` and a covariate's `b_kappa_cond`), `"sd"`
  every `sd_<group>__<par>_*`, `"cor"` every `cor_<group>__*` and
  `"subject"` every `r_<group>__<par>[*]`. A level the prior fit has no
  draw of is dropped with a message.

- generator:

  `NULL` simulates each data set with
  [`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md),
  which is exact for the intercept-only formulas above and requires
  them. A function lifts that restriction: it is called once per data
  set as `generator(draws, data)`, where `draws` is that data set's
  prior draw as a **named numeric vector under the fit's own draw
  names** — every `b_`, `sd_`, `cor_` and `r_` draw of the free
  parameters, whether or not it is ranked, because a generator handed
  less than the fit puts a prior on would simulate from a point mass
  where the fit has a prior — and `data` is `data` as given. It returns
  the data frame to fit, which must carry the model's response columns
  and, when the formula has a group term, the grouping column with the
  labels of `data` if the subject level is ranked. bmmtools never
  interprets the row, so the truth names equal the draw names by
  construction. Note that brms's `r_` draws are deviations from the
  intercept: a generator that builds a subject's value adds them to the
  population value itself.

- ...:

  Passed to the fitter, for the prior fit and for every data set fit:
  `chains`, `iter`, `backend`, `init`, `control`. `cores`,
  `sample_prior`, `file`, `file_refit` and `file_compress` are refused,
  each with a message saying where it belongs.

- seed:

  Applied around the prior-draw subsample and around
  [`SBC::generate_datasets()`](https://hyunjimoon.github.io/SBC/reference/generate_datasets.html),
  so the same seed gives the same data sets, and passed to the fitter of
  the prior fit, where it enters
  [`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md)'s
  key. `NULL` leaves the random number generator alone and is recorded
  as `NA`. The data set fits take no seed: SBC's `future.seed` handles
  them, and one seed across fits would correlate them.

- file:

  Where to cache the **prior fit**, as in
  [`prior_check()`](https://www.gfrischkorn.org/bmmtools/reference/prior_check.md).
  `NULL` uses a temporary file. The data set fits are cached by SBC
  through `cache_mode`, not by
  [`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md).

- refit:

  Passed to
  [`fit_cached()`](https://www.gfrischkorn.org/bmmtools/reference/fit_cached.md)
  for the prior fit.

- cores_per_fit:

  Cores for each data set fit. `NULL` leaves
  [`SBC::compute_SBC()`](https://hyunjimoon.github.io/SBC/reference/compute_SBC.html)'s
  own default.

- thin_ranks:

  Thinning before ranking. `NULL` leaves SBC's default for a function
  backend.

- keep_fits:

  `TRUE` keeps every fit in the result, which is hundreds of MB for a
  default run. `FALSE` still keeps the ranks and the convergence
  diagnostics.

- cache_mode, cache_location:

  Passed to
  [`SBC::compute_SBC()`](https://hyunjimoon.github.io/SBC/reference/compute_SBC.html).
  `"results"` needs a `cache_location`.

- .fitter:

  The fitting function,
  [`bmm::bmm()`](https://venpopov.com/bmm/reference/bmm.html) by
  default. Used for the prior fit and inside the backend; tests inject a
  stand-in so that nothing is compiled.

## Value

[`SBC::compute_SBC()`](https://hyunjimoon.github.io/SBC/reference/compute_SBC.html)'s
`SBC_results`, unchanged apart from one attribute, `bmmtools_sbc`,
holding `model`, `n_sims`, `level`, `variables` (the draw names that
were ranked), `seed` (`NA` when none was given), `prior`
([`brms::prior_summary()`](https://mc-stan.org/rstantools/reference/prior_summary.html)
of the prior fit, which is what the draws came from), `layout`,
`diagnostics` and `generator` (`"simulate_recovery"` or `"user"`).

## Details

`sbc()` computes no ranks. It fits the prior once, turns its draws into
data sets through
[`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md),
builds the backend that fits each of them, and hands both to
[`SBC::compute_SBC()`](https://hyunjimoon.github.io/SBC/reference/compute_SBC.html).
What comes back is SBC's own `SBC_results`, so
[`SBC::plot_rank_hist()`](https://hyunjimoon.github.io/SBC/reference/plot_rank_hist.html),
[`SBC::plot_ecdf_diff()`](https://hyunjimoon.github.io/SBC/reference/ECDF-plots.html),
[`SBC::plot_coverage()`](https://hyunjimoon.github.io/SBC/reference/plot_coverage.html)
and `results$stats` all work as SBC documents them.

Two convergence bars apply to the same fits. bmmtools warns once on its
own Rhat 1.05 (see
[`check_convergence()`](https://www.gfrischkorn.org/bmmtools/reference/check_convergence.md))
together with SBC's rank-ESS 0.5, records both counts in the attribute,
and names SBC's stricter Rhat 1.01 in the message. An ESS below half the
maximum rank skews the ranks themselves, which is why it is a bar here
and not only a diagnostic.

The subject truths are **deviations**. brms's `r_` draws are deviations
from the intercept, while
[`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md)'s
`truth$subjects` are absolute values, so the generator emits each
subject's value minus that data set's population value, on the link
scale, under the label the subject has in `data`. The subject values
themselves are not read off the prior draw: they are drawn by
[`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md)
from the drawn SDs and correlations, which is the same conditional the
prior puts on them.

## References

Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A.
(2018). Validating Bayesian inference algorithms with simulation-based
calibration.
[doi:10.48550/arXiv.1804.06788](https://doi.org/10.48550/arXiv.1804.06788)

Modrák, M., Moon, A. H., Kim, S., Bürkner, P.-C., Huurre, N.,
Faltejsková, K., Gelman, A., & Vehtari, A. (2023). Simulation-based
calibration checking for Bayesian computation: The choice of test
quantities shapes sensitivity. *Bayesian Analysis*.
[doi:10.1214/23-BA1404](https://doi.org/10.1214/23-BA1404)

## Examples

``` r
if (FALSE) { # \dontrun{
model <- bmm::mixture2p(resp_error = "y")
results <- sbc(
  model,
  recovery_formula(model),
  data.frame(id = rep(1:20, each = 50), y = 0),
  n_sims = 20,
  chains = 2, iter = 500, backend = "cmdstanr",
  seed = 1
)
SBC::plot_rank_hist(results)
attr(results, "bmmtools_sbc")$variables

# any formula, with a generator that reads the draw row itself: here a
# condition effect on kappa, so `draws` carries `b_kappa_cond` too
design <- data.frame(
  id = rep(1:20, each = 50), cond = rep(c(-0.5, 0.5), 500), y = 0
)
generate <- function(draws, data) {
  # one call per subject and condition: rmixture2p() takes scalars
  cells <- split(seq_len(nrow(data)), list(data$id, data$cond))
  for (rows in cells) {
    id <- data$id[[rows[[1]]]]
    cond <- data$cond[[rows[[1]]]]
    kappa <- exp(
      draws[["b_kappa_Intercept"]] + draws[["b_kappa_cond"]] * cond +
        draws[[paste0("r_id__kappa[", id, ",Intercept]")]]
    )
    p_mem <- plogis(
      draws[["b_thetat_Intercept"]] +
        draws[[paste0("r_id__thetat[", id, ",Intercept]")]]
    )
    data$y[rows] <- bmm::rmixture2p(
      length(rows),
      kappa = kappa, p_mem = p_mem
    )
  }
  data
}
results <- sbc(
  model,
  bmm::bmf(kappa ~ 1 + cond + (1 | id), thetat ~ 1 + (1 | id)),
  design,
  generator = generate,
  n_sims = 20,
  chains = 2, iter = 500, backend = "cmdstanr"
)
} # }
```
