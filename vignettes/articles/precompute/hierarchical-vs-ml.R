# Precompute for vignettes/articles/hierarchical-vs-ml.Rmd.
#
# Needs bmm and brms throughout, and cmdstanr with a working CmdStan for
# every hierarchical fit and for the Stan route of fit_ml(). The optim
# route needs neither. Run by hand from the package root:
#   Rscript vignettes/articles/precompute/hierarchical-vs-ml.R
# The article shows the labelled chunks below through knitr::read_chunk()
# and loads the results from assets/hierarchical-vs-ml.rds, so it renders
# without Stan. Fits go to precompute/fits/, which is not tracked.
#
# The single data set is the one vignettes/articles/bmmtools.Rmd fits, at
# the same seed, the same scale and the same budget, so the hierarchical
# side is read back from the cache that article's precompute already
# wrote rather than sampled again.

devtools::load_all()
fits_dir <- "vignettes/articles/precompute/fits"
started <- Sys.time()

## ---- simulate
library(bmmtools)

model <- bmm::mixture2p(resp_error = "y")
pars <- c(kappa = log(8), thetat = qlogis(0.75))
sds <- c(kappa = 0.3, thetat = 0.5)

sim <- simulate_recovery(
  model,
  pars = pars, sds = sds,
  n_subjects = 40, n_trials = 100,
  seed = 1
)

## ---- hierarchical
fit <- fit_cached(
  recovery_formula(model), sim$data, model,
  file = file.path(fits_dir, "mixture2p-40x100"),
  seed = 1, chains = 4, iter = 2000, cores = 4, backend = "cmdstanr"
)
check_convergence(fit)

bayes <- extract_estimates(fit, level = "subject")

## ---- ml
ml <- fit_ml(
  model, sim$data,
  file = file.path(fits_dir, "ml-mixture2p-40x100"),
  seed = 1, refresh = 0, silent = 2
)
ml

## ---- bind
both <- dplyr::bind_rows(bayes, ml)
table(both$estimator, both$term)

## ---- score-link
scored_link <- recover_subjects(both, sim$truth$subjects, scale = "link")
summary(scored_link)

## ---- spread
spread <- dplyr::summarise(
  scored_link,
  sd_estimate = sd(estimate, na.rm = TRUE),
  sd_truth = sd(true_value, na.rm = TRUE),
  ratio = sd(estimate, na.rm = TRUE) / sd(true_value, na.rm = TRUE),
  .by = c(term, estimator)
)
spread

## ---- links
par_table <- bmm::parameters(model)
links <- stats::setNames(par_table$link, par_table$parameter)
links

## ---- score-natural
scored <- recover_subjects(both, sim$truth$subjects, links = links)
summary(scored)

## ---- plot
plot_recovery(scored_link, color_by = "estimator", annotate = TRUE)

## ---- population
recover(fit, sim$truth$population, level = "population")

## ---- optim
ml_optim <- fit_ml(model, sim$data, method = "optim")

summary(recover_subjects(
  dplyr::bind_rows(bayes, ml_optim), sim$truth$subjects, scale = "link"
))

## ---- grid
design <- expand.grid(n_subjects = 40, n_trials = c(25, 50, 100))

grid <- recovery_grid(
  model,
  grid = design,
  pars = pars, sds = sds,
  dir = file.path(fits_dir, "hierarchical-vs-ml"),
  reps = 3, seed = 2026,
  ml = list(method = "optim"),
  chains = 4, iter = 2000, cores = 4, backend = "cmdstanr",
  refresh = 0, silent = 2
)

## ---- grid-link
grid_link <- recovery_grid(
  model,
  grid = design,
  pars = pars, sds = sds,
  dir = file.path(fits_dir, "hierarchical-vs-ml"),
  reps = 3, seed = 2026,
  ml = list(method = "optim"),
  scale = "link",
  chains = 4, iter = 2000, cores = 4, backend = "cmdstanr",
  refresh = 0, silent = 2
)

## ---- grid-summary
cells <- attr(grid, "cells")
design_cols <- unique(cells[c("condition", "n_subjects", "n_trials")])

by_trials <- summary(grid_link) |>
  dplyr::filter(level == "subject") |>
  dplyr::left_join(design_cols, by = "condition")

by_trials[c("term", "estimator", "n_trials", "n", "n_converged",
            "bias", "rmse", "r", "coverage", "calibration_slope")]

## ---- grid-plot
library(ggplot2)

ggplot(by_trials, aes(factor(n_trials), rmse, colour = estimator)) +
  geom_point(size = 2.5) +
  geom_line(aes(group = estimator)) +
  facet_wrap(~term, scales = "free_y") +
  labs(x = "Trials per subject", y = "Subject-level RMSE (link scale)",
       colour = "Estimator")

## ---- ml-cells
attr(grid, "ml_cells")

## ---- save
# The cells' `elapsed` is the wall time of whatever run wrote this asset:
# once the fits are cached, a rerun only re-reads them and takes seconds,
# which says nothing about fitting. So the fitting time is taken from where
# it survives. A Stan fit keeps each chain's warmup and sampling time; the
# four chains ran in parallel, so a fit took as long as its slowest chain
# (compilation excluded). The optim route keeps no fit, so it is timed
# again here, on each cell's own simulated data.
stan_secs <- function(file) {
  max(rowSums(rstan::get_elapsed_time(readRDS(file)$fit)))
}
ml_secs <- function(file) {
  data <- readRDS(sub("\\.rds$", "-sim.rds", file))$data
  system.time(fit_ml(model, data, method = "optim"))[["elapsed"]]
}

results <- list(
  convergence = check_convergence(fit),
  reused = attr(fit, "bmmtools_cache")$reused,
  ml = ml,
  ml_optim = ml_optim,
  bound = table(both$estimator, both$term),
  links = links,
  scored_link = scored_link,
  scored = scored,
  spread = spread,
  optim_summary = summary(recover_subjects(
    dplyr::bind_rows(bayes, ml_optim), sim$truth$subjects, scale = "link"
  )),
  population = recover(fit, sim$truth$population, level = "population"),
  by_trials = by_trials,
  cells = cells,
  ml_cells = attr(grid, "ml_cells"),
  grid_secs = sum(vapply(cells$file, stan_secs, numeric(1))),
  grid_ml_secs = sum(vapply(cells$file, ml_secs, numeric(1))),
  minutes = as.double(difftime(Sys.time(), started, units = "mins"))
)
saveRDS(results, "vignettes/articles/assets/hierarchical-vs-ml.rds")
