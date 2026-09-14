# Precompute for vignettes/articles/bmmtools.Rmd.
#
# Needs bmm, brms, cmdstanr and CmdStan. Run by hand from the package root:
#   Rscript vignettes/articles/precompute/bmmtools.R
# The article shows the labelled chunks below through knitr::read_chunk()
# and loads the results from assets/bmmtools.rds, so it renders without
# Stan. Fits go to precompute/fits/, which is not tracked.

devtools::load_all()
fits_dir <- "vignettes/articles/precompute/fits"
started <- Sys.time()

## ---- model
library(bmmtools)

model <- bmm::mixture2p(resp_error = "y")
pars <- c(kappa = log(8), thetat = qlogis(0.75))
sds <- c(kappa = 0.3, thetat = 0.5)

## ---- simulate
sim <- simulate_recovery(
  model,
  pars = pars, sds = sds,
  n_subjects = 40, n_trials = 100,
  seed = 1
)
sim$data
sim$truth$population

## ---- fit
fit <- fit_cached(
  recovery_formula(model), sim$data, model,
  file = file.path(fits_dir, "mixture2p-40x100"),
  seed = 1, chains = 4, iter = 2000, cores = 4, backend = "cmdstanr"
)

## ---- refit
fit <- fit_cached(
  recovery_formula(model), sim$data, model,
  file = file.path(fits_dir, "mixture2p-40x100"),
  seed = 1, chains = 4, iter = 2000, cores = 4, backend = "cmdstanr"
)
attr(fit, "bmmtools_cache")$reused

## ---- convergence
check_convergence(fit)

## ---- recover
population <- recover(fit, sim$truth$population)
population

## ---- subjects
subjects <- recover_subjects(fit, sim$truth$subjects)
summary(subjects)

## ---- plot
plot_recovery(subjects)

## ---- save
results <- list(
  data = sim$data,
  truth = sim$truth,
  reused = attr(fit, "bmmtools_cache")$reused,
  convergence = check_convergence(fit),
  population = population,
  subjects = subjects,
  minutes = as.double(difftime(Sys.time(), started, units = "mins"))
)
saveRDS(results, "vignettes/articles/assets/bmmtools.rds")
