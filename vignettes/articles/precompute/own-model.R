# Precompute for vignettes/articles/own-model.Rmd.
#
# Needs bmm, brms, cmdstanr and CmdStan. Run by hand from the package root:
#   Rscript vignettes/articles/precompute/own-model.R
# The article shows the labelled chunks below through knitr::read_chunk()
# and loads the results from assets/own-model.rds, so it renders without
# Stan. Fits go to precompute/fits/, which is not tracked.

devtools::load_all()
fits_dir <- "vignettes/articles/precompute/fits"
started <- Sys.time()

## ---- model
library(bmmtools)

model <- bmm::mixture3p(
  resp_error = "y", nt_features = c("nt1", "nt2"), set_size = 3
)
model$links

## ---- generator
generate_mixture3p <- function(pars, n_trials, model) {
  # non-target locations relative to the target, one pair per trial
  nt <- matrix(stats::runif(2 * n_trials, -pi, pi), ncol = 2)
  # thetat and thetant are softmax weights against guessing, fixed at 0
  weights <- exp(c(pars$thetat, pars$thetant, 0))
  weights <- weights / sum(weights)
  y <- vapply(seq_len(n_trials), function(i) {
    bmm::rmixture3p(
      1,
      mu = c(pars$mu1, nt[i, ]), kappa = pars$kappa,
      p_mem = weights[[1]], p_nt = weights[[2]]
    )
  }, numeric(1))
  data.frame(y = y, nt1 = nt[, 1], nt2 = nt[, 2])
}

## ---- simulate
sim <- simulate_recovery(
  model,
  pars = c(kappa = log(8), thetat = 1.5, thetant = 0),
  sds = c(kappa = 0.3, thetat = 0.5),
  n_subjects = 30, n_trials = 100,
  generator = generate_mixture3p,
  seed = 1
)
sim$data
sim$truth$subjects

## ---- fit
fit <- fit_cached(
  recovery_formula(model), sim$data, model,
  file = file.path(fits_dir, "mixture3p-30x100"),
  seed = 1, chains = 4, iter = 2000, cores = 4, backend = "cmdstanr"
)
check_convergence(fit)

## ---- score
recover(fit, sim$truth$population)
summary(recover_subjects(fit, sim$truth$subjects))

## ---- save
results <- list(
  links = model$links,
  data = sim$data,
  truth = sim$truth,
  convergence = check_convergence(fit),
  population = recover(fit, sim$truth$population),
  subjects = recover_subjects(fit, sim$truth$subjects),
  minutes = as.double(difftime(Sys.time(), started, units = "mins"))
)
saveRDS(results, "vignettes/articles/assets/own-model.rds")
