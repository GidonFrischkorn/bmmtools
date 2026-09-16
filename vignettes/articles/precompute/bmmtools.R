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

## ---- sbc-prior
sbc_prior <- brms::set_prior(
  "normal(2, 0.5)",
  class = "b", coef = "Intercept", nlpar = "kappa"
) +
  brms::set_prior("normal(0, 0.3)", class = "sd", nlpar = "kappa") +
  brms::set_prior("normal(0, 0.5)", class = "sd", nlpar = "thetat")

## ---- sbc
layout <- data.frame(id = rep(seq_len(20), each = 50), y = 0)

ranks <- sbc(
  model, recovery_formula(model), layout,
  prior = sbc_prior,
  n_sims = 20,
  level = c("population", "sd"),
  seed = 20260916,
  chains = 2, iter = 500, backend = "cmdstanr",
  cores_per_fit = 4, keep_fits = FALSE,
  file = file.path(fits_dir, "sbc-mixture2p")
)

## ---- sbc-plot
SBC::plot_rank_hist(ranks)

## ---- cross-check-fit
sdt <- bmm::sdt_yn(response = "hits", stimulus = "stim", n_trials = "n")
sdt_sim <- simulate_recovery(
  sdt,
  pars = c(d = 1, criterion = 0.3),
  sds = c(d = 0.3, criterion = 0.3),
  n_subjects = 20, n_trials = 50, seed = 2026
)
sdt_fit <- fit_cached(
  recovery_formula(sdt), sdt_sim$data, sdt,
  file = file.path(fits_dir, "sdt-yn-20x50"),
  seed = 1, chains = 4, iter = 2000, cores = 4, backend = "cmdstanr"
)

## ---- cross-check
observed <- sdt_sim$data
hit_rate <- sum(observed$hits[observed$stim == 1]) /
  sum(observed$n[observed$stim == 1])
fa_rate <- sum(observed$hits[observed$stim == 0]) /
  sum(observed$n[observed$stim == 0])

reference <- tibble::tibble(
  term = c("d", "criterion"),
  estimate = c(
    bmm::sdt_d(hit_rate = hit_rate, fa_rate = fa_rate),
    bmm::sdt_criterion(hit_rate = hit_rate, fa_rate = fa_rate)
  ),
  source = "closed form"
)

checked <- cross_check(sdt_fit, reference)
checked

## ---- cross-check-plot
plot_recovery(checked)

## ---- save
results <- list(
  data = sim$data,
  truth = sim$truth,
  reused = attr(fit, "bmmtools_cache")$reused,
  convergence = check_convergence(fit),
  population = population,
  subjects = subjects,
  ranks = ranks,
  sbc_info = attr(ranks, "bmmtools_sbc"),
  hit_rate = hit_rate,
  fa_rate = fa_rate,
  reference = reference,
  checked = checked,
  minutes = as.double(difftime(Sys.time(), started, units = "mins"))
)
saveRDS(results, "vignettes/articles/assets/bmmtools.rds")
