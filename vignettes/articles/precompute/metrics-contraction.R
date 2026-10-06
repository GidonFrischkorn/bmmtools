# Precompute for the section "Coverage is not recovery" of the article
# on recovery summary columns (metrics.Rmd in vignettes/articles).
#
# Needs bmm, brms, cmdstanr and CmdStan. Run by hand from the package
# root:
#   Rscript vignettes/articles/precompute/metrics-contraction.R
# The article shows the labelled chunks below and loads the results from
# assets/metrics-contraction.rds, so it renders without Stan. Fits go to
# precompute/fits/, which is not tracked.
#
# The design is chosen to be weak: with one or three trials per subject
# the data say little about the population values and less about the
# subject SDs, so the posteriors range from prior-dominated to
# data-dominated across the five trial counts.

devtools::load_all()
fits_dir <- "vignettes/articles/precompute/fits"
started <- Sys.time()
# The one-trial row comes last because it was appended after the first
# run: a cell's seed follows its row, so appending kept the cached fits
# of the first four rows.

## ---- run
library(bmmtools)

model <- bmm::mixture2p(resp_error = "y")
pars <- c(kappa = log(8), thetat = qlogis(0.75))
sds <- c(kappa = 0.3, thetat = 0.5)

out <- recovery_grid(
  model,
  grid = data.frame(n_subjects = 20, n_trials = c(3, 10, 30, 100, 1)),
  pars = pars, sds = sds,
  dir = file.path(fits_dir, "metrics-contraction"),
  reps = 10,
  prior_sd = "analytic",
  levels = c("population", "sd"),
  scale = "link",
  seed = 2031,
  chains = 4, iter = 2000, cores = 4,
  backend = "cmdstanr", refresh = 0, silent = 2
)

## ---- prior-sd
prior_sd <- readRDS(
  file.path(fits_dir, "metrics-contraction", "prior-sd.rds")
)
prior_sd

## ---- summary
summary(out)[c(
  "condition", "term", "level", "n_converged", "coverage", "z_mean",
  "z_sd", "contraction"
)]

## ---- save
cells <- attr(out, "cells")
# A resumed grid records read time in `elapsed`; `fit_seconds` is each
# fit's own time, kept in its cell file, so it survives a resume.
minutes_fitting <- sum(cells$fit_seconds) / 60
results <- list(
  out = out,
  summary = summary(out),
  prior_sd = prior_sd,
  pars = pars,
  sds = sds,
  cells = cells,
  minutes_fitting = minutes_fitting,
  toolchain = c(
    bmm = as.character(utils::packageVersion("bmm")),
    brms = as.character(utils::packageVersion("brms")),
    cmdstan = as.character(cmdstanr::cmdstan_version())
  ),
  minutes = as.double(difftime(Sys.time(), started, units = "mins"))
)
saveRDS(results, "vignettes/articles/assets/metrics-contraction.rds")
