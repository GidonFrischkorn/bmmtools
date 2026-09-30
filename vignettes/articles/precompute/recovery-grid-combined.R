# Precompute for the section "One grid, every estimand" of
# vignettes/articles/recovery-grid.Rmd.
#
# Needs bmm, brms, cmdstanr and CmdStan. Run by hand from the package
# root:
#   Rscript vignettes/articles/precompute/recovery-grid-combined.R
# The article shows the labelled chunks below and loads the results from
# assets/recovery-grid-combined.rds, so it renders without Stan. Fits go
# to precompute/fits/, which is not tracked.

devtools::load_all()
fits_dir <- "vignettes/articles/precompute/fits"
started <- Sys.time()

## ---- model
library(bmmtools)

model <- bmm::mixture2p(resp_error = "y")
tasks <- c("1", "2")
# cell values on the link scale: task 2 lowers precision and memory
pars <- c(
  kappa_task1 = log(8), kappa_task2 = log(4),
  thetat_task1 = qlogis(0.8), thetat_task2 = qlogis(0.6)
)
sds <- c(
  kappa_task1 = 0.3, kappa_task2 = 0.3,
  thetat_task1 = 0.5, thetat_task2 = 0.5
)
# a person precise in one task is precise in the other
cors <- diag(4)
dimnames(cors) <- rep(list(names(pars)), 2)
cors["kappa_task1", "kappa_task2"] <- 0.7
cors["kappa_task2", "kappa_task1"] <- 0.7
cors["thetat_task1", "thetat_task2"] <- 0.7
cors["thetat_task2", "thetat_task1"] <- 0.7

## ---- run
out <- recovery_grid(
  model,
  grid = data.frame(n_subjects = c(40, 100), n_trials = 100),
  pars = pars,
  dir = file.path(fits_dir, "recovery-grid-combined"),
  reps = 10,
  sds = sds, cors = cors, tasks = tasks,
  coding = "contrast", re_cor = "all",
  levels = c("population", "effect", "sd", "subject"),
  correlations = "model",
  seed = 2028,
  chains = 4, iter = 4000, warmup = 1000, cores = 4,
  threads = brms::threading(3),
  backend = "cmdstanr", refresh = 0, silent = 2
)

## ---- families
families <- split(summary(out), ~level)
names(families)
correlations <- summary(attr(out, "correlations"))

## ---- save
cells <- attr(out, "cells")
# A resumed grid records read time in `elapsed`; `fit_seconds` is each
# fit's own time, kept in its cell file, so it survives a resume.
minutes_fitting <- sum(cells$fit_seconds) / 60
results <- list(
  out = out,
  families = families,
  correlations = correlations,
  cells = cells,
  minutes_fitting = minutes_fitting,
  toolchain = c(
    bmm = as.character(utils::packageVersion("bmm")),
    brms = as.character(utils::packageVersion("brms")),
    cmdstan = as.character(cmdstanr::cmdstan_version())
  ),
  minutes = as.double(difftime(Sys.time(), started, units = "mins"))
)
saveRDS(results, "vignettes/articles/assets/recovery-grid-combined.rds")
