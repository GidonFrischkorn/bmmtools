# Precompute for vignettes/articles/recovery-ezdm.Rmd: the EZ diffusion
# model (3 parameters) on the standard grid of 20 or 50 subjects and 30,
# 100 or 300 trials, 10 replications each.
#
# Needs bmm (1.4.1 or later: its ezdm likelihood no longer assumes normal
# response times), brms, cmdstanr and CmdStan. Run by hand from the
# package root:
#   Rscript vignettes/articles/precompute/recovery-ezdm.R
# The article shows the labelled chunks below and loads the results from
# assets/recovery-ezdm.rds, so it renders without Stan. Fits go to
# precompute/fits/, which is not tracked.

devtools::load_all()
fits_dir <- "vignettes/articles/precompute/fits"
started <- Sys.time()

## ---- model
library(bmmtools)

model <- bmm::ezdm(
  mean_rt = "mean_rt", var_rt = "var_rt", n_upper = "n_upper",
  n_trials = "n_trials"
)
# population means on each parameter's link scale: drift rate 1.5
# (identity link), boundary separation exp(0) = 1 (log link),
# non-decision time 0.3 s (log link)
pars <- c(drift = 1.5, bound = 0, ndt = log(0.3))
# between-subject SDs of the person-level values, on the same scales
sds <- c(drift = 0.3, bound = 0.2, ndt = 0.2)
# between-subject correlations: people who drift faster also set wider
# boundaries (0.4); the other two pairs are uncorrelated
cors <- diag(3)
dimnames(cors) <- rep(list(names(pars)), 2)
cors["drift", "bound"] <- cors["bound", "drift"] <- 0.4

# The prior SD of each parameter under bmm's default priors, in closed
# form, for the contraction column: normal(., 0.5) on the bound and ndt
# intercepts, exponential(1) on sd(drift) and exponential(2) on the other
# two SDs. drift's intercept prior is cauchy(0, 1), which has no SD, so it
# has no row and its contraction is NA.
prior_sd <- data.frame(
  term = c("bound", "ndt", "drift", "bound", "ndt"),
  level = c("population", "population", "sd", "sd", "sd"),
  prior_sd_link = c(0.5, 0.5, 1, 0.5, 0.5)
)

## ---- run
out <- recovery_grid(
  model,
  grid = expand.grid(n_subjects = c(20, 50), n_trials = c(30, 100, 300)),
  pars = pars, sds = sds, cors = cors,
  dir = file.path(fits_dir, "recovery-ezdm"),
  reps = 10,
  re_cor = "all",
  levels = c("population", "sd", "subject"),
  correlations = "model",
  prior_sd = prior_sd,
  seed = 3131,
  chains = 4, iter = 2000, warmup = 1000, cores = 4,
  backend = "cmdstanr", refresh = 0, silent = 2
)

## ---- families
design <- attr(out, "grid")[c("n_subjects", "n_trials")]
design$condition <- sprintf("row-%d", seq_len(nrow(design)))
with_design <- function(x) {
  x <- merge(design, as.data.frame(x), by = "condition", sort = FALSE)
  x[order(x$n_subjects, x$n_trials), ]
}
families <- lapply(split(summary(out), ~level), with_design)
names(families)
correlations <- with_design(summary(attr(out, "correlations")))

## ---- convergence
cells <- attr(out, "cells")
# the package default (Rhat <= 1.05, bulk and tail ESS >= 400, at most
# 10 divergences, over every variable) and the stricter Rhat <= 1.01
cells$converged_101 <- cells$converged & cells$max_rhat <= 1.01
convergence <- with_design(aggregate(
  cbind(converged, converged_101, fit_seconds, min_ess_bulk) ~ condition,
  data = cells,
  FUN = function(x) c(mean = mean(x), min = min(x))
))

## ---- read-off
# Computed by hand from the stored rows, pending Milestone 24 (#17),
# which adds these as summary() columns and smallest_design().
reps <- length(unique(cells$replication))
# Monte Carlo SE of a rate from `reps` replications, at the observed rate
# and at the nominal 0.95 (the second is what D83's coverage leg uses)
mcse_rate <- function(p, n = reps) sqrt(p * (1 - p) / n)
coverage_mcse <- mcse_rate(0.95)
# the Wilson score interval's lower end of a rate
wilson_low <- function(p, n = reps, z = stats::qnorm(0.975)) {
  (p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
    (1 + z^2 / n)
}
population <- families$population
population$coverage_mcse_observed <- mcse_rate(population$coverage)
population$coverage_wilson_low <- wilson_low(population$coverage)
population$coverage_ok <-
  abs(population$coverage - 0.95) <= coverage_mcse + 1e-12

# Reliability of the subject estimates on the link scale, per
# replication, then averaged. With the truth: the squared correlation of
# posterior mean and generating value. From the model alone: the variance
# of the posterior means over itself plus the mean posterior variance.
# The generating value on the link scale is the posterior mean minus z
# posterior SDs.
subject_rows <- as.data.frame(out[out$level == "subject", ])
subject_rows$truth_link <- subject_rows$post_mean_link -
  subject_rows$z * subject_rows$post_sd_link
per_rep <- do.call(rbind, lapply(
  split(subject_rows, ~ condition + term + replication, drop = TRUE),
  function(d) {
    data.frame(
      condition = d$condition[[1]], term = d$term[[1]],
      true = stats::cor(d$post_mean_link, d$truth_link)^2,
      model = stats::var(d$post_mean_link) /
        (stats::var(d$post_mean_link) + mean(d$post_sd_link^2))
    )
  }
))
per_rep$gap <- per_rep$model - per_rep$true
reliability <- with_design(do.call(rbind, lapply(
  split(per_rep, ~ condition + term, drop = TRUE),
  function(d) {
    data.frame(
      condition = d$condition[[1]], term = d$term[[1]],
      reliability_true = mean(d$true),
      reliability_true_mcse = stats::sd(d$true) / sqrt(nrow(d)),
      reliability_model = mean(d$model),
      reliability_model_mcse = stats::sd(d$model) / sqrt(nrow(d)),
      reliability_gap = mean(d$gap),
      reliability_gap_mcse = stats::sd(d$gap) / sqrt(nrow(d))
    )
  }
)))

# D83: the smallest n_trials at which r_low clears the threshold for
# every parameter, and stays clear at every larger n_trials in the grid,
# with every population coverage within one MCSE of 0.95 at that design
smallest_trials <- function(threshold, subject = families$subject) {
  do.call(rbind, lapply(split(subject, ~n_subjects), function(s) {
    trials <- sort(unique(s$n_trials))
    meets <- function(col) {
      vapply(trials, function(t) {
        all(s[[col]][s$n_trials == t] >= threshold)
      }, logical(1))
    }
    coverage <- vapply(trials, function(t) {
      all(population$coverage_ok[
        population$n_subjects == s$n_subjects[[1]] &
          population$n_trials == t
      ])
    }, logical(1))
    first_stable <- function(ok) {
      stays <- rev(cumprod(rev(ok)) == 1)
      if (any(stays)) trials[which(stays)[[1]]] else NA_real_
    }
    bound_ok <- meets("r_low")
    data.frame(
      n_subjects = s$n_subjects[[1]],
      threshold = threshold,
      smallest = first_stable(bound_ok),
      smallest_point = first_stable(meets("r")),
      smallest_with_coverage = first_stable(bound_ok & coverage)
    )
  }))
}
readoff <- do.call(rbind, lapply(c(0.7, 0.8, 0.9), smallest_trials))
readoff

## ---- save
# A resumed grid records read time in `elapsed`; `fit_seconds` is each
# fit's own time, kept in its cell file, so it survives a resume.
minutes_fitting <- sum(cells$fit_seconds) / 60
# converged-only sensitivity: the same summary over the replications
# that passed the package's default convergence gate
families_converged <- lapply(
  split(summary(out[out$converged %in% TRUE, ]), ~level), with_design
)
results <- list(
  out = out,
  design = design,
  families = families,
  families_converged = families_converged,
  correlations = correlations,
  cells = cells,
  convergence = convergence,
  population = population,
  reliability = reliability,
  readoff = readoff,
  coverage_mcse = coverage_mcse,
  pars = pars, sds = sds, cors = cors, prior_sd = prior_sd,
  minutes_fitting = minutes_fitting,
  run_info = run_info(),
  toolchain = c(
    bmm = as.character(utils::packageVersion("bmm")),
    brms = as.character(utils::packageVersion("brms")),
    cmdstan = as.character(cmdstanr::cmdstan_version())
  ),
  minutes = as.double(difftime(Sys.time(), started, units = "mins"))
)
saveRDS(results, "vignettes/articles/assets/recovery-ezdm.rds")
# the variable each fit's lowest ESS and largest R-hat sit on, read from
# the fits just saved
source("vignettes/articles/precompute/recovery-ezdm-diagnostics.R")
