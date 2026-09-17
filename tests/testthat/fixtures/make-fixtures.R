# Build the saved fixtures the score-layer tests read.
#
# This script is NOT run by the test suite. It compiles Stan and takes
# minutes; the suite reads its output and compiles nothing. Run it by
# hand, from the package root, when the fixture needs rebuilding:
#
#   Rscript tests/testthat/fixtures/make-fixtures.R
#
# A `set.seed()` here is deliberate and is not a violation of the
# "no set.seed() in package code" rule: this is a generator for a
# committed artefact, and the artefact has to be reproducible.

library(bmm)
library(brms)

set.seed(2026)

fixture_dir <- "tests/testthat/fixtures"
stopifnot(dir.exists(fixture_dir))

# Deliberately tiny. The fixture exists to exercise extract_estimates()'s
# handling of the *structure* of a bmm fit --- population terms, subject
# terms, fixed constants with NA diagnostics --- not to be a good fit.
n_subjects <- 8L
n_trials <- 30L

# Generating values on the natural scale, then the link scale, because
# the truth tibble contract (local/ARCHITECTURE.md decision 3) records truth on
# the link scale under bmm's parameter names.
kappa_pop <- 8
thetat_pop <- 0.75
kappa_sd <- 0.3 # on the log scale
thetat_sd <- 0.5 # on the logit scale

subject_kappa_link <- log(kappa_pop) + stats::rnorm(n_subjects, 0, kappa_sd)
subject_thetat_link <- stats::qlogis(thetat_pop) +
  stats::rnorm(n_subjects, 0, thetat_sd)

sim <- do.call(rbind, lapply(seq_len(n_subjects), function(i) {
  data.frame(
    id = factor(i, levels = seq_len(n_subjects)),
    y = bmm::rmixture2p(
      n = n_trials,
      mu = 0,
      kappa = exp(subject_kappa_link[i]),
      p_mem = stats::plogis(subject_thetat_link[i])
    )
  )
}))

model <- bmm::mixture2p(resp_error = "y")
formula <- bmm::bmf(
  kappa ~ 1 + (1 | id),
  thetat ~ 1 + (1 | id)
)

fit <- bmm::bmm(
  formula = formula,
  data = sim,
  model = model,
  # short on purpose: enough draws for summarise_draws() to return
  # finite rhat and ESS, few enough to keep the saved object small
  chains = 2,
  iter = 1000,
  warmup = 500,
  refresh = 0,
  backend = "cmdstanr",
  cores = 2
)

# The truth that goes with the fit, on the link scale under bmm's
# parameter names --- the contract recover() expects.
truth_population <- tibble::tibble(
  term = c("kappa", "thetat"),
  true_value = c(log(kappa_pop), stats::qlogis(thetat_pop))
)

truth_subjects <- tibble::tibble(
  id = rep(as.character(seq_len(n_subjects)), times = 2L),
  term = rep(c("kappa", "thetat"), each = n_subjects),
  true_value = c(subject_kappa_link, subject_thetat_link)
)

saveRDS(fit, file.path(fixture_dir, "mixture2p-fit.rds"), compress = "xz")
saveRDS(
  list(population = truth_population, subjects = truth_subjects),
  file.path(fixture_dir, "mixture2p-truth.rds"),
  compress = "xz"
)

size_mb <- file.size(file.path(fixture_dir, "mixture2p-fit.rds")) / 1024^2
cat(sprintf("mixture2p-fit.rds: %.3f MB\n", size_mb))
if (size_mb > 1) {
  cat(
    "Fixture exceeds the 1 MB target. Thin it with",
    "posterior::thin_draws() before committing.\n"
  )
}
cat("variables:\n")
print(posterior::variables(posterior::as_draws_array(fit)))

# ---- mixture2p-cor-draws.rds (Milestone 5.2) --------------------------
#
# A separate block with its own seed, so that rebuilding it does not touch
# the random stream of the fixture above and the fixture above does not
# have to be rebuilt to get it. To run it on its own, start R in the
# package root, read this file with readLines(), find this marker line
# with grep(), and eval(parse()) the lines from the marker to the end.
#
# It saves the draws, not the fit: the "sd" and "cor" levels read only
# the draws, `fit$ranef` and the link table, and a draws array is much
# smaller than a brmsfit.

library(bmm)

set.seed(20260914)

fixture_dir <- "tests/testthat/fixtures"
stopifnot(dir.exists(fixture_dir))

cor_n_subjects <- 20L
cor_n_trials <- 50L
cor_pop <- c(kappa = log(8), thetat = stats::qlogis(0.75))
cor_sd <- c(kappa = 0.3, thetat = 0.5)
cor_rho <- 0.5

cor_z <- matrix(stats::rnorm(cor_n_subjects * 2L), nrow = cor_n_subjects)
cor_x <- cor_z %*% chol(matrix(c(1, cor_rho, cor_rho, 1), 2L))
cor_kappa <- cor_pop[["kappa"]] + cor_sd[["kappa"]] * cor_x[, 1L]
cor_thetat <- cor_pop[["thetat"]] + cor_sd[["thetat"]] * cor_x[, 2L]

cor_data <- do.call(rbind, lapply(seq_len(cor_n_subjects), function(i) {
  data.frame(
    id = factor(i, levels = seq_len(cor_n_subjects)),
    y = bmm::rmixture2p(
      n = cor_n_trials,
      mu = 0,
      kappa = exp(cor_kappa[i]),
      p_mem = stats::plogis(cor_thetat[i])
    )
  )
}))

cor_backend <- if (
  requireNamespace("cmdstanr", quietly = TRUE) &&
    !inherits(try(cmdstanr::cmdstan_version(), silent = TRUE), "try-error")
) {
  "cmdstanr"
} else {
  "rstan"
}

cor_fit <- bmm::bmm(
  formula = bmm::bmf(
    kappa ~ 1 + (1 | p | id),
    thetat ~ 1 + (1 | p | id)
  ),
  data = cor_data,
  model = bmm::mixture2p(resp_error = "y"),
  chains = 2,
  iter = 500,
  warmup = 250,
  refresh = 0,
  backend = cor_backend,
  cores = 2,
  seed = 20260914
)

# Only the population, SD and correlation draws are kept. Keeping the 40
# subject deviations as well gave 0.162 MB against 0.019 MB (measured
# 2026-09-14), and the 5.2 tests do not read them.
cor_draws <- posterior::as_draws_array(cor_fit)
cor_keep <- grep("^(b|sd|cor)_", posterior::variables(cor_draws), value = TRUE)
cor_draws <- posterior::subset_draws(cor_draws, variable = cor_keep)

cor_file <- file.path(fixture_dir, "mixture2p-cor-draws.rds")
saveRDS(
  list(
    draws = cor_draws,
    ranef = cor_fit$ranef,
    links = cor_fit$bmm$model$links
  ),
  cor_file,
  compress = "xz"
)

cat(sprintf("backend: %s\n", cor_backend))
cor_size_mb <- file.size(cor_file) / 1024^2
cat(sprintf("mixture2p-cor-draws.rds: %.3f MB\n", cor_size_mb))
cat("variables:\n")
print(posterior::variables(cor_draws))

# ---- mixture2p-effect-fit.rds (Milestone 9.3f) -------------------------
#
# The first fixture built under coding = "contrast" and the first place
# any of Milestone 9.3's design meets real brms draws --- 9.3a-e were all
# built and tested against the mock fitter. The draw names the extractor
# side assumes (b_<par>_<task_col>1, r_id__<par>[i,<task_col>1],
# sd_id__<par>_<task_col>1) come from a make_standata() measurement
# recorded in the spec, not from a fit; this fixture is what pins them.
#
# A separate block with its own seed, so rebuilding it does not touch the
# random stream of the fixtures above, and generated through
# simulate_recovery()/recovery_formula() themselves (unlike the fixtures
# above, which are hand-rolled) because what is under test here is the
# extractor's reading of real draws, not the generation math --- that is
# already covered by test-contrasts.R's hand and Monte-Carlo checks.

devtools::load_all()
library(bmmtools)
library(bayestestR)

set.seed(20260917)

effect_fixture_dir <- "tests/testthat/fixtures"
stopifnot(dir.exists(effect_fixture_dir))

effect_n_subjects <- 10L
effect_n_trials <- 20L

effect_sim <- simulate_recovery(
  bmm::mixture2p(resp_error = "y"),
  pars = c(
    kappa_task1 = log(8), kappa_task2 = log(5),
    thetat_task1 = stats::qlogis(0.75), thetat_task2 = stats::qlogis(0.6)
  ),
  sds = c(kappa = 0.25, thetat = 0.35),
  n_subjects = effect_n_subjects, n_trials = effect_n_trials,
  tasks = c("1", "2"), task_col = "task",
  coding = "contrast", contrasts = bayestestR::contr.equalprior,
  seed = 20260917
)

effect_formula <- recovery_formula(
  effect_sim$model,
  task_col = "task", coding = "contrast"
)

effect_fit <- bmm::bmm(
  formula = effect_formula,
  data = effect_sim$data,
  model = effect_sim$model,
  chains = 2,
  iter = 400,
  warmup = 200,
  refresh = 0,
  backend = "cmdstanr",
  cores = 2,
  seed = 20260917
)

saveRDS(
  effect_fit,
  file.path(effect_fixture_dir, "mixture2p-effect-fit.rds"),
  compress = "xz"
)
saveRDS(
  effect_sim,
  file.path(effect_fixture_dir, "mixture2p-effect-sim.rds"),
  compress = "xz"
)

effect_fit_path <- file.path(effect_fixture_dir, "mixture2p-effect-fit.rds")
effect_size_mb <- file.size(effect_fit_path) / 1024^2
cat(sprintf("mixture2p-effect-fit.rds: %.3f MB\n", effect_size_mb))
if (effect_size_mb > 1) {
  cat(
    "Fixture exceeds the 1 MB target. Thin it with",
    "posterior::thin_draws() before committing.\n"
  )
}
cat("variables:\n")
print(posterior::variables(posterior::as_draws_array(effect_fit)))
