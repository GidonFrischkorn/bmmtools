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
