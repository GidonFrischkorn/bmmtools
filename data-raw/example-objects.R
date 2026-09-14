# Build the two example objects shipped in data/.
#
# Needs bmm, brms, cmdstanr and a working CmdStan; run by hand from the
# package root with `Rscript data-raw/example-objects.R`. The fits land in
# data-raw/fits/, which is not tracked; only the scored results are saved.
# Paths passed to the package stay relative, so no local directory ends up
# in the shipped objects. The labelled chunks are shown in the articles
# through knitr::read_chunk().

devtools::load_all()

## ---- grid
library(bmmtools)

recovery_mixture2p <- recovery_grid(
  bmm::mixture2p(resp_error = "y"),
  grid = expand.grid(n_subjects = c(20, 50), n_trials = c(30, 100)),
  pars = c(kappa = log(8), thetat = qlogis(0.75)),
  sds = c(kappa = 0.3, thetat = 0.5),
  dir = "data-raw/fits/recovery", reps = 5, seed = 2026,
  chains = 4, iter = 1000, cores = 4, backend = "cmdstanr",
  refresh = 0, silent = 2
)

## ---- prior-check
sdt <- bmm::sdt_yn(response = "hits", stimulus = "stim", n_trials = "n")
sdt_data <- simulate_recovery(
  sdt,
  pars = c(d = 1, criterion = 0.3), sds = c(d = 0.3, criterion = 0.3),
  n_subjects = 20, n_trials = 50, seed = 2026
)$data

narrow_sd <- brms::set_prior("normal(0, 1)", class = "sd", dpar = "d") +
  brms::set_prior("normal(0, 1)", class = "sd", dpar = "criterion")

prior_check_sdt_yn <- prior_check(
  sdt, recovery_formula(sdt), sdt_data,
  prior = list(default = NULL, narrow_sd = narrow_sd),
  n_draws = 100, seed = 2026,
  file = "data-raw/fits/prior-check/sdt_yn",
  chains = 2, iter = 1000, backend = "cmdstanr", refresh = 0, silent = 2
)

## ---- save
# No absolute path may reach a shipped object.
home <- normalizePath("~")
for (obj in list(recovery_mixture2p, prior_check_sdt_yn)) {
  text <- paste(capture.output(str(attributes(obj))), collapse = "\n")
  stopifnot(!grepl(home, text, fixed = TRUE))
}

usethis::use_data(recovery_mixture2p, prior_check_sdt_yn,
                  compress = "xz", overwrite = TRUE)
