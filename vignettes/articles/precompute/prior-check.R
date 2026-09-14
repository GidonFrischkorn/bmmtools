# Precompute for vignettes/articles/prior-check.Rmd.
#
# Needs bmm, brms, cmdstanr and CmdStan. Run by hand from the package root:
#   Rscript vignettes/articles/precompute/prior-check.R
# The article shows the labelled chunks below through knitr::read_chunk()
# and loads the results from assets/prior-check.rds, so it renders without
# Stan. Fits go to precompute/fits/, which is not tracked. The sdt_yn part
# of the article uses the shipped prior_check_sdt_yn, built in
# data-raw/example-objects.R.

devtools::load_all()
fits_dir <- "vignettes/articles/precompute/fits"

## ---- mixture-check
mix <- bmm::mixture2p(resp_error = "y")
mix_data <- simulate_recovery(
  mix,
  pars = c(kappa = log(8), thetat = qlogis(0.75)),
  sds = c(kappa = 0.3, thetat = 0.5),
  n_subjects = 10, n_trials = 40, seed = 1
)$data

abs_error <- function(yrep, data) {
  tibble::tibble(
    statistic = c("abs_q50", "abs_q90"),
    value = unname(stats::quantile(abs(yrep), c(0.5, 0.9)))
  )
}
mix_check <- prior_check(
  mix, recovery_formula(mix), mix_data,
  summary = abs_error,
  n_draws = 50, seed = 1,
  file = file.path(fits_dir, "prior-mixture2p"),
  chains = 2, iter = 1000, backend = "cmdstanr"
)
summary(mix_check)

## ---- mixture-plot
plot_prior_check(mix_check, draws = 30)

## ---- save
saveRDS(list(mix_check = mix_check), "vignettes/articles/assets/prior-check.rds")
