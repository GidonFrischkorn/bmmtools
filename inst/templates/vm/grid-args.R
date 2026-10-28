# The grid to run, written once.
#
# grid-slice.R sources this file in every process: in the setup call and
# in each slice. Only `cells` differs between them, and bmmtools adds it,
# so change the grid here and nowhere else. A slice whose call differs
# from the setup call's stops (see ?recovery_grid, "Running a grid in
# several processes").
#
# `dir` is read relative to the directory launch.sh is started from.

library(bmmtools)

grid_args <- list(
  model = bmm::mixture2p(resp_error = "y"),
  grid = expand.grid(n_subjects = c(20, 50), n_trials = c(30, 100)),
  pars = c(kappa = log(8), thetat = qlogis(0.75)),
  sds = c(kappa = 0.3, thetat = 0.5),
  dir = "fits/mixture2p",
  reps = 10,
  seed = 1,
  chains = 4,
  cores = 4,
  iter = 2000,
  backend = "cmdstanr"
)
