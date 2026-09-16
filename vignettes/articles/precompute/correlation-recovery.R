# Precompute for vignettes/articles/correlation-recovery.Rmd.
#
# Needs bmm, brms, cmdstanr, CmdStan and lavaan. Run by hand from the
# package root:
#   Rscript vignettes/articles/precompute/correlation-recovery.R
# The article shows the labelled chunks below and loads the results from
# assets/correlation-recovery.rds, so it renders without Stan or lavaan.
# Fits go to precompute/fits/, which is not tracked.

devtools::load_all()
fits_dir <- "vignettes/articles/precompute/fits"
started <- Sys.time()

## ---- model
library(bmmtools)

model <- bmm::mixture2p(resp_error = "y")
pars <- c(kappa = log(8), thetat = qlogis(0.75))
sds <- c(kappa = 0.3, thetat = 0.5)

## ---- truth
loadings <- matrix(
  c(
    sqrt(0.7), sqrt(0.7), 0, 0, 0,
    0, 0, sqrt(0.7), sqrt(0.7), 0,
    0, 0, 0, 0, 1
  ),
  ncol = 3,
  dimnames = list(
    c("kappa_task1", "kappa_task2", "thetat_task1", "thetat_task2", "G"),
    c("precision", "memory", "ability")
  )
)
factor_cors <- function(rho) {
  phi <- matrix(rho, 3, 3, dimnames = rep(list(colnames(loadings)), 2))
  diag(phi) <- 1
  phi
}
cors_for <- function(row) cors_from_factors(loadings, factor_cors(row$rho))
round(cors_for(data.frame(rho = 0.5)), 2)

## ---- grid
grid <- expand.grid(
  rho = c(0, 0.5),
  re_cor = c("none", "all"),
  stringsAsFactors = FALSE
)
grid$n_subjects <- 100
grid$n_trials <- 100
grid

## ---- run
out <- recovery_grid(
  model, grid, pars,
  dir = file.path(fits_dir, "correlation-recovery"),
  reps = 20,
  sds = sds, cors = cors_for,
  covariates = list(G = c(mean = 0, sd = 1)),
  tasks = c("1", "2"),
  formula = function(row) {
    recovery_formula(model, re_cor = row$re_cor, task_col = "task")
  },
  correlations = c("model", "draws", "point"),
  seed = 2026,
  chains = 4, iter = 1000, cores = 4, threads = brms::threading(3),
  backend = "cmdstanr", refresh = 0, silent = 2
)

## ---- cells
cells <- attr(out, "cells")
# A resumed grid records read time in `elapsed`, so the fitting time is
# taken from the files: first simulation written to last fit written.
minutes_fitting <- as.double(difftime(
  max(file.mtime(cells$file)),
  min(file.mtime(sub("\\.rds$", "-sim.rds", cells$file))),
  units = "mins"
))
minutes_fitting

## ---- convergence
# The grid's own verdict uses check_convergence()'s defaults. At iter 1000
# many fits miss the usual bar (Rhat < 1.01, ESS > 400), mostly on the
# cor_ terms and a few of the 400 subject effects. The article scores the
# fits that pass the looser bar below and reports both verdicts.
bar <- list(rhat_max = 1.05, ess_bulk_min = 150, ess_tail_min = 150,
            divergent_max = 0)
strict_bar <- list(rhat_max = 1.01, ess_bulk_min = 400, ess_tail_min = 400,
                   divergent_max = 0)
convergence <- dplyr::bind_rows(lapply(seq_len(nrow(cells)), function(i) {
  fit <- readRDS(cells$file[[i]])
  verdict <- rlang::exec(check_convergence, fit, !!!bar)
  strict <- rlang::exec(check_convergence, fit, !!!strict_bar)
  tibble::tibble(
    condition = cells$condition[[i]],
    replication = cells$replication[[i]],
    verdict,
    pass_strict = strict$pass
  )
}))
kept <- convergence[convergence$pass, c("condition", "replication")]
dplyr::summarise(
  convergence,
  fits = dplyr::n(), pass = sum(pass), pass_strict = sum(pass_strict),
  .by = condition
)

## ---- convergence-groups
# Where the fits miss the strict bar: per fit and group of variables, how
# many variables do.
variable_group <- function(variable) {
  dplyr::case_when(
    startsWith(variable, "b_") ~ "population",
    startsWith(variable, "sd_") ~ "sd",
    startsWith(variable, "cor_") ~ "cor",
    startsWith(variable, "r_") ~ "subject"
  )
}
convergence_groups <- dplyr::bind_rows(lapply(
  seq_len(nrow(cells)), function(i) {
    draws <- posterior::subset_draws(
      posterior::as_draws_df(readRDS(cells$file[[i]])),
      variable = "^(b|sd|cor|r)_", regex = TRUE
    )
    diag <- posterior::summarise_draws(draws, "rhat", "ess_bulk", "ess_tail")
    # a variable that is constant across draws has no rhat
    diag <- diag[!is.na(diag$rhat), ]
    diag$group <- variable_group(diag$variable)
    dplyr::summarise(
      diag,
      condition = cells$condition[[i]],
      replication = cells$replication[[i]],
      n_variables = dplyr::n(),
      n_rhat = sum(rhat >= strict_bar$rhat_max),
      n_ess = sum(ess_bulk <= strict_bar$ess_bulk_min |
                    ess_tail <= strict_bar$ess_tail_min),
      max_rhat = max(rhat),
      min_ess_bulk = min(ess_bulk),
      .by = group
    )
  }
))
dplyr::summarise(
  convergence_groups,
  fits_rhat = sum(n_rhat > 0), fits_ess = sum(n_ess > 0),
  max_rhat = max(max_rhat), min_ess_bulk = min(min_ess_bulk),
  .by = c(condition, group)
)

## ---- correlations
design <- data.frame(condition = sprintf("row-%d", seq_len(nrow(grid))), grid)
cors_all <- attr(out, "correlations")
cors <- dplyr::semi_join(cors_all, kept, by = c("condition", "replication"))
cor_summary <- dplyr::left_join(summary(cors), design, by = "condition")
cor_summary

## ---- plot
plot_recovery(dplyr::filter(cors, condition == "row-4"), truth = "sample") +
  ggplot2::theme(strip.text = ggplot2::element_text(size = 7))

## ---- subject-table
tab <- subject_table(out, point = "median")
tab <- dplyr::semi_join(tab, kept, by = c("condition", "replication"))
tab

## ---- lavaan-grid
stopifnot(requireNamespace("lavaan", quietly = TRUE))
# Two indicators per factor are identified only through the factor's
# correlations with the rest of the model, which vanish at rho = 0, so the
# two tasks share a loading: they are parallel by construction. kappa and
# thetat of one task are estimated from the same trials, so their
# estimation errors covary; the residual covariances take that up.
cfa_model <- "
  precision =~ l1 * est_kappa_task1 + l1 * est_kappa_task2
  memory    =~ l2 * est_thetat_task1 + l2 * est_thetat_task2
  est_kappa_task1 ~~ est_thetat_task1
  est_kappa_task2 ~~ est_thetat_task2
  precision ~~ memory + G
  memory    ~~ G
"
# Every fit is kept with its status: warnings are recorded, not dropped.
factor_correlations <- function(data, model, factors) {
  warnings <- character()
  fit <- withCallingHandlers(
    tryCatch(
      lavaan::cfa(model, data = data, std.lv = TRUE),
      error = function(e) e
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  status <- if (inherits(fit, "error")) {
    "error"
  } else if (!lavaan::lavInspect(fit, "converged")) {
    "not converged"
  } else if (length(warnings) > 0) {
    "warning"
  } else {
    "ok"
  }
  if (status %in% c("error", "not converged")) {
    return(data.frame(
      pair = NA_character_, estimate = NA_real_, ci_low = NA_real_,
      ci_high = NA_real_, status = status,
      message = paste(c(conditionMessage(fit)[inherits(fit, "error")],
                        warnings), collapse = " | ")
    ))
  }
  sol <- lavaan::standardizedSolution(fit)
  sol <- sol[sol$op == "~~" & sol$lhs %in% factors & sol$lhs != sol$rhs, ]
  data.frame(
    pair = paste(sol$lhs, sol$rhs, sep = "__"),
    estimate = sol$est.std,
    ci_low = sol$ci.lower,
    ci_high = sol$ci.upper,
    status = status,
    message = paste(unique(trimws(warnings)), collapse = " | ")
  )
}
lavaan_grid <- dplyr::bind_rows(
  lapply(split(tab, ~ condition + replication, drop = TRUE), function(d) {
    fc <- factor_correlations(d, cfa_model, c("precision", "memory"))
    cbind(d[1, c("condition", "replication")], fc, row.names = NULL)
  })
)
lavaan_grid <- dplyr::left_join(lavaan_grid, design, by = "condition")
lavaan_status <- dplyr::count(
  dplyr::distinct(lavaan_grid, condition, replication, status),
  condition, status
)
lavaan_grid_summary <- dplyr::summarise(
  dplyr::filter(lavaan_grid, !is.na(pair)),
  n_fits = dplyr::n(),
  n_warning = sum(status == "warning"),
  mean_estimate = mean(estimate),
  sd_estimate = stats::sd(estimate),
  # every latent correlation is rho: G loads 1 on its own factor
  coverage = mean(ci_low <= rho & ci_high >= rho),
  .by = c(condition, rho, re_cor, pair)
)
lavaan_status
lavaan_grid_summary

## ---- components
a <- recovery_component(model, pars, n_trials = 100, sds = sds, name = "a")
b <- recovery_component(model, pars, n_trials = 100, sds = sds, name = "b")
set_loadings <- loadings
rownames(set_loadings) <- c("a_kappa", "b_kappa", "a_thetat", "b_thetat", "G")
set <- simulate_components(
  list(a, b),
  n_subjects = 100,
  cors = cors_from_factors(set_loadings, factor_cors(0.5)),
  covariates = list(G = c(mean = 0, sd = 1)),
  seed = 2026
)
set$truth$cor

## ---- fit-components
fits <- fit_components(
  set,
  dir = file.path(fits_dir, "correlation-recovery-components"),
  chains = 4, iter = 1000, cores = 4, threads = brms::threading(3),
  backend = "cmdstanr", refresh = 0, silent = 2
)
set_convergence <- dplyr::bind_rows(
  lapply(fits, function(fit) {
    verdict <- rlang::exec(check_convergence, fit, !!!bar)
    verdict$pass_strict <- rlang::exec(check_convergence, fit, !!!strict_bar)$pass
    verdict
  }),
  .id = "component"
)
set_convergence

## ---- separate-fits
set_recovery <- recover_correlations(fits, set, estimator = c("draws", "point"))
set_summary <- summary(set_recovery)
set_summary

## ---- reliability
set_tab <- subject_table(set, fits)
terms <- c("a_kappa", "b_kappa", "a_thetat", "b_thetat")
reliability <- vapply(
  terms,
  function(term) {
    true <- set_tab[[paste0("true_", term)]]
    est <- set_tab[[paste0("est_", term)]]
    stats::cor(true, est)^2
  },
  numeric(1)
)
pairs <- c("a_kappa__b_kappa", "a_thetat__b_thetat")
attenuation <- set_summary[set_summary$term %in% pairs, ]
attenuation$rel_a <- reliability[sub("__.*", "", attenuation$term)]
attenuation$rel_b <- reliability[sub(".*__", "", attenuation$term)]
attenuation$expected <- ifelse(
  attenuation$estimator == "point",
  attenuation$true_value * sqrt(attenuation$rel_a * attenuation$rel_b),
  attenuation$true_value * attenuation$rel_a * attenuation$rel_b
)
attenuation$disattenuated <- attenuation$mean_estimate / ifelse(
  attenuation$estimator == "point",
  sqrt(attenuation$rel_a * attenuation$rel_b),
  attenuation$rel_a * attenuation$rel_b
)
attenuation[, c("term", "estimator", "true_value", "mean_estimate",
                "rel_a", "rel_b", "expected", "disattenuated")]

## ---- lavaan-components
# Within each model, kappa and thetat come from the same fit and trials,
# so their estimation errors covary.
cfa_set <- "
  precision =~ est_a_kappa + est_b_kappa
  memory    =~ est_a_thetat + est_b_thetat
  est_a_kappa ~~ est_a_thetat
  est_b_kappa ~~ est_b_thetat
  precision ~~ memory + G
  memory    ~~ G
"
lavaan_set <- factor_correlations(set_tab, cfa_set, c("precision", "memory"))
lavaan_set

## ---- save
results <- list(
  grid = grid,
  design = design,
  cells = cells,
  bar = bar,
  strict_bar = strict_bar,
  convergence = convergence,
  convergence_groups = convergence_groups,
  minutes_fitting = minutes_fitting,
  cor_recovery = cors,
  cor_summary = cor_summary,
  lavaan_grid = lavaan_grid,
  lavaan_status = lavaan_status,
  lavaan_grid_summary = lavaan_grid_summary,
  set_truth = set$truth$cor,
  set_convergence = set_convergence,
  set_recovery = set_recovery,
  set_summary = set_summary,
  reliability = reliability,
  attenuation = attenuation,
  lavaan_set = lavaan_set,
  versions = c(
    bmm = as.character(packageVersion("bmm")),
    brms = as.character(packageVersion("brms")),
    cmdstanr = as.character(packageVersion("cmdstanr")),
    lavaan = as.character(packageVersion("lavaan"))
  ),
  minutes_grid = sum(cells$elapsed) / 60,
  minutes = as.double(difftime(Sys.time(), started, units = "mins"))
)
saveRDS(results, "vignettes/articles/assets/correlation-recovery.rds")
