# Adds to assets/recovery-ezdm.rds, for every fit, the variable at which
# the convergence gate's numbers sit: the lowest bulk ESS, the lowest tail
# ESS and the largest R-hat, and the parameter with the lowest bulk ESS
# (lp__ and lprior left out). The article uses the two bulk-ESS columns;
# the tail-ESS and R-hat variables are recorded for readers of the asset. It
# reads the fits that recovery-ezdm.R saved and fits nothing.
# recovery-ezdm.R sources it at its end, and it then reads the fits under
# that script's `fits_dir`, in its recovery-ezdm subfolder. To add the
# columns to a stored asset, run it alone from the package root; if the
# fits were moved, give the recovery-ezdm folder itself (the one holding
# the cell-*.rds files, not its parent `fits`) as the argument:
#   Rscript vignettes/articles/precompute/recovery-ezdm-diagnostics.R \
#     [fits/recovery-ezdm folder]
#
# Needs bmm, brms, cmdstanr, posterior and rstan installed: readRDS() of a
# stored fit (fitted with the cmdstanr backend) loads the bmm, brms,
# cmdstanr and posterior namespaces, and reading its draws loads rstan. It
# compiles and samples nothing.

local({
  asset <- "vignettes/articles/assets/recovery-ezdm.rds"
  fits <- commandArgs(trailingOnly = TRUE)
  # `fits` is the recovery-ezdm folder; `fits_dir` and `root` its parent
  if (length(fits) != 1L) {
    # sourced by recovery-ezdm.R: its fits_dir; run alone: the default
    root <- if (exists("fits_dir", mode = "character")) {
      fits_dir
    } else {
      "vignettes/articles/precompute/fits"
    }
    fits <- file.path(root, "recovery-ezdm")
  }
  res <- readRDS(asset)
  cells <- res$cells

  # The gate summarises every variable of a fit and drops those without an
  # R-hat (constants); the same set is read here, so the numbers found
  # must be the ones the gate stored.
  weakest <- function(file) {
    fit <- readRDS(file)
    d <- posterior::summarise_draws(
      posterior::as_draws_array(fit),
      posterior::default_convergence_measures()
    )
    d <- d[!is.na(d$rhat), ]
    # the log posterior and log prior are not parameters of the model
    p <- d[!d$variable %in% c("lp__", "lprior"), ]
    data.frame(
      n_variables_check = nrow(d),
      min_ess_bulk_check = min(d$ess_bulk),
      min_ess_tail_check = min(d$ess_tail, na.rm = TRUE),
      max_rhat_check = max(d$rhat),
      min_ess_bulk_variable = d$variable[which.min(d$ess_bulk)],
      min_ess_bulk_parameter = p$variable[which.min(p$ess_bulk)],
      min_ess_tail_variable = d$variable[which.min(d$ess_tail)],
      max_rhat_variable = d$variable[which.max(d$rhat)]
    )
  }
  found <- do.call(
    rbind, lapply(file.path(fits, basename(cells$file)), weakest)
  )
  # The same fits and the same variable set as the gate: the count of
  # variables exactly; the ESS minima to all.equal()'s relative tolerance,
  # since an ESS runs to the thousands; R-hat, which sits near 1, to an
  # absolute 1e-8 in every fit.
  stopifnot(
    "fits read do not match the gate's stored values: n_variables differs" =
      identical(found$n_variables_check, as.integer(cells$n_variables)),
    "fits read do not match the gate's stored values: min_ess_bulk differs" =
      isTRUE(all.equal(found$min_ess_bulk_check, cells$min_ess_bulk)),
    "fits read do not match the gate's stored values: min_ess_tail differs" =
      isTRUE(all.equal(found$min_ess_tail_check, cells$min_ess_tail)),
    "fits read do not match the gate's stored values: max_rhat differs" =
      max(abs(found$max_rhat_check - cells$max_rhat)) < 1e-8
  )
  added <- grep("_check$", names(found), invert = TRUE, value = TRUE)
  cells[added] <- found[added]

  res$cells <- cells
  saveRDS(res, asset)
})
