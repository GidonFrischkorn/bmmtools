# Adds to assets/recovery-ezdm.rds, for every fit, the variable at which
# the convergence gate's numbers sit: the lowest bulk ESS, the lowest tail
# ESS and the largest R-hat, and the parameter with the lowest bulk ESS
# (lp__ and lprior left out). It reads the fits that recovery-ezdm.R saved
# and fits nothing. recovery-ezdm.R sources it at its end; to add the
# columns to a stored asset, run it alone from the package root, with the
# fits folder as an argument if the fits were moved:
#   Rscript vignettes/articles/precompute/recovery-ezdm-diagnostics.R [fits]
#
# Needs brms (to read the draws of a fit), not Stan.

asset <- "vignettes/articles/assets/recovery-ezdm.rds"
fits <- commandArgs(trailingOnly = TRUE)
fits <- if (length(fits) == 1L) {
  fits
} else {
  "vignettes/articles/precompute/fits/recovery-ezdm"
}
res <- readRDS(asset)
cells <- res$cells

# The gate summarises every variable of a fit and drops those without an
# R-hat (constants); the same set is read here, so the minimum found must
# be the one the gate stored.
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
    min_ess_bulk_check = min(d$ess_bulk),
    min_ess_bulk_variable = d$variable[which.min(d$ess_bulk)],
    min_ess_bulk_parameter = p$variable[which.min(p$ess_bulk)],
    min_ess_tail_variable = d$variable[which.min(d$ess_tail)],
    max_rhat_variable = d$variable[which.max(d$rhat)]
  )
}
found <- do.call(rbind, lapply(file.path(fits, basename(cells$file)), weakest))
stopifnot(isTRUE(all.equal(found$min_ess_bulk_check, cells$min_ess_bulk)))
cells[names(found)[-1L]] <- found[-1L]

res$cells <- cells
saveRDS(res, asset)
