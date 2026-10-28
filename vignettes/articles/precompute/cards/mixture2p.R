# Precompute for the mixture2p validation card: the shipped
# validation_defaults("mixture2p") entry, run through every stage
# bmmtools has built, on the standard grid of 20 or 50 subjects and 30,
# 100 or 300 trials, 10 replications each, with mixtur as the agreement
# stage's reference.
#
# Needs bmm (1.4.1 or later), brms, cmdstanr, CmdStan and SBC. Run by hand
# from the package root, when the machine is otherwise idle:
#   Rscript vignettes/articles/precompute/cards/mixture2p.R
# run_validation() resumes by stage, so a rerun reads every stage whose
# inputs are unchanged and fits only what changed. Fits go to
# precompute/fits/cards/mixture2p/, which is not tracked. The result is
# stored as assets/cards/mixture2p.rds, from which the card renders
# without Stan.
#
# The reference is fitted to one of this run's cells, so the card takes
# two passes: this script, then precompute/references/mixture2p-mixtur.R,
# then this script again, which reruns the agreement stage only.

devtools::load_all()
started <- Sys.time()

## ---- plan
library(bmmtools)

# mixtur's estimates, once the reference script has stored them
reference_file <- "vignettes/articles/assets/references/mixture2p-mixtur.rds"
reference <- if (file.exists(reference_file)) readRDS(reference_file)
# every stage derives its seed from this one: recovery 2603, prior
# predictive 2604, SBC 2605
plan <- validation_plan("mixture2p", seed = 2603, reference = reference)
plan

## ---- run
result <- run_validation(
  plan,
  dir = "vignettes/articles/precompute/fits/cards/mixture2p"
)
result

## ---- store
# The run record names the machine it ran on, and the fit paths name the
# folder it ran in; the stored copy carries neither. Paths are kept
# relative to the package root, so a rerun from the root still finds the
# fits.
root <- paste0(normalizePath("."), "/")
relative <- function(path) sub(root, "", path, fixed = TRUE)
result$run_info$hostname <- NA_character_
result$run_record$hostname <- NA_character_
cells <- attr(result$recovery, "cells")
cells$file <- relative(cells$file)
attr(result$recovery, "cells") <- cells
# a recovery grid that records its run (#19) keeps its own copy
grid_info <- attr(result$recovery, "run_info")
if (!is.null(grid_info)) {
  grid_info$hostname <- NA_character_
  attr(result$recovery, "run_info") <- grid_info
}
files <- attr(result$prior_predictive, "files")
files[] <- relative(files)
attr(result$prior_predictive, "files") <- files
# under load_all() the call holds prior_check() itself, with a srcref
# that names the folder; its name is enough
call <- attr(result$prior_predictive, "call")
call[[1L]] <- quote(prior_check)
if (!is.null(call$file)) call$file <- relative(call$file)
attr(result$prior_predictive, "call") <- call
# The prior check keeps its draws, 500 data sets of every trial per prior
# set, 8 MB here. The verdict is computed from all of them already; the
# stored copy keeps the first 50, what plot() draws by default.
draws <- attr(result$prior_predictive, "draws")
draws <- lapply(draws, function(d) {
  d[seq_len(min(50L, nrow(d))), , drop = FALSE]
})
attr(result$prior_predictive, "draws") <- draws

bytes <- serialize(result, NULL)
leaks <- c(root, "/Users/", Sys.info()[["user"]], Sys.info()[["nodename"]])
found <- vapply(leaks, function(s) {
  length(grepRaw(s, bytes, fixed = TRUE, all = TRUE))
}, 1L)
if (any(found > 0L)) {
  stop("The result still names the machine or folder: ",
       paste(names(found)[found > 0L], found[found > 0L], collapse = ", "))
}
dir.create("vignettes/articles/assets/cards", showWarnings = FALSE)
saveRDS(result, "vignettes/articles/assets/cards/mixture2p.rds")

print(as.data.frame(result$run_info))
cat(
  "minutes:", format(as.double(difftime(Sys.time(), started, units = "mins"))),
  "\n"
)
