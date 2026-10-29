# Reference for the mixture2p validation card: mixtur's maximum-likelihood
# fits of the two-component mixture model, on the data of the recovery
# cell the card's agreement stage compares (the largest design's first
# fitted replication).
#
# mixtur fits each subject on its own and returns point estimates. The
# population values are the per-subject estimates averaged on the link
# scale bmm estimates on, mean(log kappa) and mean(qlogis(p_t)), which is
# what bmm's population intercepts are, and reported on the natural scale
# the card compares on; the interval is a 95 % percentile interval over
# subjects resampled with replacement.
#
# mixtur is not a dependency of bmmtools, not even in Suggests: install it
# by hand. Run from the package root after the first run of
# precompute/cards/mixture2p.R, which fits the cell:
#   Rscript vignettes/articles/precompute/references/mixture2p-mixtur.R
# then rerun the card script, which reads the stored reference and reruns
# the agreement stage only. The result is stored in assets/references/.

devtools::load_all()
library(bmmtools)

## ---- cell
# the cell stage_agreement() compares
plan <- validation_plan("mixture2p", seed = 2603)
dir <- "vignettes/articles/precompute/fits/cards/mixture2p"
stored <- readRDS(file.path(dir, "result.rds"))
cells <- attr(stored$recovery, "cells")
condition <- bmmtools:::largest_condition(plan$grid)
cells <- cells[cells$condition == condition & cells$status == "ok", ]
cell <- cells[order(cells$replication), ][1L, ]
fit <- readRDS(cell$file)
data <- fit$data
response <- plan$model$resp_vars$resp_error
cat("cell:", condition, "replication", cell$replication, "-",
    nlevels(factor(data$id)), "subjects,", nrow(data), "trials\n")

## ---- fit
input <- data.frame(id = data$id, response = data[[response]], target = 0)
subjects <- mixtur::fit_mixtur(
  input,
  model = "2_component", unit = "radians", id_var = "id",
  response_var = "response", target_var = "target"
)
degenerate <- !is.finite(log(subjects$kappa)) |
  !is.finite(stats::qlogis(subjects$p_t))
if (any(degenerate)) {
  stop(
    sum(degenerate), " subject fit(s) at a boundary (kappa 0 or p_t 0 or 1),",
    " so the link-scale mean is not defined: ",
    paste(subjects$id[degenerate], collapse = ", ")
  )
}

## ---- reduce
link_means <- function(s) {
  c(kappa = mean(log(s$kappa)), thetat = mean(stats::qlogis(s$p_t)))
}
estimate <- link_means(subjects)
# the seed is the plan's plus 3, after the three the stages use
boot <- withr::with_seed(2606, replicate(2000L, {
  link_means(subjects[sample.int(nrow(subjects), replace = TRUE), ])
}))
bounds <- apply(boot, 1L, stats::quantile, probs = c(0.025, 0.975))
# the natural scale the card compares on, as cross_check() transforms the
# fit: kappa through exp(), thetat through plogis() to the probability of
# a response from memory
natural <- function(x) {
  c(kappa = exp(x[["kappa"]]), thetat = stats::plogis(x[["thetat"]]))
}

reference <- tibble::tibble(
  term = names(estimate),
  estimate = natural(estimate),
  ci_low = natural(bounds[1L, ]),
  ci_high = natural(bounds[2L, ]),
  source = "mixtur",
  package = "mixtur",
  version = as.character(utils::packageVersion("mixtur")),
  call = paste(
    "fit_mixtur(model = '2_component', unit = 'radians') per subject;",
    "link-scale mean; 2000 subject resamples, seed 2606"
  ),
  data_hash = rlang::hash(data),
  date = as.character(Sys.Date())
)
reference

## ---- store
dir.create("vignettes/articles/assets/references", showWarnings = FALSE)
saveRDS(reference, "vignettes/articles/assets/references/mixture2p-mixtur.rds")
