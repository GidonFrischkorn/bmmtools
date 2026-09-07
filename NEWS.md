# bmmtools 0.0.0.9000

* `recover()` and `recover_subjects()` score fits against known
  generating values and return a `bmmtools_recovery` object: one row per
  fit and parameter, with the estimate, its interval, the generating
  value, the signed error and whether the interval covered. Both take a
  `brmsfit`, a list of them, or an estimates tibble, so a scoring
  pipeline runs with no fitting package installed.
* `summary()` of a recovery object gives the per-parameter metrics: bias,
  RMSE, coverage, mean interval width, the Pearson correlation with a
  Fisher-z interval, the Spearman correlation, and Lin's concordance with
  its scale and location components. Where fewer than three complete
  pairs or no spread makes a correlation undefined it is `NA`, never `0`.
  Subject-level correlations are computed within replication and combined
  on Fisher's z scale.
* `extract_estimates()` turns a fit into the estimates tibble, the
  contract between the run layer and the score layer. Subject-level rows
  are the per-draw sum of the population intercept and the group-level
  deviation; parameters the model fixed are identified by zero posterior
  variance, not by a missing rhat, and dropped by default.
* Scoring on the natural scale handles the two links whose inverse is
  not monotone across zero: a `sqrt`-link interval that spans zero on
  the link scale becomes `[0, max]`, and an `inverse`-link interval that
  spans zero has no natural-scale image, so its bounds and `covered` are
  `NA` and a warning names the term.
* A recovery object survives the `dplyr` verbs that keep its contract
  columns and becomes a plain tibble under one that drops them, so a
  reduced object prints as a tibble instead of failing inside a metric.
* `plot_recovery()` plots estimates against generating values, one panel
  per parameter, with interval bars and an identity line.
* `inverse_link()` transforms a vector from the link scale to the natural
  scale for the twelve links bmm uses.
* Package skeleton and design record (`ARCHITECTURE.md`), Milestone 0.
