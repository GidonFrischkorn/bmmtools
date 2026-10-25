# targets without a Monte Carlo bound are refused, naming one

    Code
      smallest_design(x, c(r = 0.8))
    Condition
      Error in `smallest_design()`:
      ! "r" has no Monte Carlo interval to certify a design with.
      i Use "r_low".
    Code
      smallest_design(x, c(ccc = 0.8))
    Condition
      Error in `smallest_design()`:
      ! "ccc" has no Monte Carlo interval to certify a design with.
      i Use "ccc_low".
    Code
      smallest_design(x, c(ci_width = 1))
    Condition
      Error in `smallest_design()`:
      ! "ci_width" has no Monte Carlo interval to certify a design with.
      i Use "precision_met".
    Code
      smallest_design(x, c(coverage = 0.95))
    Condition
      Error in `smallest_design()`:
      ! "coverage" is not a target: a coverage should be near its nominal mass, not above a threshold.
      i Use `coverage_check = TRUE`.
    Code
      smallest_design(x, c(bias = 0))
    Condition
      Error in `smallest_design()`:
      ! `target` names "bias", which is not a target column.
      i Target columns: "detected", "sign_recovery", "precision_met", "reliability_true", "reliability_model", "r_low", and "ccc_low".
    Code
      smallest_design(x, c(precision_met = 0.8))
    Condition
      Error in `smallest_design()`:
      ! `x` has no columns "precision_met" and "precision_met_mcse".
      i "precision_met" comes from `summary(x, precision = )`.
    Code
      smallest_design(x, 0.8)
    Condition
      Error in `smallest_design()`:
      ! `target` must be a named numeric vector of thresholds.
      i For example `target = c(r_low = 0.8)`.
    Code
      smallest_design(x, c(r_low = NA))
    Condition
      Error in `smallest_design()`:
      ! Every `target` must be a finite number.
    Code
      smallest_design(x, c(r_low = 0.8), along = "n_items")
    Condition
      Error in `smallest_design()`:
      ! `x` has no design column "n_items".
      i A grid summary carries the grid's columns after "condition": summarise a `recovery_grid()` result.
    Code
      smallest_design(x, c(r_low = 0.8), level = 1.2)
    Condition
      Error in `smallest_design()`:
      ! `level` must be a single number between 0 and 1.
    Code
      smallest_design(list(), c(r_low = 0.8))
    Condition
      Error in `smallest_design()`:
      ! `x` must be a grid summary, not an empty list.
    Code
      smallest_design(x, c(r_low = 0.8), along = c("n_trials", "n_subjects"))
    Condition
      Error in `smallest_design()`:
      ! `along` must be one column name.
    Code
      smallest_design(x, c(r_low = 0.8), by = 1)
    Condition
      Error in `smallest_design()`:
      ! `by` must be column names or `NULL`.
    Code
      smallest_design(dplyr::mutate(x, term = NA), c(r_low = 0.8))
    Condition
      Error in `smallest_design()`:
      ! `x` has missing values in "term".
    Code
      smallest_design(dplyr::mutate(x, r_low = NA_real_), c(r_low = 0.8))
    Condition
      Error in `smallest_design()`:
      ! `x` has no rows with a value for the target.

# the print states each status in a sentence

    Code
      smallest_design(x, c(r_low = 0.8, detected = 0.8))
    Output
      <bmmtools_smallest_design>
      n_subjects = 50, c (subject), r_low >= 0.8: not reached in this grid.
      n_subjects = 50, kappa (subject), r_low >= 0.8: met from n_trials = 30, the grid's smallest; the true minimum may be lower.
      n_subjects = 50, thetat (subject), r_low >= 0.8: met from n_trials = 100.
      n_subjects = 50, mu (population), detected >= 0.8: undecided: the point estimate meets it from n_trials = 100, the 95 % lower bound does not; a perfect rate needs at least 16 replications to certify it.

