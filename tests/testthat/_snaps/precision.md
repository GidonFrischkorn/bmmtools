# precision is refused in a shape it cannot use

    Code
      summary(recovery, precision = 1)
    Condition
      Error in `summary()`:
      ! A `precision` vector must name the term of every width.
      i For example `precision = c(kappa = 1)`, which sets the population-level target; a table with "term", "level" and "width" sets other levels.
    Code
      summary(recovery, precision = c(kappa = 1, 2))
    Condition
      Error in `summary()`:
      ! A `precision` vector must name the term of every width.
      i For example `precision = c(kappa = 1)`, which sets the population-level target; a table with "term", "level" and "width" sets other levels.
    Code
      summary(recovery, precision = c(kappa = 0))
    Condition
      Error in `summary()`:
      ! Every `precision` width must be a positive, finite number.
      i It has "0".
    Code
      summary(recovery, precision = c(kappa = -1))
    Condition
      Error in `summary()`:
      ! Every `precision` width must be a positive, finite number.
      i It has "-1".
    Code
      summary(recovery, precision = c(kappa = Inf))
    Condition
      Error in `summary()`:
      ! Every `precision` width must be a positive, finite number.
      i It has "Inf".
    Code
      summary(recovery, precision = c(kappa = NA))
    Condition
      Error in `summary()`:
      ! Every `precision` width must be a positive, finite number.
      i It has "NA".
    Code
      summary(recovery, precision = "1")
    Condition
      Error in `summary()`:
      ! `precision` must be a named numeric vector or a data frame, not a string.
    Code
      summary(recovery, precision = c(thetat = 1))
    Condition
      Error in `summary()`:
      ! `precision` sets a width for "thetat (population)", which the recovery has no rows for.
      i It has "kappa (population)".
    Code
      summary(recovery, precision = data.frame(term = "kappa", width = 1))
    Condition
      Error in `summary()`:
      ! A `precision` table needs the column "level".
      i It has "term" and "width".
    Code
      summary(recovery, precision = data.frame(term = "kappa", level = "subject",
        width = 1))
    Condition
      Error in `summary()`:
      ! `precision` sets a width for "kappa (subject)", which the recovery has no rows for.
      i It has "kappa (population)".
    Code
      summary(recovery, precision = data.frame(term = c("kappa", "kappa"), level = "population",
      width = 1:2))
    Condition
      Error in `summary()`:
      ! `precision` sets a width for "kappa (population)" more than once.
      i Give one row per term and level, or add a "condition" column to tell rows apart.
    Code
      summary(recovery, precision = data.frame(condition = "row-1", term = "kappa",
        level = "population", width = 1))
    Condition
      Error in `summary()`:
      ! `precision` sets a width for "row-1: kappa (population)", which the recovery has no rows for.
      i It has "kappa (population)".
      i It has no conditions; drop the "condition" column.
    Code
      summary(recovery_mixture2p, precision = data.frame(condition = "rwo-1", term = "kappa",
        level = "population", width = 1))
    Condition
      Error in `summary()`:
      ! `precision` sets a width for "rwo-1: kappa (population)", which the recovery has no rows for.
      i It has "kappa (population)", "thetat (population)", "kappa (subject)", and "thetat (subject)".
      i Its conditions are "row-1", "row-2", "row-3", and "row-4".

