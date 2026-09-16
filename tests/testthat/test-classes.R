# Tests for the recovery classes and their methods, written against
# local/dev/spec-milestone-1-score-layer.md section 6.

# Truth varies across replications, so the correlation metrics are
# estimable once there are three or more of them --- which is what makes
# the degenerate-case note appear for one fit and not for eight.
recovery_example <- function(n_replications = 5L) {
  terms <- c("kappa", "thetat")
  estimates <- list()
  truth <- list()
  for (rep in seq_len(n_replications)) {
    true_value <- c(1, 2) + rep / 10
    estimates[[rep]] <- fake_estimates(
      terms,
      estimate = true_value + c(0.05, -0.07) * rep, replication = rep
    )
    truth[[rep]] <- fake_truth(terms,
      true_value = true_value,
      replication = rep
    )
  }
  recover(dplyr::bind_rows(estimates), dplyr::bind_rows(truth),
    scale = "link"
  )
}

# the constructor -------------------------------------------------------

test_that("the constructor names the first missing contract column", {
  x <- recovery_example()
  incomplete <- tibble::as_tibble(x)[setdiff(names(x), "covered")]

  err <- expect_error(
    new_bmmtools_recovery(incomplete, scale = "link", ci_level = 0.95)
  )
  expect_match(conditionMessage(err), "covered")
})

test_that("the constructor rejects a column of the wrong type", {
  x <- tibble::as_tibble(recovery_example())
  x$covered <- as.character(x$covered)

  err <- expect_error(
    new_bmmtools_recovery(x, scale = "link", ci_level = 0.95)
  )
  expect_match(conditionMessage(err), "covered")
})

test_that("the constructor stores scale, ci_level and the call", {
  x <- tibble::as_tibble(recovery_example())
  out <- new_bmmtools_recovery(
    x,
    scale = "natural", ci_level = 0.89, call = quote(f(1))
  )
  expect_equal(attr(out, "scale"), "natural")
  expect_equal(attr(out, "ci_level"), 0.89)
  expect_equal(attr(out, "call"), quote(f(1)))
})

# printing --------------------------------------------------------------

test_that("print returns its input invisibly", {
  x <- recovery_example()
  expect_invisible(withr::with_output_sink(nullfile(), print(x)))
  returned <- withr::with_output_sink(nullfile(), print(x))
  expect_identical(returned, x)
})

test_that("format names the scale, the fits and the parameters", {
  x <- recovery_example(n_replications = 4L)
  out <- paste(format(x), collapse = "\n")

  expect_type(format(x), "character")
  expect_match(out, "link")
  expect_match(out, "4")
  expect_match(out, "kappa")
})

test_that("a single fit prints the explanation rather than NA cells", {
  # Ported from recovery_report(): a table of NA in the correlation
  # columns is not an answer, it is an unexplained blank.
  x <- recovery_example(n_replications = 1L)
  out <- paste(format(x), collapse = "\n")
  expect_match(out, "at least 3")
})

test_that("several fits print without the degenerate-case explanation", {
  x <- recovery_example(n_replications = 8L)
  out <- paste(format(x), collapse = "\n")
  expect_no_match(out, "at least 3")
})

test_that("the summary prints and returns invisibly", {
  s <- summary(recovery_example())
  expect_invisible(withr::with_output_sink(nullfile(), print(s)))
  out <- utils::capture.output(print(s))
  expect_true(any(grepl("kappa", out)))
})

# dplyr ------------------------------------------------------------------

test_that("dplyr verbs keep the object usable", {
  x <- recovery_example()
  filtered <- dplyr::filter(x, term == "kappa")

  expect_s3_class(filtered, "bmmtools_recovery")
  expect_equal(nrow(filtered), 5L)
  expect_equal(attr(filtered, "scale"), "link")
  expect_s3_class(summary(filtered), "bmmtools_recovery_summary")
})

test_that("dropping a contract column drops the class, not the data", {
  # What is left is still a tibble a user can work with. The negative
  # assertion is the point: `tbl_df` stays true whether or not the
  # recovery class was dropped, so asserting only that would pass while
  # the object was still mislabelled.
  x <- recovery_example()
  reduced <- dplyr::select(x, "term", "estimate")

  expect_false(inherits(reduced, "bmmtools_recovery"))
  expect_s3_class(reduced, "tbl_df")
  expect_equal(ncol(reduced), 2L)
  # and it prints as a tibble instead of dying inside a metric
  expect_output(print(reduced), "term")
})

test_that("subsetting rows or all contract columns keeps the class", {
  x <- recovery_example()

  expect_s3_class(x[1:3, ], "bmmtools_recovery")
  expect_s3_class(x[recovery_contract_columns()], "bmmtools_recovery")
  expect_s3_class(dplyr::mutate(x, extra = 1), "bmmtools_recovery")
  # dropping one contract column is enough to demote
  expect_false(inherits(x[setdiff(names(x), "covered")], "bmmtools_recovery"))
})

test_that("a hand-broken object is named, not crashed into", {
  # Reachable only by attaching the class to something that does not
  # satisfy the contract; the methods say which column is missing rather
  # than failing three frames deeper inside a metric.
  x <- recovery_example()
  broken <- tibble::as_tibble(x)[c("term", "estimate")]
  class(broken) <- c("bmmtools_recovery", class(broken))

  expect_error(summary(broken), "true_value")
  expect_error(format(broken), "true_value")
})

# summary ----------------------------------------------------------------

test_that("summary is a tibble subclass and is idempotent", {
  s <- summary(recovery_example())
  expect_s3_class(s, "bmmtools_recovery_summary")
  expect_s3_class(s, "tbl_df")
  expect_identical(summary(s), s)
})

test_that("summary of a mixed-level object keeps the levels apart", {
  population <- recovery_example()
  ids <- as.character(1:5)
  subject_estimates <- fake_estimates(
    rep("kappa", 5),
    estimate = c(1, 2, 3, 4, 5),
    level = "subject", id = ids
  )
  subject_truth <- tibble::tibble(
    id = ids, term = "kappa", true_value = c(1, 2, 3, 4, 5.5)
  )
  subjects <- recover_subjects(subject_estimates, subject_truth,
    scale = "link"
  )

  combined <- new_bmmtools_recovery(
    dplyr::bind_rows(
      tibble::as_tibble(population),
      tibble::as_tibble(subjects)
    ),
    scale = "link", ci_level = 0.95
  )
  s <- summary(combined)

  expect_setequal(s$level, c("population", "subject"))
  expect_equal(nrow(s), 3L)
  expect_equal(s$n_replications[s$level == "subject"], 1L)
})

# combining correlations -------------------------------------------------

test_that("fisher_z_combine reduces to the weighted average it documents", {
  # equal subject counts: the plain average of the transformed values
  out <- fisher_z_combine(r = c(0.5, 0.8), n = c(20, 20))
  expect_equal(out$r, tanh(mean(atanh(c(0.5, 0.8)))))

  # unequal counts: the larger replication pulls harder, by n - 3
  uneven <- fisher_z_combine(r = c(0.5, 0.8), n = c(5, 100))
  expect_gt(uneven$r, out$r)
  expect_equal(
    uneven$r,
    tanh(sum(c(2, 97) * atanh(c(0.5, 0.8))) / 99)
  )
})

test_that("opposite perfect correlations give NA, not NaN", {
  # Two replications recovered perfectly in opposite directions put z at
  # +Inf and -Inf; their average is undefined. NaN would print as a
  # different missing value from every other guard in the package, so
  # both the weighted and the unweighted branch return NA_real_.
  weighted <- fisher_z_combine(r = c(1, -1), n = c(20, 20))
  expect_identical(weighted$r, NA_real_)
  expect_identical(weighted$r_low, NA_real_)

  # n <= 3 everywhere floors every weight to zero, taking the other branch
  unweighted <- fisher_z_combine(r = c(1, -1), n = c(3, 3))
  expect_identical(unweighted$r, NA_real_)

  # a perfect correlation on its own is still 1, not NA
  expect_equal(fisher_z_combine(r = c(1, 1), n = c(20, 20))$r, 1)
})

# the cross-check class --------------------------------------------------

# One fit, several subjects: the shape summary() needs to have three or
# more complete pairs per term, so the NA guards can be tested on both
# sides of the boundary.
cross_check_example <- function(n_subjects = 6L,
                                terms = c("kappa", "thetat"),
                                intervals = TRUE) {
  ids <- as.character(seq_len(n_subjects))
  fit <- structure(
    list(parameters = terms, ids = ids), class = "mockfit"
  )
  reference <- tibble::tibble(
    term = rep(terms, each = n_subjects),
    estimate = rep(seq_len(n_subjects) / 10, times = length(terms)),
    id = rep(ids, times = length(terms))
  )
  if (intervals) {
    reference$ci_low <- reference$estimate - 0.05
    reference$ci_high <- reference$estimate + 0.05
  }
  cross_check(fit, reference, scale = "link", level = "subject")
}

test_that("the cross-check constructor names the first missing column", {
  x <- cross_check_example()
  incomplete <- tibble::as_tibble(x)[setdiff(names(x), "overlap")]

  err <- expect_error(
    new_bmmtools_cross_check(incomplete, scale = "link", ci_level = 0.95)
  )
  expect_match(conditionMessage(err), "overlap")
})

test_that("the cross-check constructor rejects a column of the wrong type", {
  x <- tibble::as_tibble(cross_check_example())
  x$covered <- as.character(x$covered)

  err <- expect_error(
    new_bmmtools_cross_check(x, scale = "link", ci_level = 0.95)
  )
  expect_match(conditionMessage(err), "covered")
})

test_that("a cross-check keeps its class through filter, not select", {
  x <- cross_check_example()

  kept <- dplyr::filter(x, .data$term == "kappa")
  expect_s3_class(kept, "bmmtools_cross_check")

  dropped <- dplyr::select(x, "term", "estimate")
  expect_false(inherits(dropped, "bmmtools_cross_check"))
  expect_s3_class(dropped, "tbl_df")

  # `reference` is a contract column, so dropping it demotes
  without_reference <- x[setdiff(names(x), "reference")]
  expect_false(inherits(without_reference, "bmmtools_cross_check"))
})

test_that("row subsetting a cross-check keeps the class", {
  x <- cross_check_example()
  expect_s3_class(x[1:3, ], "bmmtools_cross_check")
})

test_that("summary of a cross-check returns the documented columns", {
  s <- summary(cross_check_example())
  expect_equal(names(s), cross_check_summary_columns())
  expect_s3_class(s, "bmmtools_cross_check_summary")
  expect_equal(nrow(s), 2L)
  expect_setequal(s$term, c("kappa", "thetat"))
  expect_equal(unique(s$level), "subject")
  expect_equal(unique(s$scale), "link")
  expect_equal(s$n, c(6, 6))
  # the mock reports a converged fit on every row
  expect_equal(s$n_converged, c(6L, 6L))
})

test_that("summary of a cross-check is NA below three pairs", {
  s <- summary(cross_check_example(n_subjects = 2L))
  expect_true(all(is.na(s$r)))
  expect_true(all(is.na(s$ccc)))
  # the errors are still defined with two pairs
  expect_false(anyNA(s$bias))
  expect_false(anyNA(s$rmse))
})

test_that("share_overlap is NA without reference intervals", {
  with_intervals <- summary(cross_check_example())
  expect_false(anyNA(with_intervals$share_overlap))

  without <- summary(cross_check_example(intervals = FALSE))
  expect_true(all(is.na(without$share_overlap)))
})

test_that("summary of a cross-check summary is itself", {
  s <- summary(cross_check_example())
  expect_identical(summary(s), s)
})

test_that("an empty cross-check summarises to the empty contract", {
  x <- cross_check_example()
  empty <- x[integer(0), ]
  s <- summary(empty)
  expect_equal(names(s), cross_check_summary_columns())
  expect_equal(nrow(s), 0L)
})

test_that("printing a cross-check names the scale and the reference", {
  x <- cross_check_example()
  out <- format(x)
  expect_match(out[[1L]], "bmmtools_cross_check")
  expect_true(any(grepl("link scale", out)))
  expect_true(any(grepl("reference", out)))
  expect_output(print(x), "bmmtools_cross_check")
})

test_that("printing an empty cross-check says so", {
  x <- cross_check_example()
  expect_true(any(grepl("No parameters", format(x[integer(0), ]))))
})

test_that("a cross-check whose contract was broken by hand is refused", {
  x <- cross_check_example()
  broken <- structure(
    tibble::as_tibble(x)[setdiff(names(x), "covered")],
    class = class(x)
  )
  expect_error(summary(broken), "covered")
  expect_error(format(broken), "covered")
})

test_that("the cross-check constructor refuses a non-data-frame", {
  expect_error(
    new_bmmtools_cross_check(1:3, scale = "link", ci_level = 0.95),
    "data frame"
  )
})

test_that("a cross-check without a verdict counts n_converged as NA", {
  x <- cross_check_example()
  bare <- tibble::as_tibble(x)[setdiff(names(x), "converged")]
  rebuilt <- new_bmmtools_cross_check(bare, scale = "link", ci_level = 0.95)

  expect_true(all(is.na(rebuilt$converged)))
  expect_identical(unique(summary(rebuilt)$n_converged), NA_integer_)
})

test_that("a verb that rebuilds a cross-check without the contract demotes", {
  x <- cross_check_example()
  reduced <- dplyr::summarise(x, m = mean(.data$bias))
  expect_false(inherits(reduced, "bmmtools_cross_check"))
  expect_s3_class(reduced, "tbl_df")
})

test_that("dplyr_reconstruct on a cross-check demotes a broken contract", {
  # `[` intercepts every column-dropping verb, so this branch is the
  # belt to that method's braces and is reached directly
  x <- cross_check_example()
  out <- dplyr::dplyr_reconstruct(tibble::tibble(term = "kappa"), x)
  expect_false(inherits(out, "bmmtools_cross_check"))
  expect_s3_class(out, "tbl_df")
})
