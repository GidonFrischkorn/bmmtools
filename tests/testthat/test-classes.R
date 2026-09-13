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
