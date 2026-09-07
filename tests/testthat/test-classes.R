# Tests for the recovery classes and their methods, written against
# dev/spec-milestone-1-score-layer.md section 6.

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
  # Standard tibble-subclass behaviour. Documented rather than fought:
  # what is left is still a tibble a user can work with.
  x <- recovery_example()
  reduced <- dplyr::select(x, "term", "estimate")
  expect_s3_class(reduced, "tbl_df")
  expect_equal(ncol(reduced), 2L)
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
