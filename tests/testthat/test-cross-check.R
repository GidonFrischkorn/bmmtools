# Tests for cross_check(), written against
# local/dev/spec-milestone-6-sbc-cross-check.md section 2.
#
# Three routes, and the split is the point. The validation table and the
# arithmetic need nothing at all. The wiring runs on the mockfit whose
# extract_estimates() method helper-generate.R registers. The two claims
# that only a real fit can support --- that a natural-scale cross_check()
# agrees with recover() to the last digit, and that the subject level
# joins on brms's own ids --- run on the committed fixture. No test
# compiles Stan.

#' A mock fit with the terms and ids the registered method reads
#' @noRd
mock_cross_fit <- function(terms = c("kappa", "thetat"),
                           ids = as.character(1:6)) {
  structure(list(parameters = terms, ids = ids), class = "mockfit")
}

#' A reference the mock's estimates (0, [-10, 10]) can be scored against
#'
#' `kappa` is covered and its interval overlaps; `thetat` sits outside the
#' fit's interval and its own interval does not reach it.
#'
#' @noRd
mock_reference <- function(intervals = TRUE) {
  out <- tibble::tibble(
    term = c("kappa", "thetat"),
    estimate = c(0.5, 11),
    ci_low = c(-1, 10.5),
    ci_high = c(1, 11.5)
  )
  if (intervals) out else out[c("term", "estimate")]
}

# the fit argument -------------------------------------------------------

test_that("a list of fits is an error naming recover()", {
  fits <- list(mock_cross_fit(), mock_cross_fit())
  err <- expect_error(cross_check(fits, mock_reference(), scale = "link"))
  expect_match(conditionMessage(err), "recover")
})

test_that("something that is not a fit is an error", {
  expect_error(
    cross_check(1:3, mock_reference(), scale = "link"),
    "extract_estimates"
  )
})

# the level argument -----------------------------------------------------

test_that("level 'sd' and 'cor' are an error naming the two levels", {
  for (level in c("sd", "cor")) {
    err <- expect_error(
      cross_check(mock_cross_fit(), mock_reference(),
        scale = "link", level = level
      )
    )
    expect_match(conditionMessage(err), "population")
    expect_match(conditionMessage(err), "subject")
  }
})

test_that("an unknown level is an error", {
  expect_error(
    cross_check(mock_cross_fit(), mock_reference(),
      scale = "link", level = "nonsense"
    )
  )
})

test_that("several levels at once are an error", {
  expect_error(
    cross_check(mock_cross_fit(), mock_reference(),
      scale = "link", level = c("population", "sd")
    ),
    "single level"
  )
})

# check_reference() ------------------------------------------------------

test_that("a reference that is neither a data frame nor a fit is an error", {
  expect_error(
    cross_check(mock_cross_fit(), "kappa", scale = "link"),
    "data frame"
  )
})

test_that("a reference without term or estimate is an error", {
  expect_error(
    cross_check(mock_cross_fit(), tibble::tibble(estimate = 1),
      scale = "link"
    ),
    "term"
  )
  expect_error(
    cross_check(mock_cross_fit(), tibble::tibble(term = "kappa"),
      scale = "link"
    ),
    "estimate"
  )
})

test_that("a non-character term or non-numeric estimate is an error", {
  expect_error(
    cross_check(
      mock_cross_fit(),
      tibble::tibble(term = 1, estimate = 1),
      scale = "link"
    ),
    "character"
  )
  expect_error(
    cross_check(
      mock_cross_fit(),
      tibble::tibble(term = "kappa", estimate = "a"),
      scale = "link"
    ),
    "numeric"
  )
})

test_that("one interval bound without the other is an error", {
  reference <- tibble::tibble(term = "kappa", estimate = 1, ci_low = 0)
  err <- expect_error(
    cross_check(mock_cross_fit(), reference, scale = "link")
  )
  expect_match(conditionMessage(err), "ci_high")
})

test_that("source defaults to 'reference' and is kept when given", {
  x <- cross_check(mock_cross_fit(), mock_reference(), scale = "link")
  expect_equal(unique(x$source), "reference")

  named <- mock_reference()
  named$source <- c("sdt_d", "published")
  y <- cross_check(mock_cross_fit(), named, scale = "link")
  expect_equal(y$source, c("sdt_d", "published"))
})

# the arithmetic ---------------------------------------------------------

test_that("bias, covered and overlap are the documented arithmetic", {
  x <- cross_check(mock_cross_fit(), mock_reference(), scale = "link")

  expect_equal(x$term, c("kappa", "thetat"))
  expect_equal(x$estimate, c(0, 0))
  expect_equal(x$reference, c(0.5, 11))
  expect_equal(x$bias, c(-0.5, -11))
  expect_equal(x$covered, c(TRUE, FALSE))
  expect_equal(x$overlap, c(TRUE, FALSE))
})

test_that("a reference without intervals gives NA bounds and NA overlap", {
  x <- cross_check(mock_cross_fit(), mock_reference(intervals = FALSE),
    scale = "link"
  )
  expect_true(all(is.na(x$ref_low)))
  expect_true(all(is.na(x$ref_high)))
  expect_true(all(is.na(x$overlap)))
  # the comparison itself is unaffected
  expect_equal(x$bias, c(-0.5, -11))
  expect_equal(x$covered, c(TRUE, FALSE))
})

test_that("overlap is NA per row, not FALSE, where one bound is missing", {
  reference <- mock_reference()
  reference$ci_low[[1L]] <- NA_real_
  reference$ci_high[[1L]] <- NA_real_
  x <- cross_check(mock_cross_fit(), reference, scale = "link")
  expect_equal(x$overlap, c(NA, FALSE))
})

test_that("cross_check_rows() computes the same on hand-built tibbles", {
  estimates <- fake_estimates(
    c("kappa", "thetat"),
    estimate = c(1, 2), ci_low = c(0, 1), ci_high = c(2, 3)
  )
  estimates$converged <- TRUE
  reference <- tibble::tibble(
    term = c("kappa", "thetat"),
    reference = c(1.5, 5),
    ref_low = c(1.4, NA_real_),
    ref_high = c(1.6, NA_real_),
    source = "reference",
    id = NA_character_
  )
  out <- cross_check_rows(estimates, reference, "population", "link")

  expect_equal(out$bias, c(-0.5, -3))
  expect_equal(out$covered, c(TRUE, FALSE))
  expect_equal(out$overlap, c(TRUE, NA))
  expect_equal(unique(out$scale), "link")
})

# the join ---------------------------------------------------------------

test_that("a reference term the fit lacks is a warning and is dropped", {
  reference <- mock_reference()
  reference$term <- c("kappa", "not_a_parameter")
  expect_warning(
    x <- cross_check(mock_cross_fit(), reference, scale = "link"),
    "not_a_parameter"
  )
  expect_equal(x$term, "kappa")
})

test_that("no reference term in the fit at all is an error", {
  reference <- mock_reference()
  reference$term <- c("nope", "also_nope")
  err <- expect_error(
    cross_check(mock_cross_fit(), reference, scale = "link")
  )
  expect_match(conditionMessage(err), "reference")
  expect_match(conditionMessage(err), "kappa")
})

test_that("a fit term the reference lacks is dropped in silence", {
  reference <- mock_reference()[1L, ]
  x <- expect_silent(
    cross_check(mock_cross_fit(), reference, scale = "link")
  )
  expect_equal(x$term, "kappa")
})

# the subject level ------------------------------------------------------

test_that("a subject-level reference without id is an error", {
  err <- expect_error(
    cross_check(mock_cross_fit(), mock_reference(),
      scale = "link", level = "subject"
    )
  )
  expect_match(conditionMessage(err), "id")
})

test_that("a subject-level reference id must be character", {
  reference <- tibble::tibble(
    term = "kappa", estimate = 1, id = 1L
  )
  err <- expect_error(
    cross_check(mock_cross_fit(), reference, scale = "link", level = "subject")
  )
  expect_match(conditionMessage(err), "id")
})

test_that("reference ids the fit lacks are a warning and are dropped", {
  reference <- tibble::tibble(
    term = "kappa",
    estimate = 0.5,
    id = c("1", "2", "99")
  )
  expect_warning(
    x <- cross_check(mock_cross_fit(), reference,
      scale = "link", level = "subject"
    ),
    "99"
  )
  expect_equal(x$id, c("1", "2"))
})

test_that("no id in common is an error", {
  reference <- tibble::tibble(term = "kappa", estimate = 0.5, id = "99")
  expect_error(
    cross_check(mock_cross_fit(), reference,
      scale = "link", level = "subject"
    )
  )
})

test_that("a fit reference with no rows at the level is an error", {
  err <- expect_error(
    cross_check(
      mock_cross_fit(),
      mock_cross_fit(ids = character(0)),
      scale = "link", level = "subject"
    )
  )
  expect_match(conditionMessage(err), "reference")
})

# the committed fixture --------------------------------------------------

test_that("a natural-scale cross_check equals recover() on the same fit", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  truth <- mixture2p_truth()

  scored <- recover(fit, truth$population)
  reference <- tibble::tibble(
    term = truth$population$term,
    estimate = truth$population$true_value
  )
  crossed <- cross_check(fit, reference)

  expect_equal(crossed$term, scored$term)
  expect_equal(crossed$estimate, scored$estimate, tolerance = 1e-12)
  expect_equal(crossed$ci_low, scored$ci_low, tolerance = 1e-12)
  expect_equal(crossed$ci_high, scored$ci_high, tolerance = 1e-12)
  # the reference is on the comparison scale already and is not
  # transformed a second time
  expect_equal(crossed$reference, truth$population$true_value)
  expect_equal(crossed$bias, crossed$estimate - crossed$reference)
})

test_that("a fit against itself at subject level has bias 0 and r 1", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  x <- cross_check(fit, fit, level = "subject")

  expect_s3_class(x, "bmmtools_cross_check")
  expect_equal(unique(x$source), "fit")
  expect_equal(x$bias, rep(0, nrow(x)))
  expect_true(all(x$covered))
  expect_true(all(x$overlap))

  summarised <- summary(x)
  expect_equal(summarised$r, rep(1, nrow(summarised)))
  expect_equal(summarised$ccc, rep(1, nrow(summarised)))
  expect_equal(summarised$bias, rep(0, nrow(summarised)))
})

test_that("a fit reference is extracted on the same scale as the fit", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  natural <- cross_check(fit, fit)
  link <- cross_check(fit, fit, scale = "link")

  expect_equal(unique(natural$scale), "natural")
  expect_equal(unique(link$scale), "link")
  # kappa has a log link, so the two scales differ
  expect_false(isTRUE(all.equal(
    natural$estimate[natural$term == "kappa"],
    link$estimate[link$term == "kappa"]
  )))
  # both sides move together, so the comparison is 0 either way
  expect_equal(natural$bias, rep(0, nrow(natural)))
  expect_equal(link$bias, rep(0, nrow(link)))
})

test_that("the subject level joins on the fit's own ids", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  truth <- mixture2p_truth()

  reference <- tibble::tibble(
    term = truth$subjects$term,
    estimate = truth$subjects$true_value,
    id = truth$subjects$id
  )
  x <- cross_check(fit, reference, level = "subject", scale = "link")

  expect_equal(nrow(x), nrow(truth$subjects))
  expect_setequal(x$id, unique(truth$subjects$id))
  expect_equal(unique(x$level), "subject")
})

# the object -------------------------------------------------------------

test_that("the returned object carries the contract and its attributes", {
  x <- cross_check(mock_cross_fit(), mock_reference(),
    scale = "link", ci_level = 0.89
  )
  expect_s3_class(x, "bmmtools_cross_check")
  expect_named(x, names(cross_check_contract()))
  expect_equal(attr(x, "scale"), "link")
  expect_equal(attr(x, "ci_level"), 0.89)
  expect_false(is.null(attr(x, "call")))
})

test_that("ci_level is validated", {
  expect_error(
    cross_check(mock_cross_fit(), mock_reference(),
      scale = "link", ci_level = 2
    ),
    "between 0 and 1"
  )
})

# links ------------------------------------------------------------------

test_that("a fit with no link table falls back to the link scale and says so", {
  expect_message(
    x <- cross_check(mock_cross_fit(), mock_reference()),
    "link scale"
  )
  expect_equal(unique(x$scale), "link")
})

test_that("links given by hand transform the fit side only", {
  reference <- tibble::tibble(term = "kappa", estimate = 1)
  x <- cross_check(mock_cross_fit(), reference, links = c(kappa = "log"))
  # the mock estimates 0 on the link scale, so exp(0) = 1 naturally
  expect_equal(x$estimate, 1)
  expect_equal(x$reference, 1)
  expect_equal(x$bias, 0)
  expect_equal(unique(x$scale), "natural")
})

test_that("a fit reference is inverted with its own link table", {
  skip_if_not_installed("brms")
  # Another implementation may hold the same parameter on a different
  # link. Inverting the reference's draws with the fit's table would put
  # a silently mis-scaled number in the `reference` column --- exactly
  # the class of error this milestone exists to close --- so each side
  # goes through its own.
  fit <- mixture2p_fit()
  other <- fit
  other$bmm$model$links$kappa <- "identity"

  x <- cross_check(fit, other, level = "population")
  on_link <- cross_check(fit, other, level = "population", scale = "link")

  kappa <- x$term == "kappa"
  expect_equal(
    x$reference[kappa], on_link$reference[on_link$term == "kappa"]
  )
  expect_equal(
    x$estimate[kappa], exp(on_link$estimate[on_link$term == "kappa"])
  )
  expect_equal(x$bias[kappa], x$estimate[kappa] - x$reference[kappa])

  # thetat has the same logit link on both sides, so it is unchanged
  thetat <- x$term == "thetat"
  expect_equal(x$bias[thetat], 0)
})

test_that("a fit reference with no link table falls back to the fit's", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  stripped <- fit
  stripped$bmm <- NULL

  x <- cross_check(fit, stripped, level = "population")
  expect_equal(unique(x$scale), "natural")
  expect_equal(x$bias, rep(0, nrow(x)))
})

test_that("links given by hand override both sides", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  other <- fit
  other$bmm$model$links$kappa <- "identity"

  x <- cross_check(
    fit, other,
    level = "population", links = c(kappa = "log", thetat = "logit")
  )
  expect_equal(x$bias, rep(0, nrow(x)))
})

test_that("a non-numeric interval bound is an error, not a silent NA", {
  reference <- tibble::tibble(
    term = "kappa", estimate = 1, ci_low = "a", ci_high = "b"
  )
  err <- expect_error(
    cross_check(mock_cross_fit(), reference, scale = "link")
  )
  expect_match(conditionMessage(err), "ci_low")
  expect_match(conditionMessage(err), "numeric")
})
