# Tests for recover() and recover_subjects(), written against
# local/dev/spec-milestone-1-score-layer.md sections 4 and 5.
#
# Almost everything here runs on hand-built estimates and truth tibbles.
# That is the point of decision 1: the score layer is a contract between
# two tibbles, so it can be tested --- and installed --- with no fitting
# package present. The two tests that use the saved fixture say so and
# skip when brms is absent.

# the contract ----------------------------------------------------------

test_that("recover returns one row per replication and term", {
  terms <- c("kappa", "thetat", "a", "c")
  estimates <- dplyr::bind_rows(lapply(1:3, function(rep) {
    fake_estimates(terms,
      estimate = rep + seq_along(terms) / 10,
      replication = rep
    )
  }))
  truth <- fake_truth(terms, true_value = seq_along(terms) / 10)

  out <- recover(estimates, truth, scale = "link")

  expect_equal(nrow(out), 12L)
  expect_named(out, recovery_contract_columns())
  expect_s3_class(out, "bmmtools_recovery")
  expect_s3_class(out, "tbl_df")
  expect_setequal(out$term, terms)
  expect_setequal(out$replication, 1:3)
})

test_that("bias is estimate minus truth, row by row", {
  estimates <- fake_estimates(c("a", "b"), estimate = c(2, -1))
  truth <- fake_truth(c("a", "b"), true_value = c(1.5, -3))

  out <- recover(estimates, truth, scale = "link")

  expect_equal(out$bias, out$estimate - out$true_value)
  expect_equal(out$bias, c(0.5, 2))
})

test_that("covered is TRUE exactly when truth lies inside the interval", {
  estimates <- fake_estimates(
    c("a", "b", "c", "d"),
    estimate = c(0, 0, 0, 0),
    ci_low = c(-1, -1, -1, 2),
    ci_high = c(1, 1, 1, 3)
  )
  truth <- fake_truth(c("a", "b", "c", "d"), true_value = c(0, -1, 1.5, 5))

  out <- recover(estimates, truth, scale = "link")

  # inside, on the boundary (closed interval), outside, far outside
  expect_equal(out$covered, c(TRUE, TRUE, FALSE, FALSE))
  expect_type(out$covered, "logical")
})

test_that("the scale, ci_level and call are carried as attributes", {
  estimates <- fake_estimates("a", estimate = 1)
  truth <- fake_truth("a", true_value = 1)
  out <- recover(estimates, truth, scale = "link")

  expect_equal(attr(out, "scale"), "link")
  expect_equal(attr(out, "ci_level"), 0.95)
  expect_false(is.null(attr(out, "call")))
})

# the link transform ----------------------------------------------------

test_that("scale = natural transforms estimate, both bounds and truth", {
  estimates <- fake_estimates("kappa", estimate = 0, ci_low = -1, ci_high = 1)
  truth <- fake_truth("kappa", true_value = 0)

  out <- recover(estimates, truth, scale = "natural", links = c(kappa = "log"))

  expect_equal(out$estimate, 1)
  expect_equal(out$true_value, 1)
  expect_equal(out$bias, 0)
  expect_equal(out$ci_low, exp(-1))
  expect_equal(out$ci_high, exp(1))
  expect_equal(out$scale, "natural")
})

test_that("bias differs between the two scales for a log link", {
  # This asymmetry is the whole reason `scale` is an argument: a bias of
  # zero on the log scale is not a bias of zero on the natural scale
  # unless the estimate is exactly right.
  estimates <- fake_estimates(c("kappa", "kappa"),
    estimate = c(1, -1),
    replication = c(1L, 2L)
  )
  truth <- fake_truth("kappa", true_value = 0)

  on_link <- recover(estimates, truth, scale = "link")
  on_natural <- recover(estimates, truth,
    scale = "natural",
    links = c(kappa = "log")
  )

  expect_equal(sum(on_link$bias), 0)
  expect_false(isTRUE(all.equal(sum(on_natural$bias), 0)))
  expect_equal(on_natural$bias, exp(c(1, -1)) - 1)
})

test_that("coverage is identical on both scales for every link", {
  # A monotone transform cannot move truth from inside an interval to
  # outside it. Values are kept positive so that every link in the
  # vocabulary is monotone over the range used --- `inverse` and `sqrt`
  # are not monotone across zero, and scoring across zero on those links
  # is a modelling error, not something coverage should paper over.
  withr::local_seed(4)
  n <- 30L
  # both bounds stay strictly positive, so `inverse` and `sqrt` are
  # monotone over the whole range used
  estimate <- stats::runif(n, 0.6, 2)
  truth <- stats::runif(n, 0.2, 2)
  estimates <- fake_estimates(
    rep("p", n),
    estimate = estimate,
    ci_low = estimate - 0.3, ci_high = estimate + 0.3,
    replication = seq_len(n)
  )
  truth_tbl <- fake_truth(rep("p", n),
    true_value = truth,
    replication = seq_len(n)
  )

  on_link <- recover(estimates, truth_tbl, scale = "link")
  for (link in link_vocabulary()) {
    on_natural <- recover(
      estimates, truth_tbl,
      scale = "natural", links = c(p = link)
    )
    expect_equal(on_natural$covered, on_link$covered, info = link)
  }
})

test_that("a sqrt-link interval spanning zero becomes [0, max]", {
  # x^2 attains its minimum, 0, inside a link-scale interval that spans
  # zero, so the image of [-1, 0.5] is [0, 1], not [0.25, 1]. A truth
  # near zero must be covered.
  estimates <- fake_estimates("p", estimate = 0, ci_low = -1, ci_high = 0.5)
  truth <- fake_truth("p", true_value = 1e-4)

  out <- recover(estimates, truth, scale = "natural", links = c(p = "sqrt"))
  expect_equal(out$ci_low, 0)
  expect_equal(out$ci_high, 1)
  expect_true(out$covered)

  # an interval on one side of zero is unaffected
  one_sided <- fake_estimates("p", estimate = 1, ci_low = 0.5, ci_high = 2)
  out <- recover(one_sided, truth, scale = "natural", links = c(p = "sqrt"))
  expect_equal(c(out$ci_low, out$ci_high), c(0.25, 4))
})

test_that("an inverse-link interval spanning zero is NA and warns", {
  # 1 / x diverges inside a link-scale interval that spans zero, so the
  # image is not an interval. The bounds and coverage become NA rather
  # than a plausible-looking wrong interval, and the caller is told.
  estimates <- fake_estimates(
    c("p", "q"),
    estimate = c(0, 1), ci_low = c(-0.5, 0.5), ci_high = c(0.5, 2)
  )
  truth <- fake_truth(c("p", "q"), true_value = c(100, 1))

  expect_warning(
    out <- recover(estimates, truth,
      scale = "natural", links = c(p = "inverse", q = "inverse")
    ),
    "p"
  )
  p <- out[out$term == "p", ]
  q <- out[out$term == "q", ]
  expect_true(is.na(p$ci_low) && is.na(p$ci_high) && is.na(p$covered))
  expect_equal(c(q$ci_low, q$ci_high), c(0.5, 2))
  expect_true(q$covered)

  # the summary still runs: coverage is computed over the defined rows
  s <- summary(out)
  expect_true(is.na(s$coverage[s$term == "p"]))
  expect_equal(s$coverage[s$term == "q"], 1)
})

test_that("a malformed links argument errors", {
  estimates <- fake_estimates("p", estimate = 1)
  truth <- fake_truth("p", true_value = 1)
  expect_error(
    recover(estimates, truth, scale = "natural", links = "log"),
    "named character"
  )
  expect_error(
    recover(estimates, truth, scale = "natural", links = c(log = 1)),
    "named character"
  )
})

test_that("truth must be a data frame", {
  estimates <- fake_estimates("p", estimate = 1)
  expect_error(recover(estimates, c(p = 1)), "data frame")
})

test_that("a decreasing link leaves ci_low below ci_high", {
  # inverse and loglog reverse the order of the bounds. Reordering them
  # after the transform is what makes coverage invariant; leaving them
  # swapped would make every interval empty.
  estimates <- fake_estimates("p", estimate = 1, ci_low = 0.5, ci_high = 2)
  truth <- fake_truth("p", true_value = 1)

  for (link in c("inverse", "loglog")) {
    out <- recover(estimates, truth, scale = "natural", links = c(p = link))
    expect_lt(out$ci_low, out$ci_high, label = link)
    expect_true(out$covered)
  }
})

test_that("a term with no link entry is transformed as identity", {
  estimates <- fake_estimates(c("kappa", "other"), estimate = c(0, 3))
  truth <- fake_truth(c("kappa", "other"), true_value = c(0, 3))

  out <- recover(estimates, truth,
    scale = "natural",
    links = c(kappa = "log")
  )

  expect_equal(out$estimate, c(1, 3))
})

test_that("no links and no model falls back to the link scale, once", {
  estimates <- fake_estimates(c("a", "b"), estimate = c(1, 2))
  truth <- fake_truth(c("a", "b"), true_value = c(1, 2))

  expect_message(
    out <- recover(estimates, truth, scale = "natural"),
    "link scale"
  )
  expect_equal(out$scale, c("link", "link"))
  expect_equal(out$estimate, c(1, 2))
})

test_that("scale = link needs no links and stays silent", {
  estimates <- fake_estimates("a", estimate = 1)
  truth <- fake_truth("a", true_value = 1)
  expect_silent(recover(estimates, truth, scale = "link"))
})

# missing terms ---------------------------------------------------------

test_that("a term in truth that the fit does not have warns and is dropped", {
  estimates <- fake_estimates(c("kappa", "thetat"), estimate = c(1, 2))
  truth <- fake_truth(c("kappa", "thetat", "kappa_typo"),
    true_value = c(1, 2, 3)
  )

  warn <- expect_warning(out <- recover(estimates, truth, scale = "link"))
  expect_match(conditionMessage(warn), "kappa_typo")
  expect_match(conditionMessage(warn), "thetat")
  expect_equal(nrow(out), 2L)
  expect_setequal(out$term, c("kappa", "thetat"))
})

test_that("every term missing is an error, not an empty tibble", {
  # It almost always means a naming mismatch, and a zero-row result hides
  # it until a plot comes out blank.
  estimates <- fake_estimates(c("kappa", "thetat"), estimate = c(1, 2))
  truth <- fake_truth(c("k", "t"), true_value = c(1, 2))

  err <- expect_error(recover(estimates, truth, scale = "link"))
  expect_match(conditionMessage(err), "kappa")
})

test_that("a term the fit has but truth does not is dropped silently", {
  # Fits routinely estimate more than a simulation grid varies.
  estimates <- fake_estimates(c("kappa", "thetat", "mu1"),
    estimate = c(1, 2, 0)
  )
  truth <- fake_truth(c("kappa", "thetat"), true_value = c(1, 2))

  expect_silent(out <- recover(estimates, truth, scale = "link"))
  expect_equal(nrow(out), 2L)
})

test_that("recover rejects a truth tibble without the required columns", {
  estimates <- fake_estimates("a", estimate = 1)
  expect_error(
    recover(estimates, tibble::tibble(term = "a"), scale = "link"),
    "true_value"
  )
  expect_error(
    recover(estimates, tibble::tibble(true_value = 1), scale = "link"),
    "term"
  )
})

test_that("recover rejects an object that is neither fit nor estimates", {
  expect_error(recover(1:10, fake_truth("a", 1)), "estimates tibble")
})

# summary ---------------------------------------------------------------

test_that("summary returns one row per term and level with the contract", {
  terms <- c("kappa", "thetat")
  estimates <- dplyr::bind_rows(lapply(1:5, function(rep) {
    fake_estimates(terms, estimate = c(1, 2) + rep / 10, replication = rep)
  }))
  truth <- fake_truth(terms, true_value = c(1, 2))

  out <- summary(recover(estimates, truth, scale = "link"))

  expect_s3_class(out, "bmmtools_recovery_summary")
  expect_equal(nrow(out), 2L)
  expect_named(out, recovery_summary_columns())
  expect_setequal(out$term, terms)
  expect_true(all(out$level == "population"))
  expect_true(all(out$scale == "link"))
  expect_equal(out$n, c(5, 5))
  expect_equal(out$n_replications, c(5L, 5L))
})

test_that("summary metrics equal the metric functions on the same columns", {
  withr::local_seed(7)
  n <- 12L
  true_value <- stats::rnorm(n)
  estimate <- true_value * 0.9 + stats::rnorm(n, sd = 0.2)
  estimates <- fake_estimates(
    rep("kappa", n),
    estimate = estimate,
    ci_low = estimate - 0.5, ci_high = estimate + 0.5,
    replication = seq_len(n)
  )
  truth <- fake_truth(rep("kappa", n),
    true_value = true_value,
    replication = seq_len(n)
  )

  rec <- recover(estimates, truth, scale = "link")
  out <- summary(rec)

  expect_equal(out$bias, metric_bias(rec$estimate, rec$true_value))
  expect_equal(out$rmse, metric_rmse(rec$estimate, rec$true_value))
  expect_equal(
    out$coverage, metric_coverage(rec$true_value, rec$ci_low, rec$ci_high)
  )
  expect_equal(out$ci_width, metric_ci_width(rec$ci_low, rec$ci_high))

  r <- metric_r(rec$estimate, rec$true_value)
  expect_equal(out$r, r$r)
  expect_equal(out$r_low, r$r_low)
  expect_equal(out$r_high, r$r_high)
  expect_equal(out$rank_r, metric_rank_r(rec$estimate, rec$true_value))

  ccc <- metric_ccc(rec$estimate, rec$true_value)
  expect_equal(out$ccc, ccc$ccc)
  expect_equal(out$ccc_scale_bias, ccc$scale_bias)
  expect_equal(out$ccc_location_bias, ccc$location_bias)
})

test_that("n_converged is NA when no fit carried a convergence flag", {
  # A hand-built estimates tibble has no `converged` column; reporting
  # n_converged as n would assert something that was never measured.
  estimates <- fake_estimates("a", estimate = 1)
  truth <- fake_truth("a", true_value = 1)
  out <- recover(estimates, truth, scale = "link")
  expect_true("converged" %in% names(out))
  expect_identical(out$converged, NA)
  expect_true(is.na(summary(out)$n_converged))
})

test_that("n_converged counts the replications whose fit passed the gate", {
  estimates <- dplyr::bind_rows(lapply(1:3, function(rep) {
    out <- fake_estimates(c("a", "b"),
      estimate = c(1, 2) + rep / 10,
      replication = rep
    )
    out$converged <- c(TRUE, TRUE, FALSE)[[rep]]
    out
  }))
  truth <- fake_truth(c("a", "b"), true_value = c(1, 2))

  out <- summary(recover(estimates, truth, scale = "link"))
  expect_identical(out$n_converged, c(2L, 2L))
  expect_identical(out$n_replications, c(3L, 3L))

  # the subject level counts replications the same way
  subjects <- dplyr::bind_rows(lapply(1:3, function(rep) {
    out <- fake_estimates(rep("a", 5),
      estimate = seq_len(5) + rep / 10,
      level = "subject", id = as.character(1:5), replication = rep
    )
    out$converged <- c(TRUE, FALSE, FALSE)[[rep]]
    out
  }))
  truth_subjects <- fake_truth(rep("a", 5),
    true_value = seq_len(5),
    id = as.character(1:5)
  )
  out <- summary(recover_subjects(subjects, truth_subjects, scale = "link"))
  expect_identical(out$n_converged, 1L)
})

test_that("a single fit gives NA correlations and numeric bias", {
  estimates <- fake_estimates(c("kappa", "thetat"), estimate = c(1.2, 2.1))
  truth <- fake_truth(c("kappa", "thetat"), true_value = c(1, 2))
  out <- summary(recover(estimates, truth, scale = "link"))

  expect_true(all(is.na(out$r)))
  expect_true(all(is.na(out$rank_r)))
  expect_true(all(is.na(out$ccc)))
  expect_false(any(is.na(out$bias)))
  expect_false(any(is.na(out$rmse)))
  expect_equal(out$n, c(1, 1))
})

# recover_subjects ------------------------------------------------------

subject_case <- function(n_subjects = 20L, n_replications = 1L, seed = 12) {
  withr::local_seed(seed)
  terms <- c("kappa", "thetat")
  ids <- as.character(seq_len(n_subjects))
  estimates <- list()
  truth <- list()
  for (rep in seq_len(n_replications)) {
    for (term in terms) {
      true_value <- stats::rnorm(n_subjects)
      estimate <- true_value * 0.8 + stats::rnorm(n_subjects, sd = 0.4)
      estimates[[length(estimates) + 1L]] <- fake_estimates(
        rep(term, n_subjects),
        estimate = estimate,
        ci_low = estimate - 1, ci_high = estimate + 1,
        level = "subject", id = ids, replication = rep
      )
      truth[[length(truth) + 1L]] <- tibble::tibble(
        id = ids, term = term, true_value = true_value, replication = rep
      )
    }
  }
  list(
    estimates = dplyr::bind_rows(estimates),
    truth = dplyr::bind_rows(truth)
  )
}

test_that("recover_subjects returns one row per subject, term and rep", {
  case <- subject_case(n_subjects = 20L)
  out <- recover_subjects(case$estimates, case$truth, scale = "link")

  expect_equal(nrow(out), 40L)
  expect_named(out, recovery_contract_columns())
  expect_true(all(out$level == "subject"))
  expect_false(any(is.na(out$id)))
  expect_equal(length(unique(out$id)), 20L)
})

test_that("subject-level r and its interval match cor.test per term", {
  case <- subject_case(n_subjects = 20L)
  rec <- recover_subjects(case$estimates, case$truth, scale = "link")
  out <- summary(rec)

  expect_equal(nrow(out), 2L)
  for (term in c("kappa", "thetat")) {
    rows <- rec[rec$term == term, ]
    oracle <- stats::cor.test(rows$estimate, rows$true_value)
    got <- out[out$term == term, ]
    expect_equal(got$r, unname(oracle$estimate), tolerance = 1e-8)
    expect_equal(got$r_low, oracle$conf.int[1], tolerance = 1e-8)
    expect_equal(got$r_high, oracle$conf.int[2], tolerance = 1e-8)
    expect_equal(got$n, 20)
  }
  expect_equal(unique(out$n_replications), 1L)
})

test_that("two replications average r via Fisher-z, not by pooling", {
  case <- subject_case(n_subjects = 20L, n_replications = 2L, seed = 21)
  rec <- recover_subjects(case$estimates, case$truth, scale = "link")
  out <- summary(rec)

  term <- "kappa"
  rows <- rec[rec$term == term, ]
  within <- vapply(1:2, function(rep) {
    sub <- rows[rows$replication == rep, ]
    stats::cor(sub$estimate, sub$true_value)
  }, numeric(1))
  # equal subject counts, so the inverse-variance weights are equal and
  # the combined value is the plain Fisher-z average
  fisher_z_average <- tanh(mean(atanh(within)))
  pooled <- stats::cor(rows$estimate, rows$true_value)

  got <- out[out$term == term, ]
  expect_equal(got$r, fisher_z_average, tolerance = 1e-8)
  expect_equal(got$n, 20)
  expect_equal(got$n_replications, 2L)
  expect_false(isTRUE(all.equal(got$r, pooled)))
})

test_that("pooling and Fisher-z averaging differ when they must", {
  # Constructed so the two answers cannot coincide: within each
  # replication the estimates track truth perfectly, but the two
  # replications sit on different offsets, so pooling across them
  # destroys the correlation that each replication has.
  ids <- as.character(1:5)
  estimates <- dplyr::bind_rows(
    fake_estimates(rep("p", 5),
      estimate = c(1, 2, 3, 4, 5),
      level = "subject", id = ids, replication = 1L
    ),
    fake_estimates(rep("p", 5),
      estimate = c(5, 4, 3, 2, 1),
      level = "subject", id = ids, replication = 2L
    )
  )
  truth <- dplyr::bind_rows(
    tibble::tibble(
      id = ids, term = "p", true_value = c(1, 2, 3, 4, 5),
      replication = 1L
    ),
    tibble::tibble(
      id = ids, term = "p", true_value = c(5, 4, 3, 2, 1),
      replication = 2L
    )
  )

  out <- summary(recover_subjects(estimates, truth, scale = "link"))
  expect_equal(out$r, 1)
  expect_equal(out$n_replications, 2L)
})

test_that("subject-level coverage lies in the unit interval", {
  case <- subject_case(n_subjects = 20L, n_replications = 2L)
  out <- summary(recover_subjects(case$estimates, case$truth, scale = "link"))
  expect_true(all(out$coverage >= 0 & out$coverage <= 1))
})

test_that("a subject in truth that the fit does not have warns", {
  case <- subject_case(n_subjects = 5L)
  extra <- tibble::tibble(
    id = "99", term = "kappa", true_value = 0, replication = 1L
  )
  truth <- dplyr::bind_rows(case$truth, extra)

  warn <- expect_warning(
    out <- recover_subjects(case$estimates, truth, scale = "link")
  )
  expect_match(conditionMessage(warn), "99")
  expect_equal(nrow(out), 10L)
})

test_that("zero spread in subject truth gives NA r, not zero", {
  ids <- as.character(1:5)
  estimates <- fake_estimates(rep("p", 5),
    estimate = c(1, 2, 3, 4, 5),
    level = "subject", id = ids
  )
  truth <- tibble::tibble(id = ids, term = "p", true_value = rep(2, 5))

  out <- summary(recover_subjects(estimates, truth, scale = "link"))
  expect_true(is.na(out$r))
  expect_false(isTRUE(out$r == 0))
  expect_false(is.na(out$bias))
})

test_that("recover_subjects rejects truth without an id column", {
  case <- subject_case(n_subjects = 5L)
  expect_error(
    recover_subjects(case$estimates, case$truth[c("term", "true_value")]),
    "id"
  )
})

# fits as input ---------------------------------------------------------

test_that("recover scores a saved bmmfit and reads its links", {
  skip_if_not_installed("brms")
  skip_if_not_installed("bmm")
  fit <- mixture2p_fit()
  truth <- mixture2p_truth()$population

  out <- recover(fit, truth)

  expect_s3_class(out, "bmmtools_recovery")
  expect_setequal(out$term, c("kappa", "thetat"))
  expect_true(all(out$scale == "natural"))
  # kappa has a log link, so the natural-scale estimate is the exponential
  # of the link-scale posterior median
  link_scale <- extract_estimates(fit)
  expect_equal(
    out$estimate[out$term == "kappa"],
    exp(link_scale$estimate[link_scale$term == "kappa"])
  )
  expect_equal(
    out$true_value[out$term == "kappa"],
    exp(truth$true_value[truth$term == "kappa"])
  )
})

test_that("recover_subjects scores a saved bmmfit against subject truth", {
  skip_if_not_installed("brms")
  skip_if_not_installed("bmm")
  fit <- mixture2p_fit()
  truth <- mixture2p_truth()$subjects

  out <- recover_subjects(fit, truth, scale = "link")

  expect_equal(nrow(out), 16L)
  expect_setequal(out$id, as.character(1:8))
  summarised <- summary(out)
  expect_equal(summarised$n, c(8, 8))
  expect_false(any(is.na(summarised$r)))
})

test_that("a list of fits becomes one replication each", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  truth <- mixture2p_truth()$population

  out <- recover(list(fit, fit), truth, scale = "link")

  expect_equal(nrow(out), 4L)
  expect_setequal(out$replication, 1:2)
  # the same fit twice: identical estimates, so a perfect but degenerate
  # correlation guard applies at n = 2
  expect_true(is.na(summary(out)$r[1]))
})
