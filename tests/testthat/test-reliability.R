# Reliability of subject estimates and the MCSEs of the rates
# (Milestone 24.1, D96 and D98). Every expected value is a formula from the
# decision applied to hand-made rows, or the normal-normal identity of the
# spec's section 1b; no test compiles Stan.

#' Subject rows of one term with link-scale posterior moments
#' @noRd
moment_subjects <- function(mean, sd, truth, replication, term = "kappa",
                            estimator = "bayes") {
  id <- as.character(stats::ave(seq_along(mean), replication, FUN = seq_along))
  estimates <- fake_estimates(term,
    estimate = mean, level = "subject", id = id, replication = replication
  )
  estimates$post_mean_link <- as.double(mean)
  estimates$post_sd_link <- as.double(sd)
  estimates$estimator <- estimator
  list(
    estimates = estimates,
    truth = fake_truth(term, truth, id = id, replication = replication)
  )
}

subject_summary <- function(rows, scale = "link", links = NULL) {
  summary(recover_subjects(rows$estimates, rows$truth,
    scale = scale, links = links
  ))
}

#' True values T ~ N(0, 1), observations x = T + e with e ~ N(0, s^2), and
#' w = 1 / (1 + s^2) the design's reliability. "bayes" rows carry the
#' sigma-known posterior, mean w * x and SD sqrt(w) * s; "ml" rows carry
#' the unshrunk x and its standard error s.
#' @noRd
normal_normal <- function(n, reps, rho2, estimator = "bayes", seed) {
  withr::local_seed(seed)
  s <- sqrt((1 - rho2) / rho2)
  replication <- rep(seq_len(reps), each = n)
  truth <- stats::rnorm(n * reps)
  x <- truth + stats::rnorm(n * reps, 0, s)
  if (identical(estimator, "ml")) {
    moment_subjects(x, s, truth, replication, estimator = "ml")
  } else {
    moment_subjects(rho2 * x, sqrt(rho2) * s, truth, replication)
  }
}

# the contract ---------------------------------------------------------------

planning_columns <- c(
  "reliability_true", "reliability_true_mcse", "reliability_model",
  "reliability_model_mcse", "reliability_gap", "reliability_gap_mcse",
  "coverage_mcse", "coverage_50_mcse", "detected_mcse", "sign_recovery_mcse"
)

test_that("the reliability and rate MCSEs follow z_sd_mcse in the contract", {
  columns <- recovery_summary_columns()
  at <- match("z_sd_mcse", columns)
  expect_identical(columns[at + seq_along(planning_columns)], planning_columns)
})

test_that("an empty summary carries the new columns as doubles", {
  empty <- empty_recovery_summary()
  expect_named(empty, recovery_summary_columns())
  for (column in planning_columns) expect_type(empty[[column]], "double")
})

# reliability, by hand ---------------------------------------------------------

test_that("both reliabilities are per-replication means, by hand", {
  m <- list(c(0.1, 0.5, 0.9, 1.6), c(-0.2, 0.4, 0.3, 1.1))
  s <- list(c(0.3, 0.3, 0.4, 0.5), c(0.2, 0.2, 0.3, 0.3))
  t <- list(c(0, 0.7, 0.8, 1.9), c(-0.5, 0.2, 0.6, 1.0))
  rows <- moment_subjects(unlist(m), unlist(s), unlist(t), rep(1:2, each = 4))
  out <- subject_summary(rows)

  a <- c(stats::cor(m[[1]], t[[1]])^2, stats::cor(m[[2]], t[[2]])^2)
  b <- c(
    stats::var(m[[1]]) / (stats::var(m[[1]]) + mean(s[[1]]^2)),
    stats::var(m[[2]]) / (stats::var(m[[2]]) + mean(s[[2]]^2))
  )
  expect_equal(out$reliability_true, mean(a))
  expect_equal(out$reliability_true_mcse, stats::sd(a) / sqrt(2))
  expect_equal(out$reliability_model, mean(b))
  expect_equal(out$reliability_model_mcse, stats::sd(b) / sqrt(2))
  # the gap is paired by replication; its MCSE is not the marginal ones'
  expect_equal(out$reliability_gap, mean(b - a))
  expect_equal(out$reliability_gap_mcse, stats::sd(b - a) / sqrt(2))
})

test_that("the posterior mean is used, not the median in `estimate`", {
  rows <- moment_subjects(
    c(0.1, 0.5, 0.9, 1.6), 0.3, c(0, 0.7, 0.8, 1.9),
    replication = 1
  )
  # a median that orders the subjects backwards must not change anything
  rows$estimates$estimate <- rev(rows$estimates$estimate)
  out <- subject_summary(rows)
  expect_equal(
    out$reliability_true,
    stats::cor(c(0.1, 0.5, 0.9, 1.6), c(0, 0.7, 0.8, 1.9))^2
  )
})

test_that("an ML row takes the separation form, not clamped at 0", {
  m <- c(0, 0.1, 0.2, 0.3)
  rows <- moment_subjects(m, 1, c(0, 0.1, 0.2, 0.3), 1, estimator = "ml")
  out <- subject_summary(rows)
  expect_identical(out$estimator, "ml")
  expect_equal(out$reliability_model, 1 - 1 / stats::var(m))
  expect_lt(out$reliability_model, 0)
})

test_that("reliability is on the link scale whatever `scale` says", {
  # the truth handed to recover_subjects() is on the link scale; scored
  # natural, it is reconstructed from post_mean_link - z * post_sd_link
  rows <- moment_subjects(
    c(1.1, 1.6, 2.0, 2.4, 0.8, 1.9, 2.2, 2.9), 0.25,
    c(1.0, 1.8, 1.9, 2.6, 1.0, 1.7, 2.5, 2.7),
    replication = rep(1:2, each = 4)
  )
  linked <- subject_summary(rows)
  natural <- subject_summary(rows, scale = "natural", links = c(kappa = "log"))
  expect_identical(natural$scale, "natural")
  reliability <- planning_columns[1:6]
  expect_equal(as.list(natural[reliability]), as.list(linked[reliability]))
  # the row's own correlation is on its scale, so it moved
  expect_false(isTRUE(all.equal(natural$r, linked$r)))
})

# reliability, the normal-normal identity (spec section 1b) ----------------

test_that("normal-normal: both estimators recover the design's reliability", {
  # both are ratios of a sample variance, so each has a finite-n bias; at
  # 400 subjects it is a tenth of the MCSE for the EAP form (a 20,000-draw
  # Monte Carlo at 100 subjects put it at -0.002, at rho2 = 0.5)
  for (rho2 in c(0.5, 0.8)) {
    out <- subject_summary(normal_normal(400, 50, rho2, seed = 2401))
    expect_lt(abs(out$reliability_true - rho2), 3 * out$reliability_true_mcse)
    expect_lt(abs(out$reliability_model - rho2), 3 * out$reliability_model_mcse)
    expect_lt(abs(out$reliability_gap), 3 * out$reliability_gap_mcse)
    expect_equal(out$reliability_true, rho2, tolerance = 0.05)
    expect_equal(out$reliability_model, rho2, tolerance = 0.05)
  }
})

test_that("normal-normal: unshrunk ML estimates need the separation form", {
  rho2 <- 0.5
  n <- 400
  rows <- normal_normal(n, 50, rho2, estimator = "ml", seed = 2402)
  ml <- subject_summary(rows)
  # var(x) = (1 + s^2) chi^2_(n-1) / (n - 1), so E[1 / var(x)] carries
  # (n - 1) / (n - 3), and the separation form's expectation is exact
  expected <- 1 - (1 - rho2) * (n - 1) / (n - 3)
  expect_lt(abs(ml$reliability_model - expected), 3 * ml$reliability_model_mcse)
  expect_equal(ml$reliability_model, rho2, tolerance = 0.05)

  # the same rows read with the posterior-mean formula report about 2/3:
  # var(x) / (var(x) + s^2) = (1 + s^2) / (1 + 2 s^2) at s^2 = 1
  rows$estimates$estimator <- "bayes"
  eap <- subject_summary(rows)
  expect_equal(eap$reliability_model, 2 / 3, tolerance = 0.05)
  expect_gt(eap$reliability_model - ml$reliability_model, 0.1)
})

# reliability, NA rules -------------------------------------------------------

test_that("a replication with fewer than three complete subjects is left out", {
  m <- c(0.1, 0.5, 0.9, 1.6)
  t <- c(0, 0.7, 0.8, 1.9)
  rows <- moment_subjects(c(m, 0.2, 0.4), 0.3, c(t, 0.1, 0.6),
    replication = c(1, 1, 1, 1, 2, 2)
  )
  out <- subject_summary(rows)
  expect_equal(out$reliability_true, stats::cor(m, t)^2)
  expect_equal(
    out$reliability_model,
    stats::var(m) / (stats::var(m) + 0.09)
  )
  # one contributing replication: no Monte Carlo error to report
  expect_identical(
    c(
      out$reliability_true_mcse, out$reliability_model_mcse,
      out$reliability_gap_mcse
    ),
    rep(NA_real_, 3)
  )
})

test_that("a replication without spread contributes nothing", {
  m <- c(0.1, 0.5, 0.9, 1.6)
  t <- c(0, 0.7, 0.8, 1.9)
  # replication 2: identical posterior means; replication 3: identical truth
  rows <- moment_subjects(
    c(m, rep(0.5, 4), m), 0.3, c(t, t, rep(1, 4)),
    replication = rep(1:3, each = 4)
  )
  out <- subject_summary(rows)
  b <- stats::var(m) / (stats::var(m) + 0.09)
  expect_equal(out$reliability_true, stats::cor(m, t)^2)
  # the model form does not need the truth, so replication 3 counts
  expect_equal(out$reliability_model, b)
  expect_equal(out$reliability_model_mcse, 0)
  expect_equal(out$reliability_gap, b - stats::cor(m, t)^2)
  expect_identical(out$reliability_gap_mcse, NA_real_)
})

test_that("without a posterior SD the model form is NA", {
  m <- c(0.1, 0.5, 0.9, 1.6)
  t <- c(0, 0.7, 0.8, 1.9)
  rows <- moment_subjects(m, NA, t, replication = 1)
  linked <- subject_summary(rows)
  expect_equal(linked$reliability_true, stats::cor(m, t)^2)
  expect_identical(
    c(linked$reliability_model, linked$reliability_gap),
    rep(NA_real_, 2)
  )
  # scored natural, the link truth cannot be rebuilt without z either
  natural <- subject_summary(rows, scale = "natural", links = c(kappa = "log"))
  expect_identical(natural$reliability_true, NA_real_)
})

test_that("population and sd rows have no reliability", {
  estimates <- dplyr::bind_rows(
    fake_estimates("kappa", c(1, 2, 3), replication = 1:3),
    fake_estimates("kappa", c(0.5, 0.7, 0.9), level = "sd", replication = 1:3)
  )
  estimates$post_mean_link <- estimates$estimate
  estimates$post_sd_link <- 0.5
  truth <- list(
    population = fake_truth("kappa", c(1.1, 2.2, 2.9), replication = 1:3),
    sd = fake_truth("kappa", c(0.6, 0.6, 0.8), replication = 1:3)
  )
  out <- summary(recover(estimates, truth,
    level = c("population", "sd"), scale = "link"
  ))
  expect_identical(out$level, c("population", "sd"))
  expect_identical(
    unlist(out[planning_columns[1:6]], use.names = FALSE),
    rep(NA_real_, 12)
  )
})

test_that("a group that mixes ci_levels has no nominal coverage MCSE", {
  estimates <- fake_estimates("kappa", rep(1, 4),
    replication = 1:4, ci_level = c(0.95, 0.95, 0.9, 0.9)
  )
  out <- summary(recover(estimates, fake_truth("kappa", 1), scale = "link"))
  expect_identical(out$coverage_mcse, NA_real_)
})

test_that("the shipped example gives the spec's independent numbers", {
  # spec 24, section 1a: computed in session 24.0 by a scratch script that
  # did not use summary(); values there are rounded to three decimals
  s <- summary(recovery_mixture2p)
  at <- function(condition, term) {
    s[s$condition == condition & s$term == term & s$level == "subject", ]
  }
  row1 <- at("row-1", "kappa")
  expect_equal(row1$reliability_true, 0.333, tolerance = 5e-4 / 0.333)
  expect_equal(row1$reliability_true_mcse, 0.103, tolerance = 5e-4 / 0.103)
  expect_equal(row1$reliability_model, 0.200, tolerance = 5e-4 / 0.2)
  expect_equal(row1$reliability_gap, -0.133, tolerance = 5e-4 / 0.133)
  row4 <- at("row-4", "kappa")
  expect_equal(row4$reliability_gap, -0.034, tolerance = 5e-4 / 0.034)
  expect_equal(row4$reliability_gap_mcse, 0.007, tolerance = 5e-4 / 0.007)
  thetat <- at("row-4", "thetat")
  expect_equal(thetat$reliability_true, 0.703, tolerance = 5e-4 / 0.703)
  expect_equal(thetat$reliability_model, 0.697, tolerance = 5e-4 / 0.697)
})

# MCSEs of the rates --------------------------------------------------------

test_that("population coverage MCSEs are binomial at the nominal rate", {
  # ten replications, all covered: the MCSE at the observed rate would be 0
  estimates <- fake_estimates("kappa", rep(1, 10), replication = 1:10)
  estimates$ci_low_50 <- 0.9
  estimates$ci_high_50 <- 1.1
  out <- summary(recover(estimates, fake_truth("kappa", 1), scale = "link"))
  expect_equal(out$coverage, 1)
  expect_equal(out$coverage_mcse, sqrt(0.95 * 0.05 / 10))
  expect_equal(out$coverage_50_mcse, sqrt(0.5 * 0.5 / 10))
})

test_that("the coverage MCSE follows ci_level and counts complete rows", {
  estimates <- fake_estimates("kappa", rep(1, 8),
    replication = 1:8, ci_level = 0.9
  )
  estimates$ci_low[1:2] <- NA
  out <- summary(recover(estimates, fake_truth("kappa", 1),
    scale = "link", ci_level = 0.9
  ))
  expect_equal(out$coverage_mcse, sqrt(0.9 * 0.1 / 6))
  # no 50 % interval: no rate and no MCSE
  expect_identical(c(out$coverage_50, out$coverage_50_mcse), c(NA_real_, NA))
})

test_that("detected and sign_recovery MCSEs are binomial at the observed rate", {
  # intervals 1 wide: 3 of 5 exclude zero; truth 0.4, signs right 4 of 5
  estimate <- c(2, 1.5, -2, 0.2, 0.3)
  estimates <- fake_estimates("effect", estimate,
    ci_low = estimate - 0.5, ci_high = estimate + 0.5, replication = 1:5
  )
  out <- summary(recover(estimates, fake_truth("effect", 0.4), scale = "link"))
  expect_equal(out$detected, 0.6)
  expect_equal(out$detected_mcse, sqrt(0.6 * 0.4 / 5))
  expect_equal(out$sign_recovery, 0.8)
  expect_equal(out$sign_recovery_mcse, sqrt(0.8 * 0.2 / 5))
})

test_that("sign_recovery_mcse counts only rows whose truth has a sign", {
  estimates <- fake_estimates("effect", c(1, 2, -1, 3), replication = 1:4)
  truth <- fake_truth("effect", c(1, 1, 1, 0), replication = 1:4)
  out <- summary(recover(estimates, truth, scale = "link"))
  expect_equal(out$sign_recovery, 2 / 3)
  expect_equal(out$sign_recovery_mcse, sqrt(2 / 3 * 1 / 3 / 3))
})

test_that("subject coverage MCSEs are between replications", {
  # covered per replication: 4/4, 2/4, 3/4; truth 0, intervals +/- 1
  estimate <- c(0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, 2)
  rows <- moment_subjects(estimate, 0.5, 0, rep(1:3, each = 4))
  # 50 % intervals: covered 4/4, 4/4, 2/4 (+/- 0.1, except two subjects of
  # replication 3 whose 50 % interval sits off the truth)
  rows$estimates$ci_low_50 <- c(rep(-0.1, 10), 0.5, 0.5)
  rows$estimates$ci_high_50 <- c(rep(0.1, 10), 0.7, 0.7)
  out <- subject_summary(rows)
  shares <- c(1, 0.5, 0.75)
  shares_50 <- c(1, 1, 0.5)
  expect_equal(out$coverage, mean(shares))
  expect_equal(out$coverage_mcse, stats::sd(shares) / sqrt(3))
  expect_equal(out$coverage_50_mcse, stats::sd(shares_50) / sqrt(3))
  expect_identical(c(out$detected_mcse, out$sign_recovery_mcse), c(NA, NA_real_))
})

test_that("one replication of subjects has no coverage MCSE", {
  rows <- moment_subjects(c(0, 0, 2), 0.5, 0, 1)
  out <- subject_summary(rows)
  expect_identical(out$coverage_mcse, NA_real_)
})

# printing ------------------------------------------------------------------

test_that("the print adds both reliabilities when a subject level is there", {
  rows <- moment_subjects(c(0.1, 0.5, 0.9), 0.3, c(0, 0.7, 0.8), 1)
  columns <- summary_print_columns(subject_summary(rows))
  at <- match("ccc", columns)
  expect_identical(
    columns[at + 1:2], c("reliability_true", "reliability_model")
  )
  expect_false(any(grepl("_mcse$|reliability_gap", columns)))

  population <- summary(recover(fake_estimates("kappa", 1),
    fake_truth("kappa", 1),
    scale = "link"
  ))
  expect_false(any(grepl("^reliability", summary_print_columns(population))))
})
