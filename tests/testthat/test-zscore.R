# Posterior z-scores and contraction (Milestone 22.2, D59-D66, D86-D90).
# Every case is analytic or a hand-built posterior; no test compiles Stan.

moment_estimates <- function(term, mean, sd, level = "population",
                             id = NA_character_, replication = 1L) {
  estimates <- fake_estimates(
    term,
    estimate = mean, level = level, id = id, replication = replication
  )
  estimates$post_mean_link <- as.double(mean)
  estimates$post_sd_link <- as.double(sd)
  estimates
}

prior_table <- function(term, prior_sd_link, level = "population", ...) {
  tibble::tibble(
    term = term, level = level, prior_sd_link = as.double(prior_sd_link),
    ...
  )
}

# per-row z ----------------------------------------------------------------

test_that("z and contraction follow the moments in the contract", {
  contract <- names(recovery_contract())
  expect_identical(
    contract[match("post_sd_link", contract) + 1:2],
    c("z", "contraction")
  )
})

test_that("a known-bias estimator has a known z", {
  # the posterior mean sits 1.5 posterior SDs above the truth
  estimates <- moment_estimates(c("a", "b"), c(1.75, -1), c(0.5, 2))
  out <- recover(estimates, fake_truth(c("a", "b"), c(1, 0)), scale = "link")
  expect_equal(out$z, c(1.5, -0.5))
})

test_that("z is computed on the link scale before the natural transform", {
  estimates <- moment_estimates("kappa", 2.3, 0.2)
  truth <- fake_truth("kappa", 2)
  linked <- recover(estimates, truth, scale = "link")
  natural <- recover(estimates, truth,
    scale = "natural", links = c(kappa = "log")
  )
  expect_equal(linked$z, 1.5)
  expect_identical(natural$z, linked$z)
  # the truth was transformed for the natural row; z was not recomputed
  expect_equal(natural$true_value, exp(2))
})

test_that("a zero or missing posterior SD gives z NA, never Inf", {
  estimates <- moment_estimates(c("a", "b", "c"), c(1, 1, 1), c(0, NA, 1))
  out <- recover(estimates, fake_truth(c("a", "b", "c"), c(0, 0, 0)),
    scale = "link"
  )
  expect_identical(out$z, c(NA_real_, NA_real_, 1))
})

test_that("estimates without moments score with z NA", {
  estimates <- fake_estimates("a", estimate = 1)
  out <- recover(estimates, fake_truth("a", 0), scale = "link")
  expect_identical(out$z, NA_real_)
  expect_identical(out$contraction, NA_real_)
})

test_that("softmax terms of mixture3p get a link-scale z, nothing natural", {
  terms <- c("kappa", "thetat", "thetant")
  estimates <- moment_estimates(terms, c(2.1, 1.3, -0.4), c(0.2, 0.3, 0.2))
  truth <- fake_truth(terms, c(2, 1.5, -0.5))
  links <- c(
    mu1 = "tan_half", kappa = "log", thetat = "softmax",
    thetant = "softmax"
  )
  out <- suppressMessages(
    recover(estimates, truth, scale = "natural", links = links)
  )
  expect_equal(out$z, c(0.5, -2 / 3, 0.5))
  expect_identical(out$scale, c("natural", "link", "link"))
  expect_equal(out$true_value[2:3], c(1.5, -0.5))
})

# per-row contraction --------------------------------------------------------

test_that("without prior_sd, contraction is NA", {
  estimates <- moment_estimates("a", 1, 0.5)
  out <- recover(estimates, fake_truth("a", 1), scale = "link")
  expect_identical(out$contraction, NA_real_)
})

test_that("conjugate normal: contraction is exact", {
  # theta ~ N(0, tau^2), n observations with known sigma: the posterior
  # variance is 1 / (1 / tau^2 + n / sigma^2), and contraction is
  # n tau^2 / (sigma^2 + n tau^2)
  tau <- 2
  sigma <- 1
  n <- 4
  post_sd <- sqrt(1 / (1 / tau^2 + n / sigma^2))
  estimates <- moment_estimates("mu", 0.3, post_sd)
  out <- recover(estimates, fake_truth("mu", 0.25),
    scale = "link", prior_sd = prior_table("mu", tau)
  )
  expect_equal(out$contraction, n * tau^2 / (sigma^2 + n * tau^2))
})

test_that("a posterior wider than its prior keeps a negative contraction", {
  estimates <- moment_estimates("a", 0, 2)
  out <- recover(estimates, fake_truth("a", 0),
    scale = "link", prior_sd = prior_table("a", 1)
  )
  expect_equal(out$contraction, -3)
})

test_that("a zero, missing or infinite prior SD gives contraction NA", {
  estimates <- moment_estimates(c("a", "b", "c"), c(0, 0, 0), c(1, 1, 1))
  expect_message(
    out <- recover(estimates, fake_truth(c("a", "b", "c"), c(0, 0, 0)),
      scale = "link", prior_sd = prior_table(c("a", "b", "c"), c(0, NA, Inf))
    ),
    "No prior SD for .*b \\(population\\)"
  )
  expect_identical(out$contraction, rep(NA_real_, 3))
})

test_that("prior_sd matches by term and level, and says once what is missing", {
  estimates <- dplyr::bind_rows(
    moment_estimates(c("a", "b"), c(0, 0), c(1, 1)),
    moment_estimates("a", 0.5, 0.1, level = "sd")
  )
  truth <- list(
    population = fake_truth(c("a", "b"), c(0, 0)),
    sd = fake_truth("a", 0.5)
  )
  prior <- prior_table(c("a", "a"), c(2, 0.2), level = c("population", "sd"))
  expect_message(
    out <- recover(estimates, truth,
      level = c("population", "sd"), scale = "link", prior_sd = prior
    ),
    "No prior SD for .*b \\(population\\)"
  )
  expect_equal(out$contraction, c(1 - 1 / 4, NA, 1 - 0.01 / 0.04))
})

test_that("a prior_sd table with a condition column joins per condition", {
  estimates <- moment_estimates(c("a", "a"), c(0, 0), c(1, 1))
  estimates$condition <- c("row-1", "row-2")
  truth <- fake_truth(c("a", "a"), c(0, 0), condition = c("row-1", "row-2"))
  prior <- prior_table(c("a", "a"), c(2, 4), condition = c("row-1", "row-2"))
  out <- recover(estimates, truth, scale = "link", prior_sd = prior)
  expect_equal(out$contraction, c(1 - 1 / 4, 1 - 1 / 16))
})

test_that("ML rows get z on the standard error and no contraction (D64)", {
  estimates <- moment_estimates("a", c(1.2, 1.4), c(0.2, 0.2),
    replication = 1:2
  )
  estimates$estimator <- c("bayes", "ml")
  out <- suppressWarnings(recover(estimates, fake_truth("a", 1),
    scale = "link", prior_sd = prior_table("a", 1)
  ))
  expect_equal(out$z, c(1, 2))
  expect_equal(out$contraction, c(1 - 0.04, NA))
})

test_that("prior_sd is refused in a shape it cannot use", {
  estimates <- moment_estimates("a", 0, 1)
  truth <- fake_truth("a", 0)
  expect_error(
    recover(estimates, truth, scale = "link", prior_sd = 2),
    "must be a data frame of prior SDs"
  )
  expect_error(
    recover(estimates, truth,
      scale = "link", prior_sd = tibble::tibble(term = "a", sd = 1)
    ),
    "needs the columns .*level.*prior_sd_link"
  )
})

test_that("recover_subjects() has no prior_sd and gives subjects NA (D63)", {
  expect_false("prior_sd" %in% names(formals(recover_subjects)))
  estimates <- moment_estimates(c("a", "a"), c(0.5, 1), c(0.5, 0.5),
    level = "subject", id = c("1", "2")
  )
  truth <- fake_truth(c("a", "a"), c(0, 1), id = c("1", "2"))
  out <- recover_subjects(estimates, truth, scale = "link")
  expect_equal(out$z, c(1, 0))
  expect_identical(out$contraction, c(NA_real_, NA_real_))
})

# the summary ----------------------------------------------------------------

test_that("the summary carries z_mean, z_sd, contraction and their MCSE", {
  expect_true(all(
    c("z_mean", "z_sd", "contraction", "z_mean_mcse", "z_sd_mcse") %in%
      recovery_summary_columns()
  ))
  z <- c(-1, 0.5, 2, 0.25)
  estimates <- moment_estimates("a", z, 1, replication = 1:4)
  out <- recover(estimates, fake_truth("a", 0),
    scale = "link", prior_sd = prior_table("a", 2)
  )
  s <- summary(out)
  expect_equal(s$z_mean, mean(z))
  expect_equal(s$z_sd, sd(z))
  expect_equal(s$contraction, 0.75)
  expect_equal(s$z_mean_mcse, sd(z) / sqrt(4))
  expect_equal(s$z_sd_mcse, sd(z) / sqrt(2 * 3))
})

test_that("the summary pools z and contraction over non-NA rows", {
  estimates <- moment_estimates("a", c(1, 3, 5), c(1, NA, 1),
    replication = 1:3
  )
  out <- recover(estimates, fake_truth("a", 0),
    scale = "link", prior_sd = prior_table("a", 2)
  )
  s <- summary(out)
  expect_equal(s$z_mean, 3)
  expect_equal(s$z_sd, sd(c(1, 5)))
  expect_equal(s$z_mean_mcse, sd(c(1, 5)) / sqrt(2))
  expect_equal(s$contraction, 0.75)
})

test_that("with fewer than two z, z_sd and both MCSE are NA", {
  out <- recover(moment_estimates("a", 1, 1), fake_truth("a", 0),
    scale = "link"
  )
  s <- summary(out)
  expect_equal(s$z_mean, 1)
  expect_identical(
    c(s$z_sd, s$z_mean_mcse, s$z_sd_mcse, s$contraction),
    rep(NA_real_, 4)
  )
})

test_that("subject rows pool z across replications, with no within step", {
  estimates <- moment_estimates(
    rep("a", 4), c(1, 2, -1, 0), 1,
    level = "subject", id = c("1", "2", "1", "2"), replication = c(1, 1, 2, 2)
  )
  truth <- fake_truth(rep("a", 4), 0,
    id = c("1", "2", "1", "2"), replication = c(1, 1, 2, 2)
  )
  s <- summary(recover_subjects(estimates, truth, scale = "link"))
  expect_equal(s$z_mean, 0.5)
  expect_equal(s$z_sd, sd(c(1, 2, -1, 0)))
  expect_identical(s$contraction, NA_real_)
})

test_that("a calibrated simulation has z_sd within its MCSE of 1", {
  # theta ~ N(0, tau^2); ybar ~ N(theta, sigma^2 / n): the conjugate
  # posterior is calibrated, so z ~ N(0, 1) over replications
  withr::local_seed(2210)
  reps <- 400
  tau <- 1.5
  se <- 0.5
  theta <- stats::rnorm(reps, 0, tau)
  ybar <- stats::rnorm(reps, theta, se)
  post_var <- 1 / (1 / tau^2 + 1 / se^2)
  post_mean <- post_var * ybar / se^2
  estimates <- moment_estimates("mu", post_mean, sqrt(post_var),
    replication = seq_len(reps)
  )
  truth <- fake_truth("mu", theta, replication = seq_len(reps))
  s <- summary(recover(estimates, truth,
    scale = "link", prior_sd = prior_table("mu", tau)
  ))
  expect_lt(abs(s$z_sd - 1), 3 * s$z_sd_mcse)
  expect_lt(abs(s$z_mean), 3 * s$z_mean_mcse)
  expect_equal(s$contraction, 1 - post_var / tau^2)
})

test_that("an empty summary has the new columns with their types", {
  empty <- empty_recovery_summary()
  for (column in c(
    "z_mean", "z_sd", "contraction", "z_mean_mcse",
    "z_sd_mcse"
  )) {
    expect_type(empty[[column]], "double")
  }
  expect_true(all(recovery_summary_columns() %in% names(empty)))
})

test_that("the print shows z_mean, z_sd, contraction after coverage_50", {
  out <- recover(moment_estimates("a", c(1, 2), 1, replication = 1:2),
    fake_truth("a", 0),
    scale = "link"
  )
  columns <- summary_print_columns(summary(out))
  at <- match("coverage_50", columns)
  expect_identical(columns[at + 1:3], c("z_mean", "z_sd", "contraction"))
  expect_false(any(c("z_mean_mcse", "z_sd_mcse") %in% columns))
})

# prior_sd from a prior-only fit or a prior check (D87) -----------------------

#' A stand-in brmsfit whose extractor gives the moments of stored draws
#' @noRd
fake_prior_fit <- function(draws, sample_prior = "only") {
  structure(
    list(
      draws = draws,
      prior = structure(data.frame(), sample_prior = sample_prior)
    ),
    class = c("zpriorfit", "brmsfit")
  )
}

registerS3method(
  "extract_estimates", "zpriorfit",
  function(fit, level = "population", ...) {
    rows <- lapply(names(fit$draws), function(term) {
      d <- fit$draws[[term]]
      moment_estimates(term, mean(d), stats::sd(d))
    })
    out <- dplyr::bind_rows(rows)
    out[out$level %in% level, ]
  },
  envir = asNamespace("bmmtools")
)

test_that("no data: a posterior equal to the prior has contraction near 0", {
  withr::local_seed(2211)
  prior <- fake_prior_fit(list(mu = stats::rnorm(4000, 0, 2)))
  posterior <- stats::rnorm(4000, 0, 2)
  estimates <- moment_estimates("mu", mean(posterior), stats::sd(posterior))
  out <- recover(estimates, fake_truth("mu", 0),
    scale = "link", prior_sd = prior
  )
  # the variance ratio of two samples of 4000 has an SE near 0.03
  expect_lt(abs(out$contraction), 0.15)
})

test_that("a fit that did not sample the prior alone is refused", {
  estimates <- moment_estimates("mu", 0, 1)
  expect_error(
    recover(estimates, fake_truth("mu", 0),
      scale = "link",
      prior_sd = fake_prior_fit(list(mu = c(-1, 1)), sample_prior = "no")
    ),
    "must have sampled the prior alone"
  )
})

fake_prior_check <- function(files) {
  structure(
    tibble::tibble(),
    class = c("bmmtools_prior_check", "tbl_df", "tbl", "data.frame"),
    sets = names(files), files = files
  )
}

test_that("a one-set prior check gives the prior SDs of its cached fit", {
  dir <- withr::local_tempdir()
  stem <- file.path(dir, "prior-default")
  saveRDS(fake_prior_fit(list(mu = c(-2, 2))), paste0(stem, ".rds"))
  estimates <- moment_estimates("mu", 0, 1)
  out <- recover(estimates, fake_truth("mu", 0),
    scale = "link", prior_sd = fake_prior_check(c(default = stem))
  )
  expect_equal(out$contraction, 1 - 1 / stats::var(c(-2, 2)))
})

test_that("a prior check with several sets or no file is refused", {
  estimates <- moment_estimates("mu", 0, 1)
  truth <- fake_truth("mu", 0)
  expect_error(
    recover(estimates, truth,
      scale = "link",
      prior_sd = fake_prior_check(c(default = "a", wide = "b"))
    ),
    "must hold one prior set"
  )
  expect_error(
    recover(estimates, truth,
      scale = "link",
      prior_sd = fake_prior_check(c(default = "/no/such/prior-fit"))
    ),
    "fit is not at"
  )
})

# prior_sd_table(): closed forms from a fit's priors (D88) -------------------

#' A prior_summary() of the shape the mixture2p fixture has (measured
#' 2026-10-05), with a cell-means parameter and a flat one added
fake_priors <- function() {
  data.frame(
    prior = c(
      "", "normal(2, 1)", "", "logistic(0, 1)", "constant(-100)",
      "student_t(3, 0, 2.5)", "student_t(3, 0, 2.5)", "", "", "", "",
      "logistic(0, 1)", "normal(0, 0.5)", "", "", ""
    ),
    class = c(
      "b", "b", "b", "b", "Intercept", "sd", "sd", "sd", "sd", "sd", "sd",
      "theta2", "b", "b", "b", "b"
    ),
    coef = c(
      "", "Intercept", "", "Intercept", "", "", "", "", "Intercept", "",
      "Intercept", "", "", "tasklow", "taskhigh", "Intercept"
    ),
    group = c(
      "", "", "", "", "", "", "", "id", "id", "id", "id", "", "", "",
      "", ""
    ),
    nlpar = c(
      "kappa", "kappa", "thetat", "thetat", "", "kappa", "thetat", "kappa",
      "kappa", "thetat", "thetat", "", "c", "c", "c", "flat"
    ),
    dpar = c("", "", "", "", "kappa2", rep("", 11)),
    lb = c("", "", "", "", "", "0", "0", "", "", "", "", "", "", "", "", ""),
    ub = "",
    source = "default"
  )
}

test_that("prior_sd_table() gives the closed-form prior SD of each term", {
  out <- suppressMessages(prior_sd_table(fake_priors()))
  expect_named(out, c("term", "level", "prior", "prior_sd_link"))
  get <- function(term, level) {
    out$prior_sd_link[out$term == term & out$level == level]
  }
  expect_equal(get("kappa", "population"), 1)
  expect_equal(get("thetat", "population"), pi / sqrt(3))
  # cell means inherit the class prior of their parameter
  expect_equal(get("c_tasklow", "population"), 0.5)
  expect_equal(get("c_taskhigh", "population"), 0.5)
  # half-t(3, 0, 2.5): checked against integrate() in spec 22.2.0
  expect_equal(get("kappa", "sd"), 3.339298, tolerance = 1e-6)
  expect_equal(get("thetat", "sd"), 3.339298, tolerance = 1e-6)
  expect_identical(get("flat", "population"), NA_real_)
  # fixed parameters and brms's own classes have no row
  expect_false(any(out$term %in% c("kappa2", "theta2", "")))
})

test_that("prior_sd_table() says which terms have no closed form", {
  expect_message(
    prior_sd_table(fake_priors()),
    "No closed-form prior SD for .*flat \\(population\\)"
  )
})

test_that("the closed forms match the distributions they describe", {
  table <- function(prior, class = "b", lb = "") {
    data.frame(
      prior = prior, class = class, coef = "Intercept",
      group = if (class == "sd") "id" else "", nlpar = "a", dpar = "",
      lb = lb, ub = ""
    )
  }
  sd_of <- function(...) {
    suppressMessages(prior_sd_table(table(...)))$prior_sd_link
  }
  expect_equal(sd_of("normal(0.5, 0.3)"), 0.3)
  expect_equal(sd_of("logistic(1, 2)"), 2 * pi / sqrt(3))
  expect_equal(sd_of("student_t(5, 0, 2)"), 2 * sqrt(5 / 3))
  expect_identical(sd_of("student_t(2, 0, 2)"), NA_real_)
  expect_identical(sd_of("cauchy(0, 1)"), NA_real_)
  # a truncated b prior is not the closed form's distribution
  expect_identical(sd_of("normal(0, 1)", lb = "0"), NA_real_)
  moment <- function(k) {
    integrate(function(x) x^k * 2 * dnorm(x, 0, 1.5), 0, Inf)$value
  }
  half_normal <- moment(2) - moment(1)^2
  expect_equal(sd_of("normal(0, 1.5)", class = "sd"), sqrt(half_normal),
    tolerance = 1e-6
  )
  # a half prior is only half of a distribution centred on 0
  expect_identical(sd_of("normal(1, 1.5)", class = "sd"), NA_real_)
  # bmm 1.3.2.9000's default on a standard deviation (measured 2026-10-05)
  expect_equal(sd_of("exponential(1)", class = "sd"), 1)
  expect_equal(sd_of("exponential(4)", class = "sd"), 0.25)
  # on an unbounded coefficient it is not the prior's distribution
  expect_identical(sd_of("exponential(1)"), NA_real_)
})

test_that("prior_sd_table() refuses what is not a prior table or a fit", {
  expect_error(prior_sd_table(1:3), "must be a fit or a prior table")
  expect_error(
    prior_sd_table(data.frame(prior = "x")),
    "lacks the columns .*class.*coef.*group.*nlpar"
  )
})

test_that("prior_sd_table() of the fixture fit reads its prior_summary()", {
  skip_if_not_installed("brms")
  skip_if_not_installed("bmm")
  out <- suppressMessages(prior_sd_table(mixture2p_fit()))
  estimates <- extract_estimates(mixture2p_fit(), level = c("population", "sd"))
  key <- paste(out$term, out$level)
  expect_true(all(paste(estimates$term, estimates$level) %in% key))
  expect_equal(
    out$prior_sd_link[match(c("kappa population", "thetat sd"), key)],
    c(1, 3.339298),
    tolerance = 1e-6
  )
})

test_that("the table it returns is a prior_sd for recover()", {
  estimates <- moment_estimates("kappa", 2, 0.25)
  out <- suppressMessages(recover(estimates, fake_truth("kappa", 2),
    scale = "link", prior_sd = prior_sd_table(fake_priors())
  ))
  expect_equal(out$contraction, 1 - 0.25^2)
})

# recovery_grid(prior_sd = ) (D89) ---------------------------------------------

#' The grid mock fitter, but a prior-only call returns a prior fit with
#' known SDs and is counted apart
prior_grid_fitter <- function() {
  mock <- grid_mock_fitter()
  prior_calls <- new.env(parent = emptyenv())
  prior_calls$n <- 0L
  prior_calls$backend <- list()
  fitter <- function(formula, data, model, prior = NULL, ...) {
    dots <- list(...)
    if (identical(dots$sample_prior, "only")) {
      prior_calls$n <- prior_calls$n + 1L
      prior_calls$backend[[prior_calls$n]] <- dots$backend
      return(fake_prior_fit(list(
        kappa = c(-10, 10), thetat = c(-20, 20)
      )))
    }
    mock$fitter(formula, data, model, prior = prior, ...)
  }
  list(fitter = fitter, calls = mock$calls, prior_calls = prior_calls)
}

zgrid_run <- function(dir, fitter, prior_sd = NULL, grid = small_grid(), ...) {
  suppressMessages(recovery_grid(
    bmm::mixture2p(resp_error = "y"),
    grid = grid,
    pars = c(kappa = log(8), thetat = stats::qlogis(0.75)),
    dir = dir,
    reps = 2L,
    sds = c(kappa = 0.3, thetat = 0.4),
    seed = 100,
    prior_sd = prior_sd,
    ...,
    .fitter = fitter
  ))
}

test_that("a grid without prior_sd has z, no contraction and no extra fit", {
  skip_if_not_installed("bmm")
  mock <- prior_grid_fitter()
  out <- zgrid_run(withr::local_tempdir(), mock$fitter)
  expect_identical(mock$prior_calls$n, 0L)
  # four cells and the preflight
  expect_identical(mock$calls$n, 5L)
  # the mock's posterior: mean 0, SD 5 on the link scale
  population <- out[out$level == "population", ]
  truth_link <- c(kappa = log(8), thetat = stats::qlogis(0.75))
  expect_equal(population$z, unname(-truth_link[population$term] / 5))
  expect_false(anyNA(out$z))
  expect_true(all(is.na(out$contraction)))
})

test_that("a grid takes a prior_sd table", {
  skip_if_not_installed("bmm")
  mock <- prior_grid_fitter()
  out <- zgrid_run(withr::local_tempdir(), mock$fitter,
    prior_sd = prior_table(c("kappa", "thetat"), c(10, 25))
  )
  population <- out[out$level == "population", ]
  expected <- 1 - 25 / c(kappa = 100, thetat = 625)[population$term]
  expect_equal(population$contraction, unname(expected))
  expect_true(all(is.na(out$contraction[out$level == "subject"])))
})

test_that("prior_sd = 'fit' fits once per set, cached, with a backend", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- prior_grid_fitter()
  out <- zgrid_run(dir, mock$fitter, prior_sd = "fit")
  # two grid rows share the formula, model and prior: one prior fit
  expect_identical(mock$prior_calls$n, 1L)
  expect_false(is.null(mock$prior_calls$backend[[1L]]))
  population <- out[out$level == "population", ]
  expected <- 1 - 25 / c(kappa = 200, thetat = 800)[population$term]
  expect_equal(population$contraction, unname(expected))
  # a second run reuses it
  again <- zgrid_run(dir, mock$fitter, prior_sd = "fit")
  expect_identical(mock$prior_calls$n, 1L)
  expect_equal(again$contraction, out$contraction)
})

test_that("prior_sd = 'analytic' reads bmm's default priors, no extra fit", {
  skip_if_not_installed("bmm")
  mock <- prior_grid_fitter()
  out <- zgrid_run(withr::local_tempdir(), mock$fitter, prior_sd = "analytic")
  expect_identical(mock$prior_calls$n, 0L)
  population <- out[out$level == "population", ]
  # the measured defaults: normal(2, 1) on kappa, logistic(0, 1) on thetat
  prior <- c(kappa = 1, thetat = pi / sqrt(3))[population$term]
  expect_equal(population$contraction, unname(1 - 25 / prior^2))
})

test_that("recovery_grid() refuses a prior_sd it cannot use", {
  skip_if_not_installed("bmm")
  mock <- prior_grid_fitter()
  expect_error(
    zgrid_run(withr::local_tempdir(), mock$fitter, prior_sd = "both"),
    "must be .*NULL.*fit.*analytic"
  )
  expect_identical(mock$calls$n, 0L)
})

# collect_grid(prior_sd = ) ----------------------------------------------------

test_that("collect_grid() gives the contraction the grid gave", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- prior_grid_fitter()
  run <- zgrid_run(dir, mock$fitter, prior_sd = "fit")
  collected <- suppressMessages(collect_grid(dir))
  expect_equal(collected$contraction, run$contraction)
  expect_false(all(is.na(collected$contraction)))
  expect_identical(mock$prior_calls$n, 1L)
})

test_that("collect_grid() takes a prior_sd that overrides the grid's", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  zgrid_run(dir, prior_grid_fitter()$fitter, prior_sd = "fit")
  collected <- suppressMessages(collect_grid(dir,
    prior_sd = prior_table(c("kappa", "thetat"), c(10, 10))
  ))
  population <- collected[collected$level == "population", ]
  expect_equal(population$contraction, rep(0.75, nrow(population)))
})

test_that("a grid rerun without prior_sd leaves collect_grid() none", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- prior_grid_fitter()
  zgrid_run(dir, mock$fitter, prior_sd = prior_table("kappa", 10))
  zgrid_run(dir, mock$fitter)
  collected <- suppressMessages(collect_grid(dir))
  expect_true(all(is.na(collected$contraction)))
})

# review fixes (22.2) ----------------------------------------------------------

test_that("prior_sd = 'analytic' with a user prior is refused before fitting", {
  # bmm::default_prior(prior = ) replaces bmm's own priors instead of
  # merging them as bmm() does (measured 2026-10-05, bmm 1.3.2.9000)
  skip_if_not_installed("bmm")
  mock <- prior_grid_fitter()
  expect_error(
    zgrid_run(withr::local_tempdir(), mock$fitter,
      prior_sd = "analytic",
      prior = brms::set_prior("normal(0, 5)", class = "b", nlpar = "kappa")
    ),
    "analytic.*prior"
  )
  expect_identical(mock$calls$n, 0L)
})

test_that("prior_sd = 'fit' fits before the cells and not for subjects only", {
  skip_if_not_installed("bmm")
  log <- new.env(parent = emptyenv())
  log$order <- character()
  mock <- prior_grid_fitter()
  fitter <- function(formula, data, model, prior = NULL, ...) {
    only <- identical(list(...)$sample_prior, "only")
    log$order <- c(log$order, if (only) "prior" else "cell")
    mock$fitter(formula, data, model, prior = prior, ...)
  }
  zgrid_run(withr::local_tempdir(), fitter, prior_sd = "fit")
  order <- log$order
  expect_identical(order[[1L]] == "prior" || order[[2L]] == "prior", TRUE)
  expect_false("prior" %in% order[-(1:2)])
  subjects_only <- prior_grid_fitter()
  zgrid_run(withr::local_tempdir(), subjects_only$fitter,
    prior_sd = "fit", levels = "subject"
  )
  expect_identical(subjects_only$prior_calls$n, 0L)
})

test_that("a prior_sd table with a key twice is refused", {
  estimates <- moment_estimates("a", 0, 1)
  expect_error(
    recover(estimates, fake_truth("a", 0),
      scale = "link", prior_sd = prior_table(c("a", "a"), c(1, 2))
    ),
    "more than once"
  )
})

test_that("a Student t without a finite nu has no closed form", {
  table <- data.frame(
    prior = "student_t(Inf, 0, 1)", class = "b", coef = "Intercept",
    group = "", nlpar = "a"
  )
  expect_identical(
    suppressMessages(prior_sd_table(table))$prior_sd_link, NA_real_
  )
})
