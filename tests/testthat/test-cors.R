# Tests for correlated truths: check_cors(), draw_correlated_pars() and
# cors_from_factors(), written against
# local/dev/spec-milestone-5-correlation-recovery.md section 5.1. Pure R;
# nothing here needs bmm.

named_cors <- function(r, names = c("kappa", "thetat")) {
  m <- matrix(c(1, r, r, 1), 2, 2)
  dimnames(m) <- list(names, names)
  m
}

toy_pars <- c(kappa = 2, thetat = 1)

# check_cors -------------------------------------------------------------

test_that("NULL cors is the identity over varying terms then covariates", {
  full <- check_cors(
    NULL, toy_pars, c(kappa = 0.3, thetat = 0.5),
    covariates = list(G = c(mean = 0, sd = 1))
  )
  expect_equal(dimnames(full), list(
    c("kappa", "thetat", "G"), c("kappa", "thetat", "G")
  ))
  expect_equal(unname(full), diag(3))
})

test_that("a partial matrix is embedded and ordered like the draw", {
  given <- named_cors(0.4, c("G", "thetat"))
  full <- check_cors(
    given, toy_pars, c(kappa = 0.3, thetat = 0.5),
    covariates = list(G = c(mean = 0, sd = 1))
  )
  expect_equal(rownames(full), c("kappa", "thetat", "G"))
  expect_equal(full["thetat", "G"], 0.4)
  expect_equal(full["G", "thetat"], 0.4)
  expect_equal(full["kappa", "G"], 0)
})

test_that("a non-varying parameter may appear only with zero correlations", {
  sds <- c(kappa = 0.3, thetat = 0)
  ok <- check_cors(named_cors(0), toy_pars, sds, covariates = NULL)
  expect_equal(dimnames(ok), list("kappa", "kappa"))
  expect_error(
    check_cors(named_cors(0.5), toy_pars, sds, covariates = NULL),
    "thetat"
  )
})

test_that("check_cors refuses malformed matrices, naming the problem", {
  sds <- c(kappa = 0.3, thetat = 0.5)
  check <- function(x) check_cors(x, toy_pars, sds, covariates = NULL)

  expect_error(check("a"), "matrix")
  expect_error(check(matrix(1, 2, 3)), "square")
  expect_error(check(diag(2)), "dimnames")
  not_sym <- named_cors(0.5)
  not_sym[1, 2] <- 0.2
  expect_error(check(not_sym), "symmetric")
  bad_diag <- named_cors(0.5)
  diag(bad_diag) <- 2
  expect_error(check(bad_diag), "diagonal")
  expect_error(check(named_cors(1.5)), "between")
  expect_error(check(named_cors(0.5, c("kappa", "zeta"))), "zeta")
  with_na <- named_cors(NA)
  expect_error(check(with_na), "missing")
  expect_error(check(named_cors(0.2, c("kappa", "kappa"))), "twice")
  # nothing varies: a matrix of zero correlations reduces to 0 by 0
  none <- check_cors(
    named_cors(0), toy_pars, c(kappa = 0, thetat = 0),
    covariates = NULL
  )
  expect_equal(dim(none), c(0L, 0L))

  not_pd <- diag(3)
  dimnames(not_pd) <- rep(list(c("kappa", "thetat", "G")), 2)
  not_pd[1, 2] <- not_pd[2, 1] <- 0.9
  not_pd[1, 3] <- not_pd[3, 1] <- 0.9
  not_pd[2, 3] <- not_pd[3, 2] <- -0.9
  expect_error(
    check_cors(
      not_pd, toy_pars, sds,
      covariates = list(G = c(mean = 0, sd = 1))
    ),
    "positive definite"
  )
})

# draw_correlated_pars ---------------------------------------------------

test_that("the identity draw reproduces the Milestone 3 draw bit for bit", {
  pars <- c(a = 1, b = -1, c = 0.5, d = 2)
  sds <- c(a = 0.3, b = 0, c = 0.7, d = 1.1)
  n <- 25L
  # the Milestone 3 formula, one rnorm() per parameter in pars order
  oracle <- withr::with_seed(42, vapply(names(pars), function(p) {
    pars[[p]] + stats::rnorm(n, 0, sds[[p]])
  }, numeric(n)))
  full <- check_cors(NULL, pars, sds, covariates = NULL)
  drawn <- withr::with_seed(
    42, draw_correlated_pars(pars, sds, full, NULL, n)
  )
  expect_identical(unname(drawn[, names(pars)]), unname(oracle))
  expect_equal(rownames(drawn), as.character(seq_len(n)))
})

test_that("appended covariates never change the parameter columns", {
  sds <- c(kappa = 0.3, thetat = 0.5)
  covariates <- list(G = c(mean = 10, sd = 2))
  without <- withr::with_seed(
    7, draw_correlated_pars(
      toy_pars, sds, check_cors(named_cors(0.5), toy_pars, sds, NULL),
      NULL, 50L
    )
  )
  given <- diag(3)
  dimnames(given) <- rep(list(c("kappa", "thetat", "G")), 2)
  given["kappa", "thetat"] <- given["thetat", "kappa"] <- 0.5
  given["kappa", "G"] <- given["G", "kappa"] <- 0.6
  full <- check_cors(given, toy_pars, sds, covariates)
  with_cov <- withr::with_seed(
    7, draw_correlated_pars(toy_pars, sds, full, covariates, 50L)
  )
  expect_identical(
    with_cov[, c("kappa", "thetat")], without[, c("kappa", "thetat")]
  )
  expect_true("G" %in% colnames(with_cov))
})

test_that("large samples hit the target correlations, means and SDs", {
  sds <- c(kappa = 0.3, thetat = 0.5)
  covariates <- list(G = c(mean = 10, sd = 2))
  given <- diag(3)
  dimnames(given) <- rep(list(c("kappa", "thetat", "G")), 2)
  given["kappa", "thetat"] <- given["thetat", "kappa"] <- 0.5
  given["kappa", "G"] <- given["G", "kappa"] <- 0.3
  full <- check_cors(given, toy_pars, sds, covariates)
  drawn <- withr::with_seed(
    3, draw_correlated_pars(toy_pars, sds, full, covariates, 4000L)
  )
  r <- stats::cor(drawn)
  expect_equal(r["kappa", "thetat"], 0.5, tolerance = 0.05 / 0.5)
  expect_lt(abs(r["kappa", "G"] - 0.3), 0.05)
  expect_lt(abs(r["thetat", "G"]), 0.05)
  expect_lt(abs(mean(drawn[, "G"]) - 10), 0.1)
  expect_lt(abs(stats::sd(drawn[, "thetat"]) - 0.5), 0.02)
})

test_that("non-varying parameters stay at their population value", {
  sds <- c(kappa = 0.3, thetat = 0)
  full <- check_cors(NULL, toy_pars, sds, NULL)
  drawn <- withr::with_seed(
    1, draw_correlated_pars(toy_pars, sds, full, NULL, 5L)
  )
  expect_equal(unname(drawn[, "thetat"]), rep(1, 5))
})

# cors_from_factors ------------------------------------------------------

test_that("cors_from_factors matches a tcrossprod oracle", {
  loadings <- matrix(
    c(0.8, 0.7, 0.1, 0.2, 0.0, 0.3, 0.6, 0.75),
    ncol = 2,
    dimnames = list(c("a", "b", "c", "d"), c("F1", "F2"))
  )
  phi <- matrix(c(1, 0.4, 0.4, 1), 2, dimnames = list(
    c("F1", "F2"), c("F1", "F2")
  ))
  # phi = t(U) U, so L t(U) (L t(U))' = L phi L'
  oracle <- tcrossprod(loadings %*% t(chol(phi)))
  diag(oracle) <- 1

  out <- cors_from_factors(loadings, phi)
  expect_equal(out, oracle)
  expect_equal(dimnames(out), list(rownames(loadings), rownames(loadings)))

  single <- cors_from_factors(c(a = 0.8, b = 0.6, c = 0.5))
  expect_equal(single["a", "b"], 0.48)
  expect_equal(unname(diag(single)), rep(1, 3))
})

test_that("cors_from_factors refuses impossible loadings", {
  expect_error(cors_from_factors(c(a = 1.2, b = 0.5)), "communalit")
  expect_error(cors_from_factors(c(0.5, 0.5)), "names")
  expect_error(cors_from_factors("x"), "numeric")
  expect_error(cors_from_factors(c(a = NA, b = 0.5)), "missing")
  expect_error(
    cors_from_factors(c(a = 0.5, b = 0.5), diag(2)),
    "1 by 1"
  )
  loadings <- matrix(c(0.6, 0.6, 0.5, 0.5), 2, dimnames = list(
    c("a", "b"), c("F1", "F2")
  ))
  bad_phi <- matrix(c(1, 2, 2, 1), 2)
  expect_error(cors_from_factors(loadings, bad_phi), "factor_cors")
  wrong_names <- diag(2)
  dimnames(wrong_names) <- list(c("F1", "F3"), c("F1", "F3"))
  expect_error(cors_from_factors(loadings, wrong_names), "F3")
  # a communality of exactly 1 leaves a singular matrix
  expect_error(
    cors_from_factors(c(a = 1, b = 1, c = 0.5)),
    "positive definite"
  )
})
