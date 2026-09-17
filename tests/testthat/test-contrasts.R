# Tests for the contrast design: check_contrasts() and contrast_truth().

# the cell-means truth the transform is given, for one parameter over two
# tasks: population means 2 and 3, SD 0.3 each, correlated 0.5
cell_case <- function(rho = 0.5, sd = 0.3, tasks = c("1", "2")) {
  terms <- paste0("kappa_task", tasks)
  pars <- stats::setNames(c(2, 3), terms)
  sds <- stats::setNames(rep(sd, 2L), terms)
  cors <- matrix(c(1, rho, rho, 1), 2L, dimnames = list(terms, terms))
  values <- matrix(
    c(2, 2.5, 3, 3.5),
    nrow = 2L, dimnames = list(c("1", "2"), terms)
  )
  list(pars = pars, sds = sds, cors = cors, values = values, tasks = tasks)
}

test_that("check_contrasts takes a matrix or a function and refuses the rest", {
  expect_equal(check_contrasts(NULL, 2L), stats::contr.treatment(2))
  expect_equal(
    check_contrasts(stats::contr.sum, 3L), stats::contr.sum(3)
  )
  given <- matrix(c(-0.5, 0.5), nrow = 2L)
  expect_equal(check_contrasts(given, 2L), given)

  expect_error(check_contrasts(matrix(1, 3L, 1L), 2L), "row")
  expect_error(check_contrasts(matrix(1, 2L, 2L), 2L), "column")
  expect_error(check_contrasts("treatment", 2L), "numeric matrix")
  expect_error(check_contrasts(matrix(c(1, NA), 2L, 1L), 2L), "NA")
  # a contrast column that is the intercept again: no inverse
  expect_error(check_contrasts(matrix(c(1, 1), 2L, 1L), 2L), "invertible")
  expect_error(
    check_contrasts(function(n) stop("no"), 2L), "failed when called"
  )
})

test_that("contrast_truth transforms treatment coding by hand", {
  case <- cell_case()
  out <- contrast_truth(
    case$pars, case$sds, case$values, case$cors,
    contrasts = stats::contr.treatment(2), tasks = case$tasks,
    task_col = "task"
  )
  expect_equal(names(out$pars), c("kappa", "kappa_task1"))
  # intercept is the first cell, slope the difference
  expect_equal(unname(out$pars), c(2, 1))
  expect_equal(unname(out$values[, "kappa"]), c(2, 2.5))
  expect_equal(unname(out$values[, "kappa_task1"]), c(1, 1))

  # sd .3 / .3 with rho .5 gives sd .3 / .3 and a correlation of -.5
  expect_equal(unname(out$sds), c(0.3, 0.3))
  expect_equal(out$cors["kappa", "kappa_task1"], -0.5)
})

test_that("contrast_truth transforms equalprior coding by hand", {
  skip_if_not_installed("bayestestR")
  case <- cell_case()
  out <- contrast_truth(
    case$pars, case$sds, case$values, case$cors,
    contrasts = bayestestR::contr.equalprior(2), tasks = case$tasks,
    task_col = "task"
  )
  # the intercept is the grand mean, the coefficient the difference over
  # sqrt(2), because contr.equalprior(2) codes the levels at -+0.7071
  expect_equal(unname(out$pars), c(2.5, 1 / sqrt(2)))

  # an orthogonal contrast decorrelates the intercept from the slope
  expect_equal(out$cors["kappa", "kappa_task1"], 0)
  # measured 2026-09-17: sd * sqrt((1 + rho) / 2) and sd * sqrt(1 - rho).
  # The spec's table gave the slope as sd * sqrt((1 - rho) / 2), which is
  # the value for no contrast matrix bmm uses: contr.equalprior_pairs(2)
  # codes -+0.5, whose coefficient is the plain difference (sd 0.3 here).
  expect_equal(unname(out$sds), c(0.3 * sqrt(0.75), 0.3 * sqrt(0.5)))

  pairs <- contrast_truth(
    case$pars, case$sds, case$values, case$cors,
    contrasts = bayestestR::contr.equalprior_pairs(2), tasks = case$tasks,
    task_col = "task"
  )
  expect_equal(unname(pairs$pars), c(2.5, 1))
  expect_equal(unname(pairs$sds), c(0.3 * sqrt(0.75), 0.3))
})

test_that("contrast_truth matches the covariance of simulated values", {
  # the check that catches a transposed map: draw cell values from the
  # stated covariance, transform them by hand, and compare the empirical
  # covariance with the transform's SD and correlation truth
  skip_if_not_installed("bayestestR")
  case <- cell_case(rho = 0.3)
  contrasts <- bayestestR::contr.equalprior(2)
  n <- 20000L
  sigma <- outer(case$sds, case$sds) * case$cors
  draws <- withr::with_seed(
    42, matrix(stats::rnorm(n * 2L), ncol = 2L) %*% chol(sigma)
  )
  colnames(draws) <- names(case$pars)
  drawn <- draws + rep(case$pars, each = n)
  rownames(drawn) <- as.character(seq_len(n))

  out <- contrast_truth(
    case$pars, case$sds, drawn, case$cors,
    contrasts = contrasts, tasks = case$tasks, task_col = "task"
  )
  empirical <- stats::cov(out$values)
  expect_equal(sqrt(diag(empirical)), out$sds, tolerance = 0.03)
  expect_equal(
    stats::cov2cor(empirical)[1L, 2L],
    out$cors["kappa", "kappa_task1"],
    tolerance = 0.03
  )
  # and the population values are the mean of the transformed draws
  expect_equal(colMeans(out$values), out$pars, tolerance = 0.02)
})

test_that("contrast_truth passes covariates and non-varying cells through", {
  case <- cell_case()
  case$sds[["kappa_task2"]] <- 0
  cors <- matrix(
    c(1, 0.4, 0.4, 1), 2L,
    dimnames = list(
      c("kappa_task1", "ability"), c("kappa_task1", "ability")
    )
  )
  values <- cbind(case$values, ability = c(-1, 1))
  out <- contrast_truth(
    case$pars, case$sds, values, cors,
    contrasts = stats::contr.treatment(2), tasks = case$tasks,
    task_col = "task", covariates = list(ability = NULL)
  )
  expect_equal(unname(out$values[, "ability"]), c(-1, 1))
  expect_true("ability" %in% rownames(out$cors))
  # one cell varies, so both the intercept and the slope inherit its SD
  expect_equal(unname(out$sds), c(0.3, 0.3))
  expect_equal(out$cors["kappa", "kappa_task1"], -1)
  # the covariate keeps its correlation with the first cell, which is the
  # intercept under treatment coding
  expect_equal(out$cors["ability", "kappa"], 0.4)
})

test_that("contrast_truth handles several parameters and three tasks", {
  tasks <- c("a", "b", "c")
  terms <- c(
    paste0("kappa_task", tasks), paste0("thetat_task", tasks)
  )
  pars <- stats::setNames(c(2, 3, 4, -1, 0, 1), terms)
  sds <- stats::setNames(rep(0.2, 6L), terms)
  cors <- diag(6L)
  dimnames(cors) <- list(terms, terms)
  values <- matrix(
    rep(pars, each = 2L),
    nrow = 2L, dimnames = list(c("1", "2"), terms)
  )
  out <- contrast_truth(
    pars, sds, values, cors,
    contrasts = stats::contr.treatment(3), tasks = tasks, task_col = "task"
  )
  expect_equal(
    names(out$pars),
    c(
      "kappa", "kappa_task1", "kappa_task2",
      "thetat", "thetat_task1", "thetat_task2"
    )
  )
  expect_equal(unname(out$pars), c(2, 1, 2, -1, 1, 2))
  # the two parameters are transformed independently: no cross term
  expect_equal(out$cors["kappa", "thetat"], 0)
})
