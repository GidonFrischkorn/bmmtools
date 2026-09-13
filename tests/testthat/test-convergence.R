# Tests for check_convergence(), written against
# local/dev/spec-milestone-2-run-layer.md section 1.
#
# The gate's logic runs on two tables --- the summarise_draws() output
# and the nuts_params() output --- so it is tested on hand-built tables
# first. The fixture tests exercise the brmsfit method and skip when
# brms is absent. No Stan is compiled.

fake_diagnostics <- function(variable, rhat, ess_bulk, ess_tail) {
  tibble::tibble(
    variable = variable,
    rhat = as.double(rhat),
    ess_bulk = as.double(ess_bulk),
    ess_tail = as.double(ess_tail)
  )
}

fake_nuts <- function(divergent, treedepth) {
  n <- max(length(divergent), length(treedepth))
  divergent <- rep_len(divergent, n)
  treedepth <- rep_len(treedepth, n)
  rbind(
    data.frame(
      Chain = 1L, Iteration = seq_len(n), Parameter = "divergent__",
      Value = as.double(divergent)
    ),
    data.frame(
      Chain = 1L, Iteration = seq_len(n), Parameter = "treedepth__",
      Value = as.double(treedepth)
    )
  )
}

default_thresholds <- function(...) {
  utils::modifyList(
    list(
      rhat_max = 1.05, ess_bulk_min = 400, ess_tail_min = NULL,
      divergent_max = 10
    ),
    list(...)
  )
}

# hand-built tables ------------------------------------------------------

test_that("a fit meeting every criterion passes with the right numbers", {
  diagnostics <- fake_diagnostics(
    c("b_a", "b_b", "sd_id__a"),
    rhat = c(1.001, 1.02, 1.01),
    ess_bulk = c(900, 450, 600),
    ess_tail = c(800, 500, 410)
  )
  nuts <- fake_nuts(
    divergent = c(0, 1, 0, 0, 0),
    treedepth = c(3, 5, 10, 4, 10)
  )

  out <- convergence_from_parts(
    diagnostics, nuts,
    treedepth_max = 10, thresholds = default_thresholds()
  )

  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 1L)
  expect_named(out, c(
    "max_rhat", "min_ess_bulk", "min_ess_tail", "n_divergent",
    "n_max_treedepth", "n_variables", "pass", "failed"
  ))
  expect_equal(out$max_rhat, 1.02)
  expect_equal(out$min_ess_bulk, 450)
  expect_equal(out$min_ess_tail, 410)
  expect_identical(out$n_divergent, 1L)
  expect_identical(out$n_max_treedepth, 2L)
  expect_identical(out$n_variables, 3L)
  expect_true(out$pass)
  expect_identical(out$failed, "")
})

test_that("rows with NA rhat are dropped, not failed on", {
  # constants have NA diagnostics (measured on the fixture); they must
  # neither fail the gate nor count as assessed variables
  diagnostics <- fake_diagnostics(
    c("b_mu1", "b_a"),
    rhat = c(NA, 1.01), ess_bulk = c(NA, 500), ess_tail = c(NA, 500)
  )
  out <- convergence_from_parts(
    diagnostics, NULL,
    treedepth_max = 10, thresholds = default_thresholds()
  )
  expect_true(out$pass)
  expect_identical(out$n_variables, 1L)
  expect_equal(out$max_rhat, 1.01)
})

test_that("rhat above the threshold fails and is named", {
  diagnostics <- fake_diagnostics(
    c("b_a", "b_b"),
    rhat = c(1.01, 1.08), ess_bulk = c(500, 500), ess_tail = c(500, 500)
  )
  out <- convergence_from_parts(
    diagnostics, NULL,
    treedepth_max = 10, thresholds = default_thresholds()
  )
  expect_false(out$pass)
  expect_identical(out$failed, "rhat")
})

test_that("failed criteria are listed in a fixed order", {
  diagnostics <- fake_diagnostics(
    "b_a",
    rhat = 1.2, ess_bulk = 50, ess_tail = 20
  )
  out <- convergence_from_parts(
    diagnostics, NULL,
    treedepth_max = 10, thresholds = default_thresholds()
  )
  expect_identical(out$failed, "rhat, ess_bulk")

  # tail ESS enters only when a threshold is given
  gated <- convergence_from_parts(
    diagnostics, NULL,
    treedepth_max = 10, thresholds = default_thresholds(ess_tail_min = 100)
  )
  expect_identical(gated$failed, "rhat, ess_bulk, ess_tail")
})

test_that("tail ESS is gated only when a threshold is supplied", {
  diagnostics <- fake_diagnostics(
    "b_a",
    rhat = 1.0, ess_bulk = 1000, ess_tail = 50
  )
  free <- convergence_from_parts(
    diagnostics, NULL,
    treedepth_max = 10, thresholds = default_thresholds()
  )
  expect_true(free$pass)
  expect_equal(free$min_ess_tail, 50)

  gated <- convergence_from_parts(
    diagnostics, NULL,
    treedepth_max = 10, thresholds = default_thresholds(ess_tail_min = 100)
  )
  expect_false(gated$pass)
  expect_identical(gated$failed, "ess_tail")
})

test_that("divergences are gated at the threshold, inclusive", {
  diagnostics <- fake_diagnostics(
    "b_a",
    rhat = 1.0, ess_bulk = 1000, ess_tail = 1000
  )
  ten <- fake_nuts(divergent = c(rep(1, 10), rep(0, 5)), treedepth = 3)
  eleven <- fake_nuts(divergent = c(rep(1, 11), rep(0, 4)), treedepth = 3)

  out_ten <- convergence_from_parts(
    diagnostics, ten,
    treedepth_max = 10, thresholds = default_thresholds()
  )
  out_eleven <- convergence_from_parts(
    diagnostics, eleven,
    treedepth_max = 10, thresholds = default_thresholds()
  )
  expect_true(out_ten$pass)
  expect_identical(out_ten$n_divergent, 10L)
  expect_false(out_eleven$pass)
  expect_identical(out_eleven$failed, "divergent")
})

test_that("a missing NUTS table leaves the sampler counts NA, not failed", {
  diagnostics <- fake_diagnostics(
    "b_a",
    rhat = 1.0, ess_bulk = 1000, ess_tail = 1000
  )
  out <- convergence_from_parts(
    diagnostics, NULL,
    treedepth_max = 10, thresholds = default_thresholds()
  )
  expect_identical(out$n_divergent, NA_integer_)
  expect_identical(out$n_max_treedepth, NA_integer_)
  expect_true(out$pass)
})

test_that("tree-depth hits are counted but not gated", {
  diagnostics <- fake_diagnostics(
    "b_a",
    rhat = 1.0, ess_bulk = 1000, ess_tail = 1000
  )
  nuts <- fake_nuts(divergent = 0, treedepth = c(10, 10, 10, 9, 11))
  out <- convergence_from_parts(
    diagnostics, nuts,
    treedepth_max = 10, thresholds = default_thresholds()
  )
  expect_identical(out$n_max_treedepth, 4L)
  expect_true(out$pass)
  expect_identical(out$failed, "")
})

test_that("nothing assessable gives an NA verdict, not a pass", {
  diagnostics <- fake_diagnostics(
    c("b_mu1", "b_mu2"),
    rhat = c(NA, NA), ess_bulk = c(NA, NA), ess_tail = c(NA, NA)
  )
  out <- convergence_from_parts(
    diagnostics, NULL,
    treedepth_max = 10, thresholds = default_thresholds()
  )
  expect_identical(out$pass, NA)
  expect_identical(out$failed, NA_character_)
  expect_identical(out$n_variables, 0L)
  expect_true(is.na(out$max_rhat))
})

test_that("the thresholds used are carried as an attribute", {
  diagnostics <- fake_diagnostics(
    "b_a",
    rhat = 1.0, ess_bulk = 1000, ess_tail = 1000
  )
  thresholds <- default_thresholds(rhat_max = 1.01, divergent_max = 0)
  out <- convergence_from_parts(
    diagnostics, NULL,
    treedepth_max = 12, thresholds = thresholds
  )
  used <- attr(out, "thresholds")
  expect_equal(used$rhat_max, 1.01)
  expect_equal(used$divergent_max, 0)
  expect_equal(used$treedepth_max, 12)
})

# argument checks -------------------------------------------------------

test_that("thresholds must be single numbers", {
  expect_error(check_thresholds(rhat_max = "a"), "rhat_max")
  expect_error(check_thresholds(ess_bulk_min = c(1, 2)), "ess_bulk_min")
  expect_error(check_thresholds(divergent_max = NA_real_), "divergent_max")
  expect_no_error(check_thresholds(ess_tail_min = NULL))
})

test_that("the default method errors on a non-fit", {
  expect_error(check_convergence(1:3), "brmsfit")
  expect_error(check_convergence(list()), "brmsfit")
})

# the fixture -----------------------------------------------------------

test_that("check_convergence() on the fixture agrees with posterior and brms", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()

  out <- check_convergence(fit)

  expect_equal(nrow(out), 1L)
  expect_named(out, c(
    "max_rhat", "min_ess_bulk", "min_ess_tail", "n_divergent",
    "n_max_treedepth", "n_variables", "pass", "failed"
  ))

  # the oracle, computed here rather than remembered
  conv <- posterior::summarise_draws(
    posterior::as_draws_array(fit), posterior::default_convergence_measures()
  )
  expect_equal(out$max_rhat, max(conv$rhat, na.rm = TRUE))
  expect_equal(out$min_ess_bulk, min(conv$ess_bulk, na.rm = TRUE))
  expect_identical(out$n_variables, sum(!is.na(conv$rhat)))

  np <- brms::nuts_params(fit)
  expect_identical(
    out$n_divergent,
    as.integer(sum(np$Value[np$Parameter == "divergent__"]))
  )
  expect_type(out$pass, "logical")
})

test_that("treedepth_max is read from the fit when NULL", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  expect_equal(fit_treedepth_max(fit), 10)
  expect_equal(attr(check_convergence(fit), "thresholds")$treedepth_max, 10)
  expect_equal(
    attr(check_convergence(fit, treedepth_max = 6), "thresholds")$treedepth_max,
    6
  )
})

test_that("variables restricts what the gate looks at", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()

  one <- check_convergence(fit, variables = "b_kappa_Intercept")
  expect_identical(one$n_variables, 1L)

  expect_error(
    check_convergence(fit, variables = c("b_kappa_Intercept", "nope")),
    "nope"
  )
})

# the branches a fixture cannot reach ----------------------------------

test_that("treedepth_max falls back to 10 when the fit carries none", {
  expect_equal(fit_treedepth_max(list()), 10)
  expect_equal(fit_treedepth_max(list(fit = NULL)), 10)
})

test_that("a NUTS table without a diagnostic leaves that count NA", {
  nuts <- fake_nuts(divergent = c(0, 1), treedepth = c(3, 4))
  no_treedepth <- nuts[nuts$Parameter != "treedepth__", ]
  diagnostics <- fake_diagnostics(
    "b_a",
    rhat = 1.0, ess_bulk = 1000, ess_tail = 1000
  )
  out <- convergence_from_parts(
    diagnostics, no_treedepth,
    treedepth_max = 10, thresholds = default_thresholds()
  )
  expect_identical(out$n_divergent, 1L)
  expect_identical(out$n_max_treedepth, NA_integer_)
})

test_that("treedepth_max and variables are validated", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  expect_error(check_convergence(fit, treedepth_max = "a"), "treedepth_max")
  expect_error(check_convergence(fit, treedepth_max = c(1, 2)), "treedepth_max")
  expect_error(check_convergence(fit, variables = 1), "variables")
})
