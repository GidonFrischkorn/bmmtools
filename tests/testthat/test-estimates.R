# Tests for extract_estimates(), written against
# dev/spec-milestone-1-score-layer.md section 2.
#
# Two kinds of test live here. Those that need the *structure* of a real
# bmm fit read the saved fixture and skip when brms is absent; those that
# need a structure the fixture does not have (two grouping factors, a
# chain with a missing draw) build a draws object by hand and need no
# suggested package at all. Nothing here compiles Stan.
#
# Every expected number comes from an oracle computed in the test ---
# posterior::summarise_draws() or posterior::as_draws_df() on the same
# object --- never from a typed constant.

probs_from <- function(ci_level) {
  c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2)
}

# the contract ----------------------------------------------------------

test_that("extract_estimates returns the contract columns in order", {
  skip_if_not_installed("brms")
  out <- extract_estimates(mixture2p_fit())

  expect_s3_class(out, "tbl_df")
  expect_named(out, names(estimates_contract()))
  for (col in names(estimates_contract())) {
    expect_type(out[[col]], estimates_contract()[[col]])
  }
})

test_that("the population level has one row per estimated parameter", {
  skip_if_not_installed("brms")
  skip_if_not_installed("bmm")
  fit <- mixture2p_fit()
  pars <- bmm::parameters(fit)
  estimated <- pars$parameter[!pars$fixed]

  out <- extract_estimates(fit, level = "population")

  expect_equal(nrow(out), length(estimated))
  expect_setequal(out$term, estimated)
  expect_true(all(out$level == "population"))
  expect_true(all(is.na(out$id)))
})

test_that("ci_method and ci_level are carried through", {
  skip_if_not_installed("brms")
  out <- extract_estimates(mixture2p_fit(), ci_level = 0.89)
  expect_true(all(out$ci_method == "eti"))
  expect_true(all(out$ci_level == 0.89))
})

# the numbers -----------------------------------------------------------

test_that("population estimates match summarise_draws on the same fit", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  out <- extract_estimates(fit, level = "population")

  draws <- posterior::subset_draws(
    posterior::as_draws_array(fit),
    variable = paste0("b_", out$term, "_Intercept")
  )
  oracle <- posterior::summarise_draws(
    draws,
    "median",
    ~ posterior::quantile2(.x, probs = probs_from(0.95)),
    posterior::default_convergence_measures()
  )

  expect_equal(out$estimate, oracle$median)
  expect_equal(out$ci_low, oracle$q2.5)
  expect_equal(out$ci_high, oracle$q97.5)
  expect_equal(out$rhat, oracle$rhat)
  expect_equal(out$ess_bulk, oracle$ess_bulk)
  expect_equal(out$ess_tail, oracle$ess_tail)
})

test_that("ci_level moves the bounds to the matching quantiles", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  wide <- extract_estimates(fit, ci_level = 0.95)
  narrow <- extract_estimates(fit, ci_level = 0.89)

  expect_true(all(narrow$ci_low > wide$ci_low))
  expect_true(all(narrow$ci_high < wide$ci_high))
  # the point estimate is a median and does not depend on the interval
  expect_equal(narrow$estimate, wide$estimate)

  draws <- posterior::subset_draws(
    posterior::as_draws_array(fit),
    variable = paste0("b_", narrow$term, "_Intercept")
  )
  oracle <- posterior::summarise_draws(
    draws, ~ posterior::quantile2(.x, probs = probs_from(0.89))
  )
  expect_equal(narrow$ci_low, oracle$q5.5)
  expect_equal(narrow$ci_high, oracle$q94.5)
})

# constants -------------------------------------------------------------

test_that("drop_constants removes the parameters bmm fixes", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()

  dropped <- extract_estimates(fit, drop_constants = TRUE)
  kept <- extract_estimates(fit, drop_constants = FALSE)

  expect_false(any(c("mu1", "mu2", "kappa2") %in% dropped$term))
  expect_true(all(c("mu1", "mu2", "kappa2") %in% kept$term))
  expect_true(all(is.na(kept$rhat[kept$term %in% c("mu1", "mu2", "kappa2")])))
  # and the estimated ones survive either way
  expect_true(all(c("kappa", "thetat") %in% dropped$term))
})

test_that("a constant is zero posterior variance, not a missing rhat", {
  # A broken chain also has NA rhat. Dropping on is.na(rhat) would remove
  # it silently and report the remaining parameters as if the fit were
  # fine; dropping on zero variance keeps it visible.
  broken <- c(stats::rnorm(40), rep(NA_real_, 40))
  draws <- fake_draws(list(
    b_fixed_Intercept = 0,
    b_broken_Intercept = broken,
    b_ok_Intercept = seq_len(80) / 10
  ))

  out <- estimates_from_draws(draws, groups = character(0))

  expect_true(is.na(posterior::summarise_draws(draws, "rhat")$rhat[2]))
  expect_gt(stats::var(broken, na.rm = TRUE), 0)
  expect_setequal(out$term, c("broken", "ok"))
})

test_that("dropping every term warns and returns a zero-row tibble", {
  draws <- fake_draws(list(b_a_Intercept = 0, b_b_Intercept = -100))
  expect_warning(
    out <- estimates_from_draws(draws, groups = character(0)),
    "constant"
  )
  expect_equal(nrow(out), 0L)
  expect_named(out, names(estimates_contract()))
})

# subject level ---------------------------------------------------------

test_that("subject level returns one row per subject and parameter", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  out <- extract_estimates(fit, level = "subject")

  n_subjects <- length(unique(out$id))
  expect_equal(n_subjects, 8L)
  expect_equal(nrow(out), n_subjects * 2L)
  expect_setequal(out$term, c("kappa", "thetat"))
  expect_false(any(is.na(out$id)))
  expect_true(all(out$level == "subject"))
})

test_that("a subject value is the per-draw sum of intercept and deviation", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  out <- extract_estimates(fit, level = "subject")

  df <- posterior::as_draws_df(fit)
  for (par in c("kappa", "thetat")) {
    for (subject in c("1", "5", "8")) {
      summed <- df[[paste0("b_", par, "_Intercept")]] +
        df[[paste0("r_id__", par, "[", subject, ",Intercept]")]]
      row <- out[out$term == par & out$id == subject, ]
      expect_equal(row$estimate, stats::median(summed))
      expect_equal(
        row$ci_low, unname(stats::quantile(summed, 0.025, names = FALSE))
      )
      expect_equal(
        row$ci_high, unname(stats::quantile(summed, 0.975, names = FALSE))
      )
    }
  }
})

test_that("the sum is taken per draw, not by adding two summaries", {
  # Adding the two medians would give a different number whenever the
  # two are not perfectly aligned, and intervals that are too narrow.
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  out <- extract_estimates(fit, level = "subject")
  df <- posterior::as_draws_df(fit)

  naive_low <- unname(stats::quantile(
    df[["b_kappa_Intercept"]], 0.025,
    names = FALSE
  )) +
    unname(stats::quantile(
      df[["r_id__kappa[1,Intercept]"]], 0.025,
      names = FALSE
    ))
  row <- out[out$term == "kappa" & out$id == "1", ]
  expect_gt(row$ci_low, naive_low)
})

test_that("both levels stack and the level column separates them", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  both <- extract_estimates(fit, level = c("population", "subject"))

  expect_equal(nrow(both), 2L + 16L)
  expect_setequal(unique(both$level), c("population", "subject"))
  expect_equal(
    both[both$level == "population", ],
    extract_estimates(fit, level = "population")
  )
})

test_that("group-level SDs are not returned in Milestone 1", {
  skip_if_not_installed("brms")
  out <- extract_estimates(mixture2p_fit(), level = c("population", "subject"))
  expect_false(any(grepl("^sd_", out$term)))
  expect_false("sd" %in% out$level)
})

# grouping factors ------------------------------------------------------

test_that("an ambiguous grouping factor is an error that lists them", {
  draws <- fake_draws(list(
    b_kappa_Intercept = seq_len(80) / 10,
    `r_id__kappa[1,Intercept]` = stats::rnorm(80),
    `r_item__kappa[1,Intercept]` = stats::rnorm(80)
  ))
  err <- expect_error(
    estimates_from_draws(draws, groups = c("id", "item"), level = "subject")
  )
  expect_match(conditionMessage(err), "id")
  expect_match(conditionMessage(err), "item")
  expect_match(conditionMessage(err), "group")
})

test_that("a fit with no group-level effects errors on level = subject", {
  draws <- fake_draws(list(b_kappa_Intercept = seq_len(80) / 10))
  expect_error(
    estimates_from_draws(draws, groups = character(0), level = "subject"),
    "no group-level effects"
  )
})

test_that("an unknown group names the ones the fit has", {
  skip_if_not_installed("brms")
  expect_error(
    extract_estimates(mixture2p_fit(), level = "subject", group = "subject"),
    "id"
  )
})

test_that("a named group resolves without ambiguity", {
  draws <- fake_draws(list(
    b_kappa_Intercept = seq_len(80) / 10,
    `r_id__kappa[1,Intercept]` = stats::rnorm(80),
    `r_item__kappa[1,Intercept]` = stats::rnorm(80)
  ))
  out <- estimates_from_draws(
    draws,
    groups = c("id", "item"), level = "subject", group = "item"
  )
  expect_equal(nrow(out), 1L)
  expect_equal(out$id, "1")
})

# term naming -----------------------------------------------------------

test_that("strip_design_suffix reduces a term to the bare parameter name", {
  expect_equal(strip_design_suffix("b_kappa_Intercept"), "kappa")
  expect_equal(strip_design_suffix("b_thetat_Intercept"), "thetat")
  expect_equal(strip_design_suffix("b_kappa_setsize2"), "kappa")
  expect_equal(strip_design_suffix("b_Intercept"), "Intercept")
})

test_that("terms that do not reduce to unique names are an error", {
  # A design-structure fit: two coefficients collapse onto "kappa". A
  # silent dplyr fan-out at the join is the failure mode this prevents.
  draws <- fake_draws(list(
    b_kappa_Intercept = seq_len(80) / 10,
    b_kappa_setsize2 = stats::rnorm(80)
  ))
  expect_error(estimates_from_draws(draws, groups = character(0)), "kappa")
})

# argument checking -----------------------------------------------------

test_that("extract_estimates rejects a non-brmsfit", {
  expect_error(extract_estimates(1:10), "brmsfit")
  expect_error(extract_estimates(tibble::tibble(a = 1)), "brmsfit")
})

test_that("ci_level must lie strictly inside zero and one", {
  draws <- fake_draws(list(b_kappa_Intercept = seq_len(80) / 10))
  expect_error(
    estimates_from_draws(draws, character(0), ci_level = 0),
    "0 and 1"
  )
  expect_error(
    estimates_from_draws(draws, character(0), ci_level = 1),
    "0 and 1"
  )
  expect_error(
    estimates_from_draws(draws, character(0), ci_level = c(0.9, 0.95)),
    "single number"
  )
})

test_that("ci_method beyond eti is refused rather than ignored", {
  draws <- fake_draws(list(b_kappa_Intercept = seq_len(80) / 10))
  expect_error(
    estimates_from_draws(draws, character(0), ci_method = "hdi"),
    "eti"
  )
})

test_that("an unknown level is refused", {
  draws <- fake_draws(list(b_kappa_Intercept = seq_len(80) / 10))
  expect_error(estimates_from_draws(draws, character(0), level = "sd"))
})

# the apabayes contract -------------------------------------------------

test_that("the estimates tibble satisfies the apabayes parameters contract", {
  skip_if_not_installed("brms")
  skip_if_not_installed("apabayes")
  out <- extract_estimates(mixture2p_fit())
  expect_no_error(apabayes::apabayes_tidy(out, type = "parameters"))
})

# the converged column (spec 2, section 2) -------------------------------

test_that("extract_estimates carries the convergence gate as a column", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()

  out <- extract_estimates(fit, level = c("population", "subject"))

  expect_type(out$converged, "logical")
  expect_length(unique(out$converged), 1L)
  expect_identical(out$converged[[1L]], check_convergence(fit)$pass)
})

test_that("converged can be supplied instead of computed", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()

  expect_true(all(!extract_estimates(fit, converged = FALSE)$converged))
  expect_true(all(is.na(extract_estimates(fit, converged = NA)$converged)))
  expect_error(extract_estimates(fit, converged = "yes"), "converged")
})

test_that("estimates_from_draws fills converged without a fit", {
  draws <- fake_draws(list(b_a_Intercept = stats::rnorm(80)))
  out <- estimates_from_draws(draws, groups = character(0), converged = TRUE)
  expect_true(all(out$converged))
  default <- estimates_from_draws(draws, groups = character(0))
  expect_true(all(is.na(default$converged)))
})

# error attribution ------------------------------------------------------

test_that("a bad argument is blamed on extract_estimates, not a helper", {
  # The validation runs inside estimates_from_draws(), several frames
  # down. Without a threaded call the error points at that helper, which
  # a user never called and cannot find in their script.
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()

  bad_level <- tryCatch(
    extract_estimates(fit, ci_level = 0),
    error = function(e) e
  )
  expect_match(
    paste(deparse(conditionCall(bad_level)), collapse = " "),
    "^extract_estimates\\("
  )

  bad_group <- tryCatch(
    extract_estimates(fit, level = "subject", group = "nope"),
    error = function(e) e
  )
  expect_match(
    paste(deparse(conditionCall(bad_group)), collapse = " "),
    "^extract_estimates\\("
  )
})
