# Tests for extract_estimates(), written against
# local/dev/spec-milestone-1-score-layer.md section 2.
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

test_that("group-level SDs are absent unless requested", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  out <- extract_estimates(fit, level = c("population", "subject"))
  expect_false(any(grepl("^sd_", out$term)))
  expect_false("sd" %in% out$level)

  sds <- extract_estimates(fit, level = "sd")
  expect_setequal(sds$term, c("kappa", "thetat"))
  expect_true(all(sds$level == "sd"))
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

test_that("split_coefficient separates the parameter from the coefficient", {
  # the four names the retired strip_design_suffix() test pinned
  parts <- split_coefficient(
    c("kappa_Intercept", "thetat_Intercept", "kappa_setsize2", "Intercept")
  )
  expect_equal(parts$par, c("kappa", "thetat", "kappa", ""))
  expect_equal(parts$coef, c("Intercept", "Intercept", "setsize2", "Intercept"))
})

test_that("a lone coefficient keeps the bare parameter name", {
  # one non-Intercept coefficient reads back as the parameter, as before 5.4
  draws <- fake_draws(list(
    b_kappa_setsize2 = seq_len(80) / 10,
    b_Intercept = stats::rnorm(80)
  ))
  out <- estimates_from_draws(draws, groups = character(0))
  expect_equal(out$term, c("kappa", "Intercept"))
})

test_that("terms that do not reduce to unique names are an error", {
  # Two coefficients collapse onto "kappa". A silent dplyr fan-out at the
  # join is the failure mode this prevents.
  draws <- fake_draws(list(
    b_kappa = seq_len(80) / 10,
    b_kappa_Intercept = stats::rnorm(80)
  ))
  expect_error(estimates_from_draws(draws, groups = character(0)), "kappa")
})

test_that("an intercept with other coefficients is a contrast design", {
  # what `coding = "contrast"` fits: the intercept is the parameter's
  # population value and every other coefficient is an effect
  draws <- fake_draws(list(
    b_kappa_Intercept = seq_len(80) / 10,
    b_kappa_setsize2 = stats::rnorm(80)
  ))
  out <- estimates_from_draws(
    draws,
    groups = character(0), level = c("population", "effect")
  )
  expect_equal(out$term, c("kappa", "kappa_setsize2"))
  expect_equal(out$level, c("population", "effect"))

  # each level selects its own rows, and neither counts the other as a
  # parameter it found and dropped
  population <- estimates_from_draws(draws, groups = character(0))
  expect_equal(population$term, "kappa")
  effect <- estimates_from_draws(
    draws,
    groups = character(0), level = "effect"
  )
  expect_equal(effect$term, "kappa_setsize2")

  # a cell-means fit has no effects, and asking for them is not a warning
  # about everything having been dropped
  cells <- fake_draws(list(
    b_kappa_task1 = stats::rnorm(80), b_kappa_task2 = stats::rnorm(80)
  ))
  expect_silent(
    none <- estimates_from_draws(
      cells,
      groups = character(0), level = "effect"
    )
  )
  expect_equal(nrow(none), 0L)
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
  expect_error(estimates_from_draws(draws, character(0), level = "bogus"))
})

# the apabayes contract -------------------------------------------------

test_that("the estimates tibble satisfies the apabayes parameters contract", {
  skip_if_not_installed("brms")
  skip_if_not_installed("apabayes")
  out <- extract_estimates(mixture2p_fit())
  # Resolved at call time, not as a literal `apabayes::apabayes_tidy()`, for
  # the reason `bmm_fun()` gives in R/adapters.R: a literal is scanned
  # statically, and `skip_if_not_installed()` guards the run, not the scan.
  # apabayes is on no repository and is deliberately not a declared
  # dependency (neither package depends on the other), so R CMD check reads
  # it as an undeclared import. Measured 2026-09-18: with the CRAN and
  # Bioconductor indices reachable the check is 0/0/0, and with them
  # unreachable the same tree is 1 WARNING --- which the workflows, running
  # `error-on = "warning"`, would fail on during any repository outage.
  apabayes_tidy <- getExportedValue("apabayes", "apabayes_tidy")
  expect_no_error(apabayes_tidy(out, type = "parameters"))
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

# the sd and cor levels (spec 5, section 5.2) -----------------------------

# Draws with two SDs and one correlation. `cor_name` lets a test put the
# two coefficient names of the correlation in either order.
sd_cor_draws <- function(cor_name = "cor_id__kappa_Intercept__thetat_Intercept",
                         seed = 11) {
  withr::local_seed(seed)
  values <- list(
    b_kappa_Intercept = stats::rnorm(80),
    b_thetat_Intercept = stats::rnorm(80),
    sd_id__kappa_Intercept = abs(stats::rnorm(80)),
    sd_id__thetat_Intercept = abs(stats::rnorm(80)),
    cor = stats::runif(80, -1, 1)
  )
  names(values)[[5L]] <- cor_name
  fake_draws(values)
}

test_that("the sd level gives one row per parameter under its bare name", {
  draws <- sd_cor_draws()
  out <- estimates_from_draws(draws, "id", level = "sd", ranef = fake_ranef())

  expect_named(out, names(estimates_contract()))
  for (col in names(estimates_contract())) {
    expect_type(out[[col]], estimates_contract()[[col]])
  }
  expect_equal(out$term, c("kappa", "thetat"))
  expect_true(all(out$level == "sd"))
  expect_true(all(is.na(out$id)))

  oracle <- posterior::summarise_draws(
    posterior::subset_draws(
      draws,
      variable = c("sd_id__kappa_Intercept", "sd_id__thetat_Intercept")
    ),
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

test_that("a cor row is named a__b and matches summarise_draws", {
  draws <- sd_cor_draws()
  out <- estimates_from_draws(draws, "id", level = "cor", ranef = fake_ranef())

  expect_equal(nrow(out), 1L)
  expect_equal(out$term, "kappa__thetat")
  expect_equal(out$level, "cor")
  expect_true(is.na(out$id))
  oracle <- posterior::summarise_draws(
    posterior::subset_draws(
      draws,
      variable = "cor_id__kappa_Intercept__thetat_Intercept"
    ),
    "median",
    posterior::default_convergence_measures()
  )
  expect_equal(out$estimate, oracle$median)
  expect_equal(out$rhat, oracle$rhat)
})

test_that("both orders of the names inside cor_ give the same row", {
  forward <- sd_cor_draws("cor_id__kappa_Intercept__thetat_Intercept")
  reversed <- sd_cor_draws("cor_id__thetat_Intercept__kappa_Intercept")

  with_ranef <- lapply(list(forward, reversed), function(d) {
    estimates_from_draws(d, "id", level = "cor", ranef = fake_ranef())
  })
  expect_equal(with_ranef[[1L]]$term, "kappa__thetat")
  expect_identical(with_ranef[[1L]], with_ranef[[2L]])

  # the ranef table lists thetat first: the term is still sorted
  flipped <- estimates_from_draws(
    forward, "id",
    level = "cor", ranef = fake_ranef(c("thetat", "kappa"))
  )
  expect_identical(flipped, with_ranef[[1L]])

  # without a ranef table the names are parsed, with the same result
  parsed <- lapply(list(forward, reversed), function(d) {
    estimates_from_draws(d, "id", level = "cor")
  })
  expect_identical(parsed[[1L]], with_ranef[[1L]])
  expect_identical(parsed[[2L]], with_ranef[[1L]])
})

test_that("the pair term sorts in the C locale", {
  withr::local_seed(3)
  draws <- fake_draws(list(
    sd_id__thetat_Intercept = abs(stats::rnorm(80)),
    sd_id__Kappa_Intercept = abs(stats::rnorm(80)),
    cor_id__thetat_Intercept__Kappa_Intercept = stats::runif(80, -1, 1)
  ))
  ranef <- fake_ranef(c("thetat", "Kappa"))
  out <- estimates_from_draws(draws, "id", level = "cor", ranef = ranef)
  expect_equal(out$term, "Kappa__thetat")
})

test_that("a ranef table without cor = TRUE gives no cor rows", {
  draws <- sd_cor_draws()
  expect_no_warning(
    out <- estimates_from_draws(
      draws, "id",
      level = "cor", ranef = fake_ranef(cor = FALSE)
    )
  )
  expect_equal(nrow(out), 0L)
  expect_named(out, names(estimates_contract()))

  # coefficients in different correlation blocks are not a pair either
  apart <- estimates_from_draws(
    draws, "id",
    level = "cor", ranef = fake_ranef(id = c(1, 2))
  )
  expect_equal(nrow(apart), 0L)
})

test_that("a pair the fit does not estimate is absent, not NA", {
  withr::local_seed(5)
  draws <- fake_draws(list(
    sd_id__a_Intercept = abs(stats::rnorm(80)),
    sd_id__b_Intercept = abs(stats::rnorm(80)),
    sd_id__c_Intercept = abs(stats::rnorm(80)),
    cor_id__a_Intercept__c_Intercept = stats::runif(80, -1, 1)
  ))
  out <- estimates_from_draws(
    draws, "id",
    level = c("sd", "cor"), ranef = fake_ranef(c("a", "b", "c"))
  )
  expect_equal(out$term[out$level == "sd"], c("a", "b", "c"))
  expect_equal(out$term[out$level == "cor"], "a__c")
  expect_false(anyNA(out$estimate))
})

test_that("two grouping factors without group is an error at sd and cor", {
  withr::local_seed(8)
  draws <- fake_draws(list(
    sd_id__kappa_Intercept = abs(stats::rnorm(80)),
    sd_item__kappa_Intercept = abs(stats::rnorm(80))
  ))
  ranef <- rbind(
    fake_ranef("kappa", group = "id"),
    fake_ranef("kappa", group = "item")
  )
  for (level in c("sd", "cor")) {
    err <- expect_error(
      estimates_from_draws(draws, c("id", "item"), level = level, ranef = ranef)
    )
    expect_match(conditionMessage(err), "item")
  }

  out <- estimates_from_draws(
    draws, c("id", "item"),
    level = "sd", group = "item", ranef = ranef
  )
  expect_equal(out$term, "kappa")
  item_draws <- as.vector(draws[, , "sd_item__kappa_Intercept"])
  expect_equal(out$estimate, stats::median(item_draws))
})

test_that("a fit with no group-level effects errors at the sd level", {
  draws <- fake_draws(list(b_kappa_Intercept = seq_len(80) / 10))
  expect_error(
    estimates_from_draws(draws, character(0), level = "sd"),
    "no group-level effects"
  )
})

test_that("a group with no SD draws is an error", {
  draws <- fake_draws(list(
    b_kappa_Intercept = seq_len(80) / 10,
    `r_id__kappa[1,Intercept]` = seq_len(80) / 10
  ))
  expect_error(
    estimates_from_draws(draws, "id", level = "sd"),
    "standard deviation"
  )
})

test_that("a dpar and a missing resp column are handled", {
  withr::local_seed(9)
  draws <- fake_draws(list(
    sd_id__sigma_Intercept = abs(stats::rnorm(80)),
    sd_id__Intercept = abs(stats::rnorm(80)),
    cor_id__Intercept__sigma_Intercept = stats::runif(80, -1, 1)
  ))
  ranef <- data.frame(
    id = 1, group = "id", coef = "Intercept",
    dpar = c("sigma", ""), cor = TRUE,
    stringsAsFactors = FALSE
  )
  out <- estimates_from_draws(
    draws, "id",
    level = c("sd", "cor"), ranef = ranef
  )
  # a coefficient with neither nlpar nor dpar keeps its coefficient name,
  # as group_coefficients() does at the subject level
  expect_setequal(out$term[out$level == "sd"], c("sigma", "Intercept"))
  expect_equal(out$term[out$level == "cor"], "Intercept__sigma")
})

test_that("a resp column is part of the brms name", {
  withr::local_seed(10)
  draws <- fake_draws(list(
    sd_id__y_kappa_Intercept = abs(stats::rnorm(80))
  ))
  ranef <- fake_ranef("kappa", resp = "y")
  out <- estimates_from_draws(draws, "id", level = "sd", ranef = ranef)
  expect_equal(out$term, "kappa")
})

test_that("a contrast design names its SDs and correlations", {
  withr::local_seed(12)
  draws <- fake_draws(list(
    sd_id__kappa_Intercept = abs(stats::rnorm(80)),
    sd_id__kappa_setsize2 = abs(stats::rnorm(80)),
    cor_id__kappa_Intercept__kappa_setsize2 = stats::runif(80, -1, 1)
  ))
  ranef <- fake_ranef(
    c("kappa", "kappa"),
    coef = c("Intercept", "setsize2")
  )
  # the SD of the intercept is the SD of the parameter, the SD of the
  # contrast the SD of the effect: the terms contrast_truth() transforms to
  sds <- estimates_from_draws(draws, "id", level = "sd", ranef = ranef)
  expect_equal(sds$term, c("kappa", "kappa_setsize2"))
  expect_equal(unique(sds$level), "sd")

  cors <- estimates_from_draws(draws, "id", level = "cor", ranef = ranef)
  expect_equal(cors$term, "kappa__kappa_setsize2")

  expect_equal(
    estimates_from_draws(draws, "id", level = "sd")$term,
    c("kappa", "kappa_setsize2")
  )
})

test_that("a malformed ranef table is refused", {
  draws <- sd_cor_draws()
  expect_error(
    estimates_from_draws(draws, "id", level = "sd", ranef = "id"),
    "ranef"
  )
  expect_error(
    estimates_from_draws(
      draws, "id",
      level = "sd", ranef = data.frame(group = "id")
    ),
    "coef"
  )
})

test_that("every SD dropped as a constant warns", {
  draws <- fake_draws(list(
    sd_id__kappa_Intercept = 0.5,
    sd_id__thetat_Intercept = 0.5
  ))
  expect_warning(
    out <- estimates_from_draws(draws, "id", level = "sd"),
    "constant"
  )
  expect_equal(nrow(out), 0L)
})

test_that("a level with nothing to extract does not claim constants", {
  # an uncorrelated fit has population draws but no cor_ draws; that is
  # not "every parameter was dropped as a constant"
  withr::local_seed(13)
  draws <- fake_draws(list(
    b_kappa_Intercept = stats::rnorm(80),
    sd_id__kappa_Intercept = abs(stats::rnorm(80))
  ))
  expect_no_warning(
    out <- estimates_from_draws(
      draws, "id",
      level = "cor", ranef = fake_ranef("kappa", cor = FALSE)
    )
  )
  expect_equal(nrow(out), 0L)

  # the same holds for a population level with no b_ draws
  expect_no_warning(
    none <- estimates_from_draws(draws[, , 2L], "id", level = "population")
  )
  expect_equal(nrow(none), 0L)
})

test_that("extract_estimates reads sd rows and no cor rows off a real fit", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  expect_no_warning(out <- extract_estimates(fit, level = c("sd", "cor")))
  expect_equal(out$term, c("kappa", "thetat"))
  expect_true(all(out$level == "sd"))
})

test_that("the correlated fixture gives sd and cor rows with finite rhat", {
  record <- mixture2p_cor_draws()
  groups <- unique(as.character(record$ranef$group))
  out <- estimates_from_draws(
    record$draws, groups,
    level = c("population", "sd", "cor"), ranef = record$ranef
  )

  sds <- out[out$level == "sd", ]
  cors <- out[out$level == "cor", ]
  expect_setequal(sds$term, c("kappa", "thetat"))
  expect_equal(cors$term, "kappa__thetat")
  expect_true(all(is.finite(sds$rhat)))
  expect_true(all(is.finite(cors$rhat)))
  expect_true(cors$ci_low >= -1 && cors$ci_high <= 1)

  # parsing the names without the ranef table agrees
  parsed <- estimates_from_draws(
    record$draws, groups,
    level = c("population", "sd", "cor")
  )
  expect_identical(parsed, out)
})

# the task dimension (spec 5, section 5.4) --------------------------------

# Cell-means draws of a two-task fit with correlated random effects, as
# brms names them for `kappa ~ 0 + task + (0 + task | id)`.
task_draws <- function(seed = 21) {
  withr::local_seed(seed)
  fake_draws(list(
    b_kappa_task1 = stats::rnorm(80, 2),
    b_kappa_task2 = stats::rnorm(80, 1),
    `r_id__kappa[1,task1]` = stats::rnorm(80),
    `r_id__kappa[1,task2]` = stats::rnorm(80),
    `r_id__kappa[2,task1]` = stats::rnorm(80),
    `r_id__kappa[2,task2]` = stats::rnorm(80),
    sd_id__kappa_task1 = abs(stats::rnorm(80)),
    sd_id__kappa_task2 = abs(stats::rnorm(80)),
    cor_id__kappa_task2__kappa_task1 = stats::runif(80, -1, 1)
  ))
}

task_ranef <- function() {
  fake_ranef(c("kappa", "kappa"), coef = c("task1", "task2"))
}

test_that("cell-means coefficients get one term per task at every level", {
  draws <- task_draws()
  levels <- c("population", "subject", "sd", "cor")
  with_ranef <- estimates_from_draws(
    draws, "id",
    level = levels, ranef = task_ranef()
  )
  parsed <- estimates_from_draws(draws, "id", level = levels)

  for (out in list(with_ranef, parsed)) {
    at <- function(lv) out$term[out$level == lv]
    expect_equal(at("population"), c("kappa_task1", "kappa_task2"))
    expect_setequal(at("subject"), c("kappa_task1", "kappa_task2"))
    expect_equal(nrow(out[out$level == "subject", ]), 4L)
    expect_equal(at("sd"), c("kappa_task1", "kappa_task2"))
    # brms put task2 first; the pair term is sorted
    expect_equal(at("cor"), "kappa_task1__kappa_task2")
  }
  expect_identical(parsed, with_ranef)

  # a subject value is the task's intercept plus that task's deviation
  at <- with_ranef$level == "subject" & with_ranef$term == "kappa_task2" &
    with_ranef$id == "2"
  row <- with_ranef[at, ]
  summed <- as.vector(draws[, , "b_kappa_task2"]) +
    as.vector(draws[, , "r_id__kappa[2,task2]"])
  expect_equal(row$estimate, stats::median(summed))
})

test_that("the subject-draws array names its terms per task", {
  out <- subject_draws_from_draws(task_draws(), "id")
  expect_equal(dimnames(out)$term, c("kappa_task1", "kappa_task2"))
  expect_equal(dimnames(out)$id, c("1", "2"))

  cors <- correlations_from_parts(out, estimator = c("draws", "point"))
  expect_equal(unique(cors$term), "kappa_task1__kappa_task2")
})

test_that("an intercept together with a contrast is scored at every level", {
  withr::local_seed(22)
  draws <- fake_draws(list(
    b_kappa_Intercept = stats::rnorm(80),
    b_kappa_task2 = stats::rnorm(80),
    `r_id__kappa[1,Intercept]` = stats::rnorm(80),
    `r_id__kappa[1,task2]` = stats::rnorm(80),
    sd_id__kappa_Intercept = abs(stats::rnorm(80)),
    sd_id__kappa_task2 = abs(stats::rnorm(80)),
    cor_id__kappa_Intercept__kappa_task2 = stats::runif(80, -1, 1)
  ))
  ranef <- fake_ranef(c("kappa", "kappa"), coef = c("Intercept", "task2"))
  expected <- list(
    population = "kappa",
    effect = "kappa_task2",
    subject = c("kappa", "kappa_task2"),
    sd = c("kappa", "kappa_task2"),
    cor = "kappa__kappa_task2"
  )
  for (level in names(expected)) {
    for (rf in list(NULL, ranef)) {
      out <- estimates_from_draws(draws, "id", level = level, ranef = rf)
      expect_equal(unique(out$term), expected[[level]])
      expect_equal(unique(out$level), level)
    }
  }
  # the subject rows are the per-draw sum of each coefficient and its own
  # deviation, so the slope is a slope and not the parameter again
  expect_equal(
    dimnames(subject_draws_from_draws(draws, "id"))$term,
    c("kappa", "kappa_task2")
  )
})

test_that("the duplicate-term error says what it prevents", {
  withr::local_seed(23)
  draws <- fake_draws(list(
    sd_id__kappa = abs(stats::rnorm(80)),
    sd_id__kappa_Intercept = abs(stats::rnorm(80))
  ))
  err <- expect_error(estimates_from_draws(draws, "id", level = "sd"), "kappa")
  expect_no_match(conditionMessage(err), "Milestone")
  expect_match(conditionMessage(err), "fan out")
})

test_that("a task term takes the link of its parameter", {
  links <- c(kappa = "log", kappa2 = "log1p", thetat = "logit")
  expect_equal(
    link_of(c("kappa_task1", "kappa2_task1", "thetat_taskB"), links),
    c("log", "log1p", "logit")
  )
})
