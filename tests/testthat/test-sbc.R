# Tests for the SBC naming layer, written against
# local/dev/spec-milestone-6-sbc-cross-check.md section 1.
#
# This is the stage where a mistake is silent. SBC matches the
# generator's truth names against the fit's draws matrix; a name that
# matches nothing is not an error there, it is a variable that never
# gets ranked. Measured 2026-09-16:
# `posterior::subset_draws(variable = "^zzz$", regex = TRUE)` returns a
# zero-column draws object and says nothing. So the patterns are tested
# against the committed fixtures --- real brms names, not hand-written
# guesses --- and check_variable_names() errors rather than proceeding.
#
# No test compiles Stan.

# sbc_variables() --------------------------------------------------------

test_that("sbc_variables builds one pattern per level and free parameter", {
  patterns <- sbc_variables(
    sbc_model(), "id",
    c("population", "sd", "cor", "subject")
  )

  expect_equal(
    patterns[sbc_variable_levels(patterns) == "population"],
    c(
      "population:kappa" = "^b_kappa_Intercept$",
      "population:thetat" = "^b_thetat_Intercept$"
    )
  )
  expect_equal(
    patterns[sbc_variable_levels(patterns) == "sd"],
    c(
      "sd:kappa" = "^sd_id__kappa_Intercept$",
      "sd:thetat" = "^sd_id__thetat_Intercept$"
    )
  )
  expect_setequal(
    unname(patterns[sbc_variable_levels(patterns) == "subject"]),
    c(
      "^r_id__kappa\\[.+,Intercept\\]$",
      "^r_id__thetat\\[.+,Intercept\\]$"
    )
  )
})

test_that("both orders of a correlation share one variable name", {
  # exactly one of the two can match, so requiring every *pattern* to
  # match would refuse every correlated fit; the name is what has to be
  # covered
  patterns <- sbc_variables(sbc_model(), "id", "cor")
  cors <- patterns[sbc_variable_levels(patterns) == "cor"]

  expect_setequal(
    unname(cors),
    c(
      "^cor_id__kappa_Intercept__thetat_Intercept$",
      "^cor_id__thetat_Intercept__kappa_Intercept$"
    )
  )
  expect_equal(unique(names(cors)), "cor:kappa__thetat")
})

test_that("a parameter the model fixes is never ranked", {
  model <- sbc_model(
    free = c("kappa", "thetat"),
    fixed = list(mu1 = 0, mu2 = 0, kappa2 = 0)
  )
  patterns <- sbc_variables(model, "id", "population")
  expect_false(any(grepl("mu1|mu2|kappa2", patterns)))
  expect_length(patterns, 2L)
})

test_that("a model with one free parameter has no correlation pattern", {
  patterns <- sbc_variables(sbc_model(free = "kappa"), "id", "cor")
  expect_length(patterns, 0L)
})

test_that("without a group there are no sd, cor or subject patterns", {
  patterns <- sbc_variables(
    sbc_model(), NULL,
    c("population", "sd", "cor", "subject")
  )
  expect_equal(unique(sbc_variable_levels(patterns)), "population")
})

test_that("a dot in a name is escaped rather than matching any character", {
  patterns <- sbc_variables(sbc_model(free = "k.a"), "id", "population")
  expect_equal(unname(patterns), "^b_k\\.a_Intercept$")
  expect_false(grepl(patterns, "b_kXa_Intercept"))
  expect_true(grepl(patterns, "b_k.a_Intercept"))
})

test_that("sbc_variables validates its arguments", {
  expect_error(sbc_variables("not a model", "id", "population"), "bmmodel")
  expect_error(sbc_variables(sbc_model(), "id", "nonsense"))
})

# the patterns against the committed fixtures ----------------------------

test_that("the population and sd patterns match the fixture's own names", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  available <- posterior::variables(posterior::as_draws_matrix(fit$fit))
  patterns <- sbc_variables(fit$bmm$model, "id", c("population", "sd"))

  matched <- unlist(lapply(
    patterns, function(p) grep(p, available, value = TRUE)
  ))
  expect_setequal(
    unname(matched),
    c(
      "b_kappa_Intercept", "b_thetat_Intercept",
      "sd_id__kappa_Intercept", "sd_id__thetat_Intercept"
    )
  )
  # the fixed mixture components are in the draws and are not ranked
  expect_true(all(c("b_mu1_Intercept", "Intercept_mu1") %in% available))
  expect_false(any(c("b_mu1_Intercept", "Intercept_mu1") %in% matched))
})

test_that("the subject pattern matches the fixture's r_id__ names", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  available <- posterior::variables(posterior::as_draws_matrix(fit$fit))
  patterns <- sbc_variables(fit$bmm$model, "id", "subject")

  matched <- unlist(lapply(
    patterns, function(p) grep(p, available, value = TRUE)
  ))
  expect_true("r_id__kappa[1,Intercept]" %in% matched)
  # eight subjects x two varying parameters, measured on the fixture
  expect_length(matched, 16L)
  expect_true(all(grepl("^r_id__", matched)))
})

test_that("the cor pattern matches the committed cor fixture in both orders", {
  skip_if_not_installed("brms")
  draws <- mixture2p_cor_draws()
  available <- posterior::variables(draws$draws)
  given <- grep("^cor_", available, value = TRUE)
  expect_length(given, 1L)

  # brms writes one order; which one is not guaranteed, so both patterns
  # are built and the reversed name is checked against them too
  reversed <- "cor_id__thetat_Intercept__kappa_Intercept"
  patterns <- sbc_variables(sbc_model(), "id", "cor")

  expect_true(any(vapply(patterns, grepl, logical(1), x = given)))
  expect_true(any(vapply(patterns, grepl, logical(1), x = reversed)))
})

# check_variable_names() -------------------------------------------------

test_that("check_variable_names passes when the two sets agree", {
  names <- c("b_kappa_Intercept", "sd_id__kappa_Intercept")
  expect_silent(check_variable_names(names, rev(names)))
})

test_that("check_variable_names errors and names both sets", {
  err <- expect_error(
    check_variable_names(
      c("kappa", "thetat"),
      c("b_kappa_Intercept", "b_thetat_Intercept")
    )
  )
  message <- conditionMessage(err)
  expect_match(message, "kappa")
  expect_match(message, "b_kappa_Intercept")
})

test_that("check_variable_names reports each direction of the difference", {
  err <- expect_error(
    check_variable_names(
      c("b_kappa_Intercept", "extra"),
      c("b_kappa_Intercept", "b_thetat_Intercept")
    )
  )
  message <- conditionMessage(err)
  expect_match(message, "extra")
  expect_match(message, "b_thetat_Intercept")
})

test_that("check_variable_names never warns and never returns quietly", {
  # SBC ranks nothing and reports nothing when a name does not match, so
  # this has to be an error: a warning in a run of a hundred fits is a
  # line nobody reads.
  expect_error(check_variable_names(character(0), "b_kappa_Intercept"))
  expect_error(check_variable_names("b_kappa_Intercept", character(0)))
})

# check_sbc_formula() ----------------------------------------------------

test_that("an intercept-only formula with a group term is accepted", {
  out <- check_sbc_formula(fake_bmmformula(
    kappa = kappa ~ 1 + (1 | id),
    thetat = thetat ~ 1 + (1 | id)
  ))
  expect_equal(out$group, "id")
  expect_false(out$correlated)
})

test_that("a correlated group term is recognised", {
  out <- check_sbc_formula(fake_bmmformula(
    kappa = kappa ~ 1 + (1 | p | id),
    thetat = thetat ~ 1 + (1 | p | id)
  ))
  expect_equal(out$group, "id")
  expect_true(out$correlated)
})

test_that("an intercept-only formula without a group term is accepted", {
  out <- check_sbc_formula(fake_bmmformula(
    kappa = kappa ~ 1,
    thetat = thetat ~ 1
  ))
  expect_null(out$group)
  expect_false(out$correlated)
})

test_that("a double bar is a group term without correlation", {
  out <- check_sbc_formula(fake_bmmformula(kappa = kappa ~ 1 + (1 || id)))
  expect_equal(out$group, "id")
  expect_false(out$correlated)
})

test_that("recovery_formula() is accepted at both re_cor settings", {
  skip_if_not_installed("bmm")
  model <- sbc_model()

  plain <- check_sbc_formula(recovery_formula(model))
  expect_equal(plain$group, "id")
  expect_false(plain$correlated)

  correlated <- check_sbc_formula(recovery_formula(model, re_cor = "all"))
  expect_equal(correlated$group, "id")
  expect_true(correlated$correlated)
})

test_that("a covariate in the fixed part is an error naming the term", {
  err <- expect_error(
    check_sbc_formula(fake_bmmformula(kappa = kappa ~ 1 + cond))
  )
  expect_match(conditionMessage(err), "kappa")
  expect_match(conditionMessage(err), "generator")
})

test_that("a task formula is an error", {
  err <- expect_error(
    check_sbc_formula(fake_bmmformula(
      kappa = kappa ~ 0 + task + (0 + task || id)
    ))
  )
  expect_match(conditionMessage(err), "kappa")
})

test_that("two different group terms are an error naming both", {
  err <- expect_error(
    check_sbc_formula(fake_bmmformula(
      kappa = kappa ~ 1 + (1 | id),
      thetat = thetat ~ 1 + (1 | session)
    ))
  )
  expect_match(conditionMessage(err), "id")
  expect_match(conditionMessage(err), "session")
})

test_that("two group terms on one parameter are an error", {
  expect_error(
    check_sbc_formula(fake_bmmformula(
      kappa = kappa ~ 1 + (1 | id) + (1 | session)
    ))
  )
})

test_that("a wrapped grouping factor is refused, not silently mis-parsed", {
  # `(1 | gr(id, cor = FALSE))` is ordinary brms syntax. Reading it off
  # the deparsed text gives the grouping factor "id, cor = FALSE", and
  # sbc_variables() would then build sd_/cor_/r_ patterns that match
  # nothing on the real fit, in silence.
  for (term in c("gr(id, cor = FALSE)", "mm(id1, id2)", "id:session")) {
    formula <- fake_bmmformula(kappa = stats::as.formula(
      paste0("kappa ~ 1 + (1 | ", term, ")"), env = globalenv()
    ))
    err <- expect_error(sbc_formula_groups(formula))
    expect_match(conditionMessage(err), "kappa")
    expect_match(conditionMessage(err), "not a name|bare column name")
  }
})

test_that("an offset is refused for being an offset, not for being a group", {
  formula <- fake_bmmformula(
    kappa = kappa ~ 1 + offset(log(n)) + (1 | id)
  )
  err <- expect_error(check_sbc_formula(formula))
  expect_match(conditionMessage(err), "intercept-only")
  expect_match(conditionMessage(err), "offset")
})

test_that("check_sbc_formula refuses something that is not a bmmformula", {
  expect_error(check_sbc_formula(kappa ~ 1), "bmmformula")
})

test_that("the group name alone is readable without the restriction", {
  # what sbc(generator = ) needs in 6.4: the group, with no check on the
  # fixed part, because bmmtools never interprets the draw row there
  out <- sbc_formula_groups(fake_bmmformula(
    kappa = kappa ~ 0 + task + (0 + task | p | id)
  ))
  expect_equal(out$group, "id")
  expect_true(out$correlated)
})

# sbc_layout() -----------------------------------------------------------

test_that("sbc_layout reads the subjects and the trials per subject", {
  data <- data.frame(id = rep(c("a", "b", "c"), each = 5L), y = 0)
  out <- sbc_layout(data, "id")
  expect_equal(out$n_subjects, 3L)
  expect_equal(out$n_trials, 5L)
  expect_equal(out$group, "id")
})

test_that("sbc_layout errors on unbalanced data and names the subjects", {
  data <- data.frame(id = c(rep("a", 5L), rep("b", 5L), rep("c", 3L)), y = 0)
  err <- expect_error(sbc_layout(data, "id"))
  expect_match(conditionMessage(err), "c")
  expect_match(conditionMessage(err), "unbalanced|same number")
})

test_that("sbc_layout errors when the group column is not there", {
  err <- expect_error(sbc_layout(data.frame(y = 1:4), "id"))
  expect_match(conditionMessage(err), "id")
})

test_that("sbc_layout errors on something that is not a data frame", {
  expect_error(sbc_layout(1:4, "id"), "data frame")
})

test_that("sbc_layout errors on an empty design", {
  expect_error(
    sbc_layout(data.frame(id = character(0), y = double(0)), "id"),
    "no rows|at least one"
  )
})

test_that("a factor group column is read by its levels", {
  data <- data.frame(id = factor(rep(1:4, each = 3L)), y = 0)
  out <- sbc_layout(data, "id")
  expect_equal(out$n_subjects, 4L)
  expect_equal(out$n_trials, 3L)
})

# check_sbc_dots() -------------------------------------------------------

test_that("check_sbc_dots refuses each argument and says where it belongs", {
  expected <- c(
    cores = "cores_per_fit",
    sample_prior = "sbc",
    file = "fit_cached",
    file_refit = "fit_cached",
    file_compress = "fit_cached"
  )
  for (name in names(expected)) {
    dots <- stats::setNames(list(1), name)
    err <- expect_error(check_sbc_dots(dots))
    expect_match(conditionMessage(err), name)
    expect_match(conditionMessage(err), expected[[name]])
  }
})

test_that("check_sbc_dots lists every offender at once", {
  err <- expect_error(check_sbc_dots(list(cores = 4, file = "x")))
  expect_match(conditionMessage(err), "cores")
  expect_match(conditionMessage(err), "file")
})

test_that("check_sbc_dots passes the arguments that do belong there", {
  dots <- list(chains = 2, iter = 500, backend = "cmdstanr", init = 0.1)
  expect_silent(check_sbc_dots(dots))
  expect_silent(check_sbc_dots(list()))
})

test_that("check_sbc_dots refuses an unnamed argument", {
  expect_error(check_sbc_dots(list(1)), "named")
})

# prior_parameter_draws() ------------------------------------------------

test_that("prior_parameter_draws subsets and samples the committed fixture", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  patterns <- sbc_variables(fit$bmm$model, "id", c("population", "sd"))

  draws <- prior_parameter_draws(fit, n_sims = 20L, variables = patterns)

  expect_s3_class(draws, "draws_matrix")
  expect_equal(dim(draws), c(20L, 4L))
  expect_setequal(
    posterior::variables(draws),
    c(
      "b_kappa_Intercept", "b_thetat_Intercept",
      "sd_id__kappa_Intercept", "sd_id__thetat_Intercept"
    )
  )
})

test_that("more n_sims than draws is an error naming iter and chains", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  patterns <- sbc_variables(fit$bmm$model, "id", "population")

  err <- expect_error(
    prior_parameter_draws(fit, n_sims = 5000L, variables = patterns)
  )
  message <- conditionMessage(err)
  expect_match(message, "iter")
  expect_match(message, "chains")
  # the number available is reported, not just the number asked for
  expect_match(message, "1000")
})

test_that("the same seed gives the same rows and no seed leaves the stream", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  patterns <- sbc_variables(fit$bmm$model, "id", "population")

  one <- prior_parameter_draws(fit, 20L, patterns, seed = 7L)
  two <- prior_parameter_draws(fit, 20L, patterns, seed = 7L)
  expect_equal(as.matrix(one), as.matrix(two))

  other <- prior_parameter_draws(fit, 20L, patterns, seed = 8L)
  expect_false(isTRUE(all.equal(as.matrix(one), as.matrix(other))))
})

test_that("a pattern that matches nothing is an error, not an empty matrix", {
  # measured 2026-09-16: subset_draws() returns a zero-column object and
  # says nothing, which is the silent failure this layer exists to close
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  err <- expect_error(
    prior_parameter_draws(fit, 5L, variables = "^not_a_variable$")
  )
  expect_match(conditionMessage(err), "not_a_variable")
})

test_that("a partial mismatch errors instead of returning a narrower matrix", {
  # the guard that matters: three patterns right and one wrong used to
  # come back three columns wide with nothing said, which is a run that
  # ranks `population` while the user asked for `population` and `sd`
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  patterns <- sbc_variables(fit$bmm$model, "id", c("population", "sd"))
  patterns[["sd:kappa"]] <- "^sd_WRONG__kappa_Intercept$"

  err <- expect_error(prior_parameter_draws(fit, 20L, variables = patterns))
  expect_match(conditionMessage(err), "sd:kappa")
  expect_match(conditionMessage(err), "sd_WRONG")
  # the three that did match are not named as missing
  expect_no_match(conditionMessage(err), "population:kappa")
})

test_that("select_variables returns the matched names in pattern order", {
  available <- c("b_kappa_Intercept", "b_thetat_Intercept", "lp__")
  out <- select_variables(
    c("population:thetat" = "^b_thetat_Intercept$",
      "population:kappa" = "^b_kappa_Intercept$"),
    available
  )
  expect_equal(out, c("b_thetat_Intercept", "b_kappa_Intercept"))
})

test_that("prior_parameter_draws validates n_sims", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  patterns <- sbc_variables(fit$bmm$model, "id", "population")
  expect_error(prior_parameter_draws(fit, 0L, patterns), "n_sims")
  expect_error(prior_parameter_draws(fit, 2.5, patterns), "n_sims")
})

test_that("the generic dispatches to a registered method", {
  fit <- sbc_mock_fit()
  draws <- prior_parameter_draws(
    fit, 10L,
    variables = c("^b_kappa_Intercept$", "^b_thetat_Intercept$")
  )
  expect_equal(dim(draws), c(10L, 2L))
  expect_setequal(
    posterior::variables(draws),
    c("b_kappa_Intercept", "b_thetat_Intercept")
  )
})

# what SBC itself sees ---------------------------------------------------

test_that("SBC reads the fixture through its brmsfit method", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("brms")
  skip_on_cran()
  fit <- mixture2p_fit()
  draws <- SBC::SBC_fit_to_draws_matrix(fit)

  expect_equal(dim(draws), c(1000L, 28L))
  # every name sbc() would rank is in what SBC hands back
  patterns <- sbc_variables(fit$bmm$model, "id", c("population", "sd"))
  available <- posterior::variables(draws)
  for (pattern in patterns) {
    expect_length(grep(pattern, available), 1L)
  }
})

# sbc() -------------------------------------------------------------------
#
# The wiring, over the mock fitter of helper-sbc.R. The model is a real
# bmm model, because the generator runs `simulate_recovery()` for real;
# only the sampler is a stand-in, so `SBC::compute_SBC()` runs its own
# pipeline over mock fits rather than being stood in for.

sbc_data <- function(n_subjects = 3L, n_trials = 4L, group = "id") {
  out <- data.frame(
    g = rep(seq_len(n_subjects), each = n_trials),
    y = 0
  )
  names(out)[[1L]] <- group
  out
}

#' The mock fits carry one deterministic sequence of draws, so their Rhat
#' says nothing and the N2 warning fires on every run. Only that warning
#' is muffled, by its class, so a run that warns about anything else
#' still shows it. `quiet` drops the cache and SBC's own progress
#' messages, which every run repeats; the two tests that assert on a
#' message pass `quiet = FALSE`.
sbc_run <- function(mock, ..., n_sims = 5L, formula = NULL, data = NULL,
                    quiet = TRUE) {
  model <- bmm::mixture2p(resp_error = "y")
  hush <- if (quiet) suppressMessages else identity
  hush(withCallingHandlers(
    sbc(
      model,
      formula %||% recovery_formula(model),
      data %||% sbc_data(),
      n_sims = n_sims,
      thin_ranks = 1,
      cores_per_fit = 1,
      .fitter = mock$fitter,
      ...
    ),
    bmmtools_sbc_diagnostics = function(w) invokeRestart("muffleWarning")
  ))
}

test_that("sbc fits the prior once and every data set once", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  mock <- sbc_mock_fitter()

  results <- sbc_run(mock, n_sims = 5L)

  expect_s3_class(results, "SBC_results")
  expect_identical(mock$calls$n, 6L)
  # the prior fit is the one that asks for prior draws; no data set fit does
  expect_identical(mock$calls$log[[1L]]$sample_prior, "only")
  for (call in sbc_dataset_calls(mock$calls)) {
    expect_null(call$sample_prior)
  }
})

test_that("every data set fit gets the layout's rows and the model's columns", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  mock <- sbc_mock_fitter()

  sbc_run(mock, n_sims = 4L, data = sbc_data(n_subjects = 5L, n_trials = 6L))

  for (call in sbc_dataset_calls(mock$calls)) {
    expect_identical(call$n_rows, 30L)
    expect_true(all(c("id", "y") %in% call$columns))
    expect_length(call$ids, 5L)
  }
})

test_that("the ranked names are the fit's own draw names", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()

  population <- sbc_run(sbc_mock_fitter(), level = "population")
  expect_setequal(
    attr(population, "bmmtools_sbc")$variables,
    c("b_kappa_Intercept", "b_thetat_Intercept")
  )

  with_sd <- sbc_run(sbc_mock_fitter(), level = c("population", "sd"))
  expect_setequal(
    attr(with_sd, "bmmtools_sbc")$variables,
    c(
      "b_kappa_Intercept", "b_thetat_Intercept",
      "sd_id__kappa_Intercept", "sd_id__thetat_Intercept"
    )
  )
  # what SBC actually ranked, not only what bmmtools recorded
  expect_setequal(
    unique(with_sd$stats$variable),
    attr(with_sd, "bmmtools_sbc")$variables
  )
})

test_that("the cor level ranks the correlation under brms's own name", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  model <- bmm::mixture2p(resp_error = "y")
  # a fit of `re_cor = "all"` does carry the cor_ draw, so the mock does
  mock <- sbc_mock_fitter(
    variables = posterior::variables(sbc_prior_draws(1L))
  )

  results <- sbc_run(
    mock,
    level = c("population", "sd", "cor"),
    formula = recovery_formula(model, re_cor = "all")
  )

  expect_true(
    "cor_id__kappa_Intercept__thetat_Intercept" %in%
      attr(results, "bmmtools_sbc")$variables
  )
  expect_true(
    "cor_id__kappa_Intercept__thetat_Intercept" %in% results$stats$variable
  )
})

test_that("an uncorrelated formula never asks a fit for a cor_ draw", {
  # the default `recovery_formula()` writes one `(1 | id)` per parameter,
  # which brms fits without a correlation. Deciding what to *draw* from
  # the levels the caller asked to *rank* left `cor` in the drawn set
  # whenever `level` omitted it, so the commonest call of all --- the
  # default formula at the default level --- asked the prior fit for a
  # `cor_` draw it does not have, and only after the fit had run. The
  # mock carries no cor_ draw here, exactly as a real fit would not.
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  model <- bmm::mixture2p(resp_error = "y")
  mock <- sbc_mock_fitter()
  expect_false(
    any(grepl("^cor_", posterior::variables(mock$fitter(
      recovery_formula(model), sbc_data(), model
    )$prior_draws)))
  )

  results <- sbc_run(mock, n_sims = 3L, formula = recovery_formula(model))

  expect_setequal(
    attr(results, "bmmtools_sbc")$variables,
    c(
      "b_kappa_Intercept", "b_thetat_Intercept",
      "sd_id__kappa_Intercept", "sd_id__thetat_Intercept"
    )
  )
  expect_identical(attr(results, "bmmtools_sbc")$level, c("population", "sd"))
})

test_that("a correlated formula draws the cor_ value even when not ranked", {
  # the other half of N5: `level` chooses what is ranked, and a
  # correlation the formula does imply is still drawn and simulated from
  skip_if_not_installed("bmm")
  model <- bmm::mixture2p(resp_error = "y")
  groups <- check_sbc_formula(recovery_formula(model, re_cor = "all"))

  sets <- sbc_pattern_sets(model, groups, "id", c("population", "sd"))

  expect_true("cor:kappa__thetat" %in% names(sets$draw))
  expect_false("cor:kappa__thetat" %in% names(sets$rank))
})

test_that("the stats have one row per simulation and ranked variable", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  results <- sbc_run(sbc_mock_fitter(), n_sims = 6L)

  variables <- attr(results, "bmmtools_sbc")$variables
  expect_equal(nrow(results$stats), 6L * length(variables))
  expect_setequal(results$stats$sim_id, seq_len(6L))
})

test_that("the attribute records what bmmtools decided", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  results <- sbc_run(sbc_mock_fitter(), n_sims = 5L, seed = 11)

  info <- attr(results, "bmmtools_sbc")
  expect_named(
    info,
    c(
      "model", "n_sims", "level", "variables", "seed", "prior", "layout",
      "diagnostics", "generator"
    )
  )
  expect_identical(info$generator, "simulate_recovery")
  expect_identical(info$model, "mixture2p")
  expect_identical(info$n_sims, 5L)
  expect_identical(info$level, c("population", "sd"))
  expect_identical(info$seed, 11)
  expect_identical(
    info$layout,
    list(n_subjects = 3L, n_trials = 4L, group = "id", ids = c("1", "2", "3"))
  )
})

test_that("no seed is recorded as NA, never as a silently fixed value", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  results <- sbc_run(sbc_mock_fitter())
  expect_identical(attr(results, "bmmtools_sbc")$seed, NA_real_)
})

test_that("the same seed gives the same data sets", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  one <- sbc_mock_fitter()
  two <- sbc_mock_fitter()
  three <- sbc_mock_fitter()

  sbc_run(one, seed = 3)
  sbc_run(two, seed = 3)
  sbc_run(three, seed = 4)

  first <- vapply(sbc_dataset_calls(one$calls), function(x) x$data$y[[1L]], 0)
  second <- vapply(sbc_dataset_calls(two$calls), function(x) x$data$y[[1L]], 0)
  other <- vapply(sbc_dataset_calls(three$calls), function(x) x$data$y[[1L]], 0)

  expect_equal(first, second)
  expect_false(isTRUE(all.equal(first, other)))
})

test_that("keep_fits keeps the fits and FALSE keeps none", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()

  kept <- sbc_run(sbc_mock_fitter(), n_sims = 3L, keep_fits = TRUE)
  expect_length(kept$fits, 3L)
  expect_s3_class(kept$fits[[1L]], "sbcmockfit")

  dropped <- sbc_run(sbc_mock_fitter(), n_sims = 3L)
  expect_true(all(vapply(dropped$fits, is.null, logical(1))))
  # the diagnostics N2 reads survive keep_fits = FALSE (measured 2026-09-16)
  expect_false(is.null(dropped$default_diagnostics))
})

# what sbc() refuses ------------------------------------------------------

test_that("sbc refuses each argument that collides or is set here", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  # `file` is a formal of sbc(), the prior fit's cache, so it can never
  # reach `...`; check_sbc_dots() still refuses it and is tested above
  for (name in c("cores", "sample_prior", "file_refit", "file_compress")) {
    dots <- stats::setNames(list(1), name)
    err <- expect_error(
      do.call(sbc_run, c(list(sbc_mock_fitter()), dots))
    )
    expect_match(conditionMessage(err), name)
  }
})

test_that("file caches the prior fit and nothing else", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  path <- file.path(withr::local_tempdir(), "prior")
  mock <- sbc_mock_fitter()

  sbc_run(mock, n_sims = 2L, file = path)

  # set_path() appends the set name, as prior_check() does
  expect_true(file.exists(paste0(path, "-sbc.rds")))
  expect_length(list.files(dirname(path), pattern = "\\.rds$"), 1L)
})

test_that("a list of priors is an error pointing at prior_check", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  err <- expect_error(
    sbc_run(sbc_mock_fitter(), prior = list(a = fake_prior(), b = NULL))
  )
  expect_match(conditionMessage(err), "prior_check")
})

test_that("a level without population is an error", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  err <- expect_error(sbc_run(sbc_mock_fitter(), level = "sd"))
  expect_match(conditionMessage(err), "population")
})

test_that("unbalanced data is an error naming the odd subject", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  data <- rbind(sbc_data(), data.frame(id = 3L, y = 0))
  err <- expect_error(sbc_run(sbc_mock_fitter(), data = data))
  expect_match(conditionMessage(err), "unbalanced|same number")
})

test_that("an NA in the grouping column is refused, not counted away", {
  # table() drops NA silently, so the layout would otherwise be read off
  # the rows that happen to have an id (6.2 review, LOW, deferred to 6.3)
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  data <- sbc_data()
  data$id[[2L]] <- NA
  err <- expect_error(sbc_run(sbc_mock_fitter(), data = data))
  expect_match(conditionMessage(err), "NA")
  expect_match(conditionMessage(err), "id")
})

test_that("cache_mode results without a location is an error", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  err <- expect_error(sbc_run(sbc_mock_fitter(), cache_mode = "results"))
  expect_match(conditionMessage(err), "cache_location")
})

test_that("more simulations than prior draws errors naming iter and chains", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  mock <- sbc_mock_fitter(n_draws = 10L)
  err <- expect_error(sbc_run(mock, n_sims = 50L))
  expect_match(conditionMessage(err), "iter")
  expect_match(conditionMessage(err), "chains")
})

test_that("sbc validates model, data, n_sims, seed and .fitter", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  model <- bmm::mixture2p(resp_error = "y")
  formula <- recovery_formula(model)
  expect_error(
    sbc("not a model", formula, sbc_data(), .fitter = identity), "bmmodel"
  )
  expect_error(sbc(model, formula, 1:4, .fitter = identity), "data frame")
  expect_error(
    sbc(model, formula, sbc_data(), n_sims = 0L, .fitter = identity), "n_sims"
  )
  expect_error(
    sbc(model, formula, sbc_data(), seed = "x", .fitter = identity), "seed"
  )
  expect_error(
    sbc(model, formula, sbc_data(), .fitter = "not a function"), ".fitter"
  )
})

# levels a formula has no draws for --------------------------------------

test_that("sd and subject levels are dropped when nothing varies", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  model <- bmm::mixture2p(resp_error = "y")
  formula <- bmm::bmf(kappa ~ 1, thetat ~ 1)
  mock <- sbc_mock_fitter(
    variables = c("b_kappa_Intercept", "b_thetat_Intercept")
  )

  messages <- character(0)
  results <- withCallingHandlers(
    sbc_run(
      mock,
      level = c("population", "sd", "subject"), formula = formula,
      quiet = FALSE
    ),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(any(grepl("Dropping.*sd", messages)))
  expect_true(any(grepl("Dropping.*subject", messages)))
  expect_identical(attr(results, "bmmtools_sbc")$level, "population")
})

test_that("the cor level is dropped when the group term is uncorrelated", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  mock <- sbc_mock_fitter()

  expect_message(
    results <- sbc_run(
      mock,
      level = c("population", "sd", "cor"), quiet = FALSE
    ),
    "Dropping"
  )
  expect_identical(
    attr(results, "bmmtools_sbc")$level, c("population", "sd")
  )
})

test_that("a level the formula does imply is an error, not a drop", {
  # the difference that matters: "no group term" is a design fact and is
  # a message, while a group term whose draws are missing is a bug and
  # has to stop the run before SBC ranks a subset in silence
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  mock <- sbc_mock_fitter(
    variables = c(
      "b_kappa_Intercept", "b_thetat_Intercept", "sd_id__kappa_Intercept"
    )
  )
  err <- expect_error(sbc_run(mock, level = c("population", "sd")))
  expect_match(conditionMessage(err), "sd:thetat")
})

test_that("level selects what is ranked, not what is drawn", {
  # N5: holding the SDs at zero while the fitted model has a prior on
  # them breaks the joint p(theta) p(y | theta), and the population ranks
  # stop being calibrated too. So a group term is simulated from whether
  # or not "sd" is ranked.
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  mock <- sbc_mock_fitter()

  results <- sbc_run(mock, level = "population", n_sims = 4L)

  expect_setequal(
    attr(results, "bmmtools_sbc")$variables,
    c("b_kappa_Intercept", "b_thetat_Intercept")
  )
  # the generated subjects differ from one another, which they could not
  # do if the between-subject SDs had been left out of the simulation
  for (call in sbc_dataset_calls(mock$calls)) {
    per_subject <- tapply(call$data$y, call$data$id, mean)
    expect_gt(stats::sd(per_subject), 0)
  }
})

# sbc_generator() --------------------------------------------------------

test_that("the generator reproduces the prior rows and the layout", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("brms")
  skip_if_not_installed("bmm")
  skip_on_cran()
  fit <- mixture2p_fit()
  model <- fit$bmm$model
  patterns <- sbc_variables(model, "id", c("population", "sd"))
  draws <- prior_parameter_draws(fit, n_sims = 20L, variables = patterns)
  layout <- sbc_layout(fit$data, "id")

  generator <- sbc_generator(
    draws, posterior::variables(draws), model, layout,
    correlated = FALSE, group = "id"
  )
  datasets <- SBC::generate_datasets(generator, 20L)

  expect_length(datasets$generated, 20L)
  expect_equal(nrow(datasets$generated[[1L]]), 240L)
  expect_equal(
    as.matrix(posterior::as_draws_matrix(datasets$variables)),
    as.matrix(draws),
    ignore_attr = TRUE
  )
  expect_setequal(
    posterior::variables(datasets$variables), posterior::variables(draws)
  )
})

test_that("the generator refuses to be asked for more rows than it has", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  model <- bmm::mixture2p(resp_error = "y")
  draws <- posterior::subset_draws(
    sbc_prior_draws(3L),
    variable = c("b_kappa_Intercept", "b_thetat_Intercept")
  )
  generator <- sbc_generator(
    draws, posterior::variables(draws), model,
    list(n_subjects = 2L, n_trials = 3L, group = "id"),
    correlated = FALSE, group = "id"
  )
  expect_error(SBC::generate_datasets(generator, 4L), "3")
})

test_that("a non-default grouping factor names the generated column", {
  # simulate_recovery() always calls the column `id`; a formula that says
  # (1 | subject) would then be fitted against data that has no such
  # column, and brms would be the one to complain
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  model <- bmm::mixture2p(resp_error = "y")
  formula <- bmm::bmf(kappa ~ 1 + (1 | subject), thetat ~ 1 + (1 | subject))
  mock <- sbc_mock_fitter(group = "subject")

  sbc_run(
    mock,
    level = c("population", "sd"), n_sims = 2L, formula = formula,
    data = sbc_data(group = "subject")
  )

  for (call in sbc_dataset_calls(mock$calls)) {
    expect_true("subject" %in% call$columns)
    expect_false("id" %in% call$columns)
  }
})

# sbc_backend() ----------------------------------------------------------

test_that("the backend lambda takes cores and passes the fit's arguments", {
  skip_if_not_installed("SBC")
  skip_on_cran()
  seen <- NULL
  fitter <- function(formula, data, model, prior = NULL, ...) {
    seen <<- c(list(formula = formula, model = model, prior = prior), list(...))
    "fitted"
  }
  backend <- sbc_backend(
    "the formula", "the model", "the prior",
    dots = list(chains = 2L, iter = 500L), fitter = fitter
  )

  # cores_arg = "cores" makes SBC put cores into the call, so the lambda
  # needs the formal; a function of one argument fails inside do.call()
  expect_named(formals(backend$func), c("generated", "cores"))
  out <- SBC::SBC_fit(backend, data.frame(y = 1), cores = 3L)

  expect_identical(out, "fitted")
  expect_identical(seen$cores, 3L)
  expect_identical(seen$chains, 2L)
  expect_identical(seen$prior, "the prior")
  expect_null(seen$sample_prior)
})

test_that("SBC's own mock backend accepts the generator's data sets", {
  # SBC_backend_mock()'s SBC_fit() ignores `generated` and returns a
  # fixed result, so this is a smoke test of the datasets only; the
  # per-dataset fitter call is tested through .fitter above
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  model <- bmm::mixture2p(resp_error = "y")
  draws <- posterior::subset_draws(
    sbc_prior_draws(4L),
    variable = c("b_kappa_Intercept", "b_thetat_Intercept")
  )
  generator <- sbc_generator(
    draws, posterior::variables(draws), model,
    list(n_subjects = 2L, n_trials = 3L, group = "id"),
    correlated = FALSE, group = "id"
  )
  datasets <- SBC::generate_datasets(generator, 4L)
  # 200 draws so that SBC does not warn about too few for its checks
  backend <- SBC::SBC_backend_mock(result = sbc_prior_draws(200L))

  results <- SBC::compute_SBC(
    datasets, backend,
    keep_fits = FALSE, thin_ranks = 1, cores_per_fit = 1
  )
  expect_s3_class(results, "SBC_results")
  expect_equal(nrow(results$stats), 4L * 2L)
})

# sbc_diagnostics() (N2) -------------------------------------------------

test_that("sbc_diagnostics counts the fits past each bar and warns once", {
  results <- list(default_diagnostics = data.frame(
    sim_id = 1:4,
    max_rhat = c(1.00, 1.08, 1.02, 1.20),
    min_ess_to_rank = c(0.9, 0.8, 0.3, 0.9)
  ))

  warnings <- character()
  out <- withCallingHandlers(
    sbc_diagnostics(results),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_identical(out$n_high_rhat, 2L)
  expect_identical(out$n_low_ess_to_rank, 1L)
  expect_identical(out$rhat_max, 1.05)
  expect_identical(out$ess_to_rank_min, 0.5)
  expect_length(warnings, 1L)
  # bmmtools' bar is stated and SBC's stricter one is named
  expect_match(warnings, "1.05")
  expect_match(warnings, "1.01")
})

test_that("sbc_diagnostics is silent when every fit passes both bars", {
  results <- list(default_diagnostics = data.frame(
    sim_id = 1:3,
    max_rhat = c(1.00, 1.01, 1.02),
    min_ess_to_rank = c(0.9, 0.8, 0.7)
  ))
  expect_silent(out <- sbc_diagnostics(results))
  expect_identical(out$n_high_rhat, 0L)
  expect_identical(out$n_low_ess_to_rank, 0L)
  expect_identical(out$n_missing, 0L)
})

test_that("a fit with no Rhat at all is counted rather than passed", {
  results <- list(default_diagnostics = data.frame(
    sim_id = 1:2,
    max_rhat = c(1.00, NA_real_),
    min_ess_to_rank = c(0.9, NA_real_)
  ))
  expect_warning(out <- sbc_diagnostics(results), "1")
  expect_identical(out$n_missing, 1L)
  expect_identical(out$n_high_rhat, 0L)
})

test_that("sbc_diagnostics survives a result with no diagnostics table", {
  out <- expect_silent(sbc_diagnostics(list()))
  expect_identical(out$n_high_rhat, NA_integer_)
  expect_identical(out$n_fits, NA_integer_)
})

test_that("the warning fires from sbc() and the counts reach the attribute", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  model <- bmm::mixture2p(resp_error = "y")
  mock <- sbc_mock_fitter()

  # not through sbc_run(), which muffles this warning by class
  expect_warning(
    results <- sbc(
      model, recovery_formula(model), sbc_data(),
      n_sims = 3L, thin_ranks = 1, cores_per_fit = 1, .fitter = mock$fitter
    ),
    class = "bmmtools_sbc_diagnostics"
  )

  diagnostics <- attr(results, "bmmtools_sbc")$diagnostics
  expect_named(
    diagnostics,
    c(
      "n_fits", "rhat_max", "ess_to_rank_min", "n_high_rhat",
      "n_low_ess_to_rank", "n_missing"
    )
  )
  expect_identical(diagnostics$n_fits, 3L)
})

# sbc_pattern_sets() and its neighbours, without SBC ---------------------

test_that("the cor level is dropped when only one parameter varies", {
  expect_message(
    sets <- sbc_pattern_sets(
      sbc_model(),
      list(group = "id", correlated = TRUE, varying = "kappa"),
      "id", c("population", "sd", "cor")
    ),
    "Dropping"
  )
  expect_identical(sets$level, c("population", "sd"))
  # N5: the sd pattern is built for the parameter that varies and not
  # for the one that does not, whose sd_ draw no fit of this formula has
  expect_setequal(
    names(sets$draw),
    c("population:kappa", "population:thetat", "sd:kappa")
  )
})

test_that("what is drawn is a superset of what is ranked", {
  sets <- sbc_pattern_sets(
    sbc_model(),
    list(group = "id", correlated = TRUE, varying = c("kappa", "thetat")),
    "id", "population"
  )
  expect_identical(sets$level, "population")
  expect_setequal(names(sets$rank), c("population:kappa", "population:thetat"))
  expect_setequal(
    unique(names(sets$draw)),
    c(
      "population:kappa", "population:thetat", "sd:kappa", "sd:thetat",
      "cor:kappa__thetat"
    )
  )
})

test_that("sbc_cor_names gives nothing for fewer than two varying terms", {
  expect_identical(sbc_cor_names("kappa", "id", "anything"), list())
})

test_that("sbc_cor_names finds the pair in whichever order brms wrote it", {
  available <- "cor_id__thetat_Intercept__kappa_Intercept"
  found <- sbc_cor_names(c("kappa", "thetat"), "id", available)
  expect_length(found, 1L)
  expect_identical(found[[1L]]$name, available)
})

test_that("check_group_ids leaves a missing column to sbc_layout", {
  expect_silent(check_group_ids(data.frame(y = 1:3), "id"))
})

test_that("a draw the model cannot generate from names the draw", {
  # measured 2026-09-16 on the first real Stan run: bmm's own default
  # half-student_t(3, 0, 2.5) on the group-level SD of `kappa` draws
  # values near 9.5 on the log scale, and `rmixture2p()` then dies with
  # `node stack overflow`, naming neither the simulation nor a parameter
  skip_if_not_installed("bmm")
  model <- bmm::mixture2p(resp_error = "y")
  draws <- posterior::as_draws_matrix(cbind(
    b_kappa_Intercept = c(1.4, 1.4),
    b_thetat_Intercept = c(0.8, 0.8),
    sd_id__kappa_Intercept = c(0.3, -1)
  ))
  generator <- sbc_generator(
    draws, c("b_kappa_Intercept", "b_thetat_Intercept"), model,
    list(n_subjects = 2L, n_trials = 3L, group = "id"),
    correlated = FALSE, group = "id"
  )

  expect_silent(generator$f())
  err <- expect_error(generator$f())
  message <- conditionMessage(err)
  expect_match(message, "data set 2")
  expect_match(message, "sd_id__kappa_Intercept")
  expect_match(message, "prior")
})

# the subject level (6.4) ------------------------------------------------
#
# brms's `r_<g>__<par>[<label>,Intercept]` draws are deviations from the
# intercept --- measured on the fixture 2026-09-16: the eight
# `r_id__kappa[i,Intercept]` posterior means average -0.0278 against a
# `b_kappa_Intercept` of 1.8692 --- while `simulate_recovery()`'s
# `truth$subjects` are absolute. Emitting the absolute values would give
# ranks that look like a badly miscalibrated model and are a units error.
# The label is the other half: the prior fit names its `r_` draws after
# the values of the user's grouping column, and the simulation numbers
# its subjects 1..n, so the generated column has to be relabelled before
# the two can be matched by name.

test_that("sbc_layout records the subject labels in the order of data", {
  layout <- sbc_layout(
    data.frame(id = rep(c("s2", "s1", "s3"), each = 2L), y = 0), "id"
  )
  expect_identical(layout$ids, c("s2", "s1", "s3"))
  expect_identical(layout$n_subjects, 3L)

  # a factor is read by the labels it uses, not by its level order
  levelled <- sbc_layout(
    data.frame(
      id = factor(rep(c("b", "a"), each = 3L), levels = c("z", "a", "b")),
      y = 0
    ),
    "id"
  )
  expect_identical(levelled$ids, c("b", "a"))
  expect_identical(length(levelled$ids), levelled$n_subjects)
})

test_that("subject truths are the simulated values minus the population", {
  # pure arithmetic on a hand-built truth: no package, no simulation
  truth <- list(
    population = tibble::tibble(
      term = c("kappa", "thetat"), true_value = c(2, 0.5)
    ),
    subjects = tibble::tibble(
      id = c("1", "2", "3", "1", "2", "3"),
      term = rep(c("kappa", "thetat"), each = 3L),
      true_value = c(2.3, 1.8, 2.0, 0.1, 0.9, 0.5)
    )
  )

  truths <- sbc_subject_truths(truth, "id", ids = c("a", "b", "c"))

  expect_named(
    truths,
    c(
      "r_id__kappa[a,Intercept]", "r_id__kappa[b,Intercept]",
      "r_id__kappa[c,Intercept]", "r_id__thetat[a,Intercept]",
      "r_id__thetat[b,Intercept]", "r_id__thetat[c,Intercept]"
    )
  )
  expect_equal(unname(truths), c(0.3, -0.2, 0, -0.4, 0.4, 0))
  # the whole point: adding the population value back gives the truth
  population <- stats::setNames(
    truth$population$true_value, truth$population$term
  )
  expect_equal(
    unname(truths) + unname(population[truth$subjects$term]),
    truth$subjects$true_value
  )
})

test_that("a subject the truth has no row for is left out, not invented", {
  truth <- list(
    population = tibble::tibble(term = "kappa", true_value = 1),
    subjects = tibble::tibble(
      id = character(), term = character(), true_value = double()
    )
  )
  expect_identical(
    sbc_subject_truths(truth, "id", ids = c("a", "b")),
    stats::setNames(numeric(0), character(0))
  )
})

test_that("the generator emits deviations, under the prior fit's own labels", {
  # the regression test for the deviation trap, on the committed fixture:
  # the same seed reproduces the generator's simulation outside it, so
  # the emitted subject truths can be compared with the simulation's
  # absolute values exactly, not statistically
  skip_if_not_installed("SBC")
  skip_if_not_installed("brms")
  skip_if_not_installed("bmm")
  skip_on_cran()
  fit <- mixture2p_fit()
  model <- fit$bmm$model
  free <- model_parameters(model)$free
  patterns <- sbc_variables(model, "id", c("population", "sd", "subject"))
  draws <- prior_parameter_draws(fit, n_sims = 3L, variables = patterns)
  layout <- sbc_layout(fit$data, "id")
  rank <- select_variables(patterns, posterior::variables(draws))

  generator <- sbc_generator(
    draws, rank, model, layout, correlated = FALSE, group = "id"
  )
  out <- withr::with_seed(7, generator$f())

  # a plain named matrix, as the generator builds it: `[` on a
  # draws_matrix drops the variable names one row at a time
  values <- matrix(
    as.numeric(draws),
    nrow = posterior::ndraws(draws),
    dimnames = list(NULL, posterior::variables(draws))
  )
  row <- values[1L, ]
  sim <- withr::with_seed(7, simulate_recovery(
    model,
    pars = stats::setNames(row[paste0("b_", free, "_Intercept")], free),
    n_subjects = layout$n_subjects, n_trials = layout$n_trials,
    sds = stats::setNames(row[paste0("sd_id__", free, "_Intercept")], free)
  ))
  population <- stats::setNames(
    sim$truth$population$true_value, sim$truth$population$term
  )
  subjects <- sim$truth$subjects
  expected <- stats::setNames(
    subjects$true_value - unname(population[subjects$term]),
    paste0("r_id__", subjects$term, "[", subjects$id, ",Intercept]")
  )

  emitted <- unlist(out$variables)
  expect_setequal(names(emitted), rank)
  expect_equal(emitted[names(expected)], expected)
  # the eight labels are the fixture's own, so SBC would find every one
  expect_true(all(
    grep("^r_", rank, value = TRUE) %in%
      posterior::variables(posterior::as_draws_matrix(fit$fit))
  ))
  # and not the prior fit's r_ draws copied through: those are a
  # different random variable, and copying them would rank the simulation
  # against values it was not generated from
  expect_false(isTRUE(all.equal(
    unname(emitted[names(expected)]), unname(row[names(expected)])
  )))
  # the data are relabelled to match
  expect_setequal(as.character(unique(out$generated$id)), layout$ids)
})

test_that("a label the prior fit did not see is an error naming both sets", {
  # the trap of the 6.3 handoff: the prior fit is fitted to the user's
  # data and carries its labels, the simulation numbers subjects 1..n.
  # A layout whose labels differ from the fit's cannot be ranked, and the
  # error has to come from the name check, not from SBC ranking nothing
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  model <- bmm::mixture2p(resp_error = "y")
  draws <- sbc_prior_draws(3L, subjects = c("1", "2"))
  draws <- posterior::subset_draws(
    draws, variable = grep("^cor_", posterior::variables(draws),
      invert = TRUE, value = TRUE
    )
  )
  rank <- posterior::variables(draws)
  generator <- sbc_generator(
    draws, rank, model,
    list(n_subjects = 2L, n_trials = 3L, group = "id", ids = c("s1", "s2")),
    correlated = FALSE, group = "id"
  )

  err <- expect_error(generator$f())
  message <- conditionMessage(err)
  expect_match(message, "r_id__kappa\\[s1,Intercept\\]")
  expect_match(message, "r_id__kappa\\[1,Intercept\\]")
})

test_that("level = subject ranks every subject of every varying parameter", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  mock <- sbc_mock_fitter(subjects = c("1", "2", "3"))

  results <- sbc_run(
    mock,
    level = c("population", "sd", "subject"), n_sims = 4L
  )

  info <- attr(results, "bmmtools_sbc")
  expect_identical(info$level, c("population", "sd", "subject"))
  subject_names <- as.vector(outer(
    c("r_id__kappa[", "r_id__thetat["), c("1", "2", "3"),
    function(a, b) paste0(a, b, ",Intercept]")
  ))
  expect_setequal(
    info$variables,
    c(
      "b_kappa_Intercept", "b_thetat_Intercept",
      "sd_id__kappa_Intercept", "sd_id__thetat_Intercept",
      subject_names
    )
  )
  # what SBC ranked, not only what bmmtools recorded: n_subjects x
  # n_varying more rows per simulation, which is why the level is opt-in
  expect_setequal(unique(results$stats$variable), info$variables)
  expect_equal(nrow(results$stats), 4L * 10L)
})

test_that("the generated data carry the layout's labels, not 1..n", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  data <- sbc_data()
  data$id <- paste0("s", data$id)
  mock <- sbc_mock_fitter(subjects = c("s1", "s2", "s3"))

  results <- sbc_run(
    mock,
    level = c("population", "subject"), n_sims = 2L, data = data
  )

  for (call in sbc_dataset_calls(mock$calls)) {
    expect_setequal(call$ids, c("s1", "s2", "s3"))
  }
  expect_true(
    "r_id__kappa[s2,Intercept]" %in% attr(results, "bmmtools_sbc")$variables
  )
  expect_true("r_id__kappa[s2,Intercept]" %in% results$stats$variable)
})

test_that("subject draws the prior fit lacks are an error, not a drop", {
  # the formula implies them, so their absence is a bug in the fit, and
  # the rule of "a level the formula does imply is an error" applies
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  mock <- sbc_mock_fitter()

  err <- expect_error(
    sbc_run(mock, level = c("population", "subject"), n_sims = 2L)
  )
  expect_match(conditionMessage(err), "subject:kappa")
})

test_that("the subject level is ranked but never drawn from the prior row", {
  sets <- sbc_pattern_sets(
    sbc_model(),
    list(group = "id", correlated = FALSE, varying = c("kappa", "thetat")),
    "id", c("population", "sd", "subject")
  )
  expect_identical(sets$level, c("population", "sd", "subject"))
  expect_setequal(
    names(sets$rank),
    c(
      "population:kappa", "population:thetat", "sd:kappa", "sd:thetat",
      "subject:kappa", "subject:thetat"
    )
  )
  # the subjects are simulated from the drawn SDs, so they are never read
  # off the row; but their names are taken from the prior fit, so that
  # the generator's constructed names are checked against real ones
  expect_false(any(grepl("^subject:", names(sets$draw))))
  expect_setequal(
    names(sets$take), union(names(sets$draw), names(sets$rank))
  )
})

# the user generator (6.4, spec (d) under open question (i) = A) --------
#
# With a `generator`, bmmtools never interprets the prior draw row: the
# generator receives every b_, sd_, cor_ and r_ draw of the prior fit
# under the fit's own names and returns the data frame to fit, so any
# formula works and the truth names equal the draw names by
# construction. Each level then means every draw of that class for the
# model's free parameters, resolved against the prior fit's names
# rather than built from the formula --- which is what makes a
# cell-means formula's `b_kappa_task1` rankable.

#' Draws a fit of `kappa ~ 0 + task + (1 | id), thetat ~ 1 + (1 | id)`
#' would carry, plus two names that must never be handed over or ranked:
#' a fixed parameter's constant and the log density.
task_draws <- function(n = 60L, subjects = c("1", "2", "3")) {
  k <- seq_len(n)
  values <- cbind(
    b_kappa_task1 = log(4) + 0.01 * k,
    b_kappa_task2 = log(5) + 0.01 * k,
    b_thetat_Intercept = stats::qlogis(0.7) + 0.005 * k,
    b_mu1_Intercept = 0,
    sd_id__kappa_Intercept = 0.2 + 0.002 * k,
    sd_id__thetat_Intercept = 0.3 + 0.002 * k,
    lp__ = -k
  )
  for (i in subjects) {
    values <- cbind(values, 0.1 * sin(k + as.integer(i)))
    colnames(values)[ncol(values)] <- paste0("r_id__kappa[", i, ",Intercept]")
  }
  posterior::as_draws_matrix(values)
}

task_formula <- function() {
  bmm::bmf(kappa ~ 0 + task + (1 | id), thetat ~ 1 + (1 | id))
}

task_data <- function() {
  out <- sbc_data(n_subjects = 3L, n_trials = 4L)
  out$task <- rep(c("task1", "task2"), length.out = nrow(out))
  out
}

#' A proper prior on the cell means, which brms's default leaves flat and
#' `check_improper_priors()` would otherwise warn about. Per coefficient,
#' because that check reads the coefficient rows, which brms leaves empty
#' when a class-level prior covers them (a false positive of Milestone
#' 4's check, noted in STATE-milestone-6.md, not fixed here).
task_prior <- function() {
  brms::set_prior(
    "normal(0, 1)", class = "b", coef = "tasktask1", nlpar = "kappa"
  ) +
    brms::set_prior(
      "normal(0, 1)", class = "b", coef = "tasktask2", nlpar = "kappa"
    )
}

#' A generator that records what it was given and returns the layout
#' with fresh responses
recording_generator <- function() {
  seen <- new.env(parent = emptyenv())
  seen$draws <- list()
  seen$data <- list()
  f <- function(draws, data) {
    seen$draws[[length(seen$draws) + 1L]] <- draws
    seen$data[[length(seen$data) + 1L]] <- data
    data$y <- seq_len(nrow(data)) + length(seen$draws)
    data
  }
  list(f = f, seen = seen)
}

test_that("a generator with a task formula runs where the default errors", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  mock <- sbc_mock_fitter(draws = task_draws())

  expect_error(
    sbc_run(
      mock,
      formula = task_formula(), data = task_data(), n_sims = 2L,
      prior = task_prior()
    ),
    "intercept-only"
  )

  gen <- recording_generator()
  results <- sbc_run(
    sbc_mock_fitter(draws = task_draws()),
    formula = task_formula(), data = task_data(), n_sims = 3L,
    generator = gen$f, level = c("population", "sd"), prior = task_prior()
  )

  info <- attr(results, "bmmtools_sbc")
  expect_identical(info$generator, "user")
  expect_identical(info$level, c("population", "sd"))
  expect_setequal(
    info$variables,
    c(
      "b_kappa_task1", "b_kappa_task2", "b_thetat_Intercept",
      "sd_id__kappa_Intercept", "sd_id__thetat_Intercept"
    )
  )
  expect_setequal(unique(results$stats$variable), info$variables)
})

test_that("the generator receives the whole row of the four classes", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  gen <- recording_generator()
  mock <- sbc_mock_fitter(draws = task_draws())

  sbc_run(
    mock,
    formula = task_formula(), data = task_data(), n_sims = 3L,
    generator = gen$f, level = "population", prior = task_prior()
  )

  expect_length(gen$seen$draws, 3L)
  row <- gen$seen$draws[[1L]]
  expect_type(row, "double")
  # every b_, sd_ and r_ draw, whether or not it is ranked ...
  expect_setequal(
    names(row),
    c(
      "b_kappa_task1", "b_kappa_task2", "b_thetat_Intercept",
      "sd_id__kappa_Intercept", "sd_id__thetat_Intercept",
      paste0("r_id__kappa[", 1:3, ",Intercept]")
    )
  )
  # ... and neither a fixed parameter's constant nor the log density
  expect_false(any(c("b_mu1_Intercept", "lp__") %in% names(row)))
  # the data is the user's own frame, columns and all
  expect_identical(gen$seen$data[[1L]], task_data())
  # and what the generator returned is what was fitted
  for (k in seq_along(sbc_dataset_calls(mock$calls))) {
    call <- sbc_dataset_calls(mock$calls)[[k]]
    expect_identical(call$data$y, seq_len(12L) + k)
    expect_true("task" %in% call$columns)
  }
})

test_that("the generator's truths are the row itself, under the fit's names", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  gen <- recording_generator()

  results <- sbc_run(
    sbc_mock_fitter(draws = task_draws()),
    formula = task_formula(), data = task_data(), n_sims = 2L,
    generator = gen$f, level = "population", prior = task_prior()
  )

  first <- results$stats[results$stats$sim_id == 1L, ]
  truths <- stats::setNames(first$simulated_value, first$variable)
  ranked <- names(gen$seen$draws[[1L]])[1:3]
  expect_equal(truths[ranked], gen$seen$draws[[1L]][ranked])
})

test_that("a generator that is not a two-argument function is refused", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  mock <- sbc_mock_fitter()
  expect_error(sbc_run(mock, generator = "no"), "function")
  err <- expect_error(sbc_run(mock, generator = function(draws) draws))
  expect_match(conditionMessage(err), "generator\\(draws, data\\)")
  # dots count as accepting both
  expect_no_error(check_sbc_generator(function(...) NULL))
  expect_no_error(check_sbc_generator(NULL))
})

test_that("what the generator returns is checked, naming the data set", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  err <- expect_error(
    sbc_run(sbc_mock_fitter(), generator = function(draws, data) 1)
  )
  expect_match(conditionMessage(err), "data frame")
  expect_match(conditionMessage(err), "data set 1")

  err <- expect_error(sbc_run(
    sbc_mock_fitter(),
    generator = function(draws, data) data.frame(id = data$id)
  ))
  expect_match(conditionMessage(err), "\"y\"")

  err <- expect_error(sbc_run(
    sbc_mock_fitter(),
    generator = function(draws, data) data.frame(y = data$y)
  ))
  expect_match(conditionMessage(err), "\"id\"")
})

test_that("an error inside the generator names the data set and the draw", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  err <- expect_error(sbc_run(
    sbc_mock_fitter(),
    generator = function(draws, data) stop("boom")
  ))
  message <- conditionMessage(err)
  expect_match(message, "data set 1")
  expect_match(message, "b_kappa_Intercept")
  expect_match(conditionMessage(err$parent), "boom")
})

test_that("in generator mode a level the prior fit lacks is dropped, named", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  model <- bmm::mixture2p(resp_error = "y")
  # a formula with a group term, fitted (by the mock) without sd_ draws:
  # generator mode cannot know from the formula, so it reads the fit
  mock <- sbc_mock_fitter(
    variables = c("b_kappa_Intercept", "b_thetat_Intercept")
  )
  messages <- character(0)
  results <- withCallingHandlers(
    sbc_run(
      mock,
      level = c("population", "sd", "subject"),
      generator = function(draws, data) data, quiet = FALSE
    ),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(any(grepl("Dropping.*sd.*sd_id__", messages)))
  expect_true(any(grepl("Dropping.*subject.*r_id__", messages)))
  expect_identical(attr(results, "bmmtools_sbc")$level, "population")

  # without a group term the reason is the formula, not the fit
  messages <- character(0)
  withCallingHandlers(
    sbc_run(
      mock,
      level = c("population", "sd"),
      formula = bmm::bmf(kappa ~ 1, thetat ~ 1),
      generator = function(draws, data) data, quiet = FALSE
    ),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(any(grepl("Dropping.*sd.*group term", messages)))
})

test_that("in generator mode cor ranks every correlation of the group", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  model <- bmm::mixture2p(resp_error = "y")
  mock <- sbc_mock_fitter(
    variables = posterior::variables(sbc_prior_draws(1L))
  )

  results <- sbc_run(
    mock,
    level = c("population", "sd", "cor"),
    formula = recovery_formula(model, re_cor = "all"),
    generator = function(draws, data) data
  )

  expect_true(
    "cor_id__kappa_Intercept__thetat_Intercept" %in%
      attr(results, "bmmtools_sbc")$variables
  )
})

test_that("in generator mode subject truths must carry the layout's labels", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  mock <- sbc_mock_fitter(subjects = c("1", "2", "3"))

  # the fit's own labels: fine, and the r_ draws are ranked
  results <- sbc_run(
    mock,
    level = c("population", "subject"), n_sims = 2L,
    generator = function(draws, data) data
  )
  expect_true(
    "r_id__kappa[2,Intercept]" %in% attr(results, "bmmtools_sbc")$variables
  )

  # other labels: an error naming both sets, because SBC would rank
  # nothing for those subjects and say nothing
  err <- expect_error(sbc_run(
    sbc_mock_fitter(subjects = c("1", "2", "3")),
    level = c("population", "subject"), n_sims = 2L,
    generator = function(draws, data) {
      data$id <- data$id + 10L
      data
    }
  ))
  expect_match(conditionMessage(err), "11")
  expect_match(conditionMessage(err), "labels")

  # without the subject level the labels are the generator's business
  expect_no_error(sbc_run(
    sbc_mock_fitter(subjects = c("1", "2", "3")),
    level = "population", n_sims = 2L,
    generator = function(draws, data) {
      data$id <- data$id + 10L
      data
    }
  ))
})

test_that("an unbalanced design is fine when the generator owns the data", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  data <- sbc_data()[-1L, ]

  expect_error(sbc_run(sbc_mock_fitter(), data = data), "unbalanced")

  results <- sbc_run(
    sbc_mock_fitter(),
    data = data, n_sims = 2L, level = "population",
    generator = function(draws, data) data
  )
  layout <- attr(results, "bmmtools_sbc")$layout
  expect_identical(layout$n_subjects, 3L)
  expect_identical(layout$n_trials, 4L)
})

test_that("sbc_generator_sets resolves each level against the fit's names", {
  available <- c(
    "b_kappa_task1", "b_kappa_task2", "b_thetat_Intercept",
    "b_mu1_Intercept", "Intercept_mu1", "sd_id__kappa_Intercept",
    "cor_id__kappa_Intercept__thetat_Intercept",
    "r_id__kappa[1,Intercept]", "r_id__kappa[2,Intercept]",
    "z_1[1,1]", "lp__"
  )
  groups <- list(group = "id", correlated = TRUE, varying = "kappa")

  sets <- sbc_generator_sets(
    sbc_model(), groups, "id", c("population", "sd", "cor", "subject"),
    available
  )
  expect_identical(sets$level, c("population", "sd", "cor", "subject"))
  expect_setequal(
    sets$take,
    c(
      "b_kappa_task1", "b_kappa_task2", "b_thetat_Intercept",
      "sd_id__kappa_Intercept", "cor_id__kappa_Intercept__thetat_Intercept",
      "r_id__kappa[1,Intercept]", "r_id__kappa[2,Intercept]"
    )
  )
  expect_setequal(sets$rank, sets$take)

  # ranking less than is taken
  population <- sbc_generator_sets(
    sbc_model(), groups, "id", "population", available
  )
  expect_setequal(
    population$rank, c("b_kappa_task1", "b_kappa_task2", "b_thetat_Intercept")
  )
  expect_setequal(population$take, sets$take)

  # no population draw at all is a bug, not a drop
  expect_error(
    sbc_generator_sets(sbc_model(), groups, "id", "population", "lp__"),
    "population"
  )
})

test_that("broad patterns name every draw of a class for a parameter", {
  patterns <- sbc_variables(
    sbc_model(), "id", c("population", "sd", "cor", "subject"),
    broad = TRUE
  )
  expect_equal(
    unname(patterns[names(patterns) == "population:kappa"]), "^b_kappa_"
  )
  expect_equal(unname(patterns[names(patterns) == "cor:id"]), "^cor_id__")
  expect_true(grepl(patterns[["population:kappa"]], "b_kappa_task2"))
  expect_false(grepl(patterns[["population:kappa"]], "b_kappa2_Intercept"))
  expect_true(grepl(patterns[["sd:kappa"]], "sd_id__kappa_task1"))
  expect_true(grepl(patterns[["subject:kappa"]], "r_id__kappa[s1,task1]"))
  expect_false(grepl(patterns[["subject:kappa"]], "r_id__kappa2[1,Intercept]"))
})

test_that("the user generator refuses to be asked for more rows than it has", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran()
  model <- bmm::mixture2p(resp_error = "y")
  draws <- posterior::subset_draws(
    sbc_prior_draws(3L),
    variable = c("b_kappa_Intercept", "b_thetat_Intercept")
  )
  generator <- sbc_user_generator(
    draws, posterior::variables(draws), function(draws, data) data,
    sbc_data(), model,
    list(n_subjects = 3L, n_trials = 4L, group = "id", ids = c("1", "2", "3")),
    group = "id", needs_group = FALSE
  )
  expect_error(SBC::generate_datasets(generator, 4L), "3")
})
