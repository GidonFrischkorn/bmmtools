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
