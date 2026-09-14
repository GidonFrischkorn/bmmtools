# Tests for the correlation scorer, written against
# local/dev/spec-milestone-5-correlation-recovery.md section 5.3 (stage
# 5.3a: extract_subject_draws(), extract_correlations(),
# recover_correlations(), the bmmtools_cor_recovery class and the
# plot_recovery() generic).
#
# Everything except the tests marked "fixture" runs on hand-built subject
# draws, the same way the estimates tests run on fake_draws(): the core
# takes an array, a cor-estimates tibble and a covariate table, so no fit
# is needed.

# helpers ------------------------------------------------------------------

# Per-draw correlation of two terms across subjects, the slow way.
oracle_draw_cors <- function(x, a, b, transform_a = identity,
                             transform_b = identity) {
  apply(x, 1:2, function(slice) {
    stats::cor(transform_a(slice[, a]), transform_b(slice[, b]))
  })
}

oracle_summary <- function(per_draw, ci_level = 0.95) {
  probs <- c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2)
  list(
    estimate = stats::median(per_draw),
    ci_low = stats::quantile(per_draw, probs[[1L]], names = FALSE),
    ci_high = stats::quantile(per_draw, probs[[2L]], names = FALSE),
    rhat = posterior::rhat(per_draw)
  )
}

# A correlation-recovery input with `n_rep` replications, built by hand
# so every derived column can be computed on paper.
hand_cor_rows <- function(estimate, ci_low, ci_high,
                          term = "kappa__thetat",
                          estimator = "draws",
                          replication = seq_along(estimate),
                          scale = "link") {
  parts <- strsplit(term, "__", fixed = TRUE)[[1L]]
  tibble::tibble(
    term = term, var1 = parts[[1L]], var2 = parts[[2L]],
    estimator = estimator,
    estimate = as.double(estimate),
    ci_low = as.double(ci_low), ci_high = as.double(ci_high),
    ci_method = "eti", ci_level = 0.95,
    rhat = 1, ess_bulk = 1000, ess_tail = 1000,
    scale = scale, n = 4L, converged = TRUE,
    replication = replication
  )
}

# Subject truth for replications: id x term, long.
hand_subjects <- function(values, replication) {
  # values: a list term -> numeric vector over ids
  ids <- as.character(seq_along(values[[1L]]))
  dplyr::bind_rows(lapply(names(values), function(t) {
    tibble::tibble(
      id = ids, term = t, true_value = as.double(values[[t]]),
      replication = replication
    )
  }))
}

# A minimal bmmtools_simulation carrying only what the scorer reads.
fake_simulation <- function(subjects, cor, covariates = NULL,
                            links = list(kappa = "log", thetat = "logit")) {
  empty <- tibble::tibble(
    id = character(), term = character(), true_value = double()
  )
  structure(
    list(
      truth = list(
        subjects = subjects, cor = cor,
        covariates = covariates %||% empty
      ),
      model = list(links = links)
    ),
    class = "bmmtools_simulation"
  )
}

mock_fit <- function(parameters = c("kappa", "thetat"), ids = 1:6,
                     data = NULL) {
  structure(
    list(parameters = parameters, ids = as.character(ids), data = data),
    class = "mockfit"
  )
}

# extract_subject_draws() -------------------------------------------------

subject_level_draws <- function(n_subjects = 3L, seed = 3) {
  withr::local_seed(seed)
  n <- 80L
  values <- list(
    b_kappa_Intercept = stats::rnorm(n, 1),
    b_thetat_Intercept = stats::rnorm(n),
    b_fixed_Intercept = 0.5
  )
  for (i in seq_len(n_subjects)) {
    values[[sprintf("r_id__kappa[%d,Intercept]", i)]] <- stats::rnorm(n)
    values[[sprintf("r_id__thetat[%d,Intercept]", i)]] <- stats::rnorm(n)
    values[[sprintf("r_id__fixed[%d,Intercept]", i)]] <- 0
  }
  fake_draws(values)
}

test_that("subject draws are intercept plus deviation per draw, 4-d", {
  draws <- subject_level_draws()
  out <- subject_draws_from_draws(draws, "id")

  expect_equal(dim(out), c(40L, 2L, 3L, 2L))
  expect_named(dimnames(out), c("iteration", "chain", "id", "term"))
  expect_equal(dimnames(out)$id, c("1", "2", "3"))
  # the constant term is dropped as drop_constant_rows() would drop it
  expect_equal(dimnames(out)$term, c("kappa", "thetat"))
  expect_equal(attr(out, "group"), "id")

  plain <- unclass(draws)
  expect_equal(
    out[, , "2", "thetat"],
    plain[, , "b_thetat_Intercept"] + plain[, , "r_id__thetat[2,Intercept]"],
    ignore_attr = TRUE
  )
})

test_that("subject draws need a grouping factor with coefficients", {
  draws <- subject_level_draws()
  expect_error(subject_draws_from_draws(draws, character(0)), "no group")
  expect_error(
    subject_draws_from_draws(draws, c("id", "item")),
    "more than one grouping factor"
  )
  expect_error(
    subject_draws_from_draws(draws, c("id", "item"), group = "item"),
    "no group-level coefficients"
  )
})

test_that("the default method names the class it cannot read", {
  expect_error(extract_subject_draws(1:3), "integer")
  expect_error(extract_subject_draws(list(a = 1)), "brmsfit")
})

test_that("the mock method returns the documented shape", {
  out <- extract_subject_draws(mock_fit())
  expect_equal(dim(out), c(10L, 2L, 6L, 2L))
  expect_equal(dimnames(out)$term, c("kappa", "thetat"))
  expect_identical(out, extract_subject_draws(mock_fit()))
})

test_that("fixture: extract_subject_draws reads a real fit", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()
  out <- extract_subject_draws(fit)
  draws <- unclass(posterior::as_draws_array(fit))

  expect_equal(dim(out), c(500L, 2L, 8L, 2L))
  expect_equal(dimnames(out)$term, c("kappa", "thetat"))
  expect_equal(
    out[, , "5", "kappa"],
    draws[, , "b_kappa_Intercept"] + draws[, , "r_id__kappa[5,Intercept]"],
    ignore_attr = TRUE
  )
  expect_error(extract_subject_draws(fit, bogus = 1), "bogus")
})

# extract_correlations(): the draws and point estimators ------------------

test_that("draws is the median and ETI of the per-draw correlations", {
  x <- fake_subject_draws(n_subjects = 7L)
  out <- correlations_from_parts(x, estimator = "draws")

  expect_equal(nrow(out), 1L)
  expect_equal(out$term, "kappa__thetat")
  per_draw <- oracle_draw_cors(x, "kappa", "thetat")
  oracle <- oracle_summary(per_draw)
  expect_equal(out$estimate, oracle$estimate)
  expect_equal(out$ci_low, oracle$ci_low)
  expect_equal(out$ci_high, oracle$ci_high)
  # rhat is rank-normalised, so a 1e-16 difference between the two
  # correlation formulas can swap two ranks; hence the tolerance
  expect_equal(out$rhat, oracle$rhat, tolerance = 1e-3)
  expect_true(is.finite(out$ess_bulk) && is.finite(out$ess_tail))
  expect_equal(out$ci_method, "eti")
  expect_equal(out$n, 7L)
})

test_that("draws follows ci_level", {
  x <- fake_subject_draws()
  out <- correlations_from_parts(x, estimator = "draws", ci_level = 0.8)
  oracle <- oracle_summary(oracle_draw_cors(x, "kappa", "thetat"), 0.8)
  expect_equal(out$ci_low, oracle$ci_low)
  expect_equal(out$ci_level, 0.8)
})

test_that("point is cor.test() on the per-subject posterior means", {
  x <- fake_subject_draws(n_subjects = 9L)
  out <- correlations_from_parts(x, estimator = "point", ci_level = 0.9)

  means <- apply(x, c(3, 4), mean)
  test <- stats::cor.test(means[, "kappa"], means[, "thetat"],
    conf.level = 0.9
  )
  expect_equal(out$estimate, unname(test$estimate))
  expect_equal(c(out$ci_low, out$ci_high), as.double(test$conf.int))
  expect_equal(out$ci_method, "fisher_z")
  expect_true(is.na(out$rhat) && is.na(out$ess_bulk) && is.na(out$ess_tail))
  expect_type(out$rhat, "double")
})

test_that("the natural scale transforms each term before correlating", {
  x <- fake_subject_draws(n_subjects = 8L)
  links <- c(kappa = "log", thetat = "logit")
  out <- correlations_from_parts(
    x,
    estimator = c("draws", "point"), scale = "natural", links = links
  )
  expect_equal(out$scale, c("natural", "natural"))

  per_draw <- oracle_draw_cors(x, "kappa", "thetat", exp, stats::plogis)
  expect_equal(out$estimate[out$estimator == "draws"], stats::median(per_draw))

  # point: inverse link of the mean, not the mean of the inverse links
  means <- apply(x, c(3, 4), mean)
  expect_equal(
    out$estimate[out$estimator == "point"],
    stats::cor(exp(means[, "kappa"]), stats::plogis(means[, "thetat"]))
  )
})

test_that("the returned columns follow the spec order and types", {
  x <- fake_subject_draws()
  out <- correlations_from_parts(x, converged = TRUE)
  expect_named(out, c(
    "term", "var1", "var2", "estimator", "estimate", "ci_low", "ci_high",
    "ci_method", "ci_level", "rhat", "ess_bulk", "ess_tail", "scale", "n",
    "converged"
  ))
  expect_type(out$n, "integer")
  expect_type(out$converged, "logical")
  expect_equal(out$var1, c("kappa", "kappa"))
  expect_equal(out$var2, c("thetat", "thetat"))
  expect_true(all(out$converged))
})

test_that("every pair of terms is returned, and pairs selects in any order", {
  x <- fake_subject_draws(terms = c("thetat", "kappa", "c"))
  all_pairs <- correlations_from_parts(x, estimator = "point")
  expect_setequal(all_pairs$term, c("kappa__thetat", "c__thetat", "c__kappa"))

  one <- correlations_from_parts(
    x,
    estimator = "point", pairs = c("thetat__kappa", "kappa__thetat")
  )
  expect_equal(one$term, "kappa__thetat")

  expect_error(
    correlations_from_parts(x, pairs = "kappa__bogus"), "bogus"
  )
  expect_error(correlations_from_parts(x, pairs = "kappa"), "two different")
  expect_error(
    correlations_from_parts(x, pairs = "kappa__kappa"), "two different"
  )
  expect_error(correlations_from_parts(x, pairs = 1), "character")
})

test_that("NA, never 0, with fewer than three subjects or no spread", {
  two <- fake_subject_draws(n_subjects = 2L)
  out <- correlations_from_parts(two, estimator = c("draws", "point"))
  expect_true(all(is.na(out$estimate)))

  flat <- fake_subject_draws(n_subjects = 5L)
  covariates <- data.frame(id = as.character(1:5), G = 1)
  out <- correlations_from_parts(
    flat,
    estimator = c("draws", "point"), covariates = covariates,
    pairs = "G__kappa"
  )
  expect_true(all(is.na(out$estimate)))
  expect_false(any(out$estimate %in% 0))
})

# extract_correlations(): the model estimator ------------------------------

test_that("model rows are the 5.2 cor rows, whichever name order brms used", {
  withr::local_seed(5)
  values <- list(
    b_kappa_Intercept = stats::rnorm(80),
    b_thetat_Intercept = stats::rnorm(80),
    sd_id__kappa_Intercept = abs(stats::rnorm(80)),
    sd_id__thetat_Intercept = abs(stats::rnorm(80)),
    cor = stats::runif(80, -1, 1)
  )
  forward <- values
  names(forward)[[5L]] <- "cor_id__kappa_Intercept__thetat_Intercept"
  reversed <- values
  names(reversed)[[5L]] <- "cor_id__thetat_Intercept__kappa_Intercept"

  x <- fake_subject_draws()
  for (v in list(forward, reversed)) {
    cor_rows <- estimates_from_draws(
      fake_draws(v), "id",
      level = "cor", ranef = fake_ranef()
    )
    out <- correlations_from_parts(x, cor_rows, estimator = "model")
    expect_equal(out$term, "kappa__thetat")
    expect_equal(out$estimator, "model")
    expect_equal(out$estimate, cor_rows$estimate)
    expect_equal(out$ci_low, cor_rows$ci_low)
    expect_equal(out$rhat, cor_rows$rhat)
    expect_equal(out$ci_method, "eti")
    expect_equal(out$n, 6L)
  }
})

test_that("a hand-built cor row with reversed names is normalised", {
  cor_rows <- tibble::tibble(
    term = "thetat__kappa", estimate = 0.3, ci_low = 0.1, ci_high = 0.5,
    rhat = 1, ess_bulk = 500, ess_tail = 500
  )
  out <- correlations_from_parts(NULL, cor_rows, estimator = "model")
  expect_equal(out$term, "kappa__thetat")
  expect_true(is.na(out$n))
  expect_error(
    correlations_from_parts(NULL, cor_rows["term"], estimator = "model"),
    "estimate"
  )
})

test_that("the natural scale gives no model rows and says why", {
  x <- fake_subject_draws()
  cor_rows <- tibble::tibble(
    term = "kappa__thetat", estimate = 0.3, ci_low = 0.1, ci_high = 0.5,
    rhat = 1, ess_bulk = 500, ess_tail = 500
  )
  expect_message(
    out <- correlations_from_parts(
      x, cor_rows,
      scale = "natural", links = c(kappa = "log")
    ),
    "link scale"
  )
  expect_false("model" %in% out$estimator)
  expect_setequal(out$estimator, c("draws", "point"))

  expect_no_message(
    correlations_from_parts(x, cor_rows, estimator = "draws",
      scale = "natural", links = c(kappa = "log")
    )
  )
})

test_that("fixture: the model estimator on the correlated fixture", {
  record <- mixture2p_cor_draws()
  groups <- unique(as.character(record$ranef$group))
  cor_rows <- estimates_from_draws(
    record$draws, groups,
    level = "cor", ranef = record$ranef
  )
  out <- correlations_from_parts(
    NULL, cor_rows,
    estimator = "model", n_subjects = 20L
  )
  expect_equal(out$term, "kappa__thetat")
  expect_equal(out$estimate, cor_rows$estimate)
  expect_equal(out$ci_high, cor_rows$ci_high)
  expect_true(is.finite(out$rhat))
  expect_equal(out$n, 20L)
  expect_true(out$ci_low >= -1 && out$ci_high <= 1)
})

# covariates --------------------------------------------------------------

test_that("covariate pairs use the same vector in every draw, no model row", {
  x <- fake_subject_draws(n_subjects = 6L)
  withr::local_seed(8)
  g <- stats::rnorm(6)
  covariates <- data.frame(id = as.character(6:1), G = rev(g))
  cor_rows <- tibble::tibble(
    term = "kappa__thetat", estimate = 0.3, ci_low = 0.1, ci_high = 0.5,
    rhat = 1, ess_bulk = 500, ess_tail = 500
  )
  out <- correlations_from_parts(x, cor_rows, covariates = covariates)

  expect_setequal(
    out$term[out$estimator == "draws"],
    c("kappa__thetat", "G__kappa", "G__thetat")
  )
  expect_equal(out$term[out$estimator == "model"], "kappa__thetat")

  per_draw <- apply(x, 1:2, function(slice) stats::cor(g, slice[, "kappa"]))
  draws_row <- out[out$estimator == "draws" & out$term == "G__kappa", ]
  expect_equal(draws_row$estimate, stats::median(per_draw))
  expect_equal(draws_row$var1, "G")

  means <- apply(x, c(3, 4), mean)
  point_row <- out[out$estimator == "point" & out$term == "G__thetat", ]
  expect_equal(point_row$estimate, stats::cor(g, means[, "thetat"]))
})

test_that("a covariate repeated over trials is reduced to one per subject", {
  x <- fake_subject_draws(n_subjects = 4L)
  long <- data.frame(id = rep(1:4, each = 3), G = rep(c(2, 5, 1, 7), each = 3))
  wide <- data.frame(id = 1:4, G = c(2, 5, 1, 7))
  expect_equal(
    correlations_from_parts(x, covariates = long, estimator = "point"),
    correlations_from_parts(x, covariates = wide, estimator = "point")
  )
})

test_that("a covariate that varies within subject is an error", {
  x <- fake_subject_draws(n_subjects = 3L)
  bad <- data.frame(id = c(1, 1, 2, 3), G = c(1, 2, 3, 4))
  expect_error(
    correlations_from_parts(x, covariates = bad),
    "constant within"
  )
})

test_that("covariate tables are validated", {
  x <- fake_subject_draws(n_subjects = 3L)
  expect_error(
    correlations_from_parts(x, covariates = list(G = 1:3)), "data frame"
  )
  expect_error(
    correlations_from_parts(x, covariates = data.frame(subject = 1:3, G = 1)),
    "id"
  )
  expect_error(
    correlations_from_parts(x, covariates = data.frame(id = 1:3)),
    "no covariate"
  )
  expect_error(
    correlations_from_parts(x, covariates = data.frame(id = 1:3, G = "a")),
    "numeric"
  )
  expect_error(
    correlations_from_parts(x, covariates = data.frame(id = 1:2, G = 1:2)),
    "no value for"
  )
  expect_error(
    correlations_from_parts(
      x,
      covariates = data.frame(id = 1:3, G = c(1, NA, 3))
    ),
    "missing"
  )
  expect_error(
    correlations_from_parts(x, covariates = data.frame(id = 1:3, kappa = 1:3)),
    "subject term"
  )
  expect_error(
    correlations_from_parts(x, covariates = data.frame(id = 1:3, a__b = 1:3)),
    "__"
  )
})

test_that("the subject-draws array is validated", {
  expect_error(correlations_from_parts(matrix(1, 2, 2)), "4-d")
  bad <- fake_subject_draws()
  dimnames(bad) <- NULL
  expect_error(correlations_from_parts(bad), "dimnames")
})

# extract_correlations() on fit objects ----------------------------------

test_that("extract_correlations reads a mock fit through the generics", {
  data <- data.frame(
    id = rep(1:6, each = 2), G = rep(c(3, 1, 4, 1, 5, 9), each = 2)
  )
  fit <- mock_fit(data = data)
  expect_message(
    out <- extract_correlations(fit, covariates = "G", scale = "link"),
    NA
  )
  expect_setequal(out$estimator, c("draws", "point"))
  expect_setequal(
    unique(out$term), c("kappa__thetat", "G__kappa", "G__thetat")
  )
  # a mock fit carries no convergence verdict
  expect_true(all(is.na(out$converged)))

  expect_equal(
    out[out$estimator != "model", ],
    correlations_from_parts(
      extract_subject_draws(fit),
      covariates = data,
      estimator = c("draws", "point")
    )
  )
})

test_that("extract_correlations validates its arguments", {
  fit <- mock_fit()
  expect_error(extract_correlations(fit, estimator = "bogus"), "estimator")
  expect_error(extract_correlations(fit, scale = "bogus"), "scale")
  expect_error(extract_correlations(fit, ci_level = 2), "between 0 and 1")
  expect_error(extract_correlations(fit, ci_level = "a"), "single number")
  expect_error(extract_correlations(fit, converged = "yes"), "converged")
  expect_error(extract_correlations(fit, covariates = "G"), "no data")
  expect_error(
    extract_correlations(mock_fit(data = data.frame(id = 1:6)),
      covariates = "G"
    ),
    "G"
  )
  expect_error(extract_correlations(fit, bogus = 1), "bogus")
  expect_error(extract_correlations(1:3, estimator = "draws"), "integer")
})

test_that("natural scale without links falls back to the link scale", {
  expect_message(
    out <- extract_correlations(mock_fit(),
      estimator = "point", scale = "natural"
    ),
    "link scale"
  )
  expect_equal(out$scale, "link")
})

test_that("fixture: extract_correlations on a real fit", {
  skip_if_not_installed("brms")
  fit <- mixture2p_fit()

  out <- extract_correlations(fit)
  # the fixture has no cor_ draws, so there is no model row
  expect_setequal(out$estimator, c("draws", "point"))
  expect_equal(unique(out$term), "kappa__thetat")
  expect_equal(out$n, c(8L, 8L))
  expect_type(out$converged, "logical")
  expect_false(anyNA(out$converged))

  expect_message(
    natural <- extract_correlations(fit, scale = "natural"),
    "link scale"
  )
  expect_equal(unique(natural$scale), "natural")

  expect_error(extract_correlations(fit, covariates = "y"), "constant within")
  expect_equal(
    nrow(extract_correlations(fit, estimator = "model", converged = NA)), 0L
  )
})

# recover_correlations() --------------------------------------------------

three_reps <- function(scale = "link") {
  subjects <- dplyr::bind_rows(
    hand_subjects(list(kappa = c(1, 2, 3, 4), thetat = c(1, 3, 2, 5)), 1L),
    hand_subjects(list(kappa = c(4, 3, 2, 1), thetat = c(1, 2, 3, 5)), 2L),
    hand_subjects(list(kappa = c(1, 2, 4, 3), thetat = c(0, 0, 1, 1)), 3L)
  )
  cor <- tibble::tibble(
    term = "kappa__thetat", var1 = "kappa", var2 = "thetat", true_value = 0.5
  )
  fits <- hand_cor_rows(
    estimate = c(0.6, -0.2, 0.9),
    ci_low = c(0.2, -0.7, 0.5),
    ci_high = c(0.9, 0.3, 0.95),
    scale = scale
  )
  list(fits = fits, truth = list(cor = cor, subjects = subjects))
}

test_that("sample_value, bias, coverage and excludes_zero by hand", {
  input <- three_reps()
  out <- recover_correlations(input$fits, input$truth, scale = "link")

  expect_s3_class(out, "bmmtools_cor_recovery")
  expect_named(out, cor_recovery_contract_columns())
  sample <- c(
    stats::cor(c(1, 2, 3, 4), c(1, 3, 2, 5)),
    stats::cor(c(4, 3, 2, 1), c(1, 2, 3, 5)),
    stats::cor(c(1, 2, 4, 3), c(0, 0, 1, 1))
  )
  expect_equal(out$sample_value, sample)
  expect_equal(out$true_value, rep(0.5, 3))
  expect_equal(out$bias, c(0.1, -0.7, 0.4))
  expect_equal(out$bias_sample, c(0.6, -0.2, 0.9) - sample)
  expect_equal(out$covered, c(TRUE, FALSE, TRUE))
  expect_equal(
    out$covered_sample,
    sample >= c(0.2, -0.7, 0.5) & sample <= c(0.9, 0.3, 0.95)
  )
  expect_equal(out$excludes_zero, c(TRUE, FALSE, TRUE))
  expect_true(all(is.na(out$condition)))
})

test_that("a closed interval covers a bound equal to the value", {
  input <- three_reps()
  input$fits$ci_low[[1L]] <- 0.5
  input$fits$ci_high[[2L]] <- 0.5
  input$fits$ci_low[[3L]] <- 0.5
  out <- recover_correlations(input$fits, input$truth)
  expect_equal(out$covered, c(TRUE, TRUE, TRUE))
  expect_true(is.na(
    recover_correlations(
      dplyr::mutate(input$fits, ci_low = NA_real_), input$truth
    )$excludes_zero[[1L]]
  ))
})

test_that("natural scale: sample_value transformed, true_value 0 or NA", {
  input <- three_reps(scale = "natural")
  input$truth$cor <- tibble::tibble(
    term = c("kappa__thetat", "G__kappa"), true_value = c(0.5, 0)
  )
  input$truth$covariates <- tibble::tibble(
    id = rep(as.character(1:4), 3), term = "G",
    true_value = c(2, 1, 4, 3, 1, 1, 2, 2, 5, 4, 3, 2),
    replication = rep(1:3, each = 4)
  )
  fits <- dplyr::bind_rows(
    input$fits,
    hand_cor_rows(
      estimate = c(0.1, 0, -0.1), ci_low = -0.5, ci_high = 0.5,
      term = "G__kappa", scale = "natural"
    )
  )
  out <- recover_correlations(
    fits, input$truth,
    scale = "natural", links = c(kappa = "log", thetat = "logit")
  )
  pair <- out[out$term == "kappa__thetat", ]
  expect_true(all(is.na(pair$true_value)))
  expect_true(all(is.na(pair$covered)))
  expect_equal(
    pair$sample_value[[1L]],
    stats::cor(exp(c(1, 2, 3, 4)), stats::plogis(c(1, 3, 2, 5)))
  )
  cov_rows <- out[out$term == "G__kappa", ]
  expect_equal(cov_rows$true_value, c(0, 0, 0))
  expect_equal(
    cov_rows$sample_value[[2L]],
    stats::cor(c(1, 1, 2, 2), exp(c(4, 3, 2, 1)))
  )
  expect_equal(attr(out, "scale"), "natural")
})

test_that("a scale that disagrees with the extracted tibble is an error", {
  input <- three_reps()
  expect_error(
    recover_correlations(input$fits, input$truth,
      scale = "natural", links = c(kappa = "log")
    ),
    "extracted on"
  )
  # without scale, the tibble's own scale is used
  natural <- three_reps(scale = "natural")
  out <- recover_correlations(natural$fits, natural$truth,
    links = c(kappa = "log", thetat = "logit")
  )
  expect_equal(unique(out$scale), "natural")
})

test_that("truth without a replication column applies to every replication", {
  input <- three_reps()
  input$truth$subjects <- input$truth$subjects[
    input$truth$subjects$replication == 1L,
    c("id", "term", "true_value")
  ]
  out <- recover_correlations(input$fits, input$truth)
  first <- stats::cor(c(1, 2, 3, 4), c(1, 3, 2, 5))
  expect_equal(out$sample_value, rep(first, 3))
})

test_that("an unmatched truth pair warns and all unmatched is an error", {
  input <- three_reps()
  input$truth$cor <- tibble::tibble(
    term = c("kappa__thetat", "a__b"), true_value = c(0.5, 0.2)
  )
  expect_warning(
    out <- recover_correlations(input$fits, input$truth),
    "a__b"
  )
  expect_equal(unique(out$term), "kappa__thetat")

  input$truth$cor <- tibble::tibble(term = "a__b", true_value = 0.2)
  expect_error(recover_correlations(input$fits, input$truth), "matches")
})

test_that("pairs restricts both the estimates and the truth", {
  input <- three_reps()
  input$truth$cor <- tibble::tibble(
    term = c("kappa__thetat", "a__b"), true_value = c(0.5, 0.2)
  )
  expect_no_warning(
    out <- recover_correlations(input$fits, input$truth,
      pairs = "thetat__kappa"
    )
  )
  expect_equal(nrow(out), 3L)
})

test_that("estimator filters an extracted tibble", {
  input <- three_reps()
  fits <- dplyr::bind_rows(
    input$fits, dplyr::mutate(input$fits, estimator = "point")
  )
  out <- recover_correlations(fits, input$truth, estimator = "point")
  expect_equal(unique(out$estimator), "point")
  expect_error(
    recover_correlations(fits, input$truth, estimator = "model"),
    "No correlation estimates"
  )
})

test_that("recover_correlations scores one mock fit against a simulation", {
  fit <- mock_fit()
  withr::local_seed(4)
  subjects <- tibble::tibble(
    id = rep(as.character(1:6), 2), term = rep(c("kappa", "thetat"), each = 6),
    true_value = stats::rnorm(12)
  )
  covariates <- tibble::tibble(
    id = as.character(1:6), term = "G", true_value = stats::rnorm(6)
  )
  cor <- tibble::tibble(
    term = c("kappa__thetat", "G__kappa", "G__thetat"),
    var1 = c("kappa", "G", "G"), var2 = c("thetat", "kappa", "thetat"),
    true_value = c(0.4, 0, 0)
  )
  sim <- fake_simulation(subjects, cor, covariates)

  out <- recover_correlations(fit, sim, estimator = c("draws", "point"))
  expect_equal(nrow(out), 6L)
  expect_equal(unique(out$replication), 1L)

  # the covariate values come from the truth
  draws_g <- out[out$estimator == "draws" & out$term == "G__kappa", ]
  x <- extract_subject_draws(fit)
  per_draw <- apply(x, 1:2, function(s) {
    stats::cor(covariates$true_value, s[, "kappa"])
  })
  expect_equal(draws_g$estimate, stats::median(per_draw))

  # natural scale reads the links from the simulation's model
  expect_message(
    natural <- recover_correlations(fit, sim, scale = "natural"),
    "link scale"
  )
  expect_setequal(natural$estimator, c("draws", "point"))
  expect_equal(unique(natural$scale), "natural")
  expect_error(
    suppressMessages(
      recover_correlations(fit, sim, estimator = "model", scale = "natural")
    ),
    "nothing is left"
  )
})

test_that("a list of fits pairs with a list of simulations", {
  fits <- list(mock_fit(ids = 1:5), mock_fit(ids = 1:5))
  make_sim <- function(seed) {
    withr::local_seed(seed)
    fake_simulation(
      tibble::tibble(
        id = rep(as.character(1:5), 2),
        term = rep(c("kappa", "thetat"), each = 5),
        true_value = stats::rnorm(10)
      ),
      tibble::tibble(term = "kappa__thetat", true_value = 0)
    )
  }
  sims <- list(make_sim(1), make_sim(2))
  out <- recover_correlations(fits, sims, estimator = "point")
  expect_equal(out$replication, 1:2)
  wide <- function(sim) {
    s <- sim$truth$subjects
    stats::cor(
      s$true_value[s$term == "kappa"], s$true_value[s$term == "thetat"]
    )
  }
  expect_equal(out$sample_value, c(wide(sims[[1L]]), wide(sims[[2L]])))

  expect_error(
    recover_correlations(fits, sims[1L], estimator = "point"),
    "2 fits"
  )
  named <- recover_correlations(
    list(a = fits[[1L]], b = fits[[2L]]), sims,
    estimator = "point"
  )
  expect_equal(named$replication, c("a", "b"))
})

test_that("an extracted tibble without replication or converged is one fit", {
  input <- three_reps()
  keep <- setdiff(names(input$fits), c("replication", "converged"))
  one <- input$fits[1L, keep]
  truth <- input$truth
  truth$subjects <- truth$subjects[truth$subjects$replication == 1L, ]
  out <- recover_correlations(one, truth)
  expect_equal(out$replication, 1L)
  expect_true(is.na(out$converged))
})

test_that("a truth split by condition needs estimates that carry it", {
  input <- three_reps()
  truth <- input$truth
  truth$cor$condition <- "a"
  expect_false("condition" %in% names(input$fits))
  expect_error(
    recover_correlations(input$fits, truth),
    "condition"
  )
  keyed <- input$fits
  keyed$condition <- "a"
  out <- recover_correlations(keyed, truth)
  expect_true(all(out$condition == "a"))
})

test_that("an extracted tibble pairs with a list of simulations by label", {
  input <- three_reps()
  sims <- lapply(1:3, function(rep) {
    fake_simulation(
      input$truth$subjects[input$truth$subjects$replication == rep, 1:3],
      input$truth$cor
    )
  })
  out <- recover_correlations(input$fits, sims)
  expect_equal(
    out$sample_value,
    recover_correlations(input$fits, input$truth)$sample_value
  )
})

test_that("a model-only table with covariates keeps the covariate pairs", {
  cor_rows <- tibble::tibble(
    term = "kappa__thetat", estimate = 0.3, ci_low = 0.1, ci_high = 0.5,
    rhat = 1, ess_bulk = 500, ess_tail = 500
  )
  out <- correlations_from_parts(
    NULL, cor_rows,
    covariates = data.frame(id = 1:3, G = 1:3),
    estimator = c("model", "draws", "point"), pairs = "G__kappa"
  )
  expect_equal(nrow(out), 0L)
  expect_named(out, correlation_columns())
})

test_that("covariates in a list of simulations follow their replication", {
  fits <- list(mock_fit(ids = 1:4), mock_fit(ids = 1:4))
  make_sim <- function(g) {
    fake_simulation(
      tibble::tibble(
        id = rep(as.character(1:4), 2),
        term = rep(c("kappa", "thetat"), each = 4),
        true_value = c(1, 2, 3, 4, 2, 1, 4, 3)
      ),
      tibble::tibble(term = "G__kappa", true_value = 0),
      tibble::tibble(id = as.character(1:4), term = "G", true_value = g)
    )
  }
  g1 <- c(1, 3, 2, 4)
  g2 <- c(4, 3, 2, 1)
  out <- recover_correlations(
    fits, list(make_sim(g1), make_sim(g2)),
    estimator = "point"
  )
  x <- colMeans(extract_subject_draws(fits[[1L]]), dims = 2L)
  expect_equal(
    out$estimate,
    c(stats::cor(g1, x[, "kappa"]), stats::cor(g2, x[, "kappa"]))
  )
  expect_equal(
    out$sample_value,
    c(stats::cor(g1, 1:4), stats::cor(g2, 1:4))
  )
})

test_that("a truth table that is not a data frame is an error", {
  input <- three_reps()
  truth <- input$truth
  truth$covariates <- 1
  expect_error(
    recover_correlations(input$fits, truth),
    "must be a data frame"
  )
})

test_that("recover_correlations validates fits and truth", {
  input <- three_reps()
  expect_error(recover_correlations("x", input$truth), "fits")
  expect_error(
    recover_correlations(input$fits["term"], input$truth),
    "missing"
  )
  expect_error(recover_correlations(input$fits, 1), "truth")
  expect_error(
    recover_correlations(input$fits, list(cor = input$truth$cor)),
    "truth"
  )
  expect_error(
    recover_correlations(
      input$fits,
      list(cor = input$truth$cor["term"], subjects = input$truth$subjects)
    ),
    "true_value"
  )
  old <- structure(
    list(truth = list(population = tibble::tibble())),
    class = "bmmtools_simulation"
  )
  expect_error(recover_correlations(input$fits, old), "correlation")
})

# the class ---------------------------------------------------------------

test_that("the constructor checks the contract", {
  input <- three_reps()
  x <- tibble::as_tibble(recover_correlations(input$fits, input$truth))
  incomplete <- x[setdiff(names(x), "sample_value")]
  expect_error(
    new_bmmtools_cor_recovery(incomplete, "link", 0.95),
    "sample_value"
  )
  x$covered_sample <- as.character(x$covered_sample)
  expect_error(new_bmmtools_cor_recovery(x, "link", 0.95), "covered_sample")
  expect_error(new_bmmtools_cor_recovery(1, "link", 0.95), "data frame")
})

test_that("select() demotes the class and filter() keeps it", {
  input <- three_reps()
  out <- recover_correlations(input$fits, input$truth)

  reduced <- dplyr::select(out, "term", "estimate")
  expect_false(inherits(reduced, "bmmtools_cor_recovery"))
  expect_false(inherits(out[, 1:3], "bmmtools_cor_recovery"))

  kept <- dplyr::filter(out, replication > 1)
  expect_s3_class(kept, "bmmtools_cor_recovery")
  expect_equal(nrow(kept), 2L)

  expect_false(inherits(
    dplyr_reconstruct(tibble::as_tibble(out)["term"], out),
    "bmmtools_cor_recovery"
  ))
  expect_s3_class(
    dplyr_reconstruct(tibble::as_tibble(out), out), "bmmtools_cor_recovery"
  )

  broken <- structure(
    tibble::as_tibble(out)["term"],
    class = class(out)
  )
  expect_error(summary(broken), "missing")
})

test_that("format leads with the scale and estimators, print is invisible", {
  input <- three_reps()
  out <- recover_correlations(input$fits, input$truth)
  text <- format(out)
  expect_equal(text[[1L]], "<bmmtools_cor_recovery>")
  expect_match(text[[2L]], "link scale")
  expect_match(text[[2L]], "draws")
  expect_true(any(grepl("kappa__thetat", text, fixed = TRUE)))
  expect_output(res <- print(out), "cor_recovery")
  expect_identical(res, out)

  empty <- out[0, ]
  expect_match(format(empty)[[2L]], "No correlations")
  expect_equal(nrow(summary(empty)), 0L)
  expect_named(summary(empty), cor_recovery_summary_columns())

  two <- out[1:2, ]
  expect_true(any(grepl("not estimable", format(two))))
})

test_that("format notes the natural-scale reading of true_value", {
  natural <- three_reps(scale = "natural")
  out <- recover_correlations(natural$fits, natural$truth,
    links = c(kappa = "log", thetat = "logit")
  )
  expect_true(any(grepl("sample_value", format(out), fixed = TRUE)))
})

# summary -----------------------------------------------------------------

test_that("summary has one row per term x estimator x scale, spec columns", {
  input <- three_reps()
  fits <- dplyr::bind_rows(
    input$fits, dplyr::mutate(input$fits, estimator = "point")
  )
  out <- recover_correlations(fits, input$truth)
  s <- summary(out)
  expect_s3_class(s, "bmmtools_cor_recovery_summary")
  expect_named(s, cor_recovery_summary_columns())
  expect_equal(nrow(s), 2L)
  expect_equal(s$estimator, c("draws", "point"))

  row <- s[1L, ]
  rows <- out[out$estimator == "draws", ]
  expect_equal(row$n_replications, 3L)
  expect_equal(row$n_converged, 3L)
  expect_equal(row$true_value, 0.5)
  expect_equal(row$sample_sd, stats::sd(rows$sample_value))
  expect_equal(row$mean_estimate, mean(rows$estimate))
  expect_equal(row$bias, mean(rows$estimate - 0.5))
  expect_equal(row$rmse, sqrt(mean((rows$estimate - 0.5)^2)))
  expect_equal(row$bias_sample, mean(rows$estimate - rows$sample_value))
  expect_equal(
    row$rmse_sample, sqrt(mean((rows$estimate - rows$sample_value)^2))
  )
  expect_equal(row$coverage, mean(rows$covered))
  expect_equal(row$coverage_sample, mean(rows$covered_sample))
  expect_equal(row$ci_width, mean(rows$ci_high - rows$ci_low))
  expect_equal(row$rejection_rate, 2 / 3)
  expect_true(is.na(row$false_positive_rate))
  expect_equal(row$power, 2 / 3)

  r <- metric_r(rows$estimate, rows$sample_value)
  ccc <- recovery_ccc(rows$estimate, rows$sample_value)
  expect_equal(c(row$r, row$r_low, row$r_high), c(r$r, r$r_low, r$r_high))
  expect_equal(
    c(row$ccc, row$ccc_low, row$ccc_high),
    c(ccc$ccc, ccc$ccc_low, ccc$ccc_high)
  )
  expect_identical(summary(s), s)
  expect_output(print(s), "rejection_rate")
})

test_that("false_positive_rate needs every true_value 0, power none", {
  make <- function(true_values, excludes) {
    n <- length(true_values)
    x <- hand_cor_rows(
      estimate = rep(0.1, n),
      ci_low = ifelse(excludes, 0.05, -0.2), ci_high = 0.4
    )
    x$true_value <- as.double(true_values)
    x$sample_value <- seq_len(n) / 10
    x$bias <- x$estimate - x$true_value
    x$bias_sample <- x$estimate - x$sample_value
    x$covered <- x$covered_sample <- NA
    x$excludes_zero <- excludes
    x$converged <- NA
    new_bmmtools_cor_recovery(x, "link", 0.95)
  }
  null <- summary(make(c(0, 0, 0, 0), c(TRUE, FALSE, FALSE, FALSE)))
  expect_equal(null$false_positive_rate, 0.25)
  expect_true(is.na(null$power))
  expect_equal(null$true_value, 0)
  expect_true(is.na(null$n_converged))

  alt <- summary(make(c(0.3, 0.3, 0.3), c(TRUE, TRUE, FALSE)))
  expect_equal(alt$power, 2 / 3)
  expect_true(is.na(alt$false_positive_rate))

  # NA true_value (a nonzero correlation on the natural scale) is not 0
  natural <- summary(make(c(NA, NA), c(TRUE, TRUE)))
  expect_equal(natural$power, 1)
  expect_true(is.na(natural$false_positive_rate))
  expect_true(is.na(natural$true_value))

  mixed <- summary(make(c(0, 0.3, 0.2), c(TRUE, TRUE, FALSE)))
  expect_true(is.na(mixed$power) && is.na(mixed$false_positive_rate))
  expect_true(is.na(mixed$true_value))
  expect_equal(mixed$rejection_rate, 2 / 3)
})

test_that("summary splits by condition when there is one", {
  input <- three_reps()
  fits <- dplyr::bind_rows(
    dplyr::mutate(input$fits, condition = "a"),
    dplyr::mutate(input$fits, condition = "b")
  )
  truth <- input$truth
  truth$cor <- dplyr::bind_rows(
    dplyr::mutate(truth$cor, condition = "a"),
    dplyr::mutate(truth$cor, condition = "b", true_value = 0)
  )
  out <- recover_correlations(fits, truth)
  s <- summary(out)
  expect_named(s, c("condition", cor_recovery_summary_columns()))
  expect_equal(s$condition, c("a", "b"))
  expect_equal(s$true_value, c(0.5, 0))
})

# plot_recovery() -------------------------------------------------------

test_that("plot_recovery dispatches on both recovery classes", {
  skip_if_not_installed("ggplot2")
  input <- three_reps()
  cor_out <- recover_correlations(input$fits, input$truth)
  p <- plot_recovery(cor_out)
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  layers <- vapply(p$layers, function(l) class(l$geom)[[1L]], character(1))
  point <- which(layers == "GeomPoint")
  expect_setequal(built$data[[point]]$x, cor_out$sample_value)
  expect_setequal(built$data[[point]]$y, cor_out$estimate)
  expect_match(p$labels$x, "In-sample")

  true_plot <- plot_recovery(cor_out, truth = "true", intervals = FALSE,
    identity_line = FALSE, color_by = NULL, facet_by = NULL
  )
  expect_match(true_plot$labels$x, "Generating")
  true_layers <- vapply(
    true_plot$layers, function(l) class(l$geom)[[1L]], character(1)
  )
  expect_false("GeomLinerange" %in% true_layers)

  estimates <- fake_estimates(c("a", "b", "c"), estimate = c(1, 2, 3))
  truth <- fake_truth(c("a", "b", "c"), true_value = c(1, 2, 3))
  expect_s3_class(
    plot_recovery(recover(estimates, truth, scale = "link")), "ggplot"
  )

  expect_error(plot_recovery(cor_out, truth = "bogus"), "truth")
  expect_error(plot_recovery(cor_out, facet_by = "bogus"), "bogus")
  expect_error(plot_recovery(cor_out, bogus = 1), "bogus")
})

test_that("rows without the chosen truth are dropped with a message", {
  skip_if_not_installed("ggplot2")
  natural <- three_reps(scale = "natural")
  out <- recover_correlations(natural$fits, natural$truth,
    links = c(kappa = "log", thetat = "logit")
  )
  expect_error(plot_recovery(out, truth = "true"), "sample")

  mixed <- out
  mixed$true_value[[1L]] <- 0
  expect_message(p <- plot_recovery(mixed, truth = "true"), "2 rows")
  expect_s3_class(p, "ggplot")
})

# the trade-off property (spec 5.3; local/dev/sim/cor_estimators.R) --------

test_that("at rho = 0 with a posterior trade-off, point is biased, draws not", {
  skip_on_cran()
  withr::local_seed(20260914)
  # the cor_estimators.R cell rho = 0, rel = .5, c_err = -.5: known
  # hyperparameters, so the posterior is exact
  n <- 50L
  n_iter <- 100L
  n_chain <- 2L
  rel <- 0.5
  e2 <- (1 - rel) / rel
  psi <- e2 * matrix(c(1, -0.5, -0.5, 1), 2)
  v <- solve(diag(2) + solve(psi))
  w <- v %*% solve(psi)
  upper <- chol(v)
  n_draws <- n_iter * n_chain

  reps <- replicate(200L, {
    theta <- matrix(stats::rnorm(2 * n), n)
    y <- theta + matrix(stats::rnorm(2 * n), n) %*% chol(psi)
    post_mean <- y %*% t(w)
    noise <- matrix(stats::rnorm(n_draws * n * 2), ncol = 2) %*% upper
    x <- array(
      NA_real_,
      dim = c(n_iter, n_chain, n, 2),
      dimnames = list(
        iteration = NULL, chain = NULL,
        id = as.character(seq_len(n)), term = c("a", "b")
      )
    )
    for (k in 1:2) {
      x[, , , k] <- array(
        rep(post_mean[, k], each = n_draws) + noise[, k],
        dim = c(n_iter, n_chain, n)
      )
    }
    out <- correlations_from_parts(x, estimator = c("draws", "point"))
    stats::setNames(out$estimate, out$estimator)
  })

  expect_gt(abs(mean(reps["point", ])), 0.15)
  expect_lt(abs(mean(reps["draws", ])), 0.05)
})
