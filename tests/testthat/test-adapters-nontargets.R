# Tests for the mixture3p and imm adapters (Milestone 21.2, spec
# local/dev/spec-milestone-21-design-and-adapters.md § 9, D56). Both models
# read set size and non-target locations (and, for imm, distances) from a
# trial design. Every expected weight below is written from the brms
# formula bmm builds (configure_model.mixture3p(), configure_model.imm_*()),
# measured in the generated Stan code on bmm 1.3.2.9000 (2026-09-30 build):
# the guessing component has weight 0 on the log scale, and mixture3p's
# lures share exp(thetant) equally. Nothing here compiles Stan.

mock_bmm <- function(formula, data, model) {
  bmm::bmm(
    formula, data, model,
    backend = "mock", mock_fit = 1, rename = FALSE, silent = 2
  )
}

m3p_model <- function(set_size = "ss") {
  bmm::mixture3p(
    resp_error = "y", nt_features = paste0("nt", 1:3), set_size = set_size
  )
}

# bmm::imm() reads its own call to find missing arguments, so it cannot
# be reached through do.call()
imm_model <- function(version, set_size = "ss") {
  if (version == "abc") {
    return(bmm::imm(
      resp_error = "y", nt_features = paste0("nt", 1:3),
      set_size = set_size, version = "abc"
    ))
  }
  bmm::imm(
    resp_error = "y", nt_features = paste0("nt", 1:3),
    nt_distances = paste0("d", 1:3), set_size = set_size, version = version
  )
}

imm_pars <- list(
  full = c(kappa = log(8), a = log(0.5), c = log(2), s = log(1.5)),
  bsc = c(kappa = log(8), c = log(2), s = log(1.5)),
  abc = c(kappa = log(8), a = log(0.5), c = log(2))
)

#' Random set sizes 1 to 4, lures NA-padded beyond set size - 1
random_design <- function(distances = TRUE) {
  function(n_trials) {
    ss <- sample(1:4, n_trials, replace = TRUE)
    out <- data.frame(ss = ss)
    for (j in 1:3) {
      present <- ss >= j + 1
      locations <- stats::runif(n_trials, -pi, pi)
      out[[paste0("nt", j)]] <- ifelse(present, locations, NA)
      if (distances) {
        spread <- stats::runif(n_trials, 0.2, 2)
        out[[paste0("d", j)]] <- ifelse(present, spread, NA)
      }
    }
    out
  }
}

#' A fixed trial list: one set size, lures at fixed places
fixed_design <- function(n_trials, ss, nt = c(2, -2, 1), d = c(0.5, 1, 1.5)) {
  out <- data.frame(ss = rep(ss, n_trials))
  for (j in 1:3) {
    out[[paste0("nt", j)]] <- if (ss >= j + 1) nt[[j]] else NA_real_
    out[[paste0("d", j)]] <- if (ss >= j + 1) d[[j]] else NA_real_
  }
  out
}

#' Share of responses within `width` of `at`, on the circle
near <- function(y, at, width = 0.3) {
  mean(abs(atan2(sin(y - at), cos(y - at))) < width)
}

#' The share near `at` is `expected`, within 0.035 (about 2.7 SD of a
#' share at 1500 trials); a per-item lure weight would be 0.1 away
expect_share <- function(y, at, expected, label = NULL) {
  expect_lt(abs(near(y, at) - expected), 0.035, label = label)
}

# the tables --------------------------------------------------------------

test_that("trial_design_columns reads set size, lures, distances", {
  skip_if_not_installed("bmm")
  expect_identical(
    trial_design_columns(m3p_model()), c("ss", paste0("nt", 1:3))
  )
  # a number is the set size of every trial, not a column
  expect_identical(trial_design_columns(m3p_model(4)), paste0("nt", 1:3))
  for (v in c("full", "bsc")) {
    expect_identical(
      trial_design_columns(imm_model(v)),
      c("ss", paste0("nt", 1:3), paste0("d", 1:3)),
      label = v
    )
  }
  expect_identical(
    trial_design_columns(imm_model("abc")), c("ss", paste0("nt", 1:3))
  )
  # the seven 0.2.0 adapters still read none
  expect_identical(
    trial_design_columns(bmm::mixture2p(resp_error = "y")), character()
  )
})

test_that("generator_for knows mixture3p and every imm version", {
  skip_if_not_installed("bmm")
  expect_type(generator_for(m3p_model()), "closure")
  for (v in c("full", "bsc", "abc")) {
    expect_type(generator_for(imm_model(v)), "closure")
  }
  expect_true(all(c("mixture3p", "imm") %in% adapter_classes()))
})

test_that("fit_ml's optim route has no density for them", {
  skip_if_not_installed("bmm")
  expect_null(density_for(m3p_model()))
  expect_null(density_for(imm_model("full")))
  expect_false(any(c("mixture3p", "imm") %in% density_classes()))
  expect_setequal(
    density_classes(),
    c("sdt_yn", "sdt_mafc", "ezdm", "ddm", "cswald", "mixture2p", "sdm")
  )
})

# the softmax link (D56) --------------------------------------------------

test_that("a softmax term reaches the generator on the link scale", {
  skip_if_not_installed("bmm")
  seen <- new.env()
  recorder <- function(pars, n_trials, model, trial_design) {
    seen$pars <- pars
    data.frame(y = rep(0, n_trials))
  }
  simulate_recovery(
    m3p_model(), c(kappa = log(8), thetat = 1.5, thetant = -0.5),
    n_subjects = 1, n_trials = 4,
    trial_design = fixed_design(4, ss = 4)[c("ss", paste0("nt", 1:3))],
    generator = recorder
  )
  expect_equal(seen$pars$kappa, 8)
  expect_identical(seen$pars$thetat, 1.5)
  expect_identical(seen$pars$thetant, -0.5)
  expect_equal(seen$pars$mu1, 0)
  # the exported transform stays elementwise and keeps refusing it
  expect_error(inverse_link(1, "softmax"), "softmax")
})

test_that("mixture3p trial weights follow bmm's likelihood", {
  skip_if_not_installed("bmm")
  # set size >= 2: softmax over target, lures (shared) and guessing at 0
  w <- mixture3p_weights(thetat = log(2), thetant = 0, set_size = c(2, 3, 4))
  expect_equal(w$p_mem, rep(0.5, 3))
  expect_equal(w$p_nt, rep(0.25, 3))
  # set size 1: no lure, so no lure weight, and the target renormalised
  w1 <- mixture3p_weights(thetat = log(2), thetant = 0, set_size = 1)
  expect_equal(w1$p_mem, 2 / 3)
  expect_identical(w1$p_nt, 0)
  # extreme values neither overflow nor sum past one
  big <- mixture3p_weights(thetat = 800, thetant = 790, set_size = 3)
  expect_true(is.finite(big$p_mem) && big$p_mem + big$p_nt <= 1)
})

test_that("the mixture3p adapter's data has the shares the parameters imply", {
  skip_if_not_installed("bmm")
  # p_mem = .5, p_nt = .25 in total (not per lure), p_guess = .25;
  # a per-item lure weight would give p_mem = .4 at set size 3
  pars <- c(kappa = log(400), thetat = log(2), thetant = 0)
  n <- 1500
  sim <- simulate_recovery(
    m3p_model(), pars,
    n_subjects = 1, n_trials = n,
    trial_design = fixed_design(n, ss = 3)[c("ss", paste0("nt", 1:3))],
    seed = 21
  )
  guess <- 0.25 * 0.6 / (2 * pi)
  expect_share(sim$data$y, 0, 0.5 + guess)
  expect_share(sim$data$y, 2, 0.125 + guess)
  expect_share(sim$data$y, -2, 0.125 + guess)
  expect_share(sim$data$y, 1, guess)

  # set size 1: p_mem = 2/3, the rest guessing, nothing at a lure column
  one <- simulate_recovery(
    m3p_model(), pars,
    n_subjects = 1, n_trials = n,
    trial_design = fixed_design(n, ss = 1)[c("ss", paste0("nt", 1:3))],
    seed = 21
  )
  expect_true(all(is.na(one$data$nt1)))
  expect_share(one$data$y, 0, 2 / 3 + (1 / 3) * 0.6 / (2 * pi))
})

test_that("imm adapters give each version's weights from bmm's likelihood", {
  skip_if_not_installed("bmm")
  # natural c = 2, a = 0.5, s = 1.5, b = 1 (guessing at 0 on the log
  # scale); lure 1 at 2 with distance 0.5, lure 2 at -2 with distance 1
  n <- 1500
  design <- fixed_design(n, ss = 3)
  expected <- list(
    # target c + a, lure c * exp(-s d) + a
    full = c(2.5, 2 * exp(-0.75) + 0.5, 2 * exp(-1.5) + 0.5, 1),
    # target c, lure c * exp(-s d)
    bsc = c(2, 2 * exp(-0.75), 2 * exp(-1.5), 1),
    # target c + a, lure a
    abc = c(2.5, 0.5, 0.5, 1)
  )
  for (v in names(expected)) {
    keep <- c("ss", paste0("nt", 1:3), if (v != "abc") paste0("d", 1:3))
    sim <- simulate_recovery(
      imm_model(v), replace(imm_pars[[v]], "kappa", log(400)),
      n_subjects = 1, n_trials = n, trial_design = design[keep], seed = 7
    )
    p <- expected[[v]] / sum(expected[[v]])
    guess <- p[[4]] * 0.6 / (2 * pi)
    expect_share(sim$data$y, 0, p[[1]] + guess, label = v)
    expect_share(sim$data$y, 2, p[[2]] + guess, label = v)
    expect_share(sim$data$y, -2, p[[3]] + guess, label = v)
  }
})

test_that("the imm adapter maps each version onto rimm()'s arguments", {
  skip_if_not_installed("bmm")
  pars <- list(mu1 = 0, kappa = 8, a = 0.5, c = 2, s = 1.5)
  full <- imm_args(pars, imm_model("full"), lures = c(2, -2), dist = c(0.5, 1))
  expect_equal(full$mu, c(0, 2, -2))
  expect_equal(full$dist, c(0, 0.5, 1))
  expect_equal(
    full[c("c", "a", "s", "b", "kappa")],
    list(c = 2, a = 0.5, s = 1.5, b = 1, kappa = 8)
  )
  bsc <- imm_args(pars[c("mu1", "kappa", "c", "s")], imm_model("bsc"),
    lures = 2, dist = 0.5
  )
  expect_identical(bsc$a, 0)
  abc <- imm_args(pars[c("mu1", "kappa", "a", "c")], imm_model("abc"),
    lures = c(2, -2), dist = NULL
  )
  # no distance term: lures carry `a` alone
  expect_identical(abc$dist, c(0, Inf, Inf))
  # set size 1: the target alone
  one <- imm_args(pars, imm_model("full"), lures = numeric(), dist = numeric())
  expect_identical(one$mu, 0)
  expect_identical(one$dist, 0)
})

test_that("the location parameter is read off the model, mu or mu1", {
  model <- structure(
    list(parameters = list(mu = "", kappa = "")),
    class = c("bmmodel", "mixture3p")
  )
  expect_identical(location_parameter(model), "mu")
  model$parameters <- list(mu1 = "", kappa = "")
  expect_identical(location_parameter(model), "mu1")
  # a model listing neither keeps bmm's present name
  model$parameters <- list(kappa = "")
  expect_identical(location_parameter(model), "mu1")
})

# data bmm accepts --------------------------------------------------------

test_that("adapter data with random designs pass bmm's checks (mock backend)", {
  skip_if_not_installed("bmm")
  skip_if_not_installed("brms")
  cases <- list(
    list(
      model = m3p_model(),
      pars = c(kappa = log(8), thetat = 1, thetant = -0.5)
    ),
    list(model = imm_model("full"), pars = imm_pars$full),
    list(model = imm_model("bsc"), pars = imm_pars$bsc),
    list(model = imm_model("abc"), pars = imm_pars$abc)
  )
  for (case in cases) {
    label <- class(case$model)[[length(class(case$model))]]
    distances <- !is.null(case$model$other_vars$nt_distances)
    sim <- simulate_recovery(
      case$model, case$pars,
      n_subjects = 3, n_trials = 12,
      trial_design = random_design(distances), seed = 3
    )
    expect_identical(
      sim$generator, paste0("adapter:", adapter_name(case$model))
    )
    expect_setequal(sim$trial_design_columns, trial_design_columns(case$model))
    # the truth names exactly bmm's free parameters
    p <- bmm_parameter_info(case$model)
    expect_setequal(sim$truth$population$term, p$parameter[!p$fixed])
    expect_true(all(abs(sim$data$y) <= pi), label = label)
    fit <- mock_bmm(recovery_formula(case$model), sim$data, case$model)
    expect_s3_class(fit, "bmmfit")
  }
})

test_that("a numeric set size needs only the lure columns", {
  skip_if_not_installed("bmm")
  skip_if_not_installed("brms")
  model <- m3p_model(4)
  sim <- simulate_recovery(
    model, c(kappa = log(8), thetat = 1, thetant = -0.5),
    n_subjects = 2, n_trials = 6,
    trial_design = fixed_design(6, ss = 4)[paste0("nt", 1:3)], seed = 1
  )
  expect_named(sim$data, c("id", paste0("nt", 1:3), "y"))
  expect_s3_class(
    mock_bmm(recovery_formula(model), sim$data, model), "bmmfit"
  )
})

test_that("a set-size-1 trial list simulates and passes bmm's checks", {
  skip_if_not_installed("bmm")
  skip_if_not_installed("brms")
  design <- rbind(fixed_design(3, ss = 1), fixed_design(3, ss = 4))
  for (model in list(m3p_model(), imm_model("full"))) {
    pars <- if (inherits(model, "imm")) {
      imm_pars$full
    } else {
      c(kappa = log(8), thetat = 1, thetant = -0.5)
    }
    keep <- intersect(names(design), trial_design_columns(model))
    sim <- simulate_recovery(
      model, pars,
      n_subjects = 2, n_trials = 6, trial_design = design[keep], seed = 2
    )
    expect_equal(sim$data$ss, rep(design$ss, 2))
    expect_s3_class(
      mock_bmm(recovery_formula(model), sim$data, model), "bmmfit"
    )
  }
})

# refusals ----------------------------------------------------------------

test_that("the adapters need a design, with the columns the model names", {
  skip_if_not_installed("bmm")
  pars <- c(kappa = log(8), thetat = 1, thetant = -0.5)
  expect_error(
    simulate_recovery(m3p_model(), pars, n_subjects = 1, n_trials = 3),
    "needs `trial_design`"
  )
  design <- fixed_design(3, ss = 4)
  expect_error(
    simulate_recovery(
      m3p_model(), pars,
      n_subjects = 1, n_trials = 3, trial_design = design[c("ss", "nt1", "nt2")]
    ),
    "\"nt3\""
  )
  expect_error(
    simulate_recovery(
      imm_model("full"), imm_pars$full,
      n_subjects = 1, n_trials = 3,
      trial_design = design[c("ss", paste0("nt", 1:3))]
    ),
    "\"d1\""
  )
  # a function's result is checked per call
  expect_error(
    simulate_recovery(
      m3p_model(), pars,
      n_subjects = 1, n_trials = 3,
      trial_design = function(n) design[seq_len(n), c("ss", "nt1")]
    ),
    "\"nt2\""
  )
})

test_that("a lure the set size includes must have a location", {
  skip_if_not_installed("bmm")
  design <- fixed_design(3, ss = 3)[c("ss", paste0("nt", 1:3))]
  design$nt2[[2]] <- NA
  expect_error(
    simulate_recovery(
      m3p_model(), c(kappa = log(8), thetat = 1, thetant = -0.5),
      n_subjects = 1, n_trials = 3, trial_design = design, seed = 1
    ),
    "nt2"
  )
  # a set size beyond the lure columns cannot be laid out
  design <- fixed_design(3, ss = 4)[c("ss", paste0("nt", 1:3))]
  design$ss[[1]] <- 5
  expect_error(
    simulate_recovery(
      m3p_model(), c(kappa = log(8), thetat = 1, thetant = -0.5),
      n_subjects = 1, n_trials = 3, trial_design = design, seed = 1
    ),
    "set size"
  )
})

# prior_check -------------------------------------------------------------

test_that("prior_check knows the circular range of both models", {
  skip_if_not_installed("bmm")
  for (model in list(m3p_model(), imm_model("abc"))) {
    expect_identical(
      response_range(model, data.frame(y = 0)),
      list(floor = -pi, ceiling = pi)
    )
  }
})

# sbc ---------------------------------------------------------------------

test_that("sbc()'s default generator simulates mixture3p on data's trials", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran() # a whole sbc() run, as in every sbc() test
  model <- m3p_model()
  data <- withr::with_seed(4, {
    d <- random_design(FALSE)(12)
    data.frame(id = rep(c("b", "a", "c"), times = 4), d, y = 0)
  })
  k <- seq_len(60)
  draws <- posterior::as_draws_matrix(cbind(
    b_kappa_Intercept = log(6) + 0.01 * k,
    b_thetat_Intercept = 1 + 0.01 * k,
    b_thetant_Intercept = -0.5 + 0.01 * k
  ))
  mock <- sbc_mock_fitter(draws = draws)
  # SBC's progress and the mock's diagnostics, as in sbc_run()
  suppressMessages(withCallingHandlers(
    sbc(
      model, bmm::bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1), data,
      n_sims = 2L, thin_ranks = 1, cores_per_fit = 1, level = "population",
      .fitter = mock$fitter
    ),
    bmmtools_sbc_diagnostics = function(w) invokeRestart("muffleWarning")
  ))
  calls <- sbc_dataset_calls(mock$calls)
  expect_length(calls, 2L)
  for (call in calls) {
    generated <- call$data
    for (label in c("a", "b", "c")) {
      mine <- as.character(generated$id) == label
      for (col in trial_design_columns(model)) {
        expect_identical(generated[[col]][mine], data[[col]][data$id == label])
      }
    }
    expect_true(all(abs(generated$y) <= pi))
  }
})

test_that("the NA message says what bmm would read for each column kind", {
  skip_if_not_installed("bmm")
  design <- fixed_design(2, ss = 3)
  design$d2[[1]] <- NA
  expect_error(
    simulate_recovery(
      imm_model("full"), imm_pars$full,
      n_subjects = 1, n_trials = 2, trial_design = design, seed = 1
    ),
    "999"
  )
  design <- fixed_design(2, ss = 3)[c("ss", paste0("nt", 1:3))]
  design$nt2[[1]] <- NA
  expect_error(
    simulate_recovery(
      m3p_model(), c(kappa = log(8), thetat = 1, thetant = -0.5),
      n_subjects = 1, n_trials = 2, trial_design = design, seed = 1
    ),
    "as 0"
  )
})

test_that("set-size-1 weights survive a lure weight that would overflow", {
  w <- mixture3p_weights(thetat = 5, thetant = 800, set_size = c(1, 3))
  expect_equal(w$p_mem[[1]], stats::plogis(5))
  expect_identical(w$p_nt[[1]], 0)
  expect_equal(w$p_nt[[2]], 1)
})

test_that("a model naming its columns by regex is refused by name", {
  skip_if_not_installed("bmm")
  model <- bmm::mixture3p(
    resp_error = "y", nt_features = "nt", set_size = "ss", regex = TRUE
  )
  expect_error(
    simulate_recovery(
      model, c(kappa = log(8), thetat = 1, thetant = -0.5),
      n_subjects = 1, n_trials = 2,
      trial_design = fixed_design(2, ss = 4)[c("ss", paste0("nt", 1:3))]
    ),
    "regex"
  )
})
