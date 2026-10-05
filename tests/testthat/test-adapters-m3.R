# Tests for the m3 adapter (Milestone 21.3, spec
# local/dev/spec-milestone-21-design-and-adapters.md § 10, D57, D58). Every
# expected share below is written from bmm's activation arithmetic
# (.compute_m3_probability_vector(), read on bmm 1.3.2.9000, 2026-09-30
# build): activation per category, exp() under the softmax choice rule,
# times the category's number of options, normalised. The links depend on
# the choice rule (`simple`: c, a log, b = 0.1; `softmax`: identity, b = 0).
# Nothing here compiles Stan.

mock_bmm <- function(formula, data, model) {
  bmm::bmm(
    formula, data, model,
    backend = "mock", mock_fit = 1, rename = FALSE, silent = 2
  )
}

cats3 <- c("corr", "other", "npl")

m3_ss <- function(num_options = c(1, 4, 5), choice_rule = "simple") {
  bmm::m3(
    resp_cats = cats3, num_options = num_options,
    choice_rule = choice_rule, version = "ss"
  )
}

m3_cs <- function() {
  bmm::m3(
    resp_cats = c("corr", "dc", "other", "do", "npl"),
    num_options = c(1, 1, 4, 2, 5), choice_rule = "simple", version = "cs"
  )
}

m3_custom <- function(links = list(c = "log", a = "log")) {
  bmm::m3(
    resp_cats = cats3, num_options = c(1, 4, 5), choice_rule = "simple",
    version = "custom", links = links
  )
}

# the activations listed out of resp_cats order, then the parameters
custom_formula <- function() {
  bmm::bmf(
    npl ~ b, corr ~ b + a + c, other ~ b + a,
    c ~ 1 + (1 | id), a ~ 1 + (1 | id)
  )
}

#' The category probabilities bmm's likelihood gives, written out by hand
#' for the simple choice rule: corr ~ b + a + c, other ~ b + a, npl ~ b
expected_ss <- function(c, a, b = 0.1, num_options = c(1, 4, 5)) {
  acts <- c(b + a + c, b + a, b) * num_options
  stats::setNames(acts / sum(acts), cats3)
}

#' Pooled category shares of a simulation
shares <- function(data, cats = cats3) {
  totals <- colSums(as.data.frame(data)[cats])
  totals / sum(totals)
}

#' Every share within `width` of its expected value, on the absolute scale:
#' 0.025 is about 5 SE of a share near .3 at 8000 trials, and a category
#' swapped with another is 0.2 or more away
expect_shares <- function(data, expected, width = 0.025) {
  got <- shares(data)
  expect_identical(names(got), names(expected))
  expect_lt(max(abs(got - expected)), width)
}

#' A grid fitter for a custom m3: bmm::parameters() lists no free parameter
#' for one, so the shared mock_terms() would give a fit with no terms;
#' this one carries the terms a real fit of custom_formula() has
custom_m3_fitter <- function() {
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  fitter <- function(formula, data, model, prior = NULL, ...) {
    calls$n <- calls$n + 1L
    structure(
      list(parameters = c("c", "a"), ids = levels(data$id)),
      class = "mockfit"
    )
  }
  list(fitter = fitter, calls = calls)
}

# the table ----------------------------------------------------------------

test_that("m3 has a generator adapter and no density", {
  expect_true("m3" %in% adapter_classes())
  expect_false("m3" %in% density_classes())
})

test_that("an m3 with numeric num_options reads no trial design", {
  skip_if_not_installed("bmm")
  expect_identical(trial_design_columns(m3_ss()), character())
  expect_identical(trial_design_columns(m3_custom()), character())
  model <- m3_ss(num_options = c("n_corr", "n_other", "n_npl"))
  expect_identical(
    trial_design_columns(model), c("n_corr", "n_other", "n_npl")
  )
})

# per version, against bmm's own parameter table -----------------------------

test_that("each m3 version simulates one count row per subject", {
  skip_if_not_installed("bmm")
  cases <- list(
    list(model = m3_ss(), pars = c(c = log(3), a = log(0.5))),
    list(
      model = m3_ss(choice_rule = "softmax"), pars = c(c = 2, a = 0.5)
    ),
    list(model = m3_cs(), pars = c(c = log(3), a = log(0.5), f = 0)),
    list(
      model = m3_custom(), pars = c(c = log(3), a = log(0.5)),
      formula = custom_formula()
    )
  )
  for (case in cases) {
    cats <- case$model$resp_vars$resp_cats
    sim <- simulate_recovery(
      case$model, case$pars,
      n_subjects = 4, n_trials = 50, sds = c(c = 0.2),
      formula = case$formula, seed = 1
    )
    expect_identical(sim$generator, "adapter:m3")
    expect_named(sim$data, c("id", cats))
    expect_identical(nrow(sim$data), 4L)
    expect_true(all(rowSums(as.data.frame(sim$data)[cats]) == 50))
    # the truth names exactly the free parameters: bmm::parameters() for
    # ss and cs; for custom, bmm learns them from the formula at fit time,
    # so its table lists none and the links are what names them
    p <- bmm::parameters(case$model)
    free <- if (inherits(case$model, "m3_custom")) {
      names(case$model$links)
    } else {
      p$parameter[!p$fixed]
    }
    expect_setequal(sim$truth$population$term, free)
  }
})

test_that("tasks give one count row per subject and task", {
  skip_if_not_installed("bmm")
  sim <- simulate_recovery(
    m3_ss(), c(c = log(3), a = log(0.5)),
    n_subjects = 3, n_trials = 20, tasks = c("1", "2"), seed = 2
  )
  expect_identical(nrow(sim$data), 6L)
  expect_true(all(rowSums(as.data.frame(sim$data)[cats3]) == 20))
})

# data bmm accepts ------------------------------------------------------------

test_that("adapter data pass bmm's checks for every version (mock backend)", {
  skip_if_not_installed("bmm")
  skip_if_not_installed("brms")
  ss <- simulate_recovery(
    m3_ss(), c(c = log(3), a = log(0.5)),
    n_subjects = 3, n_trials = 30, seed = 3
  )
  expect_s3_class(
    mock_bmm(recovery_formula(ss$model), ss$data, ss$model), "bmmfit"
  )
  cs <- simulate_recovery(
    m3_cs(), c(c = log(3), a = log(0.5), f = 0),
    n_subjects = 3, n_trials = 30, seed = 3
  )
  expect_s3_class(
    mock_bmm(recovery_formula(cs$model), cs$data, cs$model), "bmmfit"
  )
  custom <- simulate_recovery(
    m3_custom(), c(c = log(3), a = log(0.5)),
    n_subjects = 3, n_trials = 30, formula = custom_formula(), seed = 3
  )
  # bmm warns that it chose priors for the custom model's parameters;
  # that is bmm's notice about fitting, not about the data. Silencing it
  # with options(bmm.default_priors = FALSE) instead trips CRAN bmm
  # 1.3.2's combine_prior() for a custom m3 (measured 2026-10-05)
  expect_warning(
    fit <- mock_bmm(custom_formula(), custom$data, custom$model),
    "Default priors"
  )
  expect_s3_class(fit, "bmmfit")
})

# category placement ----------------------------------------------------------

test_that("each category lands in the column resp_cats names", {
  skip_if_not_installed("bmm")
  # c = 5, a = 0.5 on the natural scale: shares .659 / .282 / .059. A
  # category swapped with another, or the option counts applied in the
  # formula's order (bmm's bug before the fix), is 0.2 or more away
  expected <- expected_ss(c = 5, a = 0.5)
  pars <- c(c = log(5), a = log(0.5))
  custom <- simulate_recovery(
    m3_custom(), pars,
    n_subjects = 4, n_trials = 2000, formula = custom_formula(), seed = 4
  )
  expect_shares(custom$data, expected)
  ss <- simulate_recovery(
    m3_ss(), pars,
    n_subjects = 4, n_trials = 2000, seed = 4
  )
  expect_shares(ss$data, expected)
})

test_that("num_options named after the categories in another order", {
  skip_if_not_installed("bmm")
  model <- m3_ss(num_options = c(npl = 5, corr = 1, other = 4))
  sim <- simulate_recovery(
    model, c(c = log(5), a = log(0.5)),
    n_subjects = 4, n_trials = 2000, seed = 5
  )
  expect_shares(sim$data, expected_ss(c = 5, a = 0.5))
})

test_that("the softmax choice rule exponentiates the activations", {
  skip_if_not_installed("bmm")
  model <- m3_ss(choice_rule = "softmax")
  sim <- simulate_recovery(
    model, c(c = 1, a = 0.5),
    n_subjects = 4, n_trials = 2000, seed = 6
  )
  # b = 0 under softmax
  acts <- exp(c(0 + 0.5 + 1, 0 + 0.5, 0)) * c(1, 4, 5)
  expect_shares(sim$data, stats::setNames(acts / sum(acts), cats3))
})

# the custom version and its formula (D57, D58) -------------------------------

test_that("a custom m3 without links is refused with a hint", {
  skip_if_not_installed("bmm")
  expect_error(
    simulate_recovery(
      m3_custom(links = NULL), c(c = log(3), a = log(0.5)),
      n_subjects = 2, n_trials = 10, formula = custom_formula()
    ),
    "links"
  )
})

test_that("a custom m3 needs the formula with its activations", {
  skip_if_not_installed("bmm")
  pars <- c(c = log(3), a = log(0.5))
  expect_error(
    simulate_recovery(m3_custom(), pars, n_subjects = 2, n_trials = 10),
    "needs `formula`"
  )
  # a formula without the activation of one category
  expect_error(
    simulate_recovery(
      m3_custom(), pars,
      n_subjects = 2, n_trials = 10,
      formula = bmm::bmf(corr ~ b + a + c, other ~ b + a, c ~ 1, a ~ 1)
    ),
    "\"npl\""
  )
  # an activation using a name that is no parameter of the model
  expect_error(
    simulate_recovery(
      m3_custom(), pars,
      n_subjects = 2, n_trials = 10,
      formula = bmm::bmf(npl ~ b, corr ~ b + a + c + d, other ~ b + a)
    ),
    "\"d\""
  )
  # the formula is checked before anything is drawn
  withr::local_seed(1)
  before <- .Random.seed
  expect_error(
    simulate_recovery(m3_custom(), pars, n_subjects = 2, n_trials = 10)
  )
  expect_identical(.Random.seed, before)
})

test_that("a formula reaches only an adapter that reads one", {
  skip_if_not_installed("bmm")
  expect_error(
    simulate_recovery(
      bmm::mixture2p(resp_error = "y"),
      c(kappa = log(8), thetat = 1),
      n_subjects = 2, n_trials = 10,
      formula = bmm::bmf(kappa ~ 1, thetat ~ 1)
    ),
    "does not read a formula"
  )
  generator <- function(pars, n_trials, model) {
    data.frame(y = stats::runif(n_trials, -1, 1))
  }
  expect_error(
    simulate_recovery(
      bmm::mixture2p(resp_error = "y"),
      c(kappa = log(8), thetat = 1),
      n_subjects = 2, n_trials = 10,
      formula = bmm::bmf(kappa ~ 1, thetat ~ 1), generator = generator
    ),
    "closes over"
  )
  expect_error(
    simulate_recovery(
      m3_custom(), c(c = 1, a = 0),
      n_subjects = 2, n_trials = 10, formula = "corr ~ b"
    ),
    "bmmformula"
  )
})

test_that("the simulation records the formula that reached the adapter", {
  skip_if_not_installed("bmm")
  sim <- simulate_recovery(
    m3_custom(), c(c = log(3), a = log(0.5)),
    n_subjects = 2, n_trials = 10, formula = custom_formula(), seed = 1
  )
  expect_identical(sim$formula, formula_key(custom_formula()))
  expect_output(print(sim), "activation formulas: npl, corr, other")
  plain <- simulate_recovery(
    m3_ss(), c(c = log(3), a = log(0.5)),
    n_subjects = 2, n_trials = 10, seed = 1
  )
  expect_null(plain$formula)
  expect_no_match(capture.output(print(plain)), "activation")
})

test_that("a simulation stored before the formula field is upgraded", {
  sim <- list(sds = c(kappa = 0), truth = list(sd = tibble::tibble()))
  expect_true("formula" %in% names(upgrade_simulation(sim)))
  expect_null(upgrade_simulation(sim)$formula)
})

# num_options as columns: a trial design ------------------------------------

test_that("num_options columns make one trial per design row", {
  skip_if_not_installed("bmm")
  skip_if_not_installed("brms")
  model <- m3_ss(num_options = c("n_corr", "n_other", "n_npl"))
  design <- data.frame(
    n_corr = 1, n_other = rep(c(0, 4), 10), n_npl = 5
  )
  sim <- simulate_recovery(
    model, c(c = log(3), a = log(0.5)),
    n_subjects = 3, n_trials = 20, trial_design = design, seed = 7
  )
  expect_named(sim$data, c("id", "n_corr", "n_other", "n_npl", cats3))
  expect_identical(nrow(sim$data), 60L)
  # one trial per row: one count of 1
  expect_true(all(rowSums(as.data.frame(sim$data)[cats3]) == 1))
  # no option, no response in that category
  expect_true(all(sim$data$other[sim$data$n_other == 0] == 0))
  expect_s3_class(
    mock_bmm(recovery_formula(model), sim$data, model), "bmmfit"
  )
})

test_that("num_options columns follow each row's numbers", {
  skip_if_not_installed("bmm")
  model <- m3_ss(num_options = c("n_corr", "n_other", "n_npl"))
  design <- data.frame(n_corr = 1, n_other = 4, n_npl = 5)
  design <- design[rep(1, 1000), ]
  sim <- simulate_recovery(
    model, c(c = log(5), a = log(0.5)),
    n_subjects = 2, n_trials = 1000, trial_design = design, seed = 8
  )
  expect_shares(sim$data, expected_ss(c = 5, a = 0.5), width = 0.035)
})

test_that("num_options columns need a design without NA or negatives", {
  skip_if_not_installed("bmm")
  model <- m3_ss(num_options = c("n_corr", "n_other", "n_npl"))
  pars <- c(c = log(3), a = log(0.5))
  expect_error(
    simulate_recovery(model, pars, n_subjects = 1, n_trials = 2),
    "needs `trial_design`"
  )
  design <- data.frame(n_corr = 1, n_other = c(4, NA), n_npl = 5)
  expect_error(
    simulate_recovery(
      model, pars,
      n_subjects = 1, n_trials = 2, trial_design = design, seed = 1
    ),
    "n_other"
  )
  design$n_other <- c(4, -1)
  expect_error(
    simulate_recovery(
      model, pars,
      n_subjects = 1, n_trials = 2, trial_design = design, seed = 1
    ),
    "n_other"
  )
})

# forwarding: grid, component, sbc --------------------------------------------

test_that("a grid hands its formula to the custom m3 adapter and records it", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- custom_m3_fitter()
  run <- function(formula) {
    suppressMessages(recovery_grid(
      m3_custom(),
      grid = small_grid(), reps = 1L,
      pars = c(c = log(3), a = log(0.5)),
      dir = dir, seed = 1, formula = formula, preflight = FALSE,
      .fitter = mock$fitter
    ))
  }
  run(custom_formula())
  for (row in 1:2) {
    sim <- readRDS(file.path(dir, sprintf("cell-%d-rep-1-sim.rds", row)))
    expect_named(sim$data, c("id", cats3))
    expect_identical(sim$formula, formula_key(custom_formula()))
  }
  record <- readRDS(file.path(dir, "grid.rds"))
  expect_identical(record$formula, formula_key(custom_formula()[cats3]))
})

test_that("a grid of a custom m3 without activations stops before writing", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  expect_error(
    recovery_grid(
      m3_custom(),
      grid = small_grid(), reps = 1L,
      pars = c(c = log(3), a = log(0.5)),
      dir = dir, seed = 1, .fitter = mock$fitter
    ),
    "activation"
  )
  expect_false(file.exists(file.path(dir, "grid.rds")))
  expect_identical(mock$calls$n, 0L)
})

test_that("a grid of any other model records no formula", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- grid_mock_fitter()
  suppressMessages(recovery_grid(
    bmm::mixture2p(resp_error = "y"),
    grid = small_grid(), reps = 1L,
    pars = c(kappa = log(8), thetat = 1),
    dir = dir, seed = 1, .fitter = mock$fitter
  ))
  expect_false("formula" %in% names(readRDS(file.path(dir, "grid.rds"))))
  sim <- readRDS(file.path(dir, "cell-1-rep-1-sim.rds"))
  expect_null(sim$formula)
})

test_that("a component hands its formula to the custom m3 adapter", {
  skip_if_not_installed("bmm")
  expect_error(
    recovery_component(
      m3_custom(),
      pars = c(c = log(3), a = log(0.5)), n_trials = 20, name = "mmm"
    ),
    "activation"
  )
  comp <- recovery_component(
    m3_custom(),
    pars = c(c = log(3), a = log(0.5)), n_trials = 20,
    formula = custom_formula(), name = "mmm"
  )
  set <- simulate_components(list(comp), n_subjects = 2, seed = 3)
  data <- set$components$mmm$data
  expect_named(data, c("id", cats3))
  expect_true(all(rowSums(as.data.frame(data)[cats3]) == 20))
})

test_that("sbc()'s default generator simulates a custom m3 from its formula", {
  skip_if_not_installed("SBC")
  skip_if_not_installed("bmm")
  skip_on_cran() # a whole sbc() run, as in every sbc() test
  model <- m3_custom()
  # bmm's notice that it chose priors for the custom model, as above
  withr::local_options(bmm.default_priors = FALSE)
  data <- data.frame(id = c("b", "a", "c"), corr = 20, other = 5, npl = 5)
  k <- seq_len(60)
  draws <- posterior::as_draws_matrix(cbind(
    b_c_Intercept = log(3) + 0.01 * k,
    b_a_Intercept = log(0.5) + 0.01 * k
  ))
  mock <- sbc_mock_fitter(draws = draws)
  formula <- bmm::bmf(npl ~ b, corr ~ b + a + c, other ~ b + a, c ~ 1, a ~ 1)
  suppressMessages(withCallingHandlers(
    sbc(
      model, formula, data,
      n_sims = 2L, thin_ranks = 1, cores_per_fit = 1, level = "population",
      .fitter = mock$fitter
    ),
    bmmtools_sbc_diagnostics = function(w) invokeRestart("muffleWarning")
  ))
  calls <- sbc_dataset_calls(mock$calls)
  expect_length(calls, 2L)
  for (call in calls) {
    # the activations reached the adapter: one count row per subject. The
    # trials per row are sbc_layout()'s rows per subject, not data's 30,
    # as for every count model (STATE-milestone-21, 21.3, open question)
    expect_named(call$data, c("id", cats3))
    totals <- rowSums(as.data.frame(call$data)[cats3])
    expect_length(unique(totals), 1L)
  }
})

# prior_check -----------------------------------------------------------------

test_that("prior_check's range for m3 is 0 to each row's trial count", {
  skip_if_not_installed("bmm")
  data <- data.frame(corr = c(10, 3), other = c(5, NA), npl = c(5, 1))
  expect_identical(
    response_range(m3_ss(), data),
    list(floor = 0, ceiling = c(20, 4))
  )
})

test_that("the default summary gives one block per m3 category", {
  skip_if_not_installed("bmm")
  model <- m3_ss()
  data <- data.frame(corr = c(10, 2), other = c(0, 1), npl = c(0, 1))
  # draws x rows x categories, as brms documents a multinomial prediction;
  # draw 1 puts everything in corr, draw 2 nothing
  yrep <- array(0, dim = c(2, 2, 3))
  yrep[1, , 1] <- c(10, 4)
  yrep[2, , 2] <- c(10, 4)
  out <- default_prior_summary(response_range(model, data), model)(yrep, data)
  expect_setequal(unique(out$response), cats3)
  ceiling_corr <- out$value[
    out$response == "corr" & out$statistic == "ceiling_rate"
  ]
  floor_corr <- out$value[
    out$response == "corr" & out$statistic == "floor_rate"
  ]
  expect_equal(ceiling_corr, 0.5)
  expect_equal(floor_corr, 0.5)
})

# branches added after the gate's coverage run (never RED) --------------------

test_that("counts without column names keep resp_cats order", {
  skip_if_not_installed("bmm")
  # a bmm whose rm3() returned unnamed columns: the activations are passed
  # in resp_cats order, so position is category
  local_mocked_bindings(
    rm3 = function(n, size, pars, m3_model, act_funs = NULL, ...) {
      matrix(c(3, 2, 1), nrow = 1)
    },
    .package = "bmm"
  )
  out <- generate_m3(
    list(c = 1, a = 1, b = 0.1), 6, m3_ss()
  )
  expect_identical(out, data.frame(corr = 3, other = 2, npl = 1))
})

test_that("prior_check's m3 range needs every category column", {
  skip_if_not_installed("bmm")
  expect_message(
    out <- response_range(m3_ss(), data.frame(corr = 1, other = 2)),
    "\"npl\""
  )
  expect_identical(out, list(floor = NA_real_, ceiling = NA_real_))
})

test_that("a grid records a formula function by its activations", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- custom_m3_fitter()
  formula <- function(row) custom_formula()
  suppressMessages(recovery_grid(
    m3_custom(),
    grid = small_grid(), reps = 1L,
    pars = c(c = log(3), a = log(0.5)),
    dir = dir, seed = 1, formula = formula, preflight = FALSE,
    .fitter = mock$fitter
  ))
  record <- readRDS(file.path(dir, "grid.rds"))
  expect_identical(record$formula, formula_key(custom_formula()[cats3]))
  sim <- readRDS(file.path(dir, "cell-2-rep-1-sim.rds"))
  expect_identical(sim$formula, formula_key(custom_formula()))
})

# review of 21.3 --------------------------------------------------------------

test_that("a free parameter no activation uses is refused", {
  skip_if_not_installed("bmm")
  # `a` has a link but no activation reads it: its truth would be a
  # parameter the data never depended on
  expect_error(
    simulate_recovery(
      m3_custom(), c(c = log(3), a = 0),
      n_subjects = 2, n_trials = 10,
      formula = bmm::bmf(corr ~ b + c, other ~ b, npl ~ b, c ~ 1)
    ),
    "\"a\""
  )
})

test_that("the activations reach rm3() in resp_cats order", {
  skip_if_not_installed("bmm")
  # released bmm 1.3.2 matched options to activations by position, so
  # the order handed over is what keeps categories in place there
  seen <- NULL
  local_mocked_bindings(
    rm3 = function(n, size, pars, m3_model, act_funs = NULL, ...) {
      seen <<- names(act_funs)
      matrix(c(1, 1, 1), nrow = 1, dimnames = list(NULL, cats3))
    },
    .package = "bmm"
  )
  generate_m3(
    list(c = 1, a = 1, b = 0.1), 3, m3_custom(),
    formula = custom_formula()
  )
  expect_identical(seen, cats3)
})

test_that("option counts must be finite, numeric and not all zero", {
  skip_if_not_installed("bmm")
  model <- m3_ss(num_options = c("n_corr", "n_other", "n_npl"))
  pars <- c(c = log(3), a = log(0.5))
  run <- function(design) {
    simulate_recovery(
      model, pars,
      n_subjects = 1, n_trials = 2, trial_design = design, seed = 1
    )
  }
  expect_error(
    run(data.frame(n_corr = 1, n_other = c(4, Inf), n_npl = 5)), "n_other"
  )
  expect_error(
    run(data.frame(n_corr = c(1, 0), n_other = c(4, 0), n_npl = c(5, 0))),
    "Trial 2"
  )
  # a factor column is read as the numbers it shows, not its codes
  sim <- run(data.frame(
    n_corr = 1, n_other = factor(c("0", "0")), n_npl = 5
  ))
  expect_true(all(sim$data$other == 0))
})

test_that("a grid of a custom m3 runs to the end and scores c and a", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- custom_m3_fitter()
  out <- suppressMessages(recovery_grid(
    m3_custom(),
    grid = small_grid(), reps = 1L,
    pars = c(c = log(3), a = log(0.5)),
    dir = dir, seed = 1, formula = custom_formula(), preflight = FALSE,
    .fitter = mock$fitter
  ))
  expect_setequal(out$term[out$level == "population"], c("c", "a"))
  expect_identical(attr(out, "cells")$status, c("ok", "ok"))
})

test_that("the grid records the activations, not the parameter formulas", {
  skip_if_not_installed("bmm")
  dir <- withr::local_tempdir()
  mock <- custom_m3_fitter()
  # the messages a run gives, collected rather than suppressed
  run <- function(formula) {
    said <- character()
    withCallingHandlers(
      recovery_grid(
        m3_custom(),
        grid = small_grid(), reps = 1L,
        pars = c(c = log(3), a = log(0.5)),
        dir = dir, seed = 1, formula = formula, preflight = FALSE,
        .fitter = mock$fitter
      ),
      message = function(m) {
        said <<- c(said, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
    said
  }
  run(custom_formula())
  record <- readRDS(file.path(dir, "grid.rds"))
  expect_identical(record$formula, formula_key(custom_formula()[cats3]))
  # a new random-effects structure leaves the data's activations alone
  wider <- bmm::bmf(
    npl ~ b, corr ~ b + a + c, other ~ b + a,
    c ~ 1 + (1 | id), a ~ 1
  )
  expect_false(any(grepl("Rewriting", run(wider))))
  # new activations are named
  other <- bmm::bmf(
    npl ~ b, corr ~ b + a + c, other ~ b + 2 * a,
    c ~ 1 + (1 | id), a ~ 1 + (1 | id)
  )
  expect_true(any(grepl("Rewriting.*formula", run(other))))
})
