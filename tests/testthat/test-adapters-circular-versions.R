# The circular adapters under bmm's two namings and versions (#52).
#
# CRAN bmm 1.4.1 names the location of mixture2p, mixture3p and imm `mu1`
# and gives mixture2p and mixture3p no version class. bmm's `develop`
# (1.4.2.9000, measured at c13dd2c5 on 2026-11-09) renames it `mu`, gives
# both models the versions `simple`, `slot` and `slot_averaging` with the
# class `<model>_<version>`, and adds `variable_precision = TRUE`, which
# adds `tau` and keeps the class. The adapters simulate the fixed-precision
# `simple` model (and imm's `full`, `bsc`, `abc`) only; every other
# version is refused, because simulating the simple model for it would
# score a recovery against the wrong generating process without an error.
#
# Two kinds of test below. "hand-renamed" tests take whatever model the
# installed bmm builds and rename its location to `mu` or give it a
# version class by hand, so they run on both CI legs. "real build" tests
# use the installed bmm's own objects and skip unless it names the
# location `mu` (the develop leg, or bmm 1.5.0 and later).

circular_location <- 0.5 # on the tan_half link: 2 * atan(0.5) = 0.927 rad

# the installed model with its location renamed to `mu`, wherever bmm
# keeps it: a no-op on a model that already says `mu`
mu_named <- function(model) {
  rename <- function(x) {
    if (is.null(x)) {
      return(x)
    }
    names(x)[names(x) == "mu1"] <- "mu"
    x
  }
  model$parameters <- rename(model$parameters)
  model$links <- rename(model$links)
  model$fixed_parameters <- rename(model$fixed_parameters)
  model
}

# the installed model dressed as another version: its version class
# replaced, `$version` set, as bmm's constructor sets both
as_version <- function(model, base, version) {
  cls <- class(model)
  cls <- cls[!startsWith(cls, paste0(base, "_"))]
  class(model) <- c(cls, paste0(base, "_", version))
  model$version <- version
  model
}

installed_names_mu <- function() {
  "mu" %in% names(bmm::mixture2p(resp_error = "y")$fixed_parameters)
}

skip_unless_mu_build <- function() {
  skip_if_not_installed("bmm")
  skip_if_not(
    installed_names_mu(),
    "the installed bmm names the location mu1 (bmm < 1.5.0); the \
     hand-renamed tests cover mu"
  )
}

nt_model <- function(model) {
  nt <- paste0("nt", 1:2)
  if (identical(model, "imm")) {
    return(bmm::imm(
      resp_error = "y", nt_features = nt, nt_distances = paste0("d", 1:2),
      set_size = "ss", version = "full"
    ))
  }
  bmm::mixture3p(resp_error = "y", nt_features = nt, set_size = "ss")
}

nt_design <- function(n) {
  data.frame(
    ss = rep(3, n), nt1 = rep(2, n), nt2 = rep(-2, n),
    d1 = rep(3, n), d2 = rep(3, n)
  )
}

# precise pars, almost all memory: the responses sit at the location
location_pars <- list(
  mixture2p = c(kappa = log(200), thetat = stats::qlogis(0.95)),
  mixture3p = c(kappa = log(200), thetat = 4, thetant = -4),
  imm = c(kappa = log(200), a = log(0.01), c = log(50), s = log(3)),
  sdm = c(kappa = log(50), c = log(20))
)

simulate_at_location <- function(model, name, n = 400) {
  design <- if (name %in% c("mixture3p", "imm")) {
    keep <- c("ss", "nt1", "nt2", if (name == "imm") c("d1", "d2"))
    nt_design(n)[keep]
  }
  loc <- stats::setNames(circular_location, location_of(model))
  simulate_recovery(
    model, c(location_pars[[name]], loc),
    n_subjects = 1, n_trials = n, trial_design = design, seed = 11
  )
}

expect_at_location <- function(y, label) {
  at <- 2 * atan(circular_location)
  # most responses near the location, few near 0, where an adapter that
  # dropped the location would put them
  expect_gt(mean(abs(y - at) < 0.3), 0.6, label = label)
  expect_lt(mean(abs(y) < 0.3), 0.1, label = label)
}

circular_models <- function() {
  list(
    mixture2p = bmm::mixture2p(resp_error = "y"),
    mixture3p = nt_model("mixture3p"),
    imm = nt_model("imm"),
    sdm = bmm::sdm(resp_error = "y")
  )
}

# the location -------------------------------------------------------------

test_that("the location is read off the model's own parameters, mu or mu1", {
  shape_141 <- structure(
    list(
      links = list(mu1 = "tan_half", kappa = "log", thetat = "logit"),
      fixed_parameters = list(mu1 = 0, mu2 = 0, kappa2 = -100)
    ),
    class = c("bmmodel", "circular", "mixture2p")
  )
  expect_identical(location_parameter(shape_141), "mu1")
  shape_dev <- structure(
    list(
      links = list(mu = "tan_half", kappa = "log", thetat = "logit"),
      fixed_parameters = list(mu = 0)
    ),
    class = c("bmmodel", "circular", "mixture2p", "mixture2p_simple")
  )
  expect_identical(location_parameter(shape_dev), "mu")
  # a location bmm fixes without a link is still the model's
  shape_dev$links$mu <- NULL
  expect_identical(location_parameter(shape_dev), "mu")
  # a model that names neither has no location to give: an error naming
  # it, never a guess
  shape_dev$fixed_parameters <- list(kappa2 = -100)
  expect_error(location_parameter(shape_dev), "mixture2p_simple")
})

# hand-renamed: both CI legs -----------------------------------------------

test_that("each circular adapter simulates a mu-named model (hand-renamed)", {
  skip_if_not_installed("bmm")
  models <- lapply(circular_models(), mu_named)
  for (name in names(models)) {
    sim <- simulate_at_location(models[[name]], name)
    expect_identical(sim$generator, paste0("adapter:", name), label = name)
    expect_at_location(sim$data$y, name)
  }
})

test_that("the mixture2p density reads a mu-named location (hand-renamed)", {
  skip_if_not_installed("bmm")
  model <- mu_named(bmm::mixture2p(resp_error = "y"))
  pars <- list(mu = 1, kappa = 10, thetat = 0.9)
  data <- data.frame(y = c(-1, 0, 1))
  expect_equal(
    density_for(model)(pars, data, model),
    bmm::dmixture2p(c(-1, 0, 1), mu = 1, kappa = 10, p_mem = 0.9, log = TRUE)
  )
})

# real build: the develop leg ---------------------------------------------

test_that("each circular adapter simulates bmm's own mu-named model (real build)", {
  skip_unless_mu_build()
  models <- circular_models()
  for (name in names(models)) {
    expect_identical(location_of(models[[name]]), "mu", label = name)
    sim <- simulate_at_location(models[[name]], name)
    expect_at_location(sim$data$y, name)
  }
})

# which versions an adapter takes ------------------------------------------

test_that("an adapter takes its own versions and no other", {
  model <- function(cls, version = NA_character_, vp = NULL) {
    structure(
      list(version = version, variable_precision = vp),
      class = c("bmmodel", cls)
    )
  }
  # bmm 1.4.1: mixture2p and mixture3p have one version, no version class
  # and `$version` the string "NA"; bmm's own check_data.mixture2p() reads
  # such a model as `simple`
  expect_identical(
    adapter_name(model(c("circular", "mixture2p"), "NA")), "mixture2p"
  )
  expect_identical(
    adapter_name(model(c("circular", "mixture2p"))), "mixture2p"
  )
  expect_identical(
    adapter_name(model(c("circular", "non_targets", "mixture3p"))),
    "mixture3p"
  )
  expect_identical(
    adapter_name(model(c("circular", "mixture2p", "mixture2p_simple"), "simple")),
    "mixture2p"
  )
  expect_identical(
    adapter_name(model(c("circular", "mixture3p", "mixture3p_simple"), "simple")),
    "mixture3p"
  )
  for (v in c("slot", "slot_averaging")) {
    expect_identical(
      adapter_name(model(c("circular", "mixture2p", paste0("mixture2p_", v)), v)),
      NA_character_,
      label = v
    )
    expect_identical(
      adapter_name(model(c("circular", "mixture3p", paste0("mixture3p_", v)), v)),
      NA_character_,
      label = v
    )
  }
  for (v in c("full", "bsc", "abc")) {
    expect_identical(
      adapter_name(model(c("circular", "imm", paste0("imm_", v)), v)), "imm"
    )
  }
  # a version bmm adds later gets no adapter rather than a sibling's
  expect_identical(
    adapter_name(model(c("circular", "imm", "imm_new"), "new")), NA_character_
  )
  expect_identical(
    adapter_name(model(c("circular", "sdm", "sdm_simple"), "simple")), "sdm"
  )
  # variable precision adds tau and keeps the class: no adapter draws tau
  expect_identical(
    adapter_name(model(
      c("circular", "mixture2p", "mixture2p_simple"), "simple",
      vp = TRUE
    )),
    NA_character_
  )
  expect_identical(
    adapter_name(model(c("circular", "imm", "imm_full"), "full", vp = TRUE)),
    NA_character_
  )
  expect_identical(
    adapter_name(model(
      c("circular", "mixture2p", "mixture2p_simple"), "simple",
      vp = FALSE
    )),
    "mixture2p"
  )
  # a version class without `$version` is read from the class, never as
  # the 1.4.1 single version
  expect_identical(
    adapter_name(model(c("circular", "mixture2p", "mixture2p_slot"), NULL)),
    NA_character_
  )
  expect_identical(
    adapter_name(model(c("circular", "mixture3p", "mixture3p_slot"), "NA")),
    NA_character_
  )
  # and a model of no circular class gains no class
  expect_identical(
    model_classes(model(c("ezdm", "ezdm_3par"), NULL)),
    c("bmmodel", "ezdm", "ezdm_3par")
  )
  # the ezdm rule is unchanged
  expect_identical(adapter_name(model(c("ezdm", "ezdm_3par"))), "ezdm")
})

test_that("a version the adapter cannot simulate is refused, with the way out", {
  skip_if_not_installed("bmm")
  slot <- as_version(
    mu_named(bmm::mixture2p(resp_error = "y")), "mixture2p", "slot"
  )
  expect_null(generator_for(slot))
  expect_null(density_for(slot))
  err <- expect_error(
    simulate_recovery(
      slot, c(kappa = log(10), thetat = 0),
      n_subjects = 1, n_trials = 2, seed = 1
    )
  )
  msg <- cli::ansi_strip(conditionMessage(err))
  expect_match(msg, "mixture2p")
  expect_match(msg, "\"slot\"", fixed = TRUE)
  expect_match(msg, "generator")
  expect_match(msg, "simulate_data()", fixed = TRUE)

  # the same refusal from a component
  expect_error(
    recovery_component(
      slot, c(kappa = log(10), thetat = 0),
      n_trials = 2, name = "a"
    ),
    "slot"
  )
  # and from fit_ml(), whose density would be the simple model's
  expect_error(
    fit_ml(slot, data.frame(id = 1, y = 0), method = "optim"),
    "slot"
  )

  vp <- mu_named(bmm::mixture2p(resp_error = "y"))
  vp$variable_precision <- TRUE
  err <- expect_error(
    simulate_recovery(
      vp, c(kappa = log(10), thetat = 0),
      n_subjects = 1, n_trials = 2, seed = 1
    )
  )
  expect_match(
    cli::ansi_strip(conditionMessage(err)), "variable_precision = TRUE",
    fixed = TRUE
  )
})

test_that("imm's fixed b reaches the adapter as 1 (real build)", {
  skip_unless_mu_build()
  seen <- new.env()
  recorder <- function(pars, n_trials, model, trial_design) {
    seen$pars <- pars
    data.frame(y = rep(0, n_trials))
  }
  model <- nt_model("imm")
  expect_true("b" %in% names(model$fixed_parameters))
  simulate_recovery(
    model, location_pars$imm,
    n_subjects = 1, n_trials = 2, trial_design = nt_design(2),
    generator = recorder
  )
  expect_equal(seen$pars$b, 1)
})

test_that("bmm's own slot and variable-precision models are refused (real build)", {
  skip_unless_mu_build()
  nt <- paste0("nt", 1:2)
  refused <- list(
    bmm::mixture2p(resp_error = "y", set_size = "ss", version = "slot"),
    bmm::mixture2p(
      resp_error = "y", set_size = "ss", version = "slot_averaging"
    ),
    bmm::mixture3p(
      resp_error = "y", nt_features = nt, set_size = "ss", version = "slot"
    ),
    bmm::mixture2p(resp_error = "y", variable_precision = TRUE),
    bmm::mixture3p(
      resp_error = "y", nt_features = nt, set_size = "ss",
      variable_precision = TRUE
    ),
    bmm::imm(
      resp_error = "y", nt_features = nt, nt_distances = paste0("d", 1:2),
      set_size = "ss", version = "full", variable_precision = TRUE
    )
  )
  for (model in refused) {
    expect_identical(adapter_name(model), NA_character_)
    expect_null(generator_for(model))
  }
})

test_that("a refused circular version keeps its response range", {
  skip_if_not_installed("bmm")
  slot <- as_version(bmm::mixture2p(resp_error = "y"), "mixture2p", "slot")
  expect_identical(
    response_range(slot, data.frame(y = 0)), list(floor = -pi, ceiling = pi)
  )
})

# imm's guessing weight ----------------------------------------------------

test_that("the imm adapter passes the model's b, 1 when it has none", {
  skip_if_not_installed("bmm")
  pars <- list(mu = 0, kappa = 8, a = 0.5, c = 2, s = 1.5)
  model <- mu_named(nt_model("imm"))
  # bmm 1.4.1 has no b; develop fixes it at 0 on the log link, which
  # reaches the adapter as 1
  expect_identical(imm_args(pars, model, 2, 0.5)$b, 1)
  expect_identical(imm_args(c(pars, b = 2), model, 2, 0.5)$b, 2)
})
