# Tests for the four-parameter ezdm adapter (Milestone 13.1, #8, D100).
#
# Measured 2026-10-25 on CRAN bmm 1.4.1: `ezdm(version = "4par")` has the
# class vector bmmodel/ezdm/ezdm_4par and estimates `zr` (logit link) next
# to the three-parameter version's `drift`, `bound` and `ndt`. Its data
# hold the RT summaries per boundary, named `mean_rt = c(upper, lower)`
# and `var_rt = c(upper, lower)`, which `rezdm(version = "4par")` returns
# as `mean_rt_upper`, `mean_rt_lower`, `var_rt_upper`, `var_rt_lower`, with
# `NA` at a boundary reached fewer than twice. Before this adapter a 4par
# model got none (the cap in generator_for()), and before that cap it
# would have got three-parameter data at zr = 0.5 without complaint.
# Nothing here compiles Stan.

mock_bmm <- function(formula, data, model) {
  bmm::bmm(
    formula, data, model,
    backend = "mock", mock_fit = 1, rename = FALSE, silent = 2
  )
}

ez3_model <- function() {
  bmm::ezdm(mean_rt = "mrt", var_rt = "vrt", n_upper = "nu", n_trials = "n")
}

ez4_model <- function() {
  bmm::ezdm(
    mean_rt = c("mrt_u", "mrt_l"), var_rt = c("vrt_u", "vrt_l"),
    n_upper = "nu", n_trials = "n", version = "4par"
  )
}

ez4_pars <- c(
  drift = 1, bound = log(1.5), ndt = log(0.3), zr = stats::qlogis(0.6)
)

ez4_columns <- c("mrt_u", "mrt_l", "vrt_u", "vrt_l", "nu")

# the table ------------------------------------------------------------------

test_that("each ezdm version has the class vector the dispatch reads", {
  skip_if_not_installed("bmm")
  expect_identical(class(ez3_model()), c("bmmodel", "ezdm", "ezdm_3par"))
  expect_identical(class(ez4_model()), c("bmmodel", "ezdm", "ezdm_4par"))
})

test_that("each ezdm version resolves to its own adapter", {
  skip_if_not_installed("bmm")
  expect_identical(adapter_name(ez3_model()), "ezdm")
  expect_identical(adapter_name(ez4_model()), "ezdm_4par")
  expect_identical(generator_for(ez3_model()), generate_ezdm)
  expect_identical(generator_for(ez4_model()), generate_ezdm_4par)
  expect_identical(density_for(ez3_model()), density_ezdm)
  expect_identical(density_for(ez4_model()), density_ezdm_4par)
})

test_that("an ezdm version bmmtools does not know gets no adapter", {
  # a later bmm version must not inherit a sibling's generator (D100)
  later <- structure(list(), class = c("bmmodel", "ezdm", "ezdm_5par"))
  expect_true(is.na(adapter_name(later)))
  expect_null(generator_for(later))
  expect_null(density_for(later))
})

# simulate_recovery ----------------------------------------------------------

test_that("ezdm 4par simulates one row of per-boundary summaries", {
  skip_if_not_installed("bmm")
  sim <- simulate_recovery(
    ez4_model(), ez4_pars,
    n_subjects = 3, n_trials = 200, seed = 1
  )
  expect_identical(sim$generator, "adapter:ezdm_4par")
  expect_named(sim$data, c("id", ez4_columns, "n"))
  expect_equal(sim$data$n, rep(200, 3))
  expect_true(all(sim$data$nu <= sim$data$n))
  # the truth names exactly bmm's free parameters, zr among them
  p <- bmm_parameter_info(ez4_model())
  expect_setequal(sim$truth$population$term, p$parameter[!p$fixed])
  expect_true("zr" %in% sim$truth$population$term)
})

test_that("every parameter reaches rezdm() on its natural scale", {
  skip_if_not_installed("bmm")
  seen <- NULL
  local_mocked_bindings(
    rezdm = function(n, n_trials, drift, bound, ndt, zr = 0.5, s = 1,
                     version = c("3par", "4par")) {
      seen <<- list(
        n_trials = n_trials, drift = drift, bound = bound, ndt = ndt,
        zr = zr, s = s, version = version
      )
      data.frame(
        mean_rt_upper = 0.5, mean_rt_lower = 0.6, var_rt_upper = 0.01,
        var_rt_lower = 0.02, n_upper = 10, n_trials = n_trials
      )
    },
    .package = "bmm"
  )
  simulate_recovery(ez4_model(), ez4_pars, n_subjects = 1, n_trials = 30)
  expect_equal(seen, list(
    n_trials = 30, drift = 1, bound = 1.5, ndt = 0.3, zr = 0.6, s = 1,
    version = "4par"
  ))
})

test_that("the upper and lower summaries land in the columns the model names", {
  skip_if_not_installed("bmm")
  model <- ez4_model()
  pars <- natural_pars(ez4_pars, model)
  out <- withr::with_seed(3, generate_ezdm_4par(pars, 50L, model))
  raw <- withr::with_seed(3, bmm::rezdm(
    1L, 50L,
    drift = 1, bound = 1.5, ndt = 0.3, zr = 0.6, s = 1, version = "4par"
  ))
  expect_equal(
    unlist(out[c(ez4_columns, "n")], use.names = FALSE),
    unlist(raw[c(
      "mean_rt_upper", "mean_rt_lower", "var_rt_upper", "var_rt_lower",
      "n_upper", "n_trials"
    )], use.names = FALSE)
  )
})

test_that("a 4par model naming one column per RT statistic is refused", {
  skip_if_not_installed("bmm")
  # bmm's constructor accepts it; only its data check would object
  one <- bmm::ezdm(
    mean_rt = "mrt", var_rt = "vrt", n_upper = "nu", n_trials = "n",
    version = "4par"
  )
  expect_error(
    simulate_recovery(one, ez4_pars, n_subjects = 1, n_trials = 30),
    "two column names"
  )
})

test_that("ezdm 4par data pass bmm's checks (mock backend)", {
  skip_if_not_installed("bmm")
  skip_if_not_installed("brms")
  sim <- simulate_recovery(
    ez4_model(), ez4_pars,
    n_subjects = 3, n_trials = 200, seed = 5
  )
  # at this design point both boundaries are reached often in every cell
  expect_false(anyNA(sim$data))
  fit <- mock_bmm(recovery_formula(ez4_model()), sim$data, ez4_model())
  expect_s3_class(fit, "bmmfit")
})

test_that("a boundary reached fewer than twice is NA, and bmm keeps it", {
  skip_if_not_installed("bmm")
  skip_if_not_installed("brms")
  sim <- simulate_recovery(
    ez4_model(),
    c(drift = 2, bound = log(2), ndt = log(0.3), zr = stats::qlogis(0.8)),
    n_subjects = 4, n_trials = 10, seed = 1
  )
  expect_true(anyNA(sim$data$mrt_l))
  expect_false(anyNA(sim$data[c("nu", "n")]))
  # bmm >= 1.4.1 (#430) keeps such cells and says how many there are
  expect_warning(
    fit <- mock_bmm(recovery_formula(ez4_model()), sim$data, ez4_model()),
    "fewer than two"
  )
  expect_s3_class(fit, "bmmfit")
})

# trial design ---------------------------------------------------------------

test_that("ezdm 4par simulates each row from its own number of trials", {
  skip_if_not_installed("bmm")
  expect_identical(count_design_columns(ez4_model()), "n")
  expect_identical(trial_design_columns(ez4_model()), character())
  sim <- simulate_recovery(
    ez4_model(), ez4_pars,
    n_subjects = 2, n_trials = 2,
    trial_design = data.frame(n = c(100, 3)), seed = 1
  )
  expect_named(sim$data, c("id", "n", ez4_columns))
  expect_identical(sim$data$n, c(100, 3, 100, 3))
  expect_true(all(sim$data$nu <= sim$data$n))
})

test_that("ezdm 4par refuses fewer than 3 trials in a row, before any draw", {
  skip_if_not_installed("bmm")
  expect_error(
    simulate_recovery(
      ez4_model(), ez4_pars,
      n_subjects = 1, n_trials = 2,
      trial_design = data.frame(n = c(100, 2)), seed = 1
    ),
    "at least 3"
  )
})

# density, range, validation -------------------------------------------------

test_that("the 4par density reads both boundaries and zr", {
  skip_if_not_installed("bmm")
  model <- ez4_model()
  pars <- natural_pars(ez4_pars, model)
  dat <- simulate_recovery(
    model, ez4_pars,
    n_subjects = 4, n_trials = 40, seed = 2
  )$data
  ll <- density_for(model)(pars, dat, model)
  expect_equal(ll, bmm::dezdm(
    as.matrix(dat[c("mrt_u", "mrt_l")]), as.matrix(dat[c("vrt_u", "vrt_l")]),
    dat$nu, dat$n,
    drift = 1, bound = 1.5, ndt = 0.3, zr = 0.6, s = 1,
    version = "4par", log = TRUE
  ))
  # zr is read: the same data score differently at another starting point
  other <- pars
  other$zr <- 0.4
  expect_false(isTRUE(all.equal(ll, density_for(model)(other, dat, model))))
})

test_that("prior_check knows no response range for ezdm 4par, as for 3par", {
  skip_if_not_installed("bmm")
  data <- data.frame(n = 10, nu = 5)
  expect_message(out <- response_range(ez4_model(), data), "No response range")
  expect_true(is.na(out$floor))
  expect_true(is.na(out$ceiling))
})

test_that("a 4par model is not planned under the 3par ezdm defaults", {
  skip_if_not_installed("bmm")
  expect_identical(validation_model_name(ez3_model()), "ezdm")
  expect_identical(validation_model_name(ez4_model()), "ezdm_4par")
  expect_null(validation_defaults(ez4_model()))
})
