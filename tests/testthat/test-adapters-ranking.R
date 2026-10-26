# Tests for the sdt_ranking adapter (Milestone 13.2, #8, D101).
#
# Measured 2026-10-25 on CRAN bmm 1.4.1: `sdt_ranking(response, m, dist)`
# names one rank-count column per rank and no trial-count column (bmm takes
# each row's total). `m`, the set size, is a number or a column name.
# `rsdt_ranking(n, n_trials, m, d, sdratio, dist)` takes vectors of
# `n_trials` and `m` and returns `rank1..rank<max m>`, zero beyond a row's
# `m`. The model's `sdratio` (dist = "normal") is the LOG of the SD ratio on
# an identity link, fixed at 0; `rsdt_ranking()` and `dsdt_ranking()` take
# the ratio itself. `dsdt_ranking()` refuses NA counts, while bmm's data
# check reads NA beyond a row's `m` as 0. Nothing here compiles Stan.

mock_bmm <- function(formula, data, model) {
  bmm::bmm(
    formula, data, model,
    backend = "mock", mock_fit = 1, rename = FALSE, silent = 2
  )
}

ranks <- paste0("k", 1:4)

# a number for m needs exactly m rank columns (bmm's constructor checks it)
rank_model <- function(m = 4, dist = "gumbel_min") {
  response <- if (is.numeric(m)) paste0("k", seq_len(m)) else ranks
  bmm::sdt_ranking(response = response, m = m, dist = dist)
}

# a design for a model whose m names the column `mm`
rank_design <- function(n = c(60, 40, 50), mm = c(2, 3, 4)) {
  data.frame(mm = mm, n_trials = n)
}

# the table ------------------------------------------------------------------

test_that("sdt_ranking has the class vector the dispatch reads", {
  skip_if_no_bmm_model("sdt_ranking", density = TRUE)
  expect_identical(class(rank_model()), c("bmmodel", "sdt", "sdt_ranking"))
  expect_identical(adapter_name(rank_model()), "sdt_ranking")
  expect_identical(generator_for(rank_model()), generate_sdt_ranking)
  expect_identical(density_for(rank_model()), density_sdt_ranking)
})

test_that("sdt_ranking names its design columns after m", {
  skip_if_no_bmm_model("sdt_ranking")
  expect_identical(trial_design_columns(rank_model()), character())
  expect_identical(count_design_columns(rank_model()), "n_trials")
  expect_identical(trial_design_columns(rank_model("mm")), "mm")
  expect_identical(count_design_columns(rank_model("mm")), "n_trials")
})

# a number for m -------------------------------------------------------------

test_that("sdt_ranking simulates one row of rank counts per subject", {
  skip_if_no_bmm_model("sdt_ranking")
  sim <- simulate_recovery(
    rank_model(), c(d = 1),
    n_subjects = 3, n_trials = 80, seed = 1
  )
  expect_identical(sim$generator, "adapter:sdt_ranking")
  expect_named(sim$data, c("id", ranks))
  expect_equal(unname(rowSums(as.data.frame(sim$data)[ranks])), rep(80, 3))
  p <- bmm_parameter_info(rank_model())
  expect_setequal(sim$truth$population$term, p$parameter[!p$fixed])
})

test_that("every parameter reaches rsdt_ranking(), sdratio as the ratio", {
  skip_if_no_bmm_model("sdt_ranking")
  seen <- list()
  local_mocked_bindings(
    rsdt_ranking = function(n, n_trials, m, d, sdratio = 1,
                            dist = c("gumbel_min", "normal")) {
      seen[[length(seen) + 1L]] <<- list(
        n = n, n_trials = n_trials, m = m, d = d, sdratio = sdratio,
        dist = dist
      )
      matrix(
        1, n, max(m),
        dimnames = list(NULL, paste0("rank", seq_len(max(m))))
      )
    },
    .package = "bmm"
  )
  simulate_recovery(rank_model(), c(d = 1.2), n_subjects = 1, n_trials = 30)
  # the model's sdratio is the log ratio: the free value log(1.5) is a
  # ratio of 1.5, and the fixed 0 a ratio of 1
  simulate_recovery(
    rank_model(dist = "normal"), c(d = 0.8, sdratio = log(1.5)),
    n_subjects = 1, n_trials = 30
  )
  simulate_recovery(
    rank_model(3, dist = "normal"), c(d = 0.8),
    n_subjects = 1, n_trials = 30
  )
  expect_equal(seen[[1]], list(
    n = 1L, n_trials = 30, m = 4L, d = 1.2, sdratio = 1, dist = "gumbel_min"
  ))
  expect_equal(seen[[2]]$sdratio, 1.5)
  expect_equal(seen[[2]]$dist, "normal")
  expect_equal(seen[[3]]$sdratio, 1)
  expect_equal(seen[[3]]$m, 3L)
})

test_that("the rank counts land in the columns the model names, in order", {
  skip_if_no_bmm_model("sdt_ranking")
  model <- rank_model(dist = "normal")
  pars <- natural_pars(c(d = 1, sdratio = log(2)), model)
  out <- withr::with_seed(7, generate_sdt_ranking(pars, 200L, model))
  raw <- withr::with_seed(7, bmm::rsdt_ranking(
    1L, 200L,
    m = 4, d = 1, sdratio = 2, dist = "normal"
  ))
  expect_equal(unlist(out[ranks], use.names = FALSE), as.vector(raw))
})

test_that("sdt_ranking data pass bmm's checks (mock backend)", {
  skip_if_no_bmm_model("sdt_ranking")
  skip_if_not_installed("brms")
  cases <- list(
    list(model = rank_model(), pars = c(d = 1)),
    list(
      model = rank_model(dist = "normal"),
      pars = c(d = 1, sdratio = log(1.5))
    ),
    list(model = rank_model(3), pars = c(d = 1))
  )
  for (case in cases) {
    sim <- simulate_recovery(
      case$model, case$pars,
      n_subjects = 3, n_trials = 40, seed = 5
    )
    fit <- mock_bmm(recovery_formula(case$model), sim$data, case$model)
    expect_s3_class(fit, "bmmfit")
  }
})

test_that("a count design gives each row its own number of trials", {
  skip_if_no_bmm_model("sdt_ranking")
  sim <- simulate_recovery(
    rank_model(), c(d = 1),
    n_subjects = 2, n_trials = 2,
    trial_design = data.frame(n_trials = c(100, 30)), seed = 1
  )
  expect_named(sim$data, c("id", "n_trials", ranks))
  expect_equal(
    unname(rowSums(as.data.frame(sim$data)[ranks])), c(100, 30, 100, 30)
  )
})

# m as a column --------------------------------------------------------------

test_that("m as a column needs a trial design that holds it", {
  skip_if_no_bmm_model("sdt_ranking")
  expect_error(
    simulate_recovery(
      rank_model("mm"), c(d = 1),
      n_subjects = 1, n_trials = 3, seed = 1
    ),
    "needs `trial_design`"
  )
  expect_error(
    simulate_recovery(
      rank_model("mm"), c(d = 1),
      n_subjects = 1, n_trials = 3,
      trial_design = data.frame(n_trials = c(10, 10, 10)), seed = 1
    ),
    "mm"
  )
})

test_that("m as a column draws each row at its own set size", {
  skip_if_no_bmm_model("sdt_ranking")
  sim <- simulate_recovery(
    rank_model("mm"), c(d = 1),
    n_subjects = 2, n_trials = 3, trial_design = rank_design(), seed = 2
  )
  expect_named(sim$data, c("id", "mm", "n_trials", ranks))
  expect_identical(sim$trial_design_columns, c("mm", "n_trials"))
  counts <- as.matrix(as.data.frame(sim$data)[ranks])
  expect_equal(unname(rowSums(counts)), rep(c(60, 40, 50), 2))
  # nothing is ranked beyond a row's set size
  expect_true(all(counts[col(counts) > sim$data$mm] == 0))
})

test_that("m as a column passes bmm's checks (mock backend)", {
  skip_if_no_bmm_model("sdt_ranking")
  skip_if_not_installed("brms")
  sim <- simulate_recovery(
    rank_model("mm"), c(d = 1),
    n_subjects = 3, n_trials = 3, trial_design = rank_design(), seed = 3
  )
  model <- rank_model("mm")
  expect_s3_class(
    mock_bmm(recovery_formula(model), sim$data, model), "bmmfit"
  )
})

test_that("a set size outside 2 to the number of rank columns is refused", {
  skip_if_no_bmm_model("sdt_ranking")
  bad <- list(c(2, 5, 3), c(1, 2, 3), c(2, 2.5, 3), c(2, NA, 3))
  for (mm in bad) {
    expect_error(
      simulate_recovery(
        rank_model("mm"), c(d = 1),
        n_subjects = 1, n_trials = 3,
        trial_design = rank_design(mm = mm), seed = 1
      ),
      "whole numbers from 2 to 4"
    )
  }
  # bmm wants the column numeric ("Set-size column 'mm' must be numeric"),
  # so a factor or text is refused here rather than by bmm after the draws
  for (mm in list(factor(c(2, 3, 4)), c("2", "3", "4"))) {
    expect_error(
      simulate_recovery(
        rank_model("mm"), c(d = 1),
        n_subjects = 1, n_trials = 3,
        trial_design = data.frame(mm = mm, n_trials = 10), seed = 1
      ),
      "whole numbers from 2 to 4"
    )
  }
  # a function's design is checked per call
  expect_error(
    simulate_recovery(
      rank_model("mm"), c(d = 1),
      n_subjects = 1, n_trials = 3,
      trial_design = function(n_trials) rank_design(mm = c(2, 6, 3)), seed = 1
    ),
    "whole numbers from 2 to 4"
  )
})

# density, range -------------------------------------------------------------

test_that("the density scores the counts with sdratio as the ratio", {
  skip_if_no_bmm_model("sdt_ranking", density = TRUE)
  model <- rank_model(dist = "normal")
  link <- c(d = 1, sdratio = log(1.5))
  dat <- simulate_recovery(
    model, link,
    n_subjects = 3, n_trials = 50, seed = 4
  )$data
  ll <- density_for(model)(natural_pars(link, model), dat, model)
  expect_equal(ll, bmm::dsdt_ranking(
    as.matrix(as.data.frame(dat)[ranks]),
    m = 4, d = 1, sdratio = 1.5, dist = "normal", log = TRUE
  ))
})

test_that("the density reads m from its column and NA beyond it as 0", {
  skip_if_no_bmm_model("sdt_ranking", density = TRUE)
  model <- rank_model("mm")
  dat <- simulate_recovery(
    model, c(d = 1),
    n_subjects = 1, n_trials = 3, trial_design = rank_design(), seed = 6
  )$data
  counts <- as.matrix(as.data.frame(dat)[ranks])
  expected <- bmm::dsdt_ranking(
    counts,
    m = c(2, 3, 4), d = 1, dist = "gumbel_min", log = TRUE
  )
  dat$k4[1:2] <- NA
  dat$k3[1] <- NA
  ll <- density_for(model)(natural_pars(c(d = 1), model), dat, model)
  expect_equal(ll, expected)
})

# fit_ml(method = "optim") refuses what bmm refuses ---------------------------

# bmm's rule (check_data.sdt_ranking()): an NA count beyond a row's m reads
# as 0, an NA within it is refused. The optim route never runs bmm's data
# check, so it has to apply the rule itself.

test_that("the optim route refuses an NA count where bmm does, only there", {
  skip_if_no_bmm_model("sdt_ranking", density = TRUE)
  skip_if_not_installed("brms")
  fixed <- simulate_recovery(
    rank_model(), c(d = 1),
    n_subjects = 3, n_trials = 60, seed = 1
  )$data
  by_col <- simulate_recovery(
    rank_model("mm"), c(d = 1),
    n_subjects = 2, n_trials = 3, trial_design = rank_design(), seed = 6
  )$data
  # rows 1 to 3 of `by_col` have m = 2, 3 and 4
  cases <- list(
    list(model = rank_model(), data = fixed, row = 1, col = "k3"),
    list(model = rank_model(), data = fixed, row = 2, col = "k4"),
    list(model = rank_model("mm"), data = by_col, row = 1, col = "k2"),
    list(model = rank_model("mm"), data = by_col, row = 2, col = "k3"),
    list(model = rank_model("mm"), data = by_col, row = 1, col = "k3"),
    list(model = rank_model("mm"), data = by_col, row = 2, col = "k4")
  )
  refused <- logical()
  for (case in cases) {
    # bmm takes the data as they are, so a refusal is about the NA
    expect_s3_class(
      mock_bmm(recovery_formula(case$model), case$data, case$model),
      "bmmfit"
    )
    dat <- as.data.frame(case$data)
    dat[case$row, case$col] <- NA
    bmm_refuses <- inherits(
      try(mock_bmm(recovery_formula(case$model), dat, case$model),
        silent = TRUE
      ),
      "try-error"
    )
    if (bmm_refuses) {
      expect_error(
        fit_ml(case$model, dat, method = "optim"),
        paste0("Not so in row ", case$row, " ")
      )
    } else {
      expect_s3_class(
        fit_ml(case$model, dat, method = "optim"), "bmmtools_ml"
      )
    }
    refused <- c(refused, bmm_refuses)
  }
  # both sides of the rule were exercised
  expect_true(any(refused) && !all(refused))
})

test_that("the optim route names every row with an NA count within m", {
  skip_if_no_bmm_model("sdt_ranking", density = TRUE)
  dat <- as.data.frame(simulate_recovery(
    rank_model(), c(d = 1),
    n_subjects = 4, n_trials = 60, seed = 2
  )$data)
  dat$k1[2] <- NA
  dat$k4[4] <- NA
  expect_error(
    fit_ml(rank_model(), dat, method = "optim"),
    "rows 2 and 4"
  )
})

# bmm refuses data missing the set-size column or a rank column and names it
# ("Set-size column 'mm' missing in the data", "Response columns 'k4'
# missing in the data"; mock backend, CRAN bmm 1.4.1, 2026-10-25).

test_that("the optim route names a missing set-size or rank column", {
  skip_if_no_bmm_model("sdt_ranking", density = TRUE)
  by_col <- as.data.frame(simulate_recovery(
    rank_model("mm"), c(d = 1),
    n_subjects = 2, n_trials = 3, trial_design = rank_design(), seed = 6
  )$data)
  fixed <- as.data.frame(simulate_recovery(
    rank_model(), c(d = 1),
    n_subjects = 2, n_trials = 60, seed = 1
  )$data)
  no_m <- by_col[setdiff(names(by_col), "mm")]
  expect_error(
    fit_ml(rank_model("mm"), no_m, method = "optim"),
    "no column.*mm",
    class = "rlang_error"
  )
  no_rank <- by_col[setdiff(names(by_col), "k4")]
  expect_error(
    fit_ml(rank_model("mm"), no_rank, method = "optim"),
    "no column.*k4",
    class = "rlang_error"
  )
  no_rank <- fixed[setdiff(names(fixed), "k3")]
  expect_error(
    fit_ml(rank_model(), no_rank, method = "optim"),
    "no column.*k3",
    class = "rlang_error"
  )
})

test_that("the optim route reads an NA count beyond m as the 0 it stands for", {
  skip_if_no_bmm_model("sdt_ranking", density = TRUE)
  model <- rank_model("mm")
  dat <- as.data.frame(simulate_recovery(
    model, c(d = 1),
    n_subjects = 2, n_trials = 3, trial_design = rank_design(), seed = 6
  )$data)
  # the simulated counts beyond each row's m are 0
  with_na <- dat
  with_na$k4[dat$mm < 4] <- NA
  with_na$k3[dat$mm < 3] <- NA
  expect_equal(
    fit_ml(model, with_na, method = "optim")$estimate,
    fit_ml(model, dat, method = "optim")$estimate
  )
})

test_that("prior_check's range for sdt_ranking is 0 to each row's total", {
  skip_if_no_bmm_model("sdt_ranking")
  data <- data.frame(k1 = c(5, 1), k2 = c(3, NA), k3 = 2, k4 = 0)
  expect_equal(
    response_range(rank_model(), data),
    list(floor = 0, ceiling = c(10, 3))
  )
  expect_message(
    out <- response_range(rank_model(), data.frame(k1 = 1, k2 = 2)),
    "not in"
  )
  expect_true(is.na(out$ceiling))
})

test_that("prior draws of a ranking model are named after its rank columns", {
  skip_if_no_bmm_model("sdt_ranking")
  yrep <- array(0, c(5, 2, 4))
  expect_named(yrep_list(yrep, rank_model()), ranks)
  # a model with as many response columns as slices is not renamed
  ez <- bmm::ezdm(
    mean_rt = "mrt", var_rt = "vrt", n_upper = "nu", n_trials = "n"
  )
  expect_named(yrep_list(array(0, c(5, 2, 3)), ez), paste0("response_", 1:3))
})
