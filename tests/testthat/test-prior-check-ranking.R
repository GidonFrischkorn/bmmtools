# The default prior-check summary on a ranking design whose set size is
# a column. A row of sdt_ranking(m = "mm") ranks only its first `mm`
# items; the rank columns beyond are structural zeros, or NA, and are no
# prediction at all. Counted as draws they read as floor performance
# (measured 2026-10-25: k3 and k4 at a floor rate of 0.5 and a median of
# 2.5 for draws that were all 5). Every expected value below is written
# down from the cells a row actually predicts.

# A mock fit that hands back a 3-D array whole. helper-prior-check.R's
# method subsets with two indices, so it cannot carry one, and replacing
# it would change the files that rely on it; a class of its own leaves
# that method alone.
registerS3method(
  "prior_predict_draws", "mockpriorarrayfit",
  function(fit, ndraws, seed = NULL, ...) fit$yrep,
  envir = asNamespace("bmmtools")
)

array_fitter <- function(yrep) {
  mock <- prior_mock_fitter(yrep = yrep)
  function(...) {
    fit <- mock$fitter(...)
    class(fit) <- c("mockpriorarrayfit", class(fit))
    fit
  }
}

# Four rows at set sizes 2, 4, 2, 4, every row out of 10 trials, so the
# ceiling is 10 throughout.
ranking_data <- function() {
  data.frame(
    k1 = c(6, 4, 6, 4), k2 = c(4, 3, 4, 3),
    k3 = c(0, 2, 0, 2), k4 = c(0, 1, 0, 1),
    mm = c(2, 4, 2, 4)
  )
}

# Four draws x four rows x four ranks. Ranks 1 and 2 are predicted on
# every row; ranks 3 and 4 only on rows 2 and 4, where they hold
# 0, 10, 3, 5 (row 2) and 1, 2, 3, 4 (row 4) across the draws: one floor
# cell, one ceiling cell, eight cells in all. `structural` fills rows 1
# and 3 of ranks 3 and 4.
ranking_yrep <- function(structural) {
  yrep <- array(5, dim = c(4L, 4L, 4L))
  for (rank in 3:4) {
    yrep[, 2L, rank] <- c(0, 10, 3, 5)
    yrep[, 4L, rank] <- c(1, 2, 3, 4)
    yrep[, c(1L, 3L), rank] <- structural
  }
  yrep
}

predicted_cells <- c(0, 10, 3, 5, 1, 2, 3, 4)

stat <- function(out, response, statistic) {
  out$value[out$response == response & out$statistic == statistic]
}

expect_rank_beyond_m_left_out <- function(out) {
  for (rank in c("k3", "k4")) {
    expect_equal(stat(out, rank, "floor_rate"), 1 / 8)
    expect_equal(stat(out, rank, "ceiling_rate"), 1 / 8)
    expect_equal(
      c(stat(out, rank, "q50"), stat(out, rank, "q90"), stat(out, rank, "q95")),
      stats::quantile(predicted_cells, c(0.5, 0.9, 0.95), names = FALSE)
    )
  }
  # ranks every row predicts keep all sixteen cells, every one a 5
  expect_equal(stat(out, "k1", "floor_rate"), 0)
  expect_equal(stat(out, "k1", "q50"), 5)
}

test_that("structural zeros beyond a row's set size are not floor draws", {
  skip_if_no_bmm_model("sdt_ranking")
  model <- bmm::sdt_ranking(response = paste0("k", 1:4), m = "mm")

  out <- suppressMessages(prior_check(
    model, bmm::bmf(d ~ 1), ranking_data(),
    prior = fake_prior(), n_draws = 4L,
    .fitter = array_fitter(ranking_yrep(0))
  ))

  expect_rank_beyond_m_left_out(as.data.frame(out))
})

test_that("structural NA beyond a row's set size gives the same summary", {
  skip_if_no_bmm_model("sdt_ranking")
  model <- bmm::sdt_ranking(response = paste0("k", 1:4), m = "mm")
  range <- response_range(model, ranking_data())
  summarise <- default_prior_summary(range, model)

  out <- summarise(ranking_yrep(NA), ranking_data())

  expect_rank_beyond_m_left_out(out)
  expect_equal(out, summarise(ranking_yrep(0), ranking_data()))
})

test_that("a rank no row reaches has NA statistics, not NaN", {
  skip_if_no_bmm_model("sdt_ranking")
  model <- bmm::sdt_ranking(response = paste0("k", 1:4), m = "mm")
  data <- ranking_data()
  data$mm <- 3
  summarise <- default_prior_summary(c(0, 10), model)

  out <- summarise(ranking_yrep(0), data)

  k4 <- out$value[out$response == "k4"]
  expect_true(all(is.na(k4)))
  expect_false(any(is.nan(k4)))
})

test_that("a numeric set size masks nothing", {
  skip_if_no_bmm_model("sdt_ranking")
  model <- bmm::sdt_ranking(response = paste0("k", 1:4), m = 4)
  # with every rank asked of every row, rows 1 and 3 of ranks 3 and 4 are
  # predicted zeros, and they are floor draws
  yrep <- ranking_yrep(0)
  data <- ranking_data()[paste0("k", 1:4)]
  summarise <- default_prior_summary(c(0, 10), model)

  expect_no_message(out <- summarise(yrep, data))

  all_cells <- c(predicted_cells, rep(0, 8L))
  expect_equal(stat(out, "k3", "floor_rate"), mean(all_cells == 0))
  expect_equal(stat(out, "k3", "ceiling_rate"), mean(all_cells == 10))
  expect_equal(
    stat(out, "k3", "q50"),
    stats::quantile(all_cells, 0.5, names = FALSE)
  )
})

test_that("a missing set-size column is reported and gives NA, once", {
  skip_if_no_bmm_model("sdt_ranking")
  model <- bmm::sdt_ranking(response = paste0("k", 1:4), m = "mm")
  data <- ranking_data()
  data$mm <- NULL
  summarise <- default_prior_summary(c(0, 10), model)

  expect_message(out <- summarise(ranking_yrep(0), data), "mm")
  expect_true(all(is.na(out$value)))
  # prior_check() calls the summary once per prior set
  expect_no_message(summarise(ranking_yrep(0), data))
})

test_that("sdt_yn keeps every cell of its draws", {
  skip_if_no_bmm_model("sdt_yn")
  model <- bmm::sdt_yn(response = "k", stimulus = "s", n_trials = "n")
  data <- data.frame(k = 1:3, s = c(0, 1, 0), n = 10)

  out <- default_prior_summary(c(0, 10), model)(known_yrep(), data)

  # known_yrep(): twelve values 0 to 9, 10 and 10
  expect_equal(out$value, c(1 / 12, 2 / 12, 5.5, 9.9, 10))
  expect_equal(unique(out$response), "k")
})

test_that("m3 keeps its zero counts as floor draws", {
  skip_if_not_installed("bmm")
  model <- bmm::m3(
    resp_cats = c("corr", "other", "npl"), num_options = c(1, 4, 5),
    choice_rule = "simple"
  )
  data <- data.frame(corr = c(6, 8), other = c(3, 2), npl = c(1, 0))
  # two draws x two rows x three categories; `npl` is 0 in three of its
  # four cells, and all three are predictions
  yrep <- array(
    c(6, 7, 8, 9, 3, 2, 2, 1, 0, 1, 0, 0),
    dim = c(2L, 2L, 3L)
  )

  out <- default_prior_summary(response_range(model, data), model)(yrep, data)

  expect_equal(stat(out, "npl", "floor_rate"), 3 / 4)
  expect_equal(stat(out, "npl", "q50"), 0)
  expect_equal(stat(out, "corr", "ceiling_rate"), 0)
  expect_equal(stat(out, "corr", "floor_rate"), 0)
})
