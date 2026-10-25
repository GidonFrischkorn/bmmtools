# run_validation(), card_numbers(), verdict_table(), compare_validations()
# (Milestone 26, spec 26 § Functions). Every stage runs through the
# `.fitter` mock of helper-validation.R; nothing is compiled. The mock's
# draws are sequences, so its SBC diagnostics always fail: a run against
# it is "complete, not validated" unless a test edits the result.

skip_validation_run <- function() {
  skip_if_not_installed("bmm")
  skip_if_not_installed("SBC")
  skip_on_cran() # every run fits through SBC: about 0.5 s each, 40 runs
}

#' Which stage a mock call belongs to, by its design point
#'
#' small_validation_plan() gives the prior check 3 subjects, the
#' calibration 4, and the grid 5 and 6.
#' @noRd
call_stage <- function(calls) {
  vapply(calls$log, function(entry) {
    n <- length(entry$ids)
    switch(as.character(n),
      "3" = "prior_predictive",
      "4" = "sbc",
      "recovery"
    )
  }, character(1))
}

stage_status_of <- function(result) {
  stats::setNames(result$stages$status, result$stages$stage)
}

# a full run --------------------------------------------------------------

test_that("a run gives every stage a status and keeps every result", {
  skip_validation_run()
  dir <- withr::local_tempdir()
  mock <- validation_mock_fitter()
  result <- mock_validation(small_validation_plan(), dir, mock)

  expect_s3_class(result, "bmmtools_validation")
  status <- stage_status_of(result)
  ran <- c(
    "run_record", "prior_predictive", "sbc", "recovery", "computation",
    "identifiability", "individual_differences"
  )
  expect_identical(unname(status[ran]), rep("ok", length(ran)))
  expect_identical(status[["data_requirements"]], "not built")
  expect_identical(status[["agreement"]], "not applicable")
  expect_identical(status[["model_recovery"]], "not built")
  expect_identical(status[["posterior_checks"]], "not built")
  expect_identical(
    result$stages$message[result$stages$stage == "data_requirements"],
    "not run (#17)"
  )

  # each stage's object is the stage function's, unchanged
  expect_s3_class(result$prior_predictive, "bmmtools_prior_check")
  expect_s3_class(result$sbc, "SBC_results")
  expect_s3_class(result$recovery, "bmmtools_recovery")
  expect_null(result$data_requirements)
  expect_null(result$agreement)
  expect_identical(result$run_record$bmm_version, result$run_info$bmm_version[[1L]])
})

test_that("each stage fits what the plan says, once", {
  skip_validation_run()
  mock <- validation_mock_fitter()
  mock_validation(small_validation_plan(), withr::local_tempdir(), mock)

  stage <- call_stage(mock$calls)
  kind <- mock_calls_by_kind(mock$calls)
  # two prior sets: bmm's defaults and the calibration's own prior
  expect_identical(sum(stage == "prior_predictive" & kind == "prior"), 2L)
  expect_identical(sum(stage == "prior_predictive" & kind == "fit"), 0L)
  # one prior fit and one fit per simulation
  expect_identical(sum(stage == "sbc" & kind == "prior"), 1L)
  expect_identical(sum(stage == "sbc" & kind == "fit"), 5L)
  # 2 designs x 2 replications, the preflight, one prior fit for prior_sd
  expect_identical(sum(stage == "recovery" & kind == "fit"), 5L)
  expect_identical(sum(stage == "recovery" & kind == "prior"), 1L)

  priors <- lapply(mock$calls$log[stage == "prior_predictive"], `[[`, "prior")
  expect_null(priors[[1L]])
  expect_s3_class(priors[[2L]], "brmsprior")
  expect_identical(
    priors[[2L]]$prior,
    c("normal(2, 0.5)", "normal(0, 0.3)", "normal(0, 0.5)")
  )
})

test_that("the plan's samplers reach their stages and the call's win", {
  skip_validation_run()
  mock <- validation_mock_fitter()
  plan <- small_validation_plan()
  result <- mock_validation(plan, withr::local_tempdir(), mock,
    iter = 123L, cores = 2L
  )
  stage <- call_stage(mock$calls)
  kind <- mock_calls_by_kind(mock$calls)
  cells <- mock$calls$log[stage == "recovery" & kind == "fit"]
  # the preflight runs one chain of 200 iterations by design
  cells <- Filter(function(entry) !isTRUE(entry$chains == 1), cells)
  expect_length(cells, 4L)
  expect_true(all(vapply(cells, `[[`, 0, "iter") == 123))
  expect_true(all(vapply(cells, `[[`, 0, "chains") == 4))
  prior_check_fits <- mock$calls$log[stage == "prior_predictive"]
  expect_true(all(vapply(prior_check_fits, `[[`, 0, "chains") == 2))
  expect_true(all(vapply(prior_check_fits, `[[`, 0, "iter") == 123))
  # sbc() refuses `cores` in its dots; the runner keeps it out
  expect_identical(stage_status_of(result)[["sbc"]], "ok")
})

test_that("every stage writes its status, and the run its record", {
  skip_validation_run()
  dir <- withr::local_tempdir()
  result <- mock_validation(small_validation_plan(), dir, validation_mock_fitter())

  expect_true(file.exists(file.path(dir, "plan.rds")))
  expect_true(file.exists(file.path(dir, "result.rds")))
  for (stage in validation_stage_names()) {
    paths <- validation_stage_paths(dir, stage)
    expect_true(file.exists(paths$status), label = stage)
    status <- readRDS(paths$status)
    expect_identical(status$status, stage_status_of(result)[[stage]])
    expect_identical(
      file.exists(paths$result), identical(status$status, "ok"),
      label = stage
    )
  }
  info <- readRDS(file.path(dir, "run-info.rds"))
  expect_identical(info$when, c("start", "end"))
  expect_true(info$seconds[[2L]] >= 0)
  expect_identical(result$run_info$when, c("start", "end"))
  stored <- readRDS(file.path(dir, "result.rds"))
  expect_identical(stored$headline, result$headline)
  expect_identical(readRDS(file.path(dir, "plan.rds"))$hash, result$plan$hash)
})

# resuming ----------------------------------------------------------------

test_that("a second run reads every stage and fits nothing", {
  skip_validation_run()
  dir <- withr::local_tempdir()
  plan <- small_validation_plan()
  first <- mock_validation(plan, dir, validation_mock_fitter())
  mock <- validation_mock_fitter()
  second <- mock_validation(plan, dir, mock)

  expect_identical(mock$calls$n, 0L)
  ran <- second$stages$stage[second$stages$status == "ok"]
  expect_true(all(second$stages$resumed[second$stages$stage %in% setdiff(ran, "run_record")]))
  expect_false(second$stages$resumed[second$stages$stage == "run_record"])
  expect_identical(
    second$numbers[second$numbers$metric != "seconds", ],
    first$numbers[first$numbers$metric != "seconds", ]
  )
})

test_that("a changed plan reruns only the stages whose inputs changed", {
  skip_validation_run()
  dir <- withr::local_tempdir()
  mock_validation(small_validation_plan(), dir, validation_mock_fitter())
  mock <- validation_mock_fitter()
  expect_message(
    second <- run_validation(
      small_validation_plan(reps = 3L), dir,
      .fitter = mock$fitter
    ) |>
      withCallingHandlers(
        bmmtools_sbc_diagnostics = function(w) invokeRestart("muffleWarning")
      ),
    "The plan differs"
  ) |>
    suppressMessages() # the per-stage progress messages after the one tested

  resumed <- stats::setNames(second$stages$resumed, second$stages$stage)
  expect_true(resumed[["prior_predictive"]])
  expect_true(resumed[["sbc"]])
  expect_false(resumed[["recovery"]])
  # the views read the grid, so they follow it
  expect_false(resumed[["computation"]])
  # the grid keeps its cells: only the third replication of each design
  expect_identical(mock$calls$n, 2L)
  expect_identical(unique(attr(second$recovery, "cells")$replication), 1:3)
})

test_that("a changed criterion reruns nothing and changes the verdict", {
  skip_validation_run()
  dir <- withr::local_tempdir()
  mock_validation(small_validation_plan(), dir, validation_mock_fitter())
  mock <- validation_mock_fitter()
  strict <- mock_validation(
    small_validation_plan(criteria = list(contraction_min = 0.9)), dir, mock
  )
  expect_identical(mock$calls$n, 0L)
  verdicts <- strict$verdicts
  row <- grepl("^contraction above 0.9", verdicts$criterion)
  expect_identical(verdicts$pass[row], FALSE)
})

test_that("the call's arguments are inputs: a new iter reruns the fits", {
  skip_validation_run()
  dir <- withr::local_tempdir()
  plan <- small_validation_plan()
  mock_validation(plan, dir, validation_mock_fitter())
  mock <- validation_mock_fitter()
  again <- mock_validation(plan, dir, mock, iter = 300L)
  expect_false(any(again$stages$resumed[again$stages$stage %in% c("sbc", "recovery")]))
  expect_gt(mock$calls$n, 0L)
})

test_that("a new bmm version reruns every stage", {
  skip_validation_run()
  dir <- withr::local_tempdir()
  plan <- small_validation_plan()
  mock_validation(plan, dir, validation_mock_fitter())
  real <- package_version_string
  local_mocked_bindings(package_version_string = function(pkg) {
    if (identical(pkg, "bmm")) "99.0.0" else real(pkg)
  })
  mock <- validation_mock_fitter()
  again <- mock_validation(plan, dir, mock)
  ran <- again$stages$status == "ok" & again$stages$stage != "run_record"
  expect_false(any(again$stages$resumed[ran]))
  # the fit caches key on the version too, so the grid's cells refit
  expect_gt(sum(call_stage(mock$calls) == "recovery"), 0L)
})

test_that("an error after the stages is reported and the result written", {
  skip_validation_run()
  dir <- withr::local_tempdir()
  local_mocked_bindings(card_numbers = function(result) stop("boom"))
  run <- collect_warnings(
    mock_validation(small_validation_plan(), dir, validation_mock_fitter())
  )
  expect_true(any(grepl("could not be computed", run$warnings)))
  expect_identical(run$value$headline$state, "incomplete")
  expect_match(run$value$headline$text, "boom")
  expect_true(file.exists(file.path(dir, "result.rds")))
  expect_identical(stage_status_of(run$value)[["recovery"]], "ok")
})

# failures ----------------------------------------------------------------

test_that("a stage that fails is recorded and the others still run", {
  skip_validation_run()
  dir <- withr::local_tempdir()
  # every prior-only fit at the prior check's design point fails
  fails <- function(data, dots, prior) {
    identical(dots$sample_prior, "only") && length(unique(data$id)) == 3L
  }
  mock <- validation_mock_fitter(fail_on = fails)
  run <- collect_warnings(mock_validation(small_validation_plan(), dir, mock))
  failed <- grep("^Stage prior_predictive failed", run$warnings, value = TRUE)
  expect_length(failed, 1L)
  # the stage's own condition follows on the next line
  expect_match(failed, "the other stages go on.\nx .*mock sampler failure")
  result <- run$value
  status <- stage_status_of(result)
  expect_identical(status[["prior_predictive"]], "failed")
  expect_identical(
    unname(status[c("sbc", "recovery", "computation", "identifiability")]),
    rep("ok", 4L)
  )
  expect_match(
    result$stages$message[result$stages$stage == "prior_predictive"],
    "mock sampler failure"
  )
  expect_null(result$prior_predictive)
  expect_identical(result$headline$state, "incomplete")
  expect_match(result$headline$text, "incomplete \\(stages Prior predictive\\)")
  row <- result$verdicts[result$verdicts$stage == "prior_predictive", ]
  expect_true(is.na(row$pass))
  expect_match(row$note, "mock sampler failure")

  # fixed, a rerun runs the failed stage alone
  again <- validation_mock_fitter()
  fixed <- mock_validation(small_validation_plan(), dir, again)
  expect_identical(stage_status_of(fixed)[["prior_predictive"]], "ok")
  expect_identical(unique(call_stage(again$calls)), "prior_predictive")
  expect_true(fixed$stages$resumed[fixed$stages$stage == "sbc"])

  # a stage that fails on a rerun leaves no older result behind
  collect_warnings(mock_validation(
    small_validation_plan(), dir, validation_mock_fitter(fail_on = fails),
    refit = "always"
  ))
  paths <- validation_stage_paths(dir, "prior_predictive")
  expect_false(file.exists(paths$result))
  expect_identical(readRDS(paths$status)$status, "failed")
})

test_that("an SBC result without bmmtools' record or band is not judged", {
  skip_validation_run()
  result <- mock_validation(
    small_validation_plan(), withr::local_tempdir(), validation_mock_fitter()
  )
  attr(result$sbc, "bmmtools_sbc") <- NULL
  rows <- verdict_sbc(result, result$plan$criteria)
  expect_true(is.na(rows$pass[[2L]]))
  band <- sbc_ecdf_band(list())
  expect_false(is.data.frame(band))
  expect_true(nzchar(attr(band, "reason")))
  result$sbc <- list()
  rows <- verdict_sbc(result, result$plan$criteria)
  expect_match(rows$note[[1L]], "^not judged: no ECDF band: ")
})

test_that("a stage that reads a failed stage is skipped, not failed", {
  skip_validation_run()
  # the grid's simulation cannot be generated: that stops recovery_grid()
  plan <- small_validation_plan(generator = function(pars, n_trials, model) {
    stop("no data for you")
  })
  # the prior check and SBC simulate with the same generator
  run <- collect_warnings(
    mock_validation(plan, withr::local_tempdir(), validation_mock_fitter())
  )
  failed <- grep("^Stage .* failed", run$warnings, value = TRUE)
  expect_identical(
    sub("^Stage ([a-z_]+) failed.*$", "\\1", failed),
    c("prior_predictive", "sbc", "recovery")
  )
  result <- run$value
  status <- stage_status_of(result)
  expect_identical(status[["recovery"]], "failed")
  for (stage in c("computation", "identifiability", "individual_differences")) {
    expect_identical(status[[stage]], "skipped")
    expect_identical(
      result$stages$message[result$stages$stage == stage],
      "needs the recovery stage, which produced no result"
    )
  }
})

test_that("cells that fail lower the convergence rate and fail computation", {
  skip_validation_run()
  # both replications of the second design (6 subjects) fail to fit
  fails <- function(data, dots, prior) {
    is.null(dots$sample_prior) && length(unique(data$id)) == 6L
  }
  run <- collect_warnings(mock_validation(
    small_validation_plan(), withr::local_tempdir(),
    validation_mock_fitter(fail_on = fails)
  ))
  # the grid names its failed cells; no stage fails
  expect_true(any(grepl("cell", run$warnings)))
  expect_false(any(grepl("^Stage", run$warnings)))
  result <- run$value
  expect_identical(stage_status_of(result)[["recovery"]], "ok")
  computation <- result$computation
  expect_identical(computation$pass_rate[computation$condition == "all"], 0.5)
  expect_identical(computation$n_cells[computation$condition == "6 x 20"], 2L)
  expect_identical(computation$n_ok[computation$condition == "6 x 20"], 0L)
  row <- result$verdicts[result$verdicts$stage == "computation", ]
  expect_identical(row$value, 0.5)
  expect_false(row$pass)
  expect_identical(result$headline$state, "not validated")
  expect_match(result$headline$text, "Computation")
})

test_that("a rerun retries the failed cells and reads the rest", {
  skip_validation_run()
  dir <- withr::local_tempdir()
  fails <- function(data, dots, prior) {
    is.null(dots$sample_prior) && length(unique(data$id)) == 6L
  }
  collect_warnings(mock_validation(
    small_validation_plan(), dir, validation_mock_fitter(fail_on = fails)
  ))
  mock <- validation_mock_fitter()
  fixed <- mock_validation(small_validation_plan(), dir, mock)
  resumed <- stats::setNames(fixed$stages$resumed, fixed$stages$stage)
  expect_false(resumed[["recovery"]])
  expect_true(resumed[["sbc"]])
  # the two failed cells of the second design, nothing else
  expect_identical(mock$calls$n, 2L)
  expect_identical(fixed$computation$pass_rate[fixed$computation$condition == "all"], 1)
  # a complete grid is read on the next run
  again <- validation_mock_fitter()
  mock_validation(small_validation_plan(), dir, again)
  expect_identical(again$calls$n, 0L)
})

test_that("an unreadable stored result is run again, not trusted", {
  skip_validation_run()
  dir <- withr::local_tempdir()
  mock_validation(small_validation_plan(), dir, validation_mock_fitter())
  writeLines("not an rds file", validation_stage_paths(dir, "computation")$result)
  result <- mock_validation(small_validation_plan(), dir, validation_mock_fitter())
  expect_identical(stage_status_of(result)[["computation"]], "ok")
  expect_false(result$stages$resumed[result$stages$stage == "computation"])
  expect_s3_class(result$computation, "tbl_df")
})

# refit -------------------------------------------------------------------

test_that("refit = 'never' reads what is stored and skips the rest", {
  skip_validation_run()
  empty <- withr::local_tempdir()
  mock <- validation_mock_fitter()
  none <- mock_validation(small_validation_plan(), empty, mock, refit = "never")
  expect_identical(mock$calls$n, 0L)
  status <- stage_status_of(none)
  expect_identical(status[["run_record"]], "ok")
  expect_identical(
    unname(status[c("prior_predictive", "sbc", "recovery", "computation")]),
    rep("skipped", 4L)
  )
  expect_identical(none$headline$state, "incomplete")

  dir <- withr::local_tempdir()
  first <- mock_validation(small_validation_plan(), dir, validation_mock_fitter())
  same <- mock_validation(small_validation_plan(), dir, mock, refit = "never")
  expect_identical(stage_status_of(same)[["recovery"]], "ok")
  # a different plan is read as it is stored, with no refit, and marked
  # stale: it is not the new plan's result
  read <- mock_validation(small_validation_plan(reps = 5L), dir, mock,
    refit = "never"
  )
  expect_identical(mock$calls$n, 0L)
  status <- stage_status_of(read)
  expect_identical(status[["recovery"]], "stale")
  expect_identical(status[["prior_predictive"]], "ok")
  # the views read the stale grid, so their key differs too
  expect_identical(status[["computation"]], "stale")
  expect_match(
    read$stages$message[read$stages$stage == "recovery"], "other inputs"
  )
  expect_identical(max(attr(read$recovery, "cells")$replication), 2L)
  expect_identical(read$headline$state, "incomplete")
  # the stored plan is the one the results were made from
  expect_identical(readRDS(file.path(dir, "plan.rds"))$hash, first$plan$hash)
})

test_that("refit = 'always' reruns every stage and refits the prior fits", {
  skip_validation_run()
  dir <- withr::local_tempdir()
  plan <- small_validation_plan()
  mock_validation(plan, dir, validation_mock_fitter())
  mock <- validation_mock_fitter()
  again <- mock_validation(plan, dir, mock, refit = "always")
  expect_false(any(again$stages$resumed))
  stage <- call_stage(mock$calls)
  kind <- mock_calls_by_kind(mock$calls)
  expect_identical(sum(stage == "prior_predictive" & kind == "prior"), 2L)
  expect_identical(sum(stage == "sbc" & kind == "prior"), 1L)
  # the grid's cells are kept by their cache keys
  expect_identical(sum(stage == "recovery" & kind == "fit"), 0L)
})

# applicability and agreement --------------------------------------------

test_that("one subject-varying parameter runs no correlations", {
  skip_validation_run()
  result <- mock_validation(
    small_validation_plan(sds = c(kappa = 0.3)), withr::local_tempdir(),
    validation_mock_fitter()
  )
  status <- stage_status_of(result)
  expect_identical(status[["individual_differences"]], "not applicable")
  expect_null(attr(result$recovery, "correlations"))
  expect_null(result$individual_differences)
})

test_that("excluded stages leave the card incomplete", {
  skip_validation_run()
  mock <- validation_mock_fitter()
  result <- mock_validation(
    small_validation_plan(stages = "prior_predictive"),
    withr::local_tempdir(), mock
  )
  status <- stage_status_of(result)
  expect_identical(status[["prior_predictive"]], "ok")
  expect_identical(status[["sbc"]], "excluded")
  expect_identical(unique(call_stage(mock$calls)), "prior_predictive")
  expect_identical(result$headline$state, "incomplete")
  expect_match(result$headline$text, "Implementation \\(SBC\\)")
})

reference_table <- function() {
  tibble::tibble(
    term = c("kappa", "thetat"), estimate = c(0.1, -0.1),
    package = "mixtur", version = "1.2.3", call = "mixtur::fit_mixtur()",
    data_hash = "abc", date = "2026-10-24", source = "mixtur"
  )
}

test_that("agreement compares the largest design's first fit", {
  skip_validation_run()
  result <- mock_validation(
    small_validation_plan(reference = reference_table()),
    withr::local_tempdir(), validation_mock_fitter()
  )
  expect_identical(stage_status_of(result)[["agreement"]], "ok")
  expect_s3_class(result$agreement, "bmmtools_cross_check")
  expect_setequal(result$agreement$term, c("kappa", "thetat"))
  row <- result$verdicts[result$verdicts$stage == "agreement", ]
  expect_false(row$required)
  expect_true(is.na(row$pass))
  # the mock's estimates are 0 on the link scale
  expect_equal(row$value, 0.1)
  numbers <- result$numbers[
    result$numbers$stage == "agreement" & result$numbers$metric != "seconds",
  ]
  expect_setequal(unique(numbers$metric), c("estimate", "reference", "bias", "covered"))
})

test_that("agreement fails, alone, when the fit file is gone", {
  skip_validation_run()
  dir <- withr::local_tempdir()
  plan <- small_validation_plan(reference = reference_table())
  # run the grid, then delete its fits before agreement reads one
  mock_validation(
    validation_plan(
      plan$model,
      grid = plan$grid, reps = plan$reps, sbc = plan$sbc,
      prior_check = plan$prior_check, seed = plan$seed,
      stages = "recovery"
    ),
    dir, validation_mock_fitter()
  )
  cells <- list.files(
    file.path(dir, "recovery", "cells"), "^cell-[0-9]+-rep-[0-9]+\\.rds$",
    full.names = TRUE
  )
  expect_length(cells, 4L)
  unlink(cells)
  run <- collect_warnings(mock_validation(plan, dir, validation_mock_fitter()))
  expect_true(any(grepl("^Stage agreement failed", run$warnings)))
  result <- run$value
  expect_identical(stage_status_of(result)[["agreement"]], "failed")
  expect_match(
    result$stages$message[result$stages$stage == "agreement"], "is gone"
  )
  expect_identical(stage_status_of(result)[["recovery"]], "ok")
})

test_that("the ezdm entry runs through the grid with its correlation", {
  skip_validation_run()
  plan <- validation_plan(
    "ezdm",
    grid = data.frame(n_subjects = 5L, n_trials = 10L), reps = 2L,
    stages = c("recovery", "computation", "individual_differences"),
    seed = 3
  )
  mock <- validation_mock_fitter()
  result <- mock_validation(plan, withr::local_tempdir(), mock)
  status <- stage_status_of(result)
  expect_identical(status[["recovery"]], "ok")
  expect_identical(status[["individual_differences"]], "ok")
  expect_setequal(
    unique(result$recovery$term[result$recovery$level == "population"]),
    c("drift", "bound", "ndt")
  )
  truth <- result$individual_differences
  expect_equal(
    truth$true_value[truth$term == "bound__drift"], 0.4
  )
})

# card_numbers -------------------------------------------------------------

test_that("card_numbers is one long tibble over every stage that ran", {
  skip_validation_run()
  result <- mock_validation(
    small_validation_plan(), withr::local_tempdir(), validation_mock_fitter()
  )
  numbers <- card_numbers(result)
  expect_identical(
    names(numbers),
    c("stage", "condition", "level", "term", "metric", "value", "mcse")
  )
  expect_identical(numbers, result$numbers)
  expect_setequal(
    unique(numbers$stage),
    c(
      "run_record", "prior_predictive", "sbc", "recovery", "computation",
      "identifiability", "individual_differences"
    )
  )
  keys <- paste(
    numbers$stage, numbers$condition, numbers$level,
    numbers$term, numbers$metric
  )
  expect_false(anyDuplicated(keys) > 0L)
  # designs are labelled by their size, not by the grid's row number
  expect_setequal(
    unique(numbers$condition[
      numbers$stage == "recovery" & numbers$metric != "seconds"
    ]),
    c("5 x 10", "6 x 20")
  )
  contraction <- numbers$value[
    numbers$stage == "identifiability" & numbers$metric == "contraction" &
      numbers$level == "population"
  ]
  expect_equal(unique(contraction), 0.75)
  rate <- numbers[numbers$stage == "computation" & numbers$condition == "all" &
    numbers$metric == "pass_rate", ]
  expect_identical(rate$value, 1)
  # four cells all converged: the rate (4 + 1) / (4 + 2) gives the error
  expect_equal(rate$mcse, sqrt(5 / 6 * 1 / 6 / 4))
  bias <- numbers[numbers$stage == "recovery" & numbers$metric == "bias", ]
  expect_false(anyNA(bias$mcse))
  floor <- numbers[numbers$stage == "prior_predictive" &
    numbers$metric == "floor_rate", ]
  expect_identical(floor$condition, c("default", "sbc"))
  q50 <- numbers[numbers$stage == "prior_predictive" & numbers$metric == "q50", ]
  expect_true(all(is.na(q50$mcse)))
  expect_error(card_numbers(list()), "validation result")
})

test_that("the MCSE helpers follow their closed forms", {
  expect_equal(rate_mcse(0.9, 10), sqrt(10 / 12 * 2 / 12 / 10))
  expect_gt(rate_mcse(1, 10), 0)
  expect_gt(rate_mcse(0, 10), 0)
  expect_true(is.na(rate_mcse(1.5, 10)))
  expect_true(is.na(rate_mcse(NA_real_, 10)))
  # an RMSE of 0.5 around a bias of 0.3 over 5 rows
  expect_equal(bias_mcse(0.3, 0.5, 5), sqrt((0.25 - 0.09) / 4))
  expect_true(is.na(bias_mcse(0.3, 0.5, 1)))
})

# verdicts and the headline -----------------------------------------------

#' A mock result in which SBC passes, so that the headline can say
#' "validated" when nothing else fails
passing_result <- function(dir) {
  result <- mock_validation(small_validation_plan(), dir, validation_mock_fitter())
  meta <- attr(result$sbc, "bmmtools_sbc")
  meta$diagnostics[c("n_high_rhat", "n_low_ess_to_rank", "n_missing")] <- 0L
  attr(result$sbc, "bmmtools_sbc") <- meta
  result$verdicts <- verdict_table(result)
  result
}

test_that("a run in which every criterion holds is validated at its plan", {
  skip_validation_run()
  result <- passing_result(withr::local_tempdir())
  required <- result$verdicts$required & result$verdicts$status == "ok"
  expect_true(all(result$verdicts$pass[required]))
  headline <- validation_headline(result)
  expect_identical(headline$state, "validated")
  expect_identical(
    headline$text,
    sprintf("validated at plan 0.1.0 (hash %s)", substr(result$plan$hash, 1, 8))
  )
  # the not-built data requirements stage does not count against it
  expect_false("data_requirements" %in% headline$stages)
})

test_that("the mock's SBC diagnostics leave it complete, not validated", {
  skip_validation_run()
  result <- mock_validation(
    small_validation_plan(), withr::local_tempdir(), validation_mock_fitter()
  )
  expect_identical(result$headline$state, "not validated")
  expect_identical(result$headline$stages, "sbc")
  expect_identical(
    result$headline$text, "complete, not validated (stages Implementation (SBC))"
  )
  sbc <- result$verdicts[result$verdicts$stage == "sbc", ]
  expect_identical(nrow(sbc), 2L)
  expect_identical(sbc$value[[2L]], 10)
  expect_false(sbc$pass[[2L]])
})

test_that("each criterion fails on a result that misses it", {
  skip_validation_run()
  result <- passing_result(withr::local_tempdir())

  high <- result
  pc <- high$prior_predictive
  pc$value[pc$prior == "sbc" & pc$statistic == "ceiling_rate"] <- 0.2
  high$prior_predictive <- pc
  row <- verdict_table(high)[verdict_table(high)$stage == "prior_predictive", ]
  expect_identical(row$value, 0.2)
  expect_false(row$pass)
  expect_match(row$criterion, "under the calibration's prior at most 0.1")

  low <- result
  rec <- low$recovery
  out <- rec$level == "population" & rec$term == "kappa" &
    rec$condition == "row-2"
  # coverage is computed from the interval: move it off the truth
  rec$ci_low[out] <- rec$true_value[out] + 1
  rec$ci_high[out] <- rec$true_value[out] + 2
  low$recovery <- rec
  row <- verdict_table(low)[verdict_table(low)$stage == "recovery", ]
  expect_identical(row$value, 0)
  expect_false(row$pass)
  expect_identical(row$note, "kappa at 6 x 20")

  flat <- result
  flat$identifiability$contraction <- -0.1
  rows <- verdict_table(flat)[verdict_table(flat)$stage == "identifiability", ]
  expect_identical(rows$pass, c(FALSE, TRUE))
  expect_identical(validation_headline(flat)$state, "not validated")
})

test_that("a criterion that cannot be evaluated is named, not failed", {
  skip_validation_run()
  result <- passing_result(withr::local_tempdir())
  pc <- result$prior_predictive
  pc$value[pc$statistic %in% c("floor_rate", "ceiling_rate")] <- NA_real_
  result$prior_predictive <- pc
  result$verdicts <- verdict_table(result)
  row <- result$verdicts[result$verdicts$stage == "prior_predictive", ]
  expect_true(is.na(row$pass))
  expect_match(row$note, "no floor or ceiling")
  headline <- validation_headline(result)
  # D80: "validated" needs every required criterion met, so one that
  # could not be evaluated withholds it
  expect_identical(headline$state, "not judged")
  expect_identical(headline$stages, "prior_predictive")
  expect_identical(headline$text, "complete, not judged (stages Prior predictive)")
})

test_that("the coverage bar is nominal minus one MCSE at nominal", {
  criteria <- default_validation_criteria()
  rows <- tibble::tibble(
    term = c("a", "b"), condition = "x",
    coverage = c(0.9, 0.85), n_replications = 10L
  )
  # 0.95 - sqrt(0.95 * 0.05 / 10) = 0.881
  out <- coverage_verdict("text", rows, "condition", criteria, "none")
  expect_false(out$pass)
  expect_identical(out$note, "b at x")
  expect_identical(out$value, 0.85)
  rows$coverage[[2L]] <- 0.89
  expect_true(coverage_verdict("text", rows, "condition", criteria, "none")$pass)
  empty <- coverage_verdict("text", rows[0, ], "condition", criteria, "none")
  expect_true(is.na(empty$pass))
  expect_identical(empty$note, "none")
})

# compare_validations ------------------------------------------------------

test_that("a run compared with itself moved nowhere", {
  skip_validation_run()
  dir <- withr::local_tempdir()
  result <- mock_validation(small_validation_plan(), dir, validation_mock_fitter())
  changes <- compare_validations(result, result)
  expect_identical(
    names(changes),
    c(
      "stage", "condition", "level", "term", "metric", "old", "new", "diff",
      "mcse_old", "mcse_new", "moved"
    )
  )
  expect_false(any(changes$moved, na.rm = TRUE))
  expect_identical(nrow(changes), nrow(result$numbers))
  expect_identical(attr(changes, "old_bmm"), attr(changes, "new_bmm"))
  # a path and a directory are read as the result
  from_dir <- compare_validations(dir, file.path(dir, "result.rds"))
  expect_identical(from_dir$moved, changes$moved)
})

test_that("moved uses the MCSE where both runs have one, else a tolerance", {
  numbers <- tibble::tibble(
    stage = "recovery", condition = "20 x 30", level = "population",
    term = c("kappa", "kappa", "kappa", "kappa"),
    metric = c("coverage", "rmse", "seconds", "bias"),
    value = c(0.9, 0.2, 10, 0.05), mcse = c(0.095, NA, NA, 0.01)
  )
  later <- numbers
  later$value <- c(0.8, 0.21, 99, 0.2)
  changes <- compare_validations(numbers, later)
  moved <- stats::setNames(changes$moved, changes$metric)
  # 0.1 < 2 * sqrt(2) * 0.095
  expect_false(moved[["coverage"]])
  # no MCSE and the default tolerance of 0: any change moves
  expect_true(moved[["rmse"]])
  # timings never move
  expect_false(moved[["seconds"]])
  # 0.15 > 2 * sqrt(2) * 0.01
  expect_true(moved[["bias"]])
  expect_equal(changes$diff, later$value - numbers$value)

  added <- dplyr::bind_rows(later, tibble::tibble(
    stage = "sbc", condition = NA, level = NA, term = NA,
    metric = "n_fits", value = 100, mcse = NA
  ))
  changes <- compare_validations(numbers, added)
  expect_true(is.na(changes$moved[changes$metric == "n_fits"]))
  expect_identical(attr(changes, "n_compared"), 4L)
  expect_error(compare_validations(numbers[1:3], later), "lacks the column")
  expect_error(compare_validations(1, later), "validation result")
})

test_that("previous makes the changes of a run", {
  skip_validation_run()
  dir <- withr::local_tempdir()
  plan <- small_validation_plan()
  mock_validation(plan, dir, validation_mock_fitter())
  later <- mock_validation(plan, dir, validation_mock_fitter(),
    previous = file.path(dir, "result.rds")
  )
  expect_false(is.null(later$changes))
  expect_false(any(later$changes$moved, na.rm = TRUE))
  expect_identical(attr(later$changes, "new_bmm"), as.character(utils::packageVersion("bmm")))
  first <- mock_validation(plan, withr::local_tempdir(), validation_mock_fitter())
  expect_null(first$changes)
})

# arguments and printing --------------------------------------------------

test_that("run_validation checks its arguments", {
  skip_validation_run()
  plan <- small_validation_plan()
  dir <- withr::local_tempdir()
  expect_error(run_validation(list(), dir), "must be a validation plan")
  expect_error(run_validation(plan, c("a", "b")), "single path")
  expect_error(run_validation(plan, dir, workers = 2), "reserved")
  expect_error(run_validation(plan, dir, .fitter = 1), "must be a function")
  expect_error(
    run_validation(plan, dir, NULL, "on_change", NULL, 3),
    "must be named"
  )
  expect_error(
    run_validation(plan, dir, previous = file.path(dir, "none.rds")),
    "does not exist"
  )
  other <- withr::local_tempfile(fileext = ".rds")
  saveRDS(list(), other)
  expect_error(run_validation(plan, dir, previous = other), "not a validation result")
})

test_that("print shows the headline and every stage", {
  skip_validation_run()
  result <- mock_validation(
    small_validation_plan(), withr::local_tempdir(), validation_mock_fitter()
  )
  lines <- format(result)
  expect_match(lines[[1L]], "mixture2p: complete, not validated")
  expect_match(lines[[2L]], "^Plan 0.1.0 \\(hash [0-9a-f]{8}\\), bmm ")
  expect_length(lines, 2L + 11L)
  expect_true(any(grepl("data_requirements +not built .*not run \\(#17\\)", lines)))
  expect_output(print(result), "<bmmtools validation>")
})
