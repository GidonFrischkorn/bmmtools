# Tests for benchmark() and its summary, written against
# local/dev/spec-milestone-15-benchmark.md (decisions D69-D77).
#
# Every run goes through a mock fitter and a mock `.metrics`, so the
# orchestration is tested with no model, no Stan and no brms.

two_sets <- function(n = 3L) {
  lapply(seq_len(n), function(i) data.frame(y = seq_len(10L + i)))
}

mock_implementations <- function(a_leapfrog = 1000, b_leapfrog = 2000,
                                 seconds = 0, error_on = integer()) {
  list(
    a = sleeping_fitter(seconds, leapfrog = a_leapfrog, label = "a"),
    b = sleeping_fitter(seconds,
      leapfrog = b_leapfrog, label = "b",
      error_on = error_on
    )
  )
}

# order and bookkeeping --------------------------------------------------

test_that("runs are interleaved A, B, A, B (D74)", {
  out <- benchmark(
    mock_implementations(),
    data = two_sets(3L), runs = 2L, warmup_runs = 1L,
    compile = FALSE, .metrics = mock_metrics()
  )
  expect_identical(out$implementation, c("a", "b", "a", "b", "a", "b"))
  expect_identical(out$pair, c(0L, 0L, 1L, 1L, 2L, 2L))
  expect_identical(out$warmup_run, c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_identical(out$run, c(1L, 1L, 2L, 2L, 3L, 3L))
})

test_that("the number of rows is (runs + warmup_runs) x implementations", {
  out <- benchmark(
    mock_implementations(),
    data = two_sets(4L), runs = 3L, warmup_runs = 1L,
    compile = FALSE, .metrics = mock_metrics()
  )
  expect_equal(nrow(out), 8L)
  expect_equal(sum(out$warmup_run), 2L)

  none <- benchmark(
    mock_implementations(),
    data = two_sets(2L), runs = 2L, warmup_runs = 0L,
    compile = FALSE, .metrics = mock_metrics()
  )
  expect_equal(nrow(none), 4L)
  expect_false(any(none$warmup_run))
})

test_that("one data set per pair, the same for both implementations", {
  sets <- two_sets(3L)
  out <- benchmark(
    mock_implementations(),
    data = sets, runs = 2L, warmup_runs = 1L,
    compile = FALSE, .metrics = mock_metrics()
  )
  by_pair <- split(out$data_id, out$pair)
  expect_true(all(vapply(by_pair, function(x) length(unique(x)) == 1L, logical(1))))
  expect_equal(length(unique(out$data_id)), 3L)
  expect_identical(out$data_rows, as.integer(rep(vapply(sets, nrow, integer(1)), each = 2L)))
})

test_that("a single data frame is used for every pair", {
  out <- benchmark(
    mock_implementations(),
    data = data.frame(y = 1:10), runs = 2L,
    warmup_runs = 0L, compile = FALSE, .metrics = mock_metrics()
  )
  expect_equal(length(unique(out$data_id)), 1L)
  expect_true(all(out$data_rows == 10L))
})

test_that("wall time is measured around the fitter", {
  out <- benchmark(
    mock_implementations(seconds = 0.05),
    data = data.frame(y = 1:5),
    runs = 1L, warmup_runs = 0L, compile = FALSE, .metrics = mock_metrics()
  )
  expect_true(all(out$wall_seconds >= 0.05))
})

test_that("every run_info() column is present and timestamps do not go back", {
  out <- benchmark(
    mock_implementations(),
    data = two_sets(2L), runs = 2L, warmup_runs = 0L,
    compile = FALSE, .metrics = mock_metrics()
  )
  expect_true(all(names(run_info()) %in% names(out)))
  expect_false(is.unsorted(out$timestamp))
})

test_that("a benchmark writes no fit and leaves no cache behind (D73)", {
  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  out <- benchmark(
    mock_implementations(),
    data = two_sets(2L), runs = 2L, warmup_runs = 0L,
    compile = FALSE, .metrics = mock_metrics()
  )
  expect_equal(nrow(out), 4L)
  expect_length(list.files(dir,
    recursive = TRUE, all.files = TRUE,
    pattern = "\\.rds$|\\.key$"
  ), 0L)
})

# errors in a run --------------------------------------------------------

test_that("on_error = 'record' keeps the row and goes on", {
  expect_warning(
    out <- benchmark(
      mock_implementations(error_on = 2L),
      data = two_sets(2L), runs = 2L,
      warmup_runs = 0L, compile = FALSE, .metrics = mock_metrics()
    ),
    "failed"
  )
  expect_equal(nrow(out), 4L)
  failed <- out[nzchar(out$error), ]
  expect_equal(nrow(failed), 1L)
  expect_identical(failed$implementation, "b")
  expect_match(failed$error, "mock failure")
  expect_true(is.na(failed$us_per_gradient))
  expect_true(is.na(failed$wall_seconds))
})

test_that("on_error = 'stop' aborts", {
  expect_error(
    benchmark(
      mock_implementations(error_on = 1L),
      data = two_sets(2L), runs = 2L,
      warmup_runs = 0L, compile = FALSE, on_error = "stop",
      .metrics = mock_metrics()
    ),
    "mock failure"
  )
})

# compile measurement ----------------------------------------------------

test_that("a compile hook is timed once per implementation", {
  calls <- 0L
  impls <- list(
    a = list(
      fit = sleeping_fitter(leapfrog = 1000, label = "a"),
      compile = function(data, ...) {
        calls <<- calls + 1L
        "// stan code a"
      }
    ),
    b = list(
      fit = sleeping_fitter(leapfrog = 2000, label = "b"),
      compile = function(data, ...) "// stan code b"
    )
  )
  out <- benchmark(
    impls,
    data = two_sets(2L), runs = 2L, warmup_runs = 0L, compile = TRUE,
    .metrics = mock_metrics(),
    .compiler = function(code, backend, force, ...) list(exe = NA_character_)
  )
  expect_equal(calls, 1L)
  expect_true(all(out$compile_source == "forced"))
  expect_true(all(!is.na(out$compile_seconds)))
  compile_table <- attr(out, "compile")
  expect_equal(nrow(compile_table), 2L)
  expect_identical(compile_table$implementation, c("a", "b"))
})

test_that("an implementation without a compile hook says so", {
  expect_message(
    out <- benchmark(
      mock_implementations(),
      data = two_sets(2L), runs = 1L, warmup_runs = 0L,
      compile = TRUE, .metrics = mock_metrics()
    ),
    "compile"
  )
  expect_true(all(is.na(out$compile_seconds)))
  expect_true(all(out$compile_source == "none"))
})

# input validation -------------------------------------------------------

test_that("benchmark() validates its arguments", {
  good_data <- two_sets(2L)
  expect_error(benchmark(list(), data = good_data), "implementations")
  expect_error(benchmark("a", data = good_data), "implementations")
  expect_error(
    benchmark(list(sleeping_fitter()), data = good_data), "names"
  )
  expect_error(
    benchmark(stats::setNames(list(sleeping_fitter(), sleeping_fitter()), c("a", "a")),
      data = good_data
    ),
    "unique"
  )
  expect_error(
    benchmark(mock_implementations(), data = 1), "data frame"
  )
  expect_error(
    benchmark(mock_implementations(), data = good_data, runs = 0L), "runs"
  )
  expect_error(
    benchmark(mock_implementations(), data = good_data, warmup_runs = -1L),
    "warmup_runs"
  )
})

test_that("a data list shorter than the pairs needed is an error", {
  expect_error(
    benchmark(mock_implementations(),
      data = two_sets(2L), runs = 3L,
      warmup_runs = 1L, compile = FALSE, .metrics = mock_metrics()
    ),
    "4"
  )
})

# the summary ------------------------------------------------------------

test_that("summary() excludes the warmup runs and reports spread beside the centre", {
  out <- benchmark(
    mock_implementations(),
    data = two_sets(4L), runs = 3L, warmup_runs = 1L,
    compile = FALSE, .metrics = mock_metrics()
  )
  s <- summary(out)
  expect_s3_class(s, "bmmtools_benchmark_summary")
  per <- s$per_implementation
  expect_true(all(per$n_runs == 3L))
  expect_true(all(c(
    "median", "mad", "min", "max", "n_runs",
    "n_not_converged", "n_error"
  ) %in% names(per)))
  expect_setequal(unique(per$implementation), c("a", "b"))
})

test_that("the ratio is paired, with an interval (D77)", {
  # b takes twice as many leapfrog steps for the same time: half the cost
  out <- benchmark(
    mock_implementations(a_leapfrog = 1000, b_leapfrog = 2000),
    data = two_sets(4L), runs = 3L, warmup_runs = 1L, compile = FALSE,
    .metrics = mock_metrics()
  )
  ratios <- summary(out)$ratios
  row <- ratios[ratios$implementation == "b" &
    ratios$metric == "us_per_gradient", ]
  expect_equal(nrow(row), 1L)
  expect_equal(row$ratio, 0.5)
  expect_lte(row$ci_low, 0.5)
  expect_gte(row$ci_high, 0.5)
  expect_equal(row$n_pairs, 3L)
  expect_identical(row$ci_method, "paired-t (log)")
  expect_identical(unique(ratios$reference), "a")
})

test_that("one pair gives no interval, and says why", {
  out <- benchmark(
    mock_implementations(),
    data = two_sets(2L), runs = 1L, warmup_runs = 1L,
    compile = FALSE, .metrics = mock_metrics()
  )
  ratios <- summary(out)$ratios
  expect_true(all(is.na(ratios$ci_low)))
  measured <- ratios[ratios$n_pairs > 0L, ]
  expect_true(nrow(measured) > 0L)
  expect_true(all(measured$ci_method == "single pair"))
  expect_true(all(measured$n_pairs == 1L))
  # a metric nothing measured says so rather than claiming a single pair
  unmeasured <- ratios[ratios$metric == "compile_seconds", ]
  expect_identical(unmeasured$ci_method, "not measured")
})

test_that("a single implementation has no ratios", {
  out <- benchmark(
    list(only = sleeping_fitter(leapfrog = 1000)),
    data = two_sets(2L),
    runs = 2L, warmup_runs = 0L, compile = FALSE, .metrics = mock_metrics()
  )
  s <- summary(out)
  expect_equal(nrow(s$ratios), 0L)
  expect_equal(nrow(s$per_implementation) > 0L, TRUE)
  expect_output(print(s), "one implementation")
})

test_that("runs that did not converge are kept, counted and warned about", {
  failing_metrics <- function(fit, wall_seconds = NA_real_, ...) {
    row <- mock_metrics()(fit, wall_seconds = wall_seconds, ...)
    if (identical(fit$label, "b")) row$converged <- FALSE
    row
  }
  out <- benchmark(
    mock_implementations(),
    data = two_sets(2L), runs = 2L, warmup_runs = 0L,
    compile = FALSE, .metrics = failing_metrics
  )
  expect_warning(s <- summary(out), "converge")
  per <- s$per_implementation
  expect_true(all(per$n_not_converged[per$implementation == "b"] == 2L))
  expect_true(all(per$n_runs[per$implementation == "b"] == 2L))
})

test_that("differing sampler settings between implementations warn", {
  settings_metrics <- function(fit, wall_seconds = NA_real_, ...) {
    row <- mock_metrics()(fit, wall_seconds = wall_seconds, ...)
    if (identical(fit$label, "b")) row$iter <- 4000L
    row
  }
  out <- benchmark(
    mock_implementations(),
    data = two_sets(2L), runs = 2L, warmup_runs = 0L,
    compile = FALSE, .metrics = settings_metrics
  )
  expect_warning(summary(out), "iter")
})

test_that("print() shows the implementations and the run count", {
  out <- benchmark(
    mock_implementations(),
    data = two_sets(2L), runs = 2L, warmup_runs = 0L,
    compile = FALSE, .metrics = mock_metrics()
  )
  expect_output(print(out), "benchmark")
  expect_output(print(summary(out)), "us_per_gradient")
})

# the remaining branches -------------------------------------------------

test_that("benchmark() validates compile, .metrics and a malformed spec", {
  good <- two_sets(2L)
  expect_error(
    benchmark(mock_implementations(), data = good, compile = NA), "compile"
  )
  expect_error(
    benchmark(mock_implementations(), data = good, .metrics = 1), ".metrics"
  )
  expect_error(
    benchmark(list(a = list(fit = "not a function")), data = good),
    "fit"
  )
})

test_that("settings the metrics did not know come from the call", {
  bare_metrics <- function(fit, wall_seconds = NA_real_, ...) {
    row <- mock_metrics()(fit, wall_seconds = wall_seconds, ...)
    row$backend <- NA_character_
    row$fit_cores <- NA_integer_
    row$seed <- NA_real_
    row
  }
  out <- benchmark(
    mock_implementations(), data = two_sets(2L), runs = 2L, warmup_runs = 0L,
    compile = FALSE, .metrics = bare_metrics, backend = "rstan", cores = 3,
    seed = 42
  )
  expect_true(all(out$backend == "rstan"))
  expect_true(all(out$fit_cores == 3L))
  expect_true(all(out$seed == 42))
})

test_that("the compile hook is compiled for the backend the call names", {
  seen <- NULL
  impls <- list(
    a = list(
      fit = sleeping_fitter(leapfrog = 1000, label = "a"),
      compile = function(data, ...) "// code"
    )
  )
  out <- benchmark(
    impls, data = two_sets(2L), runs = 2L, warmup_runs = 0L, compile = TRUE,
    .metrics = mock_metrics(), backend = "rstan",
    .compiler = function(code, backend, force, ...) {
      seen <<- backend
      list(exe = NA_character_)
    }
  )
  expect_identical(seen, "rstan")
  expect_identical(attr(out, "compile")$compile_source, "forced")
})

test_that("a run whose metrics carry no per-chain table still records", {
  chainless <- function(fit, wall_seconds = NA_real_, ...) {
    row <- mock_metrics()(fit, wall_seconds = wall_seconds, ...)
    attr(row, "chains") <- NULL
    row
  }
  out <- benchmark(
    mock_implementations(), data = two_sets(2L), runs = 2L, warmup_runs = 0L,
    compile = FALSE, .metrics = chainless
  )
  expect_equal(nrow(out), 4L)
  expect_equal(nrow(attr(out, "chains")), 0L)
})

test_that("print() names the warmup runs and the failures", {
  expect_warning(
    out <- benchmark(
      mock_implementations(error_on = 2L), data = two_sets(3L), runs = 2L,
      warmup_runs = 1L, compile = FALSE, .metrics = mock_metrics()
    ),
    "failed"
  )
  expect_output(print(out), "warmup run")
  expect_output(print(out), "failed")
})

test_that("summary() validates its metrics, reference and level", {
  out <- benchmark(
    mock_implementations(), data = two_sets(2L), runs = 2L, warmup_runs = 0L,
    compile = FALSE, .metrics = mock_metrics()
  )
  expect_error(summary(out, metrics = "not_a_column"), "metrics")
  expect_error(summary(out, reference = "c"), "reference")
  expect_error(summary(out, ci_level = 0), "ci_level")
  expect_error(summary(out, ci_level = 1.5), "ci_level")
})

test_that("a verb that drops a contract column drops the class", {
  out <- benchmark(
    mock_implementations(), data = two_sets(2L), runs = 2L, warmup_runs = 0L,
    compile = FALSE, .metrics = mock_metrics()
  )
  expect_s3_class(out[1:2, ], "bmmtools_benchmark")
  expect_null(attr(out[1:2, ], "chains"))
  dropped <- dplyr::select(out, "implementation", "us_per_gradient")
  expect_false(inherits(dropped, "bmmtools_benchmark"))
  expect_s3_class(dropped, "tbl_df")
  kept <- dplyr::filter(out, .data$implementation == "a")
  expect_s3_class(kept, "bmmtools_benchmark")
})

# what the review of 15.2 found ------------------------------------------

test_that("compile time is one measurement, not one per pair", {
  impls <- list(
    a = list(
      fit = sleeping_fitter(leapfrog = 1000, label = "a"),
      compile = function(data, ...) "// a"
    ),
    b = list(
      fit = sleeping_fitter(leapfrog = 2000, label = "b"),
      compile = function(data, ...) "// b"
    )
  )
  calls <- 0L
  out <- benchmark(
    impls, data = two_sets(4L), runs = 3L, warmup_runs = 1L, compile = TRUE,
    .metrics = mock_metrics(),
    .compiler = function(code, backend, force, ...) {
      calls <<- calls + 1L
      Sys.sleep(0.01 * calls)
      list(exe = NA_character_)
    }
  )
  s <- summary(out)
  compile_rows <- s$per_implementation[
    s$per_implementation$metric == "compile_seconds",
  ]
  expect_true(all(compile_rows$n_measured == 1L))
  # n_runs still counts the recorded runs; one value has no spread
  expect_true(all(compile_rows$n_runs == 3L))
  expect_true(all(is.na(compile_rows$mad)))
  ratio <- s$ratios[s$ratios$metric == "compile_seconds", ]
  expect_equal(ratio$n_pairs, 1L)
  expect_identical(ratio$ci_method, "single measurement")
  expect_true(is.na(ratio$ci_low))
  expect_true(is.finite(ratio$ratio))
  # a per-run metric still gets its paired interval
  per_run <- s$ratios[s$ratios$metric == "us_per_gradient", ]
  expect_equal(per_run$n_pairs, 3L)
})

test_that("an interval widens with the spread between pairs", {
  varying <- local({
    seen <- 0L
    function(fit, wall_seconds = NA_real_, ...) {
      seen <<- seen + 1L
      scale <- if (identical(fit$label, "b")) c(1, 2, 4)[[(seen + 1L) %/% 2L]] else 1
      row <- mock_metrics()(fit, wall_seconds = wall_seconds, ...)
      row$us_per_gradient <- row$us_per_gradient * scale
      row
    }
  })
  out <- benchmark(
    mock_implementations(), data = two_sets(3L), runs = 3L, warmup_runs = 0L,
    compile = FALSE, .metrics = varying
  )
  ratio <- summary(out)$ratios
  ratio <- ratio[ratio$metric == "us_per_gradient", ]
  expect_equal(ratio$n_pairs, 3L)
  expect_lt(ratio$ci_low, ratio$ratio)
  expect_gt(ratio$ci_high, ratio$ratio)
})

test_that("a failure in the measuring step is recorded like any other", {
  broken <- function(fit, ...) {
    if (identical(fit$label, "b")) stop("cannot measure this fit", call. = FALSE)
    mock_metrics()(fit, ...)
  }
  expect_warning(
    out <- benchmark(
      mock_implementations(), data = two_sets(2L), runs = 2L, warmup_runs = 0L,
      compile = FALSE, .metrics = broken
    ),
    "failed"
  )
  failed <- out[nzchar(out$error), ]
  expect_equal(nrow(failed), 2L)
  expect_match(failed$error[[1L]], "cannot measure")
  expect_true(all(!is.na(failed$timestamp)))

  expect_error(
    benchmark(
      mock_implementations(), data = two_sets(2L), runs = 2L, warmup_runs = 0L,
      compile = FALSE, on_error = "stop", .metrics = broken
    ),
    "Measuring"
  )
})

test_that("a failed run still carries the run_info() columns", {
  expect_warning(
    out <- benchmark(
      mock_implementations(error_on = 1L), data = two_sets(2L), runs = 2L,
      warmup_runs = 0L, compile = FALSE, .metrics = mock_metrics()
    ),
    "failed"
  )
  expect_true(all(names(run_info()) %in% names(out)))
  expect_true(all(!is.na(out$timestamp)))
  expect_true(all(!is.na(out$platform)))
})

test_that("the fitter is never handed a cache file (D73)", {
  seen <- list()
  recording <- function(data, ...) {
    seen[[length(seen) + 1L]] <<- names(rlang::list2(...))
    structure(list(label = "a", leapfrog = 1000), class = "mock_bench_fit")
  }
  out <- benchmark(
    list(a = recording), data = two_sets(2L), runs = 2L, warmup_runs = 0L,
    compile = FALSE, .metrics = mock_metrics(), chains = 2
  )
  expect_equal(nrow(out), 2L)
  passed <- unique(unlist(seen))
  expect_false(any(c("file", "file_refit", "file_compress") %in% passed))
  expect_true("chains" %in% passed)
})
