# Hand-built parts for the benchmark module, so that every rule is
# tested without a fit and without Stan (spec 15, § "Tests to write").

#' A chains x {warmup, sample} elapsed-time matrix, as
#' `rstan::get_elapsed_time()` returns it
#' @noRd
fake_elapsed <- function(warmup, sample) {
  n <- max(length(warmup), length(sample))
  out <- cbind(
    warmup = as.double(rep_len(warmup, n)),
    sample = as.double(rep_len(sample, n))
  )
  rownames(out) <- paste0("chain:", seq_len(n))
  out
}

#' Per-chain leapfrog counts, as the module reduces `nuts_params()` to
#' @noRd
fake_leapfrog <- function(leapfrog, draws = 500L) {
  tibble::tibble(
    chain = seq_along(leapfrog),
    leapfrog = as.double(leapfrog),
    draws = as.integer(rep_len(draws, length(leapfrog)))
  )
}

#' A `summarise_draws()`-shaped table; `NA` rhat marks a constant
#' @noRd
fake_bench_diagnostics <- function(rhat = c(1.00, 1.01),
                                   ess_bulk = c(800, 600),
                                   ess_tail = c(700, 500),
                                   variable = NULL) {
  variable <- variable %||% paste0("b_p", seq_along(rhat))
  tibble::tibble(
    variable = variable, rhat = as.double(rhat),
    ess_bulk = as.double(ess_bulk), ess_tail = as.double(ess_tail)
  )
}

#' The settings a run records, with sensible stand-ins
#' @noRd
fake_bench_settings <- function(...) {
  utils::modifyList(
    list(
      backend = "cmdstanr", chains = 2L, iter = 1000L, warmup = 500L,
      thin = 1L, fit_cores = 2L, fit_threads = NA_integer_, seed = 1
    ),
    list(...)
  )
}

#' A fitter that sleeps a known time and returns a marked object
#'
#' The orchestration is tested through this and through a `.metrics`
#' stand-in, so no model is compiled (COMMON § 6).
#'
#' @noRd
sleeping_fitter <- function(seconds = 0, leapfrog = 1000, label = "mock",
                            error_on = integer()) {
  calls <- 0L
  function(data, ...) {
    calls <<- calls + 1L
    if (calls %in% error_on) {
      stop("mock failure on call ", calls, call. = FALSE)
    }
    if (seconds > 0) Sys.sleep(seconds)
    structure(
      list(
        label = label, call_index = calls, rows = nrow(data),
        leapfrog = leapfrog
      ),
      class = "mock_bench_fit"
    )
  }
}

#' Metrics for a `mock_bench_fit`, in the row's own shape
#' @noRd
mock_metrics <- function(scale = 1) {
  function(fit, wall_seconds = NA_real_, compile_seconds = NA_real_,
           compile_source = NA_character_, ...) {
    benchmark_from_parts(
      elapsed = fake_elapsed(warmup = c(1, 1) * scale, sample = c(1, 1) * scale),
      leapfrog = fake_leapfrog(rep(fit$leapfrog, 2L)),
      diagnostics = fake_bench_diagnostics(),
      nuts = NULL,
      wall_seconds = wall_seconds,
      compile_seconds = compile_seconds,
      compile_source = compile_source,
      settings = fake_bench_settings()
    )
  }
}
