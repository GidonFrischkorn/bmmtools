# The benchmark module.
#
# A speed claim needs units that say what they measure. Two facts, both
# measured 2026-10-05 (spec 15 § 15.0a), shape everything here:
# `rstan::get_elapsed_time()` gives per-chain warmup and sampling
# seconds under both backends, and nothing in a fit records compile
# time. Per-chain times are chain-local work --- they barely move
# between `cores = 4` and `cores = 1`, while wall time nearly doubles ---
# so every rate divides by seconds summed over chains, never by wall
# time, which would be wrong by about the number of cores (D69, D70).
# The rules live in benchmark_from_parts() so that they are tested on
# hand-built tables, as convergence_from_parts() is.

#' The benchmark row contract
#'
#' Column names and their `typeof()`. Identity and settings first, then
#' time, work, efficiency and the convergence numbers a speed figure
#' must not be read without.
#'
#' @noRd
benchmark_contract <- function() {
  c(
    implementation = "character",
    pair = "integer",
    run = "integer",
    warmup_run = "logical",
    data_id = "character",
    data_rows = "integer",
    backend = "character",
    chains = "integer",
    iter = "integer",
    warmup = "integer",
    thin = "integer",
    fit_cores = "integer",
    fit_threads = "integer",
    seed = "double",
    wall_seconds = "double",
    compile_seconds = "double",
    compile_source = "character",
    warmup_chain_seconds = "double",
    sampling_chain_seconds = "double",
    total_chain_seconds = "double",
    max_chain_seconds = "double",
    overhead_seconds = "double",
    leapfrog_sampling = "double",
    leapfrog_warmup = "double",
    us_per_gradient = "double",
    us_per_gradient_min = "double",
    us_per_gradient_max = "double",
    us_per_gradient_warmup = "double",
    ess_bulk_min = "double",
    ess_tail_min = "double",
    ess_bulk_per_chain_sec_sampling = "double",
    ess_tail_per_chain_sec_sampling = "double",
    ess_bulk_per_chain_sec_total = "double",
    ess_tail_per_chain_sec_total = "double",
    max_rhat = "double",
    n_divergent = "integer",
    n_max_treedepth = "integer",
    n_variables = "integer",
    converged = "logical",
    error = "character"
  )
}

#' @noRd
benchmark_contract_columns <- function() {
  names(benchmark_contract())
}

#' Metrics measured once per implementation, not once per run
#'
#' The compile is timed once and its value copied into every run's row,
#' so treating it as one value per pair would present a single
#' measurement as `runs` independent ones, with an interval of width
#' zero. It is reported as the single measurement it is.
#'
#' @noRd
metrics_measured_once <- function() {
  "compile_seconds"
}

#' The metrics `summary()` reports unless told otherwise
#' @noRd
benchmark_default_metrics <- function() {
  c(
    "us_per_gradient",
    "ess_bulk_per_chain_sec_sampling",
    "ess_tail_per_chain_sec_sampling",
    "sampling_chain_seconds",
    "wall_seconds",
    "compile_seconds"
  )
}

#' The settings a run records, with every slot present
#' @noRd
benchmark_settings <- function(settings = list()) {
  defaults <- list(
    backend = NA_character_, chains = NA_integer_, iter = NA_integer_,
    warmup = NA_integer_, thin = NA_integer_, fit_cores = NA_integer_,
    fit_threads = NA_integer_, seed = NA_real_
  )
  given <- settings[!vapply(settings, is.null, logical(1))]
  out <- utils::modifyList(defaults, given)
  out$backend <- as.character(out$backend)[[1L]]
  counts <- c(
    "chains", "iter", "warmup", "thin", "fit_cores", "fit_threads"
  )
  for (nm in counts) {
    out[[nm]] <- as.integer(out[[nm]])[[1L]]
  }
  out$seed <- as.double(out$seed)[[1L]]
  out
}

#' Drop the sampler's own variables from a diagnostics table
#'
#' D75: the ESS minimum describes the **estimated parameters**. `lp__`
#' and `lprior` have a finite rhat and ESS, so the rhat rule that drops
#' constants keeps them, and on the committed fixture `lp__` *was* the
#' minimum (204.71 against 220.34 for the smallest parameter, measured
#' 2026-10-05). `lprior` is worse than cosmetic here: it is defined by
#' each implementation's own priors and parametrisation, so a rate built
#' on it would compare priors rather than speed.
#'
#' @noRd
drop_sampler_variables <- function(diagnostics) {
  if (is.null(diagnostics) || !"variable" %in% names(diagnostics)) {
    return(diagnostics)
  }
  diagnostics[!diagnostics$variable %in% c("lp__", "lprior"), , drop = FALSE]
}

#' A ratio that is `NA` rather than `Inf` or `NaN`
#'
#' A rate whose denominator is zero or unknown is unknown, and a
#' benchmark that printed `Inf` would be read as a measurement.
#'
#' @noRd
safe_ratio <- function(numerator, denominator) {
  ok <- length(numerator) == 1L && length(denominator) == 1L &&
    !is.na(numerator) && !is.na(denominator) && denominator > 0
  if (!ok) {
    return(NA_real_)
  }
  as.double(numerator / denominator)
}

#' Summaries that give `NA` for nothing rather than `Inf` or a warning
#' @noRd
min_finite <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else min(x)
}

#' @noRd
max_finite <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else max(x)
}

#' @noRd
median_finite <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else as.double(stats::median(x))
}

#' A sum over chains that is `NA` when any chain is unmeasured
#'
#' Summing the chains that did report would understate the total and
#' overstate every rate built on it, with nothing to show that a chain
#' was missing.
#'
#' @noRd
sum_or_na <- function(x) {
  if (length(x) == 0L || anyNA(x)) NA_real_ else sum(x)
}

#' The per-chain table, from an elapsed matrix and leapfrog counts
#'
#' One row per chain. `leapfrog_sampling` is `NA` when the fit carried no
#' sampler parameters, and `leapfrog_warmup` is `NA` unless warmup draws
#' were saved (measured: `nuts_params()` returns the saved draws only, so
#' `save_warmup = TRUE` is what makes the warmup phase visible).
#'
#' @noRd
benchmark_chain_table <- function(elapsed, leapfrog, leapfrog_warmup) {
  elapsed <- as.matrix(elapsed)
  n <- nrow(elapsed)
  pick <- function(table, chain) {
    if (is.null(table) || !"chain" %in% names(table)) {
      return(NA_real_)
    }
    hit <- table$leapfrog[match(chain, table$chain)]
    if (length(hit) != 1L) NA_real_ else as.double(hit)
  }
  chain <- seq_len(n)
  sampling <- vapply(chain, function(i) pick(leapfrog, i), double(1))
  warm <- vapply(chain, function(i) pick(leapfrog_warmup, i), double(1))
  sample_seconds <- as.double(elapsed[, "sample"])
  warmup_seconds <- as.double(elapsed[, "warmup"])
  tibble::tibble(
    chain = as.integer(chain),
    warmup_seconds = warmup_seconds,
    sample_seconds = sample_seconds,
    leapfrog_sampling = sampling,
    leapfrog_warmup = warm,
    us_per_gradient = vapply(
      chain, function(i) safe_ratio(1e6 * sample_seconds[[i]], sampling[[i]]),
      double(1)
    ),
    us_per_gradient_warmup = vapply(
      chain, function(i) safe_ratio(1e6 * warmup_seconds[[i]], warm[[i]]),
      double(1)
    )
  )
}

#' Refuse a gradient cost when the draws were thinned
#'
#' `brms::nuts_params()` reports the **saved** draws, so with `thin = k`
#' the leapfrog steps cover one iteration in k while the sampler's
#' `sample` time covers all of them. Measured 2026-10-05 on one
#' `mixture2p` fit, 2 chains, iter 1000: 256 µs per gradient at
#' `thin = 1`, 495 µs at `thin = 2` and 1271 µs at `thin = 5`, while the
#' sampling time barely moved (0.53-0.58 s) — thinning discards output,
#' not work. The count is what it says it is and stays; the cost, whose
#' numerator and denominator then cover different iterations, is `NA`
#' (D78).
#'
#' @noRd
blank_cost_if_thinned <- function(chain_table, thin) {
  if (is.na(thin) || thin <= 1L) {
    return(chain_table)
  }
  chain_table$us_per_gradient <- NA_real_
  chain_table$us_per_gradient_warmup <- NA_real_
  chain_table
}

#' One benchmark row from parts already extracted
#'
#' @param elapsed A chains x `c("warmup", "sample")` matrix of seconds,
#'   as `rstan::get_elapsed_time()` returns.
#' @param leapfrog A table with `chain` and `leapfrog`, the post-warmup
#'   leapfrog steps per chain, or `NULL` when the fit has none.
#' @param leapfrog_warmup The same for the warmup phase, or `NULL`.
#' @param diagnostics,nuts Passed to [convergence_from_parts()].
#' @param wall_seconds,compile_seconds,compile_source What the caller
#'   measured around the fit; no fit records either.
#' @param identity Identity columns (`implementation`, `pair`, ...).
#' @param settings Sampler settings as run.
#' @param error The message of a run that failed, `""` otherwise.
#' @noRd
benchmark_from_parts <- function(elapsed,
                                 leapfrog,
                                 diagnostics,
                                 nuts = NULL,
                                 leapfrog_warmup = NULL,
                                 wall_seconds = NA_real_,
                                 compile_seconds = NA_real_,
                                 compile_source = NA_character_,
                                 treedepth_max = 10,
                                 thresholds = NULL,
                                 identity = list(),
                                 settings = list(),
                                 error = "") {
  thresholds <- thresholds %||% check_thresholds()
  settings <- benchmark_settings(settings)
  chain_table <- benchmark_chain_table(elapsed, leapfrog, leapfrog_warmup)
  chain_table <- blank_cost_if_thinned(chain_table, settings$thin)
  convergence <- convergence_from_parts(
    drop_sampler_variables(diagnostics), nuts, treedepth_max, thresholds
  )

  warmup_chain_seconds <- sum_or_na(chain_table$warmup_seconds)
  sampling_chain_seconds <- sum_or_na(chain_table$sample_seconds)
  total_chain_seconds <- sum_or_na(
    chain_table$warmup_seconds + chain_table$sample_seconds
  )
  chain_totals <- chain_table$warmup_seconds + chain_table$sample_seconds
  max_chain_seconds <- max_finite(chain_totals)
  wall_seconds <- as.double(wall_seconds)[[1L]]
  overhead <- if (is.na(wall_seconds) || is.na(max_chain_seconds)) {
    NA_real_
  } else {
    wall_seconds - max_chain_seconds
  }

  leapfrog_sampling <- sum_or_na(chain_table$leapfrog_sampling)
  leapfrog_warmup_total <- sum_or_na(chain_table$leapfrog_warmup)
  gradient_cost <- chain_table$us_per_gradient
  gradient_cost_warmup <- chain_table$us_per_gradient_warmup

  row <- tibble::tibble(
    implementation = NA_character_, pair = NA_integer_, run = NA_integer_,
    warmup_run = NA, data_id = NA_character_, data_rows = NA_integer_,
    backend = settings$backend, chains = settings$chains,
    iter = settings$iter, warmup = settings$warmup, thin = settings$thin,
    fit_cores = settings$fit_cores, fit_threads = settings$fit_threads,
    seed = settings$seed,
    wall_seconds = wall_seconds,
    compile_seconds = as.double(compile_seconds)[[1L]],
    compile_source = as.character(compile_source)[[1L]],
    warmup_chain_seconds = warmup_chain_seconds,
    sampling_chain_seconds = sampling_chain_seconds,
    total_chain_seconds = total_chain_seconds,
    max_chain_seconds = max_chain_seconds,
    overhead_seconds = overhead,
    leapfrog_sampling = leapfrog_sampling,
    leapfrog_warmup = leapfrog_warmup_total,
    us_per_gradient = median_finite(gradient_cost),
    us_per_gradient_min = min_finite(gradient_cost),
    us_per_gradient_max = max_finite(gradient_cost),
    us_per_gradient_warmup = median_finite(gradient_cost_warmup),
    ess_bulk_min = as.double(convergence$min_ess_bulk),
    ess_tail_min = as.double(convergence$min_ess_tail),
    ess_bulk_per_chain_sec_sampling = safe_ratio(
      convergence$min_ess_bulk, sampling_chain_seconds
    ),
    ess_tail_per_chain_sec_sampling = safe_ratio(
      convergence$min_ess_tail, sampling_chain_seconds
    ),
    ess_bulk_per_chain_sec_total = safe_ratio(
      convergence$min_ess_bulk, total_chain_seconds
    ),
    ess_tail_per_chain_sec_total = safe_ratio(
      convergence$min_ess_tail, total_chain_seconds
    ),
    max_rhat = as.double(convergence$max_rhat),
    n_divergent = as.integer(convergence$n_divergent),
    n_max_treedepth = as.integer(convergence$n_max_treedepth),
    n_variables = as.integer(convergence$n_variables),
    converged = convergence$pass,
    error = as.character(error)[[1L]]
  )
  row <- fill_benchmark_identity(row, identity)
  new_bmmtools_benchmark(row, chains = chain_table)
}

#' Write the identity columns a caller knows into a row
#' @noRd
fill_benchmark_identity <- function(row, identity) {
  identity <- identity[!vapply(identity, is.null, logical(1))]
  coercers <- list(
    implementation = as.character, pair = as.integer, run = as.integer,
    warmup_run = as.logical, data_id = as.character, data_rows = as.integer
  )
  for (nm in intersect(names(identity), names(coercers))) {
    row[[nm]] <- coercers[[nm]](identity[[nm]])[[1L]]
  }
  row
}

#' A row for a run that never produced a fit
#'
#' Its measurements are `NA`, not zero: the run was not measured.
#'
#' @noRd
empty_benchmark_row <- function(identity = list(), settings = list(),
                                error = "", compile_seconds = NA_real_,
                                compile_source = NA_character_) {
  # as many NA chains as the call asked for, so that a failed 4-chain
  # run does not read as a 1-chain one in the `chains` attribute
  chains <- benchmark_settings(settings)$chains
  if (is.na(chains) || chains < 1L) {
    chains <- 1L
  }
  benchmark_from_parts(
    elapsed = matrix(
      NA_real_,
      nrow = chains, ncol = 2L,
      dimnames = list(paste0("chain:", seq_len(chains)), c("warmup", "sample"))
    ),
    leapfrog = NULL,
    diagnostics = tibble::tibble(
      variable = character(), rhat = double(), ess_bulk = double(),
      ess_tail = double()
    ),
    nuts = NULL, identity = identity, settings = settings, error = error,
    compile_seconds = compile_seconds, compile_source = compile_source
  )
}

#' Attach the benchmark class and its attributes
#' @noRd
new_bmmtools_benchmark <- function(x, chains = NULL, settings = NULL,
                                   compile = NULL) {
  out <- tibble::as_tibble(x)
  class(out) <- c("bmmtools_benchmark", class(out))
  attr(out, "chains") <- chains
  attr(out, "settings") <- settings
  attr(out, "compile") <- compile
  out
}

# reading the parts out of a fit ------------------------------------------

#' Per-chain warmup and sampling seconds, or an `NA` matrix
#'
#' `rstan::get_elapsed_time(fit$fit)` survives brms's conversion of a
#' cmdstanr run, because `fit$fit` is a `stanfit` under both backends
#' (measured 2026-10-05).
#'
#' @noRd
fit_elapsed_time <- function(fit) {
  out <- tryCatch(rstan::get_elapsed_time(fit$fit), error = function(e) NULL)
  ok <- is.matrix(out) && all(c("warmup", "sample") %in% colnames(out)) &&
    nrow(out) >= 1L
  if (!ok) {
    chains <- tryCatch(as.integer(fit$fit@sim$chains), error = function(e) 1L)
    if (length(chains) != 1L || is.na(chains) || chains < 1L) chains <- 1L
    return(matrix(
      NA_real_,
      nrow = chains, ncol = 2L,
      dimnames = list(paste0("chain:", seq_len(chains)), c("warmup", "sample"))
    ))
  }
  out[, c("warmup", "sample"), drop = FALSE]
}

#' Leapfrog steps per chain from a `nuts_params()` table
#' @noRd
leapfrog_from_nuts <- function(nuts) {
  wanted <- c("Chain", "Parameter", "Value")
  if (is.null(nuts) || !all(wanted %in% names(nuts))) {
    return(NULL)
  }
  lf <- nuts[nuts$Parameter == "n_leapfrog__", , drop = FALSE]
  if (nrow(lf) == 0L) {
    return(NULL)
  }
  chains <- sort(unique(as.integer(lf$Chain)))
  tibble::tibble(
    chain = chains,
    leapfrog = vapply(
      chains, function(i) sum(lf$Value[as.integer(lf$Chain) == i]), double(1)
    ),
    draws = vapply(
      chains, function(i) sum(as.integer(lf$Chain) == i), integer(1)
    )
  )
}

#' Whether the fit kept its warmup draws
#' @noRd
fit_saved_warmup <- function(fit) {
  value <- tryCatch(fit$fit@sim$warmup2, error = function(e) NULL)
  if (!is.numeric(unlist(value))) {
    return(0L)
  }
  as.integer(max(unlist(value), na.rm = TRUE))
}

#' The sampling and warmup leapfrog tables
#'
#' `brms::nuts_params()` reports the **saved** draws, so the warmup phase
#' is reachable only when the fit was run with `save_warmup = TRUE`; the
#' warmup counts are then the total minus the post-warmup ones.
#'
#' @noRd
fit_leapfrog_tables <- function(fit, inc_warmup) {
  nuts_or_null <- function(warmup) {
    tryCatch(
      if (warmup) {
        brms::nuts_params(fit, inc_warmup = TRUE)
      } else {
        brms::nuts_params(fit)
      },
      error = function(e) NULL
    )
  }
  sampling <- leapfrog_from_nuts(nuts_or_null(FALSE))
  warmup <- NULL
  if (isTRUE(inc_warmup)) {
    warmup <- leapfrog_warmup_from_totals(
      leapfrog_from_nuts(nuts_or_null(TRUE)), sampling
    )
  }
  list(sampling = sampling, warmup = warmup)
}

#' The warmup phase as the difference between all saved draws and the
#' post-warmup ones
#'
#' Kept separate from the fit so that the subtraction is tested on
#' hand-built tables; `nuts_params()` itself cannot be injected.
#'
#' @noRd
leapfrog_warmup_from_totals <- function(total, sampling) {
  if (is.null(total) || is.null(sampling)) {
    return(NULL)
  }
  idx <- match(total$chain, sampling$chain)
  extra <- total$draws - sampling$draws[idx]
  if (all(is.na(extra) | extra <= 0L)) {
    return(NULL)
  }
  leapfrog <- total$leapfrog - sampling$leapfrog[idx]
  # a chain whose warmup draws were not saved has an unknown warmup
  # count, never a count of zero: brms returns the post-warmup draws and
  # warns when asked for warmup it does not have
  leapfrog[is.na(extra) | extra <= 0L] <- NA_real_
  tibble::tibble(
    chain = total$chain,
    leapfrog = leapfrog,
    draws = as.integer(pmax(extra, 0L))
  )
}

#' The sampler settings a fit records
#' @noRd
fit_settings <- function(fit) {
  sim <- tryCatch(fit$fit@sim, error = function(e) NULL)
  args <- tryCatch(fit$fit@stan_args[[1L]], error = function(e) NULL)
  one_int <- function(x) {
    if (is.numeric(x) && length(x) >= 1L) as.integer(x[[1L]]) else NA_integer_
  }
  backend <- tryCatch(fit$backend, error = function(e) NULL)
  list(
    backend = if (is.character(backend) && length(backend) == 1L) {
      backend
    } else {
      NA_character_
    },
    chains = one_int(sim$chains), iter = one_int(sim$iter),
    warmup = one_int(sim$warmup), thin = one_int(sim$thin),
    fit_threads = one_int(args$num_threads),
    seed = if (is.numeric(args$seed)) as.double(args$seed[[1L]]) else NA_real_
  )
}

# benchmark_metrics() ----------------------------------------------------

#' Measure one fit: compile, per-chain time, ESS per second, gradient cost
#'
#' What a fit cost, in units that name their phase and their denominator.
#' Sampling time comes from the sampler's own clock per chain
#' ([rstan::get_elapsed_time()], which works under the cmdstanr backend
#' too), and the work from `n_leapfrog__`
#' ([brms::nuts_params()]). Every rate divides by seconds **summed over
#' chains**, never by wall time: per-chain times are chain-local work, so
#' a wall-time denominator would make the same fit look cheaper on more
#' cores.
#'
#' @param fit A fitted model, `bmmfit` or `brmsfit`.
#' @param wall_seconds The wall time of the fitting call, which the
#'   caller must measure: no fit records it.
#' @param compile_seconds Compile time from [benchmark_compile()]; no fit
#'   records that either.
#' @param compile_source How `compile_seconds` was obtained, for the
#'   reader of the row: `"forced"`, `"fit"` or `"none"`.
#' @param variables Restrict the rhat and ESS set to these variables.
#' @param inc_warmup Read the warmup phase's leapfrog steps as well.
#'   `NULL`, the default, uses them when the fit saved its warmup draws,
#'   which happens only with `save_warmup = TRUE`. A fit that did not
#'   save them reports `NA` for the warmup work whatever this is set to:
#'   asking brms for warmup draws it does not have returns the
#'   post-warmup ones, which would read as zero warmup work.
#' @param rhat_max,ess_bulk_min,ess_tail_min,divergent_max Thresholds for
#'   the `converged` column, as in [check_convergence()].
#' @param ... Not used.
#'
#' @return A one-row tibble of class `bmmtools_benchmark`; see
#'   [benchmark()] for the columns. The per-chain table is in the
#'   attribute `chains`.
#'
#' @details
#' The minimum ESS is taken over the **estimated** parameters: a
#' parameter fixed by a constant prior has no posterior variance, no rhat
#' and no ESS, and `n_variables` records how many were assessed. This is
#' [check_convergence()]'s rule, and the convergence columns come from
#' the same code.
#'
#' @examples
#' \dontrun{
#' started <- Sys.time()
#' fit <- bmm::bmm(formula, data, model, backend = "cmdstanr")
#' row <- benchmark_metrics(
#'   fit,
#'   wall_seconds = as.double(difftime(Sys.time(), started, units = "secs"))
#' )
#' row$us_per_gradient
#' }
#'
#' @seealso [benchmark()], [benchmark_compile()], [check_convergence()]
#' @export
benchmark_metrics <- function(fit, ...) {
  UseMethod("benchmark_metrics")
}

#' @rdname benchmark_metrics
#' @export
benchmark_metrics.default <- function(fit, ...) {
  cli::cli_abort(
    "{.arg fit} must be a {.cls brmsfit}, not {.obj_type_friendly {fit}}."
  )
}

#' @rdname benchmark_metrics
#' @export
benchmark_metrics.brmsfit <- function(fit,
                                      wall_seconds = NA_real_,
                                      compile_seconds = NA_real_,
                                      compile_source = NA_character_,
                                      variables = NULL,
                                      inc_warmup = NULL,
                                      rhat_max = 1.05,
                                      ess_bulk_min = 400,
                                      ess_tail_min = NULL,
                                      divergent_max = 10,
                                      ...) {
  rlang::check_installed("brms", "to measure a fit.")
  rlang::check_installed("rstan", "to read a fit's per-chain times.")
  rlang::check_dots_empty()
  thresholds <- check_thresholds(
    rhat_max, ess_bulk_min, ess_tail_min, divergent_max
  )
  draws <- select_draw_variables(posterior::as_draws_array(fit), variables)
  diagnostics <- posterior::summarise_draws(
    draws, posterior::default_convergence_measures()
  )
  nuts <- tryCatch(brms::nuts_params(fit), error = function(e) NULL)
  # asking for warmup a fit did not save makes brms return the
  # post-warmup draws with a warning, which would read as zero warmup
  # work; the fit decides, whatever `inc_warmup` says
  inc_warmup <- !isFALSE(inc_warmup) && fit_saved_warmup(fit) > 0L
  leapfrog <- fit_leapfrog_tables(fit, inc_warmup)
  benchmark_from_parts(
    elapsed = fit_elapsed_time(fit),
    leapfrog = leapfrog$sampling,
    leapfrog_warmup = leapfrog$warmup,
    diagnostics = diagnostics,
    nuts = nuts,
    wall_seconds = wall_seconds,
    compile_seconds = compile_seconds,
    compile_source = compile_source,
    treedepth_max = fit_treedepth_max(fit),
    thresholds = thresholds,
    settings = fit_settings(fit)
  )
}

# benchmark_compile() ----------------------------------------------------

#' Compile Stan code with the backend, forcing a fresh compile
#' @noRd
compile_stan_model <- function(code, backend, force, ...) {
  if (identical(backend, "cmdstanr")) {
    rlang::check_installed("cmdstanr", "to compile a model.")
    file <- cmdstanr::write_stan_file(code)
    model <- cmdstanr::cmdstan_model(file, force_recompile = force, ...)
    exe <- tryCatch(model$exe_file(), error = function(e) NA_character_)
    return(list(exe = exe))
  }
  rlang::check_installed("rstan", "to compile a model.")
  rstan::stan_model(model_code = code, ...)
  list(exe = NA_character_)
}

#' Time the compile step on its own
#'
#' No fit records how long it took to compile, so a benchmark times a
#' compile of its own. `force = TRUE` is what makes the number a compile
#' rather than a cache hit: under cmdstanr a second call to the same code
#' reuses the executable in 0.02 s against 4.7 s for a compile (measured
#' 2026-10-05), while rstan recompiles every time (23.3 s).
#'
#' @param code Stan code, as [bmm::stancode()] or
#'   `brms::make_stancode()` returns it.
#' @param backend `"cmdstanr"` or `"rstan"`.
#' @param force Recompile even when an executable is cached. Leaving it
#'   `TRUE` is the point; `FALSE` measures a reuse.
#' @param ... Passed to the backend's compiler.
#' @param .compiler The compiling function, the backend's by default.
#'   Tests inject a stand-in so that nothing is compiled.
#'
#' @return A one-row tibble with `compile_seconds`, `backend`, `forced`,
#'   `code_hash` and `exe`.
#'
#' @examples
#' \dontrun{
#' code <- bmm::stancode(formula, data, model, backend = "cmdstanr")
#' benchmark_compile(code, backend = "cmdstanr")
#' }
#'
#' @seealso [benchmark()], [benchmark_metrics()]
#' @export
benchmark_compile <- function(code, backend = "cmdstanr", force = TRUE, ...,
                              .compiler = NULL) {
  if (!is.character(code) || length(code) != 1L || is.na(code)) {
    cli::cli_abort(
      "{.arg code} must be a single string of Stan code, \\
       not {.obj_type_friendly {code}}."
    )
  }
  backend <- as.character(backend)
  if (length(backend) != 1L || !backend %in% c("cmdstanr", "rstan")) {
    cli::cli_abort(
      "{.arg backend} must be {.val cmdstanr} or {.val rstan}, \\
       not {.val {backend}}."
    )
  }
  if (!is.logical(force) || length(force) != 1L || is.na(force)) {
    cli::cli_abort("{.arg force} must be {.code TRUE} or {.code FALSE}.")
  }
  compiler <- .compiler %||% compile_stan_model
  if (!is.function(compiler)) {
    cli::cli_abort(
      "{.arg .compiler} must be a function, \\
       not {.obj_type_friendly {compiler}}."
    )
  }
  started <- Sys.time()
  out <- compiler(code, backend = backend, force = force, ...)
  seconds <- as.double(difftime(Sys.time(), started, units = "secs"))
  exe <- out$exe
  tibble::tibble(
    compile_seconds = seconds,
    backend = backend,
    forced = force,
    code_hash = rlang::hash(code),
    exe = if (is.character(exe) && length(exe) == 1L) exe else NA_character_
  )
}

# benchmark() ------------------------------------------------------------

#' Check and normalise the implementations
#' @noRd
normalise_implementations <- function(implementations,
                                      call = rlang::caller_env()) {
  if (!is.list(implementations) || length(implementations) == 0L) {
    cli::cli_abort(
      "{.arg implementations} must be a non-empty named list of \\
       functions, not {.obj_type_friendly {implementations}}.",
      call = call
    )
  }
  names <- names(implementations)
  if (is.null(names) || any(!nzchar(names)) || anyNA(names)) {
    cli::cli_abort(
      "{.arg implementations} must have names: one per implementation.",
      call = call
    )
  }
  duplicated_names <- unique(names[duplicated(names)])
  if (length(duplicated_names) > 0L) {
    cli::cli_abort(
      "{.arg implementations} must have unique names; \\
       {.val {duplicated_names}} {?is/are} repeated.",
      call = call
    )
  }
  lapply(stats::setNames(seq_along(implementations), names), function(i) {
    spec <- implementations[[i]]
    if (is.function(spec)) {
      spec <- list(fit = spec)
    }
    if (!is.list(spec) || !is.function(spec$fit)) {
      cli::cli_abort(
        "{.arg implementations[[{names[[i]]}]]} must be a function or a \\
         list with a {.field fit} function.",
        call = call
      )
    }
    list(fit = spec$fit, compile = spec$compile)
  })
}

#' The data sets, one per pair, checked for length
#' @noRd
benchmark_data_list <- function(data, needed, call = rlang::caller_env()) {
  if (is.data.frame(data)) {
    sets <- rep(list(data), needed)
    names(sets) <- rep("data", needed)
    return(sets)
  }
  all_frames <- is.list(data) && length(data) > 0L &&
    all(vapply(data, is.data.frame, logical(1)))
  if (!all_frames) {
    cli::cli_abort(
      "{.arg data} must be a data frame or a list of data frames, \\
       not {.obj_type_friendly {data}}.",
      call = call
    )
  }
  if (length(data) < needed) {
    cli::cli_abort(
      "{.arg data} has {length(data)} data set{?s} but \\
       {needed} {?is/are} needed \\
       ({.arg runs} + {.arg warmup_runs}); it is not recycled, because \\
       a reused data set would make the paired ratios dependent.",
      call = call
    )
  }
  sets <- data[seq_len(needed)]
  if (is.null(names(sets)) || any(!nzchar(names(sets)))) {
    names(sets) <- paste0("set", seq_len(needed))
  }
  sets
}

#' Fill settings columns the metrics did not know
#' @noRd
fill_na_columns <- function(row, values) {
  for (nm in names(values)) {
    value <- values[[nm]]
    if (is.null(value) || !nm %in% names(row)) {
      next
    }
    if (length(row[[nm]]) == 1L && is.na(row[[nm]])) {
      row[[nm]] <- coerce_like(value, row[[nm]])
    }
  }
  row
}

#' Coerce `value` to the type of `template`
#' @noRd
coerce_like <- function(value, template) {
  if (is.integer(template)) {
    return(as.integer(value)[[1L]])
  }
  if (is.double(template)) {
    return(as.double(value)[[1L]])
  }
  if (is.logical(template)) {
    return(as.logical(value)[[1L]])
  }
  as.character(value)[[1L]]
}

#' Compare implementations on the same data, in interleaved runs
#'
#' One row per run, with compile time and sampling time kept apart, ESS
#' per chain-second for bulk and tail, and the cost of one gradient
#' evaluation. Runs are interleaved --- A, B, A, B --- because running
#' every A before every B would confound the comparison with drift in
#' machine load, and the first run of each implementation is discarded:
#' it carries the fresh compile and a cold machine. The discarded run is
#' still returned, marked `warmup_run`, so that nothing is dropped
#' silently.
#'
#' No fit is written. The module never calls [fit_cached()], so a
#' benchmark fit can never be served to a recovery grid, and no cache key
#' changes because of it.
#'
#' @param implementations A named list, one element per implementation.
#'   Each is either a function `(data, ...) -> fit` or a list with a
#'   `fit` function and, optionally, a `compile` function
#'   `(data, ...) -> Stan code` whose compile is timed once.
#' @param data A data frame used for every pair, or a list of data
#'   frames, one per pair. A list is not recycled: it must hold
#'   `runs + warmup_runs` sets, because a reused data set would make the
#'   paired ratios dependent.
#' @param ... Passed to every implementation's fitter (`chains`, `iter`,
#'   `cores`, `seed`, `backend`, ...). Pass the same arguments to all of
#'   them; `summary()` warns when the recorded settings differ.
#' @param runs Recorded pairs per implementation.
#' @param warmup_runs Leading runs per implementation that are recorded
#'   and then excluded from `summary()`.
#' @param compile Measure a forced compile once per implementation, for
#'   those that supply a `compile` hook.
#' @param on_error `"record"` keeps a row whose measurements are `NA`,
#'   puts the message in `error` and goes on; `"stop"` aborts.
#' @param .metrics The measuring function, [benchmark_metrics()] by
#'   default. Tests inject a stand-in so that no model is fitted.
#' @param .compiler Passed to [benchmark_compile()].
#'
#' @return A tibble of class `bmmtools_benchmark`, one row per run:
#'
#'   * identity --- `implementation`, `pair` (0 for a warmup pair),
#'     `run`, `warmup_run`, `data_id`, `data_rows`;
#'   * settings as run --- `backend`, `chains`, `iter`, `warmup`,
#'     `thin`, `fit_cores`, `fit_threads`, `seed`;
#'   * time in seconds --- `wall_seconds` (the fitting call),
#'     `compile_seconds` with `compile_source`,
#'     `warmup_chain_seconds`, `sampling_chain_seconds`,
#'     `total_chain_seconds` (each summed over chains),
#'     `max_chain_seconds` and `overhead_seconds`
#'     (`wall_seconds - max_chain_seconds`, which contains a per-fit
#'     compile under a backend that recompiles, and R-side overhead
#'     under cmdstanr with a cached executable);
#'   * work --- `leapfrog_sampling`, `leapfrog_warmup` (`NA` unless
#'     warmup draws were saved), `us_per_gradient` (the median over
#'     chains of chain sampling time over chain leapfrog steps), with
#'     `us_per_gradient_min`, `us_per_gradient_max` and
#'     `us_per_gradient_warmup`;
#'   * efficiency --- `ess_bulk_min`, `ess_tail_min` and the four rates
#'     `ess_{bulk,tail}_per_chain_sec_{sampling,total}`, per
#'     chain-second, over post-warmup or over warmup plus sampling;
#'   * what the speed must be read with --- `max_rhat`, `n_divergent`,
#'     `n_max_treedepth`, `n_variables`, `converged`, `error`;
#'   * every [run_info()] column, per run.
#'
#'   The per-chain table of every run is in the attribute `chains`, the
#'   compile table in `compile`, and the call and its arguments in
#'   `settings`.
#'
#' @examples
#' \dontrun{
#' bench <- benchmark(
#'   implementations = list(
#'     current = function(data, ...) brms::brm(current_formula, data, ...),
#'     rebuilt = function(data, ...) {
#'       bmm::bmm(
#'         bmm::bmf(kappa ~ 1, thetat ~ 1), data,
#'         bmm::mixture2p(resp_error = "y"), ...
#'       )
#'     }
#'   ),
#'   data = data_sets, runs = 5,
#'   backend = "cmdstanr", chains = 4, iter = 2000, cores = 4
#' )
#' summary(bench)$ratios
#' }
#'
#' @seealso [benchmark_metrics()], [benchmark_compile()], [run_info()]
#' @export
benchmark <- function(implementations,
                      data,
                      ...,
                      runs = 5L,
                      warmup_runs = 1L,
                      compile = TRUE,
                      on_error = c("record", "stop"),
                      .metrics = NULL,
                      .compiler = NULL) {
  on_error <- rlang::arg_match(on_error)
  impls <- normalise_implementations(implementations)
  runs <- check_count(runs, "runs", 1L)
  warmup_runs <- check_count(warmup_runs, "warmup_runs", 0L)
  if (!is.logical(compile) || length(compile) != 1L || is.na(compile)) {
    cli::cli_abort("{.arg compile} must be {.code TRUE} or {.code FALSE}.")
  }
  metrics <- .metrics %||% benchmark_metrics
  if (!is.function(metrics)) {
    cli::cli_abort(
      "{.arg .metrics} must be a function, \\
       not {.obj_type_friendly {metrics}}."
    )
  }
  dots <- rlang::list2(...)
  sets <- benchmark_data_list(data, warmup_runs + runs)
  compile_table <- benchmark_compile_step(
    impls, sets[[1L]], dots, compile, .compiler
  )

  pairs <- c(rep(0L, warmup_runs), seq_len(runs))
  rows <- list()
  chains <- list()
  failures <- character()
  for (index in seq_along(pairs)) {
    set <- sets[[index]]
    for (name in names(impls)) {
      identity <- list(
        implementation = name, pair = pairs[[index]], run = index,
        warmup_run = pairs[[index]] == 0L,
        data_id = names(sets)[[index]], data_rows = nrow(set)
      )
      compile_row <- compile_table[compile_table$implementation == name, ]
      row <- benchmark_one_run(
        impls[[name]]$fit, set, dots, metrics, identity,
        compile_seconds = compile_row$compile_seconds,
        compile_source = compile_row$compile_source,
        on_error = on_error
      )
      if (nzchar(row$error)) {
        failures <- c(failures, paste0(name, " (run ", index, ")"))
      }
      chains[[length(chains) + 1L]] <- run_chain_table(row, identity)
      rows[[length(rows) + 1L]] <- row
    }
  }
  if (length(failures) > 0L) {
    cli::cli_warn(c(
      "{length(failures)} run{?s} failed and {?its/their} measurements \\
       {?is/are} {.code NA}.",
      i = "{.val {failures}}"
    ))
  }
  out <- dplyr::bind_rows(lapply(rows, strip_benchmark_class))
  new_bmmtools_benchmark(
    out,
    chains = dplyr::bind_rows(chains),
    settings = list(
      call = match.call(), runs = runs, warmup_runs = warmup_runs,
      implementations = names(impls), dots = dots
    ),
    compile = compile_table
  )
}

#' A plain tibble again, so that binding rows does not dispatch on the
#' benchmark class or carry one run's attributes into the whole table
#' @noRd
strip_benchmark_class <- function(row) {
  attr(row, "chains") <- NULL
  attr(row, "settings") <- NULL
  attr(row, "compile") <- NULL
  class(row) <- c("tbl_df", "tbl", "data.frame")
  row
}

#' The per-chain table of one run, labelled with the run's identity
#' @noRd
run_chain_table <- function(row, identity) {
  chains <- attr(row, "chains")
  if (is.null(chains) || nrow(chains) == 0L) {
    return(NULL)
  }
  dplyr::bind_cols(
    tibble::tibble(
      implementation = identity$implementation,
      pair = as.integer(identity$pair),
      run = as.integer(identity$run)
    )[rep(1L, nrow(chains)), ],
    chains
  )
}

#' Time one forced compile per implementation that offers a hook
#' @noRd
benchmark_compile_step <- function(impls, data, dots, compile, .compiler) {
  backend <- dots[["backend"]]
  backend <- if (is.character(backend) && length(backend) == 1L) {
    backend
  } else {
    "cmdstanr"
  }
  without <- character()
  rows <- lapply(names(impls), function(name) {
    hook <- impls[[name]]$compile
    if (!isTRUE(compile) || !is.function(hook)) {
      if (isTRUE(compile)) {
        without <<- c(without, name)
      }
      return(tibble::tibble(
        implementation = name, compile_seconds = NA_real_,
        compile_source = "none", code_hash = NA_character_
      ))
    }
    code <- rlang::exec(hook, data, !!!dots)
    timed <- benchmark_compile(
      code,
      backend = backend, force = TRUE, .compiler = .compiler
    )
    tibble::tibble(
      implementation = name, compile_seconds = timed$compile_seconds,
      compile_source = "forced", code_hash = timed$code_hash
    )
  })
  if (length(without) > 0L) {
    cli::cli_inform(c(
      "No compile hook for {.val {without}}: {.field compile_seconds} \\
       {?is/are} {.code NA}.",
      i = "Give the implementation a {.field compile} function to have its \\
           compile timed."
    ))
  }
  dplyr::bind_rows(rows)
}

#' Fit once, measure it, and label the row
#' @noRd
benchmark_one_run <- function(fitter, data, dots, metrics, identity,
                              compile_seconds, compile_source, on_error) {
  settings <- list(
    backend = dots[["backend"]], fit_cores = dots[["cores"]],
    seed = dots[["seed"]], chains = dots[["chains"]], iter = dots[["iter"]],
    warmup = dots[["warmup"]], thin = dots[["thin"]]
  )
  started <- Sys.time()
  fit <- tryCatch(
    rlang::exec(fitter, data, !!!dots),
    error = function(e) e
  )
  if (inherits(fit, "condition")) {
    if (identical(on_error, "stop")) {
      cli::cli_abort(
        "Implementation {.val {identity$implementation}} failed on run \\
         {identity$run}.",
        parent = fit
      )
    }
    return(dplyr::bind_cols(
      empty_benchmark_row(
        identity = identity, settings = settings,
        error = conditionMessage(fit), compile_seconds = compile_seconds,
        compile_source = compile_source
      ),
      run_info()
    ))
  }
  wall <- as.double(difftime(Sys.time(), started, units = "secs"))
  # measuring can fail where fitting did not (a fit without draws, a
  # backend whose times are unreadable); `on_error` covers that too,
  # because a benchmark that fits for hours must not be lost to it
  row <- tryCatch(
    metrics(
      fit,
      wall_seconds = wall, compile_seconds = compile_seconds,
      compile_source = compile_source
    ),
    error = function(e) e
  )
  if (inherits(row, "condition")) {
    if (identical(on_error, "stop")) {
      cli::cli_abort(
        "Measuring implementation {.val {identity$implementation}} failed \\
         on run {identity$run}.",
        parent = row
      )
    }
    return(dplyr::bind_cols(
      empty_benchmark_row(
        identity = identity, settings = settings,
        error = conditionMessage(row), compile_seconds = compile_seconds,
        compile_source = compile_source
      ),
      run_info()
    ))
  }
  row <- fill_benchmark_identity(row, identity)
  row$wall_seconds <- wall
  row$compile_seconds <- as.double(compile_seconds)[[1L]]
  row$compile_source <- as.character(compile_source)[[1L]]
  row <- fill_na_columns(row, settings)
  if (!nzchar(row$error %||% "")) {
    row$error <- ""
  }
  chains <- attr(row, "chains")
  row <- dplyr::bind_cols(row, run_info())
  attr(row, "chains") <- chains
  row
}

# summary ----------------------------------------------------------------

#' `TRUE` where `x` is `TRUE`, never where it is `NA`
#' @noRd
is_true_vec <- function(x) {
  !is.na(x) & x
}

#' Warn when the implementations were not run the same way
#'
#' A rate compared across implementations that sampled differently
#' carries a design difference as well as a speed difference.
#'
#' @noRd
warn_benchmark_settings <- function(rows) {
  columns <- c("chains", "iter", "warmup", "thin")
  differing <- columns[vapply(columns, function(nm) {
    values <- unique(rows[[nm]][!is.na(rows[[nm]])])
    length(values) > 1L
  }, logical(1))]
  if (length(differing) > 0L) {
    cli::cli_warn(c(
      "The implementations were not run with the same \\
       {.field {differing}}.",
      i = "A rate then compares a design difference as well as a speed \\
           difference."
    ))
  }
  invisible(rows)
}

#' Warn about runs that did not converge, which are kept in the summary
#' @noRd
warn_benchmark_convergence <- function(rows) {
  bad <- sum(!is.na(rows$converged) & !rows$converged)
  if (bad > 0L) {
    cli::cli_warn(c(
      "{bad} recorded run{?s} did not converge.",
      i = "{cli::qty(bad)}{?It is/They are} kept and counted in \\
           {.field n_not_converged}: \\
           dropping the slow non-converged runs of one implementation \\
           would bias the ratio in its favour."
    ))
  }
  invisible(rows)
}

#' One row per implementation and metric: centre, spread and counts
#' @noRd
benchmark_per_implementation <- function(rows, metrics) {
  implementations <- unique(rows$implementation)
  out <- lapply(implementations, function(name) {
    these <- rows[rows$implementation == name, , drop = FALSE]
    n_error <- sum(nzchar(these$error))
    n_bad <- sum(!is.na(these$converged) & !these$converged)
    dplyr::bind_rows(lapply(metrics, function(metric) {
      values <- these[[metric]]
      values <- values[is.finite(values)]
      if (metric %in% metrics_measured_once()) {
        # one compile, copied into every row: one measurement
        values <- unique(values)
      }
      tibble::tibble(
        implementation = name,
        metric = metric,
        n_runs = nrow(these),
        n_measured = length(values),
        median = median_finite(values),
        # one value has no spread; mad() would report 0
        mad = if (length(values) < 2L) {
          NA_real_
        } else {
          as.double(stats::mad(values))
        },
        min = min_finite(values),
        max = max_finite(values),
        n_not_converged = as.integer(n_bad),
        n_error = as.integer(n_error)
      )
    }))
  })
  dplyr::bind_rows(out)
}

#' The paired ratio of one metric between two implementations
#'
#' Pairs are the data sets each implementation ran, so the ratio is taken
#' within a pair and the interval is a t interval on the log ratios.
#'
#' @noRd
paired_ratio <- function(numerator, denominator, metric, ci_level) {
  if (metric %in% metrics_measured_once()) {
    return(single_measurement_ratio(numerator, denominator, metric))
  }
  pairs <- intersect(numerator$pair, denominator$pair)
  values <- vapply(pairs, function(p) {
    a <- numerator[[metric]][numerator$pair == p]
    b <- denominator[[metric]][denominator$pair == p]
    if (length(a) != 1L || length(b) != 1L) {
      return(NA_real_)
    }
    safe_ratio(a, b)
  }, double(1))
  values <- values[is.finite(values) & values > 0]
  n <- length(values)
  if (n == 0L) {
    return(tibble::tibble(
      ratio = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      n_pairs = 0L, ci_method = "not measured"
    ))
  }
  logs <- log(values)
  ratio <- exp(mean(logs))
  if (n == 1L) {
    return(tibble::tibble(
      ratio = ratio, ci_low = NA_real_, ci_high = NA_real_,
      n_pairs = 1L, ci_method = "single pair"
    ))
  }
  half <- stats::qt(1 - (1 - ci_level) / 2, df = n - 1L) *
    stats::sd(logs) / sqrt(n)
  tibble::tibble(
    ratio = ratio, ci_low = exp(mean(logs) - half),
    ci_high = exp(mean(logs) + half), n_pairs = as.integer(n),
    ci_method = "paired-t (log)"
  )
}

#' The ratio of two quantities each measured once
#'
#' No interval: one measurement over one measurement. Said so in
#' `ci_method`, rather than dressed as a paired interval of width zero.
#'
#' @noRd
single_measurement_ratio <- function(numerator, denominator, metric) {
  value <- unique(numerator[[metric]][is.finite(numerator[[metric]])])
  base <- unique(denominator[[metric]][is.finite(denominator[[metric]])])
  usable <- length(value) == 1L && length(base) == 1L
  if (!usable) {
    return(tibble::tibble(
      ratio = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      n_pairs = 0L, ci_method = "not measured"
    ))
  }
  tibble::tibble(
    ratio = safe_ratio(value, base), ci_low = NA_real_, ci_high = NA_real_,
    n_pairs = 1L, ci_method = "single measurement"
  )
}

#' Every implementation against the reference, every metric
#' @noRd
benchmark_ratios <- function(rows, metrics, reference, ci_level) {
  implementations <- unique(rows$implementation)
  reference <- reference %||% implementations[[1L]]
  if (!reference %in% implementations) {
    cli::cli_abort(
      "{.arg reference} must be one of {.val {implementations}}, \\
       not {.val {reference}}."
    )
  }
  others <- setdiff(implementations, reference)
  if (length(others) == 0L) {
    return(tibble::tibble(
      implementation = character(), reference = character(),
      metric = character(), ratio = double(), ci_low = double(),
      ci_high = double(), n_pairs = integer(), ci_method = character()
    ))
  }
  base <- rows[rows$implementation == reference, , drop = FALSE]
  dplyr::bind_rows(lapply(others, function(name) {
    these <- rows[rows$implementation == name, , drop = FALSE]
    dplyr::bind_rows(lapply(metrics, function(metric) {
      dplyr::bind_cols(
        tibble::tibble(
          implementation = name, reference = reference, metric = metric
        ),
        paired_ratio(these, base, metric, ci_level)
      )
    }))
  }))
}

#' Spread beside the centre, and the ratio with an interval
#'
#' @param object A [benchmark()] result.
#' @param metrics Columns to summarise; by default `us_per_gradient`,
#'   the two sampling ESS rates, `sampling_chain_seconds`,
#'   `wall_seconds` and `compile_seconds`.
#' @param reference The implementation every ratio is taken against, the
#'   first by default.
#' @param ci_level The interval's level.
#' @param ... Not used.
#'
#' @return An object of class `bmmtools_benchmark_summary`: a list with
#'   `per_implementation` (one row per implementation and metric, with
#'   `n_runs`, `n_measured`, `median`, `mad`, `min`, `max`,
#'   `n_not_converged` and `n_error`) and `ratios` (one row per
#'   implementation and metric against `reference`, with `ratio`,
#'   `ci_low`, `ci_high`, `n_pairs` and `ci_method`).
#'
#' @details
#' `mad` is `stats::mad()`, the median absolute deviation scaled by
#' 1.4826 so that it estimates the standard deviation of a normal
#' sample. `compile_seconds` is measured **once per implementation**, so
#' its row reports `n_measured = 1` and its ratio says
#' `"single measurement"` instead of an interval.
#'
#' Every ratio is the other implementation **over** the reference, for
#' all metrics. For a time or a cost, below 1 means the other
#' implementation is faster; for an ESS rate, where more is better,
#' above 1 does.
#'
#' Warmup runs are excluded. Runs that did not converge are **kept** and
#' counted: dropping the slow non-converged runs of one implementation
#' would bias the ratio in its favour. The ratio is the geometric mean of
#' the per-pair ratios with a t interval on the log ratios, which uses
#' the pairing the interleaved design creates; with one pair the bounds
#' are `NA` and `ci_method` says so, so a bare quotient is never
#' reported.
#'
#' @examples
#' \dontrun{
#' s <- summary(bench)
#' s$per_implementation
#' s$ratios
#' }
#'
#' @seealso [benchmark()]
#' @export
summary.bmmtools_benchmark <- function(object, metrics = NULL,
                                       reference = NULL, ci_level = 0.95,
                                       ...) {
  rlang::check_dots_empty()
  metrics <- metrics %||% benchmark_default_metrics()
  metrics <- intersect(metrics, names(object))
  if (length(metrics) == 0L) {
    cli::cli_abort("None of {.arg metrics} is a column of {.arg object}.")
  }
  bad_level <- !is.numeric(ci_level) || length(ci_level) != 1L ||
    is.na(ci_level) || ci_level <= 0 || ci_level >= 1
  if (bad_level) {
    cli::cli_abort("{.arg ci_level} must be a single number in (0, 1).")
  }
  recorded <- tibble::as_tibble(object)
  recorded <- recorded[!is_true_vec(recorded$warmup_run), , drop = FALSE]
  warn_benchmark_settings(recorded)
  warn_benchmark_convergence(recorded)
  structure(
    list(
      per_implementation = benchmark_per_implementation(recorded, metrics),
      ratios = benchmark_ratios(recorded, metrics, reference, ci_level),
      metrics = metrics,
      n_runs = nrow(recorded),
      implementations = unique(recorded$implementation)
    ),
    class = "bmmtools_benchmark_summary"
  )
}

# printing and subsetting -------------------------------------------------

#' Print a benchmark or its summary
#'
#' A header line naming the implementations, the runs and how many were
#' recorded, then the tibble itself.
#'
#' @param x A [benchmark()] result, or a
#'   [summary.bmmtools_benchmark()] one.
#' @param ... Passed to the tibble method.
#' @return `x`, invisibly.
#' @export
print.bmmtools_benchmark <- function(x, ...) {
  settings <- attr(x, "settings")
  implementations <- unique(x$implementation)
  recorded <- sum(!is_true_vec(x$warmup_run))
  cat(
    "# benchmark: ", nrow(x), " runs of ", length(implementations),
    " implementation(s) (", paste(implementations, collapse = ", "), "), ",
    recorded, " recorded\n",
    sep = ""
  )
  if (!is.null(settings$warmup_runs) && settings$warmup_runs > 0L) {
    cat(
      "# ", settings$warmup_runs,
      " warmup run(s) per implementation, marked `warmup_run` and ",
      "excluded from summary()\n",
      sep = ""
    )
  }
  failed <- sum(nzchar(x$error))
  if (failed > 0L) {
    cat("# ", failed, " run(s) failed; see `error`\n", sep = "")
  }
  print(tibble::as_tibble(x), ...)
  invisible(x)
}

#' @rdname print.bmmtools_benchmark
#' @export
print.bmmtools_benchmark_summary <- function(x, ...) {
  cat(
    "# benchmark summary: ", x$n_runs, " recorded run(s), ",
    length(x$implementations), " implementation(s)\n",
    "# metrics: ", paste(x$metrics, collapse = ", "), "\n",
    sep = ""
  )
  print(tibble::as_tibble(x$per_implementation), ...)
  if (nrow(x$ratios) == 0L) {
    cat(
      "# no ratios: one implementation only, so there is nothing to ",
      "compare it with\n",
      sep = ""
    )
  } else {
    cat(
      "# ratios against ", paste(unique(x$ratios$reference), collapse = ", "),
      "\n",
      sep = ""
    )
    print(tibble::as_tibble(x$ratios), ...)
  }
  invisible(x)
}

#' Subset a benchmark
#'
#' A verb that drops a contract column returns a plain tibble; the
#' per-chain and compile tables describe the object as returned, so they
#' are dropped by a subset rather than left to describe rows that are no
#' longer there.
#'
#' @param x A `bmmtools_benchmark` object.
#' @param ... Passed to the tibble method.
#' @return A benchmark while the contract holds, a plain tibble once it
#'   does not.
#' @export
`[.bmmtools_benchmark` <- function(x, ...) {
  out <- demote_if_incomplete(NextMethod(), benchmark_contract_columns())
  if (inherits(out, "bmmtools_benchmark")) {
    attr(out, "chains") <- NULL
    attr(out, "compile") <- NULL
  }
  out
}

#' Keep the class only while the contract holds
#'
#' @param data,template See [dplyr::dplyr_reconstruct()].
#' @return `data`, as a benchmark while it satisfies the contract and as
#'   a plain tibble otherwise.
#' @importFrom dplyr dplyr_reconstruct
#' @exportS3Method dplyr::dplyr_reconstruct
dplyr_reconstruct.bmmtools_benchmark <- function(data, template) {
  if (!all(benchmark_contract_columns() %in% names(data))) {
    return(tibble::as_tibble(data))
  }
  NextMethod()
}
