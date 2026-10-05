# Tests for the benchmark measurement functions, written against
# local/dev/spec-milestone-15-benchmark.md (decisions D69-D77).
#
# The rules run on two hand-built tables --- the per-chain elapsed times
# and the per-chain leapfrog counts --- so they are tested without a fit
# first. The fixture tests exercise the brmsfit method and skip when brms
# is absent. No Stan is compiled anywhere.

# the core, on hand-built parts -------------------------------------------

test_that("cost per gradient is per chain, and the median over chains", {
  row <- benchmark_from_parts(
    elapsed = fake_elapsed(warmup = c(2, 2, 2), sample = c(1, 2, 3)),
    leapfrog = fake_leapfrog(c(1000, 1000, 1000)),
    diagnostics = fake_bench_diagnostics(), nuts = NULL,
    settings = fake_bench_settings(chains = 3L)
  )
  # 1e6 * c(1, 2, 3) / 1000 = 1000, 2000, 3000 us per gradient
  expect_equal(row$us_per_gradient, 2000)
  expect_equal(row$us_per_gradient_min, 1000)
  expect_equal(row$us_per_gradient_max, 3000)
  chains <- attr(row, "chains")
  expect_equal(chains$us_per_gradient, c(1000, 2000, 3000))
})

test_that("cost per gradient does not involve wall time (D69, D70)", {
  parallel <- benchmark_from_parts(
    elapsed = fake_elapsed(warmup = c(1, 1), sample = c(1, 1)),
    leapfrog = fake_leapfrog(c(1000, 1000)),
    diagnostics = fake_bench_diagnostics(), nuts = NULL,
    wall_seconds = 2, settings = fake_bench_settings()
  )
  sequential <- benchmark_from_parts(
    elapsed = fake_elapsed(warmup = c(1, 1), sample = c(1, 1)),
    leapfrog = fake_leapfrog(c(1000, 1000)),
    diagnostics = fake_bench_diagnostics(), nuts = NULL,
    wall_seconds = 4, settings = fake_bench_settings()
  )
  # the same chains: two run at once, then the same two in sequence
  expect_equal(parallel$us_per_gradient, 1000)
  expect_equal(parallel$us_per_gradient, sequential$us_per_gradient)
  expect_equal(parallel$sampling_chain_seconds, sequential$sampling_chain_seconds)
  expect_equal(parallel$wall_seconds, 2)
  expect_equal(sequential$wall_seconds, 4)
})

test_that("chain seconds are summed over chains, not maximised", {
  row <- benchmark_from_parts(
    elapsed = fake_elapsed(warmup = c(2, 3), sample = c(1, 4)),
    leapfrog = fake_leapfrog(c(1000, 1000)),
    diagnostics = fake_bench_diagnostics(), nuts = NULL,
    settings = fake_bench_settings()
  )
  expect_equal(row$warmup_chain_seconds, 5)
  expect_equal(row$sampling_chain_seconds, 5)
  expect_equal(row$total_chain_seconds, 10)
  expect_equal(row$max_chain_seconds, 7)
})

test_that("overhead is wall time minus the slowest chain, and NA without wall", {
  with_wall <- benchmark_from_parts(
    elapsed = fake_elapsed(warmup = c(2, 2), sample = c(1, 1)),
    leapfrog = fake_leapfrog(c(1000, 1000)),
    diagnostics = fake_bench_diagnostics(), nuts = NULL,
    wall_seconds = 10, settings = fake_bench_settings()
  )
  expect_equal(with_wall$overhead_seconds, 7)
  without <- benchmark_from_parts(
    elapsed = fake_elapsed(warmup = c(2, 2), sample = c(1, 1)),
    leapfrog = fake_leapfrog(c(1000, 1000)),
    diagnostics = fake_bench_diagnostics(), nuts = NULL,
    settings = fake_bench_settings()
  )
  expect_true(is.na(without$overhead_seconds))
})

test_that("ESS rates divide by chain seconds, per phase (D70, D71)", {
  row <- benchmark_from_parts(
    elapsed = fake_elapsed(warmup = c(1, 1), sample = c(2, 2)),
    leapfrog = fake_leapfrog(c(1000, 1000)),
    diagnostics = fake_bench_diagnostics(
      ess_bulk = c(800, 400),
      ess_tail = c(600, 200)
    ),
    nuts = NULL, settings = fake_bench_settings()
  )
  expect_equal(row$ess_bulk_min, 400)
  expect_equal(row$ess_tail_min, 200)
  expect_equal(row$ess_bulk_per_chain_sec_sampling, 400 / 4)
  expect_equal(row$ess_tail_per_chain_sec_sampling, 200 / 4)
  expect_equal(row$ess_bulk_per_chain_sec_total, 400 / 6)
  expect_equal(row$ess_tail_per_chain_sec_total, 200 / 6)
})

test_that("a rate is NA, never Inf, when its denominator is zero or NA", {
  row <- benchmark_from_parts(
    elapsed = fake_elapsed(warmup = c(0, 0), sample = c(0, 0)),
    leapfrog = fake_leapfrog(c(0, 0)),
    diagnostics = fake_bench_diagnostics(), nuts = NULL,
    settings = fake_bench_settings()
  )
  expect_true(is.na(row$ess_bulk_per_chain_sec_sampling))
  expect_true(is.na(row$ess_bulk_per_chain_sec_total))
  expect_true(is.na(row$us_per_gradient))
})

test_that("constants are dropped from the ESS minimum (D75)", {
  row <- benchmark_from_parts(
    elapsed = fake_elapsed(warmup = c(1, 1), sample = c(1, 1)),
    leapfrog = fake_leapfrog(c(1000, 1000)),
    diagnostics = fake_bench_diagnostics(
      rhat = c(1.00, 1.01, NA, NA), ess_bulk = c(800, 600, NA, NA),
      ess_tail = c(700, 500, NA, NA)
    ),
    nuts = NULL, settings = fake_bench_settings()
  )
  expect_equal(row$ess_bulk_min, 600)
  expect_equal(row$n_variables, 2L)
})

test_that("warmup work is NA unless warmup draws were saved (D71)", {
  without <- benchmark_from_parts(
    elapsed = fake_elapsed(warmup = c(2, 2), sample = c(1, 1)),
    leapfrog = fake_leapfrog(c(1000, 1000)),
    diagnostics = fake_bench_diagnostics(), nuts = NULL,
    settings = fake_bench_settings()
  )
  expect_true(is.na(without$leapfrog_warmup))
  expect_true(is.na(without$us_per_gradient_warmup))

  with_warmup <- benchmark_from_parts(
    elapsed = fake_elapsed(warmup = c(2, 2), sample = c(1, 1)),
    leapfrog = fake_leapfrog(c(1000, 1000)),
    leapfrog_warmup = fake_leapfrog(c(4000, 4000)),
    diagnostics = fake_bench_diagnostics(), nuts = NULL,
    settings = fake_bench_settings()
  )
  expect_equal(with_warmup$leapfrog_warmup, 8000)
  expect_equal(with_warmup$us_per_gradient_warmup, 500)
})

test_that("a fit without sampler diagnostics keeps its times and reports NA work", {
  row <- benchmark_from_parts(
    elapsed = fake_elapsed(warmup = c(2, 2), sample = c(1, 1)),
    leapfrog = NULL,
    diagnostics = fake_bench_diagnostics(), nuts = NULL,
    wall_seconds = 5, settings = fake_bench_settings()
  )
  expect_equal(row$sampling_chain_seconds, 2)
  expect_equal(row$wall_seconds, 5)
  expect_true(is.na(row$leapfrog_sampling))
  expect_true(is.na(row$us_per_gradient))
  expect_true(is.na(row$n_divergent))
})

test_that("the row satisfies the contract and carries its per-chain table", {
  row <- benchmark_from_parts(
    elapsed = fake_elapsed(warmup = c(1, 1), sample = c(1, 1)),
    leapfrog = fake_leapfrog(c(1000, 1000)),
    diagnostics = fake_bench_diagnostics(), nuts = NULL,
    settings = fake_bench_settings()
  )
  expect_s3_class(row, "bmmtools_benchmark")
  expect_equal(nrow(row), 1L)
  contract <- benchmark_contract()
  expect_true(all(names(contract) %in% names(row)))
  for (nm in names(contract)) {
    expect_identical(typeof(row[[nm]]), contract[[nm]], info = nm)
  }
  chains <- attr(row, "chains")
  expect_equal(nrow(chains), 2L)
  expect_true(all(
    c(
      "chain", "warmup_seconds", "sample_seconds", "leapfrog_sampling",
      "us_per_gradient"
    ) %in% names(chains)
  ))
})

test_that("one chain claims no spread", {
  row <- benchmark_from_parts(
    elapsed = fake_elapsed(warmup = 2, sample = 1),
    leapfrog = fake_leapfrog(1000),
    diagnostics = fake_bench_diagnostics(), nuts = NULL,
    settings = fake_bench_settings(chains = 1L)
  )
  expect_equal(row$us_per_gradient, 1000)
  expect_equal(row$us_per_gradient_min, row$us_per_gradient)
  expect_equal(row$us_per_gradient_max, row$us_per_gradient)
})

test_that("compile seconds and their source are carried through", {
  row <- benchmark_from_parts(
    elapsed = fake_elapsed(warmup = 1, sample = 1),
    leapfrog = fake_leapfrog(1000),
    diagnostics = fake_bench_diagnostics(), nuts = NULL,
    compile_seconds = 4.68, compile_source = "forced",
    settings = fake_bench_settings()
  )
  expect_equal(row$compile_seconds, 4.68)
  expect_identical(row$compile_source, "forced")
})

# benchmark_metrics() on the committed fixture ----------------------------

test_that("benchmark_metrics() measures the fixture (no Stan)", {
  skip_if_not_installed("brms")
  skip_if_not_installed("rstan")
  fit <- mixture2p_fit()
  row <- benchmark_metrics(fit, wall_seconds = 9)

  expect_s3_class(row, "bmmtools_benchmark")
  # the fixture's own numbers, measured 2026-10-05
  expect_equal(row$warmup_chain_seconds, 3.248, tolerance = 1e-6)
  expect_equal(row$sampling_chain_seconds, 2.160, tolerance = 1e-6)
  expect_equal(row$leapfrog_sampling, 15232)
  # 28 variables, less 6 constants, less lp__ and lprior (D75)
  expect_equal(row$n_variables, 20L)
  expect_equal(row$ess_bulk_min, 220.34, tolerance = 1e-2)
  expect_equal(row$ess_tail_min, 193.72, tolerance = 1e-2)
  expect_identical(row$backend, "cmdstanr")
  expect_equal(row$chains, 2L)
  expect_equal(row$iter, 1000L)
  expect_equal(row$warmup, 500L)
  expect_true(is.na(row$leapfrog_warmup))
})

test_that("the fixture's gradient cost sits inside its per-chain range", {
  skip_if_not_installed("brms")
  skip_if_not_installed("rstan")
  row <- benchmark_metrics(mixture2p_fit())
  expect_true(is.finite(row$us_per_gradient))
  expect_gt(row$us_per_gradient, 0)
  expect_gte(row$us_per_gradient, row$us_per_gradient_min)
  expect_lte(row$us_per_gradient, row$us_per_gradient_max)
  chains <- attr(row, "chains")
  expect_equal(nrow(chains), 2L)
  expect_equal(
    chains$us_per_gradient,
    1e6 * chains$sample_seconds / chains$leapfrog_sampling
  )
})

test_that("benchmark_metrics() refuses what is not a fit", {
  expect_error(benchmark_metrics(1), class = "rlang_error")
  expect_error(benchmark_metrics(1), "brmsfit")
})

# benchmark_compile() ----------------------------------------------------

test_that("benchmark_compile() times the compiler it is given", {
  compiler <- function(code, backend, force, ...) {
    Sys.sleep(0.05)
    list(exe = "/tmp/fake-exe")
  }
  out <- benchmark_compile("// code", backend = "cmdstanr", .compiler = compiler)
  expect_equal(nrow(out), 1L)
  expect_gte(out$compile_seconds, 0.05)
  expect_identical(out$backend, "cmdstanr")
  expect_true(out$forced)
  expect_identical(out$code_hash, rlang::hash("// code"))
  expect_identical(out$exe, "/tmp/fake-exe")
})

test_that("benchmark_compile() validates code and backend", {
  expect_error(benchmark_compile(1, .compiler = function(...) list()), "code")
  expect_error(
    benchmark_compile("// code", backend = "jags", .compiler = function(...) list()),
    "cmdstanr"
  )
})

# the pieces that read a fit, and the validation branches ----------------

test_that("an object with no elapsed times gives an NA matrix of its chains", {
  stub <- structure(list(fit = NULL), class = c("bmmfit", "brmsfit"))
  out <- fit_elapsed_time(stub)
  expect_true(is.matrix(out))
  expect_identical(colnames(out), c("warmup", "sample"))
  expect_equal(nrow(out), 1L)
  expect_true(all(is.na(out)))
})

test_that("leapfrog_from_nuts() returns NULL for nothing usable", {
  expect_null(leapfrog_from_nuts(NULL))
  expect_null(leapfrog_from_nuts(data.frame(a = 1)))
  expect_null(leapfrog_from_nuts(data.frame(
    Chain = 1L, Parameter = "divergent__", Value = 0
  )))
  table <- leapfrog_from_nuts(data.frame(
    Chain = c(1L, 1L, 2L), Parameter = "n_leapfrog__", Value = c(3, 4, 5)
  ))
  expect_equal(table$chain, c(1L, 2L))
  expect_equal(table$leapfrog, c(7, 5))
  expect_equal(table$draws, c(2L, 1L))
})

test_that("fit_saved_warmup() is 0 when nothing says otherwise", {
  expect_equal(fit_saved_warmup(structure(list(), class = "brmsfit")), 0L)
})

test_that("the warmup split is the total minus the post-warmup draws", {
  total <- fake_leapfrog(c(3000, 3200), draws = 600L)
  sampling <- fake_leapfrog(c(1000, 1100), draws = 300L)
  warmup <- leapfrog_warmup_from_totals(total, sampling)
  expect_equal(warmup$leapfrog, c(2000, 2100))
  expect_equal(warmup$draws, c(300L, 300L))
  expect_null(leapfrog_warmup_from_totals(NULL, sampling))
  expect_null(leapfrog_warmup_from_totals(total, NULL))
  # nothing extra was saved: the warmup phase is unknown, not zero
  expect_null(leapfrog_warmup_from_totals(sampling, sampling))
})

test_that("asking a fit without warmup draws for them yields nothing", {
  stub <- structure(list(fit = NULL), class = c("bmmfit", "brmsfit"))
  out <- fit_leapfrog_tables(stub, inc_warmup = TRUE)
  expect_null(out$sampling)
  expect_null(out$warmup)
})

test_that("the fixture's warmup leapfrog steps are refused, not invented", {
  skip_if_not_installed("brms")
  skip_if_not_installed("rstan")
  fit <- mixture2p_fit()
  expect_equal(fit_saved_warmup(fit), 0L)
  # inc_warmup = TRUE cannot invent draws the fit never saved, and the
  # row says NA rather than a warmup cost of zero
  row <- benchmark_metrics(fit, inc_warmup = TRUE)
  expect_true(is.na(row$leapfrog_warmup))
  expect_true(is.na(row$us_per_gradient_warmup))
})

test_that("benchmark_compile() validates force and .compiler", {
  expect_error(
    benchmark_compile("// c", force = NA, .compiler = function(...) list()),
    "force"
  )
  expect_error(benchmark_compile("// c", .compiler = 1), ".compiler")
})

test_that("benchmark_compile() can measure a reuse instead of a compile", {
  out <- benchmark_compile(
    "// c", force = FALSE, .compiler = function(...) list()
  )
  expect_false(out$forced)
  expect_true(is.na(out$exe))
})

test_that("coerce_like() keeps a column's type", {
  expect_identical(coerce_like(2, NA_integer_), 2L)
  expect_identical(coerce_like("2", NA_real_), 2)
  expect_identical(coerce_like(1, NA), TRUE)
  expect_identical(coerce_like(2, NA_character_), "2")
})

test_that("fill_na_columns() fills only what is unknown", {
  row <- tibble::tibble(a = NA_integer_, b = 3L)
  out <- fill_na_columns(row, list(a = 2, b = 9, c = 1))
  expect_identical(out$a, 2L)
  expect_identical(out$b, 3L)
  expect_false("c" %in% names(out))
})

test_that("a paired ratio needs one value per pair", {
  rows <- function(pair, value) {
    tibble::tibble(pair = pair, us_per_gradient = value)
  }
  doubled <- paired_ratio(
    rows(c(1L, 1L), c(2, 3)), rows(1L, 1), "us_per_gradient", 0.95
  )
  expect_equal(doubled$n_pairs, 0L)
  expect_identical(doubled$ci_method, "not measured")
})

# what the review of 15.2 found ------------------------------------------

test_that("lp__ and lprior are not model parameters (D75)", {
  row <- benchmark_from_parts(
    elapsed = fake_elapsed(warmup = c(1, 1), sample = c(1, 1)),
    leapfrog = fake_leapfrog(c(1000, 1000)),
    diagnostics = fake_bench_diagnostics(
      rhat = c(1.00, 1.01, 1.00, 1.00),
      ess_bulk = c(800, 600, 100, 90),
      ess_tail = c(700, 500, 120, 110),
      variable = c("b_a", "b_b", "lprior", "lp__")
    ),
    nuts = NULL, settings = fake_bench_settings()
  )
  # the sampler's own variables have a finite rhat, so the rule that
  # drops constants keeps them; on the fixture lp__ was the minimum
  expect_equal(row$ess_bulk_min, 600)
  expect_equal(row$ess_tail_min, 500)
  expect_equal(row$n_variables, 2L)
})

test_that("thinned draws report the count but no gradient cost (D78)", {
  thinned <- benchmark_from_parts(
    elapsed = fake_elapsed(warmup = c(1, 1), sample = c(1, 1)),
    leapfrog = fake_leapfrog(c(1000, 1000), draws = 250L),
    leapfrog_warmup = fake_leapfrog(c(1000, 1000), draws = 250L),
    diagnostics = fake_bench_diagnostics(), nuts = NULL,
    settings = fake_bench_settings(thin = 2L)
  )
  expect_equal(thinned$leapfrog_sampling, 2000)
  expect_equal(thinned$sampling_chain_seconds, 2)
  expect_true(is.na(thinned$us_per_gradient))
  expect_true(is.na(thinned$us_per_gradient_min))
  expect_true(is.na(thinned$us_per_gradient_warmup))

  unthinned <- benchmark_from_parts(
    elapsed = fake_elapsed(warmup = c(1, 1), sample = c(1, 1)),
    leapfrog = fake_leapfrog(c(1000, 1000)),
    diagnostics = fake_bench_diagnostics(), nuts = NULL,
    settings = fake_bench_settings(thin = 1L)
  )
  expect_equal(unthinned$us_per_gradient, 1000)
})

test_that("a sum over chains is NA when one chain is unmeasured", {
  row <- benchmark_from_parts(
    elapsed = fake_elapsed(warmup = c(1, NA), sample = c(1, 2)),
    leapfrog = fake_leapfrog(c(1000, 1000)),
    diagnostics = fake_bench_diagnostics(), nuts = NULL,
    settings = fake_bench_settings()
  )
  expect_true(is.na(row$warmup_chain_seconds))
  expect_true(is.na(row$total_chain_seconds))
  expect_equal(row$sampling_chain_seconds, 3)
  expect_true(is.na(row$ess_bulk_per_chain_sec_total))
  expect_equal(row$ess_bulk_per_chain_sec_sampling, 600 / 3)

  partial_work <- benchmark_from_parts(
    elapsed = fake_elapsed(warmup = c(1, 1), sample = c(1, 1)),
    leapfrog = tibble::tibble(
      chain = 1:2, leapfrog = c(1000, NA), draws = c(500L, 500L)
    ),
    diagnostics = fake_bench_diagnostics(), nuts = NULL,
    settings = fake_bench_settings()
  )
  expect_true(is.na(partial_work$leapfrog_sampling))
  expect_equal(partial_work$us_per_gradient, 1000)
})

test_that("a row for a failed run has one chain per chain asked for", {
  row <- empty_benchmark_row(settings = fake_bench_settings(chains = 4L))
  expect_equal(nrow(attr(row, "chains")), 4L)
  expect_true(all(is.na(attr(row, "chains")$sample_seconds)))
  expect_equal(nrow(attr(empty_benchmark_row(), "chains")), 1L)
})
