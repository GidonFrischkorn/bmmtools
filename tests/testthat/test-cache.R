# Tests for fit_cached(), written against
# dev/spec-milestone-2-run-layer.md section 3.
#
# Every test injects the mock fitter from helper-cache.R, so nothing here
# needs bmm, brms or a compiler, and every file lands in a temporary
# directory that withr removes.

cache_call <- function(mock, file, ..., formula = NULL, data = NULL,
                       model = NULL, prior = NULL, refit = "on_change") {
  fit_cached(
    formula = formula %||% fake_bmmformula(a = a ~ 1 + (1 | id)),
    data = data %||% fake_data(),
    model = model %||% fake_model(),
    file = file,
    prior = prior,
    refit = refit,
    ...,
    .fitter = mock$fitter
  )
}

test_that("the first call fits, saves the fit and writes a key", {
  dir <- withr::local_tempdir()
  mock <- mock_fitter()
  file <- file.path(dir, "cell-01")

  expect_message(
    fit <- cache_call(mock, file, seed = 1, iter = 200),
    "Fitting"
  )

  expect_s3_class(fit, "mockfit")
  expect_identical(mock$calls$n, 1L)
  expect_true(file.exists(paste0(file, ".rds")))
  expect_true(file.exists(paste0(file, ".key")))
  cache <- attr(fit, "bmmtools_cache")
  expect_false(cache$reused)
  expect_identical(cache$file, paste0(file, ".rds"))
  expect_type(cache$key, "character")
})

test_that("an identical call reuses the cached fit without fitting", {
  dir <- withr::local_tempdir()
  mock <- mock_fitter()
  file <- file.path(dir, "cell-01")

  suppressMessages(first <- cache_call(mock, file, seed = 1, iter = 200))
  expect_message(
    second <- cache_call(mock, file, seed = 1, iter = 200),
    "Reusing"
  )

  expect_identical(mock$calls$n, 1L)
  expect_identical(second$call_number, first$call_number)
  expect_true(attr(second, "bmmtools_cache")$reused)
})

test_that("a changed sampler argument, prior, data or formula refits", {
  dir <- withr::local_tempdir()
  mock <- mock_fitter()
  file <- file.path(dir, "cell")
  suppressMessages(cache_call(mock, file, seed = 1, iter = 200))

  changes <- list(
    seed = function() cache_call(mock, file, seed = 2, iter = 200),
    iter = function() cache_call(mock, file, seed = 1, iter = 400),
    control = function() {
      cache_call(mock, file,
        seed = 1, iter = 200,
        control = list(adapt_delta = 0.99)
      )
    },
    prior = function() {
      cache_call(mock, file,
        seed = 1, iter = 200,
        prior = data.frame(prior = "normal(0, 1)", nlpar = "a")
      )
    },
    data = function() {
      cache_call(mock, file, seed = 1, iter = 200, data = fake_data(8L))
    },
    formula = function() {
      cache_call(mock, file,
        seed = 1, iter = 200,
        formula = fake_bmmformula(a = a ~ 1)
      )
    }
  )

  for (component in names(changes)) {
    before <- mock$calls$n
    expect_message(changes[[component]](), component, label = component)
    expect_identical(mock$calls$n, before + 1L, label = component)
  }
})

test_that("arguments that do not determine the fit do not refit", {
  dir <- withr::local_tempdir()
  mock <- mock_fitter()
  file <- file.path(dir, "cell")
  suppressMessages(cache_call(mock, file, seed = 1, iter = 200, cores = 2))

  suppressMessages(cache_call(mock, file, seed = 1, iter = 200, cores = 8))
  suppressMessages(cache_call(mock, file, seed = 1, iter = 200, refresh = 0))

  expect_identical(mock$calls$n, 1L)
})

test_that("refit = 'always' refits even when the key matches", {
  dir <- withr::local_tempdir()
  mock <- mock_fitter()
  file <- file.path(dir, "cell")
  suppressMessages(cache_call(mock, file, seed = 1))

  suppressMessages(fit <- cache_call(mock, file, seed = 1, refit = "always"))

  expect_identical(mock$calls$n, 2L)
  expect_identical(fit$call_number, 2L)
  expect_false(attr(fit, "bmmtools_cache")$reused)
})

test_that("refit = 'never' returns the cached fit and warns when it is stale", {
  dir <- withr::local_tempdir()
  mock <- mock_fitter()
  file <- file.path(dir, "cell")
  suppressMessages(cache_call(mock, file, seed = 1))

  expect_warning(
    fit <- cache_call(mock, file, seed = 2, refit = "never"),
    "seed"
  )

  expect_identical(mock$calls$n, 1L)
  expect_identical(fit$call_number, 1L)
  expect_true(attr(fit, "bmmtools_cache")$reused)

  # with a matching key there is nothing to warn about
  expect_no_warning(
    suppressMessages(cache_call(mock, file, seed = 1, refit = "never"))
  )
})

test_that("refit = 'never' still fits when nothing is cached", {
  dir <- withr::local_tempdir()
  mock <- mock_fitter()
  file <- file.path(dir, "fresh")

  suppressMessages(fit <- cache_call(mock, file, seed = 1, refit = "never"))

  expect_identical(mock$calls$n, 1L)
  expect_false(attr(fit, "bmmtools_cache")$reused)
})

test_that("a fit without a key file is refitted and says why", {
  dir <- withr::local_tempdir()
  mock <- mock_fitter()
  file <- file.path(dir, "cell")
  suppressMessages(cache_call(mock, file, seed = 1))
  unlink(paste0(file, ".key"))

  expect_message(cache_call(mock, file, seed = 1), "no cache key")
  expect_identical(mock$calls$n, 2L)
  expect_true(file.exists(paste0(file, ".key")))
})

test_that("brms's own cache arguments are refused", {
  dir <- withr::local_tempdir()
  mock <- mock_fitter()
  file <- file.path(dir, "cell")

  # brms's `file` cannot reach `...`: fit_cached() has a formal of that
  # name, and R refuses a call with two arguments matched to one formal
  expect_error(
    fit_cached(
      fake_bmmformula(a = a ~ 1), fake_data(), fake_model(),
      file = file, file = "other", .fitter = mock$fitter
    ),
    "matched by multiple"
  )
  expect_error(cache_call(mock, file, file_refit = "always"), "file_refit")
  expect_error(cache_call(mock, file, file_compress = TRUE), "file_compress")
  expect_identical(mock$calls$n, 0L)
})

test_that("the .rds extension is added when missing and kept when present", {
  expect_identical(cache_paths("a/b/cell")$rds, "a/b/cell.rds")
  expect_identical(cache_paths("a/b/cell.rds")$rds, "a/b/cell.rds")
  expect_identical(cache_paths("a/b/cell")$key, "a/b/cell.key")
  expect_identical(cache_paths("a/b/cell.rds")$key, "a/b/cell.key")
})

test_that("the formula key ignores the formula's environment", {
  dir <- withr::local_tempdir()
  mock <- mock_fitter()
  file <- file.path(dir, "cell")

  top_level <- fake_bmmformula(a = a ~ 1 + (1 | id), b = b ~ 1)
  build <- function() fake_bmmformula(a = a ~ 1 + (1 | id), b = b ~ 1)

  expect_false(identical(
    environment(top_level$a), environment(build()$a)
  ))
  expect_identical(formula_key(top_level), formula_key(build()))

  suppressMessages(cache_call(mock, file, formula = top_level, seed = 1))
  suppressMessages(cache_call(mock, file, formula = build(), seed = 1))
  expect_identical(mock$calls$n, 1L)
})

test_that("the key file lists one hash per component and the overall key", {
  dir <- withr::local_tempdir()
  mock <- mock_fitter()
  file <- file.path(dir, "cell")
  suppressMessages(fit <- cache_call(mock, file, seed = 1, chains = 2))

  lines <- readLines(paste0(file, ".key"))
  parsed <- read_cache_key(paste0(file, ".key"))

  expect_true(all(grepl("^[A-Za-z_]+=[0-9a-f]+$", lines)))
  expect_true(all(c(
    "formula", "data", "model", "prior", "seed", "chains", "iter",
    "warmup", "thin", "control", "init", "backend", "sample_prior",
    "algorithm", "stanvars", "bmm", "brms"
  ) %in% names(parsed$components)))
  expect_identical(parsed$key, attr(fit, "bmmtools_cache")$key)
})

test_that("a missing directory is created", {
  dir <- withr::local_tempdir()
  mock <- mock_fitter()
  file <- file.path(dir, "deeper", "still", "cell")

  suppressMessages(cache_call(mock, file, seed = 1))

  expect_true(file.exists(paste0(file, ".rds")))
})

test_that("the fitter receives the formula, data, model, prior and dots", {
  dir <- withr::local_tempdir()
  mock <- mock_fitter()
  file <- file.path(dir, "cell")
  prior <- data.frame(prior = "normal(0, 1)", nlpar = "a")

  suppressMessages(cache_call(mock, file,
    prior = prior, seed = 3, cores = 4
  ))

  last <- mock$calls$last
  expect_s3_class(last$formula, "bmmformula")
  expect_identical(last$data, fake_data())
  expect_s3_class(last$model, "bmmodel")
  expect_identical(last$prior, prior)
  expect_identical(last$dots$seed, 3)
  expect_identical(last$dots$cores, 4)
})

test_that(".fitter must be a function", {
  dir <- withr::local_tempdir()
  expect_error(
    fit_cached(
      fake_bmmformula(a = a ~ 1), fake_data(), fake_model(),
      file = file.path(dir, "cell"), .fitter = 1
    ),
    "function"
  )
})

test_that("refit must be one of the three options", {
  dir <- withr::local_tempdir()
  mock <- mock_fitter()
  expect_error(
    cache_call(mock, file.path(dir, "cell"), refit = "sometimes"),
    "refit"
  )
})

# the small helpers ------------------------------------------------------

test_that("cache_paths rejects anything but a single non-empty path", {
  expect_error(cache_paths(1), "file")
  expect_error(cache_paths(c("a", "b")), "file")
  expect_error(cache_paths(""), "file")
  expect_error(cache_paths(NA_character_), "file")
})

test_that("formula_key passes non-formula elements through", {
  expect_identical(formula_key("abc"), "abc")
  expect_identical(formula_key(3), 3)
  key <- formula_key(list(a = y ~ x, n = 2))
  expect_identical(key$a, "y ~ x")
  expect_identical(key$n, 2)
})

test_that("a package that is not installed is recorded as such", {
  expect_identical(
    package_version_string("no.such.package.bmmtools"),
    "not installed"
  )
})

test_that("a key file with no component lines reads as NULL", {
  path <- withr::local_tempfile(fileext = ".key")
  writeLines(c("", "garbage"), path)
  expect_null(read_cache_key(path))
  expect_null(read_cache_key(file.path(tempdir(), "does-not-exist.key")))
})

test_that("refit = 'never' warns when the key file is gone", {
  dir <- withr::local_tempdir()
  mock <- mock_fitter()
  file <- file.path(dir, "cell")
  suppressMessages(cache_call(mock, file, seed = 1))
  unlink(paste0(file, ".key"))

  expect_warning(
    fit <- cache_call(mock, file, seed = 1, refit = "never"),
    "no cache key"
  )
  expect_identical(mock$calls$n, 1L)
  expect_true(attr(fit, "bmmtools_cache")$reused)
})

# review findings (session 4) -------------------------------------------

test_that("an init function enters the key by its text, not its environment", {
  dir <- withr::local_tempdir()
  mock <- mock_fitter()
  file <- file.path(dir, "cell")

  make_init <- function() function() list(kappa = 0)
  suppressMessages(cache_call(mock, file, seed = 1, init = make_init()))
  # a closure built in another environment with the same text matches
  suppressMessages(cache_call(mock, file, seed = 1, init = make_init()))
  expect_identical(mock$calls$n, 1L)

  # a closure with different text refits
  expect_message(
    cache_call(mock, file, seed = 1, init = function() list(kappa = 1)),
    "init"
  )
  expect_identical(mock$calls$n, 2L)

  # a non-function init hashes as it is
  expect_identical(function_key(0), 0)
  expect_identical(function_key("random"), "random")
})

test_that("the toolchain of the backend in use is part of the key", {
  dir <- withr::local_tempdir()
  mock <- mock_fitter()
  file <- file.path(dir, "cell")
  suppressMessages(cache_call(mock, file, backend = "mock"))

  parsed <- read_cache_key(paste0(file, ".key"))
  expect_true("toolchain" %in% names(parsed$components))
  expect_identical(backend_version_string("mock"), "mock")
  expect_match(backend_version_string("rstan"), "^rstan ")
  expect_match(backend_version_string("cmdstanr"), "^cmdstanr ")
  expect_type(backend_version_string(NULL), "character")
})

test_that("the directory exists before the fitter runs", {
  dir <- withr::local_tempdir()
  file <- file.path(dir, "new", "cell")
  seen <- new.env()
  fitter <- function(formula, data, model, prior = NULL, ...) {
    seen$dir_existed <- dir.exists(file.path(dir, "new"))
    structure(list(), class = "mockfit")
  }
  suppressMessages(fit_cached(
    fake_bmmformula(a = a ~ 1), fake_data(), fake_model(),
    file = file, .fitter = fitter
  ))
  expect_true(seen$dir_existed)
})

test_that("a truncated cached fit is reported, not parsed", {
  dir <- withr::local_tempdir()
  mock <- mock_fitter()
  file <- file.path(dir, "cell")
  suppressMessages(cache_call(mock, file, seed = 1))
  writeLines("not an rds", paste0(file, ".rds"))

  expect_error(
    suppressMessages(cache_call(mock, file, seed = 1, refit = "never")),
    "could not be read"
  )
  # the key still matches, so on_change would reuse it too and must say
  # the same thing; "always" is the documented way out
  expect_error(
    suppressMessages(cache_call(mock, file, seed = 1)),
    "could not be read"
  )
  suppressMessages(fit <- cache_call(mock, file, seed = 1, refit = "always"))
  expect_s3_class(fit, "mockfit")
  expect_identical(mock$calls$n, 2L)
})

test_that("writes go through a temporary name and leave no stray file", {
  dir <- withr::local_tempdir()
  mock <- mock_fitter()
  file <- file.path(dir, "cell")
  suppressMessages(cache_call(mock, file, seed = 1))
  expect_setequal(list.files(dir), c("cell.rds", "cell.key"))
  expect_setequal(
    list.files(dir, all.files = TRUE, no.. = TRUE),
    c("cell.rds", "cell.key")
  )
})
