# Tests for run_info(), spec section 9.4.

test_that("run_info returns one row with the expected columns and types", {
  out <- run_info()
  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 1L)
  expect_named(out, c(
    "r_version", "platform", "hostname", "cores", "cpu_model",
    "memory_gb", "bmm_version", "brms_version", "cmdstanr_version",
    "posterior_version", "bmmtools_version", "cmdstan_version", "timestamp"
  ))
  expect_type(out$r_version, "character")
  expect_type(out$platform, "character")
  expect_type(out$hostname, "character")
  expect_type(out$cores, "integer")
  expect_type(out$cpu_model, "character")
  expect_type(out$memory_gb, "double")
  expect_type(out$bmm_version, "character")
  expect_type(out$brms_version, "character")
  expect_type(out$cmdstanr_version, "character")
  expect_type(out$posterior_version, "character")
  expect_type(out$bmmtools_version, "character")
  expect_type(out$cmdstan_version, "character")
  expect_s3_class(out$timestamp, "POSIXct")

  expect_equal(
    out$bmmtools_version, as.character(utils::packageVersion("bmmtools"))
  )
})

test_that("run_info never errors", {
  expect_no_error(run_info())
})

test_that("cpu_model_string and total_memory_gb degrade off Linux and macOS", {
  expect_true(is.na(cpu_model_string("SunOS")))
  expect_true(is.na(total_memory_gb("SunOS")))
})

test_that("cpu_model_string and total_memory_gb degrade when a source errors", {
  local_mocked_bindings(
    readLines = function(...) stop("no such file"),
    .package = "base"
  )
  expect_true(is.na(cpu_model_string("Linux")))
  expect_true(is.na(total_memory_gb("Linux")))
})

test_that("the Linux branches read /proc as they are written to", {
  # this machine is macOS, so the Linux paths are reachable only with the
  # file contents handed in
  local_mocked_bindings(
    readLines = function(con, ...) {
      if (identical(con, "/proc/cpuinfo")) {
        return(c("processor\t: 0", "model name\t: AMD EPYC 7763 64-Core"))
      }
      c("MemTotal:       16384 kB", "MemFree:         8192 kB")
    },
    .package = "base"
  )
  expect_equal(cpu_model_string("Linux"), "AMD EPYC 7763 64-Core")
  # 16384 kB is 16384 * 1024 bytes, which is 0.015625 GB
  expect_equal(total_memory_gb("Linux"), 0.015625)
})

test_that("the Linux branches give NA when the field is absent", {
  local_mocked_bindings(
    readLines = function(...) "nothing this function is looking for",
    .package = "base"
  )
  expect_true(is.na(cpu_model_string("Linux")))
  expect_true(is.na(total_memory_gb("Linux")))
})

test_that("cmdstan_version_string is NA when cmdstanr is not installed", {
  expect_true(is.na(cmdstan_version_string(installed = FALSE)))
})

test_that("cmdstan_version_string is NA when cmdstan_version() errors", {
  skip_if_not_installed("cmdstanr")
  local_mocked_bindings(
    cmdstan_version = function(...) stop("no CmdStan found"),
    .package = "cmdstanr"
  )
  expect_true(is.na(cmdstan_version_string(installed = TRUE)))
})
