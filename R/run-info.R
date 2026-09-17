# What ran a fit: the machine and the toolchain, for a study's runtime page.

#' The CPU model string, platform-specific and `NA` where there is none
#'
#' `/proc/cpuinfo`'s `model name` field on Linux,
#' `sysctl -n machdep.cpu.brand_string` on macOS, `NA` on every other
#' platform and whenever the source cannot be read.
#'
#' @noRd
cpu_model_string <- function(sysname = Sys.info()[["sysname"]]) {
  out <- tryCatch(
    {
      if (identical(sysname, "Linux")) {
        lines <- readLines("/proc/cpuinfo", warn = FALSE)
        hit <- grep("^model name\\s*:", lines, value = TRUE)
        if (length(hit) == 0L) {
          NA_character_
        } else {
          trimws(sub("^model name\\s*:", "", hit[[1L]]))
        }
      } else if (identical(sysname, "Darwin")) {
        system2(
          "sysctl", c("-n", "machdep.cpu.brand_string"),
          stdout = TRUE, stderr = FALSE
        )
      } else {
        NA_character_
      }
    },
    error = function(e) NA_character_,
    warning = function(w) NA_character_
  )
  if (length(out) != 1L || is.na(out) || !nzchar(trimws(out))) {
    return(NA_character_)
  }
  trimws(out)
}

#' Total system memory in GB, platform-specific and `NA` where unmeasured
#'
#' `/proc/meminfo`'s `MemTotal` field on Linux, `sysctl -n hw.memsize` on
#' macOS, `NA` on every other platform and whenever the source cannot be
#' read.
#'
#' @noRd
total_memory_gb <- function(sysname = Sys.info()[["sysname"]]) {
  bytes <- tryCatch(
    {
      if (identical(sysname, "Linux")) {
        lines <- readLines("/proc/meminfo", warn = FALSE)
        hit <- grep("^MemTotal:", lines, value = TRUE)
        if (length(hit) == 0L) {
          NA_real_
        } else {
          as.numeric(regmatches(hit, regexpr("[0-9]+", hit))[[1L]]) * 1024
        }
      } else if (identical(sysname, "Darwin")) {
        out <- system2(
          "sysctl", c("-n", "hw.memsize"),
          stdout = TRUE, stderr = FALSE
        )
        as.numeric(trimws(out))
      } else {
        NA_real_
      }
    },
    error = function(e) NA_real_,
    warning = function(w) NA_real_
  )
  if (length(bytes) != 1L || is.na(bytes)) NA_real_ else bytes / 1024^3
}

#' CmdStan's version, or `NA` when cmdstanr is not installed or errors
#' @noRd
cmdstan_version_string <- function(installed = NULL) {
  installed <- installed %||% rlang::is_installed("cmdstanr")
  if (!isTRUE(installed)) {
    return(NA_character_)
  }
  tryCatch(
    as.character(cmdstanr::cmdstan_version()),
    error = function(e) NA_character_
  )
}

#' The machine and package versions this session runs on
#'
#' One row, meant for the runtime table of a study's parameter-recovery
#' page: what ran the fits, not what they estimated. Never errors --- a
#' source it cannot read (no `/proc/cpuinfo`, no CmdStan installation,
#' `cmdstanr` not installed) degrades to `NA` rather than stopping the
#' call, because a runtime record is exactly the thing wanted after a run
#' that partly failed.
#'
#' @return A one-row tibble: `r_version`, `platform`, `hostname`, `cores`
#'   ([parallel::detectCores()]), `cpu_model`, `memory_gb`, the installed
#'   versions of bmm, brms, cmdstanr, posterior and bmmtools,
#'   `cmdstan_version` (from `cmdstanr::cmdstan_version()` when cmdstanr
#'   is installed and does not error), and `timestamp`.
#'
#' @examples
#' run_info()
#'
#' @export
run_info <- function() {
  sys <- Sys.info()
  hostname <- tryCatch(
    unname(sys[["nodename"]]),
    error = function(e) NA_character_
  )
  cores <- tryCatch(
    as.integer(parallel::detectCores()),
    error = function(e) NA_integer_
  )
  tibble::tibble(
    r_version = as.character(getRversion()),
    platform = R.version$platform,
    hostname = if (length(hostname) == 1L) hostname else NA_character_,
    cores = if (length(cores) == 1L) cores else NA_integer_,
    cpu_model = cpu_model_string(sys[["sysname"]]),
    memory_gb = total_memory_gb(sys[["sysname"]]),
    bmm_version = package_version_string("bmm"),
    brms_version = package_version_string("brms"),
    cmdstanr_version = package_version_string("cmdstanr"),
    posterior_version = package_version_string("posterior"),
    bmmtools_version = package_version_string("bmmtools"),
    cmdstan_version = cmdstan_version_string(),
    timestamp = Sys.time()
  )
}
