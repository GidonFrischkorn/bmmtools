# The record each grid call leaves in its directory: which process ran
# which cells, on what, and when.

#' The scheduler job this process runs in, read from its environment
#'
#' SLURM, PBS (Pro and Torque), LSF and SGE, in that order of detection;
#' `NA` throughout on a machine without one, such as a cloud VM. A core
#' count that does not read as a positive whole number is `NA`.
#'
#' @return A one-row tibble: `scheduler`, `job_id`, `array_job_id`,
#'   `array_task_id`, `cores_granted`.
#' @noRd
grid_job_info <- function() {
  env <- function(name) {
    value <- Sys.getenv(name, NA_character_)
    if (is.na(value) || !nzchar(value)) NA_character_ else value
  }
  first <- function(names) {
    for (name in names) {
      if (!is.na(env(name))) {
        return(env(name))
      }
    }
    NA_character_
  }
  count <- function(value) {
    # an unreadable count is NA by design, not a coercion warning
    out <- suppressWarnings(as.integer(value))
    if (length(out) != 1L || is.na(out) || out < 1L) NA_integer_ else out
  }
  # outside an array, LSF's index is 0 and SGE's is "undefined"
  task <- function(value, none) {
    if (is.na(value) || value %in% none) NA_character_ else value
  }
  out <- if (!is.na(env("SLURM_JOB_ID"))) {
    list(
      "slurm", env("SLURM_JOB_ID"), env("SLURM_ARRAY_JOB_ID"),
      env("SLURM_ARRAY_TASK_ID"), count(env("SLURM_CPUS_PER_TASK"))
    )
  } else if (!is.na(env("PBS_JOBID"))) {
    list(
      "pbs", env("PBS_JOBID"), NA_character_,
      first(c("PBS_ARRAY_INDEX", "PBS_ARRAYID")),
      count(first(c("NCPUS", "PBS_NP")))
    )
  } else if (!is.na(env("LSB_JOBID"))) {
    list(
      "lsf", env("LSB_JOBID"), NA_character_,
      task(env("LSB_JOBINDEX"), "0"), count(env("LSB_DJOB_NUMPROC"))
    )
  } else if (!is.na(env("JOB_ID")) &&
               !all(is.na(c(env("SGE_TASK_ID"), env("NSLOTS"))))) {
    list(
      "sge", env("JOB_ID"), NA_character_,
      task(env("SGE_TASK_ID"), "undefined"), count(env("NSLOTS"))
    )
  } else {
    list(
      NA_character_, NA_character_, NA_character_, NA_character_,
      NA_integer_
    )
  }
  names(out) <- c(
    "scheduler", "job_id", "array_job_id", "array_task_id", "cores_granted"
  )
  tibble::as_tibble(out)
}

#' Cell numbers as short text: runs of consecutive cells as `a-b`
#' @noRd
cells_text <- function(cells) {
  if (length(cells) == 0L) {
    return("")
  }
  cells <- sort(cells)
  starts <- c(TRUE, diff(cells) != 1L)
  runs <- split(cells, cumsum(starts))
  paste(vapply(runs, function(r) {
    if (length(r) == 1L) {
      as.character(r)
    } else {
      paste0(r[[1L]], "-", r[[length(r)]])
    }
  }, character(1)), collapse = ",")
}

#' Open the record of this call, just before its first fit (D117)
#'
#' Written at once with `ended` `NA`, so that a process stopped while it
#' fits still says that it ran and which cells it held. A call that
#' errors before its cells (a failed preflight, a refused slice) never
#' gets this far and leaves no record.
#'
#' @param started When the call started, so that the record's span
#'   includes the preflight and the prior fits.
#' @return The record, with its `path`, for `close_run_record()`.
#' @noRd
open_run_record <- function(dir, call, cells, started = Sys.time()) {
  info <- run_info()
  # counted before tibble(): there `cells` would be the column just made
  n_cells <- length(cells)
  row <- tibble::tibble(
    call = call,
    cells = cells_text(cells),
    n_cells = n_cells,
    n_reused = NA_integer_,
    n_failed = NA_integer_,
    started = started,
    ended = as.POSIXct(NA),
    pid = Sys.getpid(),
    grid_job_info(),
    info
  )
  path <- run_record_path(dir, started, info$hostname, Sys.getpid())
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write_atomic(path, function(tmp) saveRDS(row, tmp))
  list(path = path, row = row)
}

#' Where a call's record goes: its start, host and process id
#'
#' The host is reduced to bytes safe in a file name, whatever its
#' encoding.
#'
#' @noRd
run_record_path <- function(dir, started, host, pid) {
  host <- if (is.na(host)) "host" else host
  file.path(dir, "run-info", sprintf(
    "%s-%s-%d.rds", format(started, "%Y%m%d-%H%M%OS3"),
    gsub("[^A-Za-z0-9._-]", "_", host, useBytes = TRUE), pid
  ))
}

#' Close the record: the cells reused and failed, and the time it ended
#' @noRd
close_run_record <- function(record, runs = list()) {
  row <- record$row
  row$n_reused <- sum(vapply(runs, function(r) isTRUE(r$reused), logical(1)))
  row$n_failed <- sum(vapply(
    runs, function(r) identical(r$status, "error"), logical(1)
  ))
  row$ended <- Sys.time()
  write_atomic(record$path, function(tmp) saveRDS(row, tmp))
  invisible(row)
}

#' Every record in a grid's directory, in the order the calls started
#'
#' `NULL` when there is none, as in a directory written before bmmtools
#' kept them. A record that cannot be read is left out. Records that do
#' not bind into one table (a column of another type, from a record
#' edited or written by another version) give a warning and `NULL`:
#' called after the cells are fitted, this must not lose their results.
#'
#' @noRd
read_run_records <- function(dir) {
  files <- list.files(file.path(dir, "run-info"), "\\.rds$", full.names = TRUE)
  rows <- lapply(files, function(path) {
    tryCatch(readRDS(path), error = function(e) NULL)
  })
  rows <- Filter(function(x) is.data.frame(x), rows)
  if (length(rows) == 0L) {
    return(NULL)
  }
  out <- tryCatch(dplyr::bind_rows(rows), error = function(e) {
    cli::cli_warn(c(
      "The run records in {.file {file.path(dir, 'run-info')}} cannot be \
       read as one table, so {.code run_info} is left out.",
      x = conditionMessage(e)
    ))
    NULL
  })
  if (is.null(out) || !"started" %in% names(out)) {
    return(out)
  }
  out[order(out$started), , drop = FALSE]
}
