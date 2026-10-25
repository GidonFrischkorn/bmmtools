# run_validation(): execute a validation plan, stage by stage, into one
# directory (Milestone 26, D79-D84).
#
# Every stage caches under dir/<stage>/: its result (result.rds), its
# status (status.rds), and whatever the stage function caches itself
# (the prior fits, the grid's cells). A stage whose stored result was
# made from the same inputs is read, not rerun; a stage that errors is
# recorded as failed and the runner goes on. The result holds every
# stage's object unchanged, so the existing print, summary and plot
# methods apply, plus the numbers, verdicts and headline a card prints.

#' Run a validation plan
#'
#' Executes a [validation_plan()] stage by stage into `dir`, where every
#' stage keeps its result and its status. A rerun reads the stages whose
#' inputs have not changed and runs the rest, so an interrupted run is
#' resumed by calling it again.
#'
#' @param plan A `bmmtools_validation_plan`.
#' @param dir The directory every stage caches under, created if missing.
#' @param previous `NULL`, or an earlier result to compare with: a
#'   `bmmtools_validation`, its `result.rds`, or the directory it was run
#'   in. The comparison is the result's `changes`.
#' @param refit `"on_change"` (the default) reruns a stage when its
#'   inputs differ from those of its stored result, or when it has none;
#'   `"never"` reads what is stored and marks every stage without a stored
#'   result `"skipped"`; `"always"` reruns every stage, with
#'   `refit = "always"` for the prior fits of [prior_check()] and [sbc()].
#'   The recovery grid keeps its cells by their cache keys in every case
#'   (see [recovery_grid()]); delete `<dir>/recovery/cells` to refit them.
#' @param workers Reserved for parallel recovery cells (#19). Only `NULL`,
#'   which runs the cells one after another, is accepted.
#' @param ... Passed to the fitter of every stage, after the plan's
#'   sampler lists, so an argument given here wins (`backend`, `iter`,
#'   `refresh`, ...). The arguments [sbc()] refuses in its `...`, such as
#'   `cores`, do not reach the calibration.
#' @param .fitter The fitting function of every stage, `bmm::bmm()` by
#'   default. Tests inject a stand-in so that nothing is compiled.
#'
#' @return A `bmmtools_validation`: a list with `plan`; `run_info`, the
#'   [run_info()] rows written when the run started and ended, with
#'   `when` and the elapsed `seconds`; `stages`, a tibble with one row per
#'   stage (`stage`, `order`, `title`, `required`, `status`, `issue`,
#'   `reason`, `started`, `finished`, `seconds`, `resumed`, `message`,
#'   `warnings`); one element per stage holding its result object as the
#'   stage function returned it, or `NULL`; `numbers`, the
#'   [card_numbers()] tibble; `verdicts`, the [verdict_table()];
#'   `headline`, a list with `state` (`"validated"`, `"not validated"`,
#'   `"not judged"` or `"incomplete"`), `stages` (the stages the state
#'   names) and `text`; and `changes`, the
#'   [compare_validations()] of `previous` and this run, or `NULL`. The
#'   same object is written to `<dir>/result.rds`.
#'
#' @section Stages:
#' `status` is `"ok"`; `"failed"`, with the condition message in
#' `message`; `"skipped"`, when a stage it reads produced nothing or
#' `refit = "never"` found nothing stored; `"stale"`, when `refit =
#' "never"` read a result made from other inputs or software versions
#' than this plan's (the result is kept, but it is not this plan's); or
#' the plan's own status for a stage that does not run (`"not built"`,
#' `"not applicable"`, `"excluded"`). A failed stage is never an error of
#' `run_validation()`.
#'
#' The headline (D80) is `"incomplete"` while a required stage that
#' exists is not `"ok"`; `"not validated"` when a required criterion
#' failed; `"not judged"` when none failed but one could not be evaluated
#' (a response without a known range, say); and `"validated"`, printed
#' with the plan version and hash, only when every required criterion was
#' evaluated and met.
#'
#' What each stage runs:
#' * `run_record`: [run_info()], written to `<dir>/run-info.rds` before any
#'   other stage and again, with the end row, after the last.
#' * `prior_predictive`: [prior_check()] on data simulated at the plan's
#'   prior-check design point, under bmm's default priors and, when the
#'   calibration has its own, under that one too.
#' * `sbc`: [sbc()] at the plan's design point under the plan's prior.
#' * `recovery`: [recovery_grid()] over the plan's grid at the population,
#'   SD and subject levels, with the model correlations where there are at
#'   least two subject-varying parameters, every cell gated by
#'   [check_convergence()] and contracted against a prior-only fit
#'   (`prior_sd = "fit"`).
#' * `computation`, `identifiability`, `individual_differences`: tables
#'   read off the recovery grid's cells and summary.
#' * `agreement`: [cross_check()] of the fit of the largest design's first
#'   replication against the plan's reference.
#'
#' @section Resuming:
#' A stage's inputs are the plan fields it reads, its seed, the `...` of
#' the call, the fitter, the versions of bmm, brms, bmmtools, the Stan
#' toolchain (and SBC for the calibration) and, for a stage read off the
#' grid, the grid's inputs. They are hashed into the stage's status file,
#' and a stage is read rather than rerun while its hash is unchanged and
#' its status is `"ok"`, so a new bmm release reruns every stage. A grid
#' with failed cells is rerun, which retries those cells and reads the
#' others from their cache, and a stage read off the grid is computed
#' again whenever the grid was not read. The whole
#' plan is written to `<dir>/plan.rds`; a plan whose hash differs from
#' the stored one is announced and replaces it.
#'
#' @seealso [validation_plan()], [card_numbers()], [verdict_table()],
#'   [compare_validations()].
#'
#' @examples
#' \dontrun{
#' plan <- validation_plan("ezdm", seed = 20261024)
#' result <- run_validation(plan, "validation/ezdm")
#' result
#' verdict_table(result)
#' }
#'
#' @export
run_validation <- function(plan,
                           dir,
                           previous = NULL,
                           refit = c("on_change", "never", "always"),
                           workers = NULL,
                           ...,
                           .fitter = NULL) {
  check_validation_plan(plan)
  refit <- rlang::arg_match(refit)
  dir <- check_validation_dir(dir)
  if (!is.null(workers)) {
    cli::cli_abort(c(
      "{.arg workers} is reserved for parallel recovery cells (#19).",
      i = "Leave it {.code NULL}: the cells run one after another."
    ))
  }
  fitter <- .fitter
  if (!is.null(fitter) && !is.function(fitter)) {
    cli::cli_abort(
      "{.arg .fitter} must be a function, not {.obj_type_friendly {fitter}}."
    )
  }
  previous <- read_previous_validation(previous)
  dots <- rlang::list2(...)
  if (length(dots) > 0L && !rlang::is_named(dots)) {
    cli::cli_abort(c(
      "Every argument in {.arg ...} must be named.",
      i = "They are passed on to the fitter by name."
    ))
  }

  started <- Sys.time()
  info <- run_info()
  info$when <- "start"
  info$seconds <- 0
  write_validation_file(file.path(dir, "run-info.rds"), info)
  record_plan(dir, plan, refit)

  table <- validation_stage_table()
  results <- list()
  records <- list()
  for (stage in plan$stages$stage) {
    row <- plan$stages[plan$stages$stage == stage, ]
    reads <- table$reads[table$stage == stage]
    run <- run_validation_stage(
      stage, row, reads, plan, dir, dots, fitter, refit, results, records,
      info
    )
    records[[stage]] <- run$status
    results[stage] <- list(run$result)
  }

  ended <- run_info()
  ended$when <- "end"
  ended$seconds <- as.double(difftime(Sys.time(), started, units = "secs"))
  run_info_rows <- dplyr::bind_rows(info, ended)
  write_validation_file(file.path(dir, "run-info.rds"), run_info_rows)
  records$run_record$finished <- Sys.time()
  records$run_record$seconds <- ended$seconds
  write_stage_status(
    validation_stage_paths(dir, "run_record"), records$run_record
  )

  out <- c(
    list(
      plan = plan,
      run_info = run_info_rows,
      stages = validation_stage_rows(plan, records)
    ),
    results
  )
  out <- structure(out, class = "bmmtools_validation")
  # every stage has run and is cached; an error in what is read off them
  # is reported and the result still written, so it can be inspected
  derived <- tryCatch(
    {
      out$numbers <- card_numbers(out)
      out$verdicts <- verdict_table(out)
      out$headline <- validation_headline(out)
      out$changes <- if (!is.null(previous)) {
        compare_validations(previous, out)
      }
      out
    },
    error = function(e) e
  )
  if (inherits(derived, "error")) {
    message <- conditionMessage(derived)
    cli::cli_warn(c(
      "The numbers and verdicts could not be computed from the stages.",
      x = "{message}"
    ))
    out$headline <- list(
      state = "incomplete", stages = character(),
      text = paste("incomplete (no verdicts:", message, ")")
    )
  } else {
    out <- derived
  }
  write_validation_file(file.path(dir, "result.rds"), out)
  out
}

# the runner ---------------------------------------------------------------

#' @noRd
check_validation_dir <- function(dir, call = rlang::caller_env()) {
  bad <- !is.character(dir) || length(dir) != 1L || is.na(dir) ||
    !nzchar(dir)
  if (bad) {
    cli::cli_abort(
      "{.arg dir} must be a single path, not {.obj_type_friendly {dir}}.",
      call = call
    )
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(dir) || file.access(dir, 2L) != 0L) {
    cli::cli_abort("{.path {dir}} cannot be written to.", call = call)
  }
  normalizePath(dir, mustWork = TRUE)
}

#' Read the `previous` argument (or render_card()'s `result`) into a result
#' @noRd
read_previous_validation <- function(previous, arg = "previous",
                                     call = rlang::caller_env()) {
  if (is.null(previous) || inherits(previous, "bmmtools_validation")) {
    return(previous)
  }
  if (!is.character(previous) || length(previous) != 1L || is.na(previous)) {
    cli::cli_abort(
      "{.arg {arg}} must be a validation result, a path to one, or \\
       {.code NULL}, not {.obj_type_friendly {previous}}.",
      call = call
    )
  }
  path <- if (dir.exists(previous)) {
    file.path(previous, "result.rds")
  } else {
    previous
  }
  if (!file.exists(path)) {
    cli::cli_abort("{.file {path}} does not exist.", call = call)
  }
  out <- readRDS(path)
  if (!inherits(out, "bmmtools_validation")) {
    cli::cli_abort(
      "{.file {path}} holds {.obj_type_friendly {out}}, not a validation \\
       result.",
      call = call
    )
  }
  out
}

#' Write the plan to dir/plan.rds, saying so when it replaces another
#' @noRd
record_plan <- function(dir, plan, refit) {
  path <- file.path(dir, "plan.rds")
  if (file.exists(path)) {
    old <- tryCatch(readRDS(path), error = function(e) NULL)
    old_hash <- if (inherits(old, "bmmtools_validation_plan")) old$hash
    if (!identical(old_hash, plan$hash)) {
      # nolint next: object_usage_linter. Used by cli's glue interpolation.
      what <- switch(refit,
        never = paste(
          "Nothing is rerun: stages made from other inputs are read and",
          "marked stale, and the stored plan is kept."
        ),
        on_change = "Stages whose inputs changed are rerun; the rest are read.",
        always = "Every stage is rerun."
      )
      cli::cli_inform(c(
        "i" = "The plan differs from the one stored in {.path {dir}} \\
               ({substr(old_hash %||% 'none', 1, 8)} -> \\
               {substr(plan$hash, 1, 8)}).",
        " " = "{what}"
      ))
      if (identical(refit, "never")) {
        return(invisible(path))
      }
    }
  }
  write_validation_file(path, plan)
}

#' @noRd
write_validation_file <- function(path, object) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write_atomic(path, function(tmp) saveRDS(object, tmp))
  invisible(path)
}

#' @noRd
validation_stage_paths <- function(dir, stage) {
  stage_dir <- file.path(dir, gsub("_", "-", stage, fixed = TRUE))
  list(
    dir = stage_dir,
    result = file.path(stage_dir, "result.rds"),
    status = file.path(stage_dir, "status.rds")
  )
}

#' @noRd
read_stage_status <- function(paths) {
  if (!file.exists(paths$status)) {
    return(NULL)
  }
  tryCatch(readRDS(paths$status), error = function(e) NULL)
}

#' @noRd
write_stage_status <- function(paths, status) {
  write_validation_file(paths$status, status)
}

#' One stage's status record
#' @noRd
stage_status <- function(stage, status, message = NA_character_,
                         started = as.POSIXct(NA), finished = as.POSIXct(NA),
                         seconds = NA_real_, warnings = character(),
                         input_hash = NA_character_, resumed = FALSE,
                         retry = FALSE) {
  list(
    stage = stage, status = status, message = message, started = started,
    finished = finished, seconds = seconds, warnings = warnings,
    input_hash = input_hash, resumed = resumed, retry = retry
  )
}

#' Whether a stage that ran has parts worth trying again
#'
#' A grid whose cells failed is stored, so its other cells are read on a
#' rerun, but it is not reused as it is: rerunning the stage retries the
#' failed cells and reads the rest from their cache.
#'
#' @noRd
stage_needs_retry <- function(stage, result) {
  if (!identical(stage, "recovery")) {
    return(FALSE)
  }
  cells <- attr(result, "cells")
  !is.null(cells) && any(cells$status != "ok")
}

#' Run, read or skip one stage
#'
#' @return A list with `status` (a [stage_status()] record) and `result`.
#' @noRd
run_validation_stage <- function(stage, row, reads, plan, dir, dots, fitter,
                                 refit, results, records, info) {
  paths <- validation_stage_paths(dir, stage)
  if (!identical(row$status, "planned")) {
    status <- stage_status(stage, row$status, message = row$reason)
    write_stage_status(paths, status)
    return(list(status = status, result = NULL))
  }
  if (identical(stage, "run_record")) {
    status <- stage_status(
      stage, "ok",
      started = Sys.time(), finished = Sys.time(), seconds = 0
    )
    write_validation_file(paths$result, info)
    write_stage_status(paths, status)
    return(list(status = status, result = info))
  }

  key <- stage_input_hash(plan, stage, dots, fitter)
  old <- read_stage_status(paths)
  stored <- !is.null(old) && identical(old$status, "ok") &&
    file.exists(paths$result)
  fresh <- stored && identical(old$input_hash, key) && !isTRUE(old$retry)
  # a view is read only when what it reads was read too: after a rerun of
  # the grid its inputs hash the same, but the grid may have changed
  # (retried cells), and after a failure it reads nothing
  upstream <- if (!is.na(reads)) records[[reads]]
  upstream_read <- is.na(reads) ||
    (identical(upstream$status, "ok") && isTRUE(upstream$resumed))
  reuse <- switch(refit,
    never = stored,
    on_change = fresh && upstream_read,
    always = FALSE
  )
  stored_result <- if (reuse) {
    tryCatch(readRDS(paths$result), error = function(e) NULL)
  }
  if (!is.null(stored_result)) {
    old$resumed <- TRUE
    if (!fresh) {
      # refit = "never" reads a result made from other inputs; it is
      # returned, but it is not this plan's, and the card says so
      old$status <- "stale"
      old$message <- paste(
        "read under `refit = \"never\"`, made from other inputs or",
        "software versions than this plan's"
      )
    }
    cli::cli_inform("Stage {.field {stage}}: read from {.path {paths$dir}}.")
    return(list(status = old, result = stored_result))
  }
  if (identical(refit, "never")) {
    status <- stage_status(
      stage, "skipped",
      message = "`refit = \"never\"` and no stored result",
      input_hash = key
    )
    write_stage_status(paths, status)
    return(list(status = status, result = NULL))
  }
  if (!is.na(reads) && is.null(results[[reads]])) {
    status <- stage_status(
      stage, "skipped",
      message = sprintf("needs the %s stage, which produced no result", reads),
      input_hash = key
    )
    write_stage_status(paths, status)
    return(list(status = status, result = NULL))
  }

  cli::cli_inform("Stage {.field {stage}}: running.")
  dir.create(paths$dir, recursive = TRUE, showWarnings = FALSE)
  # an older result never survives next to a newer status: a failure, or
  # an interruption before the status is written, leaves no result
  unlink(paths$result)
  started <- Sys.time()
  warnings <- character()
  result <- tryCatch(
    {
      out <- withCallingHandlers(
        run_stage_body(stage, plan, paths, dots, fitter, refit, results),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
        }
      )
      write_validation_file(paths$result, out)
      out
    },
    error = function(e) e
  )
  finished <- Sys.time()
  seconds <- as.double(difftime(finished, started, units = "secs"))
  if (inherits(result, "error")) {
    message <- conditionMessage(result)
    cli::cli_warn(c(
      "Stage {.field {stage}} failed; the other stages go on.",
      x = "{message}"
    ))
    status <- stage_status(
      stage, "failed",
      message = message, started = started, finished = finished,
      seconds = seconds, warnings = warnings, input_hash = key
    )
    write_stage_status(paths, status)
    return(list(status = status, result = NULL))
  }
  status <- stage_status(
    stage, "ok",
    started = started, finished = finished, seconds = seconds,
    warnings = warnings, input_hash = key,
    retry = stage_needs_retry(stage, result)
  )
  write_stage_status(paths, status)
  list(status = status, result = result)
}

#' The fields of a plan each stage reads
#' @noRd
stage_input_fields <- function(stage) {
  truth <- c(
    "model", "pars", "sds", "cors", "formula", "re_cor", "trial_design",
    "generator"
  )
  switch(stage,
    prior_predictive = c(truth, "prior_check"),
    sbc = c(truth, "sbc"),
    recovery = c(
      truth, "grid", "reps", "sampler", "convergence", "correlations"
    ),
    agreement = "reference",
    character()
  )
}

#' The hash of what a stage's result was made from
#' @noRd
stage_input_hash <- function(plan, stage, dots, fitter = NULL) {
  table <- validation_stage_table()
  reads <- table$reads[table$stage == stage]
  fields <- stage_input_fields(stage)
  seed_name <- switch(stage,
    prior_predictive = "prior_predictive",
    sbc = "sbc",
    recovery = "recovery",
    NULL
  )
  parts <- list(
    stage = stage,
    fields = plan_key(unclass(plan)[fields]),
    seed = if (!is.null(seed_name)) plan$seeds[[seed_name]],
    dots = plan_key(dots),
    # the prior check reads the calibration's prior, nothing else of it
    sbc_prior = if (identical(stage, "prior_predictive")) {
      plan_key(plan$sbc$prior)
    },
    # a new bmm, brms, Stan toolchain or bmmtools reruns the stage: the
    # card is regenerated on each release (D81), and a stage read from
    # its result file never reaches the fit caches that key on them
    software = stage_software(stage, plan, dots),
    fitter = if (is.null(fitter)) "bmm::bmm" else function_key(fitter),
    upstream = if (!is.na(reads)) stage_input_hash(plan, reads, dots, fitter)
  )
  rlang::hash(parts)
}

#' The software versions a stage's result depends on
#' @noRd
stage_software <- function(stage, plan, dots) {
  sampler <- switch(stage,
    prior_predictive = plan$prior_check$sampler,
    sbc = plan$sbc$sampler,
    plan$sampler
  )
  backend <- merge_sampler(sampler, dots)$backend
  c(
    bmm = package_version_string("bmm"),
    brms = package_version_string("brms"),
    bmmtools = package_version_string("bmmtools"),
    SBC = if (identical(stage, "sbc")) package_version_string("SBC"),
    toolchain = backend_version_string(backend)
  )
}

#' The stage-specific work
#' @noRd
run_stage_body <- function(stage, plan, paths, dots, fitter, refit, results) {
  switch(stage,
    prior_predictive = stage_prior_predictive(plan, paths, dots, fitter, refit),
    sbc = stage_sbc(plan, paths, dots, fitter, refit),
    recovery = stage_recovery(plan, paths, dots, fitter),
    computation = stage_computation(plan, results$recovery),
    identifiability = stage_identifiability(plan, results$recovery),
    individual_differences = stage_individual_differences(
      plan, results$recovery
    ),
    agreement = stage_agreement(plan, results$recovery),
    cli::cli_abort("No stage {.val {stage}} can be run.")
  )
}

#' The formula every stage fits
#' @noRd
validation_formula <- function(plan) {
  plan$formula %||% recovery_formula(plan$model, re_cor = plan$re_cor)
}

#' One data set at a design point, simulated from the plan's truth
#'
#' The prior check and the calibration fit a data set of the model's own
#' layout; its responses do not enter a prior-only fit, and SBC replaces
#' them with its own simulations.
#'
#' @noRd
validation_data <- function(plan, n_subjects, n_trials, seed) {
  args <- list(
    model = plan$model, pars = plan$pars, n_subjects = n_subjects,
    n_trials = n_trials, sds = plan$sds, cors = plan$cors,
    trial_design = plan$trial_design, generator = plan$generator,
    seed = seed
  )
  if (adapter_reads_formula(plan$model)) args$formula <- plan$formula
  rlang::exec(simulate_recovery, !!!args)$data
}

#' A sampler list with the call's arguments laid over it
#' @noRd
merge_sampler <- function(sampler, dots) {
  out <- sampler %||% list()
  for (name in names(dots)) out[name] <- list(dots[[name]])
  out
}

#' A prior written as a data frame, turned into a brmsprior
#' @noRd
as_brms_prior <- function(prior) {
  if (is.null(prior) || inherits(prior, "brmsprior")) {
    return(prior)
  }
  rlang::check_installed("brms", "to build the plan's prior.")
  columns <- intersect(
    names(prior),
    c("prior", "class", "coef", "group", "nlpar", "dpar", "lb", "ub")
  )
  rows <- lapply(seq_len(nrow(prior)), function(i) {
    args <- lapply(prior[i, columns, drop = FALSE], function(v) {
      v <- as.character(v)
      if (is.na(v)) "" else v
    })
    rlang::exec(brms::set_prior, !!!args)
  })
  Reduce(`+`, rows)
}

#' @noRd
stage_prior_predictive <- function(plan, paths, dots, fitter, refit) {
  pc <- plan$prior_check
  data <- validation_data(
    plan, pc$n_subjects, pc$n_trials, plan$seeds[["prior_predictive"]]
  )
  sbc_prior <- as_brms_prior(plan$sbc$prior)
  prior <- if (!is.null(sbc_prior)) list(default = NULL, sbc = sbc_prior)
  args <- c(
    list(
      model = plan$model, formula = validation_formula(plan), data = data,
      prior = prior, summary = pc$summary, n_draws = pc$n_draws,
      range = pc$range, file = file.path(paths$dir, "fit"),
      refit = if (identical(refit, "always")) "always" else "on_change",
      seed = plan$seeds[["prior_predictive"]]
    ),
    merge_sampler(pc$sampler, dots),
    list(.fitter = fitter)
  )
  rlang::exec(prior_check, !!!args)
}

#' @noRd
stage_sbc <- function(plan, paths, dots, fitter, refit) {
  s <- plan$sbc
  data <- validation_data(plan, s$n_subjects, s$n_trials, plan$seeds[["sbc"]])
  sampler <- merge_sampler(s$sampler, dots)
  sampler <- sampler[setdiff(names(sampler), names(sbc_refused_dots()))]
  args <- c(
    list(
      model = plan$model, formula = validation_formula(plan), data = data,
      prior = as_brms_prior(s$prior), n_sims = s$n_sims,
      level = s$level %||% c("population", "sd"), generator = s$generator
    ),
    sampler,
    list(
      seed = plan$seeds[["sbc"]], file = file.path(paths$dir, "prior-fit"),
      refit = if (identical(refit, "always")) "always" else "on_change",
      cores_per_fit = s$cores_per_fit, thin_ranks = s$thin_ranks,
      .fitter = fitter
    )
  )
  rlang::exec(sbc, !!!args)
}

#' @noRd
stage_recovery <- function(plan, paths, dots, fitter) {
  args <- c(
    list(
      model = plan$model, grid = plan$grid, pars = plan$pars,
      dir = file.path(paths$dir, "cells"), reps = plan$reps,
      sds = plan$sds, cors = plan$cors, formula = plan$formula,
      trial_design = plan$trial_design, generator = plan$generator,
      seed = plan$seeds[["recovery"]], re_cor = plan$re_cor,
      levels = c("population", "sd", "subject"),
      correlations = plan$correlations, convergence = plan$convergence,
      prior_sd = "fit"
    ),
    merge_sampler(plan$sampler, dots),
    list(.fitter = fitter)
  )
  rlang::exec(recovery_grid, !!!args)
}

#' The design label of each grid row, keyed by the grid's condition
#' @noRd
design_labels <- function(grid) {
  stats::setNames(
    sprintf("%d x %d", grid$n_subjects, grid$n_trials),
    sprintf("row-%d", seq_len(nrow(grid)))
  )
}

#' The condition of the largest design: most observations, then most
#' subjects
#' @noRd
largest_condition <- function(grid) {
  size <- grid$n_subjects * grid$n_trials
  i <- order(-size, -grid$n_subjects)[[1L]]
  sprintf("row-%d", i)
}

#' @noRd
stage_computation <- function(plan, recovery) {
  cells <- attr(recovery, "cells")
  labels <- design_labels(plan$grid)
  one <- function(rows, condition) {
    converged <- rows$converged %in% TRUE
    tibble::tibble(
      condition = condition,
      n_cells = nrow(rows),
      n_ok = sum(rows$status == "ok"),
      n_converged = sum(converged),
      pass_rate = if (nrow(rows) > 0L) mean(converged) else NA_real_,
      max_rhat = safe_stat(rows$max_rhat, max),
      min_ess_bulk = safe_stat(rows$min_ess_bulk, min),
      min_ess_tail = safe_stat(rows$min_ess_tail, min),
      n_divergent = safe_stat(rows$n_divergent, sum),
      median_fit_seconds = safe_stat(rows$fit_seconds, stats::median)
    )
  }
  by_row <- lapply(names(labels), function(condition) {
    one(cells[cells$condition == condition, ], labels[[condition]])
  })
  dplyr::bind_rows(c(by_row, list(one(cells, "all"))))
}

#' A statistic over the non-missing values, `NA` when there are none
#' @noRd
safe_stat <- function(x, f) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA_real_ else as.double(f(x))
}

#' The recovery summary with a design label and the largest design marked
#' @noRd
labelled_summary <- function(plan, recovery) {
  out <- tibble::as_tibble(unclass(summary(recovery)))
  labels <- design_labels(plan$grid)
  out$largest <- out$condition == largest_condition(plan$grid)
  out$condition <- unname(labels[out$condition])
  out
}

#' @noRd
stage_identifiability <- function(plan, recovery) {
  rows <- labelled_summary(plan, recovery)
  rows <- rows[rows$level %in% c("population", "sd"), ]
  rows[c(
    "condition", "largest", "term", "level", "n_replications", "z_mean",
    "z_sd", "contraction", "z_mean_mcse", "z_sd_mcse"
  )]
}

#' @noRd
stage_individual_differences <- function(plan, recovery) {
  cors <- attr(recovery, "correlations")
  if (is.null(cors) || nrow(cors) == 0L) {
    return(tibble::as_tibble(empty_cor_recovery_summary()))
  }
  out <- tibble::as_tibble(unclass(summary(cors)))
  if ("condition" %in% names(out)) {
    out$condition <- unname(design_labels(plan$grid)[out$condition])
  }
  out
}

#' @noRd
stage_agreement <- function(plan, recovery) {
  cells <- attr(recovery, "cells")
  condition <- largest_condition(plan$grid)
  candidates <- cells[cells$condition == condition & cells$status == "ok", ]
  if (nrow(candidates) == 0L) {
    cli::cli_abort(
      "The largest design ({design_labels(plan$grid)[[condition]]}) has no \\
       fitted cell to compare."
    )
  }
  cell <- candidates[order(candidates$replication), ][1L, ]
  if (!file.exists(cell$file)) {
    cli::cli_abort(c(
      "The fit of the largest design's replication {cell$replication} is \\
       gone: {.file {cell$file}}.",
      i = "Agreement compares a fit; refit the cell or keep its file."
    ))
  }
  fit <- readRDS(cell$file)
  columns <- intersect(
    c("term", "estimate", "ci_low", "ci_high", "source"),
    names(plan$reference)
  )
  cross_check(fit, plan$reference[columns], level = "population")
}

#' The stage table of a result: the plan's, with what this run did
#' @noRd
validation_stage_rows <- function(plan, records) {
  stages <- plan$stages
  get <- function(field, empty) {
    lapply(stages$stage, function(s) records[[s]][[field]] %||% empty)
  }
  stages$status <- vapply(
    seq_len(nrow(stages)),
    function(i) records[[stages$stage[[i]]]]$status %||% stages$status[[i]],
    character(1)
  )
  stages$started <- do.call(c, get("started", as.POSIXct(NA)))
  stages$finished <- do.call(c, get("finished", as.POSIXct(NA)))
  stages$seconds <- vapply(get("seconds", NA_real_), as.double, double(1))
  stages$resumed <- vapply(get("resumed", FALSE), isTRUE, logical(1))
  message <- vapply(
    get("message", NA_character_),
    function(m) if (length(m) == 0L) NA_character_ else as.character(m[[1L]]),
    character(1)
  )
  stages$message <- ifelse(is.na(message), stages$reason, message)
  stages$warnings <- get("warnings", character())
  stages
}

# the numbers a card prints ------------------------------------------------

#' The numbers of a validation, one per row
#'
#' Every number a validation card prints, in one long tibble, so a card
#' can be checked against it and two runs compared with
#' [compare_validations()].
#'
#' @param result A `bmmtools_validation`.
#'
#' @return A tibble with the columns `stage`, `condition` (the design as
#'   `"<subjects> x <trials>"`, the prior set, `"all"`, or `NA`), `level`
#'   (the estimate level or the correlation estimator, or `NA`), `term`,
#'   `metric`, `value` and `mcse`, the Monte Carlo standard error where
#'   one follows from the stored results: binomial for a rate, at the rate
#'   `(x + 1) / (n + 2)` so that a rate of 0 or 1 keeps an error above 0,
#'   `sqrt((rmse^2 - bias^2) / (n - 1))` for a bias, and the summary's own
#'   for a mean z-score and its SD. `NA` elsewhere.
#'
#' @seealso [run_validation()], which stores it as `numbers`.
#' @export
card_numbers <- function(result) {
  check_validation_result(result)
  pieces <- list(
    numbers_stage_seconds(result$stages),
    numbers_prior_predictive(result$prior_predictive, result$plan),
    numbers_sbc(result$sbc, result$plan),
    numbers_recovery(result$recovery, result$plan),
    numbers_computation(result$computation),
    numbers_identifiability(result$identifiability),
    numbers_individual_differences(result$individual_differences),
    numbers_agreement(result$agreement)
  )
  out <- dplyr::bind_rows(c(list(empty_card_numbers()), pieces))
  out
}

#' @noRd
empty_card_numbers <- function() {
  tibble::tibble(
    stage = character(), condition = character(), level = character(),
    term = character(), metric = character(), value = double(),
    mcse = double()
  )
}

#' Rows of the numbers tibble from wide columns
#' @noRd
long_numbers <- function(stage, rows, metrics, mcse = list(),
                         condition = NULL, level = NULL, term = NULL) {
  if (is.null(rows) || nrow(rows) == 0L) {
    return(NULL)
  }
  pick <- function(column) {
    if (is.null(column) || !column %in% names(rows)) {
      return(rep(NA_character_, nrow(rows)))
    }
    as.character(rows[[column]])
  }
  metrics <- intersect(metrics, names(rows))
  keys <- list(
    condition = pick(condition), level = pick(level), term = pick(term)
  )
  dplyr::bind_rows(lapply(metrics, function(m) {
    # computed outside tibble(), whose columns would mask the arguments
    values <- as.double(rows[[m]])
    se <- mcse[[m]]
    se <- if (is.null(se)) rep(NA_real_, nrow(rows)) else as.double(se)
    tibble::tibble(
      stage = rep(stage, nrow(rows)),
      condition = keys$condition,
      level = keys$level,
      term = keys$term,
      metric = rep(m, nrow(rows)),
      value = values,
      mcse = se
    )
  }))
}

#' Binomial MCSE of a rate
#' @noRd
rate_mcse <- function(p, n) {
  p[!is.na(p) & (p < 0 | p > 1)] <- NA_real_
  # the rate of (x + 1) / (n + 2) successes, so that a rate of 0 or 1,
  # common at ten replications, keeps an error above zero
  p <- (p * n + 1) / (n + 2)
  out <- sqrt(p * (1 - p) / n)
  out[!is.finite(out)] <- NA_real_
  out
}

#' MCSE of a bias from the RMSE and the bias over `n` rows
#' @noRd
bias_mcse <- function(bias, rmse, n) {
  out <- sqrt(pmax(rmse^2 - bias^2, 0) / (n - 1))
  out[!is.finite(out) | n <= 1] <- NA_real_
  out
}

#' @noRd
numbers_stage_seconds <- function(stages) {
  rows <- stages[!is.na(stages$seconds), ]
  if (nrow(rows) == 0L) {
    return(NULL)
  }
  tibble::tibble(
    stage = rows$stage, condition = NA_character_, level = NA_character_,
    term = NA_character_, metric = "seconds", value = rows$seconds,
    mcse = NA_real_
  )
}

#' @noRd
numbers_prior_predictive <- function(x, plan) {
  if (is.null(x)) {
    return(NULL)
  }
  rows <- tibble::as_tibble(unclass(x))
  rate <- rows$statistic %in% c("floor_rate", "ceiling_rate")
  mcse <- rep(NA_real_, nrow(rows))
  mcse[rate] <- rate_mcse(rows$value[rate], plan$prior_check$n_draws)
  tibble::tibble(
    stage = "prior_predictive",
    condition = as.character(rows$prior),
    level = NA_character_,
    term = as.character(rows$response),
    metric = as.character(rows$statistic),
    value = as.double(rows$value),
    mcse = mcse
  )
}

#' @noRd
numbers_sbc <- function(x, plan) {
  if (is.null(x)) {
    return(NULL)
  }
  meta <- attr(x, "bmmtools_sbc")
  diagnostics <- meta$diagnostics
  counts <- c("n_fits", "n_high_rhat", "n_low_ess_to_rank", "n_missing")
  out <- tibble::tibble(
    stage = "sbc", condition = NA_character_, level = NA_character_,
    term = NA_character_, metric = counts,
    value = vapply(counts, function(k) {
      as.double(diagnostics[[k]] %||% NA_real_)
    }, double(1), USE.NAMES = FALSE),
    mcse = NA_real_
  )
  band <- sbc_ecdf_band(x, plan$criteria$sbc_ecdf_prob)
  if (is.data.frame(band)) {
    out <- dplyr::bind_rows(out, tibble::tibble(
      stage = "sbc", condition = NA_character_, level = NA_character_,
      term = band$variable, metric = "outside_band",
      value = as.double(band$outside), mcse = NA_real_
    ))
  }
  out
}

#' Which ranked variables leave SBC's simultaneous ECDF band
#'
#' The band is the one `SBC::plot_ecdf_diff()` draws, read from
#' `SBC::data_for_ecdf_plots()`: a variable whose ECDF leaves it at any
#' grid point is rejected as non-uniform at `prob`.
#'
#' @return A tibble with `variable` and `outside`, or `NULL` when SBC
#'   cannot compute the band, with the reason as its `"reason"` attribute.
#' @noRd
sbc_ecdf_band <- function(results, prob = 0.95) {
  if (!rlang::is_installed("SBC")) {
    return(structure(list(), reason = "SBC is not installed"))
  }
  data <- tryCatch(
    SBC::data_for_ecdf_plots(results, prob = prob),
    error = function(e) structure(list(), reason = conditionMessage(e))
  )
  if (!is.null(attr(data, "reason"))) {
    return(data)
  }
  limits <- data$limits_df
  # the limits are a step function: the row of each grid point is the one
  # at which the uniform value has reached it
  limits <- limits[abs(limits$x - limits$uniform_val) < 1e-12, ]
  ecdf <- data$ecdf_df
  at <- match(round(ecdf$z, 10), round(limits$x, 10))
  outside <- !is.na(at) & (
    ecdf$ecdf < limits$lower[at] - 1e-12 | ecdf$ecdf > limits$upper[at] + 1e-12
  )
  variables <- unique(as.character(ecdf$variable))
  tibble::tibble(
    variable = variables,
    outside = vapply(variables, function(v) {
      any(outside[ecdf$variable == v])
    }, logical(1), USE.NAMES = FALSE)
  )
}

#' @noRd
numbers_recovery <- function(x, plan) {
  if (is.null(x)) {
    return(NULL)
  }
  rows <- labelled_summary(plan, x)
  levels <- rows$level %in% c("population", "sd")
  population <- rows[levels, ]
  subject <- rows[rows$level == "subject", ]
  dplyr::bind_rows(
    long_numbers(
      "recovery", population,
      c(
        "n_replications", "bias", "rmse", "coverage", "coverage_50",
        "ci_width"
      ),
      mcse = list(
        bias = bias_mcse(population$bias, population$rmse, population$n),
        coverage = rate_mcse(population$coverage, population$n_replications),
        coverage_50 = rate_mcse(
          population$coverage_50, population$n_replications
        )
      ),
      condition = "condition", level = "level", term = "term"
    ),
    long_numbers(
      "recovery", subject,
      c("n_replications", "r", "r_low", "r_high", "rmse", "coverage"),
      mcse = list(
        coverage = rate_mcse(subject$coverage, subject$n_replications)
      ),
      condition = "condition", level = "level", term = "term"
    )
  )
}

#' @noRd
numbers_computation <- function(x) {
  long_numbers(
    "computation", x,
    c(
      "n_cells", "n_converged", "pass_rate", "max_rhat", "min_ess_bulk",
      "min_ess_tail", "n_divergent", "median_fit_seconds"
    ),
    mcse = if (!is.null(x)) list(pass_rate = rate_mcse(x$pass_rate, x$n_cells)),
    condition = "condition"
  )
}

#' @noRd
numbers_identifiability <- function(x) {
  long_numbers(
    "identifiability", x, c("z_mean", "z_sd", "contraction"),
    mcse = if (!is.null(x)) list(z_mean = x$z_mean_mcse, z_sd = x$z_sd_mcse),
    condition = "condition", level = "level", term = "term"
  )
}

#' @noRd
numbers_individual_differences <- function(x) {
  if (is.null(x) || nrow(x) == 0L) {
    return(NULL)
  }
  long_numbers(
    "individual_differences", x,
    c(
      "n_replications", "true_value", "mean_estimate", "bias", "rmse",
      "coverage", "ci_width", "z_mean"
    ),
    mcse = list(
      bias = bias_mcse(x$bias, x$rmse, x$n_replications),
      coverage = rate_mcse(x$coverage, x$n_replications),
      z_mean = x$z_mean_mcse
    ),
    condition = "condition", level = "estimator", term = "term"
  )
}

#' @noRd
numbers_agreement <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  rows <- tibble::as_tibble(unclass(x))
  rows$covered <- as.double(rows$covered)
  long_numbers(
    "agreement", rows, c("estimate", "reference", "bias", "covered"),
    level = "level", term = "term"
  )
}

# verdicts and the headline ---------------------------------------------

#' The verdict of every stage against the plan's criteria
#'
#' One row per criterion, with the bar written out, so a card can print
#' the criterion next to its verdict rather than a bare pass or fail.
#'
#' @param result A `bmmtools_validation`.
#'
#' @return A tibble with `stage`, `title`, `required`, `status` (the
#'   stage's), `criterion` (the bar in words), `value` (the number judged),
#'   `pass` (`TRUE`, `FALSE`, or `NA` when the stage did not run or the
#'   criterion cannot be evaluated) and `note`.
#'
#' @details
#' The criteria are the plan's (see [validation_plan()], "Criteria"):
#' * prior predictive: the largest floor or ceiling rate under the
#'   calibration's prior (or bmm's defaults when the calibration uses
#'   them) at most `prior_bound_rate_max`; not judged for a response
#'   without a known range.
#' * SBC: no ranked variable outside the `sbc_ecdf_prob` simultaneous
#'   ECDF band, and no fit with a high Rhat, a low rank ESS or missing
#'   diagnostics.
#' * computation: the share of converged cells at least
#'   `convergence_rate_min`.
#' * recovery: the coverage of every population parameter, pooled over the
#'   designs, not below `coverage_level` by a one-sided binomial test at
#'   `coverage_alpha` divided by the number of parameters. The coverage of
#'   each design is reported, not judged.
#' * identifiability: at the largest design, every population parameter's
#'   mean contraction above `contraction_min` and its mean z-score within
#'   `z_mean_max` of 0.
#' * individual differences: the coverage of every generating correlation
#'   held to the same pooled test, divided by the number of correlations.
#' * agreement: reported, not judged.
#'
#' The ECDF band is simultaneous over the ranks of one variable; it is not
#' corrected for the number of variables, so with more ranked variables a
#' calibrated model leaves some band more often.
#'
#' @export
verdict_table <- function(result) {
  check_validation_result(result)
  stages <- result$stages
  criteria <- result$plan$criteria
  rows <- lapply(seq_len(nrow(stages)), function(i) {
    stage <- stages$stage[[i]]
    status <- stages$status[[i]]
    judged <- if (identical(status, "ok")) {
      stage_verdicts(stage, result, criteria)
    }
    if (is.null(judged)) {
      if (identical(stage, "run_record") && identical(status, "ok")) {
        return(NULL)
      }
      judged <- verdict_row(
        stage_criterion_text(stage, result, criteria), NA_real_, NA,
        if (identical(status, "ok")) NA_character_ else stages$message[[i]]
      )
    }
    judged$stage <- stage
    judged$title <- stages$title[[i]]
    judged$required <- stages$required[[i]]
    judged$status <- status
    judged
  })
  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0L) {
    return(tibble::tibble(
      stage = character(), title = character(), required = logical(),
      status = character(), criterion = character(), value = double(),
      pass = logical(), note = character()
    ))
  }
  out[c(
    "stage", "title", "required", "status", "criterion", "value", "pass",
    "note"
  )]
}

#' @noRd
verdict_row <- function(criterion, value, pass, note = NA_character_) {
  tibble::tibble(
    criterion = criterion, value = as.double(value), pass = as.logical(pass),
    note = as.character(note)
  )
}

#' @noRd
fmt <- function(x) format(signif(x, 3), trim = TRUE)

#' The pooled coverage test, in words
#' @noRd
pooled_test_text <- function(criteria, of = "parameters") {
  sprintf(
    paste(
      "not below %s by a one-sided binomial test at %s divided by the",
      "number of %s"
    ),
    fmt(criteria$coverage_level), fmt(criteria$coverage_alpha), of
  )
}

#' The criterion of a stage that did not run, in words
#' @noRd
stage_criterion_text <- function(stage, result, criteria) {
  level <- fmt(100 * criteria$coverage_level)
  switch(stage,
    prior_predictive = sprintf(
      "largest floor or ceiling rate of the prior predictive draws at most %s",
      fmt(criteria$prior_bound_rate_max)
    ),
    sbc = sprintf(
      "no ranked variable outside SBC's %s%% simultaneous ECDF band",
      fmt(100 * criteria$sbc_ecdf_prob)
    ),
    recovery = sprintf(
      paste(
        "%s%% interval coverage of every population parameter, pooled",
        "across designs, %s"
      ),
      level, pooled_test_text(criteria)
    ),
    computation = sprintf(
      "at least %s of the recovery cells pass check_convergence()",
      fmt(criteria$convergence_rate_min)
    ),
    identifiability = sprintf(
      paste(
        "contraction above %s and mean z-score within +/-%s for every",
        "population parameter at the largest design"
      ),
      fmt(criteria$contraction_min), fmt(criteria$z_mean_max)
    ),
    individual_differences = sprintf(
      paste(
        "%s%% interval coverage of every generating correlation, pooled",
        "across designs, %s"
      ),
      level, pooled_test_text(criteria, "correlations")
    ),
    data_requirements = paste(
      "the smallest design with person-level r_low >= 0.80 at 50 subjects",
      "(D83)"
    ),
    agreement = "reported, not judged",
    model_recovery = "not defined yet",
    posterior_checks = "not defined yet",
    run_record = "the run is recorded",
    NA_character_
  )
}

#' The verdict rows of a stage that ran
#' @noRd
stage_verdicts <- function(stage, result, criteria) {
  switch(stage,
    prior_predictive = verdict_prior_predictive(result, criteria),
    sbc = verdict_sbc(result, criteria),
    computation = verdict_computation(result, criteria),
    recovery = verdict_recovery(result, criteria),
    identifiability = verdict_identifiability(result, criteria),
    individual_differences = verdict_individual(result, criteria),
    agreement = verdict_agreement(result),
    NULL
  )
}

#' @noRd
verdict_prior_predictive <- function(result, criteria) {
  x <- tibble::as_tibble(unclass(result$prior_predictive))
  set <- if ("sbc" %in% x$prior) "sbc" else "default"
  label <- if (identical(set, "sbc")) {
    "the calibration's prior"
  } else {
    "bmm's default priors"
  }
  bounds <- x$statistic %in% c("floor_rate", "ceiling_rate")
  rates <- x$value[x$prior == set & bounds]
  text <- sprintf(
    paste(
      "largest floor or ceiling rate of the prior predictive draws under %s",
      "at most %s"
    ),
    label, fmt(criteria$prior_bound_rate_max)
  )
  if (length(rates) == 0L || all(is.na(rates))) {
    return(verdict_row(
      text, NA_real_, NA,
      "not judged: no floor or ceiling is known for this response"
    ))
  }
  value <- max(rates, na.rm = TRUE)
  verdict_row(text, value, value <= criteria$prior_bound_rate_max)
}

#' @noRd
verdict_sbc <- function(result, criteria) {
  meta <- attr(result$sbc, "bmmtools_sbc")
  d <- meta$diagnostics
  counts <- c(d$n_high_rhat, d$n_low_ess_to_rank, d$n_missing)
  # without the record bmmtools writes, nothing was checked: NA, not 0
  bad <- if (length(counts) == 3L) sum(counts) else NA_real_
  d$rhat_max <- d$rhat_max %||% NA_real_
  d$ess_to_rank_min <- d$ess_to_rank_min %||% NA_real_
  band <- sbc_ecdf_band(result$sbc, criteria$sbc_ecdf_prob)
  band_text <- sprintf(
    "no ranked variable outside SBC's %s%% simultaneous ECDF band",
    fmt(100 * criteria$sbc_ecdf_prob)
  )
  band_row <- if (!is.data.frame(band)) {
    verdict_row(
      band_text, NA_real_, NA,
      paste("not judged: no ECDF band:", attr(band, "reason"))
    )
  } else {
    outside <- band$variable[band$outside]
    verdict_row(
      band_text, length(outside), length(outside) == 0L,
      if (length(outside) > 0L) paste(outside, collapse = ", ") else NA
    )
  }
  dplyr::bind_rows(
    band_row,
    verdict_row(
      sprintf(
        paste(
          "no fit with Rhat above %s, rank ESS ratio below %s, or missing",
          "diagnostics"
        ),
        fmt(d$rhat_max), fmt(d$ess_to_rank_min)
      ),
      bad, if (is.na(bad)) NA else bad == 0L,
      if (is.na(bad)) "not judged: no SBC diagnostics" else NA
    )
  )
}

#' @noRd
verdict_computation <- function(result, criteria) {
  x <- result$computation
  value <- x$pass_rate[x$condition == "all"]
  verdict_row(
    stage_criterion_text("computation", result, criteria),
    value, value >= criteria$convergence_rate_min
  )
}

#' Coverage pooled across designs, judged by a binomial test per term
#'
#' A term fails when `pbinom(covered, n, level)` is below `coverage_alpha`
#' over the number of terms. The per-design bar of 26.1 failed a
#' calibrated model's stage two times in three at 10 replications
#' (ARCHITECTURE § 26.0, amendment 3). `label` names the design column,
#' which pooling no longer needs; it is kept for the callers.
#' @noRd
coverage_verdict <- function(text, rows, label, criteria, empty_note) {
  if (is.null(rows) || nrow(rows) == 0L || all(is.na(rows$coverage))) {
    return(verdict_row(text, NA_real_, NA, empty_note))
  }
  rows <- rows[!is.na(rows$coverage), ]
  # coverage leaves out a replication without an interval while
  # n_replications counts it; such a replication is a failed fit, which
  # the runner retries and the computation stage judges
  covered <- round(rows$coverage * rows$n_replications)
  terms <- unique(rows$term)
  k <- tapply(covered, factor(rows$term, terms), sum)
  n <- tapply(rows$n_replications, factor(rows$term, terms), sum)
  p <- stats::pbinom(k, n, criteria$coverage_level)
  low <- p < criteria$coverage_alpha / length(terms)
  note <- if (any(low)) {
    paste(
      sprintf(
        "%s: %d of %d covered, p = %s", terms[low], k[low], n[low],
        fmt(p[low])
      ),
      collapse = "; "
    )
  } else {
    NA_character_
  }
  verdict_row(text, unname(min(k / n)), !any(low), note)
}

#' @noRd
verdict_recovery <- function(result, criteria) {
  rows <- labelled_summary(result$plan, result$recovery)
  rows <- rows[rows$level == "population", ]
  coverage_verdict(
    stage_criterion_text("recovery", result, criteria), rows, "condition",
    criteria, "not judged: no population rows were scored"
  )
}

#' @noRd
verdict_identifiability <- function(result, criteria) {
  x <- result$identifiability
  x <- x[x$largest & x$level == "population", ]
  design <- if (nrow(x) > 0L) x$condition[[1L]] else "none"
  contraction_text <- sprintf(
    paste(
      "contraction above %s for every population parameter at the largest",
      "design (%s)"
    ),
    fmt(criteria$contraction_min), design
  )
  z_text <- sprintf(
    paste(
      "mean z-score within +/-%s for every population parameter at the",
      "largest design (%s)"
    ),
    fmt(criteria$z_mean_max), design
  )
  contraction <- if (nrow(x) == 0L || all(is.na(x$contraction))) {
    verdict_row(contraction_text, NA_real_, NA, "not judged: no prior SD")
  } else {
    value <- min(x$contraction, na.rm = TRUE)
    verdict_row(contraction_text, value, value > criteria$contraction_min)
  }
  z <- if (nrow(x) == 0L || all(is.na(x$z_mean))) {
    verdict_row(z_text, NA_real_, NA, "not judged: no z-scores")
  } else {
    value <- max(abs(x$z_mean), na.rm = TRUE)
    verdict_row(z_text, value, value <= criteria$z_mean_max)
  }
  dplyr::bind_rows(contraction, z)
}

#' @noRd
verdict_individual <- function(result, criteria) {
  x <- result$individual_differences
  if (!is.null(x) && nrow(x) > 0L) x <- x[x$estimator == "model", ]
  coverage_verdict(
    stage_criterion_text("individual_differences", result, criteria), x,
    if (!is.null(x) && "condition" %in% names(x)) "condition" else "term",
    criteria, "not judged: no correlation was scored"
  )
}

#' @noRd
verdict_agreement <- function(result) {
  x <- result$agreement
  value <- if (nrow(x) > 0L) max(abs(x$bias), na.rm = TRUE) else NA_real_
  verdict_row(
    "reported, not judged: largest absolute difference from the reference",
    value, NA
  )
}

#' Complete, validated, or neither
#'
#' D80: complete when every required stage that exists ran to a stored
#' result; validated at plan V when, in addition, no required criterion
#' failed. A criterion that cannot be evaluated is named, not counted.
#'
#' @noRd
validation_headline <- function(result) {
  stages <- result$stages
  plan <- result$plan
  required <- stages[
    stages$required & !stages$status %in% c("not built", "not applicable"),
  ]
  missing <- required$stage[required$status != "ok"]
  titles <- function(s) stages$title[match(s, stages$stage)]
  if (length(missing) > 0L) {
    return(list(
      state = "incomplete", stages = missing,
      text = sprintf(
        "incomplete (stages %s)", paste(titles(missing), collapse = ", ")
      )
    ))
  }
  # recomputed, never read from the result, so an edited result cannot
  # carry a stale verdict into its headline; a stage that did not run is
  # judged by completeness above, not here
  verdicts <- verdict_table(result)
  verdicts <- verdicts[verdicts$required & verdicts$status == "ok", ]
  failed <- unique(verdicts$stage[verdicts$pass %in% FALSE])
  if (length(failed) > 0L) {
    return(list(
      state = "not validated", stages = failed,
      text = sprintf(
        "complete, not validated (stages %s)",
        paste(titles(failed), collapse = ", ")
      )
    ))
  }
  # D80: validated means every required criterion was met, so one that
  # could not be evaluated withholds the word
  open <- unique(verdicts$stage[is.na(verdicts$pass)])
  if (length(open) > 0L) {
    return(list(
      state = "not judged", stages = open,
      text = sprintf(
        "complete, not judged (stages %s)",
        paste(titles(open), collapse = ", ")
      )
    ))
  }
  list(
    state = "validated", stages = character(),
    text = sprintf(
      "validated at plan %s (hash %s)",
      format_plan_version(plan$version), substr(plan$hash, 1L, 8L)
    )
  )
}

# comparing two runs --------------------------------------------------------

#' Compare the numbers of two validations
#'
#' Joins the [card_numbers()] of an earlier and a later run and flags the
#' numbers that moved by more than their Monte Carlo error, so a card
#' regenerated after a bmm release can say what changed (D84).
#'
#' @param old,new Each a `bmmtools_validation`, a path to its
#'   `result.rds` or its directory, or a [card_numbers()] tibble.
#'
#' @return A tibble with the key columns of [card_numbers()] (`stage`,
#'   `condition`, `level`, `term`, `metric`), then `old`, `new`, `diff`
#'   (`new - old`), `mcse_old`, `mcse_new` and `moved`: where both runs
#'   have an MCSE, `TRUE` when `|diff|` exceeds
#'   `2 * sqrt(mcse_old^2 + mcse_new^2)`; otherwise when it exceeds the
#'   new plan's tolerance for the metric (its `criteria$tolerance`, 0 by
#'   default, so any change of such a number is flagged; timings never
#'   are). A number present in one run only has `moved = NA`. The
#'   attributes `old_bmm` and `new_bmm` name the bmm versions of the two
#'   runs, and `n_compared` counts the numbers present in both.
#'
#' @export
compare_validations <- function(old, new) {
  old_result <- validation_or_numbers(old, "old")
  new_result <- validation_or_numbers(new, "new")
  keys <- c("stage", "condition", "level", "term", "metric")
  a <- numbers_of(old_result)
  b <- numbers_of(new_result)
  joined <- dplyr::full_join(
    a[c(keys, "value", "mcse")], b[c(keys, "value", "mcse")],
    by = keys, suffix = c("_old", "_new"), na_matches = "na",
    relationship = "one-to-one"
  )
  tolerance <- if (inherits(new_result, "bmmtools_validation")) {
    new_result$plan$criteria$tolerance
  } else {
    default_validation_criteria()$tolerance
  }
  tol <- tolerance[joined$metric]
  tol[is.na(tol)] <- tolerance[["default"]] %||% 0
  diff <- joined$value_new - joined$value_old
  with_mcse <- !is.na(joined$mcse_old) & !is.na(joined$mcse_new)
  moved <- ifelse(
    with_mcse,
    abs(diff) > 2 * sqrt(joined$mcse_old^2 + joined$mcse_new^2),
    abs(diff) > tol
  )
  # a number present in one run only stays NA: it was added or dropped,
  # which is not a move; a number missing in both did not move
  both_missing <- is.na(joined$value_old) & is.na(joined$value_new)
  moved[both_missing] <- FALSE
  out <- tibble::tibble(
    joined[keys],
    old = joined$value_old, new = joined$value_new, diff = diff,
    mcse_old = joined$mcse_old, mcse_new = joined$mcse_new,
    moved = as.logical(moved)
  )
  structure(
    out,
    old_bmm = bmm_version_of(old_result),
    new_bmm = bmm_version_of(new_result),
    n_compared = sum(!is.na(out$moved))
  )
}

#' @noRd
validation_or_numbers <- function(x, arg, call = rlang::caller_env()) {
  if (is.data.frame(x) && !inherits(x, "bmmtools_validation")) {
    needed <- names(empty_card_numbers())
    missing <- setdiff(needed, names(x))
    if (length(missing) > 0L) {
      cli::cli_abort(
        "{.arg {arg}} lacks the column{?s} {.val {missing}} of \\
         {.fn card_numbers}.",
        call = call
      )
    }
    return(x)
  }
  if (inherits(x, "bmmtools_validation")) {
    return(x)
  }
  tryCatch(
    read_previous_validation(x),
    error = function(e) {
      cli::cli_abort(
        c(
          "{.arg {arg}} must be a validation result, a path to one, or a \\
           {.fn card_numbers} tibble.",
          x = conditionMessage(e)
        ),
        call = call
      )
    }
  )
}

#' @noRd
numbers_of <- function(x) {
  if (inherits(x, "bmmtools_validation")) {
    return(x$numbers %||% card_numbers(x))
  }
  x
}

#' @noRd
bmm_version_of <- function(x) {
  if (!inherits(x, "bmmtools_validation") || is.null(x$run_info)) {
    return(NA_character_)
  }
  as.character(x$run_info$bmm_version[[1L]])
}

#' @noRd
check_validation_result <- function(result, call = rlang::caller_env()) {
  if (!inherits(result, "bmmtools_validation")) {
    cli::cli_abort(
      c(
        "{.arg result} must be a validation result, \\
         not {.obj_type_friendly {result}}.",
        i = "{.fn run_validation} makes one."
      ),
      call = call
    )
  }
  invisible(result)
}

# printing -------------------------------------------------------------------

#' Format and print a validation result
#'
#' @param x A `bmmtools_validation`.
#' @param ... Not used.
#' @return `format()` returns the lines; `print()` writes them and returns
#'   `x` invisibly.
#' @export
format.bmmtools_validation <- function(x, ...) {
  plan <- x$plan
  stages <- x$stages
  seconds <- ifelse(
    is.na(stages$seconds), "", sprintf("%.1f s", stages$seconds)
  )
  note <- ifelse(
    stages$status == "ok",
    ifelse(stages$resumed, "read", ""),
    ifelse(is.na(stages$message), "", stages$message)
  )
  c(
    sprintf(
      "<bmmtools validation> %s: %s", plan$model_name, x$headline$text
    ),
    sprintf(
      "Plan %s (hash %s), bmm %s",
      format_plan_version(plan$version), substr(plan$hash, 1L, 8L),
      bmm_version_of(x)
    ),
    sprintf(
      "  %2d %-24s %-15s %9s  %s",
      stages$order, stages$stage, stages$status, seconds, note
    )
  )
}

#' @rdname format.bmmtools_validation
#' @export
print.bmmtools_validation <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}
