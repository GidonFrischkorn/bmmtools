# The cache.
#
# brms's own `file =` cache detects a change to the Stan code, the Stan
# data and the algorithm, and nothing else (measured on brms 2.23.0). In
# a recovery grid the cells that share a design row differ in `seed`
# and sometimes in `iter`, `warmup`, `control` or `init`, none of which
# it sees. fit_cached() hashes everything that determines the fit and
# keeps the component hashes in a sidecar so that a mismatch can be
# reported by name, following miniQmetrics/R/bayes_helpers.R.

#' The arguments of `...` that enter the cache key
#' @noRd
cache_key_args <- function() {
  c(
    "seed", "chains", "iter", "warmup", "thin", "control", "init",
    "backend", "sample_prior", "algorithm", "stanvars"
  )
}

#' The two files of one cache entry
#'
#' `.rds` is appended when `file` has no extension; the key sits next to
#' it with `.key`.
#'
#' @noRd
cache_paths <- function(file) {
  bad <- !is.character(file) || length(file) != 1L || is.na(file) ||
    !nzchar(file)
  if (bad) {
    cli::cli_abort(
      "{.arg file} must be a single path, not {.obj_type_friendly {file}}."
    )
  }
  stem <- sub("\\.rds$", "", file)
  list(
    rds = paste0(stem, ".rds"),
    key = paste0(stem, ".key"),
    meta = paste0(stem, ".meta.rds")
  )
}

#' Write and read the fit's wall time, beside the fit and outside the key
#'
#' How long a fit took is a fact about the run, not about the fit, so it
#' must not enter the cache key: `read_cache_key()` turns every
#' `name=value` line of a key file into a component and
#' `changed_components()` unions the names of the stored and the fresh
#' key, so a `seconds=` line would differ from every key ever written and
#' refit every cached fit on every run. It goes in its own file, which
#' nothing compares.
#'
#' A fit cached before this file existed has no meta file, and reports
#' `NA` seconds rather than failing.
#'
#' @noRd
write_cache_meta <- function(path, seconds) {
  meta <- list(
    seconds = as.double(seconds),
    fitted_at = Sys.time(),
    bmmtools_version = bmmtools_version()
  )
  write_atomic(path, function(tmp) saveRDS(meta, tmp))
  meta
}

#' @noRd
read_cache_seconds <- function(path) {
  if (!file.exists(path)) {
    return(NA_real_)
  }
  meta <- tryCatch(readRDS(path), error = function(e) NULL)
  seconds <- meta$seconds
  if (!is.numeric(seconds) || length(seconds) != 1L) {
    return(NA_real_)
  }
  as.double(seconds)
}

#' Deparse every formula inside an object
#'
#' A formula carries its environment, and `rlang::hash()` of a
#' `bmmformula` differs between two R sessions for that reason alone
#' (measured 2026-09-07). Deparsing each formula element to one string,
#' recursively, leaves what the formula says and drops where it was
#' made.
#'
#' @noRd
formula_key <- function(x) {
  if (inherits(x, "formula")) {
    return(paste(deparse(x, width.cutoff = 500L), collapse = ""))
  }
  if (is.list(x)) {
    out <- lapply(unclass(x), formula_key)
    attributes(out) <- NULL
    names(out) <- names(x)
    return(out)
  }
  x
}

#' Reduce a function to its text so its environment does not enter
#'
#' `init` is commonly a closure (`init = function() list(kappa = 0)`),
#' and a closure hashes differently in every environment it is built
#' in, so a grid loop that builds one per cell would refit every cell.
#' Its formals and body are what determine the starting values, so they
#' are what enters the key. Non-functions pass through and are hashed
#' as they are.
#'
#' @noRd
function_key <- function(x) {
  if (is.function(x)) {
    return(paste(deparse(x, width.cutoff = 500L), collapse = "\n"))
  }
  x
}

#' Package version as a string, or a marker when it is not installed
#' @noRd
package_version_string <- function(pkg) {
  if (rlang::is_installed(pkg)) {
    as.character(utils::packageVersion(pkg))
  } else {
    "not installed"
  }
}

#' The version of the Stan toolchain the backend in use runs on
#'
#' bmm and brms generate the code; the backend compiles and samples, and
#' a new CmdStan can change the draws. Only the backend in use enters
#' the key, so that upgrading the other one does not invalidate a grid.
#' `backend = NULL` means brms's default, which is rstan unless the
#' `brms.backend` option says otherwise.
#'
#' @noRd
backend_version_string <- function(backend) {
  if (is.null(backend)) backend <- getOption("brms.backend", "rstan")
  backend <- as.character(backend)[[1L]]
  if (identical(backend, "cmdstanr")) {
    cmdstan <- "cmdstan not found"
    if (rlang::is_installed("cmdstanr")) {
      cmdstan <- tryCatch(
        as.character(cmdstanr::cmdstan_version()),
        error = function(e) "cmdstan not found"
      )
    }
    return(paste("cmdstanr", package_version_string("cmdstanr"), cmdstan))
  }
  if (identical(backend, "rstan")) {
    return(paste("rstan", package_version_string("rstan")))
  }
  backend
}

#' Hash every component that determines a fit
#'
#' @return A list with `components`, a named character vector of hashes,
#'   and `key`, the hash of that vector.
#' @noRd
cache_key <- function(formula, data, model, prior, dots) {
  sampler <- lapply(cache_key_args(), function(arg) dots[[arg]])
  names(sampler) <- cache_key_args()
  # `[<-` with list(), not `$<-`: assigning NULL would drop the element
  sampler["init"] <- list(function_key(sampler$init))
  # `draws` is to a Laplace or variational fit what `iter` is to sampling.
  # It enters only when given: a component added unconditionally would
  # differ from every key file written before it existed and refit every
  # cached sampling fit for a change that never touched it.
  if (!is.null(dots$draws)) sampler$draws <- dots$draws
  parts <- c(
    list(
      formula = formula_key(formula),
      data = data,
      model = model,
      prior = prior
    ),
    sampler,
    list(
      bmm = package_version_string("bmm"),
      brms = package_version_string("brms"),
      toolchain = backend_version_string(dots$backend)
    )
  )
  components <- vapply(parts, rlang::hash, character(1))
  list(components = components, key = rlang::hash(components))
}

#' Write a file atomically: to a temporary name, then rename
#'
#' A process killed inside `saveRDS()` would otherwise leave a truncated
#' `.rds` under the final name, which `refit = "never"` would then try to
#' read.
#'
#' @noRd
write_atomic <- function(path, writer) {
  tmp <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir = dirname(path)
  )
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  writer(tmp)
  if (!file.rename(tmp, path)) {
    cli::cli_abort("Could not write {.file {path}}.")
  }
  invisible(path)
}

#' Write the key file: one `component=hash` line each, then `key=`
#' @noRd
write_cache_key <- function(path, key) {
  lines <- c(
    paste0(names(key$components), "=", key$components),
    paste0("key=", key$key)
  )
  write_atomic(path, function(tmp) writeLines(lines, tmp))
}

#' Read a key file back into the shape `cache_key()` returns
#'
#' A missing or unreadable file gives `NULL`.
#'
#' @noRd
read_cache_key <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  lines <- readLines(path, warn = FALSE)
  lines <- lines[grepl("=", lines, fixed = TRUE)]
  if (length(lines) == 0L) {
    return(NULL)
  }
  name <- sub("=.*$", "", lines)
  value <- sub("^[^=]*=", "", lines)
  is_key <- name == "key"
  components <- value[!is_key]
  names(components) <- name[!is_key]
  list(
    components = components,
    key = if (any(is_key)) value[is_key][[1L]] else NA_character_
  )
}

#' Which components differ between two keys
#' @noRd
changed_components <- function(old, new) {
  all_names <- union(names(old$components), names(new$components))
  changed <- vapply(all_names, function(name) {
    !identical(
      unname(old$components[name]),
      unname(new$components[name])
    )
  }, logical(1))
  all_names[changed]
}

#' Refuse brms's own cache arguments
#'
#' brms's `file` cannot reach `...` at all, because `fit_cached()` has a
#' formal of that name and R refuses two arguments matched to one
#' formal; the two companions can, and are refused here.
#'
#' @noRd
check_no_brms_cache_args <- function(dots, call = rlang::caller_env()) {
  clash <- intersect(names(dots), c("file_refit", "file_compress"))
  if (length(clash) > 0L) {
    cli::cli_abort(
      c(
        "{.arg {clash}} cannot be passed through {.fn fit_cached}.",
        i = "There is one cache, keyed on everything that determines the \\
             fit; use {.arg file} and {.arg refit} of {.fn fit_cached}."
      ),
      call = call
    )
  }
  invisible(NULL)
}

#' Emit a cache message of a catchable class
#' @noRd
cache_inform <- function(..., .envir = rlang::caller_env()) {
  cli::cli_inform(c(...), class = "bmmtools_cache_message", .envir = .envir)
}

#' Decide whether the cached fit can be reused, and say why not
#'
#' @return A list with `reuse` (logical), `cached` (whether the `.rds`
#'   exists), `old_key` (`NULL` when no key file) and `changed` (the
#'   components that differ).
#' @noRd
cache_lookup <- function(paths, key, refit) {
  cached <- file.exists(paths$rds)
  old_key <- if (cached) read_cache_key(paths$key) else NULL
  changed <- character()
  if (!is.null(old_key)) {
    changed <- changed_components(old_key, key)
  }
  matches <- cached && !is.null(old_key) && length(changed) == 0L
  reuse <- switch(refit,
    always = FALSE,
    never = cached,
    on_change = matches
  )
  list(reuse = reuse, cached = cached, old_key = old_key, changed = changed)
}

#' Read a cached fit, warning when it is stale
#' @noRd
cache_read <- function(paths, key, lookup, call = rlang::caller_env()) {
  if (is.null(lookup$old_key)) {
    cli::cli_warn(c(
      "Reusing {.file {paths$rds}} although no cache key was recorded.",
      i = "{.code refit = \"never\"} keeps whatever is cached."
    ))
  } else if (length(lookup$changed) > 0L) {
    cli::cli_warn(c(
      "Reusing {.file {paths$rds}} although {.val {lookup$changed}} \\
       changed.",
      i = "{.code refit = \"never\"} keeps whatever is cached."
    ))
  } else {
    cache_inform("Reusing cached fit {.file {paths$rds}}.")
  }
  fit <- tryCatch(
    readRDS(paths$rds),
    error = function(e) {
      cli::cli_abort(
        c(
          "The cached fit {.file {paths$rds}} could not be read.",
          i = "It may be truncated; delete it or use \\
               {.code refit = \"always\"}."
        ),
        parent = e,
        call = call
      )
    }
  )
  attr(fit, "bmmtools_cache") <- list(
    file = paths$rds, key = key$key, reused = TRUE,
    seconds = read_cache_seconds(paths$meta)
  )
  fit
}

#' Say why a fit is being run
#' @noRd
cache_announce <- function(paths, lookup, refit) {
  if (lookup$cached && refit == "on_change") {
    if (is.null(lookup$old_key)) {
      cache_inform("Refitting {.file {paths$rds}}: no cache key recorded.")
    } else {
      cache_inform(
        "Refitting {.file {paths$rds}}: {.val {lookup$changed}} changed."
      )
    }
  } else {
    cache_inform("Fitting {.file {paths$rds}}.")
  }
  invisible(NULL)
}

#' Fit a bmm model once and reuse it while nothing has changed
#'
#' A cache for the fits a recovery study produces. brms's own `file`
#' argument reuses a saved fit when the Stan code, the Stan data and the
#' algorithm are unchanged, and does not look at the seed, the number of
#' chains or iterations, the warmup, thinning, `control`, `init` or the
#' backend (measured on brms 2.23.0). In a simulation grid the cells that
#' share a design differ in exactly those, so `fit_cached()` keys on all
#' of them, and on the formula, data, model, prior, the installed
#' versions of bmm and brms and the Stan toolchain of the backend in use,
#' and records each component's hash next to the fit so that a mismatch
#' is reported by name.
#'
#' @param formula,data,model,prior Passed to the fitter as they are; the
#'   formula enters the key deparsed, so the environment it was built in
#'   does not matter.
#' @param file Path of the cache file. `.rds` is appended when the path
#'   has no extension; the key is written next to it as `<file>.key`. The
#'   directory is created before fitting, so an unwritable path fails
#'   before any sampling.
#' @param refit `"on_change"` (the default) fits when nothing is cached or
#'   when any key component changed; `"never"` returns whatever is cached
#'   and warns if it is stale; `"always"` fits and overwrites.
#' @param ... Passed to the fitter. `seed`, `chains`, `iter`, `warmup`,
#'   `thin`, `control`, `init`, `backend`, `sample_prior`, `algorithm`
#'   and `stanvars` enter the key, and so does `draws` (the number of
#'   draws of a Laplace or variational fit) when it is given; anything
#'   else (`cores`, `refresh`, `silent`, ...) changes how the fit is run,
#'   not what it is, and does not. An `init` given as a function enters
#'   the key by its text, so a
#'   closure built afresh in every iteration of a loop still matches.
#'   brms's `file_refit` and `file_compress` are refused: there is one
#'   cache, not two.
#' @param .fitter The fitting function, `bmm::bmm()` by default. Tests
#'   inject a stand-in so that no model is compiled.
#'
#' @return The fit, with an attribute `bmmtools_cache`: a list of
#'   `file`, `key`, `reused` and `seconds`, the wall time the fit took
#'   when it was run, which is the first fit's time when this call reused
#'   one and `NA` for a fit cached before bmmtools recorded it.
#'
#' @details
#' The fit and the key are each written to a temporary name and renamed
#' into place, so an interrupted write never leaves a truncated file
#' under the final name; the key is written after the fit. A cached fit
#' that cannot be read is reported as such rather than as a parse error.
#'
#' How long the fit took is written to `<file>.meta.rds`, beside the fit
#' and outside the key. It is deliberately not a key component: a
#' component that changes with every run would differ from every stored
#' key and refit everything a resumed study had already fitted.
#'
#' @examples
#' \dontrun{
#' fit <- fit_cached(
#'   bmm::bmf(kappa ~ 1 + (1 | id), thetat ~ 1 + (1 | id)),
#'   data, bmm::mixture2p(resp_error = "y"),
#'   file = "fits/cell-01", seed = 1, chains = 4, iter = 2000
#' )
#' attr(fit, "bmmtools_cache")$reused
#' }
#'
#' @export
fit_cached <- function(formula,
                       data,
                       model,
                       file,
                       prior = NULL,
                       refit = c("on_change", "never", "always"),
                       ...,
                       .fitter = NULL) {
  refit <- rlang::arg_match(refit)
  dots <- rlang::list2(...)
  check_no_brms_cache_args(dots)
  fitter <- .fitter
  if (!is.null(fitter) && !is.function(fitter)) {
    cli::cli_abort(
      "{.arg .fitter} must be a function, not {.obj_type_friendly {fitter}}."
    )
  }

  paths <- cache_paths(file)
  key <- cache_key(formula, data, model, prior, dots)
  lookup <- cache_lookup(paths, key, refit)
  if (lookup$reuse) {
    return(cache_read(paths, key, lookup))
  }

  cache_announce(paths, lookup, refit)
  dir.create(dirname(paths$rds), recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(dirname(paths$rds))) {
    cli::cli_abort("Could not create the directory for {.file {paths$rds}}.")
  }

  if (is.null(fitter)) {
    rlang::check_installed("bmm", "to fit a model.")
    fitter <- bmm::bmm
  }
  started <- Sys.time()
  fit <- rlang::exec(
    fitter,
    formula = formula, data = data, model = model, prior = prior, !!!dots
  )
  seconds <- as.double(difftime(Sys.time(), started, units = "secs"))

  write_atomic(paths$rds, function(tmp) saveRDS(fit, tmp))
  meta <- write_cache_meta(paths$meta, seconds)
  write_cache_key(paths$key, key)
  attr(fit, "bmmtools_cache") <- list(
    file = paths$rds, key = key$key, reused = FALSE, seconds = meta$seconds
  )
  fit
}
