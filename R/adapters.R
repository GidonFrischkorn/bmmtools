# The generator adapters.
#
# bmm's r<model>() functions do not take the model's parameter names on
# the model's link scale (measured 2026-09-07 and 2026-09-08), so each
# supported model gets a small function that maps one subject's
# natural-scale parameters onto the generator's arguments and the
# generator's output onto the model's column names. Column names come
# from the model object, never from literals. A user-supplied
# `generator =` always takes precedence over the built-in adapter. An
# upstream change letting r<model>() take the model's own parameter
# names would empty this table.

#' Resolve a bmm function by name, at call time
#'
#' Written when the signal-detection stack (`sdt_yn`, `sdt_mafc` and their
#' `r*()` and `d*()` functions) was in no released bmm: a literal
#' `bmm::rsdt_yn()` then made `R CMD check`'s "checking dependencies in R
#' code" report a missing object, a WARNING that fails the workflow. CRAN
#' bmm 1.4.1, the minimum bmmtools declares since 0.3.0 (D93), exports the
#' whole SDT stack and `rmpt()`. The lookup is kept (D94) because it costs
#' nothing and turns a function a `develop` build lacks into a named error
#' rather than R's bare "not an exported object".
#'
#' Used only for the four SDT names it was written for; every other bmm
#' call in the package stays a literal `bmm::` call, which is what
#' documents the dependency.
#'
#' @param name The name of an exported bmm function.
#'
#' @return The function.
#' @noRd
bmm_fun <- function(name) {
  if (!requireNamespace("bmm", quietly = TRUE)) {
    cli::cli_abort("The {.pkg bmm} package is needed for {.fun {name}}.")
  }
  if (!name %in% getNamespaceExports("bmm")) {
    cli::cli_abort(c(
      "The installed {.pkg bmm} ({packageVersion('bmm')}) does not export \\
       {.fun {name}}.",
      i = "bmmtools needs {.pkg bmm} 1.4.1 or later, which exports it."
    ))
  }
  getExportedValue("bmm", name)
}

#' The class name an adapter is registered under
#'
#' A model with versions of different parameters gets one adapter per
#' version, registered under the version's own class (D100), so that a
#' version bmm adds later gets no adapter rather than a sibling's.
#' @noRd
adapter_classes <- function() {
  c(
    "sdt_yn", "sdt_mafc", "sdt_ranking", "ezdm", "ezdm_4par", "ddm",
    "cswald", "mixture2p", "sdm", "mixture3p", "imm", "m3"
  )
}

#' The class a model must carry for an adapter to take it
#'
#' The adapter's own name, except the three-parameter `ezdm`'s, which keeps
#' its 0.2.0 name (it is part of every stored `"adapter:ezdm"` record and
#' of the `validation_defaults()` entry) and takes only the version class
#' `ezdm_3par`.
#' @noRd
adapter_class <- function(name) {
  if (identical(name, "ezdm")) "ezdm_3par" else name
}

#' @noRd
adapter_name <- function(model) {
  for (name in adapter_classes()) {
    if (inherits(model, adapter_class(name))) {
      return(name)
    }
  }
  NA_character_
}

#' Look an adapter up by the model's class
#'
#' `NULL` when the model has none.
#'
#' @noRd
generator_for <- function(model) {
  switch(adapter_name(model),
    sdt_yn = generate_sdt_yn,
    sdt_mafc = generate_sdt_mafc,
    sdt_ranking = generate_sdt_ranking,
    ezdm = generate_ezdm,
    ezdm_4par = generate_ezdm_4par,
    ddm = generate_ddm,
    cswald = generate_cswald,
    mixture2p = generate_mixture2p,
    sdm = generate_sdm,
    mixture3p = generate_mixture3p,
    imm = generate_imm,
    m3 = generate_m3,
    NULL
  )
}

#' Whether a model's built-in generator reads the fit formula (D58)
#'
#' Only the custom `m3`, whose activation formulas exist nowhere but in the
#' user's formula. The engine hands `formula` to the adapter when this is
#' `TRUE`, and to nothing else: a user generator closes over its own.
#'
#' @noRd
adapter_reads_formula <- function(model) {
  inherits(model, "m3_custom")
}

#' The columns an adapter reads from a trial design
#'
#' The design-side half of the adapter table (D55): the per-trial columns
#' an adapter takes from `trial_design` rather than drawing or writing
#' itself, named from the model object. Empty for every adapter that
#' reads none, which is all seven of the 0.2.0 table: those refuse a
#' design, because a column nobody generated from would sit in the data
#' as if it mattered, except a count model, which may take one
#' ([count_design_columns()]). [sbc()]'s default generator takes exactly these
#' columns from its `data`, so a column the adapter writes itself
#' (`sdt_yn`'s `stimulus`) is never handed back to it.
#'
#' `mixture3p` and `imm` read the set size when the model names it as a
#' column (a number is every trial's set size and no column), the
#' non-target locations, and for `imm` `full` and `bsc` the non-target
#' distances (`NULL` on `abc`). `m3` reads its numbers of options when the
#' model names them as columns, in the order of `resp_cats`, and
#' `sdt_ranking` its set size `m` when the model names it as a column
#' (D101).
#'
#' @param model A `bmmodel`.
#' @return A character vector, empty when the adapter reads no design.
#' @noRd
trial_design_columns <- function(model) {
  key <- adapter_name(model)
  if (identical(key, "m3")) {
    options <- model$other_vars$num_options
    return(if (is.character(options)) unname(options) else character())
  }
  if (identical(key, "sdt_ranking")) {
    m <- model$other_vars$m
    return(if (is.character(m)) unname(m) else character())
  }
  if (!key %in% c("mixture3p", "imm")) {
    return(character())
  }
  vars <- model$other_vars
  regex <- intersect(
    attr(model, "regex_vars"), c("set_size", "nt_features", "nt_distances")
  )
  if (isTRUE(attr(model, "regex")) && length(regex) > 0L) {
    cli::cli_abort(c(
      "The built-in generator for {.cls {key}} cannot read a model whose \
       {.arg {regex}} {?is a regex/are regexes}.",
      i = "It needs the column names themselves: build the model with \
           {.code regex = FALSE} and list the columns."
    ))
  }
  set_size <- if (is.character(vars$set_size)) vars$set_size
  c(set_size, vars$nt_features, vars$nt_distances)
}

#' The columns a count adapter reads from a trial design, if given one
#'
#' One row of a count model holds many trials. Without a design, every
#' row holds `n_trials` of them, as in 0.2.0. With one, each design row is
#' one row of the data and carries its own number of trials, so a design
#' whose rows differ, such as 100 signal and 50 noise trials, can be
#' simulated (#26). The design is optional, unlike the required columns of
#' [trial_design_columns()].
#'
#' `sdt_yn` reads the stimulus class with the count, so that a row keeps
#' its class. `m3` with numbers of options and `sdt_ranking` have no
#' trial-count column of their own, so their design names one `n_trials`;
#' `m3` with option columns is a per-trial model and reads none.
#'
#' @param model A `bmmodel`.
#' @return A character vector, the count column last; empty for every
#'   adapter that is not a count model.
#' @noRd
count_design_columns <- function(model) {
  key <- adapter_name(model)
  if (is.na(key)) {
    return(character())
  }
  vars <- model$other_vars
  switch(key,
    sdt_yn = unname(c(vars$stimulus, vars$n_trials)),
    sdt_mafc = ,
    ezdm = ,
    ezdm_4par = unname(vars$n_trials),
    m3 = if (is.character(vars$num_options)) character() else "n_trials",
    sdt_ranking = "n_trials",
    character()
  )
}

#' The response columns whose total is a row's number of trials
#'
#' bmm takes the trials of a row of `m3` with numbers of options, or of
#' `sdt_ranking`, to be its total over the response columns
#' (`check_data.m3()`, `check_data.sdt_ranking()`), so a design's
#' `n_trials` is that total. `NULL` for every other model.
#'
#' @noRd
count_total_columns <- function(model) {
  key <- adapter_name(model)
  if (identical(key, "sdt_ranking")) {
    return(unname(model$resp_vars$response))
  }
  if (identical(key, "m3") && !is.character(model$other_vars$num_options)) {
    return(model$resp_vars$resp_cats)
  }
  NULL
}

#' Check a count design's values, naming the rows that are wrong
#'
#' Run before anything is drawn, so that a design error is never reported
#' as a property of the parameters or the prior. A count is a whole number
#' of at least 1; `bmm::rezdm()` refuses 2 or fewer ("n_trials must be
#' larger than 2"), and `bmm::rsdt_mafc()` returns `NA` for 2.5 without an
#' error (both measured 2026-10-05 on the installed 1.3.2.9000).
#'
#' @param in_data `TRUE` when the design was read from [sbc()]'s `data`,
#'   whose rows it keeps: the messages then name `data` and, for `m3` and
#'   `sdt_ranking`, the row total over the response columns, which is what
#'   the design's `n_trials` is.
#' @noRd
check_count_design <- function(design, model, in_data = FALSE,
                               call = rlang::caller_env()) {
  columns <- count_design_columns(model)
  if (length(columns) == 0L) {
    return(invisible(design))
  }
  # nolint start: object_usage_linter. Used by cli's glue interpolation.
  where <- function(column) {
    if (!in_data) {
      return(cli::format_inline("{.arg trial_design}'s column {.val {column}}"))
    }
    cats <- count_total_columns(model)
    if (identical(column, "n_trials") && !is.null(cats)) {
      # one line: format_inline() keeps a backslash-newline as it is
      return(cli::format_inline(
        "The total of each row of {.arg data} over {.val {cats}}"
      ))
    }
    cli::format_inline("Column {.val {column}} of {.arg data}")
  }
  # nolint end
  count <- columns[[length(columns)]]
  n <- design[[count]]
  minimum <- if (adapter_name(model) %in% c("ezdm", "ezdm_4par")) 3 else 1
  ok <- if (is.numeric(n)) {
    is.finite(n) & n == round(n) & n >= minimum
  } else {
    rep(FALSE, length(n))
  }
  if (!all(ok)) {
    bad <- which(!ok)
    cli::cli_abort(
      c(
        "{where(count)} must hold whole numbers of trials of at least \\
         {minimum}.",
        x = "Not so in row{?s} {as.character(bad)}.",
        i = "{cli::qty(length(bad))}{?It holds/They hold} {.val {n[bad]}}."
      ),
      call = call
    )
  }
  check_ranking_design(design, model, where, call = call)
  if (length(columns) == 2L) {
    stimulus <- design[[columns[[1L]]]]
    # a factor of "0" and "1" would pass `%in%`, which compares as text,
    # and fail in the generator
    bad <- if (is.numeric(stimulus) || is.logical(stimulus)) {
      which(!stimulus %in% c(0, 1))
    } else {
      seq_along(stimulus)
    }
    if (length(bad) > 0L) {
      cli::cli_abort(
        c(
          "{where(columns[[1L]])} must be 0 or 1, as numbers.",
          x = "Not so in row{?s} {as.character(bad)}.",
          i = "{cli::qty(length(bad))}{?It holds/They hold} \\
               {.val {stimulus[bad]}}."
        ),
        call = call
      )
    }
  }
  invisible(design)
}

#' Check the set sizes a ranking design gives, naming the rows that are wrong
#'
#' Only when the model's `m` is a column (D101) and the design has it: a
#' number, as bmm's data check requires ("Set-size column ... must be
#' numeric"), and a whole one from 2, the smallest ranking, to the number of
#' rank columns. A factor or text is refused here rather than by bmm once
#' the data are drawn.
#'
#' @param where The function [check_count_design()] names a column with.
#' @noRd
check_ranking_design <- function(design, model, where,
                                 call = rlang::caller_env()) {
  if (!identical(adapter_name(model), "sdt_ranking")) {
    return(invisible(design))
  }
  column <- trial_design_columns(model)
  if (length(column) == 0L || !column %in% names(design)) {
    return(invisible(design))
  }
  m <- design[[column]]
  largest <- length(model$resp_vars$response)
  bad <- if (is.numeric(m)) {
    which(is.na(m) | m != round(m) | m < 2 | m > largest)
  } else {
    seq_along(m)
  }
  if (length(bad) > 0L) {
    cli::cli_abort(
      c(
        "{where(column)} must hold whole numbers from 2 to {largest}, the \\
         number of rank columns.",
        x = "Not so in row{?s} {as.character(bad)}.",
        i = "{cli::qty(length(bad))}{?It holds/They hold} \\
             {.val {design[[column]][bad]}}."
      ),
      call = call
    )
  }
  invisible(design)
}

#' Name generated columns after the model's own column names
#' @noRd
name_columns <- function(data, names) {
  names(data) <- unlist(names, use.names = FALSE)
  data
}

#' The criterion argument of `rsdt_yn()` and `dsdt_yn()`: `criterion` or
#' `bias`, whichever member `pars` carries
#'
#' bmm 1.4.1 offers the two as alternatives (D91), and its functions take
#' exactly one ("Give exactly one of 'criterion' ... or 'bias'"), so the
#' absent one is left out rather than passed as `NULL`.
#' @noRd
sdt_location <- function(pars) {
  if (!is.null(pars$bias)) {
    return(list(bias = pars$bias))
  }
  list(criterion = pars$criterion)
}

#' Yes/no signal detection: one row per stimulus class
#'
#' With a trial design, one row per design row, each of its own stimulus
#' class and number of trials; the generator then returns the counts only.
#' @noRd
generate_sdt_yn <- function(pars, n_trials, model, trial_design = NULL) {
  if (!is.null(trial_design)) {
    columns <- count_design_columns(model)
    counts <- rlang::exec(
      bmm_fun("rsdt_yn"),
      nrow(trial_design), trial_design[[columns[[2L]]]],
      trial_design[[columns[[1L]]]],
      d = pars$d, !!!sdt_location(pars), sdratio = pars$sdratio,
      dist = model$other_vars$dist
    )
    return(name_columns(data.frame(counts), list(model$resp_vars$response)))
  }
  stimulus <- c(1, 0)
  counts <- rlang::exec(
    bmm_fun("rsdt_yn"),
    2L, n_trials, stimulus,
    d = pars$d, !!!sdt_location(pars), sdratio = pars$sdratio,
    dist = model$other_vars$dist
  )
  name_columns(
    data.frame(counts, stimulus, n_trials),
    list(
      model$resp_vars$response, model$other_vars$stimulus,
      model$other_vars$n_trials
    )
  )
}

#' m-alternative forced choice: one row of correct counts
#'
#' With a trial design, one row per design row of its own number of trials.
#' @noRd
generate_sdt_mafc <- function(pars, n_trials, model, trial_design = NULL) {
  if (!is.null(trial_design)) {
    counts <- bmm_fun("rsdt_mafc")(
      nrow(trial_design), trial_design[[count_design_columns(model)]],
      m = model$other_vars$m, d = pars$d, dist = model$other_vars$dist
    )
    return(name_columns(data.frame(counts), list(model$resp_vars$response)))
  }
  counts <- bmm_fun("rsdt_mafc")(
    1L, n_trials,
    m = model$other_vars$m, d = pars$d, dist = model$other_vars$dist
  )
  name_columns(
    data.frame(counts, n_trials),
    list(model$resp_vars$response, model$other_vars$n_trials)
  )
}

#' The SD ratio `rsdt_ranking()` and `dsdt_ranking()` take
#'
#' The model's `sdratio` (`dist = "normal"`) is the **log** of the ratio of
#' signal to noise SDs, on an identity link and fixed at 0, while bmm's
#' generator and density take the ratio itself (`rsdt_ranking()` refuses
#' one of 0 or below; measured 2026-10-25 on CRAN bmm 1.4.1). Passing the
#' parameter through would fail at the fixed 0 and, for a free value, draw
#' data under another ratio than the one the model fits. `gumbel_min` has
#' no `sdratio`: a ratio of 1.
#' @noRd
ranking_sdratio <- function(pars) {
  exp(pars$sdratio %||% 0)
}

#' Ranking: one row of rank counts
#'
#' Without a design, one row per subject (and task) of `n_trials` trials at
#' the model's set size; with one, one row per design row, of its
#' `n_trials` trials and, when the model names `m` as a column, its own set
#' size (D101). `rsdt_ranking()` returns as many rank columns as the
#' largest set size it is given, so a design's counts are padded with zeros
#' up to the model's rank columns, as bmm's data check wants beyond a row's
#' `m`. A number for `m` needs no padding: bmm's constructor refuses one
#' that differs from the number of rank columns.
#' @noRd
generate_sdt_ranking <- function(pars, n_trials, model, trial_design = NULL) {
  vars <- model$other_vars
  response <- unname(model$resp_vars$response)
  if (is.null(trial_design)) {
    rows <- 1L
    sizes <- n_trials
  } else {
    rows <- nrow(trial_design)
    sizes <- trial_design[[count_design_columns(model)]]
  }
  m <- if (is.character(vars$m)) {
    trial_design[[vars$m]]
  } else {
    vars$m
  }
  counts <- bmm::rsdt_ranking(
    rows, sizes,
    m = m, d = pars$d, sdratio = ranking_sdratio(pars), dist = vars$dist
  )
  out <- matrix(0, rows, length(response))
  out[, seq_len(ncol(counts))] <- counts
  name_columns(as.data.frame(out), as.list(response))
}

#' EZ diffusion: one row of summary statistics
#'
#' With a trial design, one row per design row of its own number of trials.
#' `bmm::rezdm()` takes a single `n_trials` (a vector fails with "the
#' condition has length > 1", measured 2026-10-05 on 1.3.2.9000 and CRAN
#' 1.3.2), so it is called once per row.
#' @noRd
generate_ezdm <- function(pars, n_trials, model, trial_design = NULL) {
  if (!is.null(trial_design)) {
    rows <- lapply(trial_design[[count_design_columns(model)]], function(n) {
      bmm::rezdm(
        1L, n,
        drift = pars$drift, bound = pars$bound, ndt = pars$ndt, s = pars$s
      )[c("mean_rt", "var_rt", "n_upper")]
    })
    return(name_columns(
      do.call(rbind, rows),
      list(
        model$resp_vars$mean_rt, model$resp_vars$var_rt,
        model$resp_vars$n_upper
      )
    ))
  }
  out <- bmm::rezdm(
    1L, n_trials,
    drift = pars$drift, bound = pars$bound, ndt = pars$ndt, s = pars$s
  )
  name_columns(
    out[c("mean_rt", "var_rt", "n_upper", "n_trials")],
    list(
      model$resp_vars$mean_rt, model$resp_vars$var_rt,
      model$resp_vars$n_upper, model$other_vars$n_trials
    )
  )
}

#' The four-parameter EZ diffusion's column names, in `rezdm()`'s order
#'
#' `rezdm(version = "4par")` returns the RT summaries per boundary, upper
#' first; the model names them `mean_rt = c(upper, lower)` and `var_rt =
#' c(upper, lower)` (`check_data.ezdm()`; measured 2026-10-25 on CRAN bmm
#' 1.4.1). The constructor accepts one name each and only bmm's data check
#' objects, so a model like that is refused here, before any data are drawn.
#' @noRd
ezdm_4par_columns <- function(model) {
  vars <- model$resp_vars
  if (length(vars$mean_rt) != 2L || length(vars$var_rt) != 2L) {
    cli::cli_abort(c(
      "The built-in generator for {.cls ezdm_4par} needs two column names \\
       in each of {.arg mean_rt} and {.arg var_rt}.",
      i = "The four-parameter version has RT summaries per boundary: \\
           {.code mean_rt = c(upper, lower)}, \\
           {.code var_rt = c(upper, lower)}."
    ))
  }
  c(
    mean_rt_upper = vars$mean_rt[[1L]], mean_rt_lower = vars$mean_rt[[2L]],
    var_rt_upper = vars$var_rt[[1L]], var_rt_lower = vars$var_rt[[2L]],
    n_upper = vars$n_upper
  )
}

#' Four-parameter EZ diffusion: one row of per-boundary summaries
#'
#' As [generate_ezdm()], with `zr` and `version = "4par"` passed on. A
#' boundary reached fewer than twice has `NA` summaries; bmm 1.4.1 keeps
#' such a row and warns how many there are (bmm #430), so they are left
#' as `rezdm()` draws them.
#' @noRd
generate_ezdm_4par <- function(pars, n_trials, model, trial_design = NULL) {
  columns <- ezdm_4par_columns(model)
  draw <- function(n) {
    bmm::rezdm(
      1L, n,
      drift = pars$drift, bound = pars$bound, ndt = pars$ndt, zr = pars$zr,
      s = pars$s, version = "4par"
    )
  }
  if (!is.null(trial_design)) {
    rows <- lapply(trial_design[[count_design_columns(model)]], function(n) {
      draw(n)[names(columns)]
    })
    return(name_columns(do.call(rbind, rows), as.list(columns)))
  }
  out <- draw(n_trials)
  name_columns(
    out[c(names(columns), "n_trials")],
    c(as.list(columns), list(model$other_vars$n_trials))
  )
}

#' Diffusion decision model: one row per trial
#' @noRd
generate_ddm <- function(pars, n_trials, model) {
  out <- bmm::rddm(
    n_trials,
    drift = pars$drift, bound = pars$bound, ndt = pars$ndt, zr = pars$zr
  )
  name_columns(
    out[c("rt", "response")],
    list(model$resp_vars$rt, model$resp_vars$response)
  )
}

#' Censored-shifted Wald: one row per trial
#'
#' `bmm::rcswald()` takes no `version` argument --- measured 2026-09-17 on
#' the installed 1.4.1.9000, on the study's pinned `develop` and in the
#' 1.3.2 source. It draws from `rtdists::rdiffusion()` with `a = bound`
#' and `z = zr * bound`, which is the `crisk` parameterisation, so the two
#' versions differ here in the mapping and not in the generator.
#'
#' **`simple` doubles `bound`.** Its `bound` is the distance from an
#' unbiased starting point to the correct boundary, half the separation
#' `rcswald()` takes; bmm's own `?cswald` says to multiply by 2 to get the
#' full separation. Measured 2026-09-17 on 20,000 draws generated with
#' `bound = 2`: profiled over `bound`, the `simple` likelihood peaks at
#' 0.995 and the `crisk` likelihood at 1.977. Passing `bound` straight
#' through would have made every `simple` recovery report a bias of
#' `log(2)` on the log link and read as a defect in bmm rather than in
#' this adapter.
#'
#' `zr` is 0.5 for `simple`, which has no starting-point parameter and is
#' defined against an unbiased start. `crisk` fixes `zr` at 0.5 too but a
#' design may free it, so its own value is passed.
#'
#' `sndt` is passed only when the model carries it: the fork's cswald has
#' it and both the pinned `develop` and released 1.3.2 have no `sndt`
#' anywhere (measured 2026-09-17), and `rcswald()` there has no such
#' argument to take. Reading it off the model rather than off the
#' installed version keeps one adapter right on all three.
#' @noRd
generate_cswald <- function(pars, n_trials, model) {
  crisk <- inherits(model, "cswald_crisk")
  args <- list(
    n = n_trials,
    drift = pars$drift,
    bound = if (crisk) pars$bound else 2 * pars$bound,
    ndt = pars$ndt,
    zr = if (crisk) pars$zr %||% 0.5 else 0.5,
    s = pars$s %||% 1
  )
  if (!is.null(pars$sndt)) {
    args$sndt <- pars$sndt
  }
  out <- do.call(bmm::rcswald, args)
  name_columns(
    out[c("rt", "response")],
    list(model$resp_vars$rt, model$resp_vars$response)
  )
}

#' Two-parameter mixture model: one response error per trial
#' @noRd
generate_mixture2p <- function(pars, n_trials, model) {
  y <- bmm::rmixture2p(
    n_trials,
    mu = pars$mu1, kappa = pars$kappa, p_mem = pars$thetat
  )
  name_columns(data.frame(y), list(model$resp_vars$resp_error))
}

#' Signal discrimination model: one response error per trial
#' @noRd
generate_sdm <- function(pars, n_trials, model) {
  y <- bmm::rsdm(n_trials, mu = pars$mu, c = pars$c, kappa = pars$kappa)
  name_columns(data.frame(y), list(model$resp_vars$resp_error))
}

# models with non-targets ----------------------------------------------------

#' The location parameter of a circular model, read off the model
#'
#' `mu1` on the installed bmm; the rebuilt circular mixture models of bmm
#' 1.5.0 rename it `mu`.
#'
#' @noRd
location_parameter <- function(model) {
  found <- grep("^mu1?$", names(model$parameters), value = TRUE)
  if (length(found) == 0L) "mu1" else found[[1L]]
}

#' Each trial's set size, non-target locations and distances
#'
#' Read from the design the way bmm's likelihood reads the data
#' (`check_data.non_targets()`): a trial of set size `k` has its lures in
#' the first `k - 1` non-target columns, and the columns beyond are
#' ignored. A lure inside the set size with no location would be read by
#' bmm as 0, the target's own location, and one with no distance as 999,
#' no similarity at all, so both are refused here instead.
#'
#' @return A list: `set_size`, and `lures` and `dist`, one vector per
#'   trial (`dist` `NULL` when the model names no distances).
#' @noRd
nontarget_trials <- function(trial_design, n_trials, model) {
  vars <- model$other_vars
  features <- vars$nt_features
  set_size <- if (is.character(vars$set_size)) {
    # a value that is no number becomes NA, refused just below by name
    suppressWarnings(as.numeric(as.character(trial_design[[vars$set_size]])))
  } else {
    rep(as.numeric(vars$set_size), n_trials)
  }
  largest <- length(features) + 1L
  bad <- is.na(set_size) | set_size < 1 | set_size %% 1 != 0 |
    set_size > largest
  if (any(bad)) {
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    t <- which(bad)[[1L]]
    cli::cli_abort(c(
      "Trial {t} of {.arg trial_design} has set size {set_size[[t]]}.",
      i = "The model names {length(features)} non-target column{?s}, so a \\
           set size is a whole number from 1 to {largest}."
    ))
  }
  per_trial <- function(columns, read_as) {
    if (is.null(columns)) {
      return(NULL)
    }
    values <- as.matrix(trial_design[columns])
    lapply(seq_len(n_trials), function(t) {
      k <- seq_len(set_size[[t]] - 1L)
      row <- unname(as.double(values[t, k]))
      if (anyNA(row)) {
        cli::cli_abort(c(
          "Trial {t} of {.arg trial_design} has set size {set_size[[t]]}, \\
           but {.val {columns[k][is.na(row)]}} {?is/are} {.code NA}.",
          i = "A lure inside the set size needs a value; bmm would read \\
               {.code NA} there as {read_as}."
        ))
      }
      row
    })
  }
  list(
    set_size = set_size,
    lures = per_trial(features, "0, the target's location"),
    dist = per_trial(vars$nt_distances, "999, an unrelated item")
  )
}

#' One trial's mixture weights for `mixture3p`, from its link-scale weights
#'
#' bmm's likelihood gives the target the log weight `thetat`, a guessing
#' component 0, and each of the `set_size - 1` lures `thetant -
#' log(set_size - 1)`, so the lures share `exp(thetant)` (measured in the
#' generated Stan code, spec § 9.1). A set-size-1 trial has no lure and no
#' lure weight. The softmax is shifted, per trial, by the largest term
#' that trial has, so that no value overflows.
#'
#' @param thetat,thetant One subject's values, on the link scale (D56).
#' @param set_size The set size of each trial.
#' @return A list of `p_mem` and `p_nt`, one value per trial.
#' @noRd
mixture3p_weights <- function(thetat, thetant, set_size) {
  lure <- set_size > 1
  shift <- pmax(thetat, ifelse(lure, thetant, -Inf), 0)
  target <- exp(thetat - shift)
  lures <- ifelse(lure, exp(thetant - shift), 0)
  total <- target + lures + exp(-shift)
  p_mem <- target / total
  # rmixture3p() refuses p_mem + p_nt > 1, which rounding can reach when
  # guessing underflows to 0
  p_nt <- pmin(lures / total, 1 - p_mem)
  list(p_mem = p_mem, p_nt = p_nt)
}

#' Three-parameter mixture model: one response error per trial
#'
#' One `rmixture3p()` call per trial, since the call fixes the locations
#' for all its draws: 0.13 s per 500 trials (spec § 1). `thetat` and
#' `thetant` arrive on the link scale (D56).
#'
#' @noRd
generate_mixture3p <- function(pars, n_trials, model, trial_design = NULL) {
  trials <- nontarget_trials(trial_design, n_trials, model)
  weights <- mixture3p_weights(pars$thetat, pars$thetant, trials$set_size)
  mu <- pars[[location_parameter(model)]]
  y <- vapply(seq_len(n_trials), function(t) {
    bmm::rmixture3p(
      1L,
      mu = c(mu, trials$lures[[t]]), kappa = pars$kappa,
      p_mem = weights$p_mem[[t]], p_nt = weights$p_nt[[t]]
    )
  }, double(1))
  name_columns(data.frame(y), list(model$resp_vars$resp_error))
}

#' One trial's `rimm()` arguments, for each version of `imm`
#'
#' bmm's likelihood gives the target the weight `c + a` (`bsc`: `c`), a
#' lure at distance `d` the weight `c * exp(-s * d) + a` (`bsc`: `c *
#' exp(-s * d)`; `abc`: `a`), and guessing 1, all on the natural scale
#' (spec § 9.1). `dimm()` computes `c * exp(-s * dist) + a` against
#' `b`, so `b = 1`, `bsc` passes `a = 0`, and `abc` puts its lures at an
#' infinite distance, which leaves them `a` alone.
#'
#' @param lures,dist One trial's non-target locations and distances.
#' @noRd
imm_args <- function(pars, model, lures, dist) {
  abc <- inherits(model, "imm_abc")
  list(
    mu = c(pars[[location_parameter(model)]], lures),
    dist = c(0, if (abc) rep(Inf, length(lures)) else dist),
    c = pars$c,
    a = if (inherits(model, "imm_bsc")) 0 else pars$a,
    b = 1,
    s = if (abc) 1 else pars$s,
    kappa = pars$kappa
  )
}

#' Interference measurement model: one response error per trial
#' @noRd
generate_imm <- function(pars, n_trials, model, trial_design = NULL) {
  trials <- nontarget_trials(trial_design, n_trials, model)
  y <- vapply(seq_len(n_trials), function(t) {
    args <- imm_args(pars, model, trials$lures[[t]], trials$dist[[t]])
    do.call(bmm::rimm, c(list(n = 1L), args))
  }, double(1))
  name_columns(data.frame(y), list(model$resp_vars$resp_error))
}

# m3 -------------------------------------------------------------------------

#' The custom `m3`'s activation formulas, one per category, in `resp_cats`
#' order
#'
#' Read off the fit formula (D57, D58), which holds the activations next to
#' the parameter formulas. Put in `resp_cats` order here, because released
#' bmm 1.3.2's `rm3()` applies each number of options to whichever
#' activation stands at its position (fixed in bmm's development version,
#' which reorders them itself).
#'
#' Every free parameter must appear in an activation: one that does not
#' would get a truth the data never depended on.
#'
#' @param formula The `bmmformula` the custom model is fitted with.
#' @return The activation formulas, a `bmmformula`, with the names they
#'   use as the attribute `"used"`.
#' @noRd
m3_activations <- function(formula, model, call = rlang::caller_env()) {
  cats <- model$resp_vars$resp_cats
  if (is.null(formula)) {
    cli::cli_abort(
      c(
        "The built-in generator for a custom {.cls m3} needs {.arg formula}.",
        i = "Its activation formulas, one per category ({.val {cats}}), \
             exist only in the formula the model is fitted with: pass it \
             as {.arg formula}."
      ),
      call = call
    )
  }
  missing <- setdiff(cats, names(formula))
  if (length(missing) > 0L) {
    cli::cli_abort(
      c(
        "{.arg formula} has no activation formula for {.val {missing}}.",
        i = "A custom {.cls m3} needs one per category of \
             {.arg resp_cats}: {.val {cats}}."
      ),
      call = call
    )
  }
  activations <- formula[cats]
  used <- unique(unlist(lapply(activations, function(f) {
    all.vars(f[[length(f)]])
  })))
  info <- model_parameters(model)
  unused <- setdiff(info$free, used)
  if (length(unused) > 0L) {
    cli::cli_abort(
      c(
        "No activation formula uses {.val {unused}}, which the model \
         estimates.",
        i = "Its value would be a truth the data never depended on: use \
             it in an activation, or drop its link."
      ),
      call = call
    )
  }
  unknown <- setdiff(used, c(info$free, info$fixed))
  if (length(unknown) > 0L) {
    cli::cli_abort(
      c(
        "The activation formulas use {.val {unknown}}, which {?is no/are \
         no} parameter{?s} of the model.",
        i = "Its parameters are {.val {c(info$free, info$fixed)}}; give \
             each one a link with {.code m3(..., links = list(...))}."
      ),
      call = call
    )
  }
  attr(activations, "used") <- used
  activations
}

#' The numbers of options, as a numeric vector in `resp_cats` order
#'
#' Named `n_opt_<category>`, as `m3()` names them, and positional, which
#' every bmm reads the same way: numbers named after the categories in
#' another order are put in order here, because released bmm 1.3.2 ignored
#' those names (bmm #449). With column names, `row` is one trial of the
#' design and its numbers are read from it (a factor as the numbers it
#' shows, not its codes).
#'
#' @noRd
m3_options <- function(model, row = NULL, trial = NULL) {
  cats <- model$resp_vars$resp_cats
  options <- model$other_vars$num_options
  if (is.character(options)) {
    values <- suppressWarnings(as.double(vapply(
      row[unname(options)], function(v) as.character(v[[1L]]), character(1)
    )))
    bad <- !is.finite(values) | values < 0
    if (any(bad)) {
      cli::cli_abort(c(
        "Trial {trial} of {.arg trial_design} has {.val \
         {unname(options)[bad]}} = {values[bad]}.",
        i = "A number of options is a finite count of at least 0."
      ))
    }
    if (sum(values) == 0) {
      cli::cli_abort(c(
        "Trial {trial} of {.arg trial_design} has no option in any category.",
        i = "At least one of {.val {unname(options)}} must be above 0."
      ))
    }
  } else if (setequal(names(options), cats)) {
    values <- unname(options[cats])
  } else {
    values <- unname(options)
  }
  stats::setNames(as.double(values), paste0("n_opt_", cats))
}

#' One `rm3()` call: `n` rows of counts in the columns `resp_cats` names
#' @noRd
m3_counts <- function(size, pars, model, options, activations, n = 1L) {
  model$other_vars$num_options <- options
  counts <- bmm::rm3(n, size, pars, model, act_funs = activations)
  cats <- model$resp_vars$resp_cats
  if (setequal(colnames(counts), cats)) counts[, cats, drop = FALSE] else counts
}

#' Multinomial measurement model: category counts
#'
#' With numbers of options on the model, one row per subject (and task) of
#' `n_trials` trials, as `sdt_mafc`, or with a count design one row per
#' design row of its `n_trials` trials. With them named as columns of a trial
#' design, one row per trial, a single response counted in its category,
#' drawn with one `rm3(n, 1, ...)` call per distinct row of option counts
#' (one call per trial cost 0.55 ms, measured in the 21.3 review). The
#' parameters passed are those the activations use, which for the custom
#' version is every free one (`m3_activations()`) and the fixed `b` when
#' an activation reads it.
#' `ss` and `cs` take bmm's own activation formulas; the custom version
#' reads them off `formula` (D58).
#'
#' @noRd
generate_m3 <- function(pars, n_trials, model, trial_design = NULL,
                        formula = NULL) {
  activations <- if (adapter_reads_formula(model)) {
    m3_activations(formula, model)
  }
  values <- unlist(pars)
  if (!is.null(activations)) {
    values <- values[intersect(names(values), attr(activations, "used"))]
  }
  if (is.null(trial_design)) {
    counts <- m3_counts(
      n_trials, values, model, m3_options(model), activations
    )
  } else if (!is.character(model$other_vars$num_options)) {
    # numbers of options: a count design, each row its own total (#26).
    # rm3() reads only the first `size` of a vector (two rows of 100 and
    # 50 both summed to 100, measured 2026-10-05 on 1.3.2.9000 and CRAN
    # 1.3.2), so one call per distinct total
    sizes <- trial_design$n_trials
    counts <- NULL
    for (size in unique(sizes)) {
      rows <- which(sizes == size)
      drawn <- m3_counts(
        size, values, model, m3_options(model), activations,
        n = length(rows)
      )
      if (is.null(counts)) {
        counts <- matrix(0, length(sizes), ncol(drawn),
          dimnames = list(NULL, colnames(drawn))
        )
      }
      counts[rows, ] <- drawn
    }
  } else {
    options <- lapply(seq_len(n_trials), function(t) {
      m3_options(model, trial_design[t, , drop = FALSE], t)
    })
    key <- vapply(options, paste, character(1), collapse = "\r")
    counts <- NULL
    for (k in unique(key)) {
      rows <- which(key == k)
      drawn <- m3_counts(
        1L, values, model, options[[rows[[1L]]]], activations,
        n = length(rows)
      )
      if (is.null(counts)) {
        counts <- matrix(0, n_trials, ncol(drawn),
          dimnames = list(NULL, colnames(drawn))
        )
      }
      counts[rows, ] <- drawn
    }
  }
  name_columns(
    as.data.frame(counts, row.names = NULL),
    as.list(model$resp_vars$resp_cats)
  )
}
