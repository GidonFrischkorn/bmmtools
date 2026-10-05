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
#' The signal-detection stack (`sdt_yn`, `sdt_mafc` and their `r*()` and
#' `d*()` functions) is in no released bmm, so a literal `bmm::rsdt_yn()`
#' makes `R CMD check`'s "checking dependencies in R code" report a
#' missing object on any machine with a released bmm installed. That is a
#' WARNING, and `r-lib/actions` checks with `error-on = "warning"`, so it
#' fails the workflow. Resolving the name here keeps the check quiet and
#' turns an absent symbol into a named error rather than R's bare "not an
#' exported object".
#'
#' Used only for the four fork-only names; every other bmm call in the
#' package stays a literal `bmm::` call, which is what documents the
#' dependency.
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
      i = "The signal-detection models are not in a released {.pkg bmm} yet."
    ))
  }
  getExportedValue("bmm", name)
}

#' The class name an adapter is registered under
#' @noRd
adapter_classes <- function() {
  c(
    "sdt_yn", "sdt_mafc", "ezdm", "ddm", "cswald", "mixture2p", "sdm",
    "mixture3p", "imm", "m3"
  )
}

#' @noRd
adapter_name <- function(model) {
  for (cls in adapter_classes()) {
    if (inherits(model, cls)) {
      return(cls)
    }
  }
  NA_character_
}

#' Look an adapter up by the model's class
#'
#' `NULL` when the model has none, including the four-parameter EZ
#' variant, whose column layout the adapter does not know.
#'
#' @noRd
generator_for <- function(model) {
  if (inherits(model, "ezdm_4par")) {
    return(NULL)
  }
  switch(adapter_name(model),
    sdt_yn = generate_sdt_yn,
    sdt_mafc = generate_sdt_mafc,
    ezdm = generate_ezdm,
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
#' as if it mattered. [sbc()]'s default generator takes exactly these
#' columns from its `data`, so a column the adapter writes itself
#' (`sdt_yn`'s `stimulus`) is never handed back to it.
#'
#' `mixture3p` and `imm` read the set size when the model names it as a
#' column (a number is every trial's set size and no column), the
#' non-target locations, and for `imm` `full` and `bsc` the non-target
#' distances (`NULL` on `abc`). `m3` reads its numbers of options when the
#' model names them as columns, in the order of `resp_cats`.
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

#' Name generated columns after the model's own column names
#' @noRd
name_columns <- function(data, names) {
  names(data) <- unlist(names, use.names = FALSE)
  data
}

#' Yes/no signal detection: one row per stimulus class
#' @noRd
generate_sdt_yn <- function(pars, n_trials, model) {
  stimulus <- c(1, 0)
  counts <- bmm_fun("rsdt_yn")(
    2L, n_trials, stimulus,
    d = pars$d, criterion = pars$criterion, sdratio = pars$sdratio,
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
#' @noRd
generate_sdt_mafc <- function(pars, n_trials, model) {
  counts <- bmm_fun("rsdt_mafc")(
    1L, n_trials,
    m = model$other_vars$m, d = pars$d, dist = model$other_vars$dist
  )
  name_columns(
    data.frame(counts, n_trials),
    list(model$resp_vars$response, model$other_vars$n_trials)
  )
}

#' EZ diffusion: one row of summary statistics
#' @noRd
generate_ezdm <- function(pars, n_trials, model) {
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
#' `n_trials` trials, as `sdt_mafc`. With them named as columns of a trial
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
