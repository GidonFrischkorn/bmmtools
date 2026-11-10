# Every model and version the installed bmm offers has a simulation path in
# bmmtools, or is a known exception with a reason (#48).
#
# A path is a built-in adapter that was written for the version and checked to
# simulate it: its class is on `covered_classes` below, and adapter_name(), the
# dispatch simulate_recovery() uses, finds an adapter for it. adapter_name()
# alone is not enough. An adapter registered on a model's name takes every
# version of that model by inheritance, including one it cannot simulate. So
# any version an adapter takes that `covered_classes` does not record fails
# here, with or without an exception: an exception says bmmtools has no path,
# and an adapter that takes the version would simulate it anyway, wrongly if
# it was not written for it. A class on `covered_classes` that adapter_name()
# no longer takes fails too, as does a key on both lists.
#
# On bmm develop before PR #58, this test fails on the slot and
# slot_averaging versions of mixture2p and mixture3p and on every
# `+variable_precision` key: the circular adapters are registered on the
# model's name there, so they take these and would simulate a model they were
# not written for (a slot model as the `simple` one, a variable precision as
# a fixed one). PR #58 registers them on the version classes they simulate,
# and adapter_name() refuses `variable_precision = TRUE`.
#
# A simulate_data() method in bmm is *not* a path yet: bmmtools does not route
# generation through it until #22. Counting it would let this test pass over
# models bmmtools cannot simulate.
#
# A new model or version without a path turns this test red, on the develop
# leg first. That is the early warning #48 asks for. Answer it with an
# exception entry and its reason, or with a checked path; never by deleting or
# skipping the test. #22 removes the exceptions it gives a path to. The tests
# skip on CRAN: they are meant to fail when bmm adds a model, which must not
# fail a CRAN check of bmmtools.
#
# As of 2026-11, only bmm's development version exports simulate_data().
# Where bmm exports it, every key on `covered_classes` and
# `coverage_exceptions` must name a version bmm offers, so a misspelt key
# fails, and every exception must have a simulate_data() method, as its
# reason says. Where bmm does not export it, every key on
# `covered_release_only` must name a version bmm offers.
#
# The stub model this test hands adapter_name() carries only the model's and
# the version's class, so an adapter registered on a class in between (such
# as `sdt` or `circular`) would take real models the walk never sees. Every
# class the adapter table registers must therefore be a model's name or a key
# on `covered_classes` or `covered_release_only`.
#
# A model whose constructor has a `variable_precision` argument (imm,
# mixture2p and mixture3p on popov-lab/bmm@c13dd2c5, measured 2026-11-10; none
# in released bmm, 1.4.x as of 2026-11) is walked twice per version: once as
# the version's class, and once with `variable_precision = TRUE` under the key
# `<class>+variable_precision`. bmm gives both the same class, so only the key
# tells them apart. Options that change no parameters (the SDT models' `dist`,
# which their adapters read) need no rows of their own.

# The classes bmmtools' adapters were written for and checked to simulate.
# Measured with the real constructors on CRAN bmm 1.4.1 and
# popov-lab/bmm@c13dd2c5, 2026-11-10. A class the installed bmm does not offer
# is skipped, except by the unknown-key test.
covered_classes <- c(
  "cswald_simple", "cswald_crisk", "ddm", "ezdm_3par", "ezdm_4par",
  "imm_full", "imm_bsc", "imm_abc", "m3_custom", "m3_ss", "m3_cs",
  "mixture2p_simple", "mixture3p_simple", "sdm_simple",
  "sdt_mafc", "sdt_ranking", "sdt_yn"
)
# Released bmm's mixture2p and mixture3p (1.4.x, as of 2026-11) have no
# versions, so their class is the model's name. bmm develop gives them version
# classes instead, so these two are exempt from the unknown-key test and
# checked only where bmm does not export simulate_data(). Delete them once
# bmmtools' minimum bmm has versioned mixture2p and mixture3p.
covered_release_only <- c("mixture2p", "mixture3p")

# Why each model or version has no path, keyed like the walk: the version's
# class (`<model>_<version>`), the model's name when it has no versions, or
# `<class>+variable_precision`. An entry for a version the installed bmm does
# not offer is skipped, except by the unknown-key test. Both reasons measured
# on popov-lab/bmm@c13dd2c5, 2026-11-10: its simulate_data() returns data for
# every exception below.
no_path_reason <- paste(
  "No built-in adapter, and none is planned. bmm's development version",
  "simulates it with simulate_data(), which bmmtools is to route generation",
  "through with #22. Until then, pass `generator =`."
)
circular_reason <- paste(
  "The built-in circular adapters simulate a fixed precision and, for",
  "mixture2p and mixture3p, the `simple` version only: no slot capacity or",
  "set-size dependence, and no `tau`. No adapter is planned (#48). bmm's",
  "development version simulates it with simulate_data(), which bmmtools is",
  "to route generation through with #22. Until then, pass `generator =`."
)
# on bmm develop only
variable_precision_keys <- paste0(
  c(
    "imm_full", "imm_bsc", "imm_abc",
    "mixture2p_simple", "mixture2p_slot", "mixture2p_slot_averaging",
    "mixture3p_simple", "mixture3p_slot", "mixture3p_slot_averaging"
  ),
  "+variable_precision"
)
coverage_exceptions <- c(
  mpt = no_path_reason,
  sdt_cdp = no_path_reason,
  sdt_rating_standard = no_path_reason,
  sdt_rating_dpsdt = no_path_reason,
  sdt_rating_metad = no_path_reason,
  # this (from popov-lab/bmm#507) and the rest on bmm develop only
  psychometric_function = no_path_reason,
  mixture2p_slot = circular_reason,
  mixture2p_slot_averaging = circular_reason,
  mixture3p_slot = circular_reason,
  mixture3p_slot_averaging = circular_reason,
  stats::setNames(
    rep(circular_reason, length(variable_precision_keys)),
    variable_precision_keys
  )
)

# One row per model and version of the installed bmm, and one more per
# version with `variable_precision = TRUE` where the constructor takes it.
# `class` is the class a model object of that version carries last:
# `<model>_<version>`, or `<model>` for a model without versions (measured on
# CRAN 1.4.1 and popov-lab/bmm@c13dd2c5). `key` is `class`, or
# `<class>+variable_precision`.
installed_model_versions <- function() {
  models <- unclass(bmm::bmm_models())
  rows <- lapply(models, function(model) {
    constructor <- getExportedValue("bmm", model)
    args <- formals(constructor)
    versions <- if ("version" %in% names(args)) {
      eval(args$version, environment(constructor))
    } else {
      NA_character_
    }
    class <- ifelse(is.na(versions), model, paste0(model, "_", versions))
    fixed <- data.frame(
      model = model,
      version = versions,
      class = class,
      variable_precision = FALSE,
      key = class,
      stringsAsFactors = FALSE
    )
    if (!"variable_precision" %in% names(args)) {
      return(fixed)
    }
    variable <- fixed
    variable$variable_precision <- TRUE
    variable$key <- paste0(class, "+variable_precision")
    rbind(fixed, variable)
  })
  do.call(rbind, rows)
}

# Whether adapter_name() finds an adapter for row `i` of the walk. The stub
# carries only the classes the adapter table can name, the model's and the
# version's; the classes between `bmmodel` and the model's own (`circular`,
# `sdt`, ...) name no adapter.
adapter_takes <- function(offered, i) {
  object <- structure(
    list(
      version = offered$version[i],
      variable_precision = offered$variable_precision[i]
    ),
    class = unique(c("bmmodel", offered$model[i], offered$class[i]))
  )
  !is.na(adapter_name(object))
}

bmm_has_simulate_data <- function() {
  "simulate_data" %in% getNamespaceExports("bmm")
}

test_that("every bmm model and version has a path or a known exception", {
  skip_on_cran()
  skip_if_not_installed("bmm")
  offered <- installed_model_versions()
  expect_gt(nrow(offered), 0L)
  taken <- vapply(
    seq_len(nrow(offered)), adapter_takes, logical(1),
    offered = offered
  )
  recorded <- offered$key %in% c(covered_classes, covered_release_only)
  excepted <- offered$key %in% names(coverage_exceptions)
  bmm_version <- as.character(packageVersion("bmm"))

  unrecorded <- offered$key[taken & !recorded]
  expect(
    length(unrecorded) == 0L,
    paste0(
      "In bmm ", bmm_version, ", an adapter takes these model versions, but ",
      "`covered_classes` does not record them: ",
      paste(unrecorded, collapse = ", "), ". If the adapter cannot simulate ",
      "them, register it on the version classes it does simulate ",
      "(adapter_class() in R/adapters.R) and give each an exception; if it ",
      "can, check that and add them to `covered_classes`. A ",
      "`+variable_precision` key shares its version's class: make ",
      "adapter_name() refuse it (as PR #58 does) and give it an exception."
    )
  )

  registered <- unlist(lapply(adapter_classes(), adapter_class))
  in_between <- setdiff(
    registered,
    c(offered$model, covered_classes, covered_release_only)
  )
  expect(
    length(in_between) == 0L,
    paste0(
      "Adapters are registered on classes that are neither a model's name ",
      "nor a key on `covered_classes`: ", paste(in_between, collapse = ", "),
      ". This test's stub carries only the model's and the version's class, ",
      "so it cannot see which models such an adapter takes."
    )
  )

  gaps <- offered$key[!taken & !recorded & !excepted]
  expect(
    length(gaps) == 0L,
    paste0(
      "bmm ", bmm_version, " offers model versions without a simulation ",
      "path in bmmtools: ", paste(gaps, collapse = ", "), ". ",
      "Add each to `coverage_exceptions` with its reason."
    )
  )

  lost <- offered$key[recorded & !taken]
  expect(
    length(lost) == 0L,
    paste0(
      "In bmm ", bmm_version, ", no adapter takes these classes on ",
      "`covered_classes` any more: ", paste(lost, collapse = ", "), ". ",
      "Restore the adapter, or move each to `coverage_exceptions` with its ",
      "reason."
    )
  )

  both <- intersect(
    c(covered_classes, covered_release_only), names(coverage_exceptions)
  )
  expect(
    length(both) == 0L,
    paste0(
      "These keys are on both `covered_classes` and `coverage_exceptions`: ",
      paste(both, collapse = ", "), ". Keep each on one list."
    )
  )
})

test_that("where bmm exports simulate_data(), every listed key is offered", {
  skip_on_cran()
  skip_if_not_installed("bmm")
  skip_if_not(bmm_has_simulate_data(), "bmm does not export simulate_data()")
  offered <- installed_model_versions()
  unknown <- setdiff(
    c(covered_classes, names(coverage_exceptions)), offered$key
  )
  expect(
    length(unknown) == 0L,
    paste0(
      "bmm ", packageVersion("bmm"), " offers no model version with these ",
      "keys from `covered_classes` or `coverage_exceptions`: ",
      paste(unknown, collapse = ", "), ". Correct each key or remove it."
    )
  )
})

test_that("without bmm's simulate_data(), every release-only key is offered", {
  skip_on_cran()
  skip_if_not_installed("bmm")
  skip_if(bmm_has_simulate_data(), "bmm exports simulate_data()")
  offered <- installed_model_versions()
  unknown <- setdiff(covered_release_only, offered$key)
  expect(
    length(unknown) == 0L,
    paste0(
      "bmm ", packageVersion("bmm"), " offers no model version with these ",
      "keys from `covered_release_only`: ", paste(unknown, collapse = ", "),
      ". Correct each key or remove it."
    )
  )
})

test_that("where bmm exports simulate_data(), every exception has a method", {
  skip_on_cran()
  skip_if_not_installed("bmm")
  skip_if_not(bmm_has_simulate_data(), "bmm does not export simulate_data()")
  offered <- installed_model_versions()
  bmm_ns <- asNamespace("bmm")
  has_method <- function(class) {
    !is.null(utils::getS3method("simulate_data", class,
      optional = TRUE, envir = bmm_ns
    ))
  }
  for (key in intersect(names(coverage_exceptions), offered$key)) {
    row <- offered[offered$key == key, ]
    expect(
      has_method(row$class) || has_method(row$model),
      paste0(
        "bmm ", packageVersion("bmm"), " has no simulate_data() method for ",
        "the exception `", key, "`, whose reason says it has. Correct the ",
        "reason."
      )
    )
  }
})
