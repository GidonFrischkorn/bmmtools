# Every model and version the installed bmm offers has a simulation path in
# bmmtools, or is a known exception with a reason (#48).
#
# A path is a built-in adapter that accepts the version's class, decided by
# adapter_name(), the dispatch simulate_recovery() uses. A simulate_data()
# method in bmm is *not* a path yet: bmmtools does not route generation through
# it until #22. Counting it would let this test pass over models bmmtools
# cannot simulate.
#
# On the develop leg, a new model or version without a path turns this test
# red. That is the early warning #48 asks for. Answer it with an exception
# entry and its reason, or with a path; never by deleting or skipping the test.
# An exception that has gained a path fails too, so the list stays true: #22
# removes the entries it gives a path to.
#
# Known false pass on bmm develop until Milestone 41: the mixture2p and
# mixture3p adapters are registered under the model's name, so adapter_name()
# takes their `slot` and `slot_averaging` versions too, which those adapters
# cannot simulate. When 41 registers the circular adapters on their version
# class, mixture2p_slot, mixture2p_slot_averaging, mixture3p_slot and
# mixture3p_slot_averaging turn this test red and get entries below.

# Why each model or version has no path. Keyed by the version's class
# (`<model>_<version>`), or by the model's name when it has no versions. An
# entry for a model the installed bmm does not offer is ignored.
# `simulate_data = TRUE` records that bmm develop has a simulate_data() method
# for it (measured 2026-11-09 on develop c13dd2c5); the test checks this
# wherever bmm exports simulate_data().
no_path_reason <- paste(
  "No built-in adapter, and none is planned. bmm's development version has a",
  "simulate_data() method for it, which bmmtools routes generation through",
  "from #22 (Milestone 42). Until then, pass `generator =`."
)
coverage_exceptions <- list(
  mpt = list(reason = no_path_reason, simulate_data = TRUE),
  sdt_cdp = list(reason = no_path_reason, simulate_data = TRUE),
  sdt_rating_standard = list(reason = no_path_reason, simulate_data = TRUE),
  sdt_rating_dpsdt = list(reason = no_path_reason, simulate_data = TRUE),
  sdt_rating_metad = list(reason = no_path_reason, simulate_data = TRUE),
  # on bmm develop only (bmm #507)
  psychometric_function = list(reason = no_path_reason, simulate_data = TRUE)
)

# One row per model and version of the installed bmm. `class` is the class a
# model object of that version carries last: `<model>_<version>`, or `<model>`
# for a model without versions (measured on CRAN 1.4.1 and develop c13dd2c5).
installed_model_versions <- function() {
  models <- unclass(bmm::bmm_models())
  rows <- lapply(models, function(model) {
    constructor <- getExportedValue("bmm", model)
    versions <- if ("version" %in% names(formals(constructor))) {
      eval(formals(constructor)$version, environment(constructor))
    } else {
      NA_character_
    }
    data.frame(
      model = model,
      version = versions,
      class = ifelse(is.na(versions), model, paste0(model, "_", versions)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Whether adapter_name() finds an adapter for a model of this version. The
# object carries only the classes the adapter table can name; the classes
# between `bmmodel` and the model's own (`circular`, `sdt`, ...) name no
# adapter.
has_adapter <- function(model, version, class) {
  object <- structure(
    list(version = version),
    class = unique(c("bmmodel", model, class))
  )
  !is.na(adapter_name(object))
}

test_that("every bmm model and version has a path or a known exception", {
  skip_if_not_installed("bmm")
  offered <- installed_model_versions()
  expect_gt(nrow(offered), 0L)
  path <- mapply(has_adapter, offered$model, offered$version, offered$class)

  gaps <- offered$class[!path & !offered$class %in% names(coverage_exceptions)]
  expect(
    length(gaps) == 0L,
    paste0(
      "bmm ", packageVersion("bmm"), " offers model versions without a ",
      "simulation path in bmmtools: ", paste(gaps, collapse = ", "), ". ",
      "Add each to `coverage_exceptions` with its reason."
    )
  )

  stale <- offered$class[path & offered$class %in% names(coverage_exceptions)]
  expect(
    length(stale) == 0L,
    paste0(
      "These exceptions now have an adapter path; remove them from ",
      "`coverage_exceptions`: ", paste(stale, collapse = ", "), "."
    )
  )
})

test_that("every exception has a reason, and its simulate_data() claim holds", {
  skip_if_not_installed("bmm")
  for (name in names(coverage_exceptions)) {
    entry <- coverage_exceptions[[name]]
    expect_true(is.character(entry$reason) && nzchar(entry$reason),
      label = name
    )
    expect_true(is.logical(entry$simulate_data), label = name)
  }

  # checked only where bmm exports simulate_data() (not on CRAN 1.4.x)
  exported <- "simulate_data" %in% getNamespaceExports("bmm")
  offered <- installed_model_versions()
  claimed <- names(Filter(
    function(e) isTRUE(e$simulate_data), coverage_exceptions
  ))
  checked <- if (exported) intersect(claimed, offered$class) else character()
  for (cls in checked) {
    model <- offered$model[offered$class == cls]
    bmm_ns <- asNamespace("bmm")
    find <- function(class) {
      utils::getS3method("simulate_data", class,
        optional = TRUE, envir = bmm_ns
      )
    }
    method <- find(cls) %||% find(model)
    expect_false(is.null(method),
      label = paste("simulate_data() method for", cls)
    )
  }
})
