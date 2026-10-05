test_that("bmm_parameter_info() gives the table the tests read", {
  skip_if_not_installed("bmm")
  p <- bmm_parameter_info(bmm::mixture2p("y"))
  expect_named(p, c("parameter", "description", "fixed", "value", "link"))
  expect_identical(p$parameter[!p$fixed], c("kappa", "thetat"))
})

test_that("parameter_info() and parameters() agree where bmm has both", {
  skip_if_not_installed("bmm")
  skip_if_not(
    "parameter_info" %in% getNamespaceExports("bmm"),
    "installed bmm has no parameter_info()"
  )
  parameter_info <- getExportedValue("bmm", "parameter_info")
  parameters <- getExportedValue("bmm", "parameters")
  models <- list(
    bmm::mixture2p("y"),
    bmm::sdm("y"),
    bmm::m3(
      resp_cats = c("corr", "other", "npl"), num_options = c(1, 4, 5),
      choice_rule = "simple", version = "ss"
    )
  )
  for (model in models) {
    # parameters() warns that it is deprecated wherever parameter_info()
    # exists; the comparison is of what it returns
    expect_identical(
      suppressWarnings(parameters(model)), parameter_info(model)
    )
  }
  # and bmm_parameter_info() calls the one that does not warn
  expect_no_warning(bmm_parameter_info(bmm::mixture2p("y")))
})
