# Placeholder until Milestone 1 adds the score layer: keeps R CMD check
# from aborting on an empty test directory.
test_that("the package skeleton loads", {
  expect_true(is.character(utils::packageDescription("bmmtools")$Version))
})
