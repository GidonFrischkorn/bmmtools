# render_card(): the validation card of a stored result (spec 26, D82).
# Rendering needs the quarto package and the Quarto CLI; every test that
# renders skips without them, and the tests that do not render run on
# every runner.

skip_card_render <- function() {
  skip_validation_run()
  skip_on_cran()
  skip_if_not_installed("quarto")
  skip_if(is.null(quarto::quarto_path()), "the Quarto CLI is not installed")
}

#' A result with an earlier run to compare with and one model section
#' @noRd
card_result <- function(dir) {
  plan <- small_validation_plan()
  plan$sections <- list(
    list(title = "Mixture weights", numbers = "computation")
  )
  mock_validation(plan, dir, validation_mock_fitter())
  later <- mock_validation(plan, dir, validation_mock_fitter(),
    previous = file.path(dir, "result.rds")
  )
  # one moved row, so the changes section has a table
  later$changes$moved[[1L]] <- TRUE
  saveRDS(later, file.path(dir, "result.rds"))
  later
}

# rendering once per file: each render costs seconds
rendered <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      dir <- tempfile("bmmtools-card-")
      dir.create(dir)
      result <- card_result(dir)
      out <- file.path(dir, "card.html")
      value <- render_card(result, out)
      cache <<- list(
        dir = dir, result = result, out = out, value = value,
        html = paste(readLines(out, warn = FALSE), collapse = "\n"),
        text = card_visible_text(out)
      )
    }
    cache
  }
})

# arguments ----------------------------------------------------------------

test_that("render_card checks its arguments before it renders", {
  expect_error(render_card(NULL, tempfile()), "not `NULL`")
  expect_error(render_card(1, tempfile()), "validation result")
  expect_error(
    render_card(file.path(tempdir(), "no-such-result.rds"), tempfile()),
    "does not exist"
  )
  dir <- withr::local_tempdir()
  saveRDS(list(a = 1), file.path(dir, "other.rds"))
  expect_error(
    render_card(file.path(dir, "other.rds"), tempfile()), "validation result"
  )
  result <- structure(list(), class = "bmmtools_validation")
  expect_error(render_card(result, c("a", "b")), "single file path")
  expect_error(render_card(result, tempfile(), format = "pdf"), "html")
  expect_error(
    render_card(result, tempfile(), template = file.path(dir, "none.qmd")),
    "template"
  )
})

test_that("render_card names quarto and its CLI when either is missing", {
  result <- structure(list(), class = "bmmtools_validation")
  local_mocked_bindings(quarto_cli_path = function() NULL)
  skip_if_not_installed("quarto")
  expect_error(render_card(result, tempfile()), "Quarto CLI")
})

# the template ---------------------------------------------------------------

test_that("the template's prose carries no literal number", {
  template <- system.file("card", "card.qmd", package = "bmmtools")
  expect_true(nzchar(template))
  prose <- card_template_prose(template)
  expect_true(length(prose) > 0L)
  expect_identical(grep("[0-9]", prose, value = TRUE), character())
})

test_that("card_template_prose drops chunks and inline code only", {
  file <- withr::local_tempfile(fileext = ".qmd")
  writeLines(c(
    "---", "title: card", "---", "Coverage `r x[1]` of 3.", "```{r}",
    "y <- 2", "```", "Done."
  ), file)
  prose <- card_template_prose(file)
  expect_identical(prose, c("title: card", "Coverage  of 3.", "Done."))
})

# rendered -------------------------------------------------------------------

test_that("render_card writes one self-contained HTML file", {
  skip_card_render()
  card <- rendered()
  expect_identical(card$value, card$out)
  expect_true(file.exists(card$out))
  # embed-resources: no folder of assets beside the card
  expect_false(dir.exists(file.path(card$dir, "card_files")))
  # nothing is left beside the installed template
  shipped <- dirname(system.file("card", "card.qmd", package = "bmmtools"))
  expect_false(file.exists(file.path(shipped, "card.html")))
})

test_that("the card states the headline, the model and the plan", {
  skip_card_render()
  card <- rendered()
  result <- card$result
  text <- card$text
  expect_match(text, result$headline$text, fixed = TRUE)
  expect_match(text, "mixture2p", fixed = TRUE)
  expect_match(text, result$plan$version, fixed = TRUE)
  expect_match(text, substr(result$plan$hash, 1, 8), fixed = TRUE)
  expect_match(text, result$plan$owner, fixed = TRUE)
  expect_match(text, result$run_info$bmm_version[[1L]], fixed = TRUE)
  expect_match(text, result$run_info$brms_version[[1L]], fixed = TRUE)
})

test_that("every stage has a section, run or not", {
  skip_card_render()
  card <- rendered()
  for (title in card$result$stages$title) {
    expect_match(card$text, title, fixed = TRUE)
  }
  # in the stage's own section, not only in the verdicts' notes
  expect_match(card$text, "not built\\s*:\\s*not run \\(#17\\)")
  expect_match(
    card$text, "not applicable\\s*:\\s*no reference implementation"
  )
})

test_that("the verdicts print their criteria beside them", {
  skip_card_render()
  card <- rendered()
  for (criterion in unique(card$result$verdicts$criterion)) {
    expect_match(card$text, criterion, fixed = TRUE)
  }
  # one "not met" per criterion of a stage that ran and missed it
  ran <- card$result$verdicts$status == "ok"
  failed <- sum(ran & card$result$verdicts$pass %in% FALSE)
  expect_true(failed > 0L)
  shown <- lengths(gregexpr("not met", card$text, fixed = TRUE))
  expect_identical(shown, as.integer(failed))
})

test_that("a number with a Monte Carlo error prints it in brackets", {
  skip_card_render()
  card <- rendered()
  numbers <- card$result$numbers
  at <- numbers$stage == "prior_predictive" & numbers$metric == "floor_rate"
  row <- numbers[at, ][1, ]
  shown <- sprintf("%s (%s)",
    format(signif(row$value, 3), trim = TRUE),
    format(signif(row$mcse, 3), trim = TRUE)
  )
  expect_match(card$text, shown, fixed = TRUE)
})

test_that("the changes section lists the moved rows and the count", {
  skip_card_render()
  card <- rendered()
  changes <- card$result$changes
  expect_match(card$text, "Changes since bmm", fixed = TRUE)
  expect_match(
    card$text,
    sprintf(
      "%d of %d", sum(changes$moved, na.rm = TRUE), attr(changes, "n_compared")
    ),
    fixed = TRUE
  )
})

test_that("a model section prints its title and its numbers", {
  skip_card_render()
  card <- rendered()
  expect_match(card$text, "Mixture weights", fixed = TRUE)
})

test_that("every number on the card comes from the result", {
  skip_card_render()
  card <- rendered()
  expect_identical(card_unexplained_numbers(card$out, card$result), character())
  # and the card prints numbers at all
  expect_true(length(card_number_tokens(card$text)) > 20L)
})

test_that("a number typed into the template is caught", {
  skip_card_render()
  card <- rendered()
  template <- withr::local_tempfile(fileext = ".qmd")
  lines <- readLines(system.file("card", "card.qmd", package = "bmmtools"))
  writeLines(c(lines, "", "Coverage was 0.731 in the pilot."), template)
  out <- file.path(withr::local_tempdir(), "typed.html")
  render_card(file.path(card$dir, "result.rds"), out, template = template)
  expect_identical(card_unexplained_numbers(out, card$result), "0.731")
})

test_that("a template that writes elsewhere, or an unwritable output, fails", {
  skip_card_render()
  card <- rendered()
  path <- file.path(card$dir, "result.rds")
  template <- withr::local_tempfile(fileext = ".qmd")
  lines <- readLines(system.file("card", "card.qmd", package = "bmmtools"))
  lines <- sub("^format:$", "output-file: other.html\nformat:", lines)
  writeLines(lines, template)
  expect_error(
    render_card(path, tempfile(fileext = ".html"), template = template),
    "wrote no"
  )
  # an output whose directory is a file
  blocker <- withr::local_tempfile()
  writeLines("", blocker)
  # one error that carries the reason, and no warning beside it
  expect_no_warning(expect_error(
    render_card(path, file.path(blocker, "card.html")),
    "Could not write.*Not a directory"
  ))
})

test_that("a run with no earlier result says so, and a path is read", {
  skip_card_render()
  dir <- withr::local_tempdir()
  mock_validation(small_validation_plan(), dir, validation_mock_fitter())
  out <- file.path(dir, "card.html")
  # a directory is read as its result.rds
  render_card(dir, out)
  text <- card_visible_text(out)
  expect_match(text, "No earlier result was given", fixed = TRUE)
  expect_identical(
    card_unexplained_numbers(out, readRDS(file.path(dir, "result.rds"))),
    character()
  )
})

# the number check -----------------------------------------------------------

test_that("the visible text decodes entities and smart quotes", {
  file <- withr::local_tempfile(fileext = ".html")
  writeLines(c(
    "<html><head><title>9</title><style>p {margin: 3px}</style></head>",
    "<body><script>var x = 4;</script>",
    "<p>r_low &gt;= 0.80 &amp; SBC\u2019s band&nbsp;ok</p></body></html>"
  ), file)
  expect_identical(
    trimws(card_visible_text(file)), "r_low >= 0.80 & SBC's band ok"
  )
})

test_that("number tokens skip names, versions and issue numbers", {
  text <- paste(
    "mixture2p at bmm 1.4.1, see #17 and D83; coverage 0.95 (0.069),",
    "z -1.2, p = 1.07e-05, 6 x 20, 95%, r_low"
  )
  expect_identical(
    card_number_tokens(text),
    c("0.95", "0.069", "-1.2", "1.07e-05", "6", "20", "95")
  )
})

test_that("a token is explained by a value it rounds", {
  expect_true(card_token_explained("0.95", 0.9512))
  expect_true(card_token_explained("1.07e-05", 1.0678e-05))
  expect_true(card_token_explained("95", 0.95 * 100))
  expect_false(card_token_explained("0.731", c(0.7, 0.74)))
  expect_false(card_token_explained("3", numeric()))
  expect_true(card_token_explained("22000", 21987))
  expect_false(card_token_explained("22000", 22600))
  expect_false(card_token_explained("22", 21.4))
  expect_true(card_token_explained("0", 0.3))
  expect_false(card_token_explained("0", 0.6))
})
