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

#' One stage's section of a card's text: from the title's last occurrence
#' (the contents and the verdicts table come first) to the next stage's
#' @noRd
card_section <- function(text, title, next_title) {
  pieces <- strsplit(text, title, fixed = TRUE)[[1L]]
  strsplit(pieces[[length(pieces)]], next_title, fixed = TRUE)[[1L]][[1L]]
}

#' The small plan run against the mock under this calibration prior, and
#' its card
#' @noRd
render_prior_card <- function(sbc, envir = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = envir)
  design <- list(
    n_subjects = 4L, n_trials = 4L, n_sims = 5L, thin_ranks = 1,
    cores_per_fit = 1
  )
  plan <- small_validation_plan(sbc = c(design, sbc))
  result <- mock_validation(plan, dir, validation_mock_fitter())
  out <- file.path(dir, "card.html")
  render_card(result, out)
  list(result = result, out = out, text = card_visible_text(out))
}

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
  expect_match(card$text, "not built\\s*:\\s*not run \\(#6\\)")
  expect_match(
    card$text,
    "not applicable\\s*:\\s*the grid has no design with 50 subjects"
  )
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
  shown <- sprintf(
    "%s (%s)",
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
  # the rule's multiplier is the comparison's, not typed in the template
  expect_match(
    card$text,
    sprintf(
      "by more than %d times their combined", attr(changes, "mcse_multiplier")
    ),
    fixed = TRUE
  )
})

test_that("a model section prints its title and its numbers", {
  skip_card_render()
  card <- rendered()
  expect_match(card$text, "Mixture weights", fixed = TRUE)
})

test_that("the calibration's section prints the entry's prior and its scope", {
  skip_card_render()
  card <- rendered()
  plan <- card$result$plan
  # the small plan keeps the mixture2p entry's prior and reason
  expect_false(is.null(plan$sbc$prior))
  section <- card_section(card$text, "Implementation (SBC)", "Recovery")
  for (row in plan$sbc$prior$prior) {
    expect_match(section, row, fixed = TRUE)
    # listed once: the prior check points here instead of repeating them
    expect_length(gregexpr(row, card$text, fixed = TRUE)[[1L]], 1L)
  }
  for (nlpar in plan$sbc$prior$nlpar) {
    expect_match(section, nlpar, fixed = TRUE)
  }
  expect_match(section, plan$sbc$prior_reason, fixed = TRUE)
  expect_no_match(section, "prior: bmm's default priors", fixed = TRUE)
  expect_no_match(section, "No reason was given", fixed = TRUE)
  expect_match(
    section, "the recovery stage fits under bmm's default priors",
    fixed = TRUE
  )
  # the prior check is judged under the same prior, and says so
  pp <- card_section(card$text, "Prior predictive", "Implementation (SBC)")
  expect_match(
    pp, "judged on the draws under the calibration's prior",
    fixed = TRUE
  )
})

test_that("a calibration under bmm's defaults says so in both sections", {
  skip_card_render()
  card <- render_prior_card(list(prior = NULL))
  expect_null(card$result$plan$sbc$prior)
  section <- card_section(card$text, "Implementation (SBC)", "Recovery")
  expect_match(
    section, "The calibration's prior: bmm's default priors.",
    fixed = TRUE
  )
  # the entry's reason went with its prior
  expect_no_match(section, "overflows", fixed = TRUE)
  expect_no_match(section, "No reason was given", fixed = TRUE)
  expect_match(
    section, "the recovery stage fits under bmm's default priors",
    fixed = TRUE
  )
  pp <- card_section(card$text, "Prior predictive", "Implementation (SBC)")
  expect_match(
    pp, "judged on the draws under bmm's default priors.",
    fixed = TRUE
  )
  expect_no_match(pp, "calibration's prior", fixed = TRUE)
  expect_identical(card_unexplained_numbers(card$out, card$result), character())
})

test_that("an own prior without a reason is printed, and says it has none", {
  skip_card_render()
  own <- data.frame(
    prior = "normal(1.5, 0.4)", class = "b", coef = "Intercept",
    nlpar = "kappa"
  )
  card <- render_prior_card(list(prior = own))
  expect_null(card$result$plan$sbc$prior_reason)
  section <- card_section(card$text, "Implementation (SBC)", "Recovery")
  expect_match(section, "normal(1.5, 0.4)", fixed = TRUE)
  expect_match(section, "No reason was given for this prior.", fixed = TRUE)
  expect_no_match(section, "overflows", fixed = TRUE)
  expect_no_match(section, "prior: bmm's default priors", fixed = TRUE)
  expect_identical(card_unexplained_numbers(card$out, card$result), character())
})

#' The rendered card's result changed by `change`, rendered again
#' @noRd
rerender <- function(change, envir = parent.frame()) {
  result <- change(rendered()$result)
  out <- file.path(withr::local_tempdir(.local_envir = envir), "card.html")
  render_card(result, out)
  list(result = result, out = out, text = card_visible_text(out))
}

test_that("a prior on dpars names the parameter of each row", {
  # the ezdm prior has three SD rows that differ only in their dpar
  skip_card_render()
  prior <- validation_defaults("ezdm")$sbc$prior
  card <- rerender(function(result) {
    result$plan$sbc$prior <- prior
    result
  })
  section <- card_section(card$text, "Implementation (SBC)", "Recovery")
  expect_match(section, "prior class dpar", fixed = TRUE)
  rows <- paste(prior$prior, prior$class, prior$dpar)
  for (row in rows) expect_match(section, row, fixed = TRUE)
  expect_identical(card_unexplained_numbers(card$out, card$result), character())
})

test_that("a bounded brmsprior prints its bound", {
  skip_card_render()
  skip_if_not_installed("brms")
  card <- rerender(function(result) {
    result$plan$sbc$prior <- brms::set_prior(
      "normal(1.5, 0.4)",
      class = "b", coef = "Intercept", nlpar = "kappa", lb = 0.1
    )
    result
  })
  section <- card_section(card$text, "Implementation (SBC)", "Recovery")
  expect_match(section, "nlpar lb", fixed = TRUE)
  expect_match(
    section, "normal(1.5, 0.4) b Intercept kappa 0.1",
    fixed = TRUE
  )
  # brms' bookkeeping columns are not the prior
  expect_no_match(section, "source", fixed = TRUE)
  expect_identical(card_unexplained_numbers(card$out, card$result), character())
})

test_that("a prior without coef, dpar or nlpar columns prints", {
  skip_card_render()
  card <- rerender(function(result) {
    result$plan$sbc$prior <- data.frame(prior = "normal(0, 0.3)", class = "sd")
    result
  })
  section <- card_section(card$text, "Implementation (SBC)", "Recovery")
  expect_match(section, "prior class normal(0, 0.3) sd", fixed = TRUE)
  expect_no_match(section, "coef", fixed = TRUE)
})

test_that("a stale calibration does not claim the plan's prior ran", {
  skip_card_render()
  card <- rerender(function(result) {
    at <- result$stages$stage == "sbc"
    result$stages$status[at] <- "stale"
    result$stages$message[at] <- paste(
      "read under `refit = \"never\"`, made from other inputs or",
      "software versions than this plan's"
    )
    result
  })
  section <- card_section(card$text, "Implementation (SBC)", "Recovery")
  expect_match(section, "stale", fixed = TRUE)
  expect_match(section, "may not be the one it ran under", fixed = TRUE)
  expect_match(section, "The current plan's prior:", fixed = TRUE)
  expect_no_match(section, "The calibration's prior:", fixed = TRUE)
  for (row in card$result$plan$sbc$prior$prior) {
    expect_match(section, row, fixed = TRUE)
  }
})

test_that("an excluded calibration's prior is listed by the prior check", {
  skip_card_render()
  dir <- withr::local_tempdir()
  plan <- small_validation_plan(
    stages = setdiff(validation_stage_names(), "sbc")
  )
  result <- mock_validation(plan, dir, validation_mock_fitter())
  expect_identical(
    result$stages$status[result$stages$stage == "sbc"], "excluded"
  )
  out <- file.path(dir, "card.html")
  render_card(result, out)
  text <- card_visible_text(out)
  section <- card_section(text, "Implementation (SBC)", "Recovery")
  expect_match(section, "excluded", fixed = TRUE)
  expect_no_match(section, "prior:", fixed = TRUE)
  expect_no_match(section, "A calibration speaks for", fixed = TRUE)
  expect_no_match(section, plan$sbc$prior$prior[[1L]], fixed = TRUE)
  # the prior check drew and was judged under that prior, so its section
  # lists the rows and the reason, and points nowhere else
  pp <- card_section(text, "Prior predictive", "Implementation (SBC)")
  expect_match(
    pp, "judged on the draws under the calibration's prior",
    fixed = TRUE
  )
  expect_no_match(pp, "listed in", fixed = TRUE)
  expect_match(pp, "The calibration's prior:", fixed = TRUE)
  expect_match(pp, plan$sbc$prior_reason, fixed = TRUE)
  for (row in plan$sbc$prior$prior) {
    expect_match(pp, row, fixed = TRUE)
    # once on the card
    expect_length(gregexpr(row, text, fixed = TRUE)[[1L]], 1L)
  }
  expect_identical(card_unexplained_numbers(out, result), character())
})

test_that("a failed or skipped calibration still names its prior", {
  skip_card_render()
  for (status in c("failed", "skipped")) {
    card <- rerender(function(result) {
      at <- result$stages$stage == "sbc"
      result$stages$status[at] <- status
      result$stages$message[at] <- "mock sampler failure"
      result
    })
    section <- card_section(card$text, "Implementation (SBC)", "Recovery")
    expect_match(section, status, fixed = TRUE)
    expect_match(section, "The calibration's prior:", fixed = TRUE)
    for (row in card$result$plan$sbc$prior$prior) {
      expect_match(section, row, fixed = TRUE)
    }
    expect_match(section, "A calibration speaks for", fixed = TRUE)
  }
})

test_that("a prior check that was not judged does not say it was", {
  skip_card_render()
  card <- rerender(function(result) {
    at <- result$verdicts$stage == "prior_predictive"
    result$verdicts$pass[at] <- NA
    result$verdicts$value[at] <- NA_real_
    result$verdicts$note[at] <- "not judged: no floor or ceiling is known"
    result
  })
  pp <- card_section(card$text, "Prior predictive", "Implementation (SBC)")
  expect_no_match(pp, "judged on the draws under", fixed = TRUE)
  expect_match(pp, "The draws were not judged", fixed = TRUE)
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

test_that("the data requirements section prints the read-off and its bar", {
  skip_card_render()
  dir <- withr::local_tempdir()
  result <- mock_validation(
    readoff_validation_plan(), dir, validation_mock_fitter()
  )
  out <- file.path(dir, "card.html")
  render_card(result, out)
  text <- card_visible_text(out)
  # the title's last occurrence is the section's heading (the contents
  # and the verdicts table come first); the next stage ends it
  pieces <- strsplit(text, "Data requirements", fixed = TRUE)[[1L]]
  section <- sub(" Agreement .*", "", pieces[[length(pieces)]])
  expect_lt(nchar(section), nchar(text) / 2)
  # D83: the bar is stated as a convention in the section itself
  expect_match(
    section, "The read-off at r_low 0.8 is a convention, not a finding",
    fixed = TRUE
  )
  expect_match(section, "n_trials_at_r_low_0.7", fixed = TRUE)
  expect_match(section, "n_trials_at_r_low_0.9", fixed = TRUE)
  # the all-terms row: 10, 20 and 40 trials
  expect_match(section, "8 subjects\\s+subject\\s+all\\s+10\\s+20\\s+40")
  expect_identical(card_unexplained_numbers(out, result), character())
})

test_that("a bar no design reaches keeps its column of dashes", {
  skip_card_render()
  dir <- withr::local_tempdir()
  result <- mock_validation(
    readoff_validation_plan(criteria = list(read_off_r_low = 0.95)), dir,
    validation_mock_fitter()
  )
  out <- file.path(dir, "card.html")
  render_card(result, out)
  text <- card_visible_text(out)
  pieces <- strsplit(text, "Data requirements", fixed = TRUE)[[1L]]
  section <- sub(" Agreement .*", "", pieces[[length(pieces)]])
  expect_match(section, "n_trials_at_r_low_0.95", fixed = TRUE)
  expect_match(section, "subject\\s+all\\s+10\\s+40\\s+\u2013")
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

test_that("the numbers of the calibration's prior rows are explained", {
  # the card prints the rows; their numbers are typed by the plan, not by
  # the template
  file <- withr::local_tempfile(fileext = ".html")
  writeLines("<html><body><p>normal(7.31, 0.0613)</p></body></html>", file)
  result <- list(
    plan = list(sbc = list(prior = data.frame(prior = "normal(7.31, 0.0613)"))),
    headline = list(text = "")
  )
  expect_identical(card_unexplained_numbers(file, result), character())
  result$plan$sbc$prior <- NULL
  expect_setequal(card_unexplained_numbers(file, result), c("7.31", "0.0613"))
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
