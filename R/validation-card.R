# render_card(): the validation card of a stored result (Milestone 26,
# D82). The template, inst/card/card.qmd, reads the stored result and
# nothing else: Quarto renders it in a fresh R session, which loads the
# installed bmmtools rather than the one in memory, so the template calls
# no bmmtools function and every number it prints is read off the result.

#' Render the validation card of a result
#'
#' `r lifecycle::badge("experimental")`
#'
#' Renders the card of a [run_validation()] result to one self-contained
#' HTML file: the headline verdict, the run record and the plan, the
#' verdict of every stage beside its criterion, a section per stage in
#' the order of the workflow (results, or why the stage did not run), the
#' model's own sections, and the numbers that moved since the earlier
#' result the run was compared with.
#'
#' The template prints no number of its own: every figure on the card is
#' read off the result (its `numbers`, `verdicts`, `stages`, `run_info`
#' and `plan`), so a card always describes the run it was rendered from.
#'
#' @param result A `bmmtools_validation`, its `result.rds`, or the
#'   directory it was run in.
#' @param output The HTML file to write. Its directory is created if
#'   missing, and an existing file is replaced.
#' @param format `"html"`, the only format so far.
#' @param template `NULL` for the card bmmtools ships
#'   (`system.file("card", "card.qmd", package = "bmmtools")`), or the
#'   path to a Quarto document of your own that takes the parameter
#'   `result`, the path to a `result.rds`.
#' @param quiet Passed to [quarto::quarto_render()].
#'
#' @return `output`, invisibly.
#'
#' @section Requirements:
#' The quarto R package and the Quarto CLI (<https://quarto.org>); both
#' are checked before anything is rendered.
#'
#' @seealso [run_validation()], [verdict_table()], [card_numbers()],
#'   [compare_validations()].
#'
#' @examples
#' \dontrun{
#' plan <- validation_plan("ezdm", seed = 2026)
#' result <- run_validation(plan, "validation/ezdm")
#' render_card(result, "validation/ezdm/card.html")
#' }
#'
#' @export
render_card <- function(result, output, format = "html", template = NULL,
                        quiet = TRUE) {
  if (is.null(result)) {
    cli::cli_abort(
      "{.arg result} must be a validation result or a path to one, not \\
       {.code NULL}."
    )
  }
  result <- read_previous_validation(result, arg = "result")
  single <- is.character(output) && length(output) == 1L
  if (!single || is.na(output) || !nzchar(output)) {
    cli::cli_abort(
      "{.arg output} must be a single file path, not \\
       {.obj_type_friendly {output}}."
    )
  }
  format <- rlang::arg_match(format, "html")
  template <- template %||% system.file("card", "card.qmd",
    package = "bmmtools"
  )
  single <- is.character(template) && length(template) == 1L
  if (!single || !file.exists(template)) {
    cli::cli_abort("The template {.file {template}} does not exist.")
  }
  rlang::check_installed("quarto", reason = "to render a validation card.")
  if (is.null(quarto_cli_path())) {
    cli::cli_abort(c(
      "Rendering a card needs the Quarto CLI, which was not found.",
      i = "Install it from {.url https://quarto.org}; the quarto R \\
           package is installed."
    ))
  }

  # Quarto writes beside its input, which for the shipped template is the
  # installed package: render a copy in a directory of its own
  work <- tempfile("bmmtools-card-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "card.qmd")
  file.copy(template, input)
  path <- file.path(work, "result.rds")
  saveRDS(result, path)
  quarto::quarto_render(
    input,
    execute_params = list(result = path), quiet = quiet, as_job = FALSE
  )
  rendered <- file.path(work, "card.html")
  if (!file.exists(rendered)) {
    cli::cli_abort("Quarto wrote no {.file card.html}.")
  }
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  # file.copy() warns with the reason before it returns FALSE
  copied <- tryCatch(
    file.copy(rendered, output, overwrite = TRUE),
    warning = function(w) structure(FALSE, reason = conditionMessage(w))
  )
  if (!isTRUE(copied)) {
    cli::cli_abort(c(
      "Could not write {.file {output}}.",
      x = attr(copied, "reason")
    ))
  }
  invisible(output)
}

#' The Quarto CLI, or `NULL`; a seam for the tests
#' @noRd
quarto_cli_path <- function() {
  quarto::quarto_path()
}
