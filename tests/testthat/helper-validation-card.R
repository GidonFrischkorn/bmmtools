# The typed-number check of spec 26 (render_card(), 26.2 and 26.3): every
# number a rendered card prints must come from the result it was rendered
# from. 26.3 sources this file to check the article cards.

#' The text a reader sees: no head, scripts, styles or tags
#' @noRd
card_visible_text <- function(file) {
  html <- paste(readLines(file, warn = FALSE), collapse = "\n")
  html <- gsub("(?s)<head>.*?</head>", " ", html, perl = TRUE)
  html <- gsub("(?s)<script.*?</script>", " ", html, perl = TRUE)
  html <- gsub("(?s)<style.*?</style>", " ", html, perl = TRUE)
  html <- gsub("<[^>]+>", " ", html)
  entities <- c(
    "&gt;" = ">", "&lt;" = "<", "&quot;" = "\"", "&#39;" = "'", "&amp;" = "&"
  )
  for (entity in names(entities)) {
    html <- gsub(entity, entities[[entity]], html, fixed = TRUE)
  }
  html <- gsub("&(#[0-9]+|#x[0-9a-fA-F]+|[a-zA-Z]+);", " ", html)
  # Pandoc's smart quotes
  html <- gsub("\u2018|\u2019", "'", html)
  gsub("[[:space:]]+", " ", html)
}

#' The numbers in a text, leaving out those inside names, versions,
#' issue references and decision numbers
#' @noRd
card_number_tokens <- function(text) {
  pattern <- "(?<![\\w.#-])-?[0-9]+(?:\\.[0-9]+)?(?:[eE][-+]?[0-9]+)?(?![\\w.])"
  unlist(regmatches(text, gregexpr(pattern, text, perl = TRUE)))
}

#' Whether a printed token is one of `values` rounded to its precision
#' @noRd
card_token_explained <- function(token, values) {
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return(FALSE)
  }
  t <- as.numeric(token)
  parts <- strsplit(tolower(token), "e", fixed = TRUE)[[1L]]
  decimals <- if (grepl(".", parts[[1L]], fixed = TRUE)) {
    nchar(sub("^[^.]*\\.", "", parts[[1L]]))
  } else {
    0L
  }
  exponent <- if (length(parts) > 1L) as.numeric(parts[[2L]]) else 0
  # an integer's trailing zeros may be rounding: 22000 is 21987 to three
  # significant digits
  if (decimals == 0L && exponent == 0 && t != 0) {
    digits <- sub("^-", "", parts[[1L]])
    exponent <- nchar(digits) - nchar(sub("0+$", "", digits))
  }
  tol <- 0.5 * 10^(exponent - decimals) * (1 + 1e-9)
  any(abs(values - t) <= tol)
}

#' The numeric values of a list, recursively
#' @noRd
card_numeric_values <- function(x) {
  if (is.data.frame(x)) x <- as.list(x)
  if (is.list(x)) {
    return(unlist(lapply(x, card_numeric_values), use.names = FALSE))
  }
  if (is.numeric(x) || is.logical(x)) as.double(x) else numeric()
}

#' The strings of a list, recursively, with dates and times formatted
#' @noRd
card_strings <- function(x) {
  if (inherits(x, c("POSIXt", "Date"))) {
    return(format(x))
  }
  if (is.data.frame(x)) x <- as.list(x)
  if (is.list(x)) {
    return(unlist(lapply(x, card_strings), use.names = FALSE))
  }
  if (is.character(x) || is.factor(x)) as.character(x) else character()
}

#' The numbers on a rendered card that the result does not explain
#'
#' A number is explained when it rounds a stored value (the numbers, the
#' verdicts, the stage times, the run record and the plan's scalar
#' fields) or appears in a string the result stores (versions, dates,
#' criteria, notes, design labels, the calibration's prior rows).
#' @noRd
card_unexplained_numbers <- function(file, result) {
  plan <- result$plan
  changes <- result$changes
  values <- c(
    card_numeric_values(result$numbers[c("value", "mcse")]),
    card_numeric_values(result$verdicts["value"]),
    card_numeric_values(result$stages[c("seconds", "issue")]),
    card_numeric_values(result$run_info),
    card_numeric_values(plan[c(
      "pars", "sds", "cors", "grid", "reps", "seed", "seeds", "sampler",
      "sbc", "prior_check", "criteria"
    )]),
    if (!is.null(changes)) {
      c(
        card_numeric_values(
          changes[c("old", "new", "diff", "mcse_old", "mcse_new")]
        ),
        attr(changes, "n_compared"), sum(changes$moved, na.rm = TRUE)
      )
    }
  )
  strings <- c(
    card_strings(result$run_info),
    card_strings(plan[c(
      "model_name", "defaults", "owner", "version", "hash",
      "bmmtools_version", "sections"
    )]),
    plan$sbc$prior_reason, plan$prior_check$range_reason,
    # the prior rows the card prints: "normal(0, 2)" types 0 and 2
    card_strings(plan$sbc$prior),
    card_strings(result$stages[c("title", "status", "reason", "message")]),
    card_strings(result$verdicts[c("title", "criterion", "note")]),
    result$headline$text,
    card_strings(result$numbers[c("condition", "level", "term", "metric")]),
    if (!is.null(changes)) {
      c(attr(changes, "old_bmm"), attr(changes, "new_bmm"))
    }
  )
  known <- card_number_tokens(paste(strings, collapse = " "))
  tokens <- unique(card_number_tokens(card_visible_text(file)))
  tokens <- setdiff(tokens, known)
  explained <- vapply(tokens, card_token_explained, logical(1), values)
  tokens[!explained]
}

#' The prose of a Quarto template: no code chunks, no inline code, no
#' front-matter fences, no blank lines
#' @noRd
card_template_prose <- function(file) {
  lines <- readLines(file, warn = FALSE)
  fence <- grepl("^```", lines)
  inside <- (cumsum(fence) %% 2L == 1L) | fence
  lines <- lines[!inside]
  lines <- gsub("`r [^`]*`", "", lines)
  lines[nzchar(trimws(lines)) & lines != "---"]
}
