#' The link vocabulary
#'
#' The link names `inverse_link()` understands, in the order
#' `bmm:::link_transform()` lists them.
#'
#' @return A character vector of the twelve supported link names.
#' @noRd
link_vocabulary <- function() {
  c(
    "identity", "log", "softplus", "log1p", "logm1", "inverse",
    "sqrt", "logit", "probit", "tan_half", "loglog", "cloglog"
  )
}

#' Transform values from the link scale to the natural scale
#'
#' Cognitive measurement models are estimated on a link scale that keeps
#' parameters unconstrained --- a memory precision is sampled as
#' `log(kappa)`, a mixture weight as `logit(theta)` --- while readers
#' interpret them on the natural scale. `inverse_link()` is the
#' back-transform, applied by [recover()] to estimates, interval bounds
#' and generating values alike so that both sides of a recovery
#' comparison are on the same scale.
#'
#' The vocabulary and the closed forms are taken from bmm's internal
#' `link_transform()` (`bmm/R/helpers-parameters.R`), so the two cannot
#' disagree about what a link name means. bmmtools ships its own copy
#' because `link_transform()` is not exported.
#'
#' @param x A numeric vector of values on the link scale. `NA` is passed
#'   through at its own position.
#' @param link A single string naming the link: one of `"identity"`,
#'   `"log"`, `"softplus"`, `"log1p"`, `"logm1"`, `"inverse"`, `"sqrt"`,
#'   `"logit"`, `"probit"`, `"tan_half"`, `"loglog"` or `"cloglog"`.
#'   `NULL` is a synonym for `"identity"`.
#'
#' @return A double vector the same length as `x`, on the natural scale.
#'
#' @details
#' Values outside a link's natural domain are transformed as R computes
#' them rather than being replaced by `NA`: `inverse_link(0, "inverse")`
#' is `Inf` and `inverse_link(-2, "sqrt")` is `4`. In a scoring pipeline
#' a visible `Inf` is safer than a silent `NA`, which would be absorbed
#' by a metric's missing-data guard and reported as a smaller sample
#' rather than as a problem.
#'
#' @examples
#' inverse_link(c(-1, 0, 1), "log")
#' inverse_link(0, "logit")
#'
#' # a precision sampled on the log scale, back on its own scale
#' inverse_link(c(1.6, 2.3), "log")
#'
#' @export
inverse_link <- function(x, link) {
  if (is.null(link)) link <- "identity"

  if (!is.character(link) || length(link) != 1L || is.na(link)) {
    cli::cli_abort(
      "{.arg link} must be a single string, not {.obj_type_friendly {link}}."
    )
  }
  if (!link %in% link_vocabulary()) {
    cli::cli_abort(c(
      "{.arg link} must be one of {.val {link_vocabulary()}}, \\
       not {.val {link}}.",
      i = "The vocabulary follows bmm's own link table."
    ))
  }
  # logical is excluded deliberately: TRUE/FALSE reaching a link
  # transform means something upstream went wrong, and coercing it to
  # 1/0 would hide that.
  if (!is.numeric(x)) {
    cli::cli_abort(
      "{.arg x} must be a numeric vector, not {.obj_type_friendly {x}}."
    )
  }

  x <- as.double(x)

  switch(link,
    identity = x,
    log = exp(x),
    softplus = log1p(exp(x)),
    log1p = expm1(x),
    logm1 = exp(x) + 1,
    inverse = 1 / x,
    sqrt = x^2,
    logit = stats::plogis(x),
    probit = stats::pnorm(x),
    tan_half = 2 * atan(x),
    loglog = exp(-exp(x)),
    cloglog = -expm1(-exp(x))
  )
}
