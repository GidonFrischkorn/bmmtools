# Correlated truths.
#
# Subject values of the varying parameters, and any observed covariates,
# are one draw from a multivariate normal. The draw is `Z %*% chol(R)`
# over the varying parameters first and the covariates after them, scaled
# afterwards. Two measured properties follow from that order (2026-09-14):
# with an identity `R` the draw equals Milestone 3's per-parameter
# `rnorm()` bit for bit, because `rnorm(n, 0, 0)` consumes no random
# numbers and only varying terms get a column of `Z`; and because the
# Cholesky factor is upper triangular, a dimension appended at the end
# never changes the columns before it.

#' Validate a correlation matrix's shape and entries
#'
#' Symmetry and a unit diagonal are checked to `1e-8`, so a matrix built
#' by arithmetic (as `cors_from_factors()` does) is accepted.
#'
#' @noRd
check_correlation_matrix <- function(x, arg, call = rlang::caller_env()) {
  if (!is.matrix(x) || !is.numeric(x)) {
    cli::cli_abort(
      "{.arg {arg}} must be a numeric matrix, not {.obj_type_friendly {x}}.",
      call = call
    )
  }
  if (nrow(x) != ncol(x)) {
    cli::cli_abort(
      "{.arg {arg}} must be square, not {nrow(x)} by {ncol(x)}.",
      call = call
    )
  }
  if (anyNA(x)) {
    cli::cli_abort("{.arg {arg}} has missing values.", call = call)
  }
  if (max(abs(x - t(x))) > 1e-8) {
    cli::cli_abort("{.arg {arg}} must be symmetric.", call = call)
  }
  if (any(abs(diag(x) - 1) > 1e-8)) {
    cli::cli_abort("{.arg {arg}} must have a diagonal of 1.", call = call)
  }
  if (any(abs(x) > 1 + 1e-8)) {
    cli::cli_abort(
      "Every entry of {.arg {arg}} must be between -1 and 1.",
      call = call
    )
  }
  invisible(x)
}

#' Error unless a correlation matrix is positive definite
#'
#' Tested on the eigenvalues rather than by whether `chol()` happens to
#' succeed, because a singular matrix can pass `chol()` through rounding.
#'
#' @noRd
check_positive_definite <- function(x, arg, call = rlang::caller_env()) {
  if (nrow(x) == 0L) {
    return(invisible(x))
  }
  smallest <- min(eigen(x, symmetric = TRUE, only.values = TRUE)$values)
  if (smallest < sqrt(.Machine$double.eps)) {
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    shown <- signif(smallest, 3)
    cli::cli_abort(
      c(
        "{.arg {arg}} must be positive definite.",
        i = "Its smallest eigenvalue is {shown}; some correlations \\
             contradict each other."
      ),
      call = call
    )
  }
  invisible(x)
}

#' Validate `cors` and expand it to the full matrix the draw uses
#'
#' The full matrix is ordered as the draw is: varying parameters in `pars`
#' order, then covariates in list order. Terms `cors` does not name are
#' uncorrelated. A parameter that does not vary may be named, so that a
#' grid row which switches one SD off can share a matrix with the others,
#' but only with zero correlations.
#'
#' @param cors `NULL` or a correlation matrix with dimnames.
#' @param pars,sds Validated population values and SDs.
#' @param covariates Validated covariates or `NULL`.
#' @return The full correlation matrix, possibly 0 by 0.
#' @noRd
check_cors <- function(cors, pars, sds, covariates,
                       call = rlang::caller_env()) {
  varying <- names(sds)[sds > 0]
  dims <- c(varying, names(covariates))
  full <- diag(length(dims))
  dimnames(full) <- list(dims, dims)
  if (is.null(cors)) {
    return(full)
  }

  check_correlation_matrix(cors, "cors", call = call)
  given <- rownames(cors)
  if (is.null(given) || !identical(given, colnames(cors))) {
    cli::cli_abort(
      c(
        "{.arg cors} must have the same dimnames on both margins.",
        i = "They name the parameters and covariates it correlates."
      ),
      call = call
    )
  }
  if (anyDuplicated(given) > 0L) {
    cli::cli_abort(
      "{.arg cors} names {.val {given[duplicated(given)]}} twice.",
      call = call
    )
  }
  unknown <- setdiff(given, c(names(pars), names(covariates)))
  if (length(unknown) > 0L) {
    cli::cli_abort(
      c(
        "{.arg cors} names {.val {unknown}}, which is neither a parameter \\
         in {.arg pars} nor a covariate.",
        i = "Known: {.val {c(names(pars), names(covariates))}}."
      ),
      call = call
    )
  }

  fixed_here <- setdiff(given, dims)
  off_diagonal <- cors
  diag(off_diagonal) <- 0
  constant_cor <- fixed_here[
    vapply(fixed_here, function(p) any(off_diagonal[p, ] != 0), logical(1))
  ]
  if (length(constant_cor) > 0L) {
    cli::cli_abort(
      c(
        "{.arg cors} correlates {.val {constant_cor}}, which \\
         {?does/do} not vary.",
        i = "{cli::qty(length(constant_cor))}Give {?it/them} a positive \\
             entry in {.arg sds}, or set {?its/their} correlations to 0."
      ),
      call = call
    )
  }

  keep <- intersect(given, dims)
  full[keep, keep] <- cors[keep, keep]
  check_positive_definite(full, "cors", call = call)
  full
}

#' Draw subject values, and covariate values, on the link scale
#'
#' Two steps, so that [simulate_recovery()] can run the generator between
#' them: the parameter part of `Z` first, the covariate part after. The
#' covariate columns of `Z` are independent of the parameter columns, so
#' drawing them later changes nothing about the joint distribution, and
#' the generated responses then do not depend on whether covariates exist.
#'
#' @param cors The full matrix from `check_cors()`.
#' @return A matrix, subjects by parameters then covariates.
#' @noRd
draw_correlated_pars <- function(pars, sds, cors, covariates, n_subjects) {
  parts <- draw_parameter_values(pars, sds, cors, n_subjects)
  draw_covariate_values(parts, cors, covariates)
}

#' The parameter step: values for every parameter, and its part of `Z`
#' @noRd
draw_parameter_values <- function(pars, sds, cors, n_subjects) {
  ids <- as.character(seq_len(n_subjects))
  varying <- names(sds)[sds > 0]
  values <- matrix(
    rep(pars, each = n_subjects),
    nrow = n_subjects,
    dimnames = list(ids, names(pars))
  )
  z <- matrix(stats::rnorm(n_subjects * length(varying)), nrow = n_subjects)
  if (length(varying) == 0L) {
    return(list(values = values, z = z))
  }
  block <- cors[varying, varying, drop = FALSE]
  # the identity skips the product so that the Milestone 3 draw is kept
  # exactly, whatever the BLAS does with a multiplication by one; the
  # leading block of chol(cors) is chol() of the leading block
  x <- if (all(block[upper.tri(block)] == 0)) z else z %*% chol(block)
  for (j in seq_along(varying)) {
    p <- varying[[j]]
    values[, p] <- pars[[p]] + sds[[p]] * x[, j]
  }
  list(values = values, z = z)
}

#' The covariate step: append covariate columns to the parameter values
#' @noRd
draw_covariate_values <- function(parts, cors, covariates) {
  values <- parts$values
  if (length(covariates) == 0L) {
    return(values)
  }
  n_subjects <- nrow(values)
  cov_names <- names(covariates)
  z_cov <- matrix(
    stats::rnorm(n_subjects * length(cov_names)),
    nrow = n_subjects
  )
  # check_cors() orders the full matrix as the varying parameters, whose
  # part of Z was drawn first, then the covariates
  upper <- chol(cors)
  cov_columns <- ncol(parts$z) + seq_along(cov_names)
  x <- cbind(parts$z, z_cov) %*% upper[, cov_columns, drop = FALSE]
  out <- matrix(
    NA_real_,
    nrow = n_subjects, ncol = length(cov_names),
    dimnames = list(rownames(values), cov_names)
  )
  for (j in seq_along(cov_names)) {
    g <- cov_names[[j]]
    out[, g] <- covariates[[g]][["mean"]] + covariates[[g]][["sd"]] * x[, j]
  }
  cbind(values, out)
}

#' Validate the covariates
#'
#' Names become data columns and truth terms, so they must be syntactic,
#' carry no `_` (the separator of D27 terms) and collide with nothing the
#' model or the data already use.
#'
#' @return A named list of `c(mean =, sd =)`, or `NULL`.
#' @noRd
check_covariates <- function(covariates, model, call = rlang::caller_env()) {
  if (is.null(covariates)) {
    return(NULL)
  }
  nms <- names(covariates)
  bad <- !is.list(covariates) || is.data.frame(covariates) ||
    length(covariates) == 0L || is.null(nms) || anyNA(nms) ||
    !all(nzchar(nms))
  if (bad) {
    cli::cli_abort(
      "{.arg covariates} must be a named list of {.code c(mean = , sd = )}, \\
       or {.code NULL}.",
      call = call
    )
  }
  bad_shape <- nms[!vapply(covariates, function(g) {
    is.numeric(g) && length(g) == 2L && setequal(names(g), c("mean", "sd")) &&
      all(is.finite(g))
  }, logical(1))]
  if (length(bad_shape) > 0L) {
    cli::cli_abort(
      "{.arg covariates} entr{?y/ies} {.val {bad_shape}} must be \\
       {.code c(mean = , sd = )} with finite values.",
      call = call
    )
  }
  covariates <- lapply(covariates, function(g) {
    c(mean = g[["mean"]], sd = g[["sd"]])
  })
  no_spread <- nms[vapply(covariates, function(g) g[["sd"]] <= 0, logical(1))]
  if (length(no_spread) > 0L) {
    cli::cli_abort(
      "{.arg covariates} {.val {no_spread}} must have {.code sd > 0}.",
      call = call
    )
  }

  info <- model_parameters(model)
  columns <- Filter(is.character, c(model$resp_vars, model$other_vars))
  taken <- c("id", info$free, info$fixed, unlist(columns, use.names = FALSE))
  bad_name <- nms[
    make.names(nms) != nms | grepl("_", nms, fixed = TRUE) | nms %in% taken |
      duplicated(nms)
  ]
  if (length(bad_name) > 0L) {
    cli::cli_abort(
      c(
        "{.arg covariates} name{?s} {.val {unique(bad_name)}} cannot be used.",
        i = "A covariate name must be syntactic, unique, contain no \\
             {.code _}, and differ from {.val id}, the model's parameters \\
             and its data columns."
      ),
      call = call
    )
  }
  covariates
}

#' Evaluate a truth argument that may be a zero-argument function
#'
#' Called inside the seeded block, in the order `pars`, `sds`, `cors`, so
#' that random hyperparameters continue the random number stream rather
#' than reuse the part the subject draws will take.
#'
#' @noRd
resolve_truth_arg <- function(x, arg, want = c("vector", "matrix"),
                              call = rlang::caller_env()) {
  want <- rlang::arg_match(want)
  if (!is.function(x)) {
    return(x)
  }
  if (length(formals(x)) > 0L) {
    cli::cli_abort(
      c(
        "{.arg {arg}} is a function, so it must take no arguments.",
        i = "In {.fn recovery_grid}, {.code function(row)} is allowed; the \\
             grid calls it with the row."
      ),
      call = call
    )
  }
  out <- x()
  ok <- if (identical(want, "vector")) {
    is.numeric(out) && !is.matrix(out)
  } else {
    is.null(out) || (is.matrix(out) && is.numeric(out))
  }
  if (!ok) {
    # nolint next: object_usage_linter. Used by cli's glue interpolation.
    expected <- if (identical(want, "vector")) {
      "a named numeric vector"
    } else {
      "a numeric matrix or NULL"
    }
    cli::cli_abort(
      "{.arg {arg}} is a function that returned \\
       {.obj_type_friendly {out}}; it must return {expected}.",
      call = call
    )
  }
  out
}

#' The truth table of SDs: varying parameters only
#' @noRd
sd_table <- function(sds) {
  varying <- names(sds)[sds > 0]
  tibble::tibble(term = varying, true_value = as.double(unname(sds[varying])))
}

#' A pair term: the two names in C-locale order
#' @noRd
pair_term <- function(a, b) {
  sorted <- sort(c(a, b), method = "radix")
  list(
    term = paste(sorted, collapse = "__"),
    var1 = sorted[[1L]],
    var2 = sorted[[2L]]
  )
}

#' The truth table of correlations: every pair of the full matrix
#' @noRd
cor_table <- function(cors) {
  empty <- tibble::tibble(
    term = character(), var1 = character(), var2 = character(),
    true_value = double()
  )
  if (is.null(cors) || nrow(cors) < 2L) {
    return(empty)
  }
  dims <- rownames(cors)
  combos <- utils::combn(length(dims), 2L)
  pieces <- lapply(seq_len(ncol(combos)), function(k) {
    i <- combos[1L, k]
    j <- combos[2L, k]
    pair <- pair_term(dims[[i]], dims[[j]])
    tibble::tibble(
      term = pair$term, var1 = pair$var1, var2 = pair$var2,
      true_value = as.double(cors[i, j])
    )
  })
  dplyr::bind_rows(empty, pieces)
}

#' Build a correlation matrix from factor loadings
#'
#' The model-implied correlation matrix of a factor model, for the `cors`
#' argument of [simulate_recovery()] and [recovery_grid()]: tasks that
#' load on a common ability, parameters that share a latent source, or a
#' covariate that measures the factor. The result is `L Φ Lᵀ` with its
#' diagonal set to 1, so each term's unique variance is one minus its
#' communality.
#'
#' This is the **population** matrix. A simulation that instead takes the
#' empirical correlation of a finite sample of factor scores (as the
#' indSim scripts do with 1000 lavaan draws) adds sampling noise to the
#' generating values.
#'
#' @param loadings A numeric matrix, terms by factors, with the terms as
#'   row names; or a named numeric vector for a single factor.
#' @param factor_cors The factor correlation matrix. `NULL` means
#'   uncorrelated factors. When both it and `loadings` carry factor names,
#'   they must agree.
#'
#' @return A correlation matrix with the terms as dimnames.
#'
#' @examples
#' # two tasks per ability, the abilities correlated .5
#' loadings <- matrix(
#'   c(0.8, 0.7, 0, 0, 0, 0, 0.75, 0.85),
#'   ncol = 2,
#'   dimnames = list(
#'     c("kappa_task1", "kappa_task2", "thetat_task1", "thetat_task2"),
#'     c("precision", "memory")
#'   )
#' )
#' phi <- matrix(c(1, 0.5, 0.5, 1), 2,
#'   dimnames = rep(list(c("precision", "memory")), 2)
#' )
#' cors_from_factors(loadings, phi)
#'
#' # one factor
#' cors_from_factors(c(a = 0.8, b = 0.6, G = 0.5))
#'
#' @export
cors_from_factors <- function(loadings, factor_cors = NULL) {
  if (!is.numeric(loadings)) {
    cli::cli_abort(
      "{.arg loadings} must be numeric, not {.obj_type_friendly {loadings}}."
    )
  }
  if (!is.matrix(loadings)) {
    terms <- names(loadings)
    loadings <- matrix(loadings, ncol = 1L, dimnames = list(terms, NULL))
  }
  terms <- rownames(loadings)
  bad_terms <- is.null(terms) || anyNA(terms) || !all(nzchar(terms)) ||
    anyDuplicated(terms) > 0L
  if (bad_terms) {
    cli::cli_abort(
      "{.arg loadings} must carry unique term names (row names, or the \\
       names of a vector)."
    )
  }
  if (anyNA(loadings) || !all(is.finite(loadings))) {
    cli::cli_abort("{.arg loadings} has missing or infinite values.")
  }

  n_factors <- ncol(loadings)
  if (is.null(factor_cors)) {
    factor_cors <- diag(n_factors)
  } else {
    check_correlation_matrix(factor_cors, "factor_cors")
    check_positive_definite(factor_cors, "factor_cors")
    if (nrow(factor_cors) != n_factors) {
      cli::cli_abort(
        "{.arg factor_cors} must be {n_factors} by {n_factors}, one row \\
         per column of {.arg loadings}."
      )
    }
    factors <- colnames(loadings)
    given <- rownames(factor_cors)
    if (!is.null(factors) && !is.null(given) && !identical(factors, given)) {
      cli::cli_abort(
        c(
          "The factor names of {.arg factor_cors} and {.arg loadings} differ.",
          i = "{.arg loadings}: {.val {factors}}; {.arg factor_cors}: \\
               {.val {given}}."
        )
      )
    }
  }

  implied <- loadings %*% factor_cors %*% t(loadings)
  communality <- diag(implied)
  too_high <- terms[communality > 1 + 1e-8]
  if (length(too_high) > 0L) {
    cli::cli_abort(
      c(
        "The communalit{?y/ies} of {.val {too_high}} exceed{?s/} 1.",
        i = "A term cannot share more variance with the factors than it has."
      )
    )
  }
  diag(implied) <- 1
  dimnames(implied) <- list(terms, terms)
  check_positive_definite(implied, "the implied correlation matrix")
  implied
}
