# The contrast design: what the fit estimates when the data are cell means.
#
# Data generation never changes (D44). Subject values are drawn around cell
# means exactly as they always were, so `subjects = "fixed"`, the
# correlation machinery and the RNG stream all keep their behaviour. What a
# contrast design changes is the *design the fit sees* --- the task column
# carries a contrast matrix, so brms estimates an intercept and k - 1
# contrasts --- and therefore the truth those coefficients must be scored
# against. This file is that transform, and nothing else.

#' Resolve and validate the `contrasts` argument
#'
#' @param contrasts A `k` by `k - 1` matrix, or a function `(n) -> matrix`
#'   such as [stats::contr.treatment] or
#'   `bayestestR::contr.equalprior`. `NULL` is treatment coding.
#' @param k The number of task levels.
#' @return A `k` by `k - 1` matrix.
#' @noRd
check_contrasts <- function(contrasts, k, call = rlang::caller_env()) {
  if (is.null(contrasts)) {
    contrasts <- stats::contr.treatment
  }
  if (is.function(contrasts)) {
    contrasts <- tryCatch(
      contrasts(k),
      error = function(e) {
        cli::cli_abort(
          c(
            "{.arg contrasts} failed when called with {.val {k}} levels.",
            x = conditionMessage(e)
          ),
          call = call
        )
      }
    )
  }
  if (!is.matrix(contrasts) || !is.numeric(contrasts)) {
    cli::cli_abort(
      "{.arg contrasts} must be a numeric matrix or a function returning \\
       one, not {.obj_type_friendly {contrasts}}.",
      call = call
    )
  }
  if (nrow(contrasts) != k || ncol(contrasts) != k - 1L) {
    cli::cli_abort(
      c(
        "{.arg contrasts} must have {k} row{?s} and {k - 1} column{?s}, \\
         not {nrow(contrasts)} and {ncol(contrasts)}.",
        i = "One row per task level and one column per contrast."
      ),
      call = call
    )
  }
  if (anyNA(contrasts)) {
    cli::cli_abort(
      "{.arg contrasts} must not contain {.val {NA}}.",
      call = call
    )
  }
  design <- cbind(1, contrasts)
  inverse <- tryCatch(solve(design), error = function(e) NULL)
  if (is.null(inverse)) {
    cli::cli_abort(
      c(
        "{.arg contrasts} does not give an invertible design.",
        i = "Its columns must be linearly independent of each other and \\
             of the intercept."
      ),
      call = call
    )
  }
  contrasts
}

#' The terms a contrast design has, in place of one parameter's cells
#'
#' `kappa` over two tasks of `task` becomes `kappa` (the intercept) and
#' `kappa_task1` (the contrast), which is what brms names the coefficients
#' of a factor carrying a contrast matrix: measured 2026-09-17, the
#' coefficient of a two-level factor is `<task_col>1` whatever the levels
#' are called.
#'
#' @noRd
contrast_terms <- function(par, k, task_col) {
  c(par, paste0(par, "_", task_col, seq_len(k - 1L)))
}

#' The linear map from cell means to intercept and contrasts
#'
#' With `M = cbind(1, C)` the design over the `k` task levels,
#' `mu = M beta`, so `beta = L mu` with `L = solve(M)`. The whole transform
#' is that one matrix: population values, subject values, SDs and
#' correlations all follow from it.
#'
#' @noRd
contrast_map <- function(contrasts) {
  solve(cbind(1, contrasts))
}

#' The block-diagonal transform over every term of a simulation
#'
#' Each parameter's `k` cell terms are replaced by its intercept and
#' contrasts through `L`; covariates, and any term that is not a task cell,
#' pass through unchanged.
#'
#' @param terms The cell terms, parameter-major, as [task_terms()] orders
#'   them, followed by any covariate names.
#' @return A matrix whose rows are the new terms and whose columns are
#'   `terms`.
#' @noRd
contrast_transform <- function(terms, pars, contrasts, tasks, task_col) {
  k <- length(tasks)
  map <- contrast_map(contrasts)
  cells <- task_terms(pars, tasks, task_col)
  blocks <- split(cells, rep(pars, each = k))[pars]
  new_terms <- unlist(lapply(pars, contrast_terms, k = k, task_col = task_col))
  passthrough <- setdiff(terms, cells)
  out <- matrix(
    0,
    nrow = length(new_terms) + length(passthrough),
    ncol = length(terms),
    dimnames = list(c(new_terms, passthrough), terms)
  )
  for (p in pars) {
    rownames_p <- contrast_terms(p, k, task_col)
    out[rownames_p, blocks[[p]]] <- map
  }
  for (g in passthrough) out[g, g] <- 1
  out
}

#' The truth of a contrast design, from the truth of the cells it generated
#'
#' A pure function: no fit is involved and no data are read. Population and
#' subject values become `beta = L mu`; the covariance of the subject values
#' becomes `T Sigma T'`, whose diagonal's square root is the SD truth and
#' whose off-diagonal, standardised, is the correlation truth --- both on
#' the link scale, as they already are.
#'
#' The transform is exact, not an approximation: the subject values *are* a
#' linear function of the cell values, so their covariance is the
#' transformed covariance and nothing is lost. A design whose contrast
#' column is not orthogonal to the intercept simply gets a nonzero
#' intercept-contrast correlation, which is the truth for that design.
#'
#' @param pars,sds,values,cors The generating truth over cell terms, as
#'   [simulate_recovery()] holds it: a named vector of population values, a
#'   named vector of SDs, a subjects-by-terms matrix, and the correlation
#'   matrix over the varying terms and the covariates.
#' @param contrasts The `k` by `k - 1` contrast matrix.
#' @param tasks,task_col The task levels, in the order the contrast matrix
#'   codes them, and the column that holds them.
#' @param covariates The covariate specification, whose names pass through
#'   the transform unchanged.
#'
#' @return A list with `pars`, `sds`, `values` and `cors` over the contrast
#'   terms, in the shape [truth_tables()] takes.
#' @noRd
contrast_truth <- function(pars, sds, values, cors, contrasts, tasks,
                           task_col, covariates = NULL) {
  parameters <- unique(sub(
    paste0("_", task_col, "(", paste(tasks, collapse = "|"), ")$"), "",
    names(pars)
  ))
  cells <- task_terms(parameters, tasks, task_col)
  cov_names <- names(covariates)
  terms <- c(cells, cov_names)
  transform <- contrast_transform(
    terms, parameters, contrasts, tasks, task_col
  )
  new_terms <- rownames(transform)

  population <- stats::setNames(rep(0, length(terms)), terms)
  known <- intersect(terms, names(pars))
  population[known] <- pars[known]
  new_pars <- drop(transform %*% population)
  names(new_pars) <- new_terms

  # the subject values are transformed as they stand; their covariance is
  # transformed with them, so the SD and correlation truths stay exactly
  # the covariance of the values a fit sees
  new_values <- values[, terms, drop = FALSE] %*% t(transform)
  dimnames(new_values) <- list(rownames(values), new_terms)

  scale <- c(sds[cells], stats::setNames(rep(1, length(cov_names)), cov_names))
  sigma <- outer(scale, scale) * full_correlation(cors, terms)
  new_sigma <- transform %*% sigma %*% t(transform)
  new_sds <- sqrt(pmax(diag(new_sigma), 0))
  names(new_sds) <- new_terms

  correlations <- correlation_from_covariance(
    new_sigma,
    keep = c(names(new_sds)[new_sds > 0], cov_names)
  )
  list(
    # an orthogonal contrast has a true correlation of exactly zero and a
    # balanced design a true effect of exactly zero; solve() reaches them
    # as 1e-17 or so, and a truth of 1e-17 is not zero to the summary,
    # which asks whether a truth is zero before it reports a sign
    pars = zapsmall(new_pars[new_terms], digits = 12L),
    sds = new_sds[setdiff(new_terms, cov_names)],
    values = new_values,
    cors = zapsmall(correlations, digits = 12L),
    # which of the terms are effects rather than values of a parameter:
    # the same split `coefficient_kinds()` reads off a fit
    effects = setdiff(new_terms, c(parameters, cov_names))
  )
}

#' The correlation matrix over `terms`, identity where it says nothing
#'
#' `cors` names only the varying terms and the covariates (see
#' [check_cors()]); a cell whose SD is zero is absent from it and
#' contributes a zero row and column to the covariance, which the identity
#' gives once it is scaled by an SD of zero.
#' @noRd
full_correlation <- function(cors, terms) {
  out <- diag(length(terms))
  dimnames(out) <- list(terms, terms)
  if (is.null(cors) || nrow(cors) < 1L) {
    return(out)
  }
  keep <- intersect(rownames(cors), terms)
  out[keep, keep] <- cors[keep, keep]
  out
}

#' A correlation matrix from a covariance, over the terms that vary
#'
#' [stats::cov2cor()] divides by the SDs, so a term that does not vary
#' gives `NaN`. [check_cors()] names only the varying terms and the
#' covariates, and [cor_table()] enumerates every pair of what it is given,
#' so the result is restricted to `keep` for the same reason: a constant
#' has no correlation to record, and recording one would put a row in the
#' truth that the cell-means truth never had.
#' @noRd
correlation_from_covariance <- function(sigma, keep) {
  keep <- intersect(rownames(sigma), keep)
  out <- diag(length(keep))
  dimnames(out) <- list(keep, keep)
  if (length(keep) > 0L) {
    out[keep, keep] <- stats::cov2cor(sigma[keep, keep, drop = FALSE])
  }
  out
}
