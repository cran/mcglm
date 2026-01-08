#' @title Wald Tests for Dispersion Components
#'
#' @description
#' Performs Wald chi-square tests for dispersion (covariance) parameters
#' by response variable in multivariate covariance generalized linear
#' models fitted with \code{\link{mcglm}}. This function is intended for
#' joint hypothesis testing of dispersion coefficients associated with
#' categorical covariates with more than two levels. It is not designed
#' for model comparison.
#'
#' @param object An object of class \code{"mcglm"}, typically the result
#'   of a call to \code{\link{mcglm}}.
#' @param idx_list A list of integer vectors indexing dispersion
#'   parameters to be jointly tested for each response.
#' @param names_list A list of character vectors with covariate names to
#'   be displayed in the output tables.
#' @param ... Currently not used.
#'
#' @return
#' The object is a list of data frames, one per response variable. Each
#' data frame contains the following columns:
#' \describe{
#'   \item{Covariate}{Name of the covariate associated with the dispersion
#'   parameters being tested.}
#'   \item{Chi.Square}{Wald chi-square test statistic.}
#'   \item{Df}{Degrees of freedom of the test.}
#'   \item{p.value}{P-value associated with the chi-square test.}
#' }
#'
#' @seealso \code{\link{mcglm}}, \code{\link{vcov}}, \code{\link{coef}}
#'
#' @export

mc_anova_disp <- function(object, idx_list, names_list, ...) {

  ## ---- Argument checks ----
  if (!inherits(object, "mcglm")) {
    stop("'object' must be an object of class 'mcglm'.", call. = FALSE)
  }

  n_resp <- length(object$mu_list)

  if (length(idx_list) != n_resp) {
    stop("'idx_list' must have one element per response.", call. = FALSE)
  }

  if (length(names_list) != n_resp) {
    stop("'names_list' must have one element per response.", call. = FALSE)
  }

  ## ---- Extract covariance matrix ----
  vv <- vcov(object)
  n_par <- nrow(vv)

  ## ---- Identify dispersion parameter blocks ----
  n_disp <- vapply(object$matrix_pred, length, integer(1))
  idx_disp <- rep(seq_len(n_resp), n_disp)

  n_beta <- n_par - length(idx_disp)
  idx_full <- c(rep(0L, n_beta), idx_disp)

  ## ---- Split covariance matrix and parameters by response ----
  vcov_list <- vector("list", n_resp)
  disp_list <- vector("list", n_resp)

  for (i in seq_len(n_resp)) {
    sel <- idx_full == i
    vcov_list[[i]] <- vv[sel, sel, drop = FALSE]
    disp_list[[i]] <- coef(object, type = "tau", response = i)$Estimates
  }

  ## ---- Wald tests ----
  out <- vector("list", n_resp)

  for (i in seq_len(n_resp)) {

    idx <- idx_list[[i]]
    nm  <- names_list[[i]]

    if (length(idx) != length(disp_list[[i]])) {
      stop("Length mismatch between 'idx_list' and dispersion parameters.",
           call. = FALSE)
    }

    ## Remove intercept if present
    if (!is.null(nm) && nm[1] == "(Intercept)") {
      idx <- idx[-1]
      nm  <- nm[-1]
      disp_list[[i]] <- disp_list[[i]][-1]
      vcov_list[[i]] <- vcov_list[[i]][-1, -1, drop = FALSE]
    }

    terms <- sort(unique(idx))
    res_i <- vector("list", length(terms))

    for (j in seq_along(terms)) {

      sel <- idx == terms[j]
      beta <- disp_list[[i]][sel]
      V    <- vcov_list[[i]][sel, sel, drop = FALSE]

      invV <- tryCatch(
        chol2inv(chol(V)),
        error = function(e) solve(V, tol = 1e-10)
      )

      X2 <- as.numeric(t(beta) %*% invV %*% beta)
      df <- length(beta)

      res_i[[j]] <- data.frame(
        Covariate  = nm[sel][1],
        Chi.Square = X2,
        Df         = df,
        p.value    = pchisq(X2, df, lower.tail = FALSE)
      )
    }

    out[[i]] <- do.call(rbind, res_i)
    rownames(out[[i]]) <- NULL
  }
  return(out)
}


