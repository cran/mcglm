#' @title MANOVA-Type Test for Dispersion Components of mcglm Models
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description
#' Performs a MANOVA-type Wald test for the dispersion parameters
#' of multivariate covariance generalized linear models fitted
#' using \code{\link{mcglm}}. The test is based on quadratic forms
#' of the estimated dispersion parameters and their covariance
#' matrix, yielding statistics analogous to the Hotelling–Lawley trace.
#'
#' @param object An object of class \code{"mcglm"}.
#' @param idx An integer vector defining the grouping structure
#'   of dispersion parameters to be tested.
#' @param effect_names A character vector with labels for the tested
#'   dispersion effects.
#' @param ... Further arguments (currently not used).
#'
#' @return
#' A data frame containing the MANOVA-type test results for the
#' dispersion parameters with the following columns:
#' \describe{
#'   \item{Effects}{Names of the tested dispersion effects.}
#'   \item{Df}{Degrees of freedom associated with each effect.}
#'   \item{Hotelling-Lawley}{Hotelling–Lawley trace statistic.}
#'   \item{Chi-square}{Chi-square test statistic.}
#'   \item{p-value}{P-values from the chi-square approximation.}
#' }
#'
#' @seealso \code{\link{mcglm}}, \code{\link{mc_manova}},
#'   \code{\link{coef.mcglm}}, \code{\link{vcov.mcglm}}
#'
#' @export

# manova for objects of mcglm class ----------------------------------
mc_manova_disp <- function(object, idx, effect_names, ...) {

  if (!inherits(object, "mcglm")) {
    stop("object must be of class 'mcglm'.", call. = FALSE)
  }

  disp <- coef(object, type = "tau")[, 1]
  n_disp <- length(disp)

  Parameters <- coef(object, type = "tau")$Parameters
  VCOV <- vcov(object)[Parameters, Parameters, drop = FALSE]

  if (length(idx) != n_disp) {
    stop("Length of 'idx' must match the number of dispersion parameters.",
         call. = FALSE)
  }

  FF <- mc_build_F(vector = idx)
  G <- Diagonal(length(object$mu_list), 1)
  CC <- lapply(FF, function(x) kronecker(G, x))

  N <- object$n_obs
  if (N <= 0) {
    stop("Number of observations must be positive.", call. = FALSE)
  }

  k <- length(CC)
  test_W <- numeric(k)
  df <- integer(k)
  p_value <- numeric(k)

  for (i in seq_along(CC)) {
    CVC <- CC[[i]] %*% VCOV %*% t(CC[[i]])

    CVC_inv <- tryCatch(
      solve(CVC),
      error = function(e) {
        stop("Singular covariance matrix in MANOVA dispersion test.",
             call. = FALSE)
      }
    )

    test_W[i] <- as.numeric(
      t(CC[[i]] %*% disp) %*% CVC_inv %*% (CC[[i]] %*% disp)
    )

    df[i] <- nrow(CC[[i]])
    p_value[i] <- stats::pchisq(test_W[i], df = df[i], lower.tail = FALSE)
  }

  out <- data.frame(
    Effects = effect_names,
    Df = df,
    `Hotelling-Lawley` = round(test_W / N, 3),
    `Chi-square` = round(test_W, 3),
    `p-value` = round(p_value, 3),
    check.names = FALSE
  )

  return(out)
}


# mc_manova_disp <- function(object, idx, names, ...) {
#   disp <- coef(object, type = "tau")[,1]
#   n_disp <- length(disp)
#   Parameters <- coef(object, type =  "tau")$Parameters
#   VCOV <- vcov(object)[Parameters, Parameters]
#   ## Conditional variance-covariance model
#   FF <- mc_build_F(vector = idx)
#   G <- Diagonal(length(object$mu_list), 1)
#   CC <- lapply(FF, function(x, G){kronecker(G,x)}, G = G)
#   N <- object$n_obs
#   test_W <- c()
#   df <- c()
#   p_value <- c()
#   for(i in 1:length(CC)) {
#     test_W[i] <- as.numeric(t(CC[[i]]%*%disp)%*%
#                               solve(CC[[i]]%*%VCOV%*%t(CC[[i]]))
#                             %*%(CC[[i]]%*%disp))
#     df[i] <- dim(CC[[i]])[1]
#     p_value[i] <- stats::pchisq(test_W[i], df = df[i], lower.tail = FALSE)
#   }
#   out <- data.frame("Effects" = names, "Df" = df,
#                     "Hotelling-Lawley" = round(test_W/N,3),
#                     "Qui-square" = round(test_W,3),
#                     "p-value" = round(p_value, 3))
#   return(out)
# }
