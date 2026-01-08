#' @title MANOVA-Type Test for Multivariate Covariance GLMs
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description
#' Performs a MANOVA-type Wald test for multivariate covariance
#' generalized linear models fitted using \code{\link{mcglm}}.
#' The test is based on quadratic forms of the estimated regression
#' parameters and their covariance matrix, yielding statistics
#' analogous to the Hotelling–Lawley trace.
#'
#' @param object An object of class \code{"mcglm"}.
#' @param ... Further arguments (currently not used).
#'
#' @return
#' A data frame containing the MANOVA-type test results with the following
#' columns:
#' \describe{
#'   \item{Effects}{Names of the tested model effects.}
#'   \item{Df}{Degrees of freedom associated with each effect.}
#'   \item{Hotelling-Lawley}{Hotelling–Lawley trace statistic.}
#'   \item{Chi-square}{Chi-square test statistic.}
#'   \item{p-value}{P-values from the chi-square approximation.}
#' }
#'
#' @seealso \code{\link{mcglm}}, \code{\link{coef.mcglm}},
#'   \code{\link{vcov.mcglm}}
#'
#' @export

# manova for objects of mcglm class ----------------------------------
mc_manova <- function(object, ...) {

  if (!inherits(object, "mcglm")) {
    stop("object must be of class 'mcglm'.", call. = FALSE)
  }

  beta <- coef(object, type = "beta")[, 1]
  n_beta <- length(beta)

  VCOV <- vcov(object)[1:n_beta, 1:n_beta, drop = FALSE]

  FF <- mc_build_F(vector = attr(object$list_X[[1]], "assign"))

  G <- Diagonal(length(object$mu_list), 1)

  CC <- lapply(FF, function(x) kronecker(G, x))

  N <- object$n_obs - ncol(object$list_X[[1]])
  if (N <= 0) {
    stop("Residual degrees of freedom must be positive.", call. = FALSE)
  }

  test_W <- numeric(length(CC))
  df <- integer(length(CC))
  p_value <- numeric(length(CC))

  for (i in seq_along(CC)) {
    CVC <- CC[[i]] %*% VCOV %*% t(CC[[i]])

    CVC_inv <- tryCatch(
      solve(CVC),
      error = function(e) {
        stop("Singular covariance matrix in MANOVA test.", call. = FALSE)
      }
    )

    test_W[i] <- as.numeric(
      t(CC[[i]] %*% beta) %*% CVC_inv %*% (CC[[i]] %*% beta)
    )

    df[i] <- nrow(CC[[i]])
    p_value[i] <- stats::pchisq(test_W[i], df = df[i], lower.tail = FALSE)
  }

  effect_names <- c(
    "Intercept",
    attr(stats::terms(object$linear_pred[[1]]), "term.labels")
  )

  out <- data.frame(
    Effects = effect_names,
    Df = df,
    `Hotelling-Lawley` = round(test_W / N, 4),
    `Chi-square` = round(test_W, 4),
    `p-value` = round(p_value, 4),
    check.names = FALSE
  )

  return(out)
}

#
#
#
#
#
#
#
# mc_manova <- function (object, ...)
# {
#   beta <- coef(object, type = "beta")[, 1]
#   n_beta <- length(beta)
#   VCOV <- vcov(object)[1:n_beta, 1:n_beta]
#   FF <- mc_build_F(vector = attr(object$list_X[[1]], "assign"))
#   G <- Diagonal(length(object$mu_list), 1)
#   CC <- lapply(FF, function(x, G) {
#     kronecker(G, x)
#   }, G = G)
#   N <- object$n_obs - dim(object$list_X[[1]])[2]
#   test_W <- c()
#   df <- c()
#   p_value <- c()
#   for (i in 1:length(CC)) {
#     test_W[i] <- as.numeric(t(CC[[i]] %*% beta) %*% solve(CC[[i]] %*%
#                                                             VCOV %*% t(CC[[i]])) %*% (CC[[i]] %*% beta))
#     df[i] <- dim(CC[[i]])[1]
#     p_value[i] <- stats::pchisq(test_W[i], df = df[i], lower.tail = FALSE)
#   }
#   names <- c("Intercept", attr(stats::terms(object$linear_pred[[1]]),
#                                "term.labels"))
#   out <- data.frame(Effects = names, Df = df,
#                     `Hotelling-Lawley` = round(test_W/N,4),
#                     `Qui-square` = round(test_W, 4),
#                     `p-value` = round(p_value, 4))
#   return(out)
# }
