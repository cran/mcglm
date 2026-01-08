#' @title Quasi-Score Function
#' @author Wagner Hugo Bonat
#'
#' @description
#' Computes the quasi-score function for the regression parameters in
#' multivariate covariance generalized linear models, together with its
#' associated sensitivity and variability matrices. These quantities are
#' key components of the estimating function approach used to fit
#' \code{mcglm} models.
#'
#' @param D A numeric matrix corresponding to the derivative of the mean
#'   vector with respect to the regression parameters. This matrix is
#'   typically obtained from the output of
#'   \code{\link[mcglm]{mc_link_function}}.
#' @param inv_C A numeric matrix giving the inverse of the covariance
#'   matrix of the response vector, usually obtained from
#'   \code{\link[mcglm]{mc_build_C}}.
#' @param y_vec A numeric vector containing the stacked observed responses.
#' @param mu_vec A numeric vector containing the stacked fitted mean values.
#' @param W A numeric matrix of weights, typically diagonal, accounting for
#'   missing observations or differential weighting of the responses.
#'
#' @return
#' A list with the following components:
#' \describe{
#'   \item{Score}{A numeric vector containing the quasi-score values for the
#'   regression parameters.}
#'   \item{Sensitivity}{A numeric matrix giving the sensitivity matrix of
#'   the quasi-score function.}
#'   \item{Variability}{A numeric matrix giving the variability matrix of
#'   the quasi-score function.}
#' }
#'
#' @details
#' Let \eqn{y} denote the response vector, \eqn{\mu} its mean, \eqn{D} the
#' derivative of \eqn{\mu} with respect to the regression parameters, and
#' \eqn{C} the covariance matrix of \eqn{y}. The quasi-score is defined as
#' \deqn{
#'   U_\beta = D^\top C^{-1} W (y - \mu),
#' }
#' where \eqn{W} is a weight matrix. The sensitivity and variability matrices
#' are computed according to standard estimating function theory and are
#' used in the iterative fitting algorithm and for inference.
#'
#' This function is internal and not intended to be called directly by
#' end users.
#'
#' @keywords internal
#' @export

mc_quasi_score <- function(D, inv_C, y_vec, mu_vec, W) {
    res <- y_vec - mu_vec
    t_D <- t(D)
    part1 <- t_D %*% inv_C
    score <- part1 %*% W %*% res
    sensitivity <- -part1 %*% W %*% D
    variability <- part1 %*% W^2 %*% D
    output <- list(Score = score, Sensitivity = sensitivity,
                   Variability = variability)
    return(output)
}
