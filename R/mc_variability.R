#' @title Variability Matrix
#' @author Wagner Hugo Bonat
#'
#' @description
#' Computes the variability matrix associated with the Pearson estimating
#' function. This function is intended for internal use in the fitting
#' algorithm and implements Equation 8 from Bonat and Jorgensen (2016).
#'
#' @param sensitivity A numeric matrix representing the sensitivity matrix.
#'     Typically obtained from \code{\link{mc_sensitivity}}.
#'
#' @param product A list of numeric matrices, usually used to compute
#'     contributions to the variability matrix.
#'
#' @param inv_C A numeric matrix representing the inverse of the covariance
#'     matrix, usually obtained from \code{\link{mc_build_C}}.
#'
#' @param C A numeric matrix representing the covariance matrix, usually
#'     obtained from \code{\link{mc_build_C}}.
#'
#' @param res A numeric vector of residuals, defined as \eqn{y - \mu}.
#'
#' @param W A numeric matrix of weights.
#'
#' @return
#' A symmetric numeric matrix representing the variability matrix associated
#' with the Pearson estimating function. The returned object is intended for internal use only.
#'
#' @keywords internal
#' @details
#' This function implements Equation 8 from Bonat and Jorgensen (2016), which
#' defines the variability matrix used in generalized estimating equations
#' for multiple response variables.

mc_variability <- function(sensitivity, product, inv_C, C, res, W) {
    WE <- lapply(product, mc_multiply2, bord2 = inv_C)
    n_par <- length(product)
    k4 <- res^4 - 3 * diag(C)^2
    #Variability <- matrix(NA, nrow = n_par, ncol = n_par)
    #for (i in 1:n_par) {
    #    for (j in 1:n_par) {
    #        Variability[i, j] <-
    #            as.numeric(-2 * sensitivity[i, j] +
    #                            sum(k4 * diag(W[[i]]) * diag(W[[j]])))
    #    }
    #}
    Sensitivity2 <- mc_sensitivity_op(products = product, W = W^2)
    Sensitivity2 <- forceSymmetric(Sensitivity2)
    #sourceCpp("src/mc_variability_op.cpp")
    W <- as.vector(diag(W))
    Variability = mc_variability_op(sensitivity = Sensitivity2, WE = WE, k4 = k4, W = W)
    #for (i in 1:n_par) {
    #    for (j in 1:i) {
    #        Variability[i, j] <-
    #            as.numeric(-2 * sensitivity[i, j] +
    #                            sum(k4 * diag(W[[i]]) * diag(W[[j]])))
    #    }
    #}
    Variability <- forceSymmetric(Variability)
    return(Variability)
}
