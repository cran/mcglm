#' @title Pearson Estimating Function
#' @author Wagner Hugo Bonat
#'
#' @description
#' Computes the Pearson estimating function for the dispersion parameters
#' in multivariate covariance generalized linear models, together with its
#' associated sensitivity and variability matrices. This function is used
#' internally in the estimating function framework adopted by
#' \code{mcglm}.
#'
#' @param y_vec Numeric vector of observed responses stacked across
#'   response variables.
#' @param mu_vec Numeric vector of fitted means corresponding to
#'   \code{y_vec}.
#' @param Cfeatures A list containing covariance-related components,
#'   typically including the covariance matrix \code{C}, its inverse
#'   \code{inv_C}, and the derivatives of \code{C} with respect to the
#'   dispersion parameters (\code{D_C}).
#' @param inv_J_beta Optional matrix giving the inverse of the sensitivity
#'   matrix associated with the regression parameters. Required only when
#'   bias correction is requested.
#' @param D Optional matrix of derivatives of the mean vector with respect
#'   to the regression parameters. Required only when bias correction is
#'   requested.
#' @param correct Logical indicating whether the bias-corrected Pearson
#'   estimating function should be computed. Defaults to \code{FALSE}.
#' @param compute_sensitivity Logical indicating whether the sensitivity
#'   matrix of the Pearson estimating function should be computed.
#'   Defaults to \code{TRUE}.
#' @param compute_variability Logical indicating whether the variability
#'   matrix of the Pearson estimating function should be computed.
#'   Defaults to \code{FALSE}.
#' @param W Numeric vector or diagonal matrix of weights associated with
#'   the observations.
#'
#' @return
#' A list with the following components:
#' \describe{
#'   \item{Score}{A numeric vector containing the values of the Pearson
#'   estimating function for the dispersion parameters.}
#'   \item{Sensitivity}{A matrix giving the sensitivity (expected Jacobian)
#'   of the Pearson estimating function. Returned only if
#'   \code{compute_sensitivity = TRUE}.}
#'   \item{Variability}{A matrix giving the variability of the Pearson
#'   estimating function. Returned only if
#'   \code{compute_variability = TRUE}.}
#'   \item{Extra}{A list of intermediate quantities used in the
#'   computation, mainly products involving derivatives of the covariance
#'   matrix.}
#' }
#'
#' @details
#' The Pearson estimating function is based on quadratic forms of the
#' residuals and the inverse covariance matrix. When
#' \code{correct = TRUE}, a bias-corrected version is computed using the
#' correction term described in Bonat and Jørgensen (2016). The sensitivity
#' and variability matrices correspond to Equations (6), (7) and (8) of
#' that reference.
#'
#' This function is intended for internal use and is not designed to be
#' called directly by end users.
#'
#' @keywords internal
#'
#' @source
#' Bonat, W. H. and Jørgensen, B. (2016). Multivariate covariance generalized
#' linear models. Journal of the Royal Statistical Society: Series C,
#' 65, 649--675.


mc_pearson <- function(y_vec, mu_vec, Cfeatures, inv_J_beta = NULL,
                       D = NULL, correct = FALSE,
                       compute_sensitivity = TRUE,
                       compute_variability = FALSE,
                       W) {

    product <- lapply(Cfeatures$D_C, mc_multiply,
                      bord2 = Cfeatures$inv_C)
    res <- y_vec - mu_vec
    pearson_score <- unlist(lapply(product, mc_core_pearson,
                                   inv_C = Cfeatures$inv_C, res = res, W = W))

    sensitivity <- matrix(NA, length(product), length(product))
    if(compute_sensitivity == TRUE) {
      sensitivity <- mc_sensitivity(product, W = W)
    }

    output <- list(Score = pearson_score, Sensitivity = sensitivity,
                   Extra = product)
    if (correct == TRUE) {
        correction <- mc_correction(D_C = Cfeatures$D_C,
                                    inv_J_beta = inv_J_beta, D = D,
                                    inv_C = Cfeatures$inv_C)
        output <- list(Score = pearson_score + correction,
                       Sensitivity = sensitivity, Extra = product)
    }
    if (compute_variability == TRUE) {
        variability <- mc_variability(sensitivity = sensitivity,
                                      product = product,
                                      inv_C = Cfeatures$inv_C,
                                      C = Cfeatures$C, res = res, W = W)
        output$Variability <- variability
    }
    return(output)
}
