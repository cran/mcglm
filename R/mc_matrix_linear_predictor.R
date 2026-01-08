#' @title Matrix Linear Predictor
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description
#' Computes the matrix linear predictor used in multivariate covariance
#' generalized linear models. The matrix linear predictor is defined
#' as a linear combination of known matrices weighted by dispersion
#' parameters.
#'
#' @param tau A numeric vector of dispersion parameters.
#' @param Z A list of known matrices with compatible dimensions.
#'
#' @return
#' A matrix of class \code{\link[Matrix]{Matrix}} representing the
#' matrix linear predictor
#' \deqn{U = \tau_1 Z_1 + \cdots + \tau_D Z_D.}
#' The returned matrix has the same dimensions as the elements of \code{Z}.
#' The returned object is intended for internal use only.
#'
#' @details
#' Given a list of known matrices \eqn{(Z_1, \ldots, Z_D)} and a vector
#' of dispersion parameters \eqn{(\tau_1, \ldots, \tau_D)}, this function
#' computes their weighted sum. This object is typically used as a
#' component of the matrix linear predictor in covariance modeling.
#'
#' @seealso \code{mc_id}, \code{mc_dist}, \code{mc_ma}, \code{mc_rw},
#'   \code{mc_mixed}, \code{mc_car}
#'
#' @source Bonat, W. H. (2018). Multiple Response Variables Regression
#' Models in R: The mcglm Package. Journal of Statistical Software, 84(4):1--30.
#'
#' @source Bonat, W. H. and Jorgensen, B. (2016). Multivariate covariance
#' generalized linear models. Journal of the Royal Statistical Society:
#' Series C, 65:649--675.
#'
#' @examples
#' Z0 <- Matrix::Diagonal(5, 1)
#' Z1 <- Matrix::Matrix(rep(1, 5) %*% t(rep(1, 5)))
#' Z <- list(Z0, Z1)
#' mc_matrix_linear_predictor(tau = c(1, 0.8), Z = Z)
#'
#' @export
#' @import Matrix

# Matrix linear predictor ----------------------------------------------
mc_matrix_linear_predictor <- function(tau, Z) {
    if (length(Z) != length(tau)) {
        stop("Incorrect number of parameters")
    }
    output <- mapply("*", Z, tau, SIMPLIFY = FALSE)
    output <- Reduce("+", output)
    return(output)
}
