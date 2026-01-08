#' @title Robust Standard Errors for Regression Parameters
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description
#' Computes cluster-robust (sandwich-type) standard errors for the
#' regression parameters of an object of class \code{mcglm}, accounting
#' for within-cluster correlation.
#'
#' @details
#' The robust variance--covariance matrix is obtained using an empirical
#' estimator based on clustered residuals and the sensitivity matrix of
#' the estimating equations. The implementation assumes that the data are
#' correctly ordered such that observations belonging to the same cluster
#' are stored in contiguous rows.
#'
#' @param object An object of class \code{mcglm} representing a fitted
#'   marginal model.
#' @param id An integer or factor vector identifying clusters or subjects.
#'   Its length and ordering must match the number and ordering of the
#'   observations used to fit the model.
#'
#' @return
#' A list with two components:
#' \describe{
#'   \item{Std.Error}{A numeric vector containing the robust standard
#'   errors of the regression parameter estimates.}
#'   \item{vcov}{A numeric matrix giving the robust variance--covariance
#'   matrix of the regression parameter estimates.}
#' }
#' The returned objects are computed under the assumption that the data
#' are in the correct cluster order.
#'
#' @source
#' Nuamah, I. F., Qu, Y., and Aminu, S. B. (1996). A SAS macro for stepwise
#' correlated binary regression. \emph{Computer Methods and Programs in
#' Biomedicine}, 49, 199--210.
#'
#' @seealso \code{mc_bias_correct_std}
#'
#' @export

mc_robust_std <- function(object, id) {
    inv_M <- object$inv_S_beta
    temp_data <- data.frame(res = object$residuals, id)
    temp_data_group <- split(temp_data, temp_data$id)
    r_rT <- bdiag(lapply(temp_data_group,
                         function(x) {
                             tcrossprod(x[, 1])
                         }))
    D <- bdiag(lapply(object$mu_list, function(x) x$D))
    p1 <- object$inv_C %*% D
    V_robust <- inv_M %*% (t(p1) %*% r_rT %*% p1) %*% inv_M
    output <- list("Std.Error" = sqrt(diag(V_robust)),
                                      "vcov" = V_robust)
    return(output)
}
