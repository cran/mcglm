#' @title Getting information about model parameters
#' @author Wagner Hugo Bonat
#'
#' @description This computes all information required about the number
#'     of model parameters.
#'
#' @param list_initial A list of initial values.
#' @param list_power_fixed A list of logical specyfing if the power
#'     parameters should be estimated or not.
#' @param n_resp A number specyfing the nmber of response variables.
#' @keywords internal
#' @return A list containing information about the number of model parameters:
#' \describe{
#'   \item{n_betas}{A list with the number of regression parameters
#'   (\eqn{\beta}) for each response variable.}
#'   \item{n_taus}{A list with the number of dispersion parameters
#'   (\eqn{\tau}) for each response variable.}
#'   \item{n_power}{A list with the number of power parameters for each
#'   response variable, accounting for fixed parameters.}
#'   \item{n_rho}{An integer giving the number of correlation parameters.}
#'   \item{n_cov}{An integer giving the total number of covariance-related
#'   parameters, including power, dispersion, and correlation parameters.}
#' }

mc_getInformation <- function(list_initial, list_power_fixed,
                              n_resp) {
    n_betas <- lapply(list_initial$regression, length)
    n_taus <- lapply(list_initial$tau, length)
    n_power <- lapply(list_initial$power, length)
    for (i in 1:n_resp) {
        if (list_power_fixed[[i]] == TRUE) {
            n_power[i] <- 0
        }
    }
    if (n_resp == 1) {
        n_rho <- 0
    }
    if (n_resp != 1) {
        n_rho <- length(list_initial$rho)
    }
    n_cov <- sum(do.call(c, n_power)) + n_rho +
        sum(do.call(c, n_taus))
    saida <- list(n_betas = n_betas, n_taus = n_taus, n_power = n_power,
                  n_rho = n_rho, n_cov = n_cov)
    return(saida)
}
