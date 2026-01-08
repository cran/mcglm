#' @title Update Regression Parameters
#' @author Wagner Hugo Bonat
#'
#' @description
#' Updates the regression parameter vectors stored in a list used during
#' the model fitting algorithm. This function is intended for internal
#' use and is called repeatedly to redistribute the regression parameters
#' across response-specific components.
#'
#' @param list_initial
#' A list containing the current values of the model parameters. The
#' element \code{regression} is updated by this function.
#'
#' @param betas
#' A numeric vector containing the current values of the regression
#' parameters for all response variables.
#'
#' @param information
#' A list containing model dimension information, including the number
#' of regression parameters per response variable. Typically the output
#' of \code{\link{mc_getInformation}}.
#'
#' @param n_resp
#' An integer specifying the number of response variables in the model.
#'
#' @return
#' A list with the same structure as \code{list_initial}, where the
#' element \code{regression} contains the updated regression parameter
#' vectors for each response variable. The returned object is intended for internal use only.
#'
#' @keywords internal

mc_updateBeta <- function(list_initial, betas, information, n_resp) {
    cod <- rep(1:n_resp, information$n_betas)
    temp <- data.frame(beta = betas, cod)
    for (k in 1:n_resp) {
        list_initial$regression[[k]] <-
            temp[which(temp$cod == k), ]$beta
    }
    return(list_initial)
}
