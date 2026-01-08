#' @title Update Covariance Parameters
#' @author Wagner Hugo Bonat
#'
#' @description
#' Updates the covariance and power parameter vectors stored in a list used
#' during the model fitting algorithm. This function is intended for
#' internal use and distributes covariance parameters across response-specific
#' components and the correlation parameter.
#'
#' @param list_initial
#' A list containing the current values of the model parameters. The elements
#' \code{tau}, \code{power}, and \code{rho} are updated by this function.
#'
#' @param covariance
#' A numeric vector containing the current values of the covariance parameters
#' for all response variables and, if applicable, the correlation parameter.
#'
#' @param list_power_fixed
#' A list of logical values indicating for each response whether the power
#' parameter should remain fixed (\code{TRUE}) or be updated (\code{FALSE}).
#'
#' @param information
#' A list containing model dimension information, including the number of
#' \code{tau} and \code{power} parameters per response variable, and the number
#' of correlation parameters. Typically the output of \code{\link{mc_getInformation}}.
#'
#' @param n_resp
#' An integer specifying the number of response variables in the model.
#'
#' @return
#' A list with the same structure as \code{list_initial}, where the elements
#' \code{tau}, \code{power}, and \code{rho} contain the updated covariance
#' parameter values for each response variable and the correlation parameter.
#' The returned object is intended for internal use only.
#'
#' @keywords internal

mc_updateCov <- function(list_initial, covariance, list_power_fixed,
                         information, n_resp) {
    rho_cod <- rep("rho", information$n_rho)
    tau_cod <- list()
    power_cod <- list()
    for (i in 1:n_resp) {
        power_cod[[i]] <- rep(paste("power", i, sep = ""),
                              information$n_power[[i]])
        tau_cod[[i]] <- rep(paste("tau", i, sep = ""),
                            information$n_tau[[i]])
    }
    temp <- data.frame(values = covariance,
                       cod = c(rho_cod,
                               do.call(c, Map(c, power_cod, tau_cod))))
    cod.tau <- paste("tau", 1:n_resp, sep = "")
    for (i in 1:n_resp) {
        list_initial$tau[[i]] <-
            temp[which(temp$cod == cod.tau[i]), ]$values
    }
    cod.power <- paste("power", 1:n_resp, sep = "")
    for (i in 1:n_resp) {
        if (list_power_fixed[[i]] == FALSE) {
            list_initial$power[[i]] <-
                temp[which(temp$cod == cod.power[i]), ]$values
        }
    }
    if (length(information$n_betas) != 1) {
        list_initial$rho <-
            temp[which(temp$cod == "rho"), ]$values
    }
    return(list_initial)
}
