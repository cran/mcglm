#' @title Variance Functions for Generalized Linear Models
#' @author Wagner Hugo Bonat
#'
#' @description
#' Computes the variance function and its derivatives with respect to
#' regression, dispersion, and power parameters. This function supports
#' standard power variance functions as well as binomial responses.
#' Intended primarily for internal use in model fitting.
#'
#' @param mu Numeric vector of expected values. Typically obtained from
#'     \code{\link{mc_link_function}}.
#'
#' @param power Numeric value (for \code{power} and \code{binomialP})
#'     or numeric vector of length two (for \code{binomialPQ}) representing
#'     the power parameters of the variance function.
#'
#' @param Ntrial Positive integer or numeric. Number of trials for
#'     binomial response variables.
#'
#' @param variance Character string specifying the variance function type:
#'     \code{"power"}, \code{"binomialP"}, or \code{"binomialPQ"}.
#'
#' @param inverse Logical. If \code{TRUE}, computes the inverse square root
#'     of the variance function.
#'
#' @param derivative_power Logical. If \code{TRUE}, computes the derivative
#'     with respect to the power parameter.
#'
#' @param derivative_mu Logical. If \code{TRUE}, computes the derivative
#'     with respect to \code{mu}.
#'
#' @return A named list containing one or more of the following elements,
#'     depending on the combination of logical arguments:
#'     \describe{
#'       \item{V_sqrt}{Square root of the variance function.}
#'       \item{V_inv_sqrt}{Inverse square root of the variance function.}
#'       \item{D_V_sqrt_power}{Derivative of V_sqrt with respect to the power parameter.}
#'       \item{D_V_inv_sqrt_power}{Derivative of V_inv_sqrt with respect to the power parameter.}
#'       \item{D_V_sqrt_mu}{Derivative of V_sqrt with respect to \code{mu}.}
#'       \item{D_V_inv_sqrt_mu}{Derivative of V_inv_sqrt with respect to \code{mu}.}
#'     }
#'
#' @details
#' The function computes the variance function and its derivatives used in
#' the estimation of generalized linear models for multiple response variables.
#' For binomial responses, it accounts for the number of trials and supports
#' both single (\code{binomialP}) and double (\code{binomialPQ}) power specifications.
#'
#' @seealso \code{\link{mc_link_function}}
#'
#' @source Bonat, W. H. and Jorgensen, B. (2016) Multivariate covariance
#' generalized linear models. Journal of the Royal Statistical Society:
#' Series C (Applied Statistics), 65:649--675.
#'
#' @examples
#' x1 <- seq(-1, 1, length.out = 5)
#' X <- model.matrix(~x1)
#' mu <- mc_link_function(beta = c(1, 0.5), X = X, offset = NULL, link = "logit")
#' mc_variance_function(mu = mu$mu, power = c(2, 1), Ntrial = 1,
#'                      variance = "binomialPQ", inverse = FALSE,
#'                      derivative_power = TRUE, derivative_mu = TRUE)
#'
#' @export


## Generic variance function -------------------------------------------
mc_variance_function <- function(mu, power, Ntrial,
                                 variance, inverse,
                                 derivative_power,
                                 derivative_mu) {
    assert_that(is.logical(inverse))
    assert_that(is.logical(derivative_power))
    assert_that(is.logical(derivative_mu))
    switch(variance,
           power = {
               output <- mc_power(mu = mu, power = power,
                                  inverse = inverse,
                                  derivative_power = derivative_power,
                                  derivative_mu = derivative_mu)
           },
           binomialP = {
               output <- mc_binomialP(mu = mu, power = power,
                                      Ntrial = Ntrial,
                                      inverse = inverse,
                                      derivative_power =
                                          derivative_power,
                                      derivative_mu = derivative_mu)
           },
           binomialPQ = {
               output <- mc_binomialPQ(mu = mu, power = power,
                                       Ntrial = Ntrial,
                                       inverse = inverse,
                                       derivative_power =
                                           derivative_power,
                                       derivative_mu = derivative_mu)
           },
           stop(gettextf("%s variance function not recognised",
                         sQuote(variance)), domain = NA))
    return(output)
}

#' @rdname mc_variance_function
## Power variance function ---------------------------------------------
mc_power <- function(mu, power, inverse,
                     derivative_power,
                     derivative_mu) {
    ## The observed value can be zero, but not the expected value.
    assert_that(all(mu > 0))
    assert_that(is.number(power))
    mu.power <- mu^power
    sqrt.mu.power <- sqrt(mu.power)
    n <- length(mu)
    if (inverse == TRUE & derivative_power == TRUE &
            derivative_mu == FALSE) {
        output <- list(
            V_inv_sqrt = Diagonal(n = n, 1/sqrt.mu.power),
            D_V_inv_sqrt_power =
                Diagonal(n = n,
                         -(mu.power * log(mu))/(2 * (mu.power)^(1.5))))
    }
    if (inverse == TRUE & derivative_power == FALSE &
            derivative_mu == FALSE) {
        output <- list(V_inv_sqrt = Diagonal(n = n, 1/sqrt.mu.power))
    }
    if (inverse == FALSE & derivative_power == TRUE &
            derivative_mu == FALSE) {
        output <- list(
            V_sqrt = Diagonal(n = n, sqrt.mu.power),
            D_V_sqrt_power =
                Diagonal(n = n,
                         +(mu.power * log(mu))/(2 * sqrt.mu.power)))
    }
    if (inverse == FALSE & derivative_power == FALSE &
            derivative_mu == FALSE) {
        output <- list(V_sqrt = Diagonal(n = n, sqrt.mu.power))
    }
    if (inverse == TRUE & derivative_power == TRUE &
            derivative_mu == TRUE) {
        output <- list(
            V_inv_sqrt = Diagonal(n = n, 1/sqrt.mu.power),
            D_V_inv_sqrt_power =
                Diagonal(n = n,
                         -(mu.power * log(mu))/(2 * (mu.power)^(1.5))),
            D_V_inv_sqrt_mu = -(mu^(power -  1) * power)/
                                   (2 * (mu.power)^(1.5)))
    }
    if (inverse == TRUE & derivative_power == FALSE &
            derivative_mu == TRUE) {
        output <- list(
            V_inv_sqrt = Diagonal(n = n, 1/sqrt.mu.power),
            D_V_inv_sqrt_mu = -(mu^(power - 1) * power)/
                                   (2 * (mu.power)^(1.5)))
    }
    if (inverse == FALSE & derivative_power == TRUE &
            derivative_mu == TRUE) {
        output <- list(
            V_sqrt = Diagonal(n = n, sqrt.mu.power),
            D_V_sqrt_power =
                Diagonal(n = n, (mu.power * log(mu))/
                                    (2 * sqrt.mu.power)),
            D_V_sqrt_mu = (mu^(power - 1) * power)/(2 * sqrt.mu.power))
    }
    if (inverse == FALSE & derivative_power == FALSE &
            derivative_mu == TRUE) {
        output <- list(V_sqrt = Diagonal(n = n, sqrt.mu.power),
                       D_V_sqrt_mu = (mu^(power - 1) * power)/
                                         (2 * sqrt.mu.power))
    }
    return(output)
}

#' @rdname mc_variance_function
#' @usage mc_binomialP(mu, power, inverse, Ntrial,
#'                     derivative_power, derivative_mu)
## BinomialP variance function
## -----------------------------------------
mc_binomialP <- function(mu, power, inverse, Ntrial,
                         derivative_power,
                         derivative_mu) {
    ## The observed value can be 0 and 1, but not the expected value
    assert_that(all(mu > 0))
    assert_that(all(mu < 1))
    assert_that(is.number(power))
    assert_that(all(Ntrial > 0))
    constant <- (1/Ntrial)
    mu.power <- mu^power
    mu.power1 <- (1 - mu)^power
    mu1mu <- constant * (mu.power * mu.power1)
    sqrt.mu1mu <- sqrt(mu1mu)
    n <- length(mu)
    if (inverse == TRUE & derivative_power == TRUE &
            derivative_mu == FALSE) {
        output <- list(
            V_inv_sqrt = Diagonal(n = n, 1/sqrt.mu1mu),
            D_V_inv_sqrt_power =
                Diagonal(n = n, -(log(1 - mu) * mu1mu +
                                  log(mu) * mu1mu)/(2 * (mu1mu^(1.5)))))
    }
    if (inverse == TRUE & derivative_power == FALSE &
            derivative_mu == FALSE) {
        output <- list(V_inv_sqrt = Diagonal(n = n, 1/sqrt.mu1mu))
    }
    if (inverse == FALSE & derivative_power == TRUE &
            derivative_mu == FALSE) {
        output <- list(
            V_sqrt = Diagonal(n = n, sqrt.mu1mu),
            D_V_sqrt_power = Diagonal(n = n, (log(1 - mu) * mu1mu +
                                              log(mu) * mu1mu)/
                                                 (2 * sqrt.mu1mu)))
    }
    if (inverse == FALSE & derivative_power == FALSE &
            derivative_mu == FALSE) {
        output <- list(V_sqrt = Diagonal(n = n, sqrt.mu1mu))
    }
    if (inverse == TRUE & derivative_power == TRUE &
            derivative_mu == TRUE) {
        output <- list(
            V_inv_sqrt = Diagonal(n = n, 1/sqrt.mu1mu),
            D_V_inv_sqrt_power =
                Diagonal(n = n, -(log(1 - mu) * mu1mu + log(mu) *
                                  mu1mu)/(2 * (mu1mu^(1.5)))),
            D_V_inv_sqrt_mu = -(constant * (mu.power1 *
                                            (mu^(power - 1)) * power) -
                                constant * (((1 - mu)^(power - 1)) *
                                            mu.power * power))/
                                   (2 * (mu1mu^(1.5))))
    }
    if (inverse == TRUE & derivative_power == FALSE &
            derivative_mu == TRUE) {
        output <- list(
            V_inv_sqrt = Diagonal(n = n, 1/sqrt.mu1mu),
            D_V_inv_sqrt_mu = -(constant *
                                (mu.power1 * (mu^(power - 1)) * power) -
                                constant * (((1 - mu)^(power - 1)) *
                                            mu.power * power))/
                                   (2 * (mu1mu^(1.5))))
    }
    if (inverse == FALSE & derivative_power == TRUE &
            derivative_mu == TRUE) {
        output <- list(
            V_sqrt = Diagonal(n = n, sqrt.mu1mu),
            D_V_sqrt_power = Diagonal(n = n, (log(1 - mu) * mu1mu +
                                              log(mu) * mu1mu)/
                                                 (2 * sqrt.mu1mu)),
            D_V_sqrt_mu = (constant *
                           (mu.power1 * (mu^(power - 1)) * power) -
                           constant * (((1 - mu)^(power - 1)) *
                                       mu.power * power))/
                              (2 * sqrt.mu1mu))
    }
    if (inverse == FALSE & derivative_power == FALSE &
            derivative_mu == TRUE) {
        output <- list(
            V_sqrt = Diagonal(n = n, sqrt.mu1mu),
            D_V_sqrt_mu = (constant *
                           (mu.power1 * (mu^(power - 1)) * power) -
                           constant * (((1 - mu)^(power - 1)) *
                                       mu.power * power))/
                              (2 * sqrt.mu1mu))
    }
    return(output)
}

#' @rdname mc_variance_function
#' @usage mc_binomialPQ(mu, power, inverse, Ntrial,
#'                      derivative_power, derivative_mu)
## BinomialPQ variance function ----------------------------------------
mc_binomialPQ <- function(mu, power, inverse,
                          Ntrial, derivative_power,
                          derivative_mu) {
    ## The observed value can be 0 and 1, but not the expected value
    assert_that(all(mu > 0))
    assert_that(all(mu < 1))
    assert_that(length(power) == 2)
    assert_that(all(Ntrial > 0))
    constant <- (1/Ntrial)
    p <- power[1]
    q <- power[2]
    mu.p <- mu^p
    mu1.q <- (1 - mu)^q
    mu.p.mu.q <- mu.p * mu1.q
    mu1mu <- mu.p.mu.q * constant
    sqrt.mu1mu <- sqrt(mu1mu)
    n <- length(mu)
    if (inverse == TRUE & derivative_power == TRUE &
            derivative_mu == FALSE) {
        denominator <- (2 * (mu1mu^1.5) * Ntrial)
        output <- list(
            V_inv_sqrt = Diagonal(n = n, 1/sqrt.mu1mu),
            D_V_inv_sqrt_p = Diagonal(n = n,
                                      -(mu.p.mu.q * log(mu))/
                                           denominator),
            D_V_inv_sqrt_q = Diagonal(n = n,
                                      -mu.p.mu.q * log(1 - mu)/
                                           denominator))
    }
    if (inverse == TRUE & derivative_power == FALSE &
            derivative_mu == FALSE) {
        output <- list(V_inv_sqrt = Diagonal(n = n, 1/sqrt.mu1mu))
    }
    if (inverse == FALSE & derivative_power == TRUE &
            derivative_mu == FALSE) {
        denominator <- 2 * sqrt.mu1mu * Ntrial
        output <- list(
            V_sqrt = Diagonal(n = n, sqrt.mu1mu),
            D_V_sqrt_p = Diagonal(n = n,
                                  +(mu.p.mu.q * log(mu))/denominator),
            D_V_sqrt_q = Diagonal(n = n,
                                  +(mu.p.mu.q * log(1 - mu))/
                                       denominator))
    }
    if (inverse == FALSE & derivative_power == FALSE &
            derivative_mu == FALSE) {
        output <- list(V_sqrt = Diagonal(n = n, sqrt.mu1mu))
    }
    if (inverse == TRUE & derivative_power == TRUE &
            derivative_mu == TRUE) {
        denominator <- (2 * (mu1mu^1.5) * Ntrial)
        output <- list(
            V_inv_sqrt = Diagonal(n = n, 1/sqrt.mu1mu),
            D_V_inv_sqrt_p = Diagonal(n = n,
                                      -(mu.p.mu.q * log(mu))/
                                           denominator),
            D_V_inv_sqrt_q = Diagonal(n = n,
                                      -mu.p.mu.q *
                                           log(1 - mu)/denominator),
            D_V_inv_sqrt_mu = -(constant *
                                (mu1.q * (mu^(p - 1)) * p) -
                                constant * (((1 - mu)^(q - 1)) *
                                            mu.p * q))/
                                   (2 * (mu1mu^1.5)))
    }
    if (inverse == TRUE & derivative_power == FALSE &
            derivative_mu == TRUE) {
        output <- list(
            V_inv_sqrt = Diagonal(n = n, 1/sqrt.mu1mu),
            D_V_inv_sqrt_mu = -(constant * (mu1.q * (mu^(p - 1)) * p) -
                                constant * (((1 - mu)^(q - 1)) *
                                            mu.p * q))/
                                   (2 * (mu1mu^1.5)))
    }
    if (inverse == FALSE & derivative_power == TRUE &
            derivative_mu == TRUE) {
        denominator1 <- 2 * sqrt.mu1mu
        denominator2 <- denominator1 * Ntrial
        output <- list(
            V_sqrt = Diagonal(n = n, sqrt.mu1mu),
            D_V_sqrt_p = Diagonal(n = n, (mu.p.mu.q * log(mu))/
                                             denominator2),
            D_V_sqrt_q = Diagonal(n = n, (mu.p.mu.q * log(1 - mu))/
                                             denominator2),
            D_V_sqrt_mu = (constant * (mu1.q * (mu^(p - 1)) * p) -
                           constant * (((1 - mu)^(q - 1)) * mu.p * q))/
                              denominator1)
    }
    if (inverse == FALSE & derivative_power == FALSE &
            derivative_mu == TRUE) {
        output <- list(
            V_sqrt = Diagonal(n = n, sqrt.mu1mu),
            D_V_sqrt_mu = (constant * (mu1.q * (mu^(p - 1)) * p) -
                           constant * (((1 - mu)^(q - 1)) * mu.p * q))/
                              (2 * sqrt.mu1mu))
    }
    return(output)
}
