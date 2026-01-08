#' @title Measures of Goodness-of-Fit
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description Extract the pseudo Gaussian log-likelihood (plogLik),
#' pseudo Akaike Information Criterion (pAIC), pseudo Kullback-Leibler
#' Information Criterion (pKLIC) and pseudo Bayesian Information Criterion (pBIC)
#' for objects of \code{mcglm} class.
#'
#' @param object an object or a list of objects representing a model
#' of \code{mcglm} class.
#' @return A data frame with the following columns:
#' \describe{
#'   \item{plogLik}{Numeric value of the pseudo Gaussian log-likelihood.}
#'   \item{Df}{Integer giving the number of estimated parameters.}
#'   \item{pAIC}{Numeric value of the pseudo Akaike Information Criterion.}
#'   \item{pKLIC}{Numeric value of the pseudo Kullback–Leibler Information Criterion.}
#'   \item{BIC}{Numeric value of the pseudo Bayesian Information Criterion.}
#' }
#'
#' @source Bonat, W. H. (2018). Multiple Response Variables Regression
#' Models in R: The mcglm Package. Journal of Statistical Software, 84(4):1--30.
#'
#' @source Wang, M. (2014). Generalized Estimating Equations in Longitudinal Data
#' Analysis: A Review and Recent Developments. Advances in Statistics, 1(1)1--13.
#'
#' @seealso \code{plogLik}, \code{pAIC}, \code{pKLIC} and \code{pBIC}.
#' @export

gof <- function(object) {
  pl <- plogLik(object, verbose = FALSE)
  AIC <- pAIC(object, verbose = FALSE)
  KLIC <- pKLIC(object, verbose = FALSE)
  BIC <- pBIC(object, verbose = FALSE)
  output <- data.frame("plogLik" = pl$plogLik, "Df" = pl$df,
                       "pAIC" = AIC$pAIC,"pKLIC" = KLIC$pKLIC,
                       "BIC" = BIC)
  return(output)
}
