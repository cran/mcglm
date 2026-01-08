#' @title Pseudo Bayesian Information Criterion
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description
#' Computes the pseudo Bayesian information criterion (pBIC) for fitted
#' multivariate covariance generalized linear models. The pBIC is defined
#' as
#' \deqn{pBIC = -2 \, \ell_p + \text{df} \log(n),}
#' where \eqn{\ell_p} is the pseudo log-likelihood, \eqn{\text{df}} is the
#' effective number of parameters, and \eqn{n} is the total number of
#' observed responses.
#'
#' This criterion provides a more strongly penalized alternative to
#' \code{pAIC}, favoring more parsimonious models when comparing
#' \code{mcglm} fits to the same data.
#'
#' @param object An object of class \code{mcglm} or a list of such objects.
#'   When a list is supplied, the pseudo log-likelihood and the number of
#'   observations are computed by aggregating information across all
#'   models.
#' @param verbose Logical indicating whether the pBIC value should be
#'   printed to the console. Defaults to \code{TRUE}.
#'
#' @return
#' An (invisible) named list with a single element:
#' \describe{
#'   \item{pBIC}{A numeric value giving the pseudo Bayesian information
#'   criterion associated with the fitted model(s).}
#' }
#'
#' @details
#' The sample size \eqn{n} used in the penalty term corresponds to the
#' total number of observed responses, obtained from the
#' \code{observed} component of the fitted \code{mcglm} object(s).
#' As the pBIC is based on a pseudo log-likelihood, it should be used
#' cautiously and only for relative comparisons among models fitted to
#' the same data set.
#'
#' @seealso
#' \code{gof}, \code{plogLik}, \code{ESS}, \code{pAIC},
#' \code{pKLIC}, \code{GOSHO}, \code{RJC}
#'
#' @source
#' Bonat, W. H. (2018). Multiple Response Variables Regression Models in R:
#' The mcglm Package. Journal of Statistical Software, 84(4), 1--30.
#'
#' @export

pBIC <- function(object, verbose = TRUE) {
  Pseudo <- plogLik(object = object, verbose = FALSE)
  if(inherits(object, "mcglm")) {
    Y = object$observed
  }
  if(inherits(object, "list")) {
    Y <- do.call(c,lapply(object, function(x)as.numeric(x$observed)))
  }
  NS <- length(Y)
  pBIC <- Pseudo$df*log(NS) - 2*Pseudo$plogLik
  if (verbose) cat("pBIC", pBIC)
  return(invisible(list("pBIC" = pBIC)))
}
