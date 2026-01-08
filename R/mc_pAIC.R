#' @title Pseudo Akaike Information Criterion
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description
#' Computes the pseudo Akaike information criterion (pAIC) for fitted
#' multivariate covariance generalized linear models. The pAIC is defined
#' as
#' \deqn{pAIC = -2 \, \ell_p + 2 \, \text{df},}
#' where \eqn{\ell_p} is the pseudo log-likelihood and \eqn{\text{df}}
#' denotes the effective number of parameters in the model.
#'
#' This criterion is intended for model comparison within the class of
#' \code{mcglm} models fitted to the same data.
#'
#' @param object An object of class \code{mcglm} or a list of such objects.
#'   When a list is provided, the pseudo log-likelihood is computed for
#'   each model.
#' @param verbose Logical indicating whether the pAIC value should be
#'   printed to the console. Defaults to \code{TRUE}.
#'
#' @return
#' An (invisible) named list with a single element:
#' \describe{
#'   \item{pAIC}{A numeric value giving the pseudo Akaike information
#'   criterion associated with the fitted model(s).}
#' }
#'
#' @details
#' The pAIC is based on the pseudo log-likelihood returned by
#' \code{\link{plogLik}} and should be used with caution, as it does not
#' correspond to a true likelihood-based information criterion.
#' Comparisons are meaningful only for models fitted to the same response
#' data.
#'
#' @seealso
#' \code{gof}, \code{plogLik}, \code{ESS}, \code{pKLIC},
#' \code{GOSHO}, \code{RJC}
#'
#' @source
#' Bonat, W. H. (2018). Multiple Response Variables Regression Models in R:
#' The mcglm Package. Journal of Statistical Software, 84(4), 1--30.
#'
#' @export

pAIC <- function(object, verbose = TRUE) {
  Pseudo <- plogLik(object = object, verbose = FALSE)
  pAIC <- 2*Pseudo$df - 2*Pseudo$plogLik
  if (verbose) cat("pAIC", pAIC)
  return(invisible(list("pAIC" = pAIC)))
}
