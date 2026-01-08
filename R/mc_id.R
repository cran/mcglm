#' @title Independent Model Structure
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description Builds an identity matrix to be used as a component of
#' the matrix linear predictor. It is in general the first component of
#' the matrix linear predictor, a kind of intercept matrix.
#'
#' @param data the data set to be used.
#' @return A list with a single component:
#' \describe{
#'   \item{Z0}{A sparse identity matrix of class \code{ddiMatrix} with dimension
#'   equal to the number of observations, representing the intercept component
#'   of the matrix linear predictor.}
#' }
#' It is intended to be supplied to the \code{matrix_pred} argument of
#' \code{\link{mcglm}}.
#' @seealso \code{mc_dist}, \code{mc_ma}, \code{mc_rw} and \code{mc_mixed}.
#'
#' @source Bonat, W. H. (2018). Multiple Response Variables Regression
#' Models in R: The mcglm Package. Journal of Statistical Software, 84(4):1--30.
#'
#' @export

mc_id <- function(data) {
  output <- list("Z0" = Diagonal(dim(data)[1],1))
  return(output)
}
