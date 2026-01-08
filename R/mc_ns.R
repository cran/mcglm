#' @title Non-structured Covariance Model
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description
#' Constructs the components of the matrix linear predictor associated
#' with a fully non-structured covariance model in multivariate covariance
#' generalized linear models. This specification allows each pair of
#' observations within a unit to have its own covariance parameter,
#' resulting in a highly flexible but parameter-intensive model.
#'
#' Due to the quadratic growth in the number of parameters, this structure
#' is typically suitable only for datasets with a small number of repeated
#' measurements per unit.
#'
#' @param id A character string giving the name of the column in
#'   \code{data} that identifies the observational units (e.g., subjects).
#'   Each unit must have the same number of observations. For time series
#'   or spatial data without replication, the same identifier should be
#'   used for all observations.
#' @param data A \code{data.frame} containing the variables referenced by
#'   \code{id} and, optionally, \code{group}.
#' @param group An optional character string giving the name of a column in
#'   \code{data} that defines groups for which different covariance
#'   structures may be specified. If \code{NULL}, a single non-structured
#'   covariance model is used for all units.
#' @param marca An optional character string specifying the level of
#'   \code{group} for which the non-structured covariance components are
#'   excluded (i.e., set to zero). This allows selective activation of the
#'   non-structured covariance according to group membership.
#'
#' @return
#' A list of symmetric block-diagonal matrices, each representing one
#' covariance component of the non-structured matrix linear predictor.
#' The length of the list is equal to \eqn{n(n - 1) / 2}, where \eqn{n} is
#' the number of observations per unit. Each element of the list is a
#' sparse matrix of class \code{"dgCMatrix"} obtained by stacking unit-
#' specific covariance blocks along the diagonal. These matrices are used
#' internally to construct the dispersion linear predictor in
#' \code{mcglm}.
#'
#' @details
#' The function requires a balanced design, meaning that all units
#' identified by \code{id} must have the same number of observations.
#' An error is raised otherwise. When \code{group} and \code{marca} are
#' provided, covariance components are generated only for units not
#' belonging to the specified level \code{marca}; for those units, the
#' corresponding blocks are set to zero.
#'
#' @seealso
#' \code{mc_id}, \code{mc_dglm}, \code{mc_dist}, \code{mc_ma},
#' \code{mc_rw}, \code{mc_mixed}
#'
#' @source
#' Bonat, W. H. (2018). Multiple Response Variables Regression Models in R:
#' The mcglm Package. Journal of Statistical Software, 84(4), 1--30.
#'
#' @export

mc_ns <- function(id, data, group = NULL, marca = NULL) {
  mc_non_aux <- function(n.resp){
    position <- combn(n.resp,2)
    list.Derivative <- list()
    n.par <- n.resp*(n.resp-1)/2
    for(i in 1:n.par){
      Derivative <- matrix(0, ncol = n.resp, nrow = n.resp)
      Derivative[position[1,i],position[2,i]] <- Derivative[position[2,i],position[1,i]] <- 1
      list.Derivative[i][[1]] <- Derivative}
    return(list.Derivative)
  }
  data[id] <- factor(data[[id]], levels=unique(data[[id]]))
  data.id <- split(data, data[id], drop = TRUE)
  DD <- sum(abs(diff(do.call(c,lapply(data.id, function(x)dim(x)[1])))))
  if( DD != 0) {
    stop("Model requires equal number of observations by id. \n")
  }
  mat.list <- list()
  for(i in 1:length(data.id)) {
    if (!is.null(group)) {
      if (unique(data.id[[i]][group]) == marca) {
      mat.list[[i]] <- lapply(mc_non_aux(dim(data.id[[i]])[1]), function(x){0*x} )
      }
    }
    if(!is.null(group)) {
      if(unique(data.id[[i]][group]) != marca) {
        mat.list[[i]] <- mc_non_aux(dim(data.id[[i]])[1])
      }
    }
    if(is.null(group)) {
      mat.list[[i]] <- mc_non_aux(dim(data.id[[i]])[1])
    }
  }
  non_list <- list()
  for(i in 1:length(mat.list[[1]])) {
    non_list[[i]] <- bdiag(lapply(mat.list, function(x)x[[i]]))
  }
  return(non_list)
}
