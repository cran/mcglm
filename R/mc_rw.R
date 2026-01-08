#' @title Random Walk Model Structure
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description
#' Constructs the components of the matrix linear predictor associated
#' with random walk (RW) models for longitudinal or time series data.
#' The user may specify the order of the random walk process.
#'
#' @details
#' This function builds sparse precision matrix components corresponding
#' to random walk structures of a given order. It is primarily intended
#' for longitudinal data indexed by a subject identifier and a time
#' variable. For pure time series data, the same \code{id} value should
#' be used for all observations. When \code{proper = TRUE}, the precision
#' structure is decomposed into diagonal and off-diagonal components.
#'
#' @param id A character string giving the name of the column in
#'   \code{data} that identifies subjects or clusters.
#' @param time A character string giving the name of the column in
#'   \code{data} that indexes time or ordering within each subject.
#' @param data A data frame containing the variables specified in
#'   \code{id} and \code{time}.
#' @param order A positive integer specifying the order of the random
#'   walk model.
#' @param proper Logical indicating whether a proper random walk
#'   specification should be used.
#'
#' @return
#' If \code{proper = FALSE}, a list with a single component:
#' \describe{
#'   \item{Z1}{A sparse matrix of class \code{dgCMatrix} representing the
#'   random walk precision structure.}
#' }
#' If \code{proper = TRUE}, a list with two components:
#' \describe{
#'   \item{Z1}{A sparse diagonal matrix of class \code{dgCMatrix}.}
#'   \item{Z2}{A sparse off-diagonal matrix of class \code{dgCMatrix}.}
#' }
#' The matrices are ordered consistently with the original data.
#'
#' @source
#' Bonat, W. H. (2018). Multiple Response Variables Regression Models in R:
#' The mcglm Package. \emph{Journal of Statistical Software}, 84(4), 1--30.
#'
#' @seealso \code{mc_id}, \code{mc_dist}, \code{mc_car},
#'   \code{mc_ma}, \code{mc_mixed}, \code{mc_compute_rho}
#'
#' @examples
#' id <- rep(1:2, each = 4)
#' time <- rep(1:4, 2)
#' data <- data.frame(id = id, time = time)
#' mc_rw(id = "id", time = "time", data = data, order = 1, proper = FALSE)
#' mc_rw(id = "id", time = "time", data = data, order = 1, proper = TRUE)
#' mc_rw(id = "id", time = "time", data = data, order = 2, proper = TRUE)
#'
#' @export



mc_rw <- function(id, time, data, order = 1, proper = FALSE) {
    mc_rw_aux <- function(n, order, proper = TRUE) {
        U = Matrix(diff(diag(n), diff = order),sparse = TRUE)
        Z <- forceSymmetric(t(U)%*%U)
        if(proper == TRUE) {
            Z2 <- Z
            Z1 <- Diagonal(dim(Z)[1], diag(Z))
            diag(Z2) <- 0
            output <- list("Z1" = Z1, "Z2" = Z2)
        } else {
            output <- list("Z1" = Z)
        }
        return(output)
    }
    data$id2 <- 1:dim(data)[1]
    data <- data[order(data[[id]]),]
    data$id3 <- 1:dim(data)[1]
    if(proper == TRUE) {
        Z1.list <- list()
        Z2.list <- list()
    } else {
        Z1.list <- list()
    }
    data.id <- split(data, data[[id]], drop = TRUE)
    DD <- sum(abs(diff(do.call(c,lapply(data.id, function(x)dim(x)[1])))))
    if( DD != 0) {
      stop("Model requires equal number of observations by id. \n")
    }
    for(i in 1:length(data.id)) {
        n <- dim(data.id[[i]])[1]
        ordem <- as.numeric(data.id[[i]][[time]])
        if(proper == TRUE) {
            ZZ <- mc_rw_aux(n = n, order = order, proper = TRUE)
            Z1.list[[i]] <- ZZ$Z1[ordem,ordem]
            Z2.list[[i]] <- ZZ$Z2[ordem,ordem]
        } else {
            ZZ <- mc_rw_aux(n = n, order = order, proper = FALSE)
            Z1.list[[i]] <- ZZ$Z1[ordem,ordem]
        }
    }
    if(proper == TRUE) {
        ordem2 <- order(data$id2)
        output <- list("Z1" = bdiag(Z1.list)[ordem2,ordem2],
                       "Z2" = bdiag(Z2.list)[ordem2,ordem2])
    } else {
        ordem2 <- order(data$id2)
        output <- list("Z1" = bdiag(Z1.list)[ordem2,ordem2])
    }
    return(output)
}
