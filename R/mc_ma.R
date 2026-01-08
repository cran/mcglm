#' @title  Moving Average Model Structure
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description Builds components of the matrix linear predictor associated
#'   with moving average (MA) covariance structures. This function is mainly
#'   intended for longitudinal data analysis, but can also be used for
#'   time series data
#'
#' @param id name of the column (string) containing the subject index.
#' Note that this structure was designed to deal with longitudinal data.
#' For times series data use the same \code{id} for all observations
#' (one unit sample).
#' @param time name of the column (string) containing the index indicating
#' the time.
#' @param data data set.
#' @param order An integer specifying the order of the moving average process.
#'
#' @details This function was primarily designed for longitudinal data,
#'   but it can also be used for time series analysis. In this case, the
#'   \code{id} argument should contain a single identifier, representing
#'   one observational unit. Internally, the function constructs block-diagonal
#'   band matrices using \code{\link[Matrix]{bandSparse}}.
#'
#' @source Bonat, W. H. (2018). Multiple Response Variables Regression
#' Models in R: The mcglm Package. Journal of Statistical Software, 84(4):1--30.
#'
#' @return A list with the following component:
#' \describe{
#'   \item{Z1}{A sparse matrix of class \code{nsCMatrix} representing the
#'   moving average component of the matrix linear predictor. The matrix
#'   has dimension equal to the total number of observations and is
#'   constructed as a block-diagonal matrix, with one block per subject
#'   (or time series), each block encoding a moving average structure of
#'   the specified order.}
#' }
#'
#' @seealso \code{mc_id}, \code{mc_dist}, \code{mc_car},
#' \code{mc_rw} and \code{mc_mixed}.
#'
#' @examples
#' id <- rep(1:2, each = 4)
#' time <- rep(1:4, 2)
#' data <- data.frame("id" = id, "time" = time)
#' mc_ma(id = "id", time = "time", data = data, order = 1)
#' mc_ma(id = "id", time = "time", data = data, order = 2)
#'
#' @export

mc_ma <- function(id, time, data, order = 1) {
  mc_ma_aux <- function(n, order) {
    output <- bandSparse(n, n, k = order, symmetric = TRUE)
    return(output)
  }
  data$id2 <- 1:dim(data)[1]
  data <- data[order(data[[id]]),]
  data$id3 <- 1:dim(data)[1]
  Z1.list <- list()
  data.id <- split(data, data[[id]], drop = TRUE)
  DD <- sum(abs(diff(do.call(c,lapply(data.id, function(x)dim(x)[1])))))
  if( DD != 0) {
    stop("Model requires equal number of observations by id. \n")
  }
  for(i in 1:length(data.id)) {
    NN <- dim(data.id[[i]])[1]
    ordem <- as.numeric(data.id[[i]][[time]])
    Z1.list[[i]] <- mc_ma_aux(n = NN, order = order)[ordem,ordem]
  }
  ordem2 <- order(data$id2)
  output <- list("Z1" = bdiag(Z1.list)[ordem2,ordem2])
  return(output)
}
