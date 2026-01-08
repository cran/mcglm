#' @title Remove Missing Observations from Matrix Linear Predictor
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description
#' Removes rows and columns corresponding to missing observations from
#' each component of a matrix linear predictor. This function is typically
#' applied after completing the data structure (e.g., via
#' \code{\link[mcglm]{mc_complete_data}}) and before fitting the model, in
#' order to ensure compatibility between the response vector and the
#' covariance-related design matrices.
#'
#' @param matrix_pred A list of square matrices representing the components
#'   of the matrix linear predictor.
#' @param cod An integer vector giving the indices of rows and columns to
#'   be removed. These indices usually correspond to missing observations
#'   in the stacked response vector.
#'
#' @return
#' A list of square matrices of the same length as \code{matrix_pred},
#' where, for each matrix, the rows and columns indexed by \code{cod}
#' have been removed.
#'
#' @details
#' For each matrix \eqn{Z_d} in the matrix linear predictor, the function
#' returns the submatrix obtained by deleting the rows and columns
#' specified in \code{cod}. This operation preserves the symmetry and
#' relative structure of the covariance components while aligning them
#' with the reduced response vector.
#'
#' This is an internal utility function and is not intended to be called
#' directly by end users.
#'
#' @seealso \code{mc_dglm}, \code{mc_ns}, \code{mc_ma}, \code{mc_rw}
#'
#' @keywords internal
#' @export

mc_remove_na <- function(matrix_pred, cod) {
  matrix_pred_temp <- list()
  for(i in 1:length(matrix_pred)) {
    matrix_pred_temp[[i]] <- matrix_pred[[i]][-cod,-cod]
  }
  return(matrix_pred_temp)
}
