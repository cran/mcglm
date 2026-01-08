#' @title Auxiliary Function for Block-Diagonal Matrix Construction
#'
#' @description
#' Constructs block-diagonal matrices used in the computation of
#' derivatives of the covariance matrix \eqn{C}. Each input matrix is
#' placed as the non-zero block corresponding to a given response
#' variable, while all remaining blocks are set to zero.
#'
#' @param list_mat
#' A list of matrices to be inserted as non-zero blocks.
#'
#' @param mat_zero
#' A list of zero matrices defining the block-diagonal structure. This
#' object is typically obtained from \code{\link{mc_build_bdiag}}.
#'
#' @param response_number
#' An integer indicating the response variable position where each
#' matrix in \code{list_mat} will be inserted.
#'
#' @details
#' For each matrix in \code{list_mat}, a block-diagonal matrix is created
#' by replacing the block corresponding to \code{response_number} in
#' \code{mat_zero} with the given matrix. The output is a list of
#' block-diagonal matrices with identical structure.
#'
#' @return
#' A list of block-diagonal matrices, each corresponding to one element
#' of \code{list_mat}. The returned object is intended for internal use only.
#'
#' @keywords internal

mc_transform_list_bdiag <- function(list_mat, mat_zero,
                                    response_number) {
    aux.f <- function(x, mat_zero, response_number) {
        mat_zero[[response_number]] <- x
        return(bdiag(mat_zero))
    }
    output <- lapply(list_mat, aux.f, mat_zero = mat_zero,
                     response_number = response_number)
    return(output)
}
