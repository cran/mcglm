#' @title Rotnitzky--Jewell Information Criterion
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description
#' Computes the Rotnitzky--Jewell information criterion (RJC) for objects of
#' class \code{mcglm}. This criterion is based on quasi-likelihood theory
#' and is intended for model assessment in marginal models.
#'
#' @details
#' The RJC is defined using the sensitivity and variability structures of
#' the estimating equations and measures the discrepancy between them.
#' The implementation assumes that the data are correctly ordered such
#' that observations belonging to the same cluster are stored in
#' contiguous rows.
#'
#' \strong{Warning:} This function is restricted to models with a single
#' response variable.
#'
#' @param object An object of class \code{mcglm} representing a fitted
#'   marginal model.
#' @param id An integer or factor vector identifying the clusters. Its
#'   length and ordering must match the number and ordering of the
#'   observations used to fit the model.
#' @param verbose Logical. If \code{TRUE}, the value of the RJC is printed
#'   to the console.
#'
#' @return
#' An invisible list with a single component:
#' \describe{
#'   \item{RJC}{A numeric scalar giving the value of the
#'   Rotnitzky--Jewell information criterion.}
#' }
#'
#' @source
#' Wang, M. (2014). Generalized estimating equations in longitudinal data
#' analysis: A review and recent developments. \emph{Advances in Statistics},
#' 1(1), 1--13.
#'
#' @seealso \code{gof}, \code{plogLik}, \code{pAIC}, \code{pKLIC},
#'   \code{ESS}, \code{GOSHO}
#'
#' @export

RJC <- function(object, id, verbose = TRUE) {
  temp_data <- data.frame(res = object$residuals, id)
  temp_data_group <- split(temp_data, temp_data$id)
  r_rT <- bdiag(lapply(temp_data_group,
                       function(x) {
                         tcrossprod(x[, 1])
                       }))
  D <- bdiag(lapply(object$mu_list, function(x) x$D))
  p1 <- t(D)%*%object$inv_C
  Omega0 <- p1%*%r_rT%*%t(p1)
  Omega1 <- p1%*%D
  Omega <- solve(Omega0)%*%Omega1
  df <- dim(Omega)[1]
  t1 <- (1 - sum(diag(Omega))/df)^2
  t2 <- (1 - sum(diag(Omega^2))/df)^2
  RJC <- sqrt(t1 + t2)
  if (verbose) cat("RJC", RJC)
  return(invisible(list("RJC" = RJC)))
}
