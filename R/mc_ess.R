#' @title Generalized Error Sum of Squares
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description Extract the generalized error sum of squares (ESS) for
#' objects of \code{mcglm} class.
#'
#' @param object an object or a list of objects representing a model
#' of \code{mcglm} class.
#' @param verbose logical. Print or not the ESS value.
#'
#' @return Returns the value of the generalized error sum of squares (ESS).
#' @seealso \code{gof}, \code{plogLik}, \code{pAIC}, \code{pKLIC},
#' \code{GOSHO} and \code{RJC}.
#'
#' @source Bonat, W. H. (2018). Multiple Response Variables Regression
#' Models in R: The mcglm Package. Journal of Statistical Software, 84(4):1--30.
#'
#' @source Wang, M. (2014). Generalized Estimating Equations in Longitudinal Data
#' Analysis: A Review and Recent Developments. Advances in Statistics, 1(1)1--13.
#'
#' @export

ESS <- function(object, verbose = TRUE) {
  if(isa(object, "mcglm")) {
    b <- c(as.matrix(object$observed)) - c(as.matrix(object$fitted))
    ess <- as.numeric(t(b)%*%object$inv_C%*%b)
    df <- length(coef(object)$Estimates)
    df2 <- length(object$observed) - df
    ess <- ess/df2
    if (verbose) cat("ESS", ess)
    return(invisible(list("ESS" = ess)))
  }
  if(isa(object, "list")) {
    Y <- do.call(c,lapply(object, function(x)as.numeric(x$observed)))
    mu <-do.call(c,lapply(object, function(x)as.numeric(x$fitted)))
    b <- Y - mu
    C.list <- lapply(object, function(x)x$C)
    inv_C.list <- lapply(object, function(x)x$inv_C)
    inv_C <- bdiag(inv_C.list)
    ess <- as.numeric(t(b)%*%inv_C%*%b)
    df <- sum(unlist(lapply(object, function(x)length(coef(x)$Estimates))))
    df2 <- length(Y) - df
    ess <- ess/df2
    if (verbose) cat("ESS", ess)
    return(invisible(list("ESS" = ess)))
  }
}
