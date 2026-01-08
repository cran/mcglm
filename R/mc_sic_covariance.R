#' @title Score Information Criterion for Covariance Components
#'
#' @description
#' Computes the Score Information Criterion (SIC) for covariance
#' components of a fitted \code{mcglm} object. The SIC-covariance is used
#' to select components of the matrix linear predictor and can be
#' employed in stepwise selection procedures.
#'
#' @param object
#' An object of class \code{mcglm}.
#'
#' @param scope
#' A list of matrices to be tested for inclusion in the matrix linear
#' predictor.
#'
#' @param idx
#' An integer vector indicating which matrices in \code{scope} belong to
#' the same effect. This is useful when more than one matrix represents
#' a single covariance component.
#'
#' @param data
#' A data frame containing all variables involved in the model.
#'
#' @param penalty
#' A numeric penalty term applied to the SIC (default is 2).
#'
#' @param response
#' An integer indicating the response variable for which the
#' SIC-covariance is computed.
#'
#' @param weights
#' An optional numeric vector of weights used in model fitting. If not
#' provided, unit weights are assumed.
#'
#' @details
#' The SIC-covariance is computed using the Pearson estimating function.
#' For each group of matrices defined by \code{idx}, a score-based test
#' statistic is calculated to assess the contribution of the associated
#' covariance components, penalized by model complexity.
#'
#' @return
#' A data frame with the following columns:
#' \describe{
#'   \item{SIC}{Score Information Criterion value.}
#'   \item{df}{Degrees of freedom associated with the test.}
#'   \item{df_total}{Total number of covariance parameters in the extended model.}
#'   \item{Tu}{Score-based test statistic.}
#'   \item{Chisq}{Reference chi-squared quantile with 95\% confidence level.}
#' }
#'
#' @references
#' Bonat, W. H., et al. (2016). Modelling the covariance structure in
#' marginal multivariate count models: Hunting in Bioko Island.
#' \emph{Journal of Agricultural, Biological and Environmental Statistics},
#' 22(4), 446--464.
#'
#' Bonat, W. H. (2018). Multiple Response Variables Regression Models in R:
#' The mcglm Package. \emph{Journal of Statistical Software}, 84(4), 1--30.
#'
#' @seealso
#' \code{\link{mc_sic}}
#'
#' @examples
#' set.seed(123)
#' SUBJECT <- gl(10, 10)
#' y <- rnorm(100)
#' data <- data.frame(y, SUBJECT)
#'
#' Z0 <- mc_id(data)
#' Z1 <- mc_mixed(~ 0 + SUBJECT, data = data)
#'
#' fit0 <- mcglm(
#'   linear_pred = c(y ~ 1),
#'   matrix_pred = list(Z0),
#'   data = data
#' )
#'
#' mc_sic_covariance(
#'   fit0,
#'   scope = Z1,
#'   idx = 1,
#'   data = data,
#'   response = 1
#' )
#'
#' @export

mc_sic_covariance <- function(object, scope, idx, data, penalty = 2,
                              response, weights) {
  if(missing(weights)) {
    weights <- rep(1, dim(object$C)[1])
  }
  W <- Diagonal(length(weights), 1)
    SIC <- c()
    df <- c()
    df_total <- c()
    TU <- c()
    QQ <- c()
    n_terms <- length(unique(idx))
    for (j in 1:n_terms) {
        tau <- coef(object, type = "tau", response = response)$Estimates
        n_tau <- length(tau)
        n_tau_new <- length(idx[idx == j])
        list_tau_new <- list(c(tau, rep(0, n_tau_new)))
        n_tau_total <- n_tau + n_tau_new
        if (object$power_fixed[[response]]) {
            list_power <- object$list_initial$power
        } else {
            list_power <- list(coef(object, type = "power",
                                    response = response)$Estimates)
            n_tau_total <- n_tau_total + 1
            n_tau <- n_tau + 1
        }
        list_Z_new <- list(c(object$matrix_pred[[response]],
                             scope[idx == j]))
        if (length(object$mu_list) == 1) {
            rho <- 0
        } else {
            rho <- coef(object, type = "correlation")$Estimates
        }
        Cfeatures <- mc_build_C(list_mu = object$mu_list,
                                list_Ntrial = object$Ntrial, rho = rho,
                                list_tau = list_tau_new,
                                list_power = list_power,
                                list_Z = list_Z_new,
                                list_sparse = object$sparse,
                                list_variance = object$variance,
                                list_covariance = object$covariance,
                                list_power_fixed = object$power_fixed,
                                compute_C = TRUE)
        temp_score <- mc_pearson(y_vec = as.numeric(object$observed),
                                 mu_vec = as.numeric(object$mu_list[[response]]$mu),
                                 Cfeatures = Cfeatures, correct = FALSE,
                                 compute_variability = TRUE, W = W)

        J <- temp_score$Sensitivity
        Sigma <- temp_score$Variability
        Sigma22 <- Sigma[c(n_tau + 1):n_tau_total,
                         c(n_tau + 1):n_tau_total]
        J21 <- J[c(n_tau + 1):n_tau_total, 1:n_tau]
        J11 <- solve(J[1:n_tau, 1:n_tau])
        Sigma12 <- Sigma[1:n_tau, c(n_tau + 1):n_tau_total]
        Sigma21 <- Sigma[c(n_tau + 1):n_tau_total, 1:n_tau]
        J12 <- J[1:n_tau, c(n_tau + 1):n_tau_total]
        Sigma11 <- Sigma[1:n_tau, 1:n_tau]

        V2 <- Sigma22 - J21 %*% J11 %*% Sigma12 - Sigma21 %*%
            J11 %*% J12 + J21 %*% J11 %*% Sigma11 %*% J11 %*%
            J12
        TU[j] <- t(temp_score$Score[c(n_tau + 1):n_tau_total] %*%
                   solve(V2) %*%
                   temp_score$Score[c(n_tau + 1):n_tau_total])
        df[j] <- n_tau_new
        SIC[j] <- -as.numeric(TU[j]) + penalty * n_tau_total
        QQ[j] <- qchisq(0.95, df = n_tau_new)
        df_total[j] <- n_tau_total
        #print(j)
    }
    output <- data.frame(SIC = SIC, df = df, df_total = df_total,
                         Tu = TU, Chisq = QQ)
    return(output)
}
