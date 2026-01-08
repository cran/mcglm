#' @title Score Information Criterion for Regression Components
#'
#' @description
#' Computes the Score Information Criterion (SIC) for regression
#' components of a fitted \code{mcglm} object. The SIC can be used for
#' selecting covariates in the linear predictor and supports stepwise
#' selection procedures.
#'
#' @param object
#' An object of class \code{mcglm}.
#'
#' @param scope
#' A character vector with the names of covariates to be tested for
#' inclusion in the linear predictor.
#'
#' @param data
#' A data frame containing all variables involved in the model.
#'
#' @param response
#' An integer indicating the response variable for which the SIC is
#' computed.
#'
#' @param penalty
#' A numeric penalty term applied to the SIC (default is 2).
#'
#' @param weights
#' An optional numeric vector of weights used in model fitting. If not
#' provided, unit weights are assumed.
#'
#' @details
#' The SIC is computed using the quasi-score function associated with the
#' regression parameters. For each candidate covariate in \code{scope},
#' the method evaluates its contribution via a score-based test statistic
#' and applies a penalty for model complexity.
#'
#' @return
#' A data frame with the following columns:
#' \describe{
#'   \item{SIC}{Score Information Criterion value.}
#'   \item{Covariates}{Name of the candidate covariate.}
#'   \item{df}{Degrees of freedom associated with the test.}
#'   \item{df_total}{Total number of regression parameters in the extended model.}
#'   \item{Tu}{Score-based test statistic.}
#'   \item{Chisq}{Reference chi-squared quantile with 95\% confidence level.}
#' }
#'
#' @references
#' Bonat, W. H. (2018). Multiple Response Variables Regression Models in R:
#' The mcglm Package. \emph{Journal of Statistical Software}, 84(4), 1--30.
#'
#' Bonat, W. H., et al. (2016). Modelling the covariance structure in
#' marginal multivariate count models: Hunting in Bioko Island.
#' \emph{Journal of Agricultural, Biological and Environmental Statistics},
#' 22(4), 446--464.
#'
#' @seealso
#' \code{\link{mc_sic_covariance}}
#'
#' @examples
#' set.seed(123)
#' x1 <- runif(100, -1, 1)
#' x2 <- gl(2, 50)
#' beta <- c(5, 0, 3)
#' X <- model.matrix(~ x1 + x2)
#' y <- rnorm(100, mean = X %*% beta, sd = 1)
#' data <- data.frame(y, x1, x2)
#'
#' Z0 <- mc_id(data)
#' fit0 <- mcglm(
#'   linear_pred = c(y ~ 1),
#'   matrix_pred = list(Z0),
#'   data = data
#' )
#'
#' mc_sic(fit0, scope = c("x1", "x2"), data = data, response = 1)
#'
#' @export

mc_sic <- function(object, scope, data, response, penalty = 2, weights) {
  if(missing(weights)) {
    weights <- rep(1, dim(object$C)[1])
  }
  W <- Diagonal(length(weights), 1)
    SIC <- c()
    df <- c()
    df_total <- c()
    TU <- c()
    QQ <- c()
    for (i in 1:length(scope)) {
        ini_formula <- object$linear_pred[[response]]
        ext_formula <- as.formula(
            paste("~", paste(ini_formula[3], scope[i], sep = "+")))
        md <- model.frame(object$linear_pred[[response]], data = data)
        Y <- model.response(md)
        ini_beta <- coef(object, type = "beta",
                         response = response)$Estimates
        ext_X <- model.matrix(ext_formula, data = data)
        n_beta <- dim(ext_X)[2]
        n_ini_beta <- length(ini_beta)
        ext_beta <- c(ini_beta, rep(0, n_beta - n_ini_beta))
        n_total_beta <- length(ext_beta)
        mu_temp <- mc_link_function(beta = ext_beta, X = ext_X,
                                    offset = NULL,
                                    link = object$link[[response]])
        score_temp <- mc_quasi_score(D = mu_temp$D,
                                     inv_C = object$inv_C, y_vec = Y,
                                     mu_vec = mu_temp$mu, W = W)
        S11 <- score_temp$Variability[1:n_ini_beta, 1:n_ini_beta]
        S22 <- score_temp$Variability[c(n_ini_beta + 1):n_total_beta,
                                      c(n_ini_beta + 1):n_total_beta]
        S12 <- score_temp$Variability[1:n_ini_beta,
                                      c(n_ini_beta + 1):n_total_beta]
        S21 <- score_temp$Variability[c(n_ini_beta + 1):n_total_beta,
                                      1:n_ini_beta]
        VB <- S22 - S21 %*% solve(S11) %*% S12
        Tu <- t(score_temp$Score[c(n_ini_beta + 1):n_total_beta]) %*%
            solve(VB) %*%
            score_temp$Score[c(n_ini_beta + 1):n_total_beta]
        df[i] <- n_beta - n_ini_beta
        SIC[i] <- -as.numeric(Tu) + penalty * n_beta
        df_total[i] <- n_beta
        TU[i] <- as.numeric(Tu)
        QQ[i] <- qchisq(0.95, df = df[i])
    }
    output <- data.frame(SIC = SIC, Covariates = scope, df = df,
                         df_total = df_total, Tu = TU, Chisq = QQ)
    return(output)
}
