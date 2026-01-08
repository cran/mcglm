#' @title Wald Tests for Fixed Effects in mcglm Models
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description
#' Performs Wald chi-square tests for assessing the significance of
#' fixed-effect terms in the linear predictors of an \code{mcglm} model.
#' The tests are conducted separately for each response variable and
#' are particularly useful for joint hypothesis testing of regression
#' coefficients associated with categorical covariates with more than
#' two levels. This function is not intended for model comparison.
#'
#' @param object An object of class \code{mcglm}, typically the result of
#'   a call to \code{\link{mcglm}}.
#' @param ... Additional arguments. Currently ignored.
#' @param verbose
#' Logical indicating whether the Wald test results should be printed
#' to the console. If \code{FALSE}, the function silently returns the
#' results.
#' @return
#' A list of data frames, one for each response variable. Each data frame
#' contains the results of Wald chi-square tests for the fixed-effect
#' terms in the corresponding linear predictor, with the following
#' columns:
#' \describe{
#'   \item{Covariate}{Name of the covariate or model term tested.}
#'   \item{Chi.Square}{Value of the Wald chi-square statistic.}
#'   \item{Df}{Degrees of freedom associated with the test.}
#'   \item{p.value}{P-value of the Wald test.}
#' }
#' The returned object is invisible and is primarily intended for
#' programmatic use.
#'
#' @details
#' The Wald tests are computed using the observed covariance matrix of
#' the regression parameter estimates. For each response variable, joint
#' tests are performed for sets of parameters corresponding to the same
#' model term, as defined by the design matrix.
#'
#' @seealso \code{\link{summary.mcglm}}, \code{\link{coef.mcglm}},
#'   \code{\link{vcov.mcglm}}
#'
#' @examples
#' x1 <- seq(-1, 1, length.out = 100)
#' x2 <- gl(5, 20)
#' beta <- c(5, 0, -2, -1, 1, 2)
#' X <- model.matrix(~ x1 + x2)
#' set.seed(123)
#' y <- rnorm(100, mean = X %*% beta, sd = 1)
#' data <- data.frame(y = y, x1 = x1, x2 = x2)
#' fit <- mcglm(c(y ~ x1 + x2), list(mc_id(data)), data = data)
#' anova(fit)
#'
#' @method anova mcglm
#' @export

anova.mcglm <- function(object, ..., verbose = TRUE) {

  if (!inherits(object, "mcglm")) {
    stop("object must be of class 'mcglm'", call. = FALSE)
  }

  n_resp <- length(object$mu_list)
  n_beta <- vapply(object$list_X, ncol, integer(1))

  idx.list <- lapply(seq_len(n_resp), function(i) {
    rep(i, n_beta[i])
  })

  vv <- vcov(object)
  idx.vec <- c(unlist(idx.list),
               rep(0, ncol(vv) - length(unlist(idx.list))))

  temp.vcov <- vector("list", n_resp)
  temp.beta <- vector("list", n_resp)

  for (i in seq_len(n_resp)) {
    idx.id <- idx.vec == i
    temp.vcov[[i]] <- vv[idx.id, idx.id, drop = FALSE]
    temp.beta[[i]] <-
      coef(object, type = "beta", response = i)$Estimates
  }

  output <- vector("list", n_resp)

  for (i in seq_len(n_resp)) {

    assign_idx <- attr(object$list_X[[i]], "assign")
    term_names <- colnames(object$list_X[[i]])

    if (term_names[1] == "(Intercept)") {
      assign_idx <- assign_idx[-1]
      term_names <- term_names[-1]
      temp.beta[[i]] <- temp.beta[[i]][-1]
      temp.vcov[[i]] <- temp.vcov[[i]][-1, -1, drop = FALSE]
    }

    n_terms <- length(unique(assign_idx))
    res_i <- vector("list", n_terms)

    for (j in seq_len(n_terms)) {

      idx.term <- assign_idx == j
      beta_j <- temp.beta[[i]][idx.term]
      vcov_j <- temp.vcov[[i]][idx.term, idx.term, drop = FALSE]

      inv_vcov_j <- solve(vcov_j, tol = .Machine$double.eps)

      chi_sq <- as.numeric(t(beta_j) %*% inv_vcov_j %*% beta_j)
      df_j <- length(beta_j)

      res_i[[j]] <- data.frame(
        Covariate = term_names[idx.term][1],
        Chi.Square = round(chi_sq, 4),
        Df = df_j,
        p.value = round(pchisq(chi_sq, df_j, lower.tail = FALSE), 4)
      )
    }

    output[[i]] <- do.call(rbind, res_i)
  }

  if (verbose) {
    cat("Wald test for fixed effects\n")
    for (i in seq_len(n_resp)) {
      cat("Call: ")
      print(object$linear_pred[[i]])
      cat("\n")
      print(output[[i]])
      cat("\n")
    }
  }

  return(invisible(output))
}

#' @title Model Coefficients
#' @name coef.mcglm
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description
#' Extract regression, dispersion and correlation parameter estimates
#' from objects of class \code{mcglm}.
#'
#' @param object
#' An object of class \code{mcglm}.
#'
#' @param std.error
#' Logical indicating whether standard errors should be returned
#' alongside the parameter estimates. Default is \code{FALSE}.
#'
#' @param response
#' Integer vector indicating for which response variables the
#' coefficients should be returned. If \code{NA}, coefficients for
#' all response variables are returned.
#'
#' @param type
#' Character vector specifying which type of coefficients should be
#' returned. Possible values are \code{"beta"}, \code{"tau"},
#' \code{"power"} and \code{"correlation"}.
#'
#' @param ...
#' Additional arguments. Currently ignored and included for
#' compatibility with the generic \code{\link[stats]{coef}} function.
#'
#' @return
#' A \code{data.frame} with one row per parameter, containing:
#' \itemize{
#'   \item \code{Estimates}: parameter estimates;
#'   \item \code{Std.error}: standard errors (if requested);
#'   \item \code{Parameters}: parameter names;
#'   \item \code{Type}: parameter type;
#'   \item \code{Response}: response variable index.
#' }
#'
#' @method coef mcglm
#' @export

coef.mcglm <- function(object, std.error = FALSE,
                       response = c(NA, 1:length(object$beta_names)),
                       type = c("beta", "tau", "power", "correlation"),
                       ...) {
    n_resp <- length(object$beta_names)
    cod_beta <- list()
    cod_power <- list()
    cod_tau <- list()
    type_beta <- list()
    type_power <- list()
    type_tau <- list()
    resp_beta <- list()
    resp_power <- list()
    resp_tau <- list()
    response_for <- 1:n_resp
    for (i in response_for) {
        cod_beta[[i]] <- paste0(
            paste0("beta", i), 0:c(object$Information$n_betas[[i]] - 1))
        type_beta[[i]] <- rep("beta", length(cod_beta[[i]]))
        resp_beta[[i]] <- rep(response_for[i], length(cod_beta[[i]]))
        if (object$Information$n_power[[i]] != 0 |
            object$power_fixed[[i]] == FALSE) {
            cod_power[[i]] <- paste0(
                paste0("power", i), 1:object$Information$n_power[[i]])
            type_power[[i]] <- rep("power",
                                   length(cod_power[[i]]))
            resp_power[[i]] <- rep(response_for[i],
                                   length(cod_power[[i]]))
        }
        if (object$Information$n_power[[i]] == 0) {
            cod_power[[i]] <- rep(1, 0)
            type_power[[i]] <- rep(1, 0)
            resp_power[[i]] <- rep(1, 0)
        }
        cod_tau[[i]] <- paste0(
            paste0("tau", i), 1:object$Information$n_tau[[i]])
        type_tau[[i]] <- rep("tau", length(cod_tau[[i]]))
        resp_tau[[i]] <- rep(response_for[i], length(cod_tau[[i]]))
    }
    rho_names <- c()
    if (n_resp != 1) {
        combination <- combn(n_resp, 2)
        for (i in 1:dim(combination)[2]) {
            rho_names[i] <- paste0(
                paste0("rho", combination[1, i]), combination[2, i])
        }
    }
    type_rho <- rep("correlation", length(rho_names))
    resp_rho <- rep(NA, length(rho_names))
    cod <- c(do.call(c, cod_beta), rho_names,
             do.call(c, Map(c, cod_tau)))
    type_cod <- c(do.call(c, type_beta), type_rho,
                  do.call(c, Map(c, type_tau)))
    response_cod <- c(do.call(c, resp_beta), resp_rho,
                      do.call(c, Map(c, resp_tau)))

    if (length(cod_power) != 0) {
        cod <- c(do.call(c, cod_beta), rho_names,
                 do.call(c, Map(c, cod_power, cod_tau)))
        type_cod <- c(do.call(c, type_beta), type_rho,
                      do.call(c, Map(c, type_power, type_tau)))
        response_cod <- c(do.call(c, resp_beta), resp_rho,
                          do.call(c, Map(c, resp_power, resp_tau)))
    }

    Estimates <- c(object$Regression, object$Covariance)
    coef_temp <- data.frame(
        Estimates = Estimates,
        Parameters = cod,
        Type = type_cod,
        Response = response_cod)
    if (std.error == TRUE) {
        coef_temp <- data.frame(
            Estimates = Estimates,
            Std.error = sqrt(diag(object$vcov)),
            Parameters = cod, Type = type_cod,
            Response = response_cod)
    }
    output <- coef_temp[
        which(coef_temp$Response %in% response &
                  coef_temp$Type %in% type), ]
    return(output)
}


#' @title Confidence Intervals for Model Parameters
#' @name confint.mcglm
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description
#' Computes Wald-type confidence intervals for parameter estimates
#' from a fitted \code{mcglm} model, based on asymptotic normality.
#'
#' @param object
#' A fitted object of class \code{mcglm}.
#'
#' @param parm
#' Optional specification of parameters for which confidence intervals
#' are required. Can be a numeric vector of indices or a character
#' vector of parameter names. If omitted, confidence intervals for all
#' parameters are returned.
#'
#' @param level
#' Numeric value giving the confidence level. Must be between 0 and 1.
#' Default is \code{0.95}.
#'
#' @param ...
#' Additional arguments. Currently ignored and included for
#' compatibility with the generic
#' \code{\link[stats]{confint}} function.
#'
#' @return
#' A numeric matrix with two columns corresponding to the lower and
#' upper confidence limits. Rows correspond to model parameters.
#'
#' @method confint mcglm
#' @export

confint.mcglm <- function(object, parm, level = 0.95, ...) {

  if (!is.numeric(level) || length(level) != 1 ||
      level <= 0 || level >= 1) {
    stop("'level' must be a numeric value between 0 and 1.")
  }

  coef_tab <- coef(object, std.error = TRUE)

  if (is.null(coef_tab$Std.error)) {
    stop("Standard errors are required to compute confidence intervals.")
  }

  n_par <- nrow(coef_tab)

  if (missing(parm)) {
    parm_idx <- seq_len(n_par)
  } else if (is.numeric(parm)) {
    if (any(parm < 1 | parm > n_par)) {
      stop("Invalid parameter indices in 'parm'.")
    }
    parm_idx <- parm
  } else if (is.character(parm)) {
    if (any(!parm %in% coef_tab$Parameters)) {
      stop("Invalid parameter names in 'parm'.")
    }
    parm_idx <- match(parm, coef_tab$Parameters)
  } else {
    stop("'parm' must be numeric or character.")
  }

  alpha <- (1 - level) / 2
  z <- stats::qnorm(c(alpha, 1 - alpha))

  ci <- coef_tab$Estimates + coef_tab$Std.error %o% z
  colnames(ci) <- c("Lower", "Upper")
  rownames(ci) <- coef_tab$Parameters

  ci[parm_idx, , drop = FALSE]
}

#' @title Fitted Values
#' @name fitted.mcglm
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description
#' Extract fitted (mean) values from a fitted \code{mcglm} model.
#' For multivariate responses, fitted values are returned in matrix
#' form, with one column per response variable.
#'
#' @param object
#' A fitted object of class \code{mcglm}.
#'
#' @param ...
#' Additional arguments. Currently ignored and included for
#' compatibility with the generic
#' \code{\link[stats]{fitted}} function.
#'
#' @return
#' A numeric matrix of fitted values. Rows correspond to observations
#' and columns correspond to response variables.
#'
#' @method fitted mcglm
#' @export

fitted.mcglm <- function(object, ...) {

  if (is.null(object$fitted)) {
    stop("No fitted values found in 'object'.")
  }

  n_resp <- length(object$beta_names)
  n_obs  <- object$n_obs

  if (length(object$fitted) != n_obs * n_resp) {
    stop("Length of fitted values is inconsistent with model dimensions.")
  }

  fitted_mat <- Matrix::Matrix(
    object$fitted,
    nrow = n_obs,
    ncol = n_resp
  )

  colnames(fitted_mat) <- paste0("Response_", seq_len(n_resp))
  rownames(fitted_mat) <- seq_len(n_obs)

  fitted_mat
}

#' @title Diagnostic Plots for mcglm Objects
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description
#' Produces diagnostic plots for fitted \code{mcglm} objects and returns
#' the data used to generate the plots. Available diagnostics include:
#' \itemize{
#'   \item \code{"residuals"}: plots of fitted values vs. Pearson residuals
#'         and Q-Q plots of residuals;
#'   \item \code{"algorithm"}: plots of regression and covariance parameters
#'         and quasi-score iterations to inspect algorithm convergence;
#'   \item \code{"partial_residuals"}: partial residual plots for each
#'         covariate in the model.
#' }
#'
#' @param x
#' A fitted object of class \code{mcglm}.
#'
#' @param type
#' Character string specifying the type of diagnostic plot to produce.
#' One of \code{"residuals"}, \code{"algorithm"}, or \code{"partial_residuals"}.
#'
#' @param plot_graphics
#' Logical; if \code{TRUE} (default), the function produces the plots.
#' If \code{FALSE}, the function only returns the data used for plotting,
#' which can be used to produce custom plots.
#'
#' @param ...
#' Additional arguments; currently ignored. Included for compatibility
#' with the generic \code{\link[graphics]{plot}} function.
#'
#' @return
#' A list containing the data used to generate the plots:
#' \itemize{
#'   \item For \code{"residuals"}: \code{fitted} values and \code{residuals}.
#'   \item For \code{"algorithm"}: \code{IterationRegression},
#'         \code{IterationCovariance}, \code{ScoreRegression},
#'         \code{ScoreCovariance}.
#'   \item For \code{"partial_residuals"}: \code{residuals}, \code{list_beta},
#'         and \code{partial_residuals} containing x and y values for each covariate.
#' }
#'
#' @seealso
#' \code{\link[stats]{residuals}}, \code{\link[stats]{fitted}},
#' \code{\link[graphics]{plot}}
#'
#' @examples
#' library(mcglm)
#' set.seed(123)
#' mydata <- data.frame(y = rnorm(10), x1 = rnorm(10),
#'                     x2 = rbinom(10, size = 1, prob = 0.5))
#' Z0 <- mc_id(mydata)
#' fit <- mcglm(c(y ~ x1 + x2), matrix_pred = list(Z0), data = mydata)
#' # Produce plots and get data
#' diag_data <- plot(fit, type = "residuals")
#' # Only get data without plotting
#' diag_data <- plot(fit, type = "partial_residuals", plot_graphics = FALSE)
#'
#' @method plot mcglm
#' @export


plot.mcglm <- function(x,
                       type = c("residuals", "algorithm", "partial_residuals"),
                       plot_graphics = TRUE,
                       ...) {

  type <- match.arg(type)
  object <- x
  n_resp <- length(object$beta_names)

  out <- list()  # Lista que vai armazenar os dados usados nos gráficos

  if (type == "residuals") {
    res <- residuals(object, type = "pearson")
    fit <- fitted(object)

    out$residuals <- res
    out$fitted <- fit

    if (plot_graphics) {
      oldpar <- graphics::par(mfrow = c(2, n_resp), mar = c(2.6, 2.5, 0.1, 0.1),
                              mgp = c(1.6, 0.6, 0))
      on.exit(graphics::par(oldpar), add = TRUE)

      for (i in seq_len(n_resp)) {
        graphics::plot(fit[, i], res[, i],
                       xlab = "Fitted values",
                       ylab = "Pearson residuals")
        lo <- stats::loess.smooth(fit[, i], res[, i])
        graphics::lines(lo$x, lo$y)
        stats::qqnorm(res[, i])
        stats::qqline(res[, i])
      }
    }
  }

  if (type == "algorithm") {
    n_iter <- sum(!is.na(object$IterationCovariance[, 1]))
    idx <- seq_len(n_iter)

    out$IterationRegression <- object$IterationRegression[idx, , drop = FALSE]
    out$IterationCovariance <- object$IterationCovariance[idx, , drop = FALSE]
    out$ScoreRegression <- object$ScoreRegression[idx, , drop = FALSE]
    out$ScoreCovariance <- object$ScoreCovariance[idx, , drop = FALSE]

    if (plot_graphics) {
      oldpar <- graphics::par(mfrow = c(2, 2), mar = c(2.6, 2.5, 0.1, 0.1),
                              mgp = c(1.6, 0.6, 0))
      on.exit(graphics::par(oldpar), add = TRUE)

      graphics::matplot(out$IterationRegression, type = "l", lty = 2,
                        xlab = "Iterations", ylab = "Regression")
      graphics::matplot(out$IterationCovariance, type = "l", lty = 2,
                        xlab = "Iterations", ylab = "Covariance")
      graphics::matplot(out$ScoreRegression, type = "l", lty = 2,
                        xlab = "Iterations", ylab = "Quasi-score Regression")
      graphics::matplot(out$ScoreCovariance, type = "l", lty = 2,
                        xlab = "Iterations", ylab = "Quasi-score Covariance")
    }
  }

  if (type == "partial_residuals") {
    if (!exists("mc_updateBeta", mode = "function")) {
      stop("Function 'mc_updateBeta' not found.")
    }

    list_beta <- mc_updateBeta(
      list_initial = object$list_initial,
      betas = object$Regression,
      n_resp = n_resp,
      information = object$Information
    )

    res <- residuals(object, type = "pearson")
    out$residuals <- res
    out$list_beta <- list_beta

    partial_list <- vector("list", n_resp)

    for (i in seq_len(n_resp)) {
      comp_X <- as.matrix(object$list_X[[i]]) * as.numeric(list_beta$regression[[i]])
      n_cov <- ncol(comp_X)
      partial_list[[i]] <- vector("list", n_cov - 1)

      if (plot_graphics) {
        oldpar <- graphics::par(mfrow = c(1, max(1, n_cov - 1)),
                                mar = c(2.6, 2.5, 0.5, 0.5),
                                mgp = c(1.6, 0.6, 0))
        on.exit(graphics::par(oldpar), add = TRUE)
      }

      for (j in 2:n_cov) {
        p_res <- comp_X[, j] + res[, i]
        partial_list[[i]][[j - 1]] <- list(
          x = object$list_X[[i]][, j],
          y = p_res,
          xlab = object$beta_names[[i]][j],
          ylab = "Partial residuals"
        )
        if (plot_graphics) {
          graphics::plot(object$list_X[[i]][, j], p_res,
                         xlab = object$beta_names[[i]][j],
                         ylab = "Partial residuals")
        }
      }
    }
    out$partial_residuals <- partial_list
  }

  return(out)
}

#' @title Print Method for mcglm Objects
#' @name print.mcglm
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description
#' Prints a concise summary of a fitted \code{mcglm} object, including
#' the model call, link and variance functions, regression coefficients
#' and dispersion parameters for each response variable.
#'
#' @param x
#' A fitted object of class \code{mcglm}, typically returned by
#' \code{mcglm()}.
#'
#' @param ...
#' Further arguments passed to or from other methods.
#'
#' @return
#' No return value, called for its side effects.
#'
#' @seealso
#' \code{\link[base]{print}},
#' \code{\link{summary.mcglm}}
#'
#' @method print mcglm
#' @export

print.mcglm <- function(x, ...) {

  object <- x
  n_resp <- length(object$beta_names)

  if (!exists("mc_updateBeta", mode = "function")) {
    stop("Function 'mc_updateBeta' not found.")
  }

  regression <- mc_updateBeta(
    list_initial = list(),
    betas = object$Regression,
    information = object$Information,
    n_resp = n_resp
  )

  for (i in seq_len(n_resp)) {

    cat("Call:\n")
    print(object$linear_pred[[i]])
    cat("\n")

    cat("Link function: ", object$link[[i]], "\n", sep = "")
    cat("Variance function: ", object$variance[[i]], "\n", sep = "")
    cat("Covariance function: ", object$covariance[[i]], "\n\n", sep = "")

    coef_reg <- regression$regression[[i]]
    names(coef_reg) <- object$beta_names[[i]]

    cat("Regression coefficients:\n")
    print(coef_reg)
    cat("\n")

    tau_temp <- coef(object, response = i, type = "tau")$Estimates
    if (length(tau_temp) > 0) {
      cat("Dispersion parameters:\n")
      print(unname(tau_temp))
      cat("\n")
    }

    power_temp <- coef(object, response = i, type = "power")$Estimates
    if (length(power_temp) > 0) {
      cat("Power parameters:\n")
      print(unname(power_temp))
      cat("\n")
    }
  }

  invisible(x)
}

#' @title Residuals for mcglm Objects
#' @name residuals.mcglm
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description
#' Computes residuals for a fitted \code{mcglm} object. Different types
#' of residuals can be extracted, depending on the specified argument
#' \code{type}.
#'
#' @param object
#' An object of class \code{mcglm}.
#'
#' @param type
#' A character string specifying the type of residuals to be returned.
#' Options are:
#' \describe{
#'   \item{\code{"raw"}}{Raw residuals, defined as observed minus fitted values.}
#'   \item{\code{"pearson"}}{Pearson residuals, scaled by the marginal standard deviation.}
#'   \item{\code{"standardized"}}{Standardized residuals, obtained using the inverse covariance matrix.}
#' }
#'
#' @param ...
#' Further arguments passed to or from other methods. Currently ignored.
#'
#' @return
#' A numeric matrix of class \code{Matrix} with dimensions
#' \eqn{n \times r}, where \eqn{n} is the number of observations and
#' \eqn{r} is the number of response variables.
#'
#' @seealso
#' \code{\link[stats]{residuals}},
#' \code{\link{fitted.mcglm}}
#'
#' @method residuals mcglm
#' @export

residuals.mcglm <- function(object,
                            type = c("raw", "pearson", "standardized"),
                            ...) {

  type <- match.arg(type)

  n_resp <- length(object$beta_names)
  n_obs  <- object$n_obs

  res_raw <- Matrix::Matrix(
    object$residuals,
    ncol = n_resp,
    nrow = n_obs
  )

  output <- switch(
    type,

    raw = res_raw,

    pearson = {
      sd_marginal <- sqrt(diag(object$C))
      if (any(sd_marginal <= 0)) {
        stop("Non-positive marginal variances detected.")
      }
      Matrix::Matrix(
        as.numeric(object$residuals / sd_marginal),
        ncol = n_resp,
        nrow = n_obs
      )
    },

    standardized = {
      chol_invC <- tryCatch(
        chol(object$inv_C),
        error = function(e)
          stop("Cholesky decomposition of inv_C failed.")
      )
      Matrix::Matrix(
        as.numeric(object$residuals %*% chol_invC),
        ncol = n_resp,
        nrow = n_obs
      )
    }
  )

  return(output)
}

#' @title Summary for mcglm Objects
#' @name summary.mcglm
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description
#' Produces a comprehensive summary of a fitted \code{mcglm} object.
#' The summary includes estimates, standard errors, Wald Z statistics,
#' and p-values for regression, dispersion, power, and correlation parameters.
#' The function can either print the summary to the console or return it invisibly.
#'
#' @param object
#' A fitted object of class \code{mcglm}.
#'
#' @param verbose
#' Logical; if \code{TRUE} (default), prints the summary to the console.
#' If \code{FALSE}, the summary is returned invisibly.
#'
#' @param print
#' Character vector specifying which components of the summary to print.
#' Possible values are \code{"Regression"}, \code{"Power"}, \code{"Dispersion"},
#' and \code{"Correlation"}. Default prints all components.
#'
#' @param ...
#' Further arguments passed to or from other methods. Currently ignored.
#'
#' @return
#' Invisibly returns a list containing summary tables:
#' \itemize{
#'   \item For each response: \code{Regression}, \code{Power}, \code{Dispersion}.
#'   \item If applicable: \code{Correlation} summary.
#' }
#' Each table contains parameter estimates, standard errors, Wald Z values,
#' and two-sided p-values.
#'
#' @details
#' The summary also prints information about the model fitting, including
#' the link, variance, and covariance functions used, the algorithm method,
#' any correction applied, and the number of iterations completed.
#'
#' @seealso
#' \code{\link{print.mcglm}}, \code{\link{coef.mcglm}}, \code{\link{residuals.mcglm}}
#'
#' @examples
#' library(mcglm)
#' set.seed(123)
#' mydata <- data.frame(y = rnorm(10), x1 = rnorm(10),
#'                     x2 = rbinom(10, size = 1, prob = 0.5))
#' Z0 <- mc_id(mydata)
#' fit <- mcglm(c(y ~ x1 + x2), matrix_pred = list(Z0), data = mydata)
#' # Print full summary
#' summary(fit)
#' # Get summary invisibly
#' out <- summary(fit, verbose = FALSE)
#'
#'
#' @method summary mcglm
#' @export

summary.mcglm <- function(object,
                          verbose = TRUE,
                          print = c("Regression", "power",
                                    "Dispersion", "Correlation"),
                          ...) {

  print <- match.arg(print, several.ok = TRUE)
  n_resp <- length(object$beta_names)
  output <- vector("list", n_resp)

  for (i in seq_len(n_resp)) {

    resp_out <- list()

    ## Regression
    tab_beta <- coef(object, std.error = TRUE,
                     response = i, type = "beta")[, c(1,2), drop = FALSE]
    tab_beta$`Z value` <- tab_beta[, 1] / tab_beta[, 2]
    tab_beta$`Pr(>|z|)` <- 2 * stats::pnorm(-abs(tab_beta$`Z value`))
    rownames(tab_beta) <- object$beta_names[[i]]
    resp_out$Regression <- tab_beta

    ## Power
    tab_power <- coef(object, std.error = TRUE,
                      response = i, type = "power")[, c(1,2), drop = FALSE]
    if (nrow(tab_power) > 0) {
      tab_power$`Z value` <- tab_power[, 1] / tab_power[, 2]
      tab_power$`Pr(>|z|)` <- 2 * stats::pnorm(-abs(tab_power$`Z value`))
      rownames(tab_power) <- NULL
      resp_out$Power <- tab_power
    }

    ## Dispersion (tau)
    tab_tau <- coef(object, std.error = TRUE,
                    response = i, type = "tau")[, c(1,2), drop = FALSE]
    tab_tau$`Z value` <- tab_tau[, 1] / tab_tau[, 2]
    tab_tau$`Pr(>|z|)` <- 2 * stats::pnorm(-abs(tab_tau$`Z value`))
    rownames(tab_tau) <- NULL
    resp_out$Dispersion <- tab_tau

    output[[i]] <- resp_out
  }

  ## Correlation
  tab_rho <- coef(object, std.error = TRUE,
                  response = NA, type = "correlation")
  if (nrow(tab_rho) > 0) {
    tab_rho <- tab_rho[, c("Parameters", "Estimates", "Std.error")]
    tab_rho$`Z value` <- tab_rho$Estimates / tab_rho$Std.error
    tab_rho$`Pr(>|z|)` <- 2 * stats::pnorm(-abs(tab_rho$`Z value`))
    rownames(tab_rho) <- NULL
    output$Correlation <- tab_rho
  }

  ## Printing
  if (isTRUE(verbose)) {

    for (i in seq_len(n_resp)) {

      if ("Regression" %in% print) {
        cat("Call: ")
        print(object$linear_pred[[i]])
        cat("\nLink function:", object$link[[i]], "\n")
        cat("Variance function:", object$variance[[i]], "\n")
        cat("Covariance function:", object$covariance[[i]], "\n\n")
        cat("Regression:\n")
        print(output[[i]]$Regression)
        cat("\n")
      }

      if ("power" %in% print && !is.null(output[[i]]$Power)) {
        cat("Power:\n")
        print(output[[i]]$Power)
        cat("\n")
      }

      if ("Dispersion" %in% print) {
        cat("Dispersion:\n")
        print(output[[i]]$Dispersion)
        cat("\n")
      }
    }

    if ("Correlation" %in% print && !is.null(output$Correlation)) {
      cat("Correlation:\n")
      print(output$Correlation)
      cat("\n")
    }

    cat("Algorithm:", object$con$method, "\n")
    cat("Correction:", object$con$correct, "\n")
    n_iter <- length(stats::na.exclude(object$IterationCovariance[, 1]))
    cat("Number iterations:", n_iter, "\n")
  }
  return(invisible(output))
}

#' @title Variance-Covariance Matrix for mcglm Objects
#' @name vcov.mcglm
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description
#' Extracts the variance-covariance matrix of the estimated parameters
#' from a fitted \code{mcglm} object.
#'
#' @param object
#' An object of class \code{mcglm}.
#'
#' @param ...
#' Further arguments passed to or from other methods. Currently ignored.
#'
#' @return
#' A numeric matrix representing the variance-covariance matrix of all
#' estimated model parameters. Row and column names correspond to the
#' parameter identifiers.
#'
#' @seealso
#' \code{\link{coef.mcglm}},
#' \code{\link{summary.mcglm}}
#'
#' @method vcov mcglm
#' @export

vcov.mcglm <- function(object, ...) {

  vc <- object$vcov

  param_names <- coef(object)$Parameters

  if (!is.null(param_names) &&
      length(param_names) == nrow(vc)) {
    colnames(vc) <- param_names
    rownames(vc) <- param_names
  }

  return(vc)
}
