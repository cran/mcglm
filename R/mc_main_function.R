#' @title Fitting Multivariate Covariance Generalized Linear Models
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description Fits multivariate covariance generalized linear models.
#'   Models are specified through lists defining the linear predictors
#'   and matrix linear predictors. The user can choose among different
#'   link, variance, and covariance functions. Model fitting is based on
#'   an estimating function approach, combining quasi-score functions for
#'   regression parameters and Pearson estimating functions for covariance
#'   parameters.
#'
#' @param linear_pred A list of model formulas, one for each response.
#'   See \code{\link[stats]{formula}} for details.
#' @param matrix_pred A list of known matrices defining the matrix linear
#'   predictor for the covariance structure. See
#'   \code{\link[mcglm]{mc_matrix_linear_predictor}} for details.
#'
#' @param link A list of link function names, one for each response.
#'   Possible values are \code{"logit"}, \code{"probit"}, \code{"cauchit"},
#'   \code{"cloglog"}, \code{"loglog"}, \code{"identity"}, \code{"log"},
#'   \code{"sqrt"}, \code{"1/mu^2"}, and \code{"inverse"}.
#' @param variance A list of variance function names.
#'   Possible values are \code{"constant"}, \code{"tweedie"},
#'   \code{"poisson_tweedie"}, \code{"binomialP"}, and \code{"binomialPQ"}.
#' @param covariance A list of covariance link function names.
#'   Possible values are \code{"identity"}, \code{"inverse"}, and
#'   \code{"expm"}.
#' @param offset A list of numeric vectors specifying offsets for each
#'   response. Use \code{NULL} if no offset is required.
#' @param Ntrial A list of numeric vectors specifying the number of trials
#'   for binomial responses. Only used for \code{binomialP} and
#'   \code{binomialPQ} variance functions.
#' @param power_fixed A list of logical values indicating whether the power
#'   parameter should be fixed (\code{TRUE}) or estimated (\code{FALSE}).
#' @param weights A list of numeric vectors of observation weights.
#'   Each element must have length equal to the number of observations.
#'   Missing observations should be coded as \code{NA}.
#' @param control_initial A list of initial values for the fitting algorithm.
#'   If set to \code{"automatic"}, initial values are generated internally
#'   using \code{\link{mc_initial_values}}.
#' @param control_algorithm A list of control parameters passed to the
#'   fitting algorithm. See \code{\link[mcglm]{fit_mcglm}} for details.
#' @param contrasts An optional list of contrasts passed to
#'   \code{\link[stats]{model.matrix}}.
#' @param data a data frame.
#' @usage mcglm(linear_pred, matrix_pred, link, variance, covariance,
#'        offset, Ntrial, power_fixed, data, control_initial,
#'        contrasts, weights, control_algorithm)
#' @return
#' An object of class \code{"mcglm"} representing a fitted multivariate
#' covariance generalized linear model.
#'
#' The returned object is a list produced by the fitting routine
#' \code{\link{fit_mcglm}}, augmented with additional components used by
#' post–estimation methods. The main components include:
#'
#' \describe{
#'   \item{beta_names}{A list of character vectors giving the names of the
#'     regression coefficients for each response variable.}
#'   \item{power_fixed}{A list of logical values indicating whether the
#'     power parameters were fixed or estimated.}
#'   \item{list_initial}{A list of initial values used in the fitting
#'     algorithm.}
#'   \item{n_obs}{An integer giving the number of observations.}
#'   \item{link}{A list of link functions used in the model.}
#'   \item{variance}{A list of variance functions used in the model.}
#'   \item{covariance}{A list of covariance link functions used in the model.}
#'   \item{linear_pred}{A list of formulas defining the linear predictors.}
#'   \item{matrix_pred}{A list of matrices defining the matrix linear predictors.}
#'   \item{list_X}{A list of design matrices corresponding to the linear
#'     predictors.}
#'   \item{observed}{A matrix of observed response values, with rows
#'     corresponding to observations and columns to response variables.}
#'   \item{Ntrial}{A list containing the number of trials for each response
#'     variable, when applicable.}
#'   \item{offset}{A list of offset vectors used in the model.}
#'   \item{sparse}{A list of logical values indicating whether each matrix
#'     linear predictor is treated as sparse.}
#'   \item{weights}{A numeric vector of weights used in the fitting process.}
#'   \item{data}{The data frame used to fit the model.}
#'   \item{con}{A list of control parameters used by the fitting algorithm.}
#' }
#'
#' Additional components may be present for internal use by methods such as
#' \code{print}, \code{summary} and \code{predict}.
#'
#' @seealso \code{fit_mcglm}, \code{mc_link_function} and
#' \code{mc_variance_function}.
#'
#' @source Bonat, W. H. and Jorgensen, B. (2016) Multivariate
#'     covariance generalized linear models.
#'     Journal of Royal Statistical Society - Series C 65:649--675.
#'
#' @source Bonat, W. H. (2018). Multiple Response Variables Regression
#' Models in R: The mcglm Package. Journal of Statistical Software, 84(4):1--30.
#'
#' @export
#' @import Matrix

mcglm <- function(linear_pred, matrix_pred, link, variance,
                  covariance, offset, Ntrial, power_fixed,
                  data, control_initial = "automatic",
                  contrasts = NULL, weights = NULL,
                  control_algorithm = list()) {
  n_resp <- length(linear_pred)
  linear_pred <- as.list(linear_pred)
  matrix_pred <- as.list(matrix_pred)
  if (missing(link)) {
    link <- rep("identity", n_resp)
  }
  if (missing(variance)) {
    variance <- rep("constant", n_resp)
  }
  if (missing(covariance)) {
    covariance <- rep("identity", n_resp)
  }
  if (missing(offset)) {
    offset <- rep(list(NULL), n_resp)
  }
  if (missing(Ntrial)) {
    Ntrial <- rep(list(rep(1, dim(data)[1])), n_resp)
  }
  if (missing(power_fixed)) {
    power_fixed <- rep(TRUE, n_resp)
  }
  if (missing(contrasts)) {
    contrasts <- NULL
  }
  link <- as.list(link)
  variance <- as.list(variance)
  covariance <- as.list(covariance)
  offset <- as.list(offset)
  Ntrial <- as.list(Ntrial)
  power_fixed <- as.list(power_fixed)
  if (!is.list(control_initial)) {
    control_initial <-
      mc_initial_values(linear_pred = linear_pred,
                        matrix_pred = matrix_pred, link = link,
                        variance = variance,
                        covariance = covariance, offset = offset,
                        Ntrial = Ntrial, contrasts = contrasts,
                        data = data)
    message("Automatic initial values selected.", "\n")
  }
  con <- list(correct = TRUE, max_iter = 20, tol = 1e-04,
              method = "chaser", tuning = 1, verbose = FALSE)
  con[(namc <- names(control_algorithm))] <- control_algorithm
  list_model_frame <- lapply(linear_pred, model.frame, na.action = 'na.pass', data = data)

  old_na_action <- getOption("na.action")
  on.exit(options(na.action = old_na_action), add = TRUE)

  if (!is.null(contrasts)) {
    list_X <- vector("list", n_resp)
    for (i in seq_len(n_resp)) {
      options(na.action = "na.pass")
      list_X[[i]] <- model.matrix(linear_pred[[i]], contrasts = contrasts[[i]])
    }
  } else {
    options(na.action = "na.pass")
    list_X <- lapply(linear_pred, model.matrix, data = data)
  }

  list_Y <- lapply(list_model_frame, model.response)
  y_vec <- as.numeric(do.call(c, list_Y))
  if(is.null(weights)) {
    C <- rep(1, length(y_vec))
    C[is.na(y_vec)] = 0
    weights = C
    y_vec[is.na(y_vec)] <- 0
  }
  if(!is.null(weights)) {
    y_vec[is.na(y_vec)] <- 0
    if(!inherits(weights, "list")) {weights <- as.list(weights)}
    weights <- as.numeric(do.call(c, weights))
  }
  sparse <- lapply(matrix_pred, function(x) {
    if (inherits(x, "dgeMatrix")) {
      FALSE
    } else TRUE
  })
  model_fit <- try(fit_mcglm(list_initial = control_initial,
                             list_link = link,
                             list_variance = variance,
                             list_covariance = covariance,
                             list_X = list_X, list_Z = matrix_pred,
                             list_offset = offset,
                             list_Ntrial = Ntrial,
                             list_power_fixed = power_fixed,
                             list_sparse = sparse, y_vec = y_vec,
                             correct = con$correct,
                             max_iter = con$max_iter, tol = con$tol,
                             method = con$method,
                             tuning = con$tuning,
                             verbose = con$verbose,
                             weights = weights))
  if (!inherits(model_fit, "try-error")) {
    model_fit$beta_names <- lapply(list_X, colnames)
    model_fit$power_fixed <- power_fixed
    model_fit$list_initial <- control_initial
    model_fit$n_obs <- dim(data)[1]
    model_fit$link <- link
    model_fit$variance <- variance
    model_fit$covariance <- covariance
    model_fit$linear_pred <- linear_pred
    model_fit$con <- con
    model_fit$observed <- Matrix(y_vec, ncol = length(list_Y),
                                 nrow = dim(data)[1])
    model_fit$list_X <- list_X
    model_fit$matrix_pred <- matrix_pred
    model_fit$Ntrial <- Ntrial
    model_fit$offset <- offset
    model_fit$power_fixed
    model_fit$sparse <- sparse
    model_fit$data <- data
    model_fit$weights <- weights
    class(model_fit) <- "mcglm"
  }
  n_it <- length(na.exclude(model_fit$IterationCovariance[,1]))
  if(con$max_it == n_it) {warning("Maximum iterations number reached. \n", call. = FALSE)}
  return(model_fit)
}
