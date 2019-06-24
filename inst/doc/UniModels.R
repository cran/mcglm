## ----setup, include=FALSE-----------------------------------------
##----------------------------------------------------------------------

library(knitr)

opts_chunk$set(
    dev.args=list(family="Palatino"))

options(width=68)

## ---- eval=FALSE--------------------------------------------------
#  library(devtools)
#  install_git("wbonat/mcglm")

## ---- eval=FALSE, error=FALSE, message=FALSE, warning=FALSE-------
#  library(mcglm)
#  packageVersion("mcglm")

## ---- warning = FALSE, message = FALSE----------------------------
# Loading extra packages
require(mcglm)
require(Matrix)
require(mvtnorm)
require(tweedie)

# Setting the seed
set.seed(2503)

# Fixed component
x1 <- seq(-1,1, l = 100)
X <- model.matrix(~ x1)
mu1 <- mcglm::mc_link_function(beta = c(1,0.8), X = X, offset = NULL, 
                        link = "identity")
# Random component
y1 <- rnorm(100, mu1$mu, sd = 0.5)

# Data structure
data <- data.frame("y1" = y1, "x1" = x1)

# Matrix linear predictor
Z0 <- mc_id(data)

# Fit
fit1.id <- mcglm(linear_pred = c(y1 ~ x1), 
                 matrix_pred = list(Z0),
                 data = data)


## -----------------------------------------------------------------
print(methods(class = "mcglm"))

## -----------------------------------------------------------------
summary(fit1.id)

## -----------------------------------------------------------------
# Fit using inverse covariance link function 
fit1.inv <- mcglm(linear_pred = c(y1 ~ x1), 
                  matrix_pred = list(Z0),
                  covariance = "inverse", data = data)


## ---- message=FALSE, warning=FALSE--------------------------------
# Fit using expm covariance link function
fit1.expm <- mcglm(linear_pred = c(y1 ~ x1), 
                   matrix_pred = list(Z0),
                   covariance = "expm", data = data)

## -----------------------------------------------------------------
# Comparing estimates using different covariance link functions
cbind(coef(fit1.id)$Estimates,
      coef(fit1.inv)$Estimates,
      coef(fit1.expm)$Estimates)

# Applying the inverse transformation
c(coef(fit1.id)$Estimates[3],
  1/coef(fit1.inv)$Estimates[3],
  exp(coef(fit1.expm)$Estimates[3]))

## -----------------------------------------------------------------
# Mean model
set.seed(1811)
x1 <- seq(-1,1, l = 100)
X <- model.matrix(~ x1)
mu1 <- mcglm::mc_link_function(beta = c(1,0.8), X = X, offset = NULL, 
                        link = "identity")
z1 <- rnorm(100, mean = 0, sd = 0.25)
data <- data.frame("id" = 1, "x1" = x1, "z1" = z1)

# Matrix linear predictor
Z <- mc_dglm(~ z1, id = 'id', data = data)

# Covariance model
Sigma <- mcglm::mc_matrix_linear_predictor(tau = c(0.2, 0.15), Z = Z)

# Simulating the response variable
y1 <- rnorm(100, mu1$mu, sd = sqrt(diag(Sigma)))
data$y <- y1

# Fitting
fit2.id <- mcglm(linear_pred = c(y1 ~ x1), matrix_pred = list(Z), data = data)

## -----------------------------------------------------------------
# Mean model
x1 <- seq(-1,1, l = 100)
X <- model.matrix(~ x1)
mu1 <- mcglm::mc_link_function(beta = c(1,0.8), X = X, offset = NULL, 
                        link = "identity")

# Data structure
data <- data.frame("id" = as.factor(rep(1:10, each = 10)), "x1" = x1)

# Covariance model
Z0 <- mc_id(data)
Z1 <- mc_mixed(~ 0 + id, data = data)
Sigma <- mcglm::mc_matrix_linear_predictor(tau = c(0.2, 0.15), Z = c(Z0,Z1))

# Simulating the Response variable
y1 <- as.numeric(rmvnorm(1, mean = mu1$mu, sigma = as.matrix(Sigma)))
data <- data.frame("y1" = y1, "x1" = x1)

# Fit
fit3.id <- mcglm(linear_pred = c(y1 ~ x1), 
                 matrix_pred = list("resp1" = c(Z0,Z1)), data = data)

## -----------------------------------------------------------------
summary(fit3.id)

## -----------------------------------------------------------------
# Mean model
x1 <- seq(-1,1, l = 500)
X <- model.matrix(~ x1)
mu1 <- mcglm::mc_link_function(beta = c(1,0.8), X = X, offset = NULL, link = "logit")

# Data structure
data <- data.frame("x1" = x1)

# Covariance model
Z0 <- mc_id(data)

# Simulating the response variable
set.seed(123)
data$y <- rbinom(500, prob = mu1$mu, size = 10)/10

## -----------------------------------------------------------------
# Fit
fit4.logit <- mcglm(linear_pred = c(y ~ x1), 
                    matrix_pred = list(Z0),
                    link = "logit", variance = "binomialP",
                    power_fixed = TRUE,
                    Ntrial = list(rep(10,500)), data = data)

## -----------------------------------------------------------------
fit4.cauchit <- mcglm(linear_pred = c(y ~ x1), 
                      matrix_pred = list(Z0),
                      link = "cauchit", variance = "binomialP", 
                      Ntrial = list(rep(10,250)), data = data)

## -----------------------------------------------------------------
fit4.logitP <- mcglm(linear_pred = c(y ~ x1), 
                      matrix_pred = list(Z0),
                      link = "logit", variance = "binomialP",
                      power_fixed = FALSE,
                      Ntrial = list(rep(10,500)), data = data)

## -----------------------------------------------------------------
fit4.logitPQ <- mcglm(linear_pred = c(y ~ x1), 
                      matrix_pred = list(Z0),
                      link = "logit", variance = "binomialPQ",
                      power_fixed = FALSE,
                      Ntrial = list(rep(10,500)), 
                      control_algorithm = list(tuning = 0.5, max_iter = 100),
                      data = data)

## -----------------------------------------------------------------
# Mean model

x1 <- seq(-2,2, l = 200)
X <- model.matrix(~ x1)
mu <- mcglm::mc_link_function(beta = c(1,0.8), X = X, offset = NULL, link = "log")

# Data structure
data <- data.frame("x1" = x1)

# Covariance model
Z0 <- mc_id(data)

# Data structure
data$y <- rpois(200, lambda = mu$mu)

# Fit
fit.poisson <- mcglm(linear_pred = c(y ~ x1), 
                    matrix_pred = list(Z0),
                    link = "log", variance = "tweedie",
                    power_fixed = TRUE, data = data)

## -----------------------------------------------------------------
# Simulating negative binomial models
set.seed(1811)
x <- rtweedie(200, mu = mu$mu, power = 2, phi = 0.5)
y <- rpois(200, lambda = x)
data <- data.frame("y1" = y, "x1" = x1)

fit.pt <- mcglm(linear_pred = c(y ~ x1), matrix_pred = list(Z0), 
                link = "log", variance = "poisson_tweedie", 
                power_fixed = FALSE, data = data)
summary(fit.pt)

