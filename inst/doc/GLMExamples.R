## ----setup, include=FALSE-----------------------------------------
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
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))

## ---- warning = FALSE, message = FALSE----------------------------
fit.glm <- glm(counts ~ outcome + treatment, family = quasipoisson)

## ---- warning = FALSE, message = FALSE----------------------------
require(mcglm)
require(Matrix)
# Matrix linear predictor
Z0 <- mc_id(d.AD)
fit.qglm <- mcglm(linear_pred = c(counts ~ outcome + treatment),
                  matrix_pred = list("resp1" = Z0),
                  link = "log", variance = "tweedie", data = d.AD,
                  control_algorithm = list(verbose = FALSE,
                                           method = "chaser",
                                           tuning = 0.8))

## ---- warning = FALSE, message = FALSE----------------------------
cbind("GLM" = coef(fit.glm),
      "McGLM" = coef(fit.qglm, type = "beta")$Estimates)

cbind("GLM" = sqrt(diag(vcov(fit.glm))), 
      "McGLM" = coef(fit.qglm, type = "beta", std.error = TRUE)$Std.error)

## ---- warning = FALSE, message = FALSE----------------------------
# Loading the data set
utils::data(anorexia, package = "MASS")

# GLM fit
anorex.1 <- glm(Postwt ~ Prewt + Treat + offset(Prewt),
               family = gaussian, data = anorexia)

# McGLM fit
Z0 <- mc_id(anorexia)
fit.anorexia <- mcglm(linear_pred = c(Postwt ~ Prewt + Treat),
                      matrix_pred = list(Z0),
                      offset = list(anorexia$Prewt),
                      power_fixed = TRUE, data = anorexia,
                      control_algorithm = list("correct" = TRUE))

## ---- warning = FALSE, message = FALSE----------------------------
# Estimates
cbind("McGLM" = round(coef(fit.anorexia, type = "beta")$Estimates,5),
      "GLM" = round(coef(anorex.1),5))

# Standard errors
cbind("McGLM" = sqrt(diag(vcov(fit.anorexia))),
      "GLM" = c(sqrt(diag(vcov(anorex.1))),NA))

## ---- warning = FALSE, message = FALSE----------------------------
clotting <- data.frame(
  u = c(5,10,15,20,30,40,60,80,100),
  lot1 = c(118,58,42,35,27,25,21,19,18),
  lot2 = c(69,35,26,21,18,16,13,12,12))

## ---- warning = FALSE, message = FALSE----------------------------
fit.lot1 <- glm(lot1 ~ log(u), data = clotting,
                family = Gamma(link = "inverse"))
fit.lot2 <- glm(lot2 ~ log(u), data = clotting,
                family = Gamma(link = "inverse"))

## ---- warning = FALSE, message = FALSE----------------------------
list_initial = list()
list_initial$regression <- list(coef(fit.lot1))
list_initial$power <- list(c(2))
list_initial$tau <- list(summary(fit.lot1)$dispersion)
list_initial$rho = 0

## ---- warning = FALSE, message = FALSE----------------------------
Z0 <- mc_id(clotting)
fit.lot1.mcglm <- mcglm(linear_pred = c(lot1 ~ log(u)),
                        matrix_pred = list(Z0),
                        link = "inverse", variance = "tweedie",
                        data = clotting,
                        control_initial = list_initial)

## ---- warning = FALSE, message = FALSE----------------------------
# Estimates
cbind("mcglm" = round(coef(fit.lot1.mcglm, type = "beta")$Estimates,5),
      "glm" = round(coef(fit.lot1),5))
# Standard errors
cbind("mcglm" = sqrt(diag(vcov(fit.lot1.mcglm))),
      "glm" = c(sqrt(diag(vcov(fit.lot1))),NA))

## ---- warning = FALSE, message = FALSE----------------------------
list_initial$regression <- list("resp1" = coef(fit.lot2))
list_initial$tau <- list("resp1" = c(var(1/clotting$lot2)))

## ---- warning = FALSE, message = FALSE----------------------------
fit.lot2.mcglm <- mcglm(linear_pred = c(lot2 ~ log(u)),
                        matrix_pred = list(Z0),
                        link = "inverse", variance = "tweedie",
                        data = clotting,
                        control_initial = list_initial)

## ---- warning = FALSE, message = FALSE----------------------------
# Estimates
cbind("mcglm" = round(coef(fit.lot2.mcglm, type = "beta")$Estimates,5),
      "glm" = round(coef(fit.lot2),5))
# Standard errors
cbind("mcglm" = sqrt(diag(vcov(fit.lot2.mcglm))),
      "glm" = c(sqrt(diag(vcov(fit.lot2))),NA))

## ---- warning = FALSE, message = FALSE----------------------------
# Initial values
list_initial = list()
list_initial$regression <- list(coef(fit.lot1), coef(fit.lot2))
list_initial$power <- list(c(2),c(2))
list_initial$tau <- list(c(0.00149), c(0.001276))
list_initial$rho = 0.80

# Matrix linear predictor
Z0 <- mc_id(clotting)

# Fit bivariate Gamma model
fit.joint.mcglm <- mcglm(linear_pred = c(lot1 ~ log(u), lot2 ~ log(u)),
                         matrix_pred = list(Z0, Z0),
                         link = c("inverse", "inverse"),
                         variance = c("tweedie", "tweedie"),
                         data = clotting,
                         control_initial = list_initial,
                         control_algorithm = list("correct" = TRUE,
                                                 "method" = "chaser",
                                                 "tuning" = 0.1,
                                                 "max_iter" = 1000))
summary(fit.joint.mcglm)

## ---- warning = FALSE, message = FALSE----------------------------
# Initial values
list_initial = list()
list_initial$regression <- list(c(log(mean(clotting$lot1)),0),
                                c(log(mean(clotting$lot2)),0))
list_initial$power <- list(c(2), c(2))
list_initial$tau <- list(c(0.023), c(0.024))
list_initial$rho = 0

# Fit bivariate Gamma model
fit.joint.log <- mcglm(linear_pred = c(lot1 ~ log(u), lot2 ~ log(u)),
                       matrix_pred = list(Z0,Z0),
                       link = c("log", "log"),
                       variance = c("tweedie", "tweedie"),
                       data = clotting,
                       control_initial = list_initial)
summary(fit.joint.log)

## ---- warning = FALSE, message = FALSE----------------------------
require(MASS)
data(menarche)
data <- data.frame("resp" = menarche$Menarche/menarche$Total,
                   "Ntrial" = menarche$Total,
                   "Age" = menarche$Age)

## ---- warning = FALSE, message = FALSE----------------------------
glm.out = glm(cbind(Menarche, Total-Menarche) ~ Age,
              family=binomial(logit), data=menarche)

## ---- warning = FALSE, message = FALSE----------------------------
# Matrix linear predictor
Z0 <- mc_id(data)
fit.logit <- mcglm(linear_pred = c(resp ~ Age),
                   matrix_pred = list(Z0),
                   link = "logit", variance = "binomialP",
                   Ntrial = list(data$Ntrial), data = data)

## ---- warning = FALSE, message = FALSE----------------------------
# Estimates
cbind("GLM" = coef(glm.out),
      "McGLM" = coef(fit.logit, type = "beta")$Estimates)
# Standard error
cbind("GLM" = c(sqrt(diag(vcov(glm.out))),NA),
      "McGLM" =  sqrt(diag(vcov(fit.logit))))

## ---- warning = FALSE, message = FALSE----------------------------
fit.logit.power <- mcglm(linear_pred = c(resp ~ Age),
                         matrix_pred = list(Z0),
                         link = "logit", variance = "binomialP",
                         Ntrial = list(data$Ntrial),
                         power_fixed = FALSE, data = data)
summary(fit.logit.power)

