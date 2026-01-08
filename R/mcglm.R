#' @title Australian Health Survey Data
#' @name ahs
#'
#' @description
#' The Australian Health Survey (AHS) was used by Bonat and Jorgensen (2016)
#' as an example of multivariate count regression modeling. The dataset
#' contains five count response variables related to health system usage and
#' nine covariates related to social conditions in Australia for the years
#' 1987-88.
#'
#' @format A \code{data.frame} with 5190 observations and 15 variables:
#' \describe{
#'   \item{\code{sex}}{Factor with levels \code{male} and \code{female}.}
#'   \item{\code{age}}{Respondent's age in years divided by 100.}
#'   \item{\code{income}}{Respondent's annual income in Australian dollars divided by 1000.}
#'   \item{\code{levyplus}}{Factor indicating coverage by private health insurance for private patients in public hospital with doctor of choice (1) or otherwise (0).}
#'   \item{\code{freepoor}}{Factor indicating government coverage due to low income, recent immigration, or unemployment (1) or otherwise (0).}
#'   \item{\code{freerepa}}{Factor indicating government coverage due to old-age/disability pension, veteran status, or family of deceased veteran (1) or otherwise (0).}
#'   \item{\code{illnes}}{Number of illnesses in the past two weeks, capped at 5.}
#'   \item{\code{actdays}}{Number of days of reduced activity in the past two weeks due to illness or injury.}
#'   \item{\code{hscore}}{General health questionnaire score (Goldberg's method); higher scores indicate poorer health.}
#'   \item{\code{chcond}}{Factor with levels: \code{limited} (chronic condition with activity limitation), \code{nonlimited} (chronic condition without limitation), \code{otherwise} (reference level).}
#'   \item{\code{Ndoc}}{Number of consultations with a doctor or specialist (response variable).}
#'   \item{\code{Nndoc}}{Number of consultations with health professionals (response variable).}
#'   \item{\code{Nadm}}{Number of admissions to hospital, psychiatric hospital, nursing, or convalescence home in the past 12 months (response variable).}
#'   \item{\code{Nhosp}}{Number of nights in a hospital during the most recent admission.}
#'   \item{\code{Nmed}}{Total number of prescribed and non-prescribed medications used in the past two days.}
#' }
#'
#' @usage data(ahs)
#' @docType data
#' @keywords datasets
#'
#' @source Deb, P. and Trivedi, P. K. (1997) "Demand for medical care by the elderly: A finite mixture approach." \emph{Journal of Applied Econometrics}, 12(3):313--336.
#'
#' @source Bonat, W. H. and Jorgensen, B. (2016) "Multivariate covariance generalized linear models." \emph{Journal of the Royal Statistical Society: Series C}, 65:649--675.
#'
#' @examples
#' library(mcglm)
#' data(ahs, package = "mcglm")
#' form1 <- Ndoc ~ income + age
#' form2 <- Nndoc ~ income + age
#' Z0 <- mc_id(ahs)
#' fit.ahs <- mcglm(linear_pred = c(form1, form2),
#'                  matrix_pred = list(Z0, Z0),
#'                  link = c("log", "log"),
#'                  variance = c("poisson_tweedie", "poisson_tweedie"),
#'                  data = ahs)
#' summary(fit.ahs)
"ahs"

#' @title Hunting Data from Pico Basile, Bioko Island, Equatorial Guinea
#' @name Hunting
#'
#' @description
#' This dataset contains a case study analyzed in Bonat et al. (2017) regarding
#' animals hunted in the village of Basile Fang, Bioko Norte Province,
#' Bioko Island, Equatorial Guinea. Monthly counts of blue duikers and other
#' small animals shot or snared were collected from a random sample of 52
#' commercial hunters between August 2010 and September 2013. For each animal
#' caught, the species, sex, capture method, and altitude were recorded. The
#' dataset contains 1216 observations.
#'
#' @format A \code{data.frame} with 1216 observations and 11 variables:
#' \describe{
#'   \item{\code{ALT}}{Factor with five levels indicating the altitude where the animal was caught.}
#'   \item{\code{SEX}}{Factor with levels \code{Female} and \code{Male}.}
#'   \item{\code{METHOD}}{Factor with levels \code{Escopeta} and \code{Trampa} indicating the method of capture.}
#'   \item{\code{OT}}{Monthly number of other small animals hunted.}
#'   \item{\code{BD}}{Monthly number of blue duikers hunted.}
#'   \item{\code{OFFSET}}{Monthly number of hunter days.}
#'   \item{\code{HUNTER}}{Hunter index.}
#'   \item{\code{MONTH}}{Month index.}
#'   \item{\code{MONTHCALENDAR}}{Month as calendar number (1 = January, ..., 12 = December).}
#'   \item{\code{YEAR}}{Calendar year (2010–2013).}
#'   \item{\code{HUNTER.MONTH}}{Index indicating observations taken for the same hunter and month.}
#' }
#'
#' @usage data(Hunting)
#' @docType data
#' @keywords datasets
#'
#' @source Bonat, W. H., et al. (2017). "Modelling the covariance structure in marginal multivariate count models: Hunting in Bioko Island." \emph{Journal of Agricultural, Biological and Environmental Statistics}, 22(4):446–464.
#'
#' @source Bonat, W. H. (2018). "Multiple Response Variables Regression Models in R: The mcglm Package." \emph{Journal of Statistical Software}, 84(4):1–30.
#'
#' @examples
#' library(mcglm)
#' library(Matrix)
#' data(Hunting, package = "mcglm")
#' formu <- OT ~ METHOD*ALT + SEX + ALT*poly(MONTH, 4)
#' Z0 <- mc_id(Hunting)
#' Z1 <- mc_mixed(~0 + HUNTER.MONTH, data = Hunting)
#' fit <- mcglm(linear_pred = c(formu),
#'              matrix_pred = list(c(Z0, Z1)),
#'              link = c("log"),
#'              variance = c("poisson_tweedie"),
#'              power_fixed = c(FALSE),
#'              control_algorithm = list(max_iter = 100),
#'              offset = list(log(Hunting$OFFSET)),
#'              data = Hunting)
#' summary(fit)
#' anova(fit)
"Hunting"

#' @title Soil Chemistry Properties Dataset
#' @name soil
#'
#' @description
#' Soil chemistry properties measured on a regular grid of 10 × 25 points,
#' spaced by 5 meters. Each record contains the chemical composition and
#' coordinates of the soil sample.
#'
#' @format A \code{data.frame} with 250 observations and 9 variables:
#' \describe{
#'   \item{\code{COORD.X}}{X coordinate of the sampling point.}
#'   \item{\code{COORD.Y}}{Y coordinate of the sampling point.}
#'   \item{\code{SAND}}{Proportion of sand in the soil sample.}
#'   \item{\code{SILT}}{Proportion of silt in the soil sample.}
#'   \item{\code{CLAY}}{Proportion of clay in the soil sample.}
#'   \item{\code{PHWATER}}{Soil pH measured in water.}
#'   \item{\code{CA}}{Calcium content of the soil.}
#'   \item{\code{MG}}{Magnesium content of the soil.}
#'   \item{\code{K}}{Potassium content of the soil.}
#' }
#'
#' @usage data(soil)
#' @docType data
#'
#' @keywords datasets
#'
#' @source Bonat, W. H. (2018). "Multiple Response Variables Regression Models in R: The mcglm Package." \emph{Journal of Statistical Software}, 84(4):1–30.
#'
#' @examples
#' library(mcglm)
#' data(soil, package = "mcglm")
#'
#' # Spatial model (tri2nb could be used but takes long)
#' Z1 <- mc_id(soil)
#'
#' # Linear predictor example
#' form.ca <- CA ~ COORD.X*COORD.Y + SAND + SILT + CLAY + PHWATER
#' fit.ca <- mcglm(
#'   linear_pred = c(form.ca),
#'   matrix_pred = list(Z1),
#'   link = "log",
#'   variance = "tweedie",
#'   covariance = "inverse",
#'   power_fixed = TRUE,
#'   data = soil,
#'   control_algorithm = list(
#'     max_iter = 1000,
#'     tuning = 0.1,
#'     verbose = FALSE,
#'     tol = 1e-03
#'   )
#' )
"soil"

#' @title Respiratory Physiotherapy on Premature Newborns
#' @name NewBorn
#'
#' @description
#' The NewBorn dataset is from a prospective study assessing the effect of
#' respiratory physiotherapy on cardiopulmonary function in ventilated
#' preterm newborn infants with birth weight less than 1500 g. The dataset
#' was collected by the nursing team of Waldemar Monastier Hospital, Campo
#' Largo, PR, Brazil, and analyzed in Bonat and Jorgensen (2016) as an
#' example of mixed outcomes regression models.
#'
#' @format A \code{data.frame} with 270 observations and 21 variables:
#' \describe{
#'   \item{\code{Sex}}{Factor with levels \code{Female} and \code{Male}.}
#'   \item{\code{GA}}{Gestational age in weeks.}
#'   \item{\code{BW}}{Birth weight in grams.}
#'   \item{\code{APGAR1M}}{APGAR index at the first minute of life.}
#'   \item{\code{APGAR5M}}{APGAR index at the fifth minute of life.}
#'   \item{\code{PRE}}{Factor indicating prematurity (YES/NO).}
#'   \item{\code{HD}}{Factor indicating Hansen's disease (YES/NO).}
#'   \item{\code{SUR}}{Factor indicating surfactant administration (YES/NO).}
#'   \item{\code{JAU}}{Factor indicating jaundice (YES/NO).}
#'   \item{\code{PNE}}{Factor indicating pneumonia (YES/NO).}
#'   \item{\code{PDA}}{Factor indicating persistence of ductus arteriosus (YES/NO).}
#'   \item{\code{PPI}}{Factor indicating primary pulmonary infection (YES/NO).}
#'   \item{\code{OTHERS}}{Factor indicating other diseases (YES/NO).}
#'   \item{\code{DAYS}}{Age in days.}
#'   \item{\code{AUX}}{Factor indicating type of respiratory auxiliary (HOOD/OTHERS).}
#'   \item{\code{RR}}{Respiratory rate (continuous).}
#'   \item{\code{HR}}{Heart rate (continuous).}
#'   \item{\code{SPO2}}{Oxygen saturation (bounded).}
#'   \item{\code{TREAT}}{Factor with three levels: Respiratory physiotherapy, Evaluation 1, Evaluation 2, Evaluation 3.}
#'   \item{\code{NBI}}{Newborn index.}
#'   \item{\code{TIME}}{Days of treatment.}
#' }
#'
#' @usage data(NewBorn)
#' @docType data
#'
#' @keywords datasets
#'
#' @source Bonat, W. H. and Jorgensen, B. (2016). "Multivariate covariance generalized linear models." \emph{Journal of Royal Statistical Society, Series C}, 65:649–675.
#'
#' @examples
#' library(mcglm)
#' library(Matrix)
#' data(NewBorn, package = "mcglm")
#'
#' # Linear predictor example
#' formu <- SPO2 ~ Sex + APGAR1M + APGAR5M + PRE + HD + SUR
#' Z0 <- mc_id(NewBorn)
#' fit <- mcglm(
#'   linear_pred = c(formu),
#'   matrix_pred = list(Z0),
#'   link = "logit",
#'   variance = "binomialP",
#'   power_fixed = TRUE,
#'   data = NewBorn,
#'   control_algorithm = list(verbose = FALSE, tuning = 0.5)
#' )
#' summary(fit)
"NewBorn"


#' @title Soybeans Experiment Data
#' @name soya
#'
#' @description
#' Dataset from an experiment conducted in a vegetation house with soybeans.
#' Each plot contained two plants and the experiment involved three levels of
#' soil water (\code{water}) and five levels of potassium fertilization (\code{pot}),
#' arranged in five blocks (\code{block}). Three response variables are recorded:
#' grain yield, number of seeds, and number of viable peas per plant. The dataset
#' contains 75 observations and 7 variables.
#'
#' @format A \code{data.frame} with 75 observations and 7 variables:
#' \describe{
#'   \item{\code{pot}}{Factor with five levels of potassium fertilization.}
#'   \item{\code{water}}{Factor with three levels of amount of water in the soil.}
#'   \item{\code{block}}{Factor with five levels representing experimental blocks.}
#'   \item{\code{grain}}{Continuous variable representing grain yield per plant.}
#'   \item{\code{seeds}}{Count variable representing number of seeds per plant.}
#'   \item{\code{viablepeas}}{Binomial variable representing number of viable peas per plant.}
#'   \item{\code{totalpeas}}{Binomial variable representing total number of peas per plant.}
#' }
#'
#' @usage data(soya)
#' @docType data
#'
#' @keywords datasets
#'
#' @source Bonat, W. H. (2018). Multiple Response Variables Regression Models in R: The mcglm Package. \emph{Journal of Statistical Software}, 84(4):1--30.
#'
#' @examples
#' library(mcglm)
#' library(Matrix)
#' data(soya, package = "mcglm")
#'
#' # Linear predictor example
#' formu <- grain ~ block + factor(water) * factor(pot)
#' Z0 <- mc_id(soya)
#' fit <- mcglm(linear_pred = c(formu), matrix_pred = list(Z0),
#'             data = soya)
#' anova(fit)
"soya"
