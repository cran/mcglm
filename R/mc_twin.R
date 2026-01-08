#' @title Twin Model Covariance Structures
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description
#' Constructs the components of the matrix linear predictor for twin data
#' analysis under ACDE-type models. The function generates covariance
#' structures suitable for monozygotic (MZ) and dizygotic (DZ) twins and
#' supports several biologically motivated and flexible model
#' parameterizations.
#'
#' @param id
#' A string indicating the name of the column in \code{data} that
#' identifies the twin pair. The same identifier must be shared by both
#' twins in a pair.
#'
#' @param twin.id
#' A string indicating the name of the column in \code{data} that
#' identifies the twin within each pair. Typically coded as \code{1} and
#' \code{2}.
#'
#' @param type
#' A string indicating the name of the column in \code{data} that
#' identifies the zygosity type. This variable must be a factor with
#' exactly two levels: \code{"mz"} and \code{"dz"}, where \code{"mz"} is
#' taken as the reference level.
#'
#' @param replicate
#' An optional string indicating the name of the column in \code{data}
#' that identifies replicated observations within the same twin pair,
#' such as time points in longitudinal twin studies. If provided, it is
#' treated as a factor.
#'
#' @param formula
#' Internal argument used to define flexible and unstructured covariance
#' models. Not intended for direct user specification.
#'
#' @param structure
#' A string specifying the covariance structure to be constructed.
#' Available options are \code{"full"}, \code{"flex"}, \code{"uns"},
#' \code{"ACE"}, \code{"ADE"}, \code{"AE"}, \code{"CE"} and \code{"E"}.
#'
#' @param data
#' A data frame containing all variables referenced by the model.
#'
#' @details
#' For biologically motivated structures (\code{"ACE"}, \code{"ADE"},
#' \code{"AE"}, \code{"CE"}, \code{"E"}), the function builds covariance
#' matrices based on classical twin modeling assumptions. For flexible
#' and unstructured options (\code{"full"}, \code{"flex"}, \code{"uns"}),
#' the covariance structure is constructed using matrix linear predictors.
#'
#' @return
#' A list of sparse matrices of class \code{dgCMatrix}, representing the
#' components of the matrix linear predictor to be used in the
#' \code{matrix_pred} argument of \code{\link{mcglm}}.
#'
#' @source
#' Bonat, W. H. (2018). Multiple Response Variables Regression Models in R:
#' The mcglm Package. Journal of Statistical Software, 84(4), 1--30.
#'
#' @seealso
#' \code{\link{mc_id}}, \code{\link{mc_dist}}, \code{\link{mc_car}},
#' \code{\link{mc_rw}}, \code{\link{mc_ns}}, \code{\link{mc_dglm}},
#' \code{\link{mc_mixed}}.
#'
#' @importFrom stats relevel
#' @export


mc_twin <- function(id, twin.id, type, replicate = NULL, structure, data) {
  # Checking the structure
  n_levels <- length(which(levels(data[[type]]) %in% c("mz","dz")))
  if(n_levels != 2) {
    stop("Levels of type column should be mz and dz.")
  }
  classe_twin_id <- class(data[[twin.id]])
  if(classe_twin_id != "factor") {
    data[[twin.id]] <- as.factor(data[[twin.id]])
    warning("Converted twin.id class to factor.")
  }
  if(!is.null(replicate)) {
    classe_replicate <- class(data[[twin.id]])
    if(classe_replicate != "factor") {
      data[[replicate]] <- as.factor(data[[replicate]])
      warning("Converted replicate class to factor.")
    }
  }
  # mz is set up as reference level
  # Models structures
  if(structure == "ACE") {
    E <- mc_id(data)
    A <- mc_twin_bio(id = id, twin.id = twin.id, replicate = replicate,
                     type = type, structure = "A", data = data)
    C <- mc_twin_bio(id = id, twin.id = twin.id, replicate = replicate,
                     type = type, structure = "C", data = data)
    output <- c(E,A,C)
  }
  if(structure == "ADE") {
    E <- mc_id(data)
    A <- mc_twin_bio(id = id, twin.id = twin.id, replicate = replicate,
                     type = type, structure = "A", data = data)
    D <- mc_twin_bio(id = id, twin.id = twin.id, replicate = replicate,
                     type = type, structure = "D", data = data)
    output <- c(E,A,D)
  }
  if(structure == "AE") {
    E <- mc_id(data)
    A <- mc_twin_bio(id = id, twin.id = twin.id, replicate = replicate,
                     type = type, structure = "A", data = data)
    output <- c(E,A)
  }
  if(structure == "CE") {
    E <- mc_id(data)
    C <- mc_twin_bio(id = id, twin.id = twin.id, replicate = replicate,
                     type = type, structure = "C", data = data)
    output <- c(E,C)
  }
  if(structure == "E") {
    E <- mc_id(data)
    output <- E
  }
  if(structure == "full") {
    data[[type]] <- relevel(data[[type]], ref = "mz")
    warning("Type reference is mz.")
    if(is.null(replicate)) {
      formula <- as.formula(paste("~", paste(twin.id, type, sep = "*"), sep = ""))
    }
    if(!is.null(replicate)) {
      formula <- as.formula(paste("~", paste(paste(twin.id, type, sep = "*"),
                                             replicate, sep="*"), sep = ""))
    }
    output <- mc_twin_full(id = id, twin.id = twin.id, type = type,
                 replicate = replicate, formula = formula, data = data)
  }
  if(structure == "flex") {
    data[[type]] <- relevel(data[[type]], ref = "mz")
    warning("Type reference is mz.")
    if(is.null(replicate)) {
      formula <- as.formula(paste("~",type))
    }
    if(!is.null(replicate)) {
      formula <- as.formula(paste("~", paste(paste(type, sep = ""),
                                             replicate, sep="*"), sep = ""))
    }
    formula <- as.formula(paste("~", type))
    output <- mc_twin_full(id = id, twin.id = twin.id, type = type,
                           replicate = replicate, formula = formula,
                           data = data)
  }
  if(structure == "uns") {
    data[[type]] <- relevel(data[[type]], ref = "mz")
    warning("Type reference is mz.")
    output <- mc_twin_full(id = id, twin.id = twin.id, type = type,
                           replicate = replicate, formula = ~ 1, data = data)
  }
  return(output)
}

#' @rdname mc_twin
mc_twin_bio <- function(id, twin.id, type, replicate = NULL,
                        structure, data) {
  # A matrix
  if ( structure == "A" ) {
  MZ <- Matrix(c(1, 1, 1, 1), 2, 2)
  DZ <- Matrix(c(1, 0.5, 0.5, 1), 2, 2)
  }
  # C matrix
  if( structure == "C") {
  MZ <- Matrix(c(1, 1, 1, 1), 2, 2)
  DZ <- Matrix(c(1, 1, 1, 1), 2, 2)
  }
  # D matrix
  if( structure == "D") {
  MZ <- Matrix(c(1, 1, 1, 1), 2, 2)
  DZ <- Matrix(c(1, 0.25, 0.25, 1), 2, 2)
  }
  data.twin <- split(data, data[id])
  DD <- sum(abs(diff(do.call(c,lapply(data.twin, function(x)dim(x)[1])))))
  if( DD != 0) {
    stop("Model requires equal number of observations per twin pair. \n")
  }
  M_list <- list()
  for(i in 1:length(data.twin)) {
    twin_temp <- data.twin[[i]]
    twin_temp$COD_ORDER1 <- 1:dim(twin_temp)[1]
    twin_temp <- twin_temp[order(twin_temp[twin.id]),]
    twin_temp$COD_ORDER2 <- 1:dim(twin_temp)[1]
    if(!is.null(replicate)) {
      n_replicate <- length(unique(twin_temp[[replicate]]))
      R_matrix <- Diagonal(n_replicate, 1)
    }
    if(is.null(replicate)) {
      n_replicate <- 1
      R_matrix <- Diagonal(n_replicate, 1)
    }
    if(unique(twin_temp[type]) == "mz") {
      twin_temp = twin_temp[order(twin_temp$COD_ORDER1),]
      M_temp <- kronecker(MZ, R_matrix)[twin_temp$COD_ORDER2,twin_temp$COD_ORDER2]
    }
    if(unique(twin_temp[type]) == "dz") {
      twin_temp = twin_temp[order(twin_temp$COD_ORDER1),]
      M_temp <- kronecker(DZ, R_matrix)[twin_temp$COD_ORDER2,twin_temp$COD_ORDER2]
    }
    M_list[[i]] <- M_temp
  }
  output <- forceSymmetric(bdiag(M_list))
  return(list(output))
}


#' @rdname mc_twin
mc_twin_full <- function(id, twin.id, type, replicate, formula, data) {
  VV <- mc_dglm(formula, id = id, data = data)
  MAT_MZ <- mc_ns(id = id, data = data)
  MAT_DZ <- mc_ns(id = id, data = data, group = type, marca = "mz")
  output <- c(VV, MAT_MZ, MAT_DZ)
  return(output)
}




