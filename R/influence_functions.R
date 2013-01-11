leverage <- function(model, ...){
  UseMethod("leverage", model)
}

leverage.default <- function(model, ...){
  stop(paste("there is no leverage() method for models of class",
             paste(class(model), collapse=", ")))
}

covratio <- function(model, ...){
  UseMethod("covratio", model)
}

covratio.default <- function(model, ...){
  stop(paste("there is no covratio() method for models of class",
             paste(class(model), collapse=", ")))
}

covratio.lm <- function(model, ...){
  function (model, infl = lm.influence(model, do.coef = FALSE), 
            res = weighted.residuals(model)) 
  {
    n <- nrow(qr.lm(model)$qr)
    p <- model$rank
    omh <- 1 - infl$hat
    e.star <- res/(infl$sigma * sqrt(omh))
    e.star[is.infinite(e.star)] <- NaN
    1/(omh * (((n - p - 1) + e.star^2)/(n - p))^p)
  }
}

covtrace <- function(model, ...){
  UseMethod("covtrace", model)
}

covtrace.default <- function(model, ...){
  stop(paste("there is no covtrace() method for models of class",
             paste(class(object), collapse=", ")))
}

mdffits <- function(model, ...){
  UseMethod("mdffits", model)
}

mdffits.default <- function(model, ...){
  stop(paste("there is no mdffits() method for models of class",
             paste(class(model), collapse=", ")))
}

rvc <- function(model, ...){
  UseMethod("rvc", model)
}

rvc.default <- function(model, ...){
  stop(paste("there is no rvc() method for models of class",
             paste(class(model), collapse=", ")))
}


#' Leverage for mixed/hierarchical linear models
#'
#' @param model fitted model of class \code{mer}
#' @param the level at which the leverage should be calculated; either
#'   1 or 2 (\code{"both"} can be specified)
#' @references 
#'   Nobre, J. S., & Singer, J. M. (2011). 
#'   Leverage analysis for linear mixed models. 
#'   Journal of Applied Statistics, 38(5), 1063–1072.
#'   
#'   Demidenko, E., & Stukel, T. A. (2005). 
#'   Influence analysis for linear mixed-effects models. 
#'   Statistics in Medicine, 24(6), 893–909.
leverage.mer <- function(model, level = "both") {
  if(!is(model, "mer")) stop("model must be of class 'mer'")
  if(model@dims[["nest"]] == 0) {
    stop("leverage.mer has not yet been implemented for models with crossed random effects")
  }
  if(!level %in% c("both", 1, 2)) {
    stop("level can only be 1, 2, or 'both'")
  }
  
  X <- getME(model, "X")
  # Z <- BlockZ(model)
  
  n     <- nrow(X)
  nt    <- model@dims[["nt"]]  # number of random-effects terms in the model
  ngrps <- unname( summary(model)@ngrps )
  
  vc   <- VarCorr(model)
  # D  <- kronecker( Diagonal(ngrps), bdiag(vc) )
  ZDZt <- attr(vc, "sc")^2 * crossprod( getME(model, "A") )
  R    <- Diagonal( n = n, x = attr(vc, "sc")^2 )
  
  V      <- ZDZt + R
  V.chol <- chol( V )
  Vinv   <- chol2inv( V.chol )
  
  xvix.inv <- attr(vc, "sc")^2 * chol2inv(getME(model, "RX"))
  
  H1 <- X %*% xvix.inv %*% t(X) %*% Vinv
  H2 <- ZDZt %*% (Diagonal( n = n ) - H1)
  
  diag.H1 <- diag(H1)
  diag.H2 <- diag(H2)
  flist <- getME(model, "flist")
  
  if(level != 2) {
    lev1 <- data.frame(fixef = diag.H1, ranef =  diag.H2)
  }
  if(level != 1) {
    lev2 <- data.frame( fixef = aggregate(diag.H1, flist, sum)[,2], 
                        ranef = aggregate(diag.H2, flist, sum)[,2])
  }
  
  if(level == 1) return(lev1)
  if(level == 2) return(lev2)
  if(level == "both") return(list(level.1 = lev1, level.2 = lev2))
}


#' Cook's distance for mixed/hierarchical linear models
#' 
#' @param model fitted model of class \code{mer}
#' @param group variable used to define the group for which cases will be deleted.
#'   If \code{group = NULL}, then individual cases will be deleted.
#' @param delete index of individual cases to be deleted. For higher level units
#'   specified in this manner, the \code{group} parameter must also be specified.
#'   If \code{case = NULL} then all cases are iteratively deleted.
cooks.distance.mer <- function(model, group = NULL, delete = NULL) {
  if(!is(model, "mer")) stop("model must be of class 'mer'")
  if(!is.null(group)) {
    if(!group %in% names(getME(model, "flist"))) {
      stop(paste(group, "is not a valid grouping factor for this model."))
    }
  }
  if(!model@dims["LMM"]){
    stop("cooks.distance is currently not implemented for GLMMs.")
  }
  
  # Extract key pieces of the model
  mats <- .mer_matrices(model)
  
  betaHat <- with(mats, XVXinv %*% t(X) %*% Vinv %*% Y)
  
  # Obtaining the building blocks
  if(is.null(group) & is.null(delete)) {
    calc.cooksd <- .Call("cooksdObs", y_ = mats$Y, X_ = as.matrix(mats$X), 
                         Vinv_ = as.matrix(mats$Vinv), 
                         XVXinv_ = as.matrix(mats$XVXinv), 
                         beta_ = as.matrix(betaHat), PACKAGE = "HLMdiag")
    res <- calc.cooksd[[1]]
    attr(res, "beta_cdd") <- calc.cooksd[[2]]
  }
  
  else{
    e <- with(mats, Y - X %*% betaHat)
    
    if( !is.null(group) ){
      grp.names <- unique( mats$flist[, group] )
      
      if( is.null(delete) ){
        del.index <- lapply(1:mats$ngrps[group], 
                            function(x) {
                              ind <- which(mats$flist[, group] == grp.names[x]) - 1
                            })
      } else{
        del.index <- list( which(mats$flist[, group] == delete) - 1 )
      }
    } else{
      del.index <- list( delete - 1 )
    }
    
    calc.cooksd <- .Call("cooksdSubset", index = del.index, X_ = X, P_ = P, 
                         Vinv_ = as.matrix(Vinv), XVXinv_ = as.matrix(XVXinv), 
                         e_ = as.numeric(e), PACKAGE = "HLMdiag")
    
    res <- calc.cooksd[[1]]
    attr(res, "beta_cdd") <- calc.cooksd[[2]] 
  }
  
# }
return(res)
}

#' MDFFITS for mixed/hierarchical linear models
#' 
#' @param model fitted model of class \code{mer}
#' @param group variable used to define the group for which cases will be deleted.
#'   If \code{group = NULL}, then individual cases will be deleted.
#' @param delete index of individual cases to be deleted. For higher level units
#'   specified in this manner, the \code{group} parameter must also be specified.
#'   If \code{case = NULL} then all cases are iteratively deleted.
mdffits.mer <- function(model, group = NULL, delete = NULL) {
  if(!is(model, "mer")) stop("model must be of class 'mer'")
  if(!is.null(group)) {
    if(!group %in% names(getME(model, "flist"))) {
      stop(paste(group, "is not a valid grouping factor for this model."))
    }
  }
  if(!model@dims["LMM"]){
    stop("mdffits is currently not implemented for GLMMs.")
  }
  
  # Extract key pieces of the model
  mats <- .mer_matrices(model)
  
  betaHat <- with(mats, XVXinv %*% t(X) %*% Vinv %*% Y)
  e <- with(mats, Y - X %*% betaHat)
  
  if( !is.null(group) ){
    grp.names <- unique( mats$flist[, group] )
    
    if( is.null(delete) ){
      del.index <- lapply(1:mats$ngrps[group], 
                          function(x) {
                            ind <- which(mats$flist[, group] == grp.names[x]) - 1
                          })
    } else{
      del.index <- list( which(mats$flist[, group] == delete) - 1 )
    }
  } else{
    if( is.null(delete) ){
      del.index <- split(0:(mats$n-1), 0:(mats$n-1))
    } else { del.index <- list( delete - 1 ) }
  }
  
  calc.mdffits <- .Call("mdffitsSubset", index = del.index, X_ = mats$X, 
                        P_ = mats$P, Vinv_ = as.matrix(mats$Vinv), 
                        XVXinv_ = as.matrix(mats$XVXinv), 
                        e_ = as.numeric(e), PACKAGE = "HLMdiag")
  res <- calc.mdffits[[1]]
  attr(res, "beta_cdd") <- calc.mdffits[[2]] 
  
  return(res)
}

#' COVRATIO for mixed/hierarchical linear models
#' 
#' @param model fitted model of class \code{mer}
#' @param group variable used to define the group for which cases will be deleted.
#'   If \code{group = NULL}, then individual cases will be deleted.
#' @param delete index of individual cases to be deleted. For higher level units
#'   specified in this manner, the \code{group} parameter must also be specified.
#'   If \code{case = NULL} then all cases are iteratively deleted.
covratio.mer <- function(model, group = NULL, delete = NULL) {
  if(!is(model, "mer")) stop("model must be of class 'mer'")
  if(!is.null(group)) {
    if(!group %in% names(getME(model, "flist"))) {
      stop(paste(group, "is not a valid grouping factor for this model."))
    }
  }
  if(!model@dims["LMM"]){
    stop("covratio is currently not implemented for GLMMs.")
  }
  
  # Extract key pieces of the model
  mats <- .mer_matrices(model)
  
  if( !is.null(group) ){
    grp.names <- unique( mats$flist[, group] )
    
    if( is.null(delete) ){
      del.index <- lapply(1:mats$ngrps[group], 
                          function(x) {
                            ind <- which(mats$flist[, group] == grp.names[x]) - 1
                          })
    } else{
      del.index <- list( which(mats$flist[, group] == delete) - 1 )
    }
  } else{
    if( is.null(delete) ){
      del.index <- split(0:(mats$n-1), 0:(mats$n-1))
    } else { del.index <- list( delete - 1 ) }
  }
  
  res <- .Call("covratioCalc", index = del.index, X_ = mats$X, P_ = mats$P, 
               Vinv_ = as.matrix(mats$Vinv), XVXinv_ = as.matrix(mats$XVXinv), 
               PACKAGE = "HLMdiag")
  
  return(res)
}

#' COVRATIO for mixed/hierarchical linear models
#' 
#' @param model fitted model of class \code{mer}
#' @param group variable used to define the group for which cases will be deleted.
#'   If \code{group = NULL}, then individual cases will be deleted.
#' @param delete index of individual cases to be deleted. For higher level units
#'   specified in this manner, the \code{group} parameter must also be specified.
#'   If \code{case = NULL} then all cases are iteratively deleted.
covtrace.mer <- function(model, group = NULL, delete = NULL) {
  if(!is(model, "mer")) stop("model must be of class 'mer'")
  if(!is.null(group)) {
    if(!group %in% names(getME(model, "flist"))) {
      stop(paste(group, "is not a valid grouping factor for this model."))
    }
  }
  if(!model@dims["LMM"]){
    stop("covtrace is currently not implemented for GLMMs.")
  }
  
  # Extract key pieces of the model
  mats <- .mer_matrices(model)
  
  if( !is.null(group) ){
    grp.names <- unique( mats$flist[, group] )
    
    if( is.null(delete) ){
      del.index <- lapply(1:mats$ngrps[group], 
                          function(x) {
                            ind <- which(mats$flist[, group] == grp.names[x]) - 1
                          })
    } else{
      del.index <- list( which(mats$flist[, group] == delete) - 1 )
    }
  } else{
    if( is.null(delete) ){
      del.index <- split(0:(mats$n-1), 0:(mats$n-1))
    } else { del.index <- list( delete - 1 ) }
  }
  
  res <- .Call("covtraceCalc", index = del.index, X_ = mats$X, P_ = mats$P, 
               Vinv_ = as.matrix(mats$Vinv), XVXinv_ = as.matrix(mats$XVXinv), 
               PACKAGE = "HLMdiag")
  
  return(res)
}

#' Relative variance change for mer objects
#' 
#' model object of class \code{mer}
rvc.mer <- function(model, group = NULL, delete = NULL) {
    delete <- case_delete(model, group = group, type = "varcomp", delete = delete)
    return( rvc(delete) )
}