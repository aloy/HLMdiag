#' @export
leverage <- function(model, ...){
  UseMethod("leverage", model)
}

#' @export
leverage.default <- function(model, ...){
  stop(paste("there is no leverage() method for models of class",
             paste(class(model), collapse=", ")))
}

#' @export
covratio <- function(model, ...){
  UseMethod("covratio", model)
}

#' @export
covratio.default <- function(model, ...){
  stop(paste("there is no covratio() method for models of class",
             paste(class(model), collapse=", ")))
}

#' @export
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

#' @export
covtrace <- function(model, ...){
  UseMethod("covtrace", model)
}

#' @export
covtrace.default <- function(model, ...){
  stop(paste("there is no covtrace() method for models of class",
             paste(class(object), collapse=", ")))
}

#' @export
mdffits <- function(model, ...){
  UseMethod("mdffits", model)
}

#' @export
mdffits.default <- function(model, ...){
  stop(paste("there is no mdffits() method for models of class",
             paste(class(model), collapse=", ")))
}

#' @export
rvc <- function(model, ...){
  UseMethod("rvc", model)
}

#' @export
rvc.default <- function(model, ...){
  stop(paste("there is no rvc() method for models of class",
             paste(class(model), collapse=", ")))
}


#' Leverage for mixed/hierarchical linear models
#' 
#' This function calculates the leverage of
#' a mixed/hierarchical model fit by \code{lmer}. 
#' 
#' @export
#' @method leverage mer
#' @S3method leverage mer
#' @aliases leverage
#' @param model fitted model of class \code{mer}
#' @param the level at which the leverage should be calculated; either
#'   1 for observation level leverage or the name of the grouping factor 
#'   (as defined in \code{flist} of the \code{mer} object) for group level
#'   leverage. \code{leverage} assumes that the grouping factors are unique;
#'   thus, if IDs are repeated within each unit, unique IDs must be generated 
#'   by the user prior to use of \code{leverage}.
#' @references 
#'   Nobre, J. S., & Singer, J. M. (2011). 
#'   Leverage analysis for linear mixed models. 
#'   Journal of Applied Statistics, 38(5), 1063–1072.
#'   
#'   Demidenko, E., & Stukel, T. A. (2005). 
#'   Influence analysis for linear mixed-effects models. 
#'   Statistics in Medicine, 24(6), 893–909.
#' @author Adam Loy \email{aloy@@iastate.edu}
#' @export
#' @seealso \code{\link{cooks.distance.mer}}, \code{\link{mdffits.mer}},
#' \code{\link{covratio.mer}}, \code{\link{covtrace.mer}}, \code{\link{rvc.mer}}  
leverage.mer <- function(model, level) {
  if(!is(model, "mer")) stop("model must be of class 'mer'")
  if(model@dims[["nest"]] == 0) {
    stop("leverage.mer has not yet been implemented for models with 
         crossed random effects")
  }
  if(!level %in% c( 1, names(getME(model, "flist")))) {
    stop("level can only be 1 or a grouping factor from the fitted model.")
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
  
  if(level == 1) {
    lev1 <- data.frame(fixef = diag.H1, ranef =  diag.H2)
  } else {
    flist   <- data.frame( getME(model, "flist")[, level] )
    grp.lev <- data.frame( fixef = aggregate(diag.H1, flist, mean)[,2], 
                           ranef = aggregate(diag.H2, flist, mean)[,2] )
  }
  
  if(level == 1) return(lev1)
  if(level != 1) return(grp.lev)
}


#' Cook's distance for mixed/hierarchical linear models
#'
#' This function calculates Cook's distance for the fixed effects parameters
#' for a mixed/hierarchical model fit by \code{lmer}. 
#'
#'@export
#'@method cooks.distance mer
#'@S3method cooks.distance mer
#'@aliases cooks.distance
#'@param model fitted model of class \code{mer}
#'@param group variable used to define the group for which cases will be
#'deleted.  If \code{group = NULL}, then individual cases will be deleted.
#'@param delete index of individual cases to be deleted.  For higher level
#'units specified in this manner, the \code{group} parameter must also be
#'specified.  If \code{case = NULL} then all cases are iteratively deleted.
#'@author Adam Loy \email{aloy@@iastate.edu}
#'@references
#' Christensen, R., Pearson, L., & Johnson, W. (1992). 
#' Case-deletion diagnostics for mixed models. \emph{Technometrics}, 34(1), 38–45.
#'   
#' Schabenberger, O. (2004),``Mixed Model Influence Diagnostics,''
#' in \emph{Proceedings of the Twenty-Ninth SAS Users Group International Conference},
#' SAS Users Group International.
#' 
#'@keywords models regression
#' @seealso \code{\link{leverage.mer}}, \code{\link{mdffits.mer}},
#' \code{\link{covratio.mer}}, \code{\link{covtrace.mer}}, \code{\link{rvc.mer}}  
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
#' This function calculates multivariate DFFITS for the fixed effects parameters
#' for a mixed/hierarchical model fit by \code{lmer}.
#' 
#' @export
#' @method mdffits mer
#' @S3method mdffits mer
#' @aliases mdffits
#' @param model fitted model of class \code{mer}
#' @param group variable used to define the group for which cases will be deleted.
#'   If \code{group = NULL}, then individual cases will be deleted.
#' @param delete index of individual cases to be deleted. For higher level units
#'   specified in this manner, the \code{group} parameter must also be specified.
#'   If \code{case = NULL} then all cases are iteratively deleted.
#'@author Adam Loy \email{aloy@@iastate.edu}
#'@references
#'   
#' Schabenberger, O. (2004),``Mixed Model Influence Diagnostics,''
#' in \emph{Proceedings of the Twenty-Ninth SAS Users Group International Conference},
#' SAS Users Group International.
#' 
#'@keywords models regression
#' @seealso \code{\link{leverage.mer}}, \code{\link{mdffits.mer}},
#' \code{\link{covratio.mer}}, \code{\link{covtrace.mer}}, \code{\link{rvc.mer}}
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



#'COVRATIO for mixed/hierarchical linear models
#'
#' This function calculates COVRATIO for the fixed effects parameters
#' for a mixed/hierarchical model fit by \code{lmer}.
#'
#'@export
#'@method covratio mer
#'@S3method covratio mer
#'@aliases covratio
#'@param model fitted model of class \code{mer}
#'@param group variable used to define the group for which cases will be
#'deleted.  If \code{group = NULL}, then individual cases will be deleted.
#'@param delete index of individual cases to be deleted.  For higher level
#'units specified in this manner, the \code{group} parameter must also be
#'specified.  If \code{case = NULL} then all cases are iteratively deleted.
#' @return If \code{delete = NULL} then a vector corresponding to each deleted
#' observation/group is returned.
#' 
#' If \code{delete} is specified then a single value is returned corresponding
#' to the deleted subset specified.
#'@author Adam Loy \email{aloy@@iastate.edu}
#'@references
#' Christensen, R., Pearson, L., & Johnson, W. (1992). 
#' Case-deletion diagnostics for mixed models. \emph{Technometrics}, 34(1), 38–45.
#'   
#' Schabenberger, O. (2004),``Mixed Model Influence Diagnostics,''
#' in \emph{Proceedings of the Twenty-Ninth SAS Users Group International Conference},
#' SAS Users Group International.
#' 
#'@keywords models regression
#' @seealso \code{\link{leverage.mer}}, \code{\link{cooks.distance.mer}}
#' \code{\link{mdffits.mer}},
#'  \code{\link{covtrace.mer}}, \code{\link{rvc.mer}}
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



#'COVTRACE for mixed/hierarchical linear models
#'
#' This function calculates COVTRACE for the fixed effects parameters
#' for a mixed/hierarchical model fit by \code{lmer}.
#'
#'@export
#'@method covtrace mer
#'@S3method covtrace mer
#'@aliases covtrace
#'@param model fitted model of class \code{mer}
#'@param group variable used to define the group for which cases will be
#'deleted.  If \code{group = NULL}, then individual cases will be deleted.
#'@param delete index of individual cases to be deleted.  For higher level
#'units specified in this manner, the \code{group} parameter must also be
#'specified.  If \code{case = NULL} then all cases are iteratively deleted.
#' @return If \code{delete = NULL} then a vector corresponding to each deleted
#' observation/group is returned.
#' 
#' If \code{delete} is specified then a single value is returned corresponding
#' to the deleted subset specified.
#' 
#'@author Adam Loy \email{aloy@@iastate.edu}
#'@references
#' Christensen, R., Pearson, L., & Johnson, W. (1992). 
#' Case-deletion diagnostics for mixed models. \emph{Technometrics}, 
#' 34(1), 38–45.
#'   
#' Schabenberger, O. (2004),``Mixed Model Influence Diagnostics,''
#' in \emph{Proceedings of the Twenty-Ninth SAS Users Group International Conference},
#' SAS Users Group International.
#' 
#'@keywords models regression
#' @seealso \code{\link{leverage.mer}}, \code{\link{cooks.distance.mer}}, 
#' \code{\link{mdffits.mer}},
#' \code{\link{covratio.mer}}, \code{\link{rvc.mer}}
covtrace.mer <- function(model, group = NULL, delete = NULL) {
  if(!is(model, "mer")) stop("model must be of class 'mer'")
  if(!is.null(group)) {
    if(!group %in% names(getME(model, "flist"))) {
      stop(paste(group, "is not a valid grouping factor for this model."))
    }
  }
  if(!model@dims["LMM"]){
    stop("covtrace is currently not implemented for GLMMs or NLMMs.")
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

#' Relative variance change for mixed/hierarchical linear models
#' 
#' This function calculates the relative variance change (RVC) for
#' mixed/hierarchical linear models fit via \code{lmer}.
#' 
#' @export
#' @method rvc mer
#' @S3method rvc mer
#' @aliases rvc
#'@param model fitted model of class \code{mer}
#'@param group variable used to define the group for which cases will be
#'deleted.  If \code{group = NULL}, then individual cases will be deleted.
#'@param delete index of individual cases to be deleted.  For higher level
#'units specified in this manner, the \code{group} parameter must also be
#'specified.  If \code{case = NULL} then all cases are iteratively deleted.
#' @return If \code{delete = NULL} a matrix with columns corresponding to the variance 
#' components of the model and rows corresponding to the deleted 
#' observation/group is returned. 
#' 
#' If \code{delete} is specified then a named vector is returned.
#' 
#' The residual variance is named \code{sigma2} and the other variance 
#' componenets are named \code{D**} where the trailing digits give the
#' position in the covariance matrix of the random effects.
#' 
#'@author Adam Loy \email{aloy@@iastate.edu}
#'@references
#' Dillane, D. (2005). ``Deletion Diagnostics for the Linear Mixed Model.'' 
#' Ph.D. thesis, Trinity College Dublin
#' 
#' @keywords models regression
#' @seealso \code{\link{leverage.mer}}, 
#' \code{\link{cooks.distance.mer}}, \code{\link{mdffits.mer}},
#' \code{\link{covratio.mer}}, \code{\link{covtrace.mer}}
rvc.mer <- function(model, group = NULL, delete = NULL) {
    delete <- case_delete(model, group = group, type = "varcomp", delete = delete)
    return( rvc(delete) )
}
