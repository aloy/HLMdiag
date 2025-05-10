#' @export
leverage <- function(model, ...){
  UseMethod("leverage", model)
}

#' @export
#' @rdname leverage.mer
#' @method leverage default
leverage.default <- function(model, ...){
  stop(paste("there is no leverage() method for models of class",
             paste(class(model), collapse=", ")))
}

#' @export
covratio <- function(model, ...){
  UseMethod("covratio", model)
}

#' @export
#' @rdname covratio
#' @method covratio default
covratio.default <- function(model, ...){
  stop(paste("there is no covratio() method for models of class",
             paste(class(model), collapse=", ")))
}

#' @export
covtrace <- function(model, ...){
  UseMethod("covtrace", model)
}

#' @export
#' @rdname covratio
#' @method covtrace default
covtrace.default <- function(model, ...){
  stop(paste("there is no covtrace() method for models of class",
             paste(class(model), collapse=", ")))
}

#' @export
mdffits <- function(model, ...){
  UseMethod("mdffits", model)
}

#' @export
#' @rdname cooks.distance
#' @method mdffits default
mdffits.default <- function(model, ...){
  stop(paste("there is no mdffits() method for models of class",
             paste(class(model), collapse=", ")))
}

#' @export
rvc <- function(model, ...){
  UseMethod("rvc", model)
}

#' @export
#' @rdname rvc.mer
#' @method rvc default
rvc.default <- function(model, ...){
  stop(paste("there is no rvc() method for models of class",
             paste(class(model), collapse=", ")))
}


#' Leverage for HLMs
#' 
#' This function calculates the leverage of
#' a hierarchical linear model fit by \code{lmer}. 
#' 
#' @export
#' @method leverage mer
#' @aliases leverage
#' @param model fitted model object of class \code{mer} of \code{lmerMod}
#' @param level the level at which the leverage should be calculated: either
#'   1 for observation level leverage (default) or the name of the grouping factor 
#'   (as defined in \code{flist} of the \code{mer} object) for group level
#'   leverage. \code{leverage} assumes that the grouping factors are unique;
#'   thus, if IDs are repeated within each unit, unique IDs must be generated 
#'   by the user prior to use of \code{leverage}.
#' @param ... do not use
#' @details Demidenko and Stukel (2005) describe leverage for mixed (hierarchical)
#' linear models as being the sum of two components, a leverage associated with the 
#' fixed (\eqn{H_1}) and a leverage associated with the random effects (\eqn{H_2}) where
#' \deqn{H_1 = X (X^\prime V^{-1} X)^{-1} X^\prime V^{-1}}
#' and
#' \deqn{H_2 = ZDZ^{\prime} V^{-1} (I - H_1)}
#' Nobre and Singer (2011) propose using
#' \deqn{H_2^* = ZDZ^{\prime}}
#' as the random effects leverage as it does not rely on the fixed effects.
#' 
#' For individual observations \code{leverage} uses the diagonal elements of the 
#' above matrices as the measure of leverage. For higher-level units, 
#' \code{leverage} uses the mean trace of the above matrices associated with each
#' higher-level unit.
#' @return \code{leverage} returns a data frame with the following columns:
#' \describe{
#'   \item{\code{overall}}{The overall leverage, i.e. \eqn{H = H_1 + H_2}.}
#'   \item{\code{fixef}}{The leverage corresponding to the fixed effects.}
#'   \item{\code{ranef}}{The leverage corresponding to the random effects 
#'     proposed by Demidenko and Stukel (2005).}
#'   \item{\code{ranef.uc}}{The (unconfounded) leverage corresponding to the 
#'     random effects proposed by Nobre and Singer (2011).}
#' }
#' @references 
#'   Demidenko, E., & Stukel, T. A. (2005) 
#'   Influence analysis for linear mixed-effects models. 
#'   \emph{Statistics in Medicine}, \bold{24}(6), 893--909.
#'   
#'   Nobre, J. S., & Singer, J. M. (2011) 
#'   Leverage analysis for linear mixed  models. 
#'   \emph{Journal of Applied Statistics}, \bold{38}(5), 1063--1072.
#' @author Adam Loy \email{loyad01@@gmail.com}
#' @keywords models regression
#' @export
#' @seealso \code{\link{cooks.distance.mer}}, \code{\link{mdffits.mer}},
#' \code{\link{covratio.mer}}, \code{\link{covtrace.mer}}, \code{\link{rvc.mer}}  
#' 
#' @examples
#' data(sleepstudy, package = 'lme4')
#' fm <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' 
#' # Observation level leverage
#' lev1 <- leverage(fm, level = 1)
#' head(lev1)
#' 
#' # Group level leverage
#' lev2 <- leverage(fm, level = "Subject")
#' head(lev2)
leverage.mer <- function(model, level = 1, ...) {
  .Deprecated("leverage.lmerMod")
  if(!is(model, "mer")) stop("model must be of class 'mer'")
  if(model@dims[["nest"]] == 0) {
    stop("leverage.mer has not yet been implemented for models with 
         crossed random effects")
  }
  if(!level %in% c( 1, names(lme4::getME(model, "flist")))) {
    stop("level can only be 1 or a grouping factor from the fitted model.")
  }
  
  if(!model@dims["LMM"]){
    stop("leverage is currently not implemented for GLMMs or NLMMs.")
  }
  
  X <- lme4::getME(model, "X")
  # Z <- BlockZ(model)
  
  n     <- nrow(X)
  nt    <- model@dims[["nt"]]  # number of random-effects terms in the model
  ngrps <- unname( summary(model)@ngrps )
  
  vc   <- lme4::VarCorr(model)
  # D  <- kronecker( Diagonal(ngrps), bdiag(vc) )
  ZDZt <- attr(vc, "sc")^2 * crossprod( lme4::getME(model, "A") )
  R    <- Diagonal( n = n, x = attr(vc, "sc")^2 )
  
  V      <- ZDZt + R
  V.chol <- chol( V )
  Vinv   <- chol2inv( V.chol )
  
  xvix.inv <- attr(vc, "sc")^2 * chol2inv(lme4::getME(model, "RX"))
  
  H1 <- X %*% xvix.inv %*% t(X) %*% Vinv
  H2 <- ZDZt %*% Vinv %*% (Diagonal( n = n ) - H1)
  
  diag.H1 <- diag(H1)
  diag.H2 <- diag(H2)
  diag.H2.uc <- diag(ZDZt) / attr(vc, "sc")^2
  
  if(level == 1) {
    lev1 <- data.frame(overall = diag.H1 + diag.H2, fixef = diag.H1, 
                       ranef =  diag.H2, ranef.uc = diag.H2.uc)
#     class(lev1) <- "leverage"
  } else {
    flist   <- data.frame( lme4::getME(model, "flist")[, level] )
    
    grp.lev.fixef <- aggregate(diag.H1, flist, mean)[,2]
    grp.lev.ranef <- aggregate(diag.H2, flist, mean)[,2]
    grp.lev.ranef.uc <- aggregate(diag.H2.uc, flist, mean)[,2]
    
    grp.lev <- data.frame( overall = grp.lev.fixef + grp.lev.ranef,
                           fixef = grp.lev.fixef, 
                           ranef = grp.lev.ranef,
                           ranef.uc = grp.lev.ranef.uc)
#     class(grp.lev) <- "leverage"
  }
  
  if(level == 1) return(lev1)
  if(level != 1) return(grp.lev)
}

#' @export
#' @rdname leverage.mer
#' @method leverage lmerMod
leverage.lmerMod <- function(model, level = 1, ...) {
  if(!isNestedModel(model)) {
    stop("leverage.mer has not yet been implemented for models with 
         crossed random effects")
  }
  if(!level %in% c( 1, names(lme4::getME(model, "flist")))) {
    stop("level can only be 1 or a grouping factor from the fitted model.")
  }
  
  if(!lme4::isLMM(model)){
    stop("leverage is currently not implemented for GLMMs or NLMMs.")
  }
  
  X <- lme4::getME(model, "X")
  # Z <- BlockZ(model)
  
  n     <- nrow(X)
#  nt    <- model@dims[["nt"]]  # number of random-effects terms in the model
  ngrps <- unname( summary(model)$ngrps )
  
  vc   <- lme4::VarCorr(model)
  # D  <- kronecker( Diagonal(ngrps), bdiag(vc) )
  ZDZt <- attr(vc, "sc")^2 * crossprod( lme4::getME(model, "A") )
  R    <- Diagonal( n = n, x = attr(vc, "sc")^2 )
  
  V      <- ZDZt + R
  V.chol <- chol( V )
  Vinv   <- chol2inv( V.chol )
  
  xvix.inv <- attr(vc, "sc")^2 * chol2inv(lme4::getME(model, "RX"))
  
  H1 <- X %*% xvix.inv %*% t(X) %*% Vinv
  H2 <- ZDZt %*% Vinv %*% (Diagonal( n = n ) - H1)
  
  diag.H1 <- diag(H1)
  diag.H2 <- diag(H2)
  diag.H2.uc <- diag(ZDZt) / attr(vc, "sc")^2
  
  if(level == 1) {
    lev1 <- data.frame(overall = diag.H1 + diag.H2, fixef = diag.H1, 
                       ranef =  diag.H2, ranef.uc = diag.H2.uc)
    #     class(lev1) <- "leverage"
  } else {
    flist   <- data.frame( lme4::getME(model, "flist")[[level]] )
    
    grp.lev.fixef <- aggregate(diag.H1, flist, mean)[,2]
    grp.lev.ranef <- aggregate(diag.H2, flist, mean)[,2]
    grp.lev.ranef.uc <- aggregate(diag.H2.uc, flist, mean)[,2]
    
    grp.lev <- data.frame( overall = grp.lev.fixef + grp.lev.ranef,
                           fixef = grp.lev.fixef, 
                           ranef = grp.lev.ranef,
                           ranef.uc = grp.lev.ranef.uc)
    #     class(grp.lev) <- "leverage"
  }
  
  if(level == 1) return(lev1)
  if(level != 1) return(grp.lev)
  }

#' Influence on fixed effects of HLMs
#'
#' These functions calculate measures of the change in the fixed effects
#' estimates based on the deletion of an observation, or group of 
#' observations, for a hierarchical linear model fit using \code{lmer}.
#' 
#' @details
#' Both Cook's distance and MDFFITS measure the change in the 
#' fixed effects estimates based on the deletion of a subset of observations. 
#' The key difference between the two diagnostics is that Cook's distance uses
#' the covariance matrix for the fixed effects from the original
#' model while MDFFITS uses the covariance matrix from the deleted 
#' model. 
#' 
#' @note
#' Because MDFFITS requires the calculation of the covariance matrix
#' for the fixed effects for every model, it will be slower.
#' 
#' @return Both functions return a numeric vector (or single value if 
#' \code{delete} has been specified) as the default. If \code{include.attr = TRUE}, 
#' then a tibble is returned. The first column consists of the Cook's distance or
#' MDFFITS values, and the later columns capture the difference between the full
#' and deleted parameter estimates. 
#'
#'@export
#'@rdname cooks.distance
#'@method cooks.distance mer
#'@aliases cooks.distance
#'@param model fitted model of class \code{mer} or \code{lmerMod}
#'@param level variable used to define the group for which cases will be
#'deleted.  If \code{level = 1} (default), then individual cases will be deleted.
#'@param delete index of individual cases to be deleted. To delete specific 
#' observations the row number must be specified. To delete higher level
#'units the group ID and \code{group} parameter must be specified.
#' If \code{delete = NULL} then all cases are iteratively deleted.
#'@param include.attr logical value determining whether the difference between 
#'the full and deleted parameter estimates should be included. If \code{FALSE} 
#'(default), a numeric vector of Cook's distance or MDFFITS is returned. 
#'If \code{TRUE}, a tibble with the Cook's distance or MDFFITS values in the 
#'first column and the parameter differences in the remaining columns is returned. 
#' @param ... do not use
#'@author Adam Loy \email{loyad01@@gmail.com}
#'@references
#' Christensen, R., Pearson, L., & Johnson, W. (1992) 
#' Case-deletion diagnostics for mixed models. \emph{Technometrics}, \bold{34}, 
#' 38--45.
#'   
#' Schabenberger, O. (2004) Mixed Model Influence Diagnostics,
#' in \emph{Proceedings of the Twenty-Ninth SAS Users Group International Conference},
#' SAS Users Group International.
#' 
#'@keywords models regression
#' @seealso \code{\link{leverage.mer}},
#' \code{\link{covratio.mer}}, \code{\link{covtrace.mer}}, \code{\link{rvc.mer}}
#' @examples 
#' data(sleepstudy, package = 'lme4')
#' ss <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' 
#' # Cook's distance for individual observations
#' ss.cd.lev1 <- cooks.distance(ss)
#' 
#' # Cook's distance for each Subject
#' ss.cd.subject <- cooks.distance(ss, level = "Subject")
#' 
#' \dontrun{
#' data(Exam, package = 'mlmRev')
#' fm <- lme4::lmer(normexam ~ standLRT * schavg + (standLRT | school), Exam)
#' 
#' # Cook's distance for individual observations
#' cd.lev1 <- cooks.distance(fm)
#' 
#' # Cook's distance for each school
#' cd.school <- cooks.distance(fm, level = "school")
#' 
#' # Cook's distance when school 1 is deleted
#' cd.school1 <- cooks.distance(fm, level = "school", delete = 1)
#' 
#' }
#' 
cooks.distance.mer <- function(model, level = 1, delete = NULL, ...) {
  .Deprecated("cooks.distance.lmerMod")
  if (hasArg(group)) {
    group <- NULL
    warning("group is not a valid argument for this function. As of version 0.4.0, group has been replaced by level.")
  }
  
  if(!is(model, "mer")) stop("model must be of class 'mer'")
  if(level != 1) {
    if(!level %in% names(lme4::getME(model, "flist"))) {
      stop(paste(level, "is not a valid grouping factor for this model."))
    }
  }
  if(!model@dims["LMM"]){
    stop("cooks.distance is currently not implemented for GLMMs or NLMMs.")
  }
  
  XVXinv <- X <- Vinv <- Y <- NULL # Make codetools happy
  
  # Extract key pieces of the model
  mats <- .mer_matrices(model)
  
  betaHat <- with(mats, XVXinv %*% t(X) %*% Vinv %*% Y)
  
  # Obtaining the building blocks
  if(level == 1 & is.null(delete)) {
    calc.cooksd <- .Call("cooksdObs", y_ = mats$Y, X_ = as.matrix(mats$X), 
                         Vinv_ = as.matrix(mats$Vinv), 
                         XVXinv_ = as.matrix(mats$XVXinv), 
                         beta_ = as.matrix(betaHat), PACKAGE = "HLMdiag")
    res <- calc.cooksd[[1]]
    attr(res, "beta_cdd") <- calc.cooksd[[2]]
  }
  
  else{
    e <- with(mats, Y - X %*% betaHat)
    
    if( level != 1){
      grp.names <- unique( mats$flist[, level] )
      
      if( is.null(delete) ){
        del.index <- lapply(1:mats$ngrps[level], 
                            function(x) {
                              ind <- which(mats$flist[, level] == grp.names[x]) - 1
                            })
      } else{
        del.index <- list( which(mats$flist[, level] %in% delete) - 1 )
      }
    } else{
      del.index <- list( delete - 1 )
    }
    
    calc.cooksd <- .Call("cooksdSubset", index = del.index, 
                         X_ = as.matrix(mats$X), 
                         P_ = as.matrix(mats$P), 
                         Vinv_ = as.matrix(mats$Vinv), 
                         XVXinv_ = as.matrix(mats$XVXinv), 
                         e_ = as.numeric(e), PACKAGE = "HLMdiag")
    
    res <- calc.cooksd[[1]]
    attr(res, "beta_cdd") <- calc.cooksd[[2]] 
  }
  
  class(res) <- "fixef.dd"
  return(res)
}

#' @export
#' @rdname cooks.distance
#' @method cooks.distance lmerMod
cooks.distance.lmerMod <- function(model, level = 1, delete = NULL, include.attr = FALSE, ...) {
  if (hasArg(group)) {
    group <- NULL
    warning("group is not a valid argument for this function. As of version 0.4.0, group has been replaced by level.")
  }
  
  if(level != 1) {
    if(!level %in% names(lme4::getME(model, "flist"))) {
      stop(paste(level, "is not a valid grouping factor for this model."))
    }
  }
  if(!lme4::isLMM(model)){
    stop("cooks.distance is currently not implemented for GLMMs or NLMMs.")
  }
  
  XVXinv <- X <- Vinv <- Y <- NULL # Make codetools happy
  
  # Extract key pieces of the model
  mats <- .lmerMod_matrices(model)
  
  betaHat <- with(mats, XVXinv %*% t(X) %*% Vinv %*% Y)
  beta_names <- rownames(betaHat)
  
  # Obtaining the building blocks
  if(level == 1 & is.null(delete)) {
    calc.cooksd <- .Call("cooksdObs", y_ = mats$Y, X_ = as.matrix(mats$X), 
                         Vinv_ = as.matrix(mats$Vinv), 
                         XVXinv_ = as.matrix(mats$XVXinv), 
                         beta_ = as.matrix(betaHat), PACKAGE = "HLMdiag")
    
    cooksd <- calc.cooksd[[1]]
    betas <- calc.cooksd[[2]]
    
    if(!include.attr) { 
      return(cooksd)
    }
    else {
      beta_matrix <- matrix(unlist(betas), nrow = length(betas), byrow = TRUE) 
      colnames(beta_matrix) <- beta_names
      cook.tbl <- tibble::as_tibble(cbind(cooksd, beta_matrix))
      # nbetas <- ncol(cook.tbl) - 1 
      # for (i in 1:nbetas) {
      #   names(cook.tbl)[i+1] <- stringr::str_c("beta", i, sep = "_")
      # }
      return(cook.tbl)
    }
  }
  
  else{
    e <- with(mats, Y - X %*% betaHat)
    
    if( level != 1 ){
      grp.names <- unique( mats$flist[[level]] )
      
      if( is.null(delete) ){
        del.index <- lapply(1:mats$ngrps[level], 
                            function(x) {
                              ind <- which(mats$flist[[level]] == grp.names[x]) - 1
                            })
      } else{
        del.index <- list( which(mats$flist[[level]] %in% delete) - 1 )
      }
    } else{
      del.index <- list( delete - 1 )
    }
    
    calc.cooksd <- .Call("cooksdSubset", index = del.index, 
                         X_ = as.matrix(mats$X), 
                         P_ = as.matrix(mats$P), 
                         Vinv_ = as.matrix(mats$Vinv), 
                         XVXinv_ = as.matrix(mats$XVXinv), 
                         e_ = as.numeric(e), PACKAGE = "HLMdiag")
    
    cooksd <- calc.cooksd[[1]]
    betas <- calc.cooksd[[2]]
    
    if(!include.attr) { 
      return(cooksd)
    }
    else { 
      beta_matrix <- matrix(unlist(betas), nrow = length(betas), byrow = TRUE) 
      colnames(beta_matrix) <- beta_names
      cook.tbl <- tibble::as_tibble(cbind(cooksd, beta_matrix))
      # nbetas <- ncol(cook.tbl) - 1 
      # for (i in 1:nbetas) {
      #   names(cook.tbl)[i+1] <- stringr::str_c("beta", i, sep = "_")
      # }
      return(cook.tbl)
    }
  }
}

#' @export
#' @rdname cooks.distance
#' @method cooks.distance lme
cooks.distance.lme <- function(model, level = 1, delete = NULL, include.attr = FALSE, ...) {
  if (hasArg(group)) {
    group <- NULL
    warning("group is not a valid argument for this function. As of version 0.4.0, group has been replaced by level.")
  }
  
  if(level != 1) {
    if(!level %in% names(model$groups)) {
      stop(paste(level, "is not a valid grouping factor for this model."))
    }
  }
  if (any("nlme" == class(model))) 
    stop("not implemented for \"nlme\" objects")
  
  XVXinv <- X <- Vinv <- Y <- NULL # Make codetools happy
  
  # Extract key pieces of the model
  mats <- .lme_matrices(model)
  
  betaHat <- with(mats, XVXinv %*% t(X) %*% Vinv %*% Y)
  beta_names <- rownames(betaHat)
  
  # Obtaining the building blocks
  if(level == 1 & is.null(delete)) {
    calc.cooksd <- .Call("cooksdObs", y_ = mats$Y, X_ = as.matrix(mats$X), 
                         Vinv_ = as.matrix(mats$Vinv), 
                         XVXinv_ = as.matrix(mats$XVXinv), 
                         beta_ = as.matrix(betaHat), PACKAGE = "HLMdiag")
    
    cooksd <- calc.cooksd[[1]]
    betas <- calc.cooksd[[2]]
    
    if(!include.attr) { 
      return(cooksd)
    }
    else {
      beta_matrix <- matrix(unlist(betas), nrow = length(betas), byrow = TRUE) 
      colnames(beta_matrix) <- beta_names
      cook.tbl <- tibble::as_tibble(cbind(cooksd, beta_matrix))
      # nbetas <- ncol(cook.tbl) - 1 
      # for (i in 1:nbetas) {
      #   names(cook.tbl)[i+1] <- stringr::str_c("beta", i, sep = "_")
      # }
      return(cook.tbl)
    }
  }

  else{
    e <- with(mats, Y - X %*% betaHat)
    
    if( level != 1 ){
      grp.names <- unique( mats$flist[[level]] )
      
      if( is.null(delete) ){
        del.index <- lapply(1:mats$ngrps[level], 
                            function(x) {
                              ind <- which(mats$flist[[level]] == grp.names[x]) - 1
                            })
      } else{
        del.index <- list( which(mats$flist[[level]] %in% delete) - 1 )
      }
    } else{
      del.index <- list( delete - 1 )
    }
    
    calc.cooksd <- .Call("cooksdSubset", index = del.index, 
                         X_ = as.matrix(mats$X), 
                         P_ = as.matrix(mats$P), 
                         Vinv_ = as.matrix(mats$Vinv), 
                         XVXinv_ = as.matrix(mats$XVXinv), 
                         e_ = as.numeric(e), PACKAGE = "HLMdiag")
    
    cooksd <- calc.cooksd[[1]]
    betas <- calc.cooksd[[2]]
    
    if(!include.attr) { 
      return(cooksd)
    }
    else { 
      beta_matrix <- matrix(unlist(betas), nrow = length(betas), byrow = TRUE) 
      colnames(beta_matrix) <- beta_names
      cook.tbl <- tibble::as_tibble(cbind(cooksd, beta_matrix))
      # nbetas <- ncol(cook.tbl) - 1 
      # for (i in 1:nbetas) {
      #   names(cook.tbl)[i+1] <- stringr::str_c("beta", i, sep = "_")
      # }
      return(cook.tbl)
    }
  }
}

print.fixef.dd <- function(x, ...) {
  attributes(x) <- NULL
  print(x, ...)
}

print.vcov.dd <- function(x, ...) { print(unclass(x), ...) }


#' @export
#' @rdname cooks.distance
#' @method mdffits mer
#' @aliases mdffits
#' @param object fitted object of class \code{mer} or \code{lmerMod}
#' @examples
#' 
#' # MDFFITS  for individual observations
#' ss.m1 <- mdffits(ss)
#' 
#' # MDFFITS for each Subject
#' ss.m.subject <- mdffits(ss, level = "Subject")
#' 
#' \dontrun{  
#' 
#' # MDFFITS  for individual observations
#' m1 <- mdffits(fm)
#' 
#' # MDFFITS for each school
#' m.school <- mdffits(fm, level = "school")
#' }
mdffits.mer <- function(object, level = 1, delete = NULL, ...) {
  .Deprecated("mdffits.lmerMod")
  if (hasArg(group)) {
    group <- NULL
    warning("group is not a valid argument for this function. As of version 0.4.0, group has been replaced by level.")
  }
  
  if(!is(object, "mer")) stop("object must be of class 'mer'")
  if(level != 1) {
    if(!level %in% names(lme4::getME(object, "flist"))) {
      stop(paste(level, "is not a valid grouping factor for this model."))
    }
  }
  if(!object@dims["LMM"]){
    stop("mdffits is currently not implemented for GLMMs or NLMMs.")
  }
  
  XVXinv <- X <- Vinv <- Y <- NULL # Make codetools happy
  
  # Extract key pieces of the model
  mats <- .mer_matrices(object)
  
  betaHat <- with(mats, XVXinv %*% t(X) %*% Vinv %*% Y)
  e <- with(mats, Y - X %*% betaHat)
  
  if( level != 1 ){
    grp.names <- unique( mats$flist[, level] )
    
    if( is.null(delete) ){
      del.index <- lapply(1:mats$ngrps[level], 
                          function(x) {
                            ind <- which(mats$flist[, level] == grp.names[x]) - 1
                          })
    } else{
      del.index <- list( which(mats$flist[, level] %in% delete) - 1 )
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
  
  class(res) <- "fixef.dd"
  return(res)
}


#' @export
#' @rdname cooks.distance
#' @method mdffits lmerMod
mdffits.lmerMod <- function(model, level = 1, delete = NULL, include.attr = FALSE, ...) {
  if (hasArg(group)) {
    group <- NULL
    warning("group is not a valid argument for this function. As of version 0.4.0, group has been replaced by level.")
  }
  
  if(level != 1) {
    if(!level %in% names(lme4::getME(model, "flist"))) {
      stop(paste(level, "is not a valid grouping factor for this model."))
    }
  }
  if(!lme4::isLMM(model)){
    stop("mdffits is currently not implemented for GLMMs or NLMMs.")
  }
  
  XVXinv <- X <- Vinv <- Y <- NULL # Make codetools happy
  
  # Extract key pieces of the model
  mats <- .lmerMod_matrices(model)
  
  betaHat <- with(mats, XVXinv %*% t(X) %*% Vinv %*% Y)
  beta_names <- rownames(betaHat)
  e <- with(mats, Y - X %*% betaHat)
  
  if( level != 1 ){
    grp.names <- unique( mats$flist[[level]] )
    
    if( is.null(delete) ){
      del.index <- lapply(1:mats$ngrps[level], 
                          function(x) {
                            ind <- which(mats$flist[[level]] == grp.names[x]) - 1
                          })
    } else{
      del.index <- list( which(mats$flist[[level]] %in% delete) - 1 )
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
  mdffits.vals <- calc.mdffits[[1]]
  betas <- calc.mdffits[[2]]
  
  if(!include.attr) { 
    return(mdffits.vals)
  }
  else { 
    beta_matrix <- matrix(unlist(betas), nrow = length(betas), byrow = TRUE) 
    colnames(beta_matrix) <- beta_names
    mdffits.tbl <- tibble::as_tibble(cbind(mdffits = mdffits.vals, beta_matrix))
    # names(mdffits.tbl)[1] <- "mdffits"
    # nbetas <- ncol(mdffits.tbl) - 1 
    # for (i in 1:nbetas) {
    #   names(mdffits.tbl)[i+1] <- stringr::str_c("beta", i, sep = "_")
    # }
    return(mdffits.tbl)
  }
}


#' @export
#' @rdname cooks.distance
#' @method mdffits lme
mdffits.lme <- function(model, level = 1, delete = NULL, include.attr = FALSE, ...) {
  if (hasArg(group)) {
    group <- NULL
    warning("group is not a valid argument for this function. As of version 0.4.0, group has been replaced by level.")
  }
  
  if(level != 1) {
    if(!level %in% names(model$groups)) {
      stop(paste(level, "is not a valid grouping factor for this model."))
    }
  }
  if (any("nlme" == class(model))) 
    stop("not implemented for \"nlme\" objects")
  
  
  XVXinv <- X <- Vinv <- Y <- NULL # Make codetools happy
  
  # Extract key pieces of the model
  mats <- .lme_matrices(model)
  
  betaHat <- with(mats, XVXinv %*% t(X) %*% Vinv %*% Y)
  beta_names <- rownames(betaHat)
  e <- with(mats, Y - X %*% betaHat)
  
  if( level != 1 ){
    grp.names <- unique( mats$flist[[level]] )
    
    if( is.null(delete) ){
      del.index <- lapply(1:mats$ngrps[level], 
                          function(x) {
                            ind <- which(mats$flist[[level]] == grp.names[x]) - 1
                          })
    } else{
      del.index <- list( which(mats$flist[[level]] %in% delete) - 1 )
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
  mdffits.vals <- calc.mdffits[[1]]
  betas <- calc.mdffits[[2]]
  
  if(!include.attr) { 
    return(mdffits.vals)
  }
  else { 
    beta_matrix <- matrix(unlist(betas), nrow = length(betas), byrow = TRUE) 
    colnames(beta_matrix) <- beta_names
    mdffits.tbl <- tibble::as_tibble(cbind(mdffits = mdffits.vals, beta_matrix))
    # names(mdffits.tbl)[1] <- "mdffits"
    # nbetas <- ncol(mdffits.tbl) - 1 
    # for (i in 1:nbetas) {
      # names(mdffits.tbl)[i+1] <- stringr::str_c("beta", i, sep = "_")
    # }
    return(mdffits.tbl)
  }
}

#' Influence on precision of fixed effects in HLMs
#'
#' These functions calculate measures of the change in the covariance
#' matrices for the fixed effects based on the deletion of an
#' observation, or group of observations, for a hierarchical 
#' linear model fit using \code{lmer}.
#' 
#' @details
#'  Both the covariance ratio (\code{covratio}) and the covariance trace
#'  (\code{covtrace}) measure the change in the covariance matrix
#'  of the fixed effects based on the deletion of a subset of observations.
#'  The key difference is how the variance covariance matrices are compared:
#'  \code{covratio} compares the ratio of the determinants while \code{covtrace}
#'  compares the trace of the ratio. 
#'  
#'@export
#' @rdname covratio
#'@method covratio mer
#'@aliases covratio
#'@param model fitted model object of class \code{mer} or \code{lmerMod}
#'@param level variable used to define the group for which cases will be
#'deleted.  If \code{level = 1} (default), then individual cases will be deleted.
#'@param delete index of individual cases to be deleted. To delete specific 
#' observations the row number must be specified. To delete higher level
#'units the group ID and \code{group} parameter must be specified.
#' If \code{delete = NULL} then all cases are iteratively deleted.
#'@param ... do not use
#' @return If \code{delete = NULL} then a vector corresponding to each deleted
#' observation/group is returned.
#' 
#' If \code{delete} is specified then a single value is returned corresponding
#' to the deleted subset specified.
#'@author Adam Loy \email{loyad01@@gmail.com}
#'@references
#' Christensen, R., Pearson, L., & Johnson, W. (1992) 
#' Case-deletion diagnostics for mixed models. \emph{Technometrics}, \bold{34}(1), 
#' 38--45.
#'   
#' Schabenberger, O. (2004) Mixed Model Influence Diagnostics,
#' in \emph{Proceedings of the Twenty-Ninth SAS Users Group International Conference},
#' SAS Users Group International.
#' 
#'@keywords models regression
#' @seealso \code{\link{leverage.mer}}, \code{\link{cooks.distance.mer}}
#' \code{\link{mdffits.mer}}, \code{\link{rvc.mer}}
#'  
#' @examples
#' 
#' data(sleepstudy, package = 'lme4')
#' ss <- lme4::lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
#' 
#' # covratio for individual observations
#' ss.cr1 <- covratio(ss)
#' 
#' # covratio for subject-level deletion
#' ss.cr2 <- covratio(ss, level = "Subject")
#' 
#' \dontrun{
#' ## A larger example
#' data(Exam, package = 'mlmRev')
#' fm <- lme4::lmer(normexam ~ standLRT * schavg + (standLRT | school), data = Exam)
#' 
#' # covratio for individual observations
#' cr1 <- covratio(fm)
#' 
#' # covratio for school-level deletion
#' cr2 <- covratio(fm, level = "school")
#' }
covratio.mer <- function(model, level = 1, delete = NULL, ...) {
  .Deprecated("covratio.lmerMod")
  if(!is(model, "mer")) stop("model must be of class 'mer'")
  
  if (hasArg(group)) {
    group <- NULL
    warning("group is not a valid argument for this function. As of version 0.4.0, group has been replaced by level.")
  }
  
  if(level != 1) {
    if(!level %in% names(lme4::getME(model, "flist"))) {
      stop(paste(level, "is not a valid grouping factor for this model."))
    }
  }
  if(!model@dims["LMM"]){
    stop("covratio is currently not implemented for GLMMs or NLMMs.")
  }
  
  # Extract key pieces of the model
  mats <- .mer_matrices(model)
  
  if( level != 1 ){
    grp.names <- unique( mats$flist[, level] )
    
    if( is.null(delete) ){
      del.index <- lapply(1:mats$ngrps[level], 
                          function(x) {
                            ind <- which(mats$flist[, level] == grp.names[x]) - 1
                          })
    } else{
      del.index <- list( which(mats$flist[, level] %in% delete) - 1 )
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


#'@export
#'@rdname covratio
#'@method covratio lmerMod
covratio.lmerMod <- function(model, level = 1, delete = NULL, ...) {
  if (hasArg(group)) {
    group <- NULL
    warning("group is not a valid argument for this function. As of version 0.4.0, group has been replaced by level.")
  }
  
  if(level != 1) {
    if(!level %in% names(lme4::getME(model, "flist"))) {
      stop(paste(level, "is not a valid grouping factor for this model."))
    }
  }
  if(!lme4::isLMM(model)){
    stop("covratio is currently not implemented for GLMMs or NLMMs.")
  }
  
  # Extract key pieces of the model
  mats <- .lmerMod_matrices(model)
  
  if( level != 1 ){
    grp.names <- unique( mats$flist[[level]] )
    
    if( is.null(delete) ){
      del.index <- lapply(1:mats$ngrps[level], 
                          function(x) {
                            ind <- which(mats$flist[[level]] == grp.names[x]) - 1
                          })
    } else{
      del.index <- list( which(mats$flist[[level]] %in% delete) - 1 )
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


#'@export
#'@rdname covratio
#'@method covratio lme
covratio.lme <- function(model, level = 1, delete = NULL, ...) {
  if (hasArg(group)) {
    group <- NULL
    warning("group is not a valid argument for this function. As of version 0.4.0, group has been replaced by level.")
  }
  
  if(level != 1) {
    if(!level %in% names(model$groups)) {
      stop(paste(level, "is not a valid grouping factor for this model."))
    }
  }
  if (any("nlme" == class(model))) 
    stop("not implemented for \"nlme\" objects")
  
  # Extract key pieces of the model
  mats <- .lme_matrices(model)
  
  if( level != 1 ){
    grp.names <- unique( mats$flist[[level]] )
    
    if( is.null(delete) ){
      del.index <- lapply(1:mats$ngrps[level], 
                          function(x) {
                            ind <- which(mats$flist[[level]] == grp.names[x]) - 1
                          })
    } else{
      del.index <- list( which(mats$flist[[level]] %in% delete) - 1 )
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


#'@export
#'@rdname covratio
#'@method covtrace mer
#'@aliases covtrace
#' @examples
#' 
#' # covtrace for individual observations
#' ss.ct1 <- covtrace(ss)
#' 
#' # covtrace for subject-level deletion
#' ss.ct2 <- covtrace(ss, level = "Subject")
#' 
#' \dontrun{
#' ## Returning to the larger example
#' # covtrace for individual observations
#' ct1 <- covtrace(fm)
#' 
#' # covtrace for school-level deletion
#' ct2 <- covtrace(fm, level = "school")
#' }
covtrace.mer <- function(model, level = 1, delete = NULL, ...) {
  .Deprecated("covtrace.lmerMod")
  if(!is(model, "mer")) stop("model must be of class 'mer'")
  
  if (hasArg(group)) {
    group <- NULL
    warning("group is not a valid argument for this function. As of version 0.4.0, group has been replaced by level.")
  }
  if(level != 1) {
    if(!level %in% names(lme4::getME(model, "flist"))) {
      stop(paste(level, "is not a valid grouping factor for this model."))
    }
  }
  if(!model@dims["LMM"]){
    stop("covtrace is currently not implemented for GLMMs or NLMMs.")
  }
  
  # Extract key pieces of the model
  mats <- .mer_matrices(model)
  
  if( level != 1 ){
    grp.names <- unique( mats$flist[, level] )
    
    if( is.null(delete) ){
      del.index <- lapply(1:mats$ngrps[level], 
                          function(x) {
                            ind <- which(mats$flist[, level] == grp.names[x]) - 1
                          })
    } else{
      del.index <- list( which(mats$flist[, level] %in% delete) - 1 )
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


#'@export
#'@rdname covratio
#'@method covtrace lmerMod
covtrace.lmerMod <- function(model, level = 1, delete = NULL, ...) {
  if (hasArg(group)) {
    group <- NULL
    warning("group is not a valid argument for this function. As of version 0.4.0, group has been replaced by level.")
  }
  
  if(level != 1) {
    if(!level %in% names(lme4::getME(model, "flist"))) {
      stop(paste(level, "is not a valid grouping factor for this model."))
    }
  }
  if(!lme4::isLMM(model)){
    stop("covtrace is currently not implemented for GLMMs or NLMMs.")
  }
  
  # Extract key pieces of the model
  mats <- .lmerMod_matrices(model)
  
  if( level != 1 ){
    grp.names <- unique( mats$flist[[level]] )
    
    if( is.null(delete) ){
      del.index <- lapply(1:mats$ngrps[level], 
                          function(x) {
                            ind <- which(mats$flist[[level]] == grp.names[x]) - 1
                          })
    } else{
      del.index <- list( which(mats$flist[[level]] %in% delete) - 1 )
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



#'@export
#'@rdname covratio
#'@method covtrace lme
covtrace.lme <- function(model, level = 1, delete = NULL, ...) {
  if (hasArg(group)) {
    group <- NULL
    warning("group is not a valid argument for this function. As of version 0.4.0, group has been replaced by level.")
  }
  
    if(level != 1) {
      if(!level %in% names(model$groups)) {
        stop(paste(level, "is not a valid grouping factor for this model."))
      }
    }
    if (any("nlme" == class(model))) 
      stop("not implemented for \"nlme\" objects")
  
  # Extract key pieces of the model
  mats <- .lme_matrices(model)
  
  if( level != 1 ){
    grp.names <- unique( mats$flist[[level]] )
    
    if( is.null(delete) ){
      del.index <- lapply(1:mats$ngrps[level], 
                          function(x) {
                            ind <- which(mats$flist[[level]] == grp.names[x]) - 1
                          })
    } else{
      del.index <- list( which(mats$flist[[level]] %in% delete) - 1 )
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


#' Relative variance change for HLMs
#' 
#' This function calculates the relative variance change (RVC) of
#' hierarchical linear models fit via \code{lmer}.
#' 
#' @export
#' @method rvc mer
#' @aliases rvc
#'@param model fitted model object of class \code{mer} or \code{lmerMod}
#'@param level variable used to define the group for which cases will be
#'deleted.  If \code{level = 1} (default), then individual cases will be deleted.
#'@param delete index of individual cases to be deleted. To delete specific 
#' observations the row number must be specified. To delete higher level
#'units the group ID and \code{group} parameter must be specified.
#' If \code{delete = NULL} then all cases are iteratively deleted.
#'@param ... do not use
#' @return If \code{delete = NULL} a matrix with columns corresponding to the variance 
#' components of the model and rows corresponding to the deleted 
#' observation/group is returned. 
#' 
#' If \code{delete} is specified then a named vector is returned.
#' 
#' The residual variance is named \code{sigma2} and the other variance 
#' components are named \code{D**} where the trailing digits give the
#' position in the covariance matrix of the random effects.
#' 
#'@author Adam Loy \email{loyad01@@gmail.com}
#'@references
#' Dillane, D. (2005) Deletion Diagnostics for the Linear Mixed Model. 
#' Ph.D. thesis, Trinity College Dublin
#' 
#' @keywords models regression
#' @seealso \code{\link{leverage.mer}}, 
#' \code{\link{cooks.distance.mer}}, \code{\link{mdffits.mer}},
#' \code{\link{covratio.mer}}, \code{\link{covtrace.mer}}
rvc.mer <- function(model, level = 1, delete = NULL, ...) {
  .Deprecated("rvc.lmerMod")
    delete <- case_delete(model, level = 1, type = "varcomp", delete = delete)
    return( rvc(delete) )
}


#' @export
#' @rdname rvc.mer
#' @method rvc lmerMod
rvc.lmerMod <- function(model, level = 1, delete = NULL, ...) {
  if (hasArg(group)) {
    group <- NULL
    warning("group is not a valid argument for this function. As of version 0.4.0, group has been replaced by level.")
  }
  
  delete <- case_delete(model, level = level, type = "varcomp", delete = delete)
  return( rvc(delete) )
}

#' @export
#' @rdname rvc.mer
#' @method rvc lme
rvc.lme <- function(model, level = 1, delete = NULL, ...) {
  if (hasArg(group)) {
    group <- NULL
    warning("group is not a valid argument for this function. As of version 0.4.0, group has been replaced by level.")
  }
  delete <- case_delete(model, level = level, type = "varcomp", delete = delete)
  return( rvc(delete) )
}



#' @export
#' @rdname leverage.mer
#' @method leverage lme
leverage.lme <- function(model, level = 1, ...) {
  if(!isNestedModel(model)) {
    stop("leverage.mer has not yet been implemented for models with 
         crossed random effects")
  }
  if(!level %in% c( 1, names(model$groups) )) {
    stop("level can only be 1 or a grouping factor from the fitted model.")
  }
  
  mats <- .lme_matrices(model)
  
  X <- mats$X
  # Z <- BlockZ(model)
  
  n     <- nrow(X)
  #  nt    <- model@dims[["nt"]]  # number of random-effects terms in the model
  ngrps <- unname( summary(model)$ngrps )
  
  D <- mats$D
  Z <- mats$Z
  ZDZt <- model$sigma^2 * Z %*% D %*% t(Z)
  
  Vinv   <- mats$Vinv
  
  xvix.inv <- mats$XVXinv
  
  H1 <- X %*% xvix.inv %*% t(X) %*% Vinv
  H2 <- ZDZt %*% Vinv %*% (Diagonal( n = n ) - H1)
  
  diag.H1 <- diag(H1)
  diag.H2 <- diag(H2)
  diag.H2.uc <- diag(ZDZt) / model$sigma^2
  
  if(level == 1) {
    lev1 <- data.frame(overall = diag.H1 + diag.H2, fixef = diag.H1, 
                       ranef =  diag.H2, ranef.uc = diag.H2.uc)
    #     class(lev1) <- "leverage"
  } else {
    flist   <- data.frame( model$groups[[level]] )
    
    grp.lev.fixef <- aggregate(diag.H1, flist, mean)[,2]
    grp.lev.ranef <- aggregate(diag.H2, flist, mean)[,2]
    grp.lev.ranef.uc <- aggregate(diag.H2.uc, flist, mean)[,2]
    
    grp.lev <- data.frame( overall = grp.lev.fixef + grp.lev.ranef,
                           fixef = grp.lev.fixef, 
                           ranef = grp.lev.ranef,
                           ranef.uc = grp.lev.ranef.uc)
    #     class(grp.lev) <- "leverage"
  }
  
  if(level == 1) return(lev1)
  if(level != 1) return(grp.lev)
}
