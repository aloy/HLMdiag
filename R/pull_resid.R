#' @export
pull_resid <- function(object, ...){
  UseMethod("pull_resid", object)
}

#' @export
#' @rdname pull_resid.lmerMod
#' @method pull_resid default
pull_resid.default <- function(object, ...){
  stop(paste("there is no pull_resid() method for objects of class",
             paste(class(object), collapse=", ")))
}

#' Computationally Efficient HLM Residuals
#'
#' \code{pull_resid} takes a hierarchical linear model fit as a \code{lmerMod}
#' or \code{lme} object and returns various types of level-1 residuals as a
#' vector. Because the \code{pull_resid} only calculates one type of residual,
#' it is more efficient than using \code{\link{hlm_resid}} and indexing the
#' resulting tibble. \code{pull_resid} is designed to be used with methods that
#' take a long time to run, such as the resampling methods found in the
#' \code{lmeresampler} package.
#' @param object an object of class \code{lmerMod} or \code{lme}.
#' @param type which residuals should be returned. Can be either 'ls', 'eb', or
#'   'marginal'
#' @param standardize a logical indicating if residuals should be standardized
#' @param ...  not in use
#' @details \describe{
#' \item{\code{type = "ls"}}{Residuals calculated by fitting separate LS
#' regression models for each group. LS residuals are unconfounded by higher
#' level residuals, but unreliable for small within-group sample sizes. When
#' \code{standardize = TRUE}, residuals are standardized by sigma components of
#' the model object.}
#' \item{\code{type = "eb"}}{Residuals calculated using the empirical Bayes (EB)
#' method using maximum likelihood. EB residuals are interrelated with higher
#' level residuals. When \code{standardize = TRUE}, residuals are standardized
#' by sigma components of the model object.}
#' \item{\code{type = "marginal"}}{Marginal residuals only consider the fixed
#' effect portion of the estimates. When \code{standardize = TRUE}, Cholesky
#' residuals are returned.}
#' }
#' @seealso \link[HLMdiag]{hlm_resid}
#'
#' @export
#' @method pull_resid lmerMod
#' @aliases pull_resid

pull_resid.lmerMod <- function(object, type = "ls", standardize = FALSE, ...) {
  
  if(!is.null(standardize) && !standardize %in% c(TRUE, FALSE)) {
    stop("standardize can only be specified to be TRUE or FALSE.")
  }
  if(!type %in% c("ls", "eb", "marginal")) {
    stop("type must be either 'ls', 'eb', or 'marginal'.")
  }
  
  if(type == "ls") {
    ls.resid <- LSresids(object, level = 1, standardize = standardize)
    ls.resid <- ls.resid[order(as.numeric(rownames(ls.resid))),]
    
    return(ls.resid[,1])
  }
  
  if(type == "eb") {
    if (standardize == TRUE) {
      eb.resid <- data.frame(.std.resid = resid(object, scale = TRUE))
    } else {
      eb.resid <- data.frame(.resid = resid(object))
    }
    return(eb.resid[,1])
  }
  
  if(type == "marginal") {
    mr <- object@resp$y - lme4::getME(object, "X") %*% lme4::fixef(object)
    if (standardize == TRUE) {
      sig0 <- lme4::getME(object, "sigma")
      ZDZt <- sig0^2 * crossprod( lme4::getME(object, "A") )
      n    <- nrow(ZDZt)
      R      <- Diagonal( n = n, x = sig0^2 )
      V      <- R + ZDZt
      V.chol <- chol( V )
      
      Lt <- solve(t(V.chol))
      mar.resid <- data.frame(.chol.mar.resid = (Lt %*% mr)[,1])
      
    } else {
      mar.resid <- data.frame(.mar.resid = mr[,1])
    }
    return(mar.resid[,1])
  }
}


#' @export
#' @rdname pull_resid.lmerMod
#' @method pull_resid lme
pull_resid.lme <- function(object, type = "ls", standardize = FALSE, ...) {
  
  if(!is.null(standardize) && !standardize %in% c(TRUE, FALSE)) {
    stop("standardize can only be specified to be TRUE or FALSE.")
  }
  if(!type %in% c("ls", "eb", "marginal")) {
    stop("type must be either 'ls', 'eb', or 'marginal'.")
  }
  
  if(type == "ls") {
    ls.resid <- LSresids(object, level = 1, standardize = standardize)
    ls.resid <- ls.resid[order(as.numeric(rownames(ls.resid))),]
    
    return(ls.resid[,1])
  }
  
  if(type == "eb") {
    if(standardize == TRUE) {
      eb.resid <- data.frame(.std.resid = resid(object, type = "normalized"))
    } else { 
      eb.resid <- data.frame(.resid = resid(object, type = "response"))
    }
    return(eb.resid[,1])
  }
  
  if(type == "marginal") {
    mr <- resid(object, type="response", level=0)
    if (standardize == TRUE) {
      V      <- extract_design(object)$V
      V.chol <- chol( V )
      
      Lt <- solve(t(V.chol))
      mar.resid <- data.frame(.chol.mar.resid = (Lt %*% mr)[,1])
      
    } else {
      mar.resid <- data.frame(.mar.resid = mr)
    }
    return(mar.resid[,1])
  }
}
