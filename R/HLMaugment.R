#' @export
HLMaugment <- function(object, ...){
  UseMethod("HLMaugment", object)
}


#' @export
#' @rdname HLMaugment.lmerMod
#' @method HLMaugment default
#' @S3method HLMaugment default
HLMaugment.default <- function(object, ...){
  stop(paste("there is no HLMaugment() method for objects of class",
             paste(class(object), collapse=", ")))
}


#' Calculating residuals from HLMs
#'
#' \code{HLMaugment} takes a hierarchical linear model fit as a
#' \code{lmerMod} object and adds information about each observation's 
#' residuals and predicted values.
#' 
#' This function extract residuals and predicted values from the model, using
#' least squares (LS) and Empirical Bayes (EB) methods, and appends them to the
#' model data. This unified framework enables the analyst to more easily conduct
#' an upward residual analysis during model exploration/checking.
#'
#' @export
#' @method HLMaugment lmerMod
#' @S3method HLMaugment lmerMod
#' @aliases HLMaugment
#' @param object an object of class \code{lmerMod}.
#' @param level which residuals should be extracted: 1 for within-group
#'   (case-level) residuals, the name of a grouping factor (as defined in
#'   \code{flist} of the \code{lmerMod} object) for between-group residuals
#' @param standardize if \code{standardize = TRUE} the standardized residuals
#'   will be returned; if \code{standardize = "semi"} then the semi-standardized
#'   level-1 residuals will be returned
#' @param sim optional argument giving the data frame used for LS residuals.
#'   This is used mainly for dealing with simulations.
#' @param ... do not use
#' @details The \code{HLMaugment} function provides a wrapper that will extract
#' residuals and predicted values from a fitted \code{lmerMod} object. 
#' The function provides access to 
#' residual quantities already made available by the functions \code{resid},
#' \code{predict}, and \code{ranef}, but adds additional functionality. Below is
#' a list of types of residuals and predicted values that are extracted and
#' appended to the model data.
#' \describe{
#' \item{raw level-1 LS residuals}{These are equivalent to the residuals extracted
#' by \code{resid} if \code{level = 1}, \code{type = "EB"}, and 
#' \code{standardize = FALSE} is specified. }
#' \item{level-1 LS fitted values}{The predicted values }
#' }
#' Note that \code{standardize = "semi"} is only implemented for level-1 LS residuals.
HLMaugment.lmerMod <- function(object, level = 1, standardize = FALSE, sim = NULL, ...) {
  # LS Residuals
  ls.resid <- LSresids(object, level = 1, stand = standardize, sim = sim)
  ls.resid <- ls.resid[order(as.numeric(rownames(ls.resid))),]
  
  if (standardize == FALSE) {
    ls.resid <- data.frame(ls.resid["LS.resid"], ls.resid["fitted"])
    names(ls.resid) <- c("LS.resid", "LS.fitted")
    
  } else if (standardize == TRUE) {
    ls.resid <- data.frame(ls.resid["std.resid"], ls.resid["fitted"])
    names(ls.resid) <- c("std.LS.resid", "LS.fitted")
    
  } else {
    ls.resid <- data.frame(ls.resid["semi.std.resid"], ls.resid["fitted"])
    names(ls.resid) <- c("semi.LS.resid", "LS.fitted")
  }
  
  # we should refine what LSresids returns to match EB method
  # I am unsure if the above code works when a sim argument is passed in
  
  # EB Residuals
  if (standardize == TRUE) {
    mats <- .lmerMod_matrices(object)
    p_diag <- diag(mats$P)
    EB.resid <- data.frame(std.EB.resid = 
                             resid(object) / ( lme4::getME(object, "sigma") * sqrt(p_diag) ))
    
  } else {
    EB.resid <- data.frame(EB.resid = resid(object))
  }
  # note that "semi" is not implemented for EB
  
  # Fitted Values
  Xbeta  <- data.frame(Xbeta = predict(object, re.form = ~0))
  XbetaZb <- data.frame(XbetaZb = lme4::getME(object, "mu"))
  
  # Marginal Residuals
  mr <- object@resp$y - lme4::getME(object, "X") %*% lme4::fixef(object)
  if (standardize == TRUE) {
    sig0 <- lme4::getME(object, "sigma")
    ZDZt <- sig0^2 * crossprod( lme4::getME(object, "A") )
    n    <- nrow(ZDZt)
    
    R      <- Diagonal( n = n, x = sig0^2 )
    V      <- R + ZDZt
    V.chol <- chol( V )
    
    Lt <- solve(t(V.chol))
    mar.resid <- data.frame(std.mar.resid = (Lt %*% mr)[,1])
    
  } else {
    mar.resid <- data.frame(mar.resid = mr[,1])
  }
  
  return.tbl <- tibble::tibble(object@frame,
                               ls.resid,
                               EB.resid,
                               Xbeta,
                               XbetaZb,
                               mar.resid)
  # It might make sense to use tibbles earlier as well, I had issues with
  # renaming the LS columns.
  
  return(return.tbl)
}
