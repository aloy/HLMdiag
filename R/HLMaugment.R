#' Calculating residuals from HLMs
#'
#' \code{HLMaugment} is a function that extracts residuals
#' from a hierarchical linear model fit
#' using \code{lmer}. That is, it is a unified framework that
#' extracts/calculates residuals from \code{mer} or \code{lmerMod} objects.
#' 
#' This function extracts residuals from the model, 
#' using least squares (LS) and Empirical 
#' Bayes (EB) methods. This unified framework
#' enables the analyst to more easily conduct
#' an upward residual analysis during model
#' exploration/checking.
#'
#' @export
#' @method HLMresid mer
#' @S3method HLMresid mer
#' @aliases HLMresid
#' @param object an object of class \code{lmerMod}.
#' @param sim optional argument giving the data frame used for LS residuals. This
#'  is used mainly for dealing with simulations.
#' @param standardize if \code{standardize = TRUE} the standardized
#' residuals will be returned; if \code{standardize = "semi"} then
#' the semi-standardized level-1 residuals will be returned. 
#' @param ... do not use
#' @details The \code{HLMaugment} function provides a wrapper that will extract
#' residuals from a fitted \code{mer} or \code{lmerMod} object.  

HLMaugment <- function(object, level = 1, standardize = FALSE, sim = NULL) {
  # LS Residuals
  ls.resid <- LSresids(object, level = 1, stand = standardize, sim = sim)
  
  if (standardize == FALSE) {
    ls.resid <- data.frame(LS.resid = ls.resid["LS.resid"], 
                           LS.fitted = ls.resid["fitted"])
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
  XbetaZb <- data.frame(XbetaZb = getME(object, "mu"))
  
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
  
  return.tbl <- tibble(object@frame,
                       ls.resid,
                       EB.resid,
                       Xbeta,
                       XbetaZb,
                       mar.resid)
  # It might make sense to use tibbles earlier as well, I had issues with
  # renaming the LS columns.
  
  return(return.tbl)
}