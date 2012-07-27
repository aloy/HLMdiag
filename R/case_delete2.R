#' Implementing Case Deletion.
#' 
#' @param model the original hierarchical model fit using \code{lmer()}
#' @param group a variable used to define the group for which cases will be deleted.
#'   If this is left \code{FALSE}, then the function will delete individual observations.
#' @param type the part of the model for which you are obtaining deletion diagnostics:
#'   the fixed effects (\code{fixef}), random effects (\code{ranef}), 
#'   variance components (\code{varcomp}), or \code{all}.
#'@return a list with the following compontents:
#'@return fixef.original the original fixed effects
#'@return ranef.original the origingal random effects
#'@return vcov.original the original variance-covariance parameters
#'@return varcomp.original the original variance components
#'@return fixef.delete a list of the fixed effects obtained through case
#'deletion
#'@return ranef.delete a list of the changes in random effects obtained through case
#'deletion
#'@return vcov.delete a list of the variance-covariance parameters obtained
#'through case deletion -- currently NULL
#'@return fitted.delete a list of the fitted values obtained through case
#'deletion -- Currently NULL
#'@return varcomp.delete a list of the variance components obtained through
#'case deletion -- Currently NULL
#' @example
#' library(mlmRev)
#' fm <- lmer(normexam ~ standLRT + I(standLRT^2) + I(standLRT^3) + schgend + (schgend | school), data = Exam)
#' fm.del <- case_delete2(fm)
case_delete2 <- function(model, group = FALSE, type = c("both", "fixef", "varcomp")){  
  
  # Extract key pieces of the model
  Y <- model@y
  X <- getME(model, "X")
  betaHat <- matrix(fixef(model), ncol = 1)
  
  n <- length(Y)
  ngrps <- summary(model)@ngrps
  
  flist <- model@flist
  
  # Constructing V = Cov(Y)
  fmVC <- VarCorr(model)
  sig0 <- unname(attr(fmVC, "sc"))

  D <- kronecker(Diagonal(ngrps), fmVC[[1]]) / sig0^2
  
  A <- model@A
  V <- (Diagonal(n) + crossprod(model@A))
  
  # Inverting V
  Vinv <- (Diagonal(n) - t(A) %*% solve(Diagonal(nrow(A)) + tcrossprod(A)) %*% A)
  
  XVXinv <- solve(crossprod(getME(model, "RX"))) 
  
  P <- Vinv - (Vinv %*% X %*% XVXinv %*% t(X) %*% Vinv)
  
  if(group == FALSE){
    # Obtain building blocks -- individual case deletion
    fixef.delete <- .Call("bbFixef1", y_ = Y, X_ = as.matrix(X), Vinv_ = as.matrix(Vinv), 
                       XVXinv_ = as.matrix(XVXinv), beta_ = betaHat, PACKAGE = "HLMdiag")
  }
  
  if(group == TRUE){
    # Obtain building blocks -- group deletion
    groups <- unlist(unique(flist))
    groups.index <- list()
    for(i in groups){
      groups.index[[i]] <- which(flist == i) - 1
    }
    
    fixef.delete <- .Call("bbFixefGroup", groupIndex = groups.index, X_ = X, 
                       Vinv_ = as.matrix(Vinv), XVinv_ = as.matrix(XVXinv),
                       P_ = as.matrix(P), e_ = as.matrix(P %*% Y),
                       PACKAGE = "HLMdiag")
  }
  
  ranef.delete <- .Call("bbRanef", Zt_ = as.matrix(getME(model, "Zt")), 
                        D_ = as.matrix(D), P_ = as.matrix(P), 
                        e_ = as.numeric(P %*% Y), PACKAGE = "HLMdiag")
  
  sigma.original <- attr(VarCorr(model), "sc")
  varcomp.original <- c(sigma2 = sigma.original^2, diag(VarCorr(model)[[names(model@flist)]]))
  
  reslist <- list(fixef.original = betaHat, 
                  ranef.original = ranef(model)[[names(model@flist)]],
                  vcov.original = as.matrix(vcov(model)), 
                  varcomp.original = varcomp.original,
                  fixef.delete = fixef.delete, 
                  ranef.delete = ranef.delete,
                  vcov.delete = NULL, 
                  fitted.delete = NULL,
                  varcomp.delete = NULL)
  
  return(reslist)
}

#------------------------------------------------------------------------------
#' Cook's distance for fixed effects
#'
#' @example
#' library(mlmRev)
#' fm <- lmer(normexam ~ standLRT + I(standLRT^2) + I(standLRT^3) + schgend + (schgend | school), data = Exam)
#' fm.del <- case_delete2(fm)
#' fm.cd  <- cook_fixef(fm.del$fixef.delete, fixef(fm), sqrt(fm.del$varcomp.original[1]))
cook_fixef <- function(bb, fixef, sig0) {
  sapply(bb, 
         function(x){
           as.numeric(t(x$ybb - t(x$xbb) %*% fixef) %*% x$hbb %*% 
             (x$ybb - t(x$xbb) %*% fixef) / (sig0^2 * length(x$beta_cdd) * (x$sbb - x$hbb)^2))
                     })
}

