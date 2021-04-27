#' @export
resid_marginal <- function(object, type){
  UseMethod("resid_marginal", object)
}

#' @export
#' @rdname resid_marginal
#' @method resid_marginal default
resid_marginal.default <- function(object, type){
  stop(paste("there is no resid_marginal() method for objects of class",
             paste(class(object), collapse=", ")))
}


#' @export
resid_conditional <- function(object, type){
  UseMethod("resid_conditional", object)
}

#' @export
#' @rdname resid_conditional
#' @method resid_conditional default
resid_conditional.default <- function(object, type){
  stop(paste("there is no resid_conditional() method for objects of class",
             paste(class(object), collapse=", ")))
}


#' @title Marginal residuals
#' 
#' @description 
#' Calculates marginal residuals of \code{lmerMod} and \code{lme} model objects.
#' 
#' @param object an object of class \code{lmerMod} or \code{lme}.
#' @param type a character string specifying what type of residuals should be calculated.
#'   It is set to \code{"raw"} (observed - fitted) by default. Other options include
#'   \code{"pearson"}, \code{"studentized"}, and \code{"cholesky"}. Partial matching of arguments is used, 
#'   so only the first character needs to be provided.
#' 
#' @return
#' A vector of marginal residuals.
#' 
#' @details
#' For a model of the form \eqn{Y = X \beta + Z b + \epsilon},
#' four types of marginal residuals can be calculated:
#' 
#' \describe{
#'   \item{\code{raw}}{\eqn{r = Y - X \hat{beta}}}
#'   \item{\code{pearson}}{\eqn{r / \sqrt{ diag(\hat{Var}(Y)})}}
#'   \item{\code{studentized}}{\eqn{r / \sqrt{ diag(\hat{Var}(r)})}}
#'   \item{\code{cholesky}}{\eqn{\hat{C}^{-1} r} where \eqn{\hat{C}\hat{C}^\prime = \hat{Var}(Y)}}
#' }
#' 
#' @references 
#' Singer, J. M., Rocha, F. M. M., & Nobre, J. S. (2017). 
#' Graphical Tools for Detecting Departures from Linear Mixed Model 
#' Assumptions and Some Remedial Measures. 
#' \emph{International Statistical Review}, \bold{85}, 290--324.
#' 
#' Schabenberger, O. (2004) Mixed Model Influence Diagnostics,
#' in \emph{Proceedings of the Twenty-Ninth SAS Users Group International Conference},
#' SAS Users Group International.
#' 
#' @export
#' @method resid_marginal lmerMod
#' @aliases resid_marginal
#' @rdname resid_marginal
resid_marginal.lmerMod <- function(object, type = c("raw", "pearson", "studentized", "cholesky")){
  type <- match.arg(type)
  
  res <- lme4::getME(object, "y") - stats::predict(object, re.form = NA)
  
  if(type != "raw") {
    res_names <- names(res)
    mats <- .lmerMod_matrices(object)
  }
  if(type == "cholesky"){
    Lt <- solve(t(mats$V.chol))
    res <- as.numeric(Lt %*% res)
    names(res) <- res_names
  } 
  else if(type == "pearson"){
    res <- res / sqrt(diag(mats$V))
  }
  else if(type == "studentized"){
    mar_var <- mats$V - mats$X %*% tcrossprod(mats$XVXinv, mats$X)
    res <- res / sqrt(diag(mar_var))
  }
  
  res
}

#' @export
#' @method resid_marginal lme
#' @aliases resid_marginal
#' @rdname resid_marginal
resid_marginal.lme <- function(object, type = c("raw", "pearson", "studentized", "cholesky")){
  type <- match.arg(type)
  
  res <- resid(object, type = "response", level = 0)
  
  if(type != "raw") {
    res_names <- names(res)
    mats <- .lme_matrices(object)
  }
  if(type == "cholesky"){
    V.chol <- chol(mats$V)
    Lt <- solve(t(V.chol))
    res <- as.numeric(Lt %*% res)
    names(res) <- res_names
  } 
  else if(type == "pearson"){
    res <- res / sqrt(diag(mats$V))
  }
  else if(type == "studentized"){
    mar_var <- mats$V - mats$X %*% tcrossprod(mats$XVXinv, mats$X)
    res <- res / sqrt(diag(mar_var))
  }
  
  res
}


#' @title Conditional residuals
#' 
#' @description 
#' Calculates conditional residuals of \code{lmerMod} and \code{lme} model objects.
#' 
#' @param object an object of class \code{lmerMod} or \code{lme}.
#' @param type a character string specifying what type of residuals should be calculated.
#'   It is set to \code{"raw"} (observed - fitted) by default. Other options include
#'   \code{"pearson"}, \code{"studentized"}, and \code{"cholesky"}. 
#'   Partial matching of arguments is used, so only the first character needs to be provided.
#' 
#' @return
#' A vector of conditional residuals.
#' 
#' @details
#' For a model of the form \eqn{Y = X \beta + Z b + \epsilon},
#' four types of marginal residuals can be calculated:
#' 
#' \describe{
#'   \item{\code{raw}}{\eqn{e = Y - X \hat{beta} - Z \hat{b}}}
#'   \item{\code{pearson}}{\eqn{e / \sqrt{diag(\hat{Var}(Y|b)})}}
#'   \item{\code{studentized}}{\eqn{e / \sqrt{diag(\hat{Var}(e))}}}
#'   \item{\code{cholesky}}{\eqn{\hat{C}^{-1} e} where \eqn{\hat{C}\hat{C}^\prime = \hat{Var}(e)}}
#' }
#' 
#' @references 
#' Singer, J. M., Rocha, F. M. M., & Nobre, J. S. (2017). 
#' Graphical Tools for Detecting Departures from Linear Mixed Model 
#' Assumptions and Some Remedial Measures. 
#' \emph{International Statistical Review}, \bold{85}, 290--324.
#' 
#' Schabenberger, O. (2004) Mixed Model Influence Diagnostics,
#' in \emph{Proceedings of the Twenty-Ninth SAS Users Group International Conference},
#' SAS Users Group International.
#' 
#' @export
#' @method resid_conditional lmerMod
#' @aliases resid_conditional
#' @rdname resid_conditional
resid_conditional.lmerMod <- function(object, type = c("raw", "pearson", "studentized", "cholesky")){
  type <- match.arg(type)
  
  # Raw residuals
  res <- resid(object)
  
  if(type == "cholesky"){
    # For Diagonal R, this will just be the Pearson resids...
    res <- res / sqrt(sigma(object)^2)
  } 
  else if(type == "studentized") {
    mats <- .lmerMod_matrices(object)
    sig0 <- lme4::getME(object, "sigma")
    R <- Diagonal(n = mats$n, x = sig0^2)
    cond_var <- R %*% mats$P %*% R
    res <- res / sqrt(diag(cond_var))
  }
  else if(type == "pearson"){
    res <- resid(object, scaled = TRUE)
  }

  res
}

#' @export
#' @method resid_conditional lme
#' @aliases resid_conditional
#' @rdname resid_conditional
resid_conditional.lme <- function(object, type = c("raw", "pearson", "studentized", "cholesky")){
  type <- match.arg(type)
  
  # Raw residuals
  res <- resid(object)
  
  if(type == "cholesky"){
    # For Diagonal R, this will just be the Pearson resids...
    res <- resid(object, type = "normalized")
  } 
  else if(type == "studentized") {
    mats <- .lme_matrices(object)
    sig0 <- object$sigma
    R <- Diagonal(n = mats$n, x = sig0^2)
    cond_var <- R %*% mats$P %*% R
    res <- res / sqrt(diag(cond_var))
  }
  else if(type == "pearson"){
    res <- resid(object, type = "pearson")
  }
  
  res
}

#' @title Random effects residuals
#' 
#' @description 
#' Calculates Random effects  residuals of \code{lmerMod} model objects.
#' 
#' @param object an object of class \code{lmerMod}.
#' @param level DESCRIPTION
#' @param which DESCRIPTION
#' @param standardize DESCRIPTION
#' 
#' @return
#' A vector of conditional residuals.
resid_ranef <- function(object, level, which, standardize){
  UseMethod("resid_ranef", object)
}


#' @method resid_ranef lmerMod
#' @aliases resid_ranef
resid_ranef.lmerMod <- function(object, level, which, standardize = FALSE){
  if(!is.logical(standardize)) {
    stop("standardize must be logical (TRUE or FALSE).")
  }
  
  # Allow level to be numeric or character...
  
  # Allow which to follow lme example...
  
  flist <- lme4::getME(object, "flist")

  if(standardize){
    re <- lme4::ranef(object, condVar = TRUE)
    vc <- lme4::VarCorr(object)
    for(i in names(flist)) {
      diag_var <- diag(as.matrix(vc[[i]])) - diag(as.matrix(attr(re[[i]], "postVar")[,,1]))
      re[[i]] <- sweep(re[[i]], 2, sqrt(diag_var), FUN = "/")
    }
  } else {
    re <- lme4::ranef(object)
  }
  
  re
}



#' @importFrom diagonals split_vector fatdiag
mahalanobis_ranef.lmerMod <- function(object){
  ngrps <- 0
  mats <- .lmerMod_matrices(object)
  
  n_lev <- length(lme4::getME(object, "flist"))
  
  if(n_lev == 1) {
    Z <- lme4::getME(object, "Z")
    vc <- lme4::VarCorr(object)
    D  <- kronecker(Diagonal(mats$ngrps), bdiag(vc))
    
    eblup <- tcrossprod(D, Z) %*% mats$Vinv %*% resid_marginal(object)
    # vcov_eblup <- D - tcrossprod(D, Z) %*% mats$P %*% tcrossprod(Z, D)
    vcov_eblup <- tcrossprod(D, Z) %*% mats$P %*% tcrossprod(Z, D)
    
    eblup_lst      <- diagonals::split_vector(eblup, size = 2)
    vcov_eblup_lst <- diagonals::fatdiag(vcov_eblup, steps = ngrps(object)) %>%
      diagonals::split_vector(size = 4) %>%
      map(~matrix(.x, nrow = 2, byrow = TRUE))
    
    mah_dist_eblup <- purrr::map2_dbl(eblup_lst, vcov_eblup_lst, ~t(.x) %*% MASS::ginv(.y) %*% .x)
  } else{
    ### Need to check for higher-level models.... D will fail...
  } 
  return(mah_dist_eblup)
  }