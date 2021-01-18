#' @export
resid_marginal <- function(object, ...){
  UseMethod("resid_marginal", object)
}

#' @export
#' @rdname resid_marginal
#' @method resid_marginal default
resid_marginal.default <- function(object, ...){
  stop(paste("there is no resid_marginal() method for objects of class",
             paste(class(object), collapse=", ")))
}


#' @export
resid_conditional <- function(object, ...){
  UseMethod("resid_conditional", object)
}

#' @export
#' @rdname resid_conditional
#' @method resid_conditional default
resid_conditional.default <- function(object, ...){
  stop(paste("there is no resid_conditional() method for objects of class",
             paste(class(object), collapse=", ")))
}


#' @export
#' @method resid_marginal lmerMod
#' @aliases resid_marginal
resid_marginal.lmerMod <- function(object, standardize =  FALSE, ...){
  if(!is.logical(standardize)) {
    stop("standardize must be logical (TRUE or FALSE).")
  }
  
  res <- lme4::getME(object, "y") - stats::predict(object, re.form = NA)
  
  if(standardize == TRUE){
    res_names <- names(res)
    mats <- .lmerMod_matrices(object)
    Lt <- solve(t(mats$V.chol))
    res <- as.numeric(Lt %*% res)
    names(res) <- res_names
  } 
  
  res
}

#' @export
#' @method resid_marginal lme
#' @aliases resid_marginal
resid_marginal.lme <- function(object, standardize =  FALSE, ...){
  if(!is.logical(standardize)) {
    stop("standardize must be logical (TRUE or FALSE).")
  }
  
  res <- resid(object, type = "response", level = 0)
  
  if(standardize == TRUE){
    res_names <- names(res)
    V      <- extract_design(object)$V
    V.chol <- chol(V)
    Lt <- solve(t(V.chol))
    res <- as.numeric(Lt %*% res)
    names(res) <- res_names
  } 
  
  res
}


#' @export
#' @method resid_conditional lmerMod
#' @aliases resid_conditional
resid_conditional.lmerMod <- function(object, standardize = FALSE, ...){
  if(!is.logical(standardize)) {
    stop("standardize must be logical (TRUE or FALSE).")
  }
  
  if(standardize){
    resid(object, scaled = TRUE)
  } else {
    resid(object)
  }
}

#' @export
#' @method resid_conditional lme
#' @aliases resid_conditional
resid_conditional.lme <- function(object, standardize = FALSE, ...){
  if(!is.logical(standardize)) {
    stop("standardize must be logical (TRUE or FALSE).")
  }
  
  if(standardize){
    resid(object, type = "normalized")
  } else {
    resid(object)
  }
}

#' @export
resid_ranef <- function(object, ...){
  UseMethod("resid_ranef", object)
}


#' @export
#' @method resid_ranef lmerMod
#' @aliases resid_ranef
resid_ranef.lmerMod<- function(object, level, which, standardize = FALSE, ...){
  if(!is.logical(standardize)) {
    stop("standardize must be logical (TRUE or FALSE).")
  }
  
  # Allow level to be numeric or character...
  
  # Allow which to follow lme example...
  
  flist <- getME(object, "flist")

  if(standardize){
    re <- ranef(object, condVar = TRUE)
    vc <- VarCorr(object)
    for(i in names(flist)) {
      diag_var <- diag(as.matrix(vc[[i]])) - diag(as.matrix(attr(re[[i]], "postVar")[,,1]))
      re[[i]] <- sweep(re[[i]], 2, sqrt(diag_var), FUN = "/")
    }
  } else {
    re <- ranef(object)
  }
  
  re
}


#' @export
#' @import diagonals
mahalanobis_ranef.lmerMod <- function(object, ...){
  mats <- HLMdiag:::.lmerMod_matrices(object)
  
  n_lev <- length(getME(object, "flist"))
  
  if(n_lev == 1) {
    Z <- getME(object, "Z")
    vc <- VarCorr(object)
    D  <- kronecker(Diagonal(mats$ngrps), bdiag(vc))
    
    eblup <- tcrossprod(D, Z) %*% mats$Vinv %*% resid_marginal(object)
    # vcov_eblup <- D - tcrossprod(D, Z) %*% mats$P %*% tcrossprod(Z, D)
    vcov_eblup <- tcrossprod(D, Z) %*% mats$P %*% tcrossprod(Z, D)
    
    eblup_lst      <- diagonals::split_vector(eblup, size = 2)
    vcov_eblup_lst <- diagonals::fatdiag(vcov_eblup, steps = ngrps(object)) %>%
      diagonals::split_vector(size = 4) %>%
      map(~matrix(.x, nrow = 2, byrow = TRUE))
    
    mah_dist_eblup <- map2_dbl(eblup_lst, vcov_eblup_lst, ~t(.x) %*% MASS::ginv(.y) %*% .x)
  } else{
    ### Need to check for higher-level models.... D will fail...
  } 
  }