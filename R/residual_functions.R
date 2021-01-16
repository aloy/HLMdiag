#' @export
resid_marginal.lmerMod <- function(object, standardize, ...){
  res <- lme4::getME(object, "y") - stats::predict(object, re.form = NA)
  
  if(standardize == TRUE){
    mats <- .lmerMod_matrices(object)
    Lt <- solve(t(mats$V.chol))
    res <- as.numeric(Lt %*% res)
  } 
  
  res
}


#' @export
resid_conditional.lmerMod <- function(object, standardize, ...){
  if(standardize){
    res <- resid(object, scaled = TRUE)
  } else {
    resid(object)
  }
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
    vcov_eblup <- D - tcrossprod(D, Z) %*% mats$P %*% tcrossprod(Z, D)
    # vcov_eblup <- tcrossprod(D, Z) %*% mats$P %*% tcrossprod(Z, D)
    
    eblup_lst      <- diagonals::split_vector(eblup, size = 2)
    vcov_eblup_lst <- diagonals::fatdiag(vcov_eblup, steps = ngrps(object)) %>%
      diagonals::split_vector(size = 4) %>%
      map(~matrix(.x, nrow = 2, byrow = TRUE))
    
    mah_dist_eblup <- map2_dbl(eblup_lst, vcov_eblup_lst, ~t(.x) %*% MASS::ginv(.y) %*% .x)
  } else{
    ### Need to check for higher-level models.... D will fail...
  } 
  }