# Reorganizing Z matrix
# This function reorganizes the Z matrix obtained from an \code{mer} object
# into the block form discussed by Demidenko (2004). Currently, this function
# assumes there is only one level and that units are nested.
BlockZ <- function(object) {
  Z <- getME(object, "Z")
  
  grp.size <- table(object@flist)
  ngrps <- length(grp.size)
  nranef <- dim(ranef(object)[[1]])[2]
  
  base.ord <- seq(from = 1, by = ngrps, length.out = nranef)
  ord <- base.ord + rep(0:(ngrps - 1), each = nranef)
  
  perm.mat <- t(as(ord, "pMatrix"))
  
  return(Z %*% perm.mat)
}

#' Extracting variance components
#' 
#' This function extracts the variance components from a mixed/hierarchical
#' linear model fit using \code{lmer}. 
#' 
#' @return A named vector is returned. \code{sigma2} denotes the residual
#' variance. The other variance components are names \code{D**} where the
#' trailing digits specify the of that variance component in the covariance
#' matrix of the random effects.
#' 
#' @param object a fitted model object of class \code{mer} or \code{lmerMod}.
#' @author Adam Loy \email{loyad01@@gmail.com}
#' @keywords models regression
#' @export
#' @examples
#' data(sleepstudy, package = "lme4") 
#' fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
#' varcomp.mer(fm1)
varcomp.mer <- function(object) {
  vc  <- VarCorr(object)
  sig <- attr(vc, "sc")
  vc.mat <- bdiag(vc)
  
  if(isDiagonal(vc.mat)) {
    vc.vec   <- diag(vc.mat)
    vc.names <- paste("D", 1:length(vc.vec), 1:length(vc.vec), sep="")
  } else{
    vc.vec <- as.matrix(vc.mat)[!upper.tri(vc.mat)]
    vc.index <- which(!upper.tri(vc.mat) == TRUE, arr.ind = TRUE)
    vc.names <- paste("D", vc.index[,1], vc.index[,2], sep="")
  }

  res <- c(sig^2, vc.vec)
  names(res) <- c("sigma2", vc.names)
  return(res)
}

# Checking if matrix is diagonal 
# @param mat a matrix
isDiagonal <- function(mat, tol = 1e-10) {
  if( !isSymmetric(mat) ) return( FALSE )
  else { 
    diag(mat) <- 0
    return(all(abs(mat) < tol))
  }
}

# Extracting/calculating key matrices from mer object 
# @param model an mer object
.mer_matrices <- function(model) {
  Y <- model@y
  X <- getME(model, "X")
  
  n <- length(Y)
  
  flist <- getME(model, "flist")
  ngrps <- sapply(flist, function(x) length(levels(x)))
  
  # Constructing V = Cov(Y)
  sig0 <- sigma(model)
  
  ZDZt <- sig0^2 * crossprod( getME(model, "A") )
  R    <- Diagonal( n = n, x = sig0^2 )
  V    <- Diagonal(n) + ZDZt
  
  # Inverting V
  V.chol <- chol( V )
  Vinv   <- chol2inv( V.chol )
  
  # Calculating P
  XVXinv <- solve( t(X) %*% Vinv %*% X )
  VinvX  <- Vinv %*% X
  M      <- VinvX %*% XVXinv %*% t(VinvX)
  P      <- .Call("cxxmatsub", as.matrix(Vinv), as.matrix(M), 
                  PACKAGE = "HLMdiag")
  
  return( list(Y = Y, X = X, n = n, ngrps = ngrps, flist = flist,
               sig0 = sig0, V = V, Vinv = Vinv, XVXinv = XVXinv,
               M = M, P = P) )
}

# Extracting/calculating key matrices from mer object 
# @param model an lmerMod object
.lmerMod_matrices <- function(model) {
  Y <- model@resp$y
  X <- getME(model, "X")
  
  n <- length(Y)
  
  flist <- getME(model, "flist")
  ngrps <- sapply(flist, function(x) length(levels(x)))
  
  # Constructing V = Cov(Y)
  sig0 <- sigma(model)
  
  ZDZt <- sig0^2 * crossprod( getME(model, "A") )
  R    <- Diagonal( n = n, x = sig0^2 )
  V    <- Diagonal(n) + ZDZt
  
  # Inverting V
  V.chol <- chol( V )
  Vinv   <- chol2inv( V.chol )
  
  # Calculating P
  XVXinv <- solve( t(X) %*% Vinv %*% X )
  VinvX  <- Vinv %*% X
  M      <- VinvX %*% XVXinv %*% t(VinvX)
  P      <- .Call("cxxmatsub", as.matrix(Vinv), as.matrix(M), 
                  PACKAGE = "HLMdiag")
  
  return( list(Y = Y, X = X, n = n, ngrps = ngrps, flist = flist,
               sig0 = sig0, V = V, Vinv = Vinv, XVXinv = XVXinv,
               M = M, P = P) )
  
}

# 'se.ranef' is a copy of function in arm package. This is copied to ensure
# that is available to all users. This should not be exported.
se.ranef <- function (object) 
{
  se.bygroup <- ranef(object, postVar = TRUE)
  n.groupings <- length(se.bygroup)
  for (m in 1:n.groupings) {
    vars.m <- attr(se.bygroup[[m]], "postVar")
    K <- dim(vars.m)[1]
    J <- dim(vars.m)[3]
    se.bygroup[[m]] <- array(NA, c(J, K))
    for (j in 1:J) {
      se.bygroup[[m]][j, ] <- sqrt(diag(as.matrix(vars.m[, 
                                                         , j])))
    }
    names.full <- dimnames(se.bygroup)
    dimnames(se.bygroup[[m]]) <- list(names.full[[1]], names.full[[2]])
  }
  return(se.bygroup)
}

# Checking whether an LMM is nested
isNestedModel <- function(object) {
  fl   <- object@flist
  fnms <- names(fl) 
  RVAL <- all(sapply(seq_along(fl)[-1], function(i) lme4::isNested(fl[[i-1]], fl[[i]])))
  
  return(RVAL)
}