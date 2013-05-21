#' @export
rotate_ranef <- function(object, ...){
  UseMethod("rotate_ranef", object)
}

#' @export
#' @rdname rotate_ranef.mer
#' @method rotate_ranef default
rotate_ranef.default <- function(object, ...){
  stop(paste("there is no rotate_ranef() method for objects of class",
             paste(class(object), collapse=", ")))
}

#' Calculate s-dimensional rotated random effects
#' 
#' This function calculates reduced dimensional rotated random effects. 
#' The rotation reduces the influence of the residuals from other levels
#' of the model so that distributional assessment of the resulting
#' random effects is possible.
#' 
#' @export
#' @method rotate_ranef mer
#' @S3method rotate_ranef mer
#' @aliases rotate_ranef
#' @param .mod
#' @param .L
#' @param s
#' @param .varimax
rotate_ranef.mer <- function(.mod, .L, s = NULL, .varimax = FALSE) {
  y <- .mod@y
  X <- getME(.mod, "X")
  Z <- BlockZ(.mod)
  
  n <- nrow(X)
  p <- ncol(X)
  ngrps <- unname( summary(.mod)@ngrps )
  
  vc <- VarCorr(.mod)
  Di <- bdiag( VarCorr(.mod) ) / (unname(attr(vc, "sc")))^2
  D  <- kronecker( Diagonal(ngrps), Di )
  
  Aslot <- .mod@A # ZDZ'
  zdzt <- crossprod( .mod@A )
  V  <- Diagonal( n ) + zdzt
  V.chol <- chol( V )
  Vinv  <- chol2inv( V.chol ) 
  
  XVXinv <- solve( t(X) %*% Vinv %*% X )
  VinvX  <- Vinv %*% X
  M      <- VinvX %*% XVXinv %*% t(VinvX)
  P      <- cxxmatsub(as.matrix(Vinv), as.matrix(M))
  
  betahat <- solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv %*% y
  mr <- y - X %*% betahat
  
  bvec <- D %*% t(Z) %*% Vinv %*% mr
    
  pzdl <- P %*% Z %*% D %*% .L
  A <- crossprod( pzdl )
  B <- t(.L) %*% D %*% t(Z) %*% P %*% Z %*% D %*% .L ## diagnostic se
  W <- try( mcrotate(A, B, s) )
  if( class(W) == "try-error") {W <- NA} else {W <- as.matrix(W)}
    
  if( .varimax == TRUE) {
    W <- try( varimax(W, normalize = FALSE)$loadings )
    if( class(W) == "try-error" ) W <- NA 
  }
    
  return( as.numeric( t(W) %*% as.numeric( t(.L) %*% bvec ) ) )
}


#' @export
#' @rdname rotate_ranef.mer
#' @method rotate_ranef lmerMod
#' @S3method rotate_ranef lmerMod
rotate_ranef.lmerMod <- function(.mod, .L, s = NULL, .varimax = FALSE) {
  y <- .mod@resp$y
  X <- getME(.mod, "X")
  Z <- getME(.mod, "Z")
  
  n <- nrow(X)
  p <- ncol(X)
  ngrps <- unname( summary(.mod)@ngrps )
  
  vc <- VarCorr(.mod)
  Di <- bdiag( VarCorr(.mod) ) / (unname(attr(vc, "sc")))^2
  D  <- kronecker( Diagonal(ngrps), Di )
  
  zdzt <- crossprod( getME(.mod, "A") )
  V  <- Diagonal( n ) + zdzt
  V.chol <- chol( V )
  Vinv  <- chol2inv( V.chol ) 
  
  XVXinv <- solve( t(X) %*% Vinv %*% X )
  VinvX  <- Vinv %*% X
  M      <- VinvX %*% XVXinv %*% t(VinvX)
  P      <- cxxmatsub(as.matrix(Vinv), as.matrix(M))
  
  betahat <- solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv %*% y
  mr <- y - X %*% betahat
  
  bvec <- D %*% t(Z) %*% Vinv %*% mr
  
  pzdl <- P %*% Z %*% D %*% .L
  A <- crossprod( pzdl )
  B <- t(.L) %*% D %*% t(Z) %*% P %*% Z %*% D %*% .L ## diagnostic se
  W <- try( mcrotate(A, B, s) )
  if( class(W) == "try-error") {W <- NA} else {W <- as.matrix(W)}
  
  if( .varimax == TRUE) {
    W <- try( varimax(W, normalize = FALSE)$loadings )
    if( class(W) == "try-error" ) W <- NA 
  }
  
  return( as.numeric( t(W) %*% as.numeric( t(.L) %*% bvec ) ) )
}


mcrotate <- function(A, B, s) {
  r <- rankMatrix(B)
  if(is.null(s)) s <- r
  
  B.svd <- svd(B)
  Cr.diag <- B.svd$d[1:r]
  Tr <- B.svd$u[, 1:r]
  
  A.star <- Diagonal( x = 1 / sqrt(Cr.diag) ) %*% t(Tr) %*% 
    A %*% 
    Tr %*% Diagonal( x = 1 / sqrt(Cr.diag) )
  
  A.star.svd <- svd( A.star )
  
  index <- seq(r, length.out = s, by = -1)
  index <- sort(index[index >= 0])
  W <- Tr %*% Diagonal( x = 1 / sqrt( Cr.diag ) ) %*% A.star.svd$u[,index]
  
  return(W)
}