#' @export
rotate_ranef <- function(.mod, ...){
  UseMethod("rotate_ranef", .mod)
}

#' @export
#' @rdname rotate_ranef.mer
#' @method rotate_ranef default
rotate_ranef.default <- function(.mod, ...){
  stop(paste("there is no rotate_ranef() method for objects of class",
             paste(class(.mod), collapse=", ")))
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
#' @param .mod an object of class \code{mer} or \code{lmerMod}.
#' @param .L a matrix defining which combination of random effects are of interest.
#' @param s the dimension of the subspace of interest.
#' @param .varimax if \code{.varimax = TRUE} than the raw varimax rotation 
#'   will be applied to the resulting rotation.
#' @param ... do not use
#' @author Adam Loy \email{loyad01@@gmail.com}
rotate_ranef.mer <- function(.mod, .L, s = NULL, .varimax = FALSE, ...) {
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
  P      <- .Call("cxxmatsub", BB = as.matrix(Vinv), CC = as.matrix(M), 
                  PACKAGE = "HLMdiag")
  
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
rotate_ranef.lmerMod <- function(.mod, .L, s = NULL, .varimax = FALSE, ...) {
  y <- .mod@resp$y
  X <- getME(.mod, "X")
  Z <- getME(.mod, "Z")
  
  n <- nrow(X)
  p <- ncol(X)
  ngrps <- unname( summary(.mod)$ngrps )
  
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
  P      <- .Call("cxxmatsub", BB = as.matrix(Vinv), CC = as.matrix(M), 
                  PACKAGE = "HLMdiag")
  
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


# #' @export
# #' @rdname rotate_ranef.mer
# #' @method rotate_ranef lme
# #' @S3method rotate_ranef lme
# rotate_ranef.lme <- function(.mod, .L, s = NULL, .varimax = FALSE, ...) {
  # design.info <- extract.lmeDesign(.mod)
  
  # y <- design.info$y
  # X <- design.info$X
  # Z <- Matrix( design.info$Z )
  # D <- Matrix( design.info$Vr )
  
  # V  <- .extractV.lme( .mod )
  # V.chol <- chol( V )
  # Vinv  <- chol2inv( V.chol ) 
  
  # XVXinv <- solve( t(X) %*% Vinv %*% X )
  # VinvX  <- Vinv %*% X
  # M      <- VinvX %*% XVXinv %*% t(VinvX)
  # P      <- .Call("cxxmatsub", BB = as.matrix(Vinv), CC = as.matrix(M), 
                  # PACKAGE = "HLMdiag")
  
  # betahat <- solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv %*% y
  # mr <- y - X %*% betahat
  
  # bvec <- D %*% t(Z) %*% Vinv %*% mr
  
  # pzdl <- P %*% Z %*% D %*% .L
  # A <- crossprod( pzdl )
  # B <- t(.L) %*% D %*% t(Z) %*% P %*% Z %*% D %*% .L ## diagnostic se
  # W <- try( mcrotate(A, B, s) )
  # if( class(W) == "try-error") {W <- NA} else {W <- as.matrix(W)}
  
  # if( .varimax == TRUE) {
    # W <- try( varimax(W, normalize = FALSE)$loadings )
    # if( class(W) == "try-error" ) W <- NA 
  # }
  
  # return( as.numeric( t(W) %*% as.numeric( t(.L) %*% bvec ) ) )
# }