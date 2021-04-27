# Reorganizing Z matrix
# This function reorganizes the Z matrix obtained from an \code{mer} object
# into the block form discussed by Demidenko (2004). Currently, this function
# assumes there is only one level and that units are nested.
BlockZ <- function(object) {
  Z <- lme4::getME(object, "Z")
  
  grp.size <- table(object@flist)
  ngrps <- length(grp.size)
  nranef <- dim(lme4::ranef(object)[[1]])[2]
  
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
#' fm1 <- lme4::lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
#' varcomp.mer(fm1)
varcomp.mer <- function(object) {
  vc  <- lme4::VarCorr(object)
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
  
  #remove 0s 
  res <- res[res != 0]
  return(res)
}

varcomp.lme <- function(model) {
  vc <- model$modelStruct$reStruct
  sig2 <- sigma(model)^2
  vc.list <- lapply(vc, function(x) {
    vc <- as.matrix(x) * sig2
  })
  vc.mat <- bdiag(vc.list)
  
  if(isDiagonal(vc.mat)) {
    vc.vec   <- diag(vc.mat)
    vc.names <- paste("D", 1:length(vc.vec), 1:length(vc.vec), sep="")
  } else{
    vc.vec <- as.matrix(vc.mat)[!upper.tri(vc.mat)]
    vc.index <- which(!upper.tri(vc.mat) == TRUE, arr.ind = TRUE)
    vc.names <- paste("D", vc.index[,1], vc.index[,2], sep="")
  }
  
  res <- c(sig2, vc.vec)
  names(res) <- c("sigma2", vc.names)
  
  #remove 0s 
  res <- res[res != 0]
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
  X <- lme4::getME(model, "X")
  
  n <- length(Y)
  
  flist <- lme4::getME(model, "flist")
  ngrps <- sapply(flist, function(x) length(levels(x)))
  
  # Constructing V = Cov(Y)
  sig0 <- lme4::getME(model, "sigma")
  
  ZDZt <- sig0^2 * crossprod( lme4::getME(model, "A") )
  R    <- Diagonal( n = n, x = sig0^2 )
  V    <- R + ZDZt
  
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
  X <- lme4::getME(model, "X")
  
  n <- length(Y)
  
  flist <- lme4::getME(model, "flist")
  ngrps <- sapply(flist, function(x) length(levels(x)))
  
  # Constructing V = Cov(Y)
  sig0 <- lme4::getME(model, "sigma")
  
  ZDZt <- sig0^2 * crossprod( lme4::getME(model, "A") )
  R    <- Diagonal( n = n, x = sig0^2 )
  V    <- R + ZDZt
  
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
               sig0 = sig0, V = V, V.chol = V.chol, Vinv = Vinv, XVXinv = XVXinv,
               M = M, P = P) )
  
}

# Extracting/calculating key matrices from lme object 
# @param model an lme object
.lme_matrices <- function(model) {
  design.info <- suppressWarnings(extract_design(model)) 
  
  Y <- design.info$Y
  Y <- Y[!is.na(Y)]
  
  X <- design.info$X
  Z <- Matrix( design.info$Z )
  
  D <- Matrix( design.info$D )
  
  n <- length(Y)
  
  flist <- model$groups
  ngrps <- sapply(flist, function(x) length(levels(x)))
  
  # Constructing V = Cov(Y)
  sig0 <- model$sigma
  V <-  Matrix(design.info$V)
  
  # Inverting V
  V.chol <- chol( V )
  Vinv   <- chol2inv( V.chol )
  
  # Calculating P
  XVXinv <- solve( t(X) %*% Vinv %*% X )
  VinvX  <- Vinv %*% X
  M      <- VinvX %*% XVXinv %*% t(VinvX)
  P      <- .Call("cxxmatsub", as.matrix(Vinv), as.matrix(M), 
                  PACKAGE = "HLMdiag")
  
  return( list(Y = Y, X = X, Z = Z, n = n, ngrps = ngrps, flist = flist,
               sig0 = sig0, V = V, Vinv = Vinv, XVXinv = XVXinv,
               M = M, P = P, D = D) )
  
}

# 'se.ranef' is a copy of function in arm package. This is copied to ensure
# that is available to all users. This should not be exported.
se.ranef <- function (object) 
{
  se.bygroup <- lme4::ranef(object, condVar = TRUE)
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
  cl <- class(object)
  if(cl == "lme") {
    RVAL <- TRUE
  } else {
    fl   <- object@flist
    fnms <- names(fl) 
    RVAL <- all(sapply(seq_along(fl)[-1], function(i) lme4::isNested(fl[[i-1]], fl[[i]])))
  }
  
  return(RVAL)
}


# Extract the residual covariance matrix from an lme object
.extractR.lme <- function(lme.fit) {
  n <- length( nlme::getResponse(lme.fit) )
  if (length(lme.fit$group) > 1) {
    stop("not implemented for multiple levels of nesting")
  } 
  else{
    ugroups <- unique(lme.fit$groups[[1]])
    if (!is.null(lme.fit$modelStruct$corStruct)) {
      V <- Matrix( nlme::corMatrix(lme.fit$modelStruct$corStruct) )
    }
    else V <- Diagonal(n)
  }
  if (!is.null(lme.fit$modelStruct$varStruct)) 
    sds <- 1/nlme::varWeights(lme.fit$modelStruct$varStruct)
  else sds <- rep(1, n)
  sds <- lme.fit$sigma * sds
  cond.var <- t(V * sds) * sds
  
  return(cond.var / lme.fit$sigma^2)
}

# Extract the marginal covariance matrix, V, for an lme object
.extractV.lme <- function(lme.fit) {
  n <- length( nlme::getResponse(lme.fit) )
  if (length(lme.fit$group) > 1) {
    stop("not implemented for multiple levels of nesting")
  } 
  else{
    ugroups <- unique(lme.fit$groups[[1]])
    if (!is.null(lme.fit$modelStruct$corStruct)) {
      cmat <- nlme::corMatrix(lme.fit$modelStruct$corStruct)
      V <- bdiag( nlme::corMatrix(lme.fit$modelStruct$corStruct) )
    }
    else V <- Diagonal(n)
  }
  if (!is.null(lme.fit$modelStruct$varStruct)) 
    sds <- 1/nlme::varWeights(lme.fit$modelStruct$varStruct)
  else sds <- rep(1, n)
  #   sds <- lme.fit$sigma * sds
  cond.var <- t(V * sds) * sds
  
  mod.mats <- .extract.lmeDesign(lme.fit)
  D <- Matrix( mod.mats$Vr )
  Z <- Matrix( mod.mats$Z )
  
  RES <- cond.var + Z %*% D %*% t(Z)
  
  return(RES)
}


#' Extracting covariance matrices from lme
#' 
#' This function extracts the full covariance matrices from a mixed/hierarchical
#' linear model fit using \code{lme}.
#' 
#' @export
#' @rdname extract_design
#' @aliases extract_design
#' @return A list of matrices is returned.
#' \itemize{
#' \item{\code{D} contains the covariance matrix of the random effects.}
#' \item{\code{V} contains the covariance matrix of the response.}
#' \item{\code{X} contains the fixed-effect model matrix.}
#' \item{\code{Z} contains the random-effect model matrix.}}
#' @param b a fitted model object of class \code{lme}.
#' @author Adam Loy \email{loyad01@@gmail.com}
#' @references This method has been adapted from the method
#'   \code{mgcv::extract.lme.cov} in the \code{mgcv} package, written by Simon
#'   N. Wood \email{simon.wood@@r-project.org}.
extract_design <- function (b){
  if (!inherits(b, "lme")) 
    stop("object does not appear to be of class lme")
  data <- b$data
  grps <- nlme::getGroups(b)
  n <- length(grps)
  if (is.null(b$modelStruct$varStruct)) 
    w <- rep(b$sigma, n)
  else {
    w <- 1/nlme::varWeights(b$modelStruct$varStruct)
    group.name <- names(b$groups)
    order.txt <- paste("ind<-order(data[[\"", group.name[1], 
                       "\"]]", sep = "")
    if (length(b$groups) > 1) 
      for (i in 2:length(b$groups)) order.txt <- paste(order.txt, 
                                                       ",data[[\"", group.name[i], "\"]]", sep = "")
    order.txt <- paste(order.txt, ")")
    eval(parse(text = order.txt))
    w[ind] <- w
    w <- w * b$sigma
  }
  if (is.null(b$modelStruct$corStruct)) 
    V <- diag(n)
  else {
    c.m <- nlme::corMatrix(b$modelStruct$corStruct)
    if (!is.list(c.m)) 
      V <- c.m
    else {
      V <- matrix(0, n, n)
      gr.name <- names(c.m)
      n.g <- length(c.m)
      j0 <- 1
      ind <- ii <- 1:n
      for (i in 1:n.g) {
        j1 <- j0 + nrow(c.m[[i]]) - 1
        V[j0:j1, j0:j1] <- c.m[[i]]
        ind[j0:j1] <- ii[grps == gr.name[i]]
        j0 <- j1 + 1
      }
      V[ind, ] <- V
      V[, ind] <- V
    }
  }
  V <- as.vector(w) * t(as.vector(w) * V)
  X <- list()
  grp.dims <- b$dims$ncol
  Zt <- model.matrix(b$modelStruct$reStruct, data)
  cov <- as.matrix(b$modelStruct$reStruct)
  i.col <- 1
  n.levels <- length(b$groups)
  Z <- matrix(0, n, 0)
  
  for (i in 1:(n.levels)) {
    if (length(levels(b$groups[[n.levels - i + 1]])) == 
        1) {
      X[[1]] <- matrix(rep(1, nrow(b$groups)))
    }
    else {
      clist <- list(`b$groups[[n.levels - i + 1]]` = c("contr.treatment", 
                                                       "contr.treatment"))
      X[[1]] <- model.matrix(~b$groups[[n.levels - 
                                          i + 1]] - 1, contrasts.arg = clist)
    }
    X[[2]] <- Zt[, i.col:(i.col + grp.dims[i] - 1), drop = FALSE]
    i.col <- i.col + grp.dims[i]
    Z <- cbind(Z, tensor.prod.model.matrix(X))
  }
  Vr <- matrix(0, ncol(Z), ncol(Z))
  start <- 1
  for (i in 1:(n.levels)) {
    k <- n.levels - i + 1
    for (j in 1:b$dims$ngrps[i]) {
      stop <- start + ncol(cov[[k]]) - 1
      Vr[start:stop, start:stop] <- cov[[k]]
      start <- stop + 1
    }
  }
  Vr <- Vr * b$sigma^2
  V <- V + Z %*% Vr %*% t(Z)
  
  return(
    list(
      D = Vr / b$sigma^2,
      V = V,
      X = model.matrix(b, data = b$data),
      Z = Z,
      Y = nlme::getResponse(b)
    )
  )
}

#adding rows with NAs into returned tibble

.lmerMod_add_NArows <- function(model, frame, na.action, data) { 
  rownums <- NULL
  for (i in 1:length(na.action)) {
    rownums[i] <- na.action[[i]]  
  }
  df <- data %>%
    dplyr::anti_join(model@frame, by = colnames(model@frame)) %>%
    dplyr::select(colnames(model@frame))
  
  for (i in 1:nrow(df)) {
    frame <- tibble::add_row(frame, df[i,], .before = as.numeric(rownums[i]))
  }
  
  return(frame)
}

.lme_add_NArows <- function(model, frame, na.action, org.data, fixed.data) {
  rownums <- NULL
  for (i in 1:length(na.action)) {
    rownums[i] <- na.action[[i]]
  }
  
  df <- org.data %>%
    dplyr::anti_join(fixed.data, by = colnames(fixed.data)) %>%
    dplyr::select(colnames(fixed.data))
  
  for (i in 1:length(na.action)) {
    frame <- tibble::add_row(frame, df[i,], .before = as.numeric(rownums[i]))
  }
  
  return(frame)
}
