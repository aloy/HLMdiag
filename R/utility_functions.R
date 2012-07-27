#' Reorganizing Z matrix
#' 
#' This function reorganizes the Z matrix obtained from an \code{mer} object
#' into the block form discussed by Demidenko (2004). Currently, this
#' function assumes there is only one level and that units are nested.
#' 
#' @param object a fitted model object of class \code{mer}.
#' 
#' @references Demidenko, E. (2004). Mixed Models: Theory and Applications.
#'             New York, Wiley.
#' @author Adam Loy \email{aloy@@iastate.edu}
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