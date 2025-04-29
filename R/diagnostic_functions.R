#' @export
#' @rdname hlm_influence
#' @param model an object containing the output returned by \code{case_delete()}.
#' This is only named differently to agree with the generic.
#' @method cooks.distance case_delete
cooks.distance.case_delete <- function(model, ...){
  p <- length(model$fixef.original)
  
  if (colnames(model$fixef.delete)[1] == "deleted") {
    model$fixef.delete <- model$fixef.delete[,-1]
  }

  if(is(model$fixef.delete, "matrix")) {
    groups <- rownames(model$fixef.delete, do.NULL = FALSE, prefix = "")
    cook <- NULL
    for(i in 1:length(groups)){
      change.fixef <- as.matrix(model$fixef.original - model$fixef.delete[i,])
      cook <- c(cook, t(change.fixef) %*% solve( as.matrix( model$vcov.original ) ) %*% change.fixef / p)
    }
  }
  else{
    change.fixef <- as.matrix(model$fixef.original - model$fixef.delete)
    cook <- t(change.fixef) %*% solve( as.matrix( model$vcov.original ) ) %*% change.fixef / p
  }

  return(cook)
}


#' @export
#' @rdname hlm_influence
#' @method mdffits case_delete
mdffits.case_delete <- function(model, ...){
  p <- length(model$fixef.original)
  
  if (colnames(model$fixef.delete)[1] == "deleted") {
    model$fixef.delete <- model$fixef.delete[,-1]
  }

  if(is(model$fixef.delete, "matrix")) {
    groups <- rownames(model$fixef.delete, do.NULL = FALSE, prefix = "")
    MDFFITS <- NULL
    for(i in 1:length(groups)){
      change.fixef <- as.matrix(model$fixef.original - model$fixef.delete[i,])
      MDFFITS <- c(MDFFITS, t(change.fixef) %*% solve( as.matrix(model$vcov.delete[[i]]) ) %*% change.fixef / p)
    }
  }
  else{
    change.fixef <- as.matrix(model$fixef.original - model$fixef.delete)
    MDFFITS <- t(change.fixef) %*% solve( as.matrix(model$vcov.delete) ) %*% change.fixef / p
  }

  return(MDFFITS)
}


#' @export
#' @rdname hlm_influence
#' @method covtrace case_delete
covtrace.case_delete <- function(model, ...){
  p <- length(model$fixef.original)

  if(is(model$vcov.delete, "list")) {
    groups <- rownames(model$fixef.delete, do.NULL = FALSE, prefix = "")
    COVTRACE <- NULL
    for(i in 1:length(groups)){
      V.original <- as.matrix(model$vcov.original)
      V.delete <- as.matrix(model$vcov.delete[[i]])
      COVTRACE <- c(COVTRACE, abs(sum(diag(solve(V.original) %*% V.delete)) - p))
    }
  }
  else{
    V.original <- as.matrix(model$vcov.original)
    V.delete <- as.matrix(model$vcov.delete)
    COVTRACE <- abs(sum(diag(solve(V.original) %*% V.delete)) - p)
  }

  return(COVTRACE)	
}



#' @export
#' @rdname hlm_influence
#' @method covratio case_delete
covratio.case_delete <- function(model, ...){
  if(is(model$vcov.delete, "list")) {
    groups <- rownames(model$fixef.delete, do.NULL = FALSE, prefix = "")
    COVRATIO <- NULL
    for(i in 1:length(groups)){
      V.original <- as.matrix(model$vcov.original)
      V.delete <- as.matrix(model$vcov.delete[[i]])
      COVRATIO <- c(COVRATIO, det(V.delete) / det(V.original))
    }
  }
  else{
    V.original <- as.matrix(model$vcov.original)
    V.delete <- as.matrix(model$vcov.delete)
    COVRATIO <- det(V.delete) / det(V.original)
  }
  
  return(COVRATIO)
}



#' @export
#' @rdname hlm_influence
#' @param ... do not use
#' @method rvc case_delete
rvc.case_delete <- function(model, ...){
	if(inherits(model$varcomp.delete, "list")) {
	  res <- do.call('rbind', lapply(model$varcomp.delete, function(x){ (x / model$varcomp.original) - 1}))
	}
  else{
    res <- (model$varcomp.delete / model$varcomp.original) - 1
  }
	return(res)
}
