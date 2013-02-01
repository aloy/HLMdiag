#' @export
cooks.distance.case_delete <- function(model, ...){
  p <- length(model$fixef.original)

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

#--------------------------------------
#' @export
mdffits.case_delete <- function(object, ...){
  p <- length(object$fixef.original)

  if(is(object$fixef.delete, "matrix")) {
    groups <- rownames(object$fixef.delete, do.NULL = FALSE, prefix = "")
    MDFFITS <- NULL
    for(i in 1:length(groups)){
      change.fixef <- as.matrix(object$fixef.original - object$fixef.delete[i,])
      MDFFITS <- c(MDFFITS, t(change.fixef) %*% solve( as.matrix(object$vcov.delete[[i]]) ) %*% change.fixef / p)
    }
  }
  else{
    change.fixef <- as.matrix(object$fixef.original - object$fixef.delete)
    MDFFITS <- t(change.fixef) %*% solve( as.matrix(object$vcov.delete) ) %*% change.fixef / p
  }

  return(MDFFITS)
}

#--------------------------------------
#' @export
covtrace.case_delete <- function(object, ...){
  p <- length(object$fixef.original)

  if(is(object$vcov.delete, "list")) {
    groups <- rownames(object$fixef.delete, do.NULL = FALSE, prefix = "")
    COVTRACE <- NULL
    for(i in 1:length(groups)){
      V.original <- as.matrix(object$vcov.original)
      V.delete <- as.matrix(object$vcov.delete[[i]])
      COVTRACE <- c(COVTRACE, abs(sum(diag(solve(V.original) %*% V.delete)) - p))
    }
  }
  else{
    V.original <- as.matrix(object$vcov.original)
    V.delete <- as.matrix(object$vcov.delete)
    COVTRACE <- abs(sum(diag(solve(V.original) %*% V.delete)) - p)
  }

  return(COVTRACE)	
}

#--------------------------------------
#' @export
covratio.case_delete <- function(object, ...){
  if(is(object$vcov.delete, "list")) {
    groups <- rownames(object$fixef.delete, do.NULL = FALSE, prefix = "")
    COVRATIO <- NULL
    for(i in 1:length(groups)){
      V.original <- as.matrix(object$vcov.original)
      V.delete <- as.matrix(object$vcov.delete[[i]])
      COVRATIO <- c(COVRATIO, det(V.delete) / det(V.original))
    }
  }
  else{
    V.original <- as.matrix(object$vcov.original)
    V.delete <- as.matrix(object$vcov.delete)
    COVRATIO <- det(V.delete) / det(V.original)
  }
  
  return(COVRATIO)
}

#--------------------------------------
#' @export
rvc.case_delete <- function(object, ...){
	if(class(object$varcomp.delete) == "list") {
	  res <- do.call('rbind', lapply(object$varcomp.delete, function(x){ (x / object$varcomp.original) - 1}))
	}
  else{
    res <- (object$varcomp.delete / object$varcomp.original) - 1
  }
	return(res)
}

#' Calculating diagnostics for two-level hierarchical linear models.
#'
#' This group of functions is used to compute deletion diagnostics for a
#' two-level normal hierarchical model at both levels of the model.
#'
#' The primary function if \code{diagnostics} which returns either a
#' list or data frame of influence measures depending on whether
#' \code{type = "both"} or if only one aspect of the model is selected.
#' If \code{type = "both"}, then a list with Cook's distance, MDFFITS,
#' COVTRACE, and COVRATIO are returned for the fixed effects and
#' relative variance change (RVC) is returned for the variance components.
#'
#' The methods \code{cooks.distance}, \code{mdffits}, \code{covtrace},
#' \code{covratio}, and \code{rvc} can be used for direct computation
#' of the corresponding diagnostic quantities from an object of class
#' \code{case_delete}.
#'
#' @aliases cooks.distance.case_delete mdffits.case_delete covtrace.case_delete covratio.case_delete rvc.case_delete
#' @param delete an object containing the output returned by \code{case_delete()}
#' @author Adam Loy \email{aloy@@iastate.edu}
#' @references Christensen, R., Pearson, L.M., and Johnson, W. (1992),
#' ``Case-Deletion Diagnostics for Mixed Models,'' \emph{Technometrics},
#' 34, 38 -- 45.
#'
#' Dillane, D. (2005), ``Deletion Diagnostics for the Linear Mixed Model,''
#' Ph.D. thesis, Trinity College Dublin.
#'
#' Schabenberger, O. (2004),``Mixed Model Influence Diagnostics,''
#' in \emph{Proceedings of the Twenty-Ninth SAS Users Group International Conference},
#' SAS Users Group International.
#'
#' @examples
#' data(Oxboys, package = 'mlmRev')
#' fm <- lmer(formula = height ~ age + I(age^2) + (age + I(age^2)| Subject), data = Oxboys)
#' fmDel <- case_delete(model = fm, group = TRUE, type = "both")
#' fmDiag <- diagnostics(model = fm, delete = fmDel)
#' 
#' \dontrun{
#' library(mlmRev)
#' exm1 <- lmer(normexam ~ standLRT + sex + schgend + (1 | school), data = Exam)
#' exm1DEL <- case_delete(model = exm1, group = TRUE, type = "both")
#' exm1DIAG <- diagnostics(model = exm1, delete = exm1DEL)
#' }
#' @export
#' @keywords models regression
diagnostics <- function(delete){
  type <- attributes(delete)$type
  if(type %in% c("fixef", "both")){
    ids <- as.vector(rownames(delete$fixef.delete, do.NULL = FALSE, prefix = ""))
  }
  else{
    ids <- as.vector(names(delete$varcomp.delete))
    if(is.null(ids)) ids <- 1:length(delete$varcomp.delete)
  }
  if(type  %in% c("fixef", "both")){
    if(!is(delete$fixef.delete, "matrix")) {
      res1 <- data.frame(COOKSD = cooks.distance(delete),
                        MDFFITS = mdffits(delete),
                        COVTRACE = covtrace(delete),
                        COVRATIO = covratio(delete))
    }
    else {
      res1 <- data.frame(IDS = ids, COOKSD = cooks.distance(delete),
                         MDFFITS = mdffits(delete),
                         COVTRACE = covtrace(delete),
                         COVRATIO = covratio(delete))
    }
    if(type == "fixef") return(res1)
  }

  if(type %in% c("varcomp", "both")){
    if(!is(delete$varcomp.delete, "list")) { res2 <- data.frame(rvc(delete)) }
    else res2 <- data.frame(IDS = ids, rvc(delete))
    
    if(type == "varcomp") return(res2)
  }

  if(type == "both"){
    res <- list(fixef_diag = res1, varcomp_diag = res2)
    return(res)
  } 
}
