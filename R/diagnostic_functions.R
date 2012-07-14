#' @export
cooksd_hlm <- function(model, delete){
  groups <- rownames(delete$fixef.delete, do.NULL = FALSE, prefix = "")
  rank.X <- qr(model.matrix(model))$rank

  cook <- NULL
  for(i in 1:length(groups)){
    change.fixef <- as.matrix(delete$fixef.original - delete$fixef.delete[i,])
    cook <- c(cook, t(change.fixef) %*% ginv(as.matrix(vcov(model))) %*% change.fixef / rank.X)
  }

  return(cook)
}

#--------------------------------------
#' @export
mdffits_hlm <- function(model, delete){
  groups <- rownames(delete$fixef.delete, do.NULL = FALSE, prefix = "")
  rank.X <- qr(model.matrix(model))$rank

  MDFFITS <- NULL
  for(i in 1:length(groups)){
    change.fixef <- as.matrix(delete$fixef.original - delete$fixef.delete[i,])
    MDFFITS <- c(MDFFITS, t(change.fixef) %*% ginv(as.matrix(delete$vcov.delete[[i]])) %*% change.fixef / rank.X)
  }	

  return(MDFFITS)
}

#--------------------------------------
#' @export
covtrace_hlm <- function(model, delete){
  groups <- rownames(delete$fixef.delete, do.NULL = FALSE, prefix = "")
  rank.X <- qr(model.matrix(model))$rank

  COVTRACE <- NULL
  for(i in 1:length(groups)){
    V.original <- as.matrix(vcov(model))
    V.delete <- as.matrix(delete$vcov.delete[[i]])
    COVTRACE <- c(COVTRACE, abs(sum(diag(ginv(V.original) %*% V.delete)) - rank.X))
  }

  return(COVTRACE)	
}

#--------------------------------------
#' @export
covratio_hlm <- function(model, delete){
  groups <- rownames(delete$fixef.delete, do.NULL = FALSE, prefix = "")

  COVRATIO <- NULL
  for(i in 1:length(groups)){
    V.original <- as.matrix(vcov(model))
    V.delete <- as.matrix(delete$vcov.delete[[i]])
    COVRATIO <- c(COVRATIO, det(V.delete) / det(V.original))
  }

  return(COVRATIO)
}

#--------------------------------------
#' @export
rvc <- function(delete){
	res <- do.call('rbind', lapply(delete$varcomp.delete, function(x){ (x / delete$varcomp.original) - 1}))
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
#' If \code{type = "both}, then a list with Cook's distance, MDFFITS,
#' COVTRACE, and COVRATIO are returned for the fixed effects and
#' relative variance change (RVC) is returned for the variance components.
#'
#' The functions \code{cooksd_hlm}, \code{mdffits_hlm}, \code{covtrace_hlm},
#' \code{covratio_hlm}, and \code{rvc} can be used for direct computation
#' of the corresponding diagnostic quantities.
#'
#' @aliases cooksd_hlm mdffits_hlm covtrace_hlm covratio_hlm rvc
#' @param model an object contatining the original hierarchical model fit using \code{lmer()}
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
diagnostics <- function(model, delete){
  type <- attributes(delete)$type
  if(type %in% c("fixef", "both")){
    ids <- as.vector(rownames(delete$fixef.delete, do.NULL = FALSE, prefix = ""))
  }
  else{
    ids <- as.vector(names(delete$varcomp.delete))
    if(is.null(ids)) ids <- 1:length(delete$varcomp.delete)
  }
  if(type  %in% c("fixef", "both")){
    res1 <- data.frame(IDS = ids, COOKSD = cooksd_hlm(model, delete),
                       MDFFITS = mdffits_hlm(model, delete),
                       COVTRACE = covtrace_hlm(model, delete),
                       COVRATIO = covratio_hlm(model, delete))
    if(type == "fixef") return(res1)
  }

  if(type %in% c("varcomp", "both")){
    res2 <- data.frame(IDS = ids, rvc(delete))
    if(type == "varcomp") return(res2)
  }

  if(type == "both"){
    res <- list(fixef_diag = res1, varcomp_diag = res2)
    return(res)
  } 
}
