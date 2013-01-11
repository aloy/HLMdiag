#' Implementing Case Deletion.
#'
#' This function is used to iteratively delete groups corresponding to the levels of a 
#' two-level hierarchical linear model. It uses \code{lmer()} to fit the models for each
#' deleted case (i.e. uses brute force). To investigate numerous levels of
#' the model, the function will need to be called multiple times, specifying
#' the group (level) of interest each time.
#'
#' @param model the original hierarchical model fit using \code{lmer()}
#' @param group a variable used to define the group for which cases will be deleted.
#'   If this is left \code{FALSE}, then the function will delete individual observations.
#' @param type the part of the model for which you are obtaining deletion diagnostics:
#'   the fixed effects (\code{fixef}), variance components (\code{varcomp}), or \code{both}
#' @return a list with the following compontents: 
#'   \item{fixef.original}{the original fixed effects}
#'   \item{ranef.original}{the origingal random effects}
#'   \item{vcov.original}{the original variance-covariance parameters}
#'   \item{varcomp.original}{the original variance components}
#'   \item{fixef.delete}{a list of the fixed effects obtained through case deletion}
#'   \item{ranef.delete}{a list of the random effects obtained through case deletion}
#'   \item{vcov.delete}{a list of the variance-covariance parameters obtained 
#'   through case deletion}
#'   \item{fitted.delete}{a list of the fitted values obtained through case deletion}
#'   \item{varcomp.delete}{a list of the variance components obtained through case
#'   deletion}
#' @author Adam Loy \email{aloy@@istate.edu}
#' @references Christensen, R., Pearson, L.M., and Johnson, W. (1992),
#' ``Case-Deletion Diagnostics for Mixed Models,'' \emph{Technometrics},
#' 34, 38 -- 45.
#'
#' Schabenberger, O. (2004),``Mixed Model Influence Diagnostics,''
#' in \emph{Proceedings of the Twenty-Ninth SAS Users Group International Conference},
#' SAS Users Group International.
#' @examples
#' data(Oxboys, package = 'mlmRev')
#' fm <- lmer(formula = height ~ age + I(age^2) + (age + I(age^2)| Subject), data = Oxboys)
#' fmDel <- case_delete(model = fm, group = TRUE, type = "both")
#'
#' \dontrun{library(mlmRev)
#' exm1 <- lmer(normexam ~ standLRT + sex + schgend + (1 | school), data = Exam)
#' exm1DEL <- case_delete(model = exm1, group = TRUE, type = "both")}
#' @export
case_delete <- function(model, group = NULL, type = c("both", "fixef", "varcomp"), 
                        delete = NULL){
  if(!is(model, "mer")) stop("model must be of class 'mer'")
  if(!model@dims["LMM"]){
    stop("case_delete is currently only implemented for mixed/hierarchical models.")
  }
  
  fixef.delete   <- NULL
  vcov.delete    <- NULL
  varcomp.delete <- NULL
  ranef.delete   <- NULL
  fitted.delete  <- NULL
  
  type <- match.arg(type) #default is "both"
  if( is.null(group) ){ # SINGLE CASE DELETION DIAGNOSTICS
    n <- model@dims[["n"]]
    modframe <- model@frame

    if( is.null(delete) ) {
      for(i in 1:n){
        model.delete <- lmer(formula = formula(model), data = model@frame[-i,])
        
        if(type %in% c("both", "varcomp")){
          if(length(getME(model.delete, "flist")) == 1) {
            ranef.delete[[i]] <- data.frame(deleted = i, 
                                            id = rownames(ranef(model.delete)[[1]]), 
                                            ranef(model.delete)[[1]])
          }
          else{
            ranef.delete[[i]] <- ranef(model.delete)
            ranef.delete[[i]] <- lapply(ranef.delete[[i]], function(x){
              x$id <- rownames(x)
              x$deleted <- i
              return(x)
            })
          }
          
          varcomp.delete[[i]] <- varcomp.mer(model.delete)
        }
        
        if(type %in% c("both", "fixef")){
          fixef.delete[[i]] <- c(deleted = i, fixef(model.delete))
          vcov.delete[[i]]  <- as.matrix(vcov(model.delete))
        }
        
        fitted.delete[[i]] <- data.frame(deleted = i, model.delete@frame, fitted(model.delete))
        
      }
    }
    else {
      model.delete   <- lmer(formula = formula(model), data = model@frame[-delete,])
      
      if(type %in% c("both", "fixef")) {
        fixef.delete   <- fixef(model.delete)
        vcov.delete    <- as.matrix(vcov(model.delete))
      }
      
      if(type %in% c("both", "varcomp")) {
        varcomp.delete <- varcomp.mer(model.delete)
        ranef.delete   <- ranef(model.delete)
        if( length(flist) == 1 ) ranef.delete <- ranef.delete[[1]]
      }
      fitted.delete  <- fitted(model.delete)
    }
  }

  else{ # MULTIPLE CASE DELETION DIAGNOSTICS
    flist <- model@flist
    
    if(!group %in% names(flist)) {
      stop(paste(group, "is not a valid grouping factor for this model."))
    }
    
    
    if( is.null(delete) ){
      data.delete <- split(model@frame, model@frame[, group])
      data.delete <- lapply(data.delete, function(df){
        data.delete[[ unique( df[, group ] ) ]] <- NULL
        do.call('rbind', data.delete)
      })
      
      model.delete <- lapply(data.delete, lmer, formula = formula(model))
      
      
      if(length(flist) == 1) {
        ranef.delete <- lapply(model.delete, function(x){
          data.frame(deleted = setdiff(model@frame[, group], x@frame[, group]),
                     id = rownames(ranef(x)[[1]]), ranef(x)[[1]])
        })
      }
      else{
        ranef.delete  <- lapply(model.delete, ranef)
        deleted.group <- rownames(ranef(model)[[group]])
        
        ranef.delete <- lapply(1:length(ranef.delete), function(x){
          ranef.list <- ranef.delete[[x]]
          lapply(ranef.list, function(y) {
            y$id <- rownames(y)
            y$deleted <- deleted.group[x]
            return(y)
          })
        })
      }
      
      varcomp.delete <- lapply(model.delete, varcomp.mer)
      
      if(type %in% c("both", "fixef")){
        fixef.delete <- lapply(model.delete, fixef)
        
        vcov.delete <- lapply(model.delete, vcov)
        vcov.delete <- lapply(vcov.delete, as.matrix)
      }
      
      
      fitted.delete <- lapply(model.delete, function(x){
        data.frame(deleted = setdiff(model@frame[, group], x@frame[, group]),
                   x@frame, fitted(x))
      })
    }
    else{
      index <- !model@frame[,group] %in% delete
      model.delete   <- lmer(formula = formula(model), data = model@frame[index,])
      
      if(type %in% c("both", "fixef")) {
        fixef.delete   <- fixef(model.delete)
        vcov.delete    <-  as.matrix(vcov(model.delete))
      }
      
      if(type %in% c("both", "varcomp")) {
        varcomp.delete <- varcomp.mer(model.delete)
        ranef.delete   <- ranef(model.delete)
        if( length(flist) == 1 ) ranef.delete <- ranef.delete[[1]]
      }
      fitted.delete  <- fitted(model.delete)
    }
    
    
  }
  
  # Organizing results
  if(is.null(delete)) {
    if(type %in% c("both", "fixef")){
      fitted.delete <- do.call('rbind', fitted.delete)
      #if(model@dims[["p"]] > 1) 
      fixef.delete  <- do.call('rbind', fixef.delete)
    }
    
    
    if(type %in% c("both", "varcomp")){
      if(length(getME(model, "flist")) == 1) {
        ranef.delete <- do.call('rbind', ranef.delete)
      }
      else {
        flist <- names(model@flist)
        temp  <- NULL
        for(i in 1:length(flist)) {
          temp[[i]] <- ldply(ranef.delete, function(x) x[[i]])
        }
        ranef.delete <- temp
        names(ranef.delete) <- names(ranef(model))
      }
    }
  }
  
  fixef.original <- model@fixef
  ranef.original <- ranef(model)
  if(length(ranef.original) == 1) ranef.original <- ranef.original[[1]]

  vcov.original <- as.matrix(vcov(model))
  varcomp.original <- varcomp.mer(model)

  val <- list(fixef.original = fixef.original, ranef.original = ranef.original,
              vcov.original = vcov.original, varcomp.original = varcomp.original,
              fixef.delete = fixef.delete, ranef.delete = ranef.delete,
              vcov.delete = vcov.delete, fitted.delete = fitted.delete,
              varcomp.delete = varcomp.delete)

  attr(val, "type") <- type
  class(val) <- "case_delete"
  return(val)
}