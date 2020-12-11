#' @export
case_delete <- function(model, ...){
  UseMethod("case_delete", model)
}

#' @export
#' @rdname case_delete.mer
#' @method case_delete default
case_delete.default <- function(model, ...){
  stop(paste("there is no case_delete() method for objects of class",
             paste(class(model), collapse=", ")))
}

#' Case Deletion for \code{mer}/\code{lmerMod} objects
#'
#'This function is used to iteratively delete groups corresponding to the
#'levels of a hierarchical linear model. It uses \code{lmer()} to fit
#'the models for each deleted case (i.e. uses brute force). To investigate
#'numerous levels of the model, the function will need to be called multiple
#'times, specifying the group (level) of interest each time.
#'
#' @export
#' @method case_delete mer
#' @aliases case_delete
#'@param model the original hierarchical model fit using \code{lmer()}
#'@param level a variable used to define the group for which cases will be
#'deleted.  If \code{level = 1} (default), then the function will delete
#'individual observations.
#'@param type the part of the model for which you are obtaining deletion
#'diagnostics: the fixed effects (\code{"fixef"}), variance components
#'(\code{"varcomp"}), or \code{"both"} (default).
#'@param delete numeric index of individual cases to be deleted. If the \code{level} parameter 
#'is specified, \code{delete} may also take the form of a character vector consisting of group 
#'names as they appear in \code{flist}. It is possible to set \code{level} and delete individual
#'cases from different groups using \code{delete}, so numeric indices should be double checked 
#'to confirm that they encompass entire groups. If \code{delete = NULL} then all cases are iteratively deleted.
#' @param ... do not use
#'@return a list with the following components:
#' \describe{
#'   \item{\code{fixef.original}}{the original fixed effects estimates}
#'   \item{\code{ranef.original}}{the original predicted random effects}
#'   \item{\code{vcov.original}}{the original variance-covariance matrix for the fixed effects}
#'   \item{\code{varcomp.original}}{the original estimated variance components}
#'   \item{\code{fixef.delete}}{a list of the fixed effects estimated after case deletion}
#'   \item{\code{ranef.delete}}{a list of the random effects predicted after case deletion}
#'   \item{\code{vcov.delete}}{a list of the variance-covariance matrices for the fixed 
#'      effects obtained after case deletion}
#'   \item{\code{fitted.delete}}{a list of the fitted values obtained after case
#'      deletion}
#' \item{\code{varcomp.delete}}{a list of the estimated variance components obtained after
#'      case deletion}
#' }
#'@author Adam Loy \email{loyad01@@gmail.com}
#' @keywords models regression
#'@references Christensen, R., Pearson, L.M., and Johnson, W. (1992)
#'Case-Deletion Diagnostics for Mixed Models, \emph{Technometrics}, \bold{34}, 38
#'-- 45.
#'
#'Schabenberger, O. (2004) Mixed Model Influence Diagnostics, in
#'\emph{Proceedings of the Twenty-Ninth SAS Users Group International
#'Conference}, SAS Users Group International.
#'@examples
#'data(sleepstudy, package = 'lme4')
#'fm <- lme4::lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
#'
#' # Deleting every Subject
#' fmDel <- case_delete(model = fm, level = "Subject", type = "both")
#'
#' # Deleting only subject 308
#' del308 <- case_delete(model = fm, level = "Subject", type = "both", delete = 308)
#' 
#' # Deleting a subset of subjects
#' delSubset <- case_delete(model = fm, level = "Subject", type = "both", delete = 308:310)
#'
case_delete.mer <- function(model, level = 1, type = c("both", "fixef", "varcomp"), 
                        delete = NULL, ...){
  if(!is(model, "mer")) stop("model must be of class 'mer'")
  if(!model@dims["LMM"]){
    stop("case_delete is currently only implemented for mixed/hierarchical models.")
  }
  
  if (hasArg(group)) {
    group <- NULL
    warning("group is not a valid argument for this function. As of version 0.4.0, group has been replaced by level. See ?hlm_influence for more information.")
  }
  
  flist <- model@flist
  
  fixef.delete   <- NULL
  vcov.delete    <- NULL
  varcomp.delete <- NULL
  ranef.delete   <- NULL
  fitted.delete  <- NULL
  
  type <- match.arg(type) #default is "both"
  if( level == 1 ){ # SINGLE CASE DELETION DIAGNOSTICS
    n <- model@dims[["n"]]
    modframe <- model@frame

    if( is.null(delete) ) {
      for(i in 1:n){
        model.delete <- lme4::lmer(formula = formula(model), data = model@frame[-i,])
        
        if(type %in% c("both", "varcomp")){
          if(length(lme4::getME(model.delete, "flist")) == 1) {
            ranef.delete[[i]] <- data.frame(deleted = i, 
                                            id = rownames(lme4::ranef(model.delete)[[1]]), 
                                            lme4::ranef(model.delete)[[1]])
          }
          else{
            ranef.delete[[i]] <- lme4::ranef(model.delete)
            ranef.delete[[i]] <- lapply(ranef.delete[[i]], function(x){
              x$id <- rownames(x)
              x$deleted <- i
              return(x)
            })
          }
          
          varcomp.delete[[i]] <- varcomp.mer(model.delete)
        }
        
        if(type %in% c("both", "fixef")){
          fixef.delete[[i]] <- c(deleted = i, lme4::fixef(model.delete))
          vcov.delete[[i]]  <- as.matrix(vcov(model.delete))
        }
        
        fitted.delete[[i]] <- data.frame(deleted = i, model.delete@frame, fitted(model.delete))
        
      }
    }
    else {
      model.delete   <- lme4::lmer(formula = formula(model), data = model@frame[-delete,])
      
      if(type %in% c("both", "fixef")) {
        fixef.delete   <- lme4::fixef(model.delete)
        vcov.delete    <- as.matrix(vcov(model.delete))
      }
      
      if(type %in% c("both", "varcomp")) {
        varcomp.delete <- varcomp.mer(model.delete)
        ranef.delete   <- lme4::ranef(model.delete)
        if( length(flist) == 1 ) ranef.delete <- ranef.delete[[1]]
      }
      fitted.delete  <- fitted(model.delete)
    }
  }

  else{ # MULTIPLE CASE DELETION DIAGNOSTICS
    
    if(!level %in% names(flist)) {
      stop(paste(level, "is not a valid grouping factor for this model."))
    }
    
    
    if( is.null(delete) ){
      data.delete <- split(model@frame, model@frame[, level])
      data.delete <- lapply(data.delete, function(df){
        index <- unique( df[, level ] )
        if(class(index) != "character") index <- as.character(index)
        data.delete[[ index ]] <- NULL
        do.call('rbind', data.delete)
      })
      
      model.delete <- lapply(data.delete, lme4::lmer, formula = formula(model))
      
      
      if(length(flist) == 1) {
        ranef.delete <- lapply(model.delete, function(x){
          data.frame(deleted = setdiff(model@frame[, level], x@frame[, level]),
                     id = rownames(lme4::ranef(x)[[1]]), lme4::ranef(x)[[1]])
        })
      }
      else{
        ranef.delete  <- lapply(model.delete, lme4::ranef)
        deleted.group <- rownames(lme4::ranef(model)[[level]])
        
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
        fixef.delete <- lapply(model.delete, lme4::fixef)
        
        vcov.delete <- lapply(model.delete, vcov)
        vcov.delete <- lapply(vcov.delete, as.matrix)
      }
      
      
      fitted.delete <- lapply(model.delete, function(x){
        data.frame(deleted = setdiff(model@frame[, level], x@frame[, level]),
                   x@frame, fitted(x))
      })
    }
    else{
      index <- !model@frame[,level] %in% delete
      model.delete   <- lme4::lmer(formula = formula(model), data = model@frame[index,])
      
      if(type %in% c("both", "fixef")) {
        fixef.delete   <- lme4::fixef(model.delete)
        vcov.delete    <-  as.matrix(vcov(model.delete))
      }
      
      if(type %in% c("both", "varcomp")) {
        varcomp.delete <- varcomp.mer(model.delete)
        ranef.delete   <- lme4::ranef(model.delete)
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
      if(length(lme4::getME(model, "flist")) == 1) {
        ranef.delete <- do.call('rbind', ranef.delete)
      }
      else {
        flist <- names(model@flist)
        temp  <- NULL
        for(i in 1:length(flist)) {
          temp[[i]] <- plyr::ldply(ranef.delete, function(x) x[[i]])
        }
        ranef.delete <- temp
        names(ranef.delete) <- names(lme4::ranef(model))
      }
    }
  }
  
  fixef.original <- model@fixef
  ranef.original <- lme4::ranef(model)
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



#' @export
#' @rdname case_delete.mer
#' @method case_delete lmerMod
case_delete.lmerMod <- function(model, level = 1, type = c("both", "fixef", "varcomp"), 
                            delete = NULL, ...){
  if(!isNestedModel(model)){
    stop("case_delete is currently only implemented for mixed/hierarchical models.")
  }
  
  if (hasArg(group)) {
    group <- NULL
    warning("group is not a valid argument for this function. As of version 0.4.0, group has been replaced by level. See ?hlm_influence for more information.")
  }
  
  flist <- model@flist
  
  fixef.delete   <- NULL
  vcov.delete    <- NULL
  varcomp.delete <- NULL
  ranef.delete   <- NULL
  fitted.delete  <- NULL
  
  type <- match.arg(type) #default is "both"
  if( level == 1 ){ # SINGLE CASE DELETION DIAGNOSTICS
    n <- lme4::getME(model, "n")
    modframe <- model@frame
    
    if( is.null(delete) ) {
      for(i in 1:n){
        model.delete <- lme4::lmer(formula = formula(model), data = model@frame[-i,])
        
        if(type %in% c("both", "varcomp")){
          if(length(lme4::getME(model.delete, "flist")) == 1) {
            ranef.delete[[i]] <- data.frame(deleted = i, 
                                            id = rownames(lme4::ranef(model.delete)[[1]]), 
                                            lme4::ranef(model.delete)[[1]])
          }
          else{
            ranef.delete[[i]] <- lme4::ranef(model.delete)
            ranef.delete[[i]] <- lapply(ranef.delete[[i]], function(x){
              x$id <- rownames(x)
              x$deleted <- i
              return(x)
            })
          }
          
          varcomp.delete[[i]] <- varcomp.mer(model.delete)
        }
        
        if(type %in% c("both", "fixef")){
          fixef.delete[[i]] <- c(deleted = i, lme4::fixef(model.delete))
          vcov.delete[[i]]  <- as.matrix(vcov(model.delete))
        }
        
        fitted.delete[[i]] <- data.frame(deleted = i, model.delete@frame, fitted(model.delete))
        
      }
    }
    else { 
      if (!is.numeric(delete)) {
        stop("For individual case deletion, the delete parameter should be a numeric vector")
      }
      model.delete   <- lme4::lmer(formula = formula(model), data = model@frame[-delete,])
      
      if(type %in% c("both", "fixef")) {
        fixef.delete   <- lme4::fixef(model.delete)
        vcov.delete    <- as.matrix(vcov(model.delete))
      }
      
      if(type %in% c("both", "varcomp")) {
        varcomp.delete <- varcomp.mer(model.delete)
        ranef.delete   <- lme4::ranef(model.delete)
        if( length(flist) == 1 ) ranef.delete <- ranef.delete[[1]]
      }
      fitted.delete  <- fitted(model.delete)
    }
  }
  
  else{ # MULTIPLE CASE DELETION DIAGNOSTICS
    
    if(!level %in% names(flist)) {
      stop(paste(level, "is not a valid grouping factor for this model."))
    }
    
    
    if( is.null(delete) ){
    
      data.delete <- split(model@frame, flist[level])
      data.delete <- lapply(data.delete, function(df) {
        df <- dplyr::anti_join(model@frame, df, by = names(model@frame))
      })
      
      
      
      model.delete <- lapply(data.delete, lme4::lmer, formula = formula(model)) #original 
    
     
      if(length(flist) == 1) {
        ranef.delete <- lapply(model.delete, function(x){
          data.frame(deleted = setdiff(model@frame[, level], x@frame[, level]),
                     id = rownames(lme4::ranef(x)[[1]]), lme4::ranef(x)[[1]])
        }) 
      }
      
      else{
        ranef.delete  <- lapply(model.delete, lme4::ranef)
        deleted.group <- rownames(lme4::ranef(model)[[level]])
        
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
        fixef.delete <- lapply(model.delete, lme4::fixef)
        
        vcov.delete <- lapply(model.delete, vcov)
        vcov.delete <- lapply(vcov.delete, as.matrix)
      }
     
      fitted.delete <- lapply(model.delete, function(x) {
        data.frame(deleted = setdiff(model@flist[[level]], x@flist[[level]]), x@frame, fitted(x))
      })
    }
    
    else{
      deleted_levels <- NULL
      
      if(is.numeric(delete)) {
        index <- rep(TRUE, length(flist[[level]]))
        index[delete] <- FALSE
        deleted_levels <- as.vector(flist[[level]][delete])
      }
      else if (is.character(delete)) {
        for (i in 1:length(delete)) {
          if(!delete[i] %in% flist[[level]]) {
            stop(paste(delete[i], "is not a valid group name to delete for this model. Names should follow the same format as in model@flist[[group]]. An example of an acceptable group name is:", model@flist[[level]][1]))
          }
        }
        index <- unlist(lapply(list(delete), function(s) {
          s <- !as.vector(model@flist[[level]]) %in% s
        }))
      }
      else {
        stop("Delete must either be a numeric or character vector.")
      }
     
      model.delete   <- lme4::lmer(formula = formula(model), data = model@frame[index,])
      
      if (sum(deleted_levels %in% model.delete@flist[[level]]) > 0) {
        warning("Level parameter is specified, but deleted cases do not encompass entire groups")
      }
      
      if(type %in% c("both", "fixef")) {
        fixef.delete   <- lme4::fixef(model.delete)
        vcov.delete    <-  as.matrix(vcov(model.delete))
      }
      
      if(type %in% c("both", "varcomp")) {
        varcomp.delete <- varcomp.mer(model.delete)
        ranef.delete   <- lme4::ranef(model.delete)
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
      if(length(lme4::getME(model, "flist")) == 1) {
        ranef.delete <- do.call('rbind', ranef.delete)
      }
      else {
        flist <- names(model@flist)
        temp  <- NULL
        for(i in 1:length(flist)) {
          temp[[i]] <- plyr::ldply(ranef.delete, function(x) x[[i]])
        }
        ranef.delete <- temp
        names(ranef.delete) <- names(lme4::ranef(model))
      }
    }
  }
  
  fixef.original <- lme4::fixef(model)
  ranef.original <- lme4::ranef(model)
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


#' @export
#' @rdname case_delete.mer
#' @method case_delete lme
#' @aliases case_delete
case_delete.lme <- function(model, level = 1, type = c("both", "fixef", "varcomp"), delete = NULL, ...){
  if(!isNestedModel(model)){
    stop("case_delete is currently only implemented for mixed/hierarchical models.")
  }
  
  flist <- model$groups
  
  fixef.delete   <- NULL
  vcov.delete    <- NULL
  varcomp.delete <- NULL
  ranef.delete   <- NULL
  fitted.delete  <- NULL
  
  type <- match.arg(type) #default is "both"
  if(level == 1) { # SINGLE CASE DELETION DIAGNOSTICS
    n <- model$dims$N
    modframe <- model$data
    dataformula <- formula(modframe)
    randcall <- as.formula(model$call$random)
    
    if( is.null(delete) ) {
      for(i in 1:n){
        if(is.null(randcall)) {
          model.delete <- nlme::lme(formula(model), data = modframe[-i,])
        } 
        else{
          model.delete <- nlme::lme(formula(model), random = randcall, data = modframe[-i,])
        } 
        
        
        if(type %in% c("both", "varcomp")){
          if(length(flist) == 1) {
            ranef.delete[[i]] <- data.frame(deleted = i, 
                                            id = rownames(nlme::ranef(model.delete)), 
                                            nlme::ranef(model.delete))
          }
          else{
            ranef.delete[[i]] <- nlme::ranef(model.delete)
            
            if(is.list(ranef.delete[[i]])) { 
              ranef.delete[[i]] <- lapply(ranef.delete[[i]], function(x){
                x$id <- rownames(x)
                x$deleted <- i
                return(x)
              })
            } else{
              ranef.delete[[i]]$id <- rownames(ranef.delete[[i]])
              ranef.delete[[i]]$deleted <- i
            }
            
          }
          
          varcomp.delete[[i]] <- varcomp.lme(model.delete) 
          
          
        } 
        
        if(type %in% c("both", "fixef")){
          fixef.delete[[i]] <- c(deleted = i, nlme::fixef(model.delete))
          vcov.delete[[i]]  <- as.matrix(vcov(model.delete))
        }
        
        
        fixed <- formula(model)
        dataform <- paste(fixed[2], "~", fixed[3], " + ",
                          paste(names(model.delete$groups), collapse = " + "))
        data <- model.delete$data %>%
          dplyr::mutate(across(where(is.character), ~ as.factor(.x))) %>%
          as.data.frame()
        new.frame <- model.frame(formula(dataform), data)
        
        
       
        fitted.delete[[i]] <- data.frame(deleted = i, new.frame, fitted(model.delete))
      }
    }  
    else {
      if (!is.numeric(delete)) {
        stop("For individual case deletion, the delete parameter should be a numeric vector")
      }
      
      if(is.null(randcall)) {
        model.delete <- nlme::lme(formula(model), data = modframe[-delete,])
      } else{
        model.delete <- nlme::lme(formula(model), data = modframe[-delete,], random = randcall)
      }
      
      if(type %in% c("both", "fixef")) {
        fixef.delete   <- nlme::fixef(model.delete)
        vcov.delete    <- as.matrix(vcov(model.delete))
      }
      
      if(type %in% c("both", "varcomp")) {
        varcomp.delete <- varcomp.lme(model.delete)
        ranef.delete   <- nlme::ranef(model.delete)
        if( length(flist) == 1 ) ranef.delete <- ranef.delete[[1]]
      }
      fitted.delete  <- fitted(model.delete)
    }
  }
  
  else{ # MULTIPLE CASE DELETION DIAGNOSTICS
    modframe <- model$data
    dataformula <- formula(modframe)
    modframe <- as.data.frame(modframe)
    randcall <- as.formula(model$call$random)
    
    if(!level %in% names(flist)) {
      stop(paste(level, "is not a valid grouping factor for this model."))
    }
    
    
    if( is.null(delete) ){
      
      data.delete <- split(modframe, modframe[, level])[1:length(unique(model$groups[[level]]))] #extract non-empty ones 
      data.delete <- lapply(data.delete, function(df) {
        df <- dplyr::anti_join(model$data, df, by = names(model$data))
      })
      
  
      if(is.null(randcall)) {
        model.delete <- lapply(data.delete, nlme::lme, fixed = formula(model))
      } else{
        model.delete <- lapply(data.delete, nlme::lme, fixed = formula(model), random = randcall)
      }
      
      
      
      if(length(flist) == 1) {
        ranef.delete <- lapply(model.delete, function(x){
          data.frame(deleted = setdiff(modframe[, level], x$data[, level]),
                     id = rownames(nlme::ranef(x)), nlme::ranef(x))
        })
      }
      else{
        ranef.delete  <- lapply(model.delete, nlme::ranef)
        deleted.group <- rownames(nlme::ranef(model))
        
        ranef.delete <- lapply(1:length(ranef.delete), function(x){
          ranef.list <- ranef.delete[[x]]
          lapply(ranef.list, function(y) {
            y$id <- rownames(y)
            y$deleted <- deleted.group[x]
            return(y)
          })
        })
      }
      
      varcomp.delete <- lapply(model.delete, varcomp.lme)
      
      if(type %in% c("both", "fixef")){
        fixef.delete <- lapply(model.delete, nlme::fixef)
        
        vcov.delete <- lapply(model.delete, vcov)
        vcov.delete <- lapply(vcov.delete, as.matrix)
      }
      
      #create new frame of data with variables used in the model 
      fixed <- formula(model)
      for (i in 1:length(model.delete)) {
        dataform <- paste(fixed[2], "~", fixed[3], " + ",
                          paste(names(model.delete[[i]]$groups), collapse = " + "))
        data <- model.delete[[i]]$data %>%
          dplyr::mutate(across(where(is.character), ~ as.factor(.x))) %>%
          as.data.frame()
        model.delete[[i]]$newframe <- model.frame(formula(dataform), data)
      }
      
      
      fitted.delete <- lapply(model.delete, function(x){
        data.frame(deleted = setdiff(modframe[, level], x$data[, level]),
        x$newframe, fitted(x))
        })
      
    }
    else{
      deleted_levels <- NULL
      
      if(is.numeric(delete)) {
        index <- rep(TRUE, length(model$groups[[level]]))
        index[delete] <- FALSE
        deleted_levels <- as.vector(model$groups[[level]][delete])
      }
      else if (is.character(delete)) {
        for (i in 1:length(delete)) {
          if(!delete[i] %in% model$groups[[level]]) {
            stop(paste(delete[i], "is not a valid group name to delete for this model. 
                       Names should follow the same format as in model$groups[[level]]. 
                       An example of an acceptable group name is:", model$groups[[level]][[1]]))
            }
          }
        index <- unlist(lapply(list(delete), function(s) {
          s <- !as.vector(model$groups[[level]]) %in% s
        }))
      }
      else {
        stop("Delete must either be a numeric or character vector.")
      }
      
      
      if(is.null(randcall)) {
        model.delete   <- nlme::lme(formula(model), data = modframe[index,])
      } else{
        model.delete   <- nlme::lme(formula(model), data = modframe[index,], random = randcall)
      }
      
      if (sum(deleted_levels %in% model.delete$groups[[level]]) > 0) {
        warning("Level parameter is specified, but deleted cases do not encompass entire groups")
      }
      
      if(type %in% c("both", "fixef")) {
        fixef.delete   <- nlme::fixef(model.delete)
        vcov.delete    <-  as.matrix(vcov(model.delete))
      }
      
      if(type %in% c("both", "varcomp")) {
        varcomp.delete <- varcomp.lme(model.delete)
        ranef.delete   <- nlme::ranef(model.delete)
        #if( length(flist) == 1 ) ranef.delete <- ranef.delete[[1]] #it fixes the issue, but might cause other problems 
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
      if(ncol(flist) == 1) {
        ranef.delete <- do.call('rbind', ranef.delete)
      }
      else {
        flist <- colnames(flist)
        temp  <- NULL
        for (i in 1:length(flist)) {
          temp[[i]] <- plyr::ldply(ranef.delete, function(x) x[[i]]) 
        }
        ranef.delete <- temp
        names(ranef.delete) <- names(nlme::ranef(model))
      }
    }
  }
  
  fixef.original <- nlme::fixef(model)
  ranef.original <- nlme::ranef(model)
  if(length(ranef.original) == 1) ranef.original <- ranef.original[[1]] 
  
  vcov.original <- as.matrix(vcov(model))
  varcomp.original <- varcomp.lme(model)
  
  val <- list(fixef.original = fixef.original, ranef.original = ranef.original,
              vcov.original = vcov.original, varcomp.original = varcomp.original,
              fixef.delete = fixef.delete, ranef.delete = ranef.delete,
              vcov.delete = vcov.delete, fitted.delete = fitted.delete,
              varcomp.delete = varcomp.delete)
  
  attr(val, "type") <- type
  class(val) <- "case_delete"
  return(val)
}



