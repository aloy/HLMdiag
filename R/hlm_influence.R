#'@export
hlm_influence <- function(model, ...) {
  UseMethod("hlm_influence", model)
}

#' @export
#' @rdname hlm_influence.lmerMod
#' @method hlm_influence default
hlm_influence.default <- function(model, ...){
  stop(paste("there is no hlm_influence() method for objects of class",
             paste(class(model), collapse=", ")))
}

#'Calculating influence diagnostics for HLMs
#'
#'@description
#'This function is used to compute influence diagnostics for a hierarchical linear model.
#'It takes a model fit as a \code{lmerMod} object or as a \code{lme} object and returns a tibble with Cook's
#'distance, MDFFITS, covtrace, covratio, and leverage.
#'
#'@export
#'@method hlm_influence lmerMod
#'@aliases hlm_influence
#'@param model an object of class \code{lmerMod} or \code{lme}
#'@param level used to define the group for which cases are deleted and influence
#'diagnostics are calculated. If \code{level = 1} (default), then influence diagnostics are
#'calculated for individual observations. Otherwise, \code{level} should be the name of a grouping
#'factor as defined in \code{flist} for a \code{lmerMod} object or as in \code{groups} for a \code{lme} object.
#'@param delete numeric index of individual cases to be deleted. If the \code{level} parameter 
#'is specified, \code{delete} may also take the form of a character vector consisting of group 
#'names as they appear in \code{flist} for \code{lme4} models or as in \code{groups} for \code{nlme} models. 
#'If \code{delete = NULL} then all cases are iteratively deleted.
#'@param approx logical parameter used to determine how the influence diagnostics are calculated.
#'If \code{FALSE} (default), influence diagnostics are calculated using a one step approximation.
#'If \code{TRUE}, influence diagnostics are calculated by iteratively deleting groups and refitting
#'the model using \code{lmer}. This method is more accurate, but slower than the one step approximation.
#'If \code{approx = FALSE}, the returned tibble also contains columns for relative variance change (RVC).
#'@param leverage a character vector to determine which types of leverage should be included in the
#'returned tibble. There are four options: 'overall' (default), 'fixef', 'ranef', or 'ranef.uc'.
#'One or more types may be specified. For additional information about the types of leverage, see
#'\code{?leverage}.
#'@param data (optional) the data frame used to fit the model. This is only necessary for \code{lmerMod} models if
#'\code{na.action = "na.exclude"} was set. 
#'@param ... not in use
#'
#'@details
#'The \code{hlm_influence} function provides a wrapper that appends influence diagnostics
#'to the original data. The approximated influence diagnostics returned by this
#'function are equivalent to those returned by \code{cooks.distance}, \code{mdffits}, \code{covtrace},
#'\code{covratio}, and \code{leverage}. The exact influence diagnostics obtained through a full
#'refit of the data are also available through \code{case_delete} and the accompanying functions
#'\code{cooks.distance}, \code{mdffits}, \code{covtrace}, and \code{covratio} that can be called
#'directly on the \code{case_delete} object.
#'@note 
#'It is possible to set \code{level} and delete individual cases from different groups using 
#'\code{delete}, so numeric indices should be double checked to confirm that they encompass entire groups.
#'Additionally, if \code{delete} is specified, leverage values are not returned in the resulting tibble. 
hlm_influence.lmerMod <- function(model, level = 1, delete = NULL, approx = TRUE, leverage = "overall", data = NULL, ...) {
  
  if (!level %in% names(model@flist) & level != 1) {
    stop(paste(level, "is not a valid level for this model"))
  }
  
  for (i in 1:length(leverage)) {
    if (!leverage[i] %in% c("overall", "fixef", "ranef", "ranef.uc")) {
      stop(paste(leverage[i], "is not a valid option for a type of leverage. Valid options are limited to: 'overall', 'fixef', 'ranef', or 'ranef.uc'."))
    }
  }
  
  if (hasArg(group)) {
    group <- NULL
    warning("group is not a valid argument for this function. As of version 0.4.0, group has been replaced by level. See ?hlm_influence for more information.")
  }
  
  na.action <- attr(model@frame, "na.action")
  if(class(na.action) == "exclude" & is.null(data)) {
    stop("Please provide the data frame used to fit the model. This is necessary when the na.action is set to na.exclude.")
  }
  
  if(!isNestedModel(model) & approx == FALSE) {
    stop("Full refits of the model are currently not implemented for models with crossed random effects.")
  }
  
  if(!is.null(delete)) { 
    warning("If the delete argument is specified, leverage cannot be returned. See ?hlm_influence for more information.")
  }
  
  if (approx) { #one step approximations
    infl.tbl <- tibble::tibble(cooksd = as.vector(cooks.distance(model, level = level, delete = delete)),
                               mdffits = as.vector(mdffits(model, level = level, delete = delete)),
                               covtrace = covtrace(model, level = level, delete = delete),
                               covratio = as.numeric(covratio(model, level = level, delete = delete)))
    
    if(!is.null(delete)) {
      return(infl.tbl)
    }
    
    if (isNestedModel(model)) {
      leverage.df <- as.data.frame(leverage(model, level = level)[,leverage])
      colnames(leverage.df) <- purrr::map_chr(leverage, function(s) stringr::str_c("leverage", s, sep = "."))
      infl.tbl <- tibble::add_column(infl.tbl, leverage.df)
    }
    else{
      warning("Leverage is currently not implemented for models with crossed random effects.")
    }
    
    if (level == 1) {
      infl.tbl <- tibble::add_column(infl.tbl, model@frame, .before = 1)
      if (class(na.action) == "exclude") {
        infl.tbl <- .lmerMod_add_NArows(model, infl.tbl, na.action, data)
      }
      infl.tbl <- tibble::add_column(infl.tbl, id = 1:nrow(infl.tbl), .before = 1)
    }
    else {
      infl.tbl <- tibble::add_column(infl.tbl, unique(model@flist[[level]]), .before = 1)
      names(infl.tbl)[1] <- level
    }
  }
  else { #full refits
    case <- case_delete(model, level = level, delete = delete)
    
    infl.tbl <- tibble::tibble(cooksd = as.vector(cooks.distance(case)),
                               mdffits = as.vector(mdffits(case)),
                               covtrace = covtrace(case),
                               covratio = covratio(case))
    
    if (!is.null(delete)) {
      rvc.df <- as.data.frame(t(rvc(case))) #need to take transpose of rvc output when delete isn't null
      colnames(rvc.df) <- purrr::map_chr(names(rvc.df), function(s) stringr::str_c("rvc", s, sep = "."))
      infl.tbl <- tibble::add_column(infl.tbl, rvc.df)
      return(infl.tbl)
    }
    else{
      rvc.df <- as.data.frame(rvc(case))
      colnames(rvc.df) <- purrr::map_chr(names(rvc.df), function(s) stringr::str_c("rvc", s, sep = "."))
      infl.tbl <- tibble::add_column(infl.tbl, rvc.df)
      
      leverage.df <- as.data.frame(leverage(model, level = level)[,leverage])
      colnames(leverage.df) <- purrr::map_chr(leverage, function(s) stringr::str_c("leverage", s, sep = "."))
      infl.tbl <- tibble::add_column(infl.tbl, leverage.df)
    }
    
    if (level == 1) {
      infl.tbl <- tibble::add_column(infl.tbl, model@frame, .before = 1)  
      if (class(na.action) == "exclude") {
        infl.tbl <- .lmerMod_add_NArows(model, infl.tbl, na.action, data)
      }
      infl.tbl <- tibble::add_column(infl.tbl, id = 1:nrow(infl.tbl), .before = 1)
    }
    else {
      infl.tbl <- tibble::add_column(infl.tbl, Group = unique(model@flist[[level]]), .before = 1)
      names(infl.tbl)[1] <- level
    }
  }
  return(infl.tbl)
}

#' @export
#' @rdname hlm_influence.lmerMod
#' @method hlm_influence lme
#' @aliases hlm_influence
hlm_influence.lme <- function(model, level = 1, delete = NULL, approx = TRUE, leverage = "overall", ...) {
  
  if (!level %in% names(model$groups) & level != 1) {
    stop(paste(level, "is not a valid level for this model"))
  }
  
  for (i in 1:length(leverage)) {
    if (!leverage[i] %in% c("overall", "fixef", "ranef", "ranef.uc")) {
      stop(paste(leverage[i], "is not a valid option for a type of leverage. Valid options are limited to: 'overall', 'fixef', 'ranef', or 'ranef.uc'."))
    }
  }
  
  if(!isNestedModel(model) & approx == FALSE) {
    stop("Full refits of the model are currently not implemented for models with crossed random effects.")
  }
  
  if (hasArg(group)) {
    group <- NULL
    warning("group is not a valid argument for this function. As of version 0.4.0, group has been replaced by level. See ?hlm_influence for more information.")
  }
  
  if(!is.null(delete)) { 
    warning("If the delete argument is specified, leverage cannot be returned. See ?hlm_influence for more information.")
  }

  
  na.action <- model$na.action
  
  if (approx) { #one step approximations
    infl.tbl <- tibble::tibble(cooksd = as.vector(cooks.distance(model, level = level, delete = delete)),
                               mdffits = as.vector(mdffits(model, level = level, delete = delete)),
                               covtrace = covtrace(model, level = level, delete = delete),
                               covratio = as.numeric(covratio(model, level = level, delete = delete)))
    
    if(!is.null(delete)) {
      return(infl.tbl)
    }
    
    if (isNestedModel(model)) {
      leverage.df <- as.data.frame(leverage(model, level = level)[,leverage])
      colnames(leverage.df) <- purrr::map_chr(leverage, function(s) stringr::str_c("leverage", s, sep = "."))
      infl.tbl <- tibble::add_column(infl.tbl, leverage.df)
    }
    else{
      warning("Leverage is currently not implemented for models with crossed random effects.")
    }
    
    if (level == 1) {
      
      fixed <- formula(model)
      dataform <- paste(fixed[2], "~", fixed[3], " + ",
                        paste(names(model$groups), collapse = " + ")) 
      data.fixed <- model$data %>%
        dplyr::mutate(across(where(is.character), ~ as.factor(.x))) %>%
        as.data.frame()
      new.data <- model.frame(formula(dataform), data.fixed)
      
      infl.tbl <- tibble::add_column(infl.tbl, new.data, .before = 1) 
      
      if (class(na.action) == "exclude") {
        infl.tbl <- .lme_add_NArows(model, infl.tbl, na.action, org.data = model$data, fixed.data = new.data)
      }
      infl.tbl <- tibble::add_column(infl.tbl, id = 1:nrow(infl.tbl), .before = 1)
    }
    else {
      infl.tbl <- tibble::add_column(infl.tbl, unique(model$groups[[level]]), .before = 1)
      names(infl.tbl)[1] <- level
    }
  }
  else { #full refits
    case <- case_delete(model, level = level, delete = delete)
    
    infl.tbl <- tibble::tibble(cooksd = as.vector(cooks.distance(case)),
                               mdffits = as.vector(mdffits(case)),
                               covtrace = covtrace(case),
                               covratio = covratio(case))
    
    if (!is.null(delete)) {
      rvc.df <- as.data.frame(t(rvc(case))) #need to take transpose of rvc output when delete isn't null
      colnames(rvc.df) <- purrr::map_chr(names(rvc.df), function(s) stringr::str_c("rvc", s, sep = "."))
      infl.tbl <- tibble::add_column(infl.tbl, rvc.df)
      return(infl.tbl)
    }
    else{
      rvc.df <- as.data.frame(rvc(case))
      colnames(rvc.df) <- purrr::map_chr(names(rvc.df), function(s) stringr::str_c("rvc", s, sep = "."))
      infl.tbl <- tibble::add_column(infl.tbl, rvc.df)
      
      leverage.df <- as.data.frame(leverage(model, level = level)[,leverage])
      colnames(leverage.df) <- purrr::map_chr(leverage, function(s) stringr::str_c("leverage", s, sep = "."))
      infl.tbl <- tibble::add_column(infl.tbl, leverage.df)
    }
    
    if (level == 1) {
      fixed <- formula(model)
      dataform <- paste(fixed[2], "~", fixed[3], " + ",
                        paste(names(model$groups), collapse = " + ")) 
      data.fixed <- model$data %>%
        dplyr::mutate(across(where(is.character), ~ as.factor(.x))) %>%
        as.data.frame()
      new.data <- model.frame(formula(dataform), data.fixed)
      infl.tbl <- tibble::add_column(infl.tbl, new.data, .before = 1)  
      if (class(na.action) == "exclude") {
        infl.tbl <- .lme_add_NArows(model, infl.tbl, na.action, org.data = model$data, fixed.data = new.data)
      }
      infl.tbl <- tibble::add_column(infl.tbl, id = 1:nrow(infl.tbl), .before = 1)
    }
    else {
      infl.tbl <- tibble::add_column(infl.tbl, Group = unique(model$groups[[level]]), .before = 1)
      names(infl.tbl)[1] <- level
    }
  }
  return(infl.tbl)
}
