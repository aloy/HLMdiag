#'@export
hlm_influence <- function(model, ...) {
  UseMethod("hlm_influence", model)
}

#' @export
#' @rdname hlm_influence.lmerMod
#' @method hlm_influence default
#' @S3method hlm_influence default
hlm_influence.default <- function(model, ...){
  stop(paste("there is no hlm_influence() method for objects of class",
             paste(class(model), collapse=", ")))
}

#'Calculating influence diagnostics for HLMs
#'
#'@description
#'This function is used to compute influence diagnostics for a hierarchical linear model. 
#'It takes a model fit as a \code{lmerMod} object and returns a tibble with Cook's 
#'distance, MDFFITS, covtrace, covratio, and leverage. 
#'
#'@export
#'@method hlm_influence lmerMod
#'@S3method hlm_influence lmerMod
#'@aliases hlm_influence
#'@param model an object of class \code{lmerMod}
#'@param level a variable used to define the group for which cases are deleted and influence 
#'diagnostics are calculated. If \code{level} equals 1 (default), then influence diagnostics are 
#'calculated for individual observations. Otherwise, \code{level} should be the name of a grouping
#'factor as defined in \code{flist} of the \code{lmerMod} object. 
#'@param approx logical parameter used to determine how the influence diagnostics are calculated.
#'If \code{FALSE} (default), influence diagnostics are calculated using a one step approximation. 
#'If \code{TRUE}, influence diagnostics are caclulated by iteratively deleting groups and refitting
#'the model using \code{lmer}. This method is more accurate, but slower than the one step approximation. 
#'@param leverage a character vector to determine which types of leverage should be included in the 
#'returned tibble. There are five options: 'overall' (default), 'fixef', 'ranef', 'ranef.uc', or 'all', 
#'which will return all four types. One or more types may be specified. For additional information 
#'about the types of leverage, see \code{?leverage}. 
#'
#'@details 
#'The \code{hlm_influence} function provides a wrapper that appends influence diagnostics 
#'to the original data. The approximated influence diagnostics returned by this 
#'function are equivalent to those returned by \code{cooks.distance}, \code{mdffits}, \code{covtrace}, 
#'\code{covratio}, and \code{leverage}. The exact influence diagnostics obtained through a full 
#'refit of the data are also avaliable through \code{case_delete} and the accompaning functions 
#'\code{cooks.distance}, \code{mdffits}, \code{covtrace}, and \code{covratio} that can be called 
#'directly on the \code{case_delete} object. 


hlm_influence.lmerMod <- function(model, level = 1, approx = TRUE, leverage = "overall", ...) {
  if (hasArg(group)) {
    warning("group is not a valid argument for this function. As of version 0.4.0, group has been replaced by level. See ?hlm_influence for more information.")
  }
  
  for (i in 1:length(leverage)) {
    if (!leverage[i] %in% c("overall", "fixef", "ranef", "ranef.uc", "all")) {
      stop(paste(leverage[i], "is not a valid option for a type of leverage. Valid options are limited to: 'overall', 'fixef', 'ranef', 'ranef.uc', or 'all'."))
    }
  }
  
  if (leverage == "all") {
    leverage <- c("overall", "fixef", "ranef", "ranef.uc")
  }
  
  if (!level %in% names(model@flist) & level != 1) {
    stop(paste(level, "is not a valid level for this model"))
  }
  
  if (approx) { #approximations
    infl.tbl <- tibble::tibble(cooksd = as.vector(cooks.distance(model, level = level)),
                    mdffits = as.vector(mdffits(model, level = level)),
                    covtrace = covtrace(model, level = level),
                    covratio = covratio(model, level = level))
        
    leverage.df <- leverage(model, level = level)[,leverage]
    names(leverage.df) <- purrr::map_chr(names(leverage.df), function(s) stringr::str_c("leverage", s, sep = "."))
    infl.tbl <- tibble::add_column(infl.tbl, leverage.df)
        
    if (level == 1) {
      infl.tbl <- tibble::add_column(infl.tbl, model@frame, .before = 1)
    }
    else {
      infl.tbl <- tibble::add_column(infl.tbl, unique(model@flist[[level]]), .before = 1)
      names(infl.tbl)[1] <- level 
    }
  }
  else { #full refits 
    case <- case_delete(model, level = level)
  
    infl.tbl <- tibble::tibble(cooksd = as.vector(cooks.distance(case)),
                            mdffits = as.vector(mdffits(case)),
                            covtrace = covtrace(case),
                            covratio = covratio(case))
    leverage.df <- leverage(model, level = level)[,leverage]
    names(leverage.df) <- purrr::map_chr(names(leverage.df), function(s) stringr::str_c("leverage", s, sep = "."))
    infl.tbl <- tibble::add_column(infl.tbl, leverage.df)
    
    
    if (level == 1) {
      infl.tbl <- tibble::add_column(infl.tbl, model@frame, .before = 1)
    }
    else {
      infl.tbl <- tibble::add_column(infl.tbl, Group = unique(model@flist[[level]]), .before = 1)
      names(infl.tbl)[1] <- level
    }
  }
  return(infl.tbl)
}

