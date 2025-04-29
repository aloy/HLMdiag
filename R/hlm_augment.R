#' Calculating residuals and influence diagnostics for HLMs
#'
#' 
#' This function is used to compute residuals, fitted values, and influence diagnostics for a 
#' hierarchical linear model. The residuals and fitted values are computed using Least Squares(LS)
#' and Empirical Bayes (EB) methods. The influence diagnostics are computed through one step 
#' approximations. 
#'
#' @export
#' @param object an object of class \code{lmerMod} or \code{lme}.
#' @param level which residuals should be extracted and what cases should be deleted for influence diagnostics.
#'If \code{level = 1} (default), then within-group (case-level) residuals are returned and influence diagnostics
#'are calculated for individual observations. Otherwise, \code{level} should be the name of a grouping
#'factor as defined in \code{flist} for a \code{lmerMod} object or as in \code{groups} for a \code{lme} object.
#'This will return between-group residuals and influence diagnostics calculated for each group. 
#' @param include.ls a logical indicating if LS residuals should be included in the
#'return tibble. \code{include.ls = FALSE} decreases runtime substantially.
#' @param data the original data frame passed to `lmer`. This is only necessary for `lmerMod` models where
#'`na.action = "na.exclude"`
#' @param ... currently not used
#' @details The \code{hlm_augment} function combines functionality from \code{hlm_resid}
#'and \code{hlm_influence} for a simpler way of obtaining residuals and influence 
#'diagnostics. Please see \code{?hlm_resid} and \code{?hlm_influence} for additional information 
#'about the returned values.
#'@note \code{hlm_augment} does not allow for the deletion of specific cases, the specification of other
#'types of leverage, or the use of full refits of the model instead of one step approximations for influence
#'diagnostics. If this additional functionality is desired, \code{hlm_influence} should be used instead. The additional
#'parameter \code{standardize} is available in \code{hlm_resid}; if this are desired, \code{hlm_resid}
#'should be used instead. 
hlm_augment <- function(object, ...){
  UseMethod("hlm_augment", object)
}

#' @export
#' @rdname hlm_augment
#' @method hlm_augment default
hlm_augment.default <- function(object, ...){
  stop(paste("there is no hlm_augment() method for objects of class",
             paste(class(object), collapse=", ")))
}



#' @export
#' @method hlm_augment lmerMod
#' @rdname hlm_augment
hlm_augment.lmerMod <- function(object, level = 1, include.ls = TRUE, data = NULL, ...) {
  residuals <- hlm_resid(object, level = level, include.ls = include.ls, data = data)
  infl <- hlm_influence(object, level = level, data = data)
  if (level == 1) {
    infl <- infl[,-c(1:(1+ncol(object@frame)))]
  }
  else {
    infl <- infl[,-1]
  }
  aug.tibble <- tibble::add_column(residuals, infl)
  return(aug.tibble)
}

#' @export
#' @rdname hlm_augment
#' @method hlm_augment lme
#' @aliases hlm_augment
hlm_augment.lme <- function(object, level = 1, include.ls = TRUE, ...) {
  residuals <- hlm_resid(object, level = level, include.ls = include.ls)
  infl <- hlm_influence(object, level = level)
  
  #getting correct model frame, without extra variables
  fixed <- as.character(formula(object))
  dataform <- paste(fixed[2], fixed[1], fixed[3], " + ", 
                    paste(names(object$groups), collapse = " + "))
  data <- object$data %>%
    dplyr::mutate(across(where(is.character), ~ as.factor(.x))) %>%
    as.data.frame()
  newdata <- model.frame(formula(dataform), data)  
  
  if (level == 1) {
    infl <- infl[,-c(1:(1+ncol(newdata)))]
  }
  else {
    infl <- infl[,-1]
  }
  aug.tibble <- tibble::add_column(residuals, infl)
  return(aug.tibble)
}
