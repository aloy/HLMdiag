#' Calculating residuals from two-level HLMs
#'
#' \code{HLMresid} is a function that extracts residuals
#' from two-level mixed/hierarchical linear models fit
#' using \code{lmer}. That is, it is a function that
#' will extract residuals from \code{mer} objects
#' in a unified framework.
#' 
#' This function can extract residuals from either
#' level of the model, and can extract residuals
#' estimated using least squares (LS), Empirical 
#' Bayes (EB), or both. This unified framework
#' enables the analyst to more easily conduct
#' an upward residual analysis during model
#' exploration or checking.
#'
#' @param object an object of class \code{mer}.
#' @param level which residuals should be plotted: 1 for within-group
#' residuals or the name of the grouping factor (as defined in \code{flist} of the 
#' \code{mer} object) for between-group residuals, or \code{marginal}.
#' @param type how are the residuals predicted: either \code{EB} or \code{LS}. 
#' @param sim data optional argument giving the data frame used for LS residuals. This
#'  is used mainly for when dealing with simulations.
#' @param standardize if \code{TRUE} the standardized level-1
#' residuals will also be returned (if \code{level = 1}); if \code{"semi"} then
#' the semi-standardized level-1 residuals will be returned.
#' @author Adam Loy \email{aloy@@istate.edu}
#' @export
#' @keywords models regression
HLMresid <- function(object, level, type = "EB", sim = NULL, standardize = TRUE){
  if(!is(object, "mer")) stop("object must be of class 'mer'")
  if(!level %in% c(1, names(object@flist), "marginal")) {
    stop("level can only be 1, a grouping factor from the fitted model,
         or marginal.")
  }
  if(!type %in% c("EB", "LS")) stop("type must be either 'EB' or 'LS'.")
  if(!is.null(standardize) && !standardize %in% c(TRUE, "semi")) {
    stop("standardize can only be specified to be TRUE or 'semi'.")
  }
	
	if(level == "marginal"){
		return(object@y - object@X %*% fixef(object))
	}
	
	if(level == 1){
		if(type == "LS"){
			return(LSresids(object = object, level = level, sim = sim, standardize = standardize))
		}
		if(type == "EB"){
			return(resid(object))
		}
	}
	
	if(level %in% names(object@flist)){
		if(type != "EB"){
			return(LSresids(object = object, level = level, sim = sim, standardize = standardize))
		}
		if(type != "LS"){
			return(ranef(object)[[level]])
		}
	}
}