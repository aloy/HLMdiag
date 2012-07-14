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
#' @param level specification of the residual of interest.
#' @param how the residuals should be estimated; default is \code{"all"}.
#' @param if \code{type = "all"} or \code{"LS"}, a linear formula
#'  used by \code{adjust_lmList} (y ~ x1 + ... + xn | g where g is a grouping factor)
#'  must be specified.
#' @param data optional argument giving the data frame used for LS residuals. This
#'  is used mainly for when dealing with simulations.
HLMresid <- function(object, level = c(1, 2), type = c("both", "LS", "EB", "marginal"), sim = NULL, semi.standardize = TRUE){
	type <- match.arg(type)
	if(is.null(data)){data <- object@frame}
	
	if(type == "marginal"){
		return(object@y - object@X %*% fixef(object))
	}
	
	if(1 %in% level){
		if(type != "EB"){
			LS1 <- LSresids(object = object, level = 1, sim = sim, semi.standardize = semi.standardize)
		}
		if(type != "LS"){
			EB1 <- resid(object)
		}
	}
	
	if(2 %in% level){
		if(type != "EB"){
			LS2 <- LSresids(object = object, level = 2, sim = sim, semi.standardize = semi.standardize)
		}
		if(type != "LS"){
			EB2 <- ranef(object)[[1]]
		}
	}
	
	level.1 <- level.2 <- NULL

	if(type == "EB"){
		if(1 %in% level){
			level.1 <- cbind(EB.resid = EB1)
		}
		if(2 %in% level){
			level.2 <- cbind(EB = EB2)
		}	
	}
	
	if(type == "LS"){
		if(1 %in% level){
			level.1 <- LS1
		}
		if(2 %in% level){
			level.2 <- cbind(LS = LS2)
		}
	}
	
	if(type == "both"){
		if(1 %in% level){
			level.1 <- cbind(LS1, EB.resid = EB1)
		}
		if(2 %in% level){
			level.2 <- cbind(LS2, EB2)
			nc <- dim(level.2)[2]/2
			colnames(level.2) <- paste(rep(c("LS", "EB"), each = nc), colnames(level.2), sep = ".")
		}
	}
	
	res <- list(level.1 = level.1, level.2 = level.2)
	
	if(sum(level %in% c(1,2)) == 2){ return(res) }
	else{return(res[[level]])}
}

