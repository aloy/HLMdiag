#'Extracting residuals from two-level HLMs
#'
#'This is a wrapper function to extract the residuals from a two-level normal
#'hierarchical linear model (HLM) fit using \code{lmer}. The EB residuals are
#'extracted using using the functions \code{resid} and \code{ranef} for levels
#'1 and 2, respectively.  To extract the LS residuals, \code{LSresids} is used.
#'
#'
#'@param object a fitted model object of class \code{mer}.
#'@param level the level from which the residuals should be extracted; either
#'\code{1} or \code{2} (both can be specified).
#'@param type the type of residual to be extracted; \code{"LS"}, \code{"EB"},
#'\code{"both"} (i.e. LS and EB), or \code{"marginal"}.
#'@param sim an optional data vector containing a (simulated) response
#'variable.
#'@param semi.standardize if \code{TRUE} the semi-standardized residuals will
#'also be returned. This argument is only used when \code{level = 1}.
#'@return If \code{type = "marginal"} is specified, a numeric vector of the
#'marginal residuals is returned.  Otherwise, a data frame is returned
#'containing the following: \itemize{ \item When \code{level=1}: the model
#'frame, the residuals (\code{LS.resid} and/or \code{EB.resid}), the fitted
#'values (\code{fitted}), and, if \code{semi.standardize = TRUE}, the diagonal
#'elements of the hat matrix (\code{hat}), the semi-standardized residuals
#'(\code{semi.std.resid})
#'
#'\item When \code{level=2}: the LS and/or EB residuals corresponding to each
#'random effect.
#'
#'\item When \code{level = c(1, 2)}: a list with elements \code{level.1} and
#'\code{level.2}. Each element is specified as stated above.  }
#'@author Adam Loy \email{aloy@@istate.edu}
#'@examples
#'
#'data(Oxboys, package = 'mlmRev')
#'fm <- lmer(formula = height ~ age + I(age^2) + (age + I(age^2)| Subject), data = Oxboys)
#'level1Resids <- HLMresid(object = fm, level = c(1,2), type = "both", semi.standardize = TRUE)
#'
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

