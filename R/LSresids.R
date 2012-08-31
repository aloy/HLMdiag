#' Calculating case-level (level-1) least squares residuals
#'
#' This function calculates the case-level least squares residuals
#' found by fitting separate LS regression models to each case.
#'
#' @param formula a linear formula that is used by \code{adjust_lmList}' (y ~ x1 + ... + xn | g where g is a grouping factor)
#' @param lme.model  an object contatining the original hierarchical model fit using lmer()
#' @param semi.standardize if \code{TRUE} the semi-standardized residuals will also be returned
#' @return a data frame with the following columns: id, residual, fitted
#' @author Adam Loy \email{aloy@@istate.edu}
#' @examples
#' data(Oxboys, package = 'mlmRev')
#' fm <- lmer(formula = height ~ age + I(age^2) + (age + I(age^2)| Subject), data = Oxboys)
#' level1Resids <- LSresids(formula = height ~ age + I(age^2) | Subject, lme.model = fm, semi.standardize = TRUE)
#' 
#' \dontrun{wages.fm1 <- lmer(lnw ~ exper + (exper | id), data = wages)
#' LSresids(formula = lnw ~ exper | id, lme.model = wages.fm1)}
#' @export
LSresids <- function(object, level, sim = NULL, semi.standardize = TRUE){
	pform <- .parse_f(formula(object))
	response <- as.character(formula(object)[[2]])
	fixed <- subset(pform, type == "fixed" & included == TRUE)
	random <- subset(pform, type == "random" & included == TRUE)
	grp <- subset(pform, type == "group" & included == TRUE)
	
	data <- object@frame
	if(!is.null(sim)){data[,response] <- sim}
	
	if(level == 1){
		# fitting a separate LS regression model to each group
			
		form <- paste(paste(response, "~"), 
			paste(as.character(fixed[,"name"]), collapse = " + "), 
			" | ", as.character(grp[, "name"]))
	
		ls.models <- adjust_lmList(formula = formula(form), data = data)
		
#	if(level == 1){
		# obtaining the residuals, fitted values, and model frame
		# from that regression for each group
		ls.residuals <- lapply(ls.models, resid)
		ls.fitted <- lapply(ls.models, fitted)
	
		if(semi.standardize == TRUE){
			ls.influence <- lapply(ls.models, lm.influence)
			ls.hat <- lapply(ls.influence, function(x) x$hat)
		}
	
		# creating a data frame of the residuals, fitted values, and model frames
		ls.data <- lapply(ls.models, model.frame)
		
		row.order <- as.numeric(unlist(lapply(ls.data, function(x) row.names(x))))
		
		return.df <- data.frame(LS.resid = unlist(ls.residuals), fitted = unlist(ls.fitted), 
			hat = unlist(ls.hat))
	
		if(semi.standardize == TRUE){
			semi.std.resid  <- with(return.df, LS.resid / sqrt(1 - hat))
			semi.std.resid[is.infinite(semi.std.resid)] <- NaN
	
			return.df <- cbind(return.df, semi.std.resid = semi.std.resid)
		}
		
		return.df <- return.df[order(row.order),]
		return.df <- cbind(data, return.df)
		
		return(return.df)
	}
	
	if(level == 2){
		ranef_names <- names(ranef(object)[[1]])
		
		form <- paste(paste(formula(object)[[2]], "~"), 
			paste(paste(as.character(fixed[,"name"]), collapse = " + "), 
				paste(as.character(random[,"name"]), collapse = " + "),
				sep = " + "),
			" | ", as.character(grp[, "name"]))
		
		gform <- paste(response, "~", deparse(formula(form)[[3]][[2]]))
		
		ls.models <- adjust_lmList(formula = formula(form), data = data)
		global.model <- lm(formula = formula(gform), data = data)
		
		ls.resid <- coef(ls.models)[,ranef_names] - coef(global.model)[ranef_names]
		ls.resid <- as.data.frame(ls.resid)
		colnames(ls.resid) <- ranef_names
		
#		ls.resid <- Matrix(unlist(ls.residuals))
#		raw_resid <- tcrossprod(crossprod(solve(tcrossprod(object@Zt)), object@Zt), t(ls.resid))
		
#		ls.resid2 <- .format_rr(rr = raw_resid, object = object)

		return(ls.resid)
	}
}

#.format_rr <- function(rr, object){
#	# Formatting LS level-2 residuals
#	cn <- lapply(object@ST, colnames)
#	nt <- 1 # there is only 1 grouping factor in a two-level model
#	lterm <- lapply(lme4:::plist(lme4:::reinds(object@Gp), cn),
#					function(el){
#						cni <- el[[2]]
#						matrix(rr[ el[[1]] ], nc = length(cni),
#								dimnames = list(NULL, cni))
#					})
#	return(lterm[[1]])
#}