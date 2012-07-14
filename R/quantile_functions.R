#' Calculating quantiles
#'
#' This function will summarize the quantiles that are found within the
#' data vector provided.
#'
#' @param x a numeric vector
#' @author Adam Loy \email{aloy@@istate.edu}
.sampleQuantiles <- function(x){
	i <- 1:length(x)
	p <- (i - .5)/length(x)
	return(p)
}

#' Constructing a normal quantile-quantile plot
#'
#' This function will construct a normal quantile-quantile plot within
#' the \code{ggplot} framework.
#'
#' @param x a numeric vector
#' @param line the method used to fit a reference line. If no reference line is desired,
#' leave the value as \code{NULL}. \code{line = "rlm"} will use robust regression to fit a
#' reference line. \code{line = "quantile"} will fit a line through the first and third quartiles.
#' @param ... other arguments to be passed to \code{qplot()}
#' @author Adam Loy \email{aloy@@istate.edu}
#' @export
ggplot_qqnorm <- function(x, line = NULL, ...){
	p <- .sampleQuantiles(x)
	theory <- qnorm(p = p)
	yp <- sort(x)
		ret <- qplot(x = theory, y = yp, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", ...)
	if(!is.null(line)){
		if(line == "quantile"){
			line.info <- qqlineInfo(x)
			ret <- ret + geom_abline(intercept = as.numeric(line.info[1]), slope = as.numeric(line.info[2]))
		}
		if(line == "rlm"){
			line.info <- coef(rlm(yp ~ theory))
			ret <- ret + geom_abline(intercept = as.numeric(line.info[1]), slope = as.numeric(line.info[2]))
		}
	}
	return(ret)
}

#' Constructing a normal probability plot
#'
#' This function will construct a normal probability plot within
#' the \code{ggplot} framework.
#'
#' @param x a numeric vector
#' @param line the method used to fit a reference line. If no reference line is desired, leave the
#' value as \code{NULL}. \code{line = "rlm"} will use robust regression to fit a reference line.
#' \code{line = "quantile"} will fit a line through the first and third quartiles.
#' @param ... other arguments to be passed to \code{qplot()}
#' @author Adam Loy \email{aloy@@istate.edu}
#' @export
ggplot_ppnorm <- function(x, line = NULL, ...){
	p <- .sampleQuantiles(x)
	yp <- sort(x)
	ret <- qplot(x = yp, y = qnorm(p), xlab = "x", ylab = "Theoretical Quantiles", ...)
	if(!is.null(line)){
		if(line == "quantile"){
			line.info <- pplineInfo(x)
			ret <- ret + geom_abline(intercept = as.numeric(line.info[1]), slope = as.numeric(line.info[2]))
		}
		if(line == "rlm"){
			line.info <- coef(rlm(qnorm(p) ~ yp))
			ret <- ret + geom_abline(intercept = as.numeric(line.info[1]), slope = as.numeric(line.info[2]))
		}
	}
	return(ret)
}

#' Adding a line to a normal probability plot
#'
#' This function will calculate the equation of the line that can be
#' used to help read the normal probability plot of interest.
#'
#' @param x a numeric vector
#' @author Adam Loy \email{aloy@@istate.edu}
pplineInfo <- function(x){
	yp <- quantile(x, c(0.25, 0.75))
	theory <- qnorm(p = c(0.25, 0.75))
	slope <- diff(theory)/diff(yp)
	intercept <- theory[1L] - slope * yp[1L]
	return(c(intercept = as.numeric(intercept), slope = as.numeric(slope)))
}

#' Adding simultaneous confidence bands to normal probability plots
#'
#' This function calculates simultaneous confidence bands for the estimated cdf using the
#' logistic-transform normal-approximation method. 
#'
#' Information on how to find the approximate factors \code{e} can be found in Meeker and Escobar (1998)
#' on page 61. They provide a table which summarizes some common choices for the factor.
#'
#' @param x a vector of probabilities
#' @param e factors for EP simultaneous approx. confidence bands (see Nair 1984 or Meeker and Escobar 1998)
#' @author Adam Loy \email{aloy@@iastate.edu}
#' @references Meeker, W. Q. and Escobar, L. A. (1998), \emph{Statistical Methods for Reliability Data},
#' New York, NY: Wiley.
#' @export
ppnormBand <- function(x, e){
	F <- x
	se <- sqrt(F * (1 - F) / length(F))
	w <- exp((e * se) / (F * (1 - F)))
	band <- list(lower = F / (F + (1 - F) * w), upper = F / (F + (1 - F) / w))
	return(band)
}

#' Adding a line to a normal quantile-quantile plot
#'
#' This function will calculate the equation of the line that can be
#' used to help read the normal quantile plot of interest.
#'
#' @param x a numeric vector
#' @author Adam Loy \email{aloy@@istate.edu}
qqlineInfo <- function(x){
	yp <- quantile(x, c(0.25, 0.75))
	theory <- qnorm(p = c(0.25, 0.75))
	slope <- diff(yp)/diff(theory)
	intercept <- yp[1L] - slope * theory[1L]
	return(c(intercept = as.numeric(intercept), slope = as.numeric(slope)))
}

#' Plotting numerous normal q-q plots on same plot
#'
#' This function will plot multiple normal q-q plots on the same plot. This
#' will be particulary useful when comparing the distrubtion between groups.
#' As differing slopes would indicate the normal distributions for the groups
#' do not share a common slope.
#'
#' @param x a numeric vector from which quantiles will be calculated
#' @param group a vector indicating group membership for each value in \code{x}.
#' @param line the method used to fit reference lines. If no reference lines are desired,
#' leave the value as \code{NULL}. \code{line = "rlm"} will use robust regression to fit
#' reference lines. \code{line = "quantile"} will fit lines through the first and third quartiles.
#' @param ... other arguments to be passed to ggplot
#' @param alpha_point alpha value specified for the points
#' @param alpha_line alpha value specified for the lines
#' @author Adam Loy \email{aloy@@istate.edu}
#' @references Hilden-Minton, J. A. (1995), ``Mulilevel Diagnostics for Mixed and Hierarchical Linear Models,''
#' Ph.D. thesis, University of California Los Angeles. 
#' @export
group_qqnorm <- function(x, group, line = NULL, alpha_point = 1, alpha_line = 1, ...){
	# Finding the slope and intercept for each group
	qq.list <- split(x, group)
	if(is.null(line) == FALSE){
		if(line == "quantile"){
			qq.coefs <- lapply(qq.list, qqlineInfo)
		}
		if(line == "rlm"){
			qq.coefs <- lapply(qq.list, function(j){
				p <- .sampleQuantiles(j)
				theory <- qnorm(p = p)
				yp <- sort(j)
				coef(rlm(yp ~ theory))
				})
		}
	qq.coefs <- do.call('rbind', qq.coefs)
	colnames(qq.coefs) <- c("intercept", "slope")
	qq.coefs <- data.frame(qq.coefs, group = row.names(qq.coefs))
	}
	
	# Defining the quantiles of interest for each group
	group.quant <- data.frame(x = x, group = group)
	group.quant <- ddply(group.quant, .(group), transform, p = .sampleQuantiles(x), yp = sort(x))
	group.quant <- ddply(group.quant, .(group), transform, theory = qnorm(p = p))
	#llply(qq.list, transform, p=sample.quantiles(df), yp= sort(df), theory = qnorm(p = p))
	
	# Plotting
	qq <- ggplot(data = group.quant, mapping = aes(x = theory, y = yp), ...) + geom_point(alpha = alpha_point)
	if(is.null(line) == FALSE){
		qq <- qq + geom_abline(aes(intercept = intercept, slope = slope), alpha = alpha_line, data = qq.coefs)
	}
	qq <- qq + xlab("Theoretical Quantiles") + ylab("Sample Quantiles")
	
	return(qq)
}
