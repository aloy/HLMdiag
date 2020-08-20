#' @export
LSresids <- function(object, ...){
  UseMethod("LSresids", object)
}

#' @export
#' @rdname LSresids.mer
#' @method LSresids default
LSresids.default <- function(object, ...){
  stop(paste("there is no LSresids() method for objects of class",
             paste(class(object), collapse=", ")))
}

#' Calculating least squares residuals
#'
#' This function calculates least squares (LS) residuals
#' found by fitting separate LS regression models to each case.
#' For examples see the documentation for \code{HLMresid}.
#'
#' @export
#' @method LSresids mer
#' @aliases LSresids
#' @param object an object of class \code{mer} or \code{lmerMod}.
#' @param level which residuals should be extracted: 1 for case-level
#' residuals or the name of a grouping factor (as defined in \code{flist} of the 
#' \code{mer} object) for between-group residuals.
#' @param sim optional argument giving the data frame used for LS residuals. This
#'  is used mainly when dealing with simulations. Removed in version 0.3.2.
#' @param standardize if \code{TRUE} the standardized level-1
#' residuals will also be returned (if \code{level = 1}); if \code{"semi"} then
#' the semi-standardized level-1 residuals will be returned.
#' @param ... do not use
#' @author Adam Loy \email{loyad01@@gmail.com}
#' @references 
#' Hilden-Minton, J. (1995) Multilevel diagnostics for mixed and hierarchical 
#' linear models. University of California Los Angeles.
#' @export
#' @seealso \code{\link{HLMresid}}
#' @keywords models regression
LSresids.mer <- function(object, level, sim = NULL, standardize = FALSE, ...){
  .Deprecated(new = 'LSresids.lmerMod')
  if(!object@dims["LMM"]){
    stop("LSresids is currently not implemented for GLMMs or NLMMs.")
  }
  if(!level %in% c(1, names(object@flist))) {
    stop("level can only be 1 or a grouping factor from the fitted model.")
  }
  if(object@dims[["nest"]] == 0) {
    stop("LSresids has not yet been implemented for models with 
          crossed random effects")
  }
  if(!is.null(standardize) && !standardize %in% c(FALSE, TRUE, "semi")) {
    stop("standardize can only be specified to be logical or 'semi' .")
  }
  
  LS.resid <- NULL # Make codetools happy
  
  fixed <- as.character( fixform( formula(object) ) )
	
	data <- object@frame
	if(!is.null(sim)){data[,fixed[2]] <- sim}
	
	if(level == 1){
		# fitting a separate LS regression model to each group
		form <- paste(fixed[2], fixed[1], fixed[3], "|", names(object@flist)[1])
	
		ls.models <- adjust_lmList(object = formula(form), data = data)
		
		ls.residuals <- lapply(ls.models, resid)
		ls.fitted <- lapply(ls.models, fitted)
	
		# creating a data frame of the residuals, fitted values, and model frames
		ls.data <- lapply(ls.models, model.frame)
		res.data <- do.call('rbind', ls.data)
		
		row.order <- unlist(lapply(ls.data, function(x) row.names(x)))

# 		row.order <- as.numeric(unlist(lapply(ls.data, function(x) row.names(x))))
		
		return.df <- data.frame(LS.resid = unlist(ls.residuals), 
                            fitted = unlist(ls.fitted))
	
		if(!is.null(standardize) && standardize == "semi"){
		  ls.influence <- lapply(ls.models, lm.influence)
		  ls.hat <- lapply(ls.influence, function(x) x$hat)
      
		  h <- unlist(ls.hat)
			semi.std.resid  <- with(return.df, LS.resid / sqrt(1 - h))
			semi.std.resid[is.infinite(semi.std.resid)] <- NA
	
			return.df <- cbind(return.df, semi.std.resid = semi.std.resid)
		}
    
		if(!is.null(standardize) && standardize == TRUE){
		  ls.rstandard <- unlist(lapply(ls.models, rstandard))
		  ls.rstandard[is.infinite(ls.rstandard)] <- NA
      
      return.df <- cbind(return.df, std.resid = ls.rstandard)
		}
		
# 		return.df <- return.df[order(row.order),]
		return.df <- cbind(data, return.df)
    rownames(return.df) <- row.order

		return(return.df)
	}
	
	if(level != 1){
    n.ranefs <- length(names(object@flist))
		ranef_names <- names( lme4::ranef(object)[[level]] )
		
		form <- paste(fixed[2], fixed[1], fixed[3], "|", level)
    
    if( level == names(object@flist)[n.ranefs] ) { # For highest level unit
      gform <- fixform( formula(object) )
      global.model <- lm(formula = formula(gform), data = data)
      
      ls.models <- adjust_lmList(object = formula(form), data = data)
      ls.resid <- coef(ls.models)[,ranef_names] - coef(global.model)[ranef_names]
      
    } else{ # For 'intermediate' level unit
      higher.level <- names(object@flist)[which(names(object@flist) == level) + 1]
      gform <- paste(fixed[2], fixed[1], fixed[3], "|", higher.level)
      global.model <- adjust_lmList(object = formula(gform), data = data)
      
      ls.models <- adjust_lmList(object = formula(form), data = data)
      
      # matching lower level units to higer level units
      units.mat <- unique(object@flist)
      n.higher.level <- table(units.mat[higher.level])
      
      global.coefs <- rep(unlist(coef(global.model)[ranef_names]), 
                          times = n.higher.level)
      ls.resid <- coef(ls.models)[,ranef_names] - global.coefs
      
    }
		
		ls.resid <- as.data.frame(ls.resid)
		colnames(ls.resid) <- ranef_names

		return(ls.resid)
	}
}

# 'fixform' is a copy of the 'nobars' function in the lme4 package, 
# renamed so it doesn't cause any conflicts.  This should not be exported.
fixform <- function (term) 
{
  if (!("|" %in% all.names(term))) 
    return(term)
  if (is.call(term) && term[[1]] == as.name("|")) 
    return(NULL)
  if (length(term) == 2) {
    nb <- fixform(term[[2]])
    if (is.null(nb)) 
      return(NULL)
    term[[2]] <- nb
    return(term)
  }
  nb2 <- fixform(term[[2]])
  nb3 <- fixform(term[[3]])
  if (is.null(nb2)) 
    return(nb3)
  if (is.null(nb3)) 
    return(nb2)
  term[[2]] <- nb2
  term[[3]] <- nb3
  term
}



#' @export
#' @rdname LSresids.mer
#' @method LSresids lmerMod
LSresids.lmerMod <- function(object, level, standardize = FALSE, ...){
  if(!lme4::isLMM(object)){
    stop("LSresids is currently not implemented for GLMMs or NLMMs.")
  }
  if(!level %in% c(1, names(object@flist))) {
    stop("level can only be 1 or a grouping factor from the fitted model.")
  }
  if(!isNestedModel(object)) {
    stop("LSresids has not yet been implemented for models with 
         crossed random effects")
  }
  if(!is.null(standardize) && !standardize %in% c(FALSE, TRUE, "semi")) {
    stop("standardize can only be specified to be logical or 'semi' .")
  }

  if(level == 1){
      y <- lme4::getME(object, "y")
      X <- lme4::getME(object, "X")
      g <- as.data.frame(object@flist[1])
      names(g) <- "group"
      
      frame <- cbind(y = y, X, g)                  #bind y, Xs, group
      frame <- frame[,-2]                          #remove intercept column
      
      frame.split <- split(frame, frame[,"group"]) #split on groups
      frame.lgl <- 
        purrr::map(frame.split, function(split) {  #check if column elements are equal
          purrr::map_lgl(split, ~all(.x == .x[1]))
        })
      
      for (i in 1:length(frame.split)) {           #remove the all equal columns
        frame.lgl[[i]][1] <- FALSE                 #make sure we never remove response
        frame.split[[i]] <- frame.split[[i]][!frame.lgl[[i]]]
      }
      
      ls.models <-                                 #fit the lm
        purrr::map(frame.split, function(split){
          lm(y ~ ., data = split)
        })
      
      ls.rank <- purrr::map_dbl(ls.models, ~qr(.x)$rank) #check rank deficiency
      deficient.groups <- names(which(ls.rank != max(ls.rank)))
      if(length(deficient.groups) != 0) {
        warning("LS residuals might be inaccurate as one or more groups are rank deficient.", 
                "\nUse the 'include.ls = FALSE' parameter to get EB residuals only.")
      }
      
      ls.residuals <- purrr::map(ls.models, resid) #calculate residuals within group
      ls.fitted <- purrr::map(ls.models, fitted)   #calculate fitted values
      
      row.order <- unlist(purrr::map(frame.split, function(x) row.names(x)))
      return.df <- data.frame(.ls.resid = unlist(ls.residuals),
                              .ls.fitted = unlist(ls.fitted))
      
      if(!is.null(standardize) && standardize == "semi"){
        ls.hat <- unlist(purrr::map(ls.models, ~lm.influence(.x)$hat))
        semi.std.resid  <- unlist(ls.residuals) / sqrt(1 - ls.hat)
        semi.std.resid[is.infinite(semi.std.resid)] <- NA
        
        return.df <- data.frame(.semi.ls.resid = semi.std.resid, 
                                .ls.fitted = unlist(ls.fitted))
      }

      if(!is.null(standardize) && standardize == TRUE){
        ls.rstandard <- unlist(purrr::map(ls.models, rstandard))
        ls.rstandard[is.infinite(ls.rstandard)] <- NA

        return.df <- data.frame(.std.ls.resid = ls.rstandard, 
                                .ls.fitted = unlist(ls.fitted))
      }
      
      rownames(return.df) <- row.order
      return(return.df)
  }
  
  if(level %in% names(object@flist)){
    fixed <- as.character(lme4::nobars( formula(object)))
    data <- object@frame
    names(data)[which(names(data) == fixed[2])] <- "y"
    fixed[2] <- "y"
    
    n.ranefs <- length(names(object@flist))
    ranef_names <- names( lme4::ranef(object)[[level]] )
    
    if(level == names(object@flist)[n.ranefs]) { #highest level
      form <- paste(fixed[2], fixed[1], fixed[3], "|", level)
      
      ls.models <- suppressWarnings(lme4::lmList(formula = formula(form), data = data))
      if (sum(purrr::map_lgl(ls.models, is.null)) == length(ls.models)) {
        stop("The model matrix is rank deficient. LS residuals cannot be calculated.", 
             "\nUse the 'include.ls = FALSE' parameter to get EB residuals only.")
      }
      if (!is.null(attr(ls.models, which = "warningMsg"))) {
        warning("The model matrix is likely rank deficient. Some LS residuals cannot be calculated.",
                "\nIt is recommended to use EB (.ranef) group level residuals for this model.")
      }
      if(standardize == "semi") standardize <- FALSE
      ls.ranef <- lme4::ranef(ls.models, standard = standardize)[ranef_names]
      ls.resid <- purrr::map_dfc(ls.ranef, ~.x)
      ls.resid <- tibble::tibble(group = row.names(coef(ls.models)),
                     ls.resid)
      
      if(ncol(ls.resid) != length(ranef_names) + 1) {
        warning("The model contains a random effect term for variables with no fixed effect term.",
                "\nSome LS group level residuals cannot be calculated.")
      }
      
      return(ls.resid)
      
    } else { #middle level
      # group by (specified level + 1)
      higher.level <- names(object@flist)[which(names(object@flist) == level) + 1]
      split_data <- split(data, data[,higher.level])
      # remove any empty groups
      split_data <- split_data[which(purrr::map_lgl(split_data, ~(nrow(.x)) != 0))]
      
      # for each group use ranef as above
      vars <- stringr::str_split(level, ":")[[1]]
      g <- vars[which(!vars %in% names(object@flist))]
      form <- paste(fixed[2], fixed[1], fixed[3], "|", g) 
      
      # recombine data frame
      ls.models <- suppressWarnings(purrr::map(split_data,
                        ~lme4::lmList(formula = formula(form), data = .x)))
      ls.warnings <- purrr::map_lgl(ls.models, ~!is.null(attr(.x, which = "warningMsg")))
      if (sum(purrr::map_lgl(ls.models, ~all(purrr::map_lgl(.x, is.null)))) == length(ls.models)) {
        stop("The model matrix is rank deficient. LS residuals cannot be calculated.", 
             "\nUse the 'include.ls = FALSE' parameter to get EB residuals only.")
        }
      if (sum(ls.warnings) > 0) {
        warning("The model matrix is likely rank deficient. Some LS residuals cannot be calculated.",
                "\nIt is recommended to use EB (.ranef) group level residuals for this model.")
      }
      
      # to fix the order 
      row.order <- purrr::map(ls.models, ~row.names(coef(.x)))
      row.order2 <- c()
      for (i in names(ls.models)) {
        row.order2 <- append(row.order2, 
                             stringr::str_c(row.order[[i]], i, sep = ":"))
      }
      
      # get the ranef
      ls.ranef <- purrr::map(ls.models, 
                        ~lme4::ranef(.x, standard = standardize)[ranef_names])
     
      # reconstruct the data frame with group id
      ls.resid <- tibble::tibble(
        group = row.order2,
        purrr::map_dfr(ls.ranef, ~dplyr::bind_cols(.x)))
      
      if(ncol(ls.resid) != length(ranef_names) + 1) {
        warning("The model contains a random effect term for variables with no fixed effect term.",
                "\nSome LS group level residuals cannot be calculated.")
      }

      return(ls.resid)
    }
  }
}

# Solving issue for R CMD check with tidyselect::where 
# This is currently the recommended practice
utils::globalVariables("where")

#' @export
#' @rdname LSresids.mer
#' @method LSresids lme
LSresids.lme <- function(object, level, standardize = FALSE, ...){
  #GLM OR NRESTED MODEL CHECK
  if(!level %in% c(1, names(object$groups))) {
    stop("level can only be 1 or a grouping factor from the fitted model.")
  }
  #NESTED MODEL CHECK
  if(!is.null(standardize) && !standardize %in% c(FALSE, TRUE, "semi")) {
    stop("standardize can only be specified to be logical or 'semi' .")
  }
  
  if(level == 1){
    y <- nlme::getResponse(object)
    y <- na.exclude(y)
    X <- model.matrix(object, data = object$data)
    g <- object$groups[length(object$groups)]
    names(g) <- "group"
    
    frame <- cbind(y = y, X, g)                  #bind y, Xs, group
    frame <- frame[,-2]                          #remove intercept column
    
    frame.split <- split(frame, frame[,"group"]) #split on groups
    frame.lgl <- 
      purrr::map(frame.split, function(split) {  #check if column elements are equal
        purrr::map_lgl(split, ~all(.x == .x[1]))
      })
    
    for (i in 1:length(frame.split)) {           #remove the all equal columns
      frame.lgl[[i]][1] <- FALSE                 #make sure we never remove response
      frame.split[[i]] <- frame.split[[i]][!frame.lgl[[i]]]
    }
    
    ls.models <-                                 #fit the lm
      purrr::map(frame.split, function(split){
        lm(y ~ ., data = split)
      })
    
    ls.rank <- purrr::map_dbl(ls.models, ~qr(.x)$rank) #check rank deficiency
    deficient.groups <- names(which(ls.rank != max(ls.rank)))
    if(length(deficient.groups) != 0) {
      warning("LS residuals might be inaccurate as one or more groups are rank deficient.", 
              "\nUse the 'include.ls = FALSE' parameter to get EB residuals only.")
    }    
    
    ls.residuals <- purrr::map(ls.models, resid) #calculate residuals within group
    ls.fitted <- purrr::map(ls.models, fitted)   #calculate fitted values
    
    row.order <- unlist(purrr::map(frame.split, function(x) row.names(x)))
    return.df <- data.frame(.ls.resid = unlist(ls.residuals),
                            .ls.fitted = unlist(ls.fitted))
    
    if(!is.null(standardize) && standardize == "semi"){
      ls.hat <- unlist(purrr::map(ls.models, ~lm.influence(.x)$hat))
      semi.std.resid  <- unlist(ls.residuals) / sqrt(1 - ls.hat)
      semi.std.resid[is.infinite(semi.std.resid)] <- NA
      
      return.df <- data.frame(.semi.ls.resid = semi.std.resid, 
                              .ls.fitted = unlist(ls.fitted))
    }
    
    if(!is.null(standardize) && standardize == TRUE){
      ls.rstandard <- unlist(purrr::map(ls.models, rstandard))
      ls.rstandard[is.infinite(ls.rstandard)] <- NA
      
      return.df <- data.frame(.std.ls.resid = ls.rstandard, 
                              .ls.fitted = unlist(ls.fitted))
    }
    
    rownames(return.df) <- row.order
    return(return.df)
  }
  
  if(level != 1){
    fixed <- as.character(formula(object))
    
    data <- object$data %>%
      dplyr::mutate(dplyr::across(where(is.character), ~ as.factor(.x))) %>%
      as.data.frame()
    
    n.ranefs <- length(names(object$groups))
    if (n.ranefs == 1){
      ranef_names <- names( nlme::ranef(object) )
    } else {
      ranef_names <- names( nlme::ranef(object)[[level]] )
    }
    
    if(level == names(object$groups)[1]) { #highest level
      form <- paste(fixed[2], fixed[1], fixed[3], "|", level)
      
      ls.models <- suppressWarnings(lme4::lmList(formula = formula(form), data = data))
      if (sum(purrr::map_lgl(ls.models, is.null)) == length(ls.models)) {
        stop("The model matrix is rank deficient. LS residuals cannot be calculated.",
             "\nUse the 'include.ls = FALSE' parameter to get EB residuals only.")
      }
      if (!is.null(attr(ls.models, which = "warningMsg"))) {
        warning("The model matrix is likely rank deficient. Some LS residuals cannot be calculated.",
                "\nIt is recommended to use EB (.ranef) group level residuals for this model.")
      }
      if(standardize == "semi") standardize <- FALSE
      ls.ranef <- lme4::ranef(ls.models, standard = standardize)[ranef_names]
      ls.resid <- purrr::map_dfc(ls.ranef, ~.x)
      ls.resid <- tibble::tibble(group = row.names(coef(ls.models)),
                                 ls.resid)
      
      if(ncol(ls.resid) != length(ranef_names) + 1) {
        warning("The model contains a random effect term for variables with no fixed effect term.",
                "\nSome LS group level residuals cannot be calculated.")
      }
      
      return(ls.resid)
      
    } else { #middle level
      # group by (specified level + 1)
      higher.level <- names(object$groups)[which(names(object$groups) == level) - 1]
      split_data <- split(data, data[,higher.level])
      # remove any empty groups
      split_data <- split_data[which(purrr::map_lgl(split_data, ~(nrow(.x)) != 0))]
      
      # for each group use ranef as above
      form <- paste(fixed[2], fixed[1], fixed[3], "|", level) 
      
      # recombine data frame
      ls.models <- suppressWarnings(purrr::map(split_data,
                                               ~lme4::lmList(formula = formula(form), data = .x)))
      ls.warnings <- purrr::map_lgl(ls.models, ~!is.null(attr(.x, which = "warningMsg")))
      if (sum(purrr::map_lgl(ls.models, ~all(purrr::map_lgl(.x, is.null)))) == length(ls.models)) {
        stop("The model matrix is rank deficient. LS residuals cannot be calculated.",
             "\nUse the 'include.ls = FALSE' parameter to get EB residuals only.")
      }
      if (sum(ls.warnings) > 0) {
        warning("The model matrix is likely rank deficient. Some LS residuals cannot be calculated.",
                "\nIt is recommended to use EB (.ranef) group level residuals for this model.")
      }
      
      # to fix the order 
      row.order <- purrr::map(ls.models, ~row.names(coef(.x)))
      row.order2 <- c()
      for (i in names(ls.models)) {
        row.order2 <- append(row.order2, 
                             stringr::str_c(i, row.order[[i]], sep = "/"))
      }
      
      # get the ranef
      ls.ranef <- purrr::map(ls.models, 
                             ~lme4::ranef(.x, standard = standardize)[ranef_names])
      
      # reconstruct the data frame with group id
      ls.resid <- tibble::tibble(
        group = row.order2,
        purrr::map_dfr(ls.ranef, ~dplyr::bind_cols(.x)))
      
      if(ncol(ls.resid) != length(ranef_names) + 1) {
        warning("The model contains a random effect term for variables with no fixed effect term.",
                "\nSome LS group level residuals cannot be calculated.")
      }
      
      return(ls.resid)
    }
  }
}
