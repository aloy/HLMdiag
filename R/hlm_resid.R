#' @export
hlm_resid <- function(object, ...){
  UseMethod("hlm_resid", object)
}

#' @export
#' @rdname hlm_resid.lmerMod
#' @method hlm_resid default
#' @S3method hlm_resid default
hlm_resid.default <- function(object, ...){
  stop(paste("there is no hlm_resid() method for objects of class",
             paste(class(object), collapse=", ")))
}

#' Calculating residuals from HLMs
#'
#' \code{hlm_resid} takes a hierarchical linear model fit as a
#' \code{lmerMod} object and adds information about each observation's 
#' residuals and predicted values.
#' 
#' This function extract residuals and predicted values from the model, using
#' least squares (LS) and Empirical Bayes (EB) methods, and appends them to the
#' model data. This unified framework enables the analyst to more easily conduct
#' an upward residual analysis during model exploration/checking.
#'
#' @export
#' @method hlm_resid lmerMod
#' @S3method hlm_resid lmerMod
#' @aliases hlm_resid
#' @param object an object of class \code{lmerMod}.
#' @param level which residuals should be extracted: 1 for within-group
#'   (case-level) residuals, the name of a grouping factor (as defined in
#'   \code{flist} of the \code{lmerMod} object) for between-group residuals
#' @param standardize if \code{standardize = TRUE} the standardized residuals
#'   will be returned; if \code{standardize = "semi"} then the semi-standardized
#'   level-1 residuals will be returned
#' @param sim optional argument giving the data frame used for LS residuals.
#'   This is used mainly for dealing with simulations.
#' @param ... do not use
#' @details The \code{hlm_resid} function provides a wrapper that will extract
#' residuals and predicted values from a fitted \code{lmerMod} object. 
#' The function provides access to 
#' residual quantities already made available by the functions \code{resid},
#' \code{predict}, and \code{ranef}, but adds additional functionality. Below is
#' a list of types of residuals and predicted values that are extracted and
#' appended to the model data.
#' \describe{
#' \item{raw level-1 LS residuals}{These are equivalent to the residuals extracted
#' by \code{resid} if \code{level = 1}, \code{type = "EB"}, and 
#' \code{standardize = FALSE} is specified. }
#' \item{level-1 LS fitted values}{The predicted values }
#' }
#' Note that \code{standardize = "semi"} is only implemented for level-1 LS residuals.
hlm_resid.lmerMod <- function(object, level = 1, standardize = FALSE, sim = NULL, ...) {
  
  if(!level %in% c(1, names(object@flist))) {
    stop("level can only be 1 or a grouping factor from the fitted model.")
  }
    if(!is.null(standardize) && !standardize %in% c(FALSE, TRUE, "semi")) {
    stop("standardize can only be specified to be logical or 'semi'.")
  }
  
  if(level == 1) { 
    # LS Residuals and Fitted
    ls.resid <- LSresids(object, level = 1, stand = standardize, sim = sim)
    ls.resid <- ls.resid[order(as.numeric(rownames(ls.resid))),]
    
    if (standardize == FALSE) {
      ls.resid <- data.frame(ls.resid["LS.resid"], ls.resid["fitted"])
      names(ls.resid) <- c(".ls.resid", ".ls.fitted")
      
    } else if (standardize == TRUE) {
      ls.resid <- data.frame(ls.resid["std.resid"], ls.resid["fitted"])
      names(ls.resid) <- c(".std.ls.resid", ".ls.fitted")
      
    } else {
      ls.resid <- data.frame(ls.resid["semi.std.resid"], ls.resid["fitted"])
      names(ls.resid) <- c(".semi.ls.resid", ".ls.fitted")
    }
    
    # EB Residuals
    if (standardize == TRUE) {
      mats <- .lmerMod_matrices(object)
      p_diag <- diag(mats$P)
      eb.resid <- data.frame(.std.resid = 
                               resid(object) / ( lme4::getME(object, "sigma") * sqrt(p_diag) ))
      
    } else {
      eb.resid <- data.frame(.resid = resid(object))
    }
    # EB Fitted
    eb.fitted <- data.frame(.fitted = lme4::getME(object, "mu"))
    
    # Marginal Residuals
    mr <- object@resp$y - lme4::getME(object, "X") %*% lme4::fixef(object)
    if (standardize == TRUE) {
      sig0 <- lme4::getME(object, "sigma")
      ZDZt <- sig0^2 * crossprod( lme4::getME(object, "A") )
      n    <- nrow(ZDZt)
      
      R      <- Diagonal( n = n, x = sig0^2 )
      V      <- R + ZDZt
      V.chol <- chol( V )
      
      Lt <- solve(t(V.chol))
      mar.resid <- data.frame(.std.mar.resid = (Lt %*% mr)[,1])
      
    } else {
      mar.resid <- data.frame(.mar.resid = mr[,1])
    }
    # Marginal Fitted
    mar.fitted  <- data.frame(.mar.fitted = predict(object, re.form = ~0))
    
    # Assemble Tibble
    return.tbl <- tibble::tibble(object@frame,
                                 eb.resid,
                                 eb.fitted,
                                 ls.resid,
                                 mar.resid,
                                 mar.fitted)
    return(return.tbl)
  }
  
  if (level %in% names(object@flist)) {
    # EB Residuals
    eb.resid <- lme4::ranef(object)[[level]]
    eb.resid <- janitor::clean_names(eb.resid)
    groups <- rownames(eb.resid)
    if (standardize == TRUE) {
      se.re <- se.ranef(object)[[level]]
      eb.resid <- eb.resid/se.re
      eb.names <- paste0(".std.ranef.", names(eb.resid))
    } else {
      eb.names <- paste0(".ranef.", names(eb.resid))
    }
    
    # LS Residuals
    ls.resid <- LSresids(object, level = level, stand = standardize, sim = sim)
    ls.resid <- ls.resid[match(groups, ls.resid$group),] # fix order
    ls.resid <- janitor::clean_names(ls.resid)
    if (standardize == TRUE) {
      ls.names <- paste0(".std.ls.", names(ls.resid))
    } else {
      ls.names <- paste0(".ls.", names(ls.resid))
    }
    
    # Grab level 2 variables
    # adjust_lmList method
    #g <- level
    #if(stringr::str_detect(level, ":")) {
    #  vars <- stringr::str_split(level, ":")[[1]]
    #  g <- vars[which(!vars %in% names(object@flist))]
    #}
    #
    #lvl1_vars <- NULL
    #fixed <- as.character(lme4::nobars( formula(object)))
    #form <- paste(fixed[2], fixed[1], fixed[3], "|", level)
    #try(lvl1_vars <- adjust_formula_lmList(formula(form), object@frame),
    #    silent = TRUE)
    #
    #if(is.null(lvl1_vars)){
    #  # model is too simple, adjust_formula fails
    #  suppressMessages(return.tbl <- tibble::tibble(
    #    groups, eb.resid, ls.resid[,-1], .name_repair = "universal"))
    #  names(return.tbl) <- c(level, eb.names, ls.names)
    #  
    #  return(return.tbl)
    #  
    #} else {
    #  lvl1_vars <- unique(unlist(purrr::map(lvl1_vars, all.names)))
    #  index <- which(!names(object@frame) %in% lvl1_vars)
    #  #use select in dplyr
    #  group_vars <- object@frame %>%
    #    dplyr::select(index)
    #  if(!is.character(group_vars[,g])) {
    #    group_vars[,g] <- as.character(group_vars[,g])
    #  }
    #  suppressMessages(return.tbl <- tibble::tibble(
    #    groups, eb.resid, ls.resid[,-1], .name_repair = "universal"))
    #  names(return.tbl) <- c(level, eb.names, ls.names)
    #  return.tbl <- tibble::tibble(
    #    unique(dplyr::left_join(group_vars, return.tbl)))
    #  
    #  return(return.tbl)
    #}

      # Grab level specific variables
      # lmList method
      fixed <- as.character(lme4::nobars( formula(object)))
      n.ranefs <- length(names(object@flist))
      ranef_names <- names( lme4::ranef(object)[[level]] )
      
      
      if(level == names(object@flist)[n.ranefs]){ # highest level
        form <- paste(fixed[2], fixed[1], fixed[3], "|", level)
        
        # Use lmList
        g.list <- lme4::lmList(formula(form), data = object@frame)
       
        # Checking if all of the values for a coef are NAs
        g.index <- which(purrr::map_lgl(coef(g.list), ~all(is.na(.x))))
        g.names <- names(g.index)
        
        # Match that index back to object@frame
        g.exp <- stringr::str_c(g.names, collapse = "|")
        g.index.frame <- which( 
          stringr::str_detect(g.exp, names(object@frame)))
        g.vars <- object@frame %>%
          dplyr::select(level, g.index.frame)
        g.vars <- unique(g.vars)
        
        # Assemble data frame
        return.tbl <- tibble::tibble(
          groups, eb.resid, ls.resid[,-1], .name_repair = "universal")
        names(return.tbl) <- c(level, eb.names, ls.names[-1])
        return.tbl <- tibble::tibble(
          unique(dplyr::left_join(g.vars, return.tbl)))
        
        return(return.tbl)
        
      } else { # inner level
        # Extract correct grouping variable
        level.var <- stringr::str_split(level, ":")[[1]]
        level.var <- level.var[which(!level.var %in% names(object@flist))]
        form <- paste(fixed[2], fixed[1], fixed[3], "|", level.var)
        
        # Use lmList
        g.list <- lme4::lmList(formula(form), data = object@frame)
        
        # Checking if all of the values for a coef are NAs
        g.index <- which(purrr::map_lgl(coef(g.list), ~all(is.na(.x))))
        g.names <- names(g.index)
        
        # Match that index back to object@frame
        higher.level <- names(object@flist[which(names(object@flist) == level) + 1])
        g.exp <- stringr::str_c(g.names, collapse = "|")    
        g.index.frame <- which( 
          stringr::str_detect(g.exp, names(object@frame)))
        g.vars <- object@frame %>%
          dplyr::select(level.var, higher.level, g.index.frame)
        g.vars <- unique(g.vars)
        
        # HERE I want to add a column to the above data frame that concatenates
        # the string stored in g.var[level] with the string stored in
        # g.var[higher.level] seperated by a colon.
        
      }
    
  }
}
