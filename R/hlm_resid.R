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
#' \code{hlm_resid} takes a hierarchical linear model fit as a \code{lmerMod} or
#' \code{lme} object and adds residuals and fitted values for individual
#' observations or higher level clusters.
#' 
#' This function extracts residuals and predicted values from the model, using
#' Least Squares (LS) and Empirical Bayes (EB) methods, and appends them to the
#' model data. This unified framework enables the analyst to more easily conduct
#' an upward residual analysis during model exploration/checking.
#'
#' @export
#' @method hlm_resid lmerMod
#' @S3method hlm_resid lmerMod
#' @aliases hlm_resid
#' @param object an object of class \code{lmerMod} or \code{lme}.
#' @param level which residuals should be extracted: 1 for within-group
#'   (case-level) residuals, the name of a grouping factor for between-group
#'   residuals (as defined in \code{flist} in \code{lmerMod} objects or in
#'   \code{groups}} in \code{lme} objects)
#' @param standardize for any level, if \code{standardize = TRUE} the
#'   standardized residuals will be returned for any group; for level-1 only, if
#'   \code{standardize = "semi"} then the semi-standardized level-1 residuals
#'   will be returned
#' @param ls.include a logical indicating if LS residuals be included in the
#'   return tibble. \code{include.ls = FALSE} decreases runtime substantially.
#' @param sim optional argument giving the data frame used for LS residuals.
#'   This is used mainly for dealing with simulations.
#' @param ... do not use
#' @details The \code{hlm_resid} function provides a wrapper that will extract
#' residuals and predicted values from a fitted \code{lmerMod} or \code{lme}
#' object.
#' The function provides access to residual quantities already made available by
#' the functions \code{resid}, \code{predict}, and \code{ranef}, but adds
#' additional functionality. Below is a list of types of residuals and predicted
#' values that are extracted and appended to the model data.
#' \describe{
#' \item{\strong{level-1 residuals}}{}
#' \item{\code{.resid} and \code{.fitted}}{Residuals calculated using
#' the EB method (using maximum likelihood). Level-1 EB residuals are interrelated
#' with higher level residuals. Equivalent to the residuals extracted by
#' \code{resid(object)} and \code{lme4::getME(object, "mu")} respectively. When
#' \code{standardize = TRUE}, residuals are standardized by sigma components of
#' the model object.}
#' \item{\code{.ls.resid} and \code{.ls.fitted}}{Residuals calculated calculated
#' by fitting separate LS regression models for each group. Level-1 LS residuals
#' are unconfounded by higher level residuals, but unreliable for small
#' within-group sample sizes. Uses \code{lme4::lmList} and \code{adjust_lmList} in
#' fitting individual LS models.}
#' \item{\code{.mar.resid} and \code{.mar.fitted}}{Marginal residuals only
#' consider the fixed effect portion of the estimates. They are calculated by 
#' \code{object@@resp$y - lme4::getME(object, "X") \%*\% lme4::fixef(object)}.
#' When \code{standardize = TRUE}, cholskey residuals are returned.}
#' \item{\strong{higher-level residuals} (random effects)}{}
#' \item{\code{.ranef.*}}{The group level random effects using the EB method of
#' estimating parameters. Equivalent to \code{lme4::ranef} on the specified
#' level. EB residuals are prefered at higher levels LS residuals are dependent
#' on a large sample size.}
#' \item{\code{.ls.*}}{The group level random effects using the LS method of
#' estimating parameters. Calculated using \code{ranef} on a \code{lmList4}
#' object to compare the random effects of individual models to the global
#' model.}
#' }
#' Note that \code{standardize = "semi"} is only implemented for level-1 LS residuals.
hlm_resid.lmerMod <- function(object, level = 1, standardize = FALSE, include.ls = TRUE, sim = NULL, ...) {
  
  if(!level %in% c(1, names(object@flist))) {
    stop("level can only be 1 or the following grouping factors from the fitted model: \n", 
         stringr::str_c(names(object@flist), collapse = ", "))
  }
  if(!is.null(standardize) && !standardize %in% c(FALSE, TRUE, "semi")) {
    stop("standardize can only be specified to be logical or 'semi'.")
  }
  
  if(level == 1) { 
    # LS Residuals and Fitted
    if(include.ls == TRUE) {
      ls.resid <- LSresids(object, level = 1, standardize = standardize, sim = sim)
      ls.resid <- ls.resid[order(as.numeric(rownames(ls.resid))),]
      # 
      # if (standardize == FALSE) {
      #   ls.resid <- ls.resid %>% 
      #     select(.ls.resid = LS.resid, .ls.fitted = fitted)
      #   
      # } else if (standardize == TRUE) {
      #   ls.resid <- ls.resid %>% 
      #     select(.std.ls.resid = std.resid, .ls.fitted = fitted)
      #   
      # } else {
      #   ls.resid <- ls.resid %>% 
      #     select(.semi.ls.resid = semi.std.resid, .ls.fitted = fitted)
      # }
    }
    
    # EB Residuals
    if (standardize == TRUE) {
      eb.resid <- data.frame(.std.resid = resid(object, scale = TRUE))
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
      mar.resid <- data.frame(.chol.mar.resid = (Lt %*% mr)[,1])
      
    } else {
      mar.resid <- data.frame(.mar.resid = mr[,1])
    }
    # Marginal Fitted
    mar.fitted  <- data.frame(.mar.fitted = predict(object, re.form = ~0))
    
    # Assemble Tibble
    if (include.ls == TRUE) {
      return.tbl <- tibble::tibble(object@frame,
                                   eb.resid,
                                   eb.fitted,
                                   ls.resid,
                                   mar.resid,
                                   mar.fitted)
    } else { 
      return.tbl <- tibble::tibble(object@frame,
                                   eb.resid,
                                   eb.fitted,
                                   mar.resid,
                                   mar.fitted)
    }
    
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
    if (include.ls == TRUE) {
      ls.resid <- LSresids(object, level = level, stand = standardize, sim = sim)
      ls.resid <- ls.resid[match(groups, ls.resid$group),] # fix order
      ls.resid <- janitor::clean_names(ls.resid)
      if (standardize == TRUE) {
        ls.names <- paste0(".std.ls.", names(ls.resid))
      } else {
        ls.names <- paste0(".ls.", names(ls.resid))
      }
    }
    
    # Grab level specific variables
    fixed <- as.character(lme4::nobars( formula(object)))
    n.ranefs <- length(names(object@flist))
    ranef_names <- names( lme4::ranef(object)[[level]] )
    
    
    if(level == names(object@flist)[n.ranefs]){ # highest level
      form <- paste(fixed[2], fixed[1], fixed[3], "|", level)
      
      # Use lmList
      g.list <- suppressWarnings(lme4::lmList(formula(form), data = object@frame))
     
      # Checking if all of the values for a coef are NAs
      g.index <- which(purrr::map_lgl(coef(g.list), ~all(is.na(.x))))
      g.names <- names(g.index)
      # Get rid of interaction terms
      interaction.index <- stringr::str_detect(g.names, ":")
      g.names <- g.names[!interaction.index]
      
      # Match that index back to object@frame
      g.exp <- stringr::str_c(g.names, collapse = "|")
      g.index.frame <- which( 
        stringr::str_detect(g.exp, names(object@frame)))
      g.vars <- object@frame %>%
        dplyr::select(all_of(level), all_of(g.index.frame))
      g.vars <- unique(g.vars)
      
      # Assemble data frame
      if (include.ls == TRUE) {
        return.tbl <- suppressMessages(tibble::tibble(
          groups, eb.resid, ls.resid[,-1], .name_repair = "universal"))
        names(return.tbl) <- c(level, eb.names, ls.names[-1])
        if(class(g.vars[level][[1]]) != class(return.tbl[level][[1]])){
          g.vars[level][[1]] <- as.character(g.vars[level][[1]])
          return.tbl[level][[1]] <- as.character(return.tbl[level][[1]])
        }
        return.tbl <- tibble::tibble(
          unique(suppressMessages(dplyr::left_join(g.vars, return.tbl))))
      
      } else {
        return.tbl <- tibble::tibble(groups, eb.resid)
        names(return.tbl) <- c(level, eb.names)
        if(class(g.vars[level][[1]]) != class(return.tbl[level][[1]])){
          g.vars[level][[1]] <- as.character(g.vars[level][[1]])
          return.tbl[level][[1]] <- as.character(return.tbl[level][[1]])
        }
        return.tbl <- tibble::tibble(
          unique(suppressMessages(dplyr::left_join(g.vars, return.tbl))))
      }
      
      
      return(return.tbl)
      
    } else { # inner level
      # Extract correct grouping variable
      level.var <- stringr::str_split(level, ":")[[1]]
      level.var <- level.var[which(!level.var %in% names(object@flist))]
      form <- paste(fixed[2], fixed[1], fixed[3], "|", level.var)
      
      # Use lmList
      g.list <- suppressWarnings(lme4::lmList(formula(form), data = object@frame))
      
      # Checking if all of the values for a coef are NAs
      g.index <- which(purrr::map_lgl(coef(g.list), ~all(is.na(.x))))
      g.names <- names(g.index)
      # Remove interaction terms
      interaction.index <- stringr::str_detect(g.names, ":")
      g.names <- g.names[!interaction.index]
      
      # Match that index back to object@frame
      higher.level <- names(object@flist[which(names(object@flist) == level) + 1])
      g.exp <- stringr::str_c(g.names, collapse = "|")    
      g.index.frame <- which( 
        stringr::str_detect(g.exp, names(object@frame)))
      g.vars <- object@frame %>%
        dplyr::select(all_of(level.var), all_of(higher.level), all_of(g.index.frame))
      g.vars <- unique(g.vars)
      
      # Add the group variable
      g.vars$group <- rep(NA, nrow(g.vars))
      for (i in 1:nrow(g.vars)){
        g.vars$group[i] <- stringr::str_c(g.vars[level.var][i,], 
                                      g.vars[higher.level][i,], sep = ":")
      }
      g.vars <- g.vars %>%
        select(ncol(g.vars), 1:(ncol(g.vars)-1))
      
      # Assemble data frame
      if (include.ls == TRUE) {
        return.tbl <- suppressMessages(tibble::tibble(
          groups, eb.resid, ls.resid[,-1], .name_repair = "universal"))
        names(return.tbl) <- c("group", eb.names, ls.names[-1])
        return.tbl <- tibble::tibble(
          unique(suppressMessages(dplyr::left_join(g.vars, return.tbl))))
        
      } else {
        return.tbl <- tibble::tibble(groups, eb.resid)
        names(return.tbl) <- c("group", eb.names)
        return.tbl <- tibble::tibble(
          unique(suppressMessages(dplyr::left_join(g.vars, return.tbl))))
      }
        
      return(return.tbl)
      
    }
  }
}

#' @export
#' @rdname hlm_resid.lmerMod
#' @method hlm_resid lme
#' @S3method hlm_resid lme
hlm_resid.lme <- function(object, level = 1, standardize = FALSE, include.ls = TRUE, sim = NULL, ...) {
  if(!level %in% c(1, names(object$groups))) {
    stop("level can only be 1 or the following grouping factors from the fitted model: \n", 
         stringr::str_c(names(object$groups), collapse = ", "))
  }
  if(!is.null(standardize) && !standardize %in% c(FALSE, TRUE, "semi")) {
    stop("standardize can only be specified to be logical or 'semi'.")
  }
  
  if(level == 1) { 
    # LS Residuals and Fitted
    if(include.ls == TRUE) {
      ls.resid <- LSresids(object, level = 1, stand = standardize, sim = sim)
      ls.resid <- ls.resid[order(as.numeric(rownames(ls.resid))),]
      
      if (standardize == FALSE) {
        ls.resid <- ls.resid %>% 
          select(.ls.resid = LS.resid, .ls.fitted = fitted)
        
      } else if (standardize == TRUE) {
        ls.resid <- ls.resid %>% 
          select(.std.ls.resid = std.resid, .ls.fitted = fitted)
        
      } else {
        ls.resid <- ls.resid %>% 
          select(.semi.ls.resid = semi.std.resid, .ls.fitted = fitted)
      }
    }
    
    # EB Residuals
    if(standardize == TRUE) {
      eb.resid <- data.frame(.std.resid = resid(object, type = "normalized"))
    } else { 
      eb.resid <- data.frame(.resid = resid(object, type = "response"))
    }

    # EB Fitted
    eb.fitted <- data.frame(.fitted = fitted(object))
    # EB Fitted
    
    # Marginal Residuals
    mr <- resid(object, type="response", level=0)
    if (standardize == TRUE) {
      V      <- extract_design(object)$V
      V.chol <- chol( V )
      
      Lt <- solve(t(V.chol))
      mar.resid <- data.frame(.chol.mar.resid = (Lt %*% mr)[,1])
    } else {
      mar.resid <- data.frame(.mar.resid = mr)
    }
    # Marginal Fitted
    mar.fitted  <- data.frame(.mar.fitted = fitted(object, level=0))
    
    # Assemble Tibble
    fixed <- as.character(formula(object))
    dataform <- paste(fixed[2], fixed[1], fixed[3], " + ", 
                      paste(names(object$groups), collapse = " + "))
    data <- object$data %>%
      dplyr::mutate(across(where(is.character), ~ as.factor(.x))) %>%
      as.data.frame()
    if (include.ls == TRUE) {
      return.tbl <- tibble::tibble(model.frame(formula(dataform), data),
                                   eb.resid,
                                   eb.fitted,
                                   ls.resid,
                                   mar.resid,
                                   mar.fitted)
    } else { 
      return.tbl <- tibble::tibble(model.frame(formula(dataform), data),
                                   eb.resid,
                                   eb.fitted,
                                   mar.resid,
                                   mar.fitted)
    }
    
    return(return.tbl)
  }
  
  if (level %in% names(object$groups)) {
    # EB Residuals
    if(standardize == "semi") standardize <- FALSE
    if (length(object$groups) != 1) {
      eb.resid <- ranef(object, standard = standardize)[[level]]
    } else { 
      eb.resid <- ranef(object, standard = standardize)
    }

    eb.resid <- janitor::clean_names(eb.resid)
    groups <- rownames(eb.resid)
    if (standardize == TRUE) {
      eb.names <- paste0(".std.ranef.", names(eb.resid))
    } else {
      eb.names <- paste0(".ranef.", names(eb.resid))
    }
    
    # LS Residuals
    if (include.ls == TRUE) {
      ls.resid <- LSresids(object, level = level, stand = standardize, sim = sim)
      ls.resid <- ls.resid[match(groups, ls.resid$group),] # fix order
      ls.resid <- janitor::clean_names(ls.resid)
      if (standardize == TRUE) {
        ls.names <- paste0(".std.ls.", names(ls.resid))
      } else {
        ls.names <- paste0(".ls.", names(ls.resid))
      }
    }
    
    # Grab level specific variables
    fixed <- as.character(formula(object))    
    n.ranefs <- length(names(object$groups))
    if (n.ranefs == 1){
      ranef_names <- names( nlme::ranef(object) )
    } else {
      ranef_names <- names( nlme::ranef(object)[[level]] )
    }
    data <- object$data %>%
      dplyr::mutate(across(where(is.character), ~ as.factor(.x))) %>%
      as.data.frame()
    
    form <- paste(fixed[2], fixed[1], fixed[3], "|", level)
    
    # Use lmList
    g.list <- suppressWarnings(lme4::lmList(formula(form), data = data))
    
    # Checking if all of the values for a coef are NAs
    g.index <- which(purrr::map_lgl(coef(g.list), ~all(is.na(.x))))
    g.names <- names(g.index)
    # remove interaction terms
    interaction.index <- stringr::str_detect(g.names, ":")
    g.names <- g.names[!interaction.index]
    
    # Match that index back to data
    g.exp <- stringr::str_c(g.names, collapse = "|")    
    g.index.frame <- which( 
      stringr::str_detect(g.exp, names(data)))
    
    if(level == names(object$groups)[1]){ # highest level
      g.vars <- data %>%
        dplyr::select(all_of(level), all_of(g.index.frame))
      g.vars <- unique(g.vars)
  
      # Assemble data frame
      if (include.ls == TRUE) {
        return.tbl <- suppressMessages(tibble::tibble(
          groups, eb.resid, ls.resid[,-1], .name_repair = "universal"))
        names(return.tbl) <- c(level, eb.names, ls.names[-1])
        if(class(g.vars[level][[1]]) != class(return.tbl[level][[1]])){
          g.vars[level][[1]] <- as.character(g.vars[level][[1]])
          return.tbl[level][[1]] <- as.character(return.tbl[level][[1]])
        }
        return.tbl <- tibble::tibble(
          unique(suppressMessages(dplyr::left_join(g.vars, return.tbl))))
        
      } else {
        return.tbl <- tibble::tibble(groups, eb.resid)
        names(return.tbl) <- c(level, eb.names)
        if(class(g.vars[level][[1]]) != class(return.tbl[level][[1]])){
          g.vars[level][[1]] <- as.character(g.vars[level][[1]])
          return.tbl[level][[1]] <- as.character(return.tbl[level][[1]])
        }
        return.tbl <- tibble::tibble(
          unique(suppressMessages(dplyr::left_join(g.vars, return.tbl))))
      }
      
      return(return.tbl)
      
    } else { # inner level
      # Match that index back to data
      higher.level <- names(object$groups[which(names(object$groups) == level) -1])
      g.vars <- data %>%
        dplyr::select(all_of(higher.level), all_of(level), all_of(g.index.frame))
      g.vars <- unique(g.vars)
      
      # Add group variable
      g.vars$group <- rep(NA, nrow(g.vars))
      for (i in 1:nrow(g.vars)){
        g.vars$group[i] <- stringr::str_c(g.vars[higher.level][i,], 
                                          g.vars[level][i,], sep = "/")
      }
      g.vars <- g.vars %>%
        select(ncol(g.vars), 1:(ncol(g.vars)-1))
      
      # Assemble data frame
      if (include.ls == TRUE) {
        return.tbl <- suppressMessages(tibble::tibble(
          groups, eb.resid, ls.resid[,-1], .name_repair = "universal"))
        names(return.tbl) <- c("group", eb.names, ls.names[-1])
        return.tbl <- tibble::tibble(
          unique(suppressMessages(dplyr::left_join(g.vars, return.tbl))))
        
      } else {
        return.tbl <- tibble::tibble(groups, eb.resid)
        names(return.tbl) <- c("group", eb.names)
        return.tbl <- tibble::tibble(
          unique(suppressMessages(dplyr::left_join(g.vars, return.tbl))))
      }
      
      return(return.tbl)
    }
  }
}