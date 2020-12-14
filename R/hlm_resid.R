#' @export
hlm_resid <- function(object, ...){
  UseMethod("hlm_resid", object)
}

#' @export
#' @rdname hlm_resid.lmerMod
#' @method hlm_resid default
hlm_resid.default <- function(object, ...){
  stop(paste("there is no hlm_resid() method for objects of class",
             paste(class(object), collapse=", ")))
}

#' Calculating residuals from HLMs
#'
#' \code{hlm_resid} takes a hierarchical linear model fit as a \code{lmerMod} or
#' \code{lme} object and extracts residuals and predicted values from the model,
#' using Least Squares (LS) and Empirical Bayes (EB) methods. It then appends
#' them to the model data frame in the form of a tibble inspired by the \code{augment}
#' function in \code{broom}. This unified framework enables the analyst to more
#' easily conduct an upward residual analysis during model exploration/checking.
#'
#' @export
#' @method hlm_resid lmerMod
#' @aliases hlm_resid
#' @param object an object of class \code{lmerMod} or \code{lme}.
#' @param level which residuals should be extracted: 1 for within-group
#'   (case-level) residuals, the name of a grouping factor for between-group
#'   residuals (as defined in \code{flist} in \code{lmerMod} objects or in
#'   \code{groups} in \code{lme} objects)
#' @param standardize for any level, if \code{standardize = TRUE} the
#'   standardized residuals will be returned for any group; for level-1 only, if
#'   \code{standardize = "semi"} then the semi-standardized level-1 residuals
#'   will be returned
#' @param include.ls a logical indicating if LS residuals be included in the
#'   return tibble. \code{include.ls = FALSE} decreases runtime substantially.
#' @param data if \code{na.action = na.exclude}, the user must provide the data
#'   set used to fit the model, otherwise \code{NULL}.
#' @param ... do not use
#' @details The \code{hlm_resid} function provides a wrapper that will extract
#' residuals and predicted values from a fitted \code{lmerMod} or \code{lme}
#' object. The function provides access to residual quantities already made available by
#' the functions \code{resid}, \code{predict}, and \code{ranef}, but adds
#' additional functionality. Below is a list of types of residuals and predicted
#' values that are extracted and appended to the model data.
#' \describe{
#' \item{\strong{level-1 residuals}}{}
#' \item{\code{.resid} and \code{.fitted}}{Residuals calculated using
#' the EB method (using maximum likelihood). Level-1 EB residuals are interrelated
#' with higher level residuals. Equivalent to the residuals extracted by
#' \code{resid(object)} and \code{predict(object)} respectively. When
#' \code{standardize = TRUE}, residuals are standardized by sigma components of
#' the model object.}
#' \item{\code{.ls.resid} and \code{.ls.fitted}}{Residuals calculated calculated
#' by fitting separate LS regression models for each group. Level-1 LS residuals
#' are unconfounded by higher level residuals, but unreliable for small
#' within-group sample sizes.}
#' \item{\code{.mar.resid} and \code{.mar.fitted}}{Marginal residuals only
#' consider the fixed effect portion of the estimates. Equivalent to
#' \code{resid(object, level=0)} in \code{nlme}, not currently implemented
#' within the \code{lme4::resid} function. When \code{standardize = TRUE},
#' Cholesky marginal residuals are returned.}
#' \item{\strong{higher-level residuals} (random effects)}{}
#' \item{\code{.ranef.var_name}}{The group level random effects using the EB method of
#' estimating parameters. Equivalent to \code{ranef(object)} on the specified
#' level. EB residuals are preferred at higher levels as LS residuals are dependent
#' on a large sample size.}
#' \item{\code{.ls.var_name}}{The group level random effects using the LS method of
#' estimating parameters. Calculated using \code{ranef} on a \code{lmList}
#' object to compare the random effects of individual models to the global
#' model.}
#' }
#' Note that \code{standardize = "semi"} is only implemented for level-1 LS residuals.
#' @author  Adam Loy \email{loyad01@@gmail.com}, Jack Moran, Jaylin Lowe
#' @keywords models regression
#' @seealso \code{\link{hlm_augment}}, \code{\link{resid}}, \code{\link{ranef}}
#' @references 
#' Hilden-Minton, J. (1995) Multilevel diagnostics for mixed and hierarchical 
#' linear models. University of California Los Angeles.
#' 
#' Houseman, E. A., Ryan, L. M., & Coull, B. A. (2004) 
#' Cholesky Residuals for Assessing Normal Errors in a Linear 
#' Model With Correlated Outcomes. 
#' \emph{Journal of the American Statistical Association}, 99(466), 383--394.
#' 
#' David Robinson and Alex Hayes (2020). broom: Convert Statistical Analysis
#' Objects into Tidy Tibbles. R package version 0.5.6.
#' \url{https://CRAN.R-project.org/package=broom}
#' @examples
#' data(sleepstudy, package = "lme4")
#' fm.lmer <- lme4::lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
#' fm.lme <- nlme::lme(Reaction ~ Days, random = ~Days|Subject, sleepstudy)
#' 
#' # level-1 and marginal residuals
#' fm.lmer.res1 <- hlm_resid(fm.lmer) ## raw level-1 and mar resids
#' fm.lmer.res1
#' fm.lme.std1 <- hlm_resid(fm.lme, standardize = TRUE) ## std level-1 and mar resids
#' fm.lme.std1
#' 
#' # level-2 residuals
#' fm.lmer.res2 <- hlm_resid(fm.lmer, level = "Subject") ## level-2 ranefs
#' fm.lmer.res2
#' fm.lme.res2 <- hlm_resid(fm.lme, level = "Subject", include.ls = FALSE) ##level-2 ranef, no LS
#' fm.lme.res2
hlm_resid.lmerMod <- function(object, level = 1, standardize = FALSE, include.ls = TRUE, data = NULL, ...) {
  
  if(!isNestedModel(object)){
    stop("hlm_resid is not currently implemented for non-nested models.")
  }
  if(!level %in% c(1, names(object@flist))) {
    stop("level can only be 1 or the following grouping factors from the fitted model: \n", 
         stringr::str_c(names(object@flist), collapse = ", "))
  }
  if(!is.null(standardize) && !standardize %in% c(FALSE, TRUE, "semi")) {
    stop("standardize can only be specified to be logical or 'semi'.")
  }
  if(class(attr(object@frame, "na.action")) == "exclude" && is.null(data) && level == 1){
    stop("Please provide the data frame used to fit the model. This is necessary when the na.action is set to na.exclude")
  }
  if(!is.null(getCall(object)$correlation)){
    warning("LS residuals for non-random correlation are not yet implemented")
    include.ls <- FALSE
  }

  # NA action
  if(class(attr(object@frame, "na.action")) == "exclude"){         #if na.exclude
    na.index <- which(!rownames(data) %in% rownames(object@frame))
    col.index <- which(colnames(data) %in% colnames(object@frame))
    data <- data[col.index]
  } else {                                                         #if na.omit
    data <- object@frame
  }
  
  if(level == 1) { 
    # LS Residuals and Fitted
    if(include.ls == TRUE) {
      ls.resid <- LSresids(object, level = 1, standardize = standardize)
      ls.resid <- ls.resid[order(as.numeric(rownames(ls.resid))),]
    }
    
    # EB Residuals
    if (standardize == TRUE) {
      eb.resid <- data.frame(.std.resid = resid(object, scale = TRUE))
    } else {
      eb.resid <- data.frame(.resid = resid(object))
    }
    # EB Fitted
    # eb.fitted <- data.frame(.fitted = lme4::getME(object, "mu"))
    eb.fitted <- data.frame(.fitted = predict(object))
    
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
    
    # NA Action
    if(class(attr(object@frame, "na.action")) == "exclude"){
      if(include.ls == TRUE){
        problem_dfs <- cbind(ls.resid, mar.resid)
        na.fix <- data.frame(LSR = rep(NA, length(na.index)), 
                             LSF = rep(NA, length(na.index)),
                             MR = rep(NA, length(na.index)))
        rownames(na.fix) <- na.index
        names(na.fix) <- c(names(ls.resid), names(mar.resid))
        problem_dfs <- rbind(problem_dfs, na.fix)
        problem_dfs <- problem_dfs[order(as.numeric(rownames(problem_dfs))),]
      } else {
        na.fix <- data.frame(MR = rep(NA, length(na.index)))
        rownames(na.fix) <- na.index
        names(na.fix) <- names(mar.resid)
        problem_dfs <- rbind(mar.resid, na.fix)
        problem_dfs <- data.frame(MR = problem_dfs[order(as.numeric(rownames(problem_dfs))),])
        names(problem_dfs) <- names(mar.resid)
      }
    } else {
      if(include.ls == TRUE){
        problem_dfs <- cbind(ls.resid, mar.resid)
      } else {
        problem_dfs <- mar.resid
      }
    }
    
    # Assemble Tibble
    if (include.ls == TRUE) {
      return.tbl <- tibble::tibble("id" = as.numeric(rownames(data)),
                                   data,
                                   eb.resid,
                                   eb.fitted,
                                   problem_dfs,
                                   mar.fitted)
    } else { 
      return.tbl <- tibble::tibble("id" = as.numeric(rownames(data)),
                                   data,
                                   eb.resid,
                                   eb.fitted,
                                   problem_dfs,
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
      ls.resid <- LSresids(object, level = level, stand = standardize)
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
    data <- object@frame
    names(data)[which(names(data) == fixed[2])] <- "y"
    fixed[2] <- "y"
    
    n.ranefs <- length(names(object@flist))
    ranef_names <- names( lme4::ranef(object)[[level]] )
    
    
    if(level == names(object@flist)[n.ranefs]){ # highest level
      form <- paste(fixed[2], fixed[1], fixed[3], "|", level)
      
      # Use lmList
      g.list <- suppressWarnings(lme4::lmList(formula(form), data = data))
     
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
        if(sum(class(g.vars[level][[1]]) != class(return.tbl[level][[1]])) != 0){
          g.vars[level][[1]] <- as.character(g.vars[level][[1]])
          return.tbl[level][[1]] <- as.character(return.tbl[level][[1]])
        }
        return.tbl <- tibble::tibble(
          unique(suppressMessages(dplyr::left_join(g.vars, return.tbl))))
      
      } else {
        return.tbl <- tibble::tibble(groups, eb.resid)
        names(return.tbl) <- c(level, eb.names)
        if(sum(class(g.vars[level][[1]]) != class(return.tbl[level][[1]])) != 0){
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
      g.list <- suppressWarnings(lme4::lmList(formula(form), data = data))
      
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
        dplyr::select(ncol(g.vars), 1:(ncol(g.vars)-1))
      
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
hlm_resid.lme <- function(object, level = 1, standardize = FALSE, include.ls = TRUE, data = NULL, ...) {
  
  if(!isNestedModel(object)){
    stop("hlm_resid is not currently implemented for non-nested models.")
  }
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
      ls.resid <- LSresids(object, level = 1, stand = standardize)
      ls.resid <- ls.resid[order(as.numeric(rownames(ls.resid))),]
    }
    
    # EB Residuals
    if(standardize == TRUE) {
      eb.resid <- data.frame(.std.resid = resid(object, type = "normalized"))
    } else { 
      eb.resid <- data.frame(.resid = resid(object, type = "response"))
    }

    # EB Fitted
    eb.fitted <- data.frame(.fitted = fitted(object))
    
    # Marginal Residuals
    mr <- resid(object, type="response", level=0)
    if (standardize == TRUE) {
      V      <- Matrix(extract_design(object)$V)
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
      dplyr::mutate(dplyr::across(where(is.character), ~ as.factor(.x))) %>%
      as.data.frame()
    model.data <- model.frame(formula(dataform), data)
    
    # NA action
    if(class(object$na.action) == "exclude"){ # if na.exclude
      # fix data frame
      na.index <- which(!rownames(data) %in% rownames(model.data))
      na.fix.data <- data[which(rownames(data) %in% na.index),] %>% 
        dplyr::select(all_of(names(model.data)))
      model.data <- rbind(model.data, na.fix.data)
      model.data <- model.data[order(as.numeric(rownames(model.data))),]
      
      # fix ls residuals
      if(include.ls == TRUE){
        na.fix.ls <- data.frame(LSR = rep(NA, length(na.index)),
                                LSF = rep(NA, length(na.index)))
        rownames(na.fix.ls) <- na.index
        names(na.fix.ls) <- c(names(ls.resid))
        ls.resid <- rbind(ls.resid, na.fix.ls)
        ls.resid <- ls.resid[order(as.numeric(rownames(ls.resid))),]
      }
    }
    
    # Continue to Assemble Tibble  
    if (include.ls == TRUE) {
      return.tbl <- tibble::tibble("id" = as.numeric(rownames(model.data)),
                                   model.data,
                                   eb.resid,
                                   eb.fitted,
                                   ls.resid,
                                   mar.resid,
                                   mar.fitted)
    } else { 
      return.tbl <- tibble::tibble("id" = as.numeric(rownames(model.data)),
                                   model.data,
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
      eb.resid <- nlme::ranef(object, standard = standardize)[[level]]
    } else { 
      eb.resid <- nlme::ranef(object, standard = standardize)
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
      ls.resid <- LSresids(object, level = level, stand = standardize)
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
      dplyr::mutate(dplyr::across(where(is.character), ~ as.factor(.x))) %>%
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
        if(sum(class(g.vars[level][[1]]) != class(return.tbl[level][[1]])) != 0){
          g.vars[level][[1]] <- as.character(g.vars[level][[1]])
          return.tbl[level][[1]] <- as.character(return.tbl[level][[1]])
        }
        return.tbl <- tibble::tibble(
          unique(suppressMessages(dplyr::left_join(g.vars, return.tbl))))
        
      } else {
        return.tbl <- tibble::tibble(groups, eb.resid)
        names(return.tbl) <- c(level, eb.names)
        if(sum(class(g.vars[level][[1]]) != class(return.tbl[level][[1]])) != 0){
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
        dplyr::select(ncol(g.vars), 1:(ncol(g.vars)-1))
      
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
