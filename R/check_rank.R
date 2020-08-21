#' @export
check_rank <- function(object, ...){
  UseMethod("check_rank", object)
}

#' @export
#' @rdname check_rank.lmerMod
#' @method check_rank default
check_rank.default <- function(object, ...){
  stop(paste("there is no check_rank() method for objects of class",
             paste(class(object), collapse=", ")))
}

#' Check Rank of Groups in HLM
#'
#' \code{check_rank} takes a hierarchical linear model fit as a \code{lmerMod}
#' or \code{lme} object and checks the rank of all lowest level groups, as they would
#' be used in calculating LS residuals. Because the LS residuals rely on rank
#' sufficient groups, a large number of rank deficient groups likely indicates
#' that the LS residuals are unreliable. By defauly, returns a named vector
#' containing a numeric of the rank of each group, but if \code{logical = TRUE},
#' returns TRUE/FALSE if a group is rank deficient.
#' @param object an object of class \code{lmerMod} or \code{lme}.
#' @param logical By default, \code{logical = FALSE}, and the method returns the
#'   numeric rank of each group. If \code{lgl = TRUE},the method returns a
#'   TRUE/FALSE indicating if a given group is rank deficient.
#' @param ... do not use.
#' @seealso \link[HLMdiag]{hlm_resid}
#'
#' @export
#' @method check_rank lmerMod
#' @aliases check_rank

check_rank.lmerMod <- function(object, logical = FALSE, ...) {
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
  
  if(logical == TRUE) {
    max.rank <- max(ls.rank)
    ls.rank.lgl <- ls.rank == max.rank
    return(ls.rank.lgl)
    
  } else {
    return(ls.rank)
  }
  
}


#' @export
#' @rdname check_rank.lmerMod
#' @method check_rank lme
check_rank.lme <- function(object, logical = FALSE, ...) {
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
  
  if(logical == TRUE) {
    max.rank <- max(ls.rank)
    ls.rank.lgl <- ls.rank == max.rank
    return(ls.rank.lgl)
    
  } else {
    return(ls.rank)
  }
}
