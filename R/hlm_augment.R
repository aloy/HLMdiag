#' @export
hlm_augment <- function(object, ...){
  UseMethod("hlm_augment", object)
}

#'Calculating residuals and influence diagnostics for HLMs
#'
#' @export
#' @rdname hlm_augment.lmerMod
#' @method hlm_augment default
#' @S3method hlm_augment default
hlm_augment.default <- function(object, ...){
  stop(paste("there is no hlm_augment() method for objects of class",
             paste(class(object), collapse=", ")))
}

#hlm_augment.lmerMod <- function(object, level = 1, )
