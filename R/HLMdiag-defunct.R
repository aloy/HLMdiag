#' Defunct functions in package HLMdiag
#' 
#' These functions are defunct and no longer available.
#' 
#' @param ... arguments passed to defunct functions
#' 
#' @details 
#' \code{HLMresid} is replaced by \code{\link{hlm_resid}}
#' 
#' \code{diagnostics} is replaced by \code{\link{hlm_influence}}
#' 
#' \code{group_qqnorm} and \code{group_qqnorm} are replaced by functions in \pkg{qqplotr}. 
#'    See \code{\link[qqplotr]{stat_qq_point}}, \code{\link[qqplotr]{stat_qq_line}}, and 
#'    \code{\link[qqplotr]{stat_qq_band}}.
#' 
#' @name HLMdiag-defunct
NULL

#' @rdname HLMdiag-defunct
HLMresid <- function (...) {
  .Defunct("hlm_resid", package = "HLMresid")
}

#' @rdname HLMdiag-defunct
HLMresid.default <- function (...) {
  .Defunct("hlm_resid", package = "HLMresid")
}

#' @rdname HLMdiag-defunct
HLMresid.lmerMod <- function (...) {
  .Defunct("hlm_resid", package = "HLMresid")
}

#' @rdname HLMdiag-defunct
HLMresid.mer <- function (...) {
  .Defunct("hlm_resid", package = "HLMresid")
}

#' @rdname HLMdiag-defunct
diagnostics <- function (...) {
  .Defunct("hlm_influence", package = "HLMresid")
}

#' @rdname HLMdiag-defunct
group_qqnorm <- function(...) {
  .Defunct("stat_qq_point", package = "qqplotr")
}

#' @rdname HLMdiag-defunct
ggplot_qqnorm <- function(...) {
  .Defunct("stat_qq_point", package = "qqplotr")
}