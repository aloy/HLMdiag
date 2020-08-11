## usethis namespace: start
#' @useDynLib HLMdiag, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("HLMdiag", libpath)
}