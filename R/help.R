#' Diagnostic tools for hierarchical (multilevel) linear models
#' 
#' HLMdiag provides a suite of diagnostic tools for hierarchical 
#' (multilevel) linear models fit using the \code{lme4} or \code{nlme}
#' packages. These tools are grouped below by purpose. 
#' See the help documentation for additional information
#' about each function.
#' 
#' \bold{Residual analysis}
#' 
#' HLMdiag's \code{\link{hlm_resid}} function provides a wrapper that 
#' extracts residuals and fitted values for individual observations
#' or groups of observations. In addition to being a wrapper function for functions
#' implemented in the \code{lme4} and \code{nlme} packages,
#' \code{\link{hlm_resid}} provides access to the marginal and least squares
#' residuals.
#' 
#' \bold{Influence analysis}
#' 
#' HLMdiag's \code{\link{hlm_influence}} function provides a convenient wrapper 
#' to obtain influence diagnostics for each observation or group of observations 
#' appended to the data used to fit the model. The diagnostics returned by 
#' \code{\link{hlm_influence}} include Cook's distance, MDFFITS, covariance trace (covtrace),
#' covariance ratio (covratio), leverage, and relative variance change (RVC). 
#' HLMdiag also contains functions to calculate these diagnostics individually, as discussed below. 
#' 
#' Influence on fitted values
#' 
#' HLMdiag provides \code{\link{leverage}} that calculates the influence
#' that observations/groups have on the fitted values (leverage). 
#' For mixed/hierarchical models leverage can be decomposed into two parts: the 
#' fixed part and the random part. We refer the user to the references
#' cited in the help documentation for additional explanation.
#' 
#' Influence on fixed effects estimates
#' 
#' HLMdiag provides \code{\link{cooks.distance}} and \code{\link{mdffits}}
#' to assess the influence of subsets of observations on the fixed effects.
#' 
#' Influence on precision of fixed effects
#' 
#' HLMdiag provides \code{\link{covratio}} and \code{\link{covtrace}}
#' to assess the influence of subsets of observations on the precision of
#' the fixed effects.
#' 
#' Influence on variance components
#' 
#' HLMdiag's \code{\link{rvc}} calculates the relative variance change to
#' assess the influence of subsets of observations on the variance
#' components.
#' 
#' \bold{Graphics}
#' 
#' HLMdiag also strives to make graphical assessment easier in the 
#' \code{ggplot2} framework by providing dotplots for influence diagnostics
#' (\code{\link{dotplot_diag}}), grouped Q-Q plots (\code{\link{group_qqnorm}}),
#' and Q-Q plots that combine the functionality of \code{\link{qqnorm}} and
#' \code{\link{qqline}} (\code{\link{ggplot_qqnorm}}).
#' 
#' @useDynLib HLMdiag, .registration = TRUE
#' @importFrom magrittr %>%
#' @importFrom reshape2 melt dcast
#' @importFrom plyr adply ddply
#' @importFrom MASS rlm
#' @importFrom mgcv tensor.prod.model.matrix
#' @importFrom dplyr select left_join mutate across bind_cols filter arrange desc
#' @importFrom stringr str_c str_detect str_split
#' @importFrom purrr map map_lgl map_df map_dfc map_dfr map_dbl
#' @importFrom tibble tibble
#' @importFrom tidyselect all_of
#' @importFrom janitor clean_names
#' @import Matrix
#' @import methods
#' @import ggplot2
#' @importFrom grDevices devAskNewPage
#' @importFrom stats coef confint IQR aggregate 
#'  complete.cases fitted formula lm lm.influence model.frame
#'  model.matrix ppoints qnorm qt quantile reorder resid rstandard
#'  varimax vcov cooks.distance covratio as.formula getCall na.exclude
#'  predict sigma
#' @name HLMdiag
#' @aliases HLMdiag package-HLMdiag
#' @keywords package
"_PACKAGE"