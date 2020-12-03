## Test environments
* local OS X install, R 4.0.3
* local OS X install, R-devel
* Ubuntu 16.04.6 LTS (on travis-ci), R 4.0.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs in R 3.2.2 or R-devel.

## Resubmission
This is a resubmission. In this version I have 
added the following to the NAMESPACE to eliminate notes:

importFrom("grDevices", "devAskNewPage")
importFrom("methods", "as", "is")
importFrom("stats", "IQR", "aggregate", "complete.cases", "fitted",
           "formula", "lm", "lm.influence", "model.frame",
           "model.matrix", "ppoints", "qnorm", "qt", "quantile",
           "reorder", "resid", "rstandard", "varimax", "vcov",
           "cooks.distance", "covratio", "coef", "confint")
