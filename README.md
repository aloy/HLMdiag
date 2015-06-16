# HLMdiag [![Build Status](https://travis-ci.org/aloy/HLMdiag.svg?branch=master)](https://travis-ci.org/aloy/HLMdiag) [![](http://cranlogs.r-pkg.org/badges/HLMdiag)](http://cran.rstudio.com/web/packages/HLMdiag/index.html)

Up to now diagnostics for mixed and hierarchical models have required much programming by 
the analyst, especially if one desires influence diagnostics. 
To help fill this need, `HLMdiag`:

* Provides convenience functions for residual analysis.
  * Allows the analyst to obtain residuals estimated by least squares (LS) or empirical Bayes (EB).
  * Allows the analyst to obtain different residual quantities (e.g. marginal, conditional, BLUPs for the two-level model).

* Implements influence analysis.
  * Leverage
  * Deletion diagnostics -- Cook's distance, MDFFITS, covariance ratio & trace

`HLMdiag` strives to provide an easy to use interface for models fit using `lmer` from the package `lme4` that is draws from the ideas of `influence.measures` for regression diagnostics.

Development may be quite slow for a while, as I am finishing my degree, so please bear with me!

## Development version

If you would like to download the development version of `HLMdiag`, I would recommend using Hadley Wickham's `devtools` package:

    # install.packages("devtools")
    library(devtools)
    install_github("aloy/HLMdiag")
