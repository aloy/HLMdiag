# HLMdiag <img src='man/figures/logo.png' align="right" height="139" />

[![R build
status](https://github.com/aloy/HLMdiag/workflows/R-CMD-check/badge.svg)](https://github.com/aloy/HLMdiag/actions)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/HLMdiag)](https://cran.r-project.org/package=HLMdiag)
![CRAN\_Downloads\_Badge](http://cranlogs.r-pkg.org/badges/HLMdiag)

The package `HLMdiag` was created in order to provide a unified
framework for analysts to diagnose hierarchical linear models (HLMs).
When `HLMdiag` was created in 2014, it made diagnostic procedures
available for HLMs that had not yet been implemented in statistical
software. Over the past 6 years, other packages have gradually
implemented some of these procedures; however, diagnosing a model still
often requires multiple packages, each of which has its own syntax and
quirks. `HLMdiag` provides diagnostic tools targeting all aspects and
levels of hierarchical linear models in a single package. `HLMdiag`
provides wrapper functions to all types of residuals implemented in
`lme4` and `nlme` as well as providing access to the marginal and least
squares residuals. For influence diagnostics, `HLMdiag` provides
functions to calculate Cook’s distance, MDFFITS, covariance trace and
ratio, relative variance change, and leverage.

## Installation

If you would like to install the development version of `HLMdiag`, you
may do so using `devtools`:

``` r
#install.packages("devtools")
library(devtools)
devtools::install_github("aloy/HLMdiag")
```

To instead download the stable CRAN version instead, use:

``` r
install.packages("HLMdiag")
```

## Details

The functions provided by this package can be separated into three
groups: residual analysis, influence analysis, and graphical tools.

## Residual Analysis

The residual functions in `HLMdiag` allow the analyst to estimate all
types of residuals defined for a hierarchical linear model. They provide
access to level-1, higher-level, and marginal residuals and use both
Least Squares and Empirical Bayes estimation methods to provide the
analyst with more choices in evaluating a model. The `hlm_resid` method
is inspired by the `augment()` function in `broom` and appends all types
of residuals and fitted values for a given level to the model frame;
however, individual types of residuals can be calculated with the
`pull_resid` method.

The functions available for residual analysis in `HLMdiag` are:
\*`hlm_resid()` calculates all residual diagnostics for a given level,
returning a tibble with the residuals and fitted values appended to the
original model frame.

\*`pull_resid()` calculates a specified type of residual, returning a
vector and prioritizing computational efficiency.

\*`hlm_augment()` combines `hlm_influence()` and `hlm_resid()` to return
a tibble with residual values and influence diagnostics appended to the
original model frame.

### Influence Analysis

The influence analysis functions provide functionality to calculate
Cook’s distance, MDFFITS, covariance ratio, covariance trace, relative
variance change, and leverage. Additionally, two functions to calculate
Cook’s distance, MDFFITS, covariance ratio, and covariance trace are
provided: a one step approximation, and a full refit method that refits
the model and recalculates the fixed and random effects. This
functionality is available through individual functions for each
diagnostic; however, the `hlm_influence` function can be used to
calculate all available diagnostics for each observation or group of
observations.

The functions available for influence analysis in `HLMdiag` are:

  - `cooks.distance()` calculates Cook’s distance values, which measures
    the difference between the original fixed effects and the deleted
    ones.
  - `mdffits()` calculates MDFFITS, a multivariate version of the DFFITS
    statistic, which is also a measure of the difference in fixed
    effects.
  - `covtrace()` calculates covariance trace, the ratio between the
    covariance matrices with and without unit *i* to the identity
    matrix.
  - `covratio()` calculate covariance ratio, a comparison of the two
    covariance matrices with and without unit *i* using their
    determinants.
  - `rvc()` calculates relative variance change, a measurement of the
    ratio of estimates of the *l* th variance component with and without
    unit *i*.
  - `leverage()` calculates leverage, he rate of change in the predicted
    response with respect to the observed response.
  - `case_delete()` iteratively deletes observations or groups of
    observations, returning a list of fixed and random components from
    the original model and the models created by deletion.
  - `hlm_influence()` calculates all of the influence diagnostics,
    returning a tibble with the influence values appended to the
    original model frame.
  - `hlm_augment()` combines `hlm_influence()` and `hlm_resid()` to
    return a tibble with residual values and influence diagnostics
    appended to the original model frame.

### Graphical Tools

`HLMdiag` provides the function `dotplot_diag()`, which creates dotplots
to visually represent influence diagnostics. It is especially useful
when used with the values returned by `hlm_influence()`. `HLMdiag` also
provides grouped Q-Q plots (`group_qqnorm()`), and Q-Q plots that
combine the functionality of qqnorm and qqline (`ggplot_qqnorm()`).

## Usage

### Residual Analysis

We will use the `sleepstudy` data set from the `lme4` package.

``` r
library(lme4)
#> Loading required package: Matrix
library(HLMdiag)
#> 
#> Attaching package: 'HLMdiag'
#> The following object is masked from 'package:stats':
#> 
#>     covratio
data(sleepstudy, package = "lme4")
sleep.lmer <- lme4::lmer(Reaction ~ Days + (Days|Subject), data = sleepstudy) 
```

We calculate the unstandardized level-1 and marginal residuals for each
observation below.

``` r
hlm_resid(sleep.lmer)
#> # A tibble: 180 x 10
#>       id Reaction  Days Subject  .resid .fitted .ls.resid .ls.fitted .mar.resid
#>    <dbl>    <dbl> <dbl> <fct>     <dbl>   <dbl>     <dbl>      <dbl>      <dbl>
#>  1     1     250.     0 308       -4.10    254.      5.37       244.      -1.85
#>  2     2     259.     1 308      -14.6     273.     -7.25       266.      -3.17
#>  3     3     251.     2 308      -42.2     293.    -36.9        288.     -21.5 
#>  4     4     321.     3 308        8.78    313.     12.0        309.      38.6 
#>  5     5     357.     4 308       24.5     332.     25.6        331.      63.6 
#>  6     6     415.     5 308       62.7     352.     61.7        353.     111.  
#>  7     7     382.     6 308       10.5     372.      7.42       375.      68.0 
#>  8     8     290.     7 308     -101.      391.   -106.         397.     -34.5 
#>  9     9     431.     8 308       19.6     411.     12.3        418.      95.4 
#> 10    10     466.     9 308       35.7     431.     26.3        440.     121.  
#> # … with 170 more rows, and 1 more variable: .mar.fitted <dbl>
```

For more information and examples of the functionality of `hlm_resid()`,
see the residual diagnostics vignette.

### Influence Analysis

We calculate influence diagnostics for each observation with the
following line:

``` r
hlm_influence(sleep.lmer)
#> # A tibble: 180 x 9
#>       id Reaction  Days Subject  cooksd mdffits covtrace covratio
#>    <int>    <dbl> <dbl> <fct>     <dbl>   <dbl>    <dbl>    <dbl>
#>  1     1     250.     0 308     1.48e-4 1.47e-4 0.00887      1.01
#>  2     2     259.     1 308     1.10e-3 1.09e-3 0.00558      1.01
#>  3     3     251.     2 308     5.13e-3 5.11e-3 0.00330      1.00
#>  4     4     321.     3 308     1.14e-4 1.14e-4 0.00175      1.00
#>  5     5     357.     4 308     3.93e-4 3.92e-4 0.000778     1.00
#>  6     6     415.     5 308     1.07e-3 1.07e-3 0.000321     1.00
#>  7     7     382.     6 308     3.49e-5 3.49e-5 0.000361     1.00
#>  8     8     290.     7 308     8.81e-3 8.80e-3 0.000944     1.00
#>  9     9     431.     8 308     8.23e-4 8.21e-4 0.00219      1.00
#> 10    10     466.     9 308     5.99e-3 5.96e-3 0.00435      1.00
#> # … with 170 more rows, and 1 more variable: leverage.overall <dbl>
```

For more information and examples of the functionality of
`hlm_influence()`, see the influence diagnostics vignette.
