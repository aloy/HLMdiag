# HLMdiag 0.5.1

* Added a `NEWS.md` file to track changes to the package.
* `cooks.distance()` and `mdffits()` with `include.attr = TRUE` now returns a tibble with column names for the fixed effects in agreement with the fixef() output.
* First columns of `hlm_resid()` output is now `.id` to avoid conflicts with `id` columns

# Version 0.5.0

## NEW FEATURES

Added separate functions to calculate residuals for LMEs. 

* `resid_marginal()` calculates (raw, studentized, Pearson, or Cholesky) marginal residuals.
* `resid_conditional` calculates (raw, studentized, Pearson, or Cholesky) conditional residuals (i.e., error terms).
* `resid_ranef` calculates raw and standardized predicted random effects.

    
## USER-VISIBLE CHANGES

The following functions are now defunct

* `HLMresid`
* `diagnostics`
* `group_qqnorm`
* `ggplot_qqnorm`


# Version 0.4.0

## DEVELOPMENT

* Influence diagnostics are now avaliable through the `case_delete` function and accompaning functions for `lme` objects from the `nlme` package. These diagnostics are also avaliable for three-level models. 
* Residual diagnostics are also now available via `hlm_resid` for `lme` model objects and three-level models.

## USER-VISIBLE CHANGES

* The `group` parameter for influence diagnostics has been changed to `level` in order to match the residual functions. 
`level` defaults to NULL, which will delete individual observations iteratively as `group = NULL` did. 
* The `cooks.distance` and `mdffits` functions now only return the values as a numeric vector, instead of also returning the 
    beta values as attributes. If these attributes are desired, the user can now set `include.attr = TRUE`, and a tibble
    will be returned instead with the influence diagnostics and the beta attributes. 
* The `dotplot_diag` function has been updated to be more efficient. Additionally, it no longer places labels on the y-axis 
    and only labels the top five observations in order to improve visibility.
* `LSresid` now returns only the residual values, excluding the model frame
* The `sim` argument for `LSresid` has been removed

## NEW FEATURES

* The `hlm_influence` function has been added. This function returns influence diagnostics appended to the model frame.
* The `hlm_resid` function has been added. This function returns residual diagnostics appended to the model frame.
* The `hlm_augment` function has been added, which combines `hlm_influence` and `hlm_resid` to return influence diagnostics and residuals appended to the model frame. 
* The `pull_resid` function has been added. This funciton returns a vector of a specified type of residual prioritizing computational efficiency.
* The `delete` parameter in `case_delete`, `hlm_influence`, and `hlm_augment` now also accepts character vectors at the second or third level. Observations or groups to be deleted can be specified by row indices in a numeric vector (as previously), or as character vectors of group level names found in `model@flist` (lmerMod models) or `model$groups` (lme models). 

## BUG FIXES

* `cooks.distance` and `mdffits` functions were fixed to solve an issue with the number of columns. 
* Fixed an issue with `case_delete` so that it works with three level models. 
* `hlm_resid`, `hlm_influence`, and `hlm_augment` properly respect `na.action` and work with models fit in `nlme`
* `LSresids` doesn't break with three-level models, or with models containing tranformed variables


# Version 0.3.1


## BUG FIXES
* `.extractV.lme` (and thus `.lme_matrices`) was fixed to work with more complex covariance structures fit via nlme.
* Updated package to work with the most recent version of `ggplot2`

# Version 0.3.0

## DEVELOPMENT

* Influence diagnostics in HLMdiag 0.3.0 are available for two-level models fit
    using the `lmer` function in lme4` or the `lme` function in `nlme`. I am still
    working to implement these methods for higher-level models using `lme`.
* HLMdiag no longer loads lme4 automatically (see above for the reason).

## NEW FEATURES

* The `rotate_ranef` function has been added. This function rotates the random 
    effects in an effort to find the least confounded residuals for distributional 
    assessment.

## BUG FIXES

* `LSresids` was fixed for an issue with the order of the resulting data frame.
* `case_delete` was fixed so that numeric group labels work properly, which fixes an issue
with `rvc`.
* Fixed an issue with `group_qqnorm`, by using `ppoints` rather than `.SampleQuantiles`.
* Fixed an issue with `case_delete.lmerMod`, to use the `getME()` function to extract `n`.
* A bug in the calculation of the Cholesky residuals was fixed (thanks to Harry Hiemstra for reporting the bug and the fix)


# Version 0.2.5

* Fixed a compatibility issue with Rcpp

# Version 0.2.4

* Added citation for the JSS paper
* Fixed a bug with the calculation of Cook's distance

# Version 0.2.3

* Added a function to calculate rotated random effects
* Added new data sets
* Maintenance for compatibility with lme4 1.0 and R 3.0.2

# Version 0.2.2 

* Changed the standardization of the EB level-1 residuals in `HLMresid` 
    to a more appropriate definition: e / var(\hat{e}).
* Added functions to add compatibility with the development version of lme4.


# Version 0.2.1 

* Fixed a bug in `group_qqnorm`
* Checked compatibility with R 2.15.3


# Version 0.2.0 

## DEVELOPMENT

* Influence diagnostics in HLMdiag 0.2.0 are compatible with hierarchical 
(multilevel) linear models of any size and with models with crossed factors.
* HLMdiag 0.2.0 offers significantly faster computation of the deletion
    diagnostics for fixed effects that are based on one-step approximations.
* S3 methods have been created for `cooks.distance`, `mdffits`, `covratio`,
`covtrace`, `rvc`, and `leverage` for objects of class `mer`.
* Full deletions are still available using the `case_delete` function, with
    corresponding S3 methods for objects of class `case_delete`.

## NEW FEATURES

* A `leverage` function has been added.
* `case_delete` and other deletion functions now allow for the user to 
    manually specify a subset to delete.

## OTHER USER-VISIBLE CHANGES

* changes to the arguments of `dotplot_diag` to accomodate a more general
    usage.
* `diagnostics` no longer requires a `model` parameter to be specified.
* A `delete` parameter has been added to `case_delete` to allow for manual
    specification of a subset that should be deleted.
* For observation-level deletion using `case_delete` the user should specify
    `group = NULL` rather than `group = FALSE`.


# Version 0.1.6 

* Updates for compatibility with ggplot2 >= 0.9.2
* Updates to NAMESPACE

# Version 0.1.5 

* Updates to NAMESPACE to fix compatibility issues

# Version 0.1.4 

* Updates for compatibility with ggplot2 0.9.0

# Version 0.1.3 

* Added "marginal" residuals to the type argument for HLMresid.

# Version 0.1.2 


## USER-VISIBLE CHANGES:

* Removed the formula argument from LSresids. The formula is now obtained automatically from the mer object. Note: we are still working on automatic recognition of math operators such as log(), but anything in I() is recognized.
* Added the function HLMresid, a wrapper that will extract both the LS and EB residuals given an mer object.
* Added a level argument to LSresids, so the function can extract LS residuals from either level of the model.

## BUG FIXES: 

* Fixed the ordering of output from 'LSresids' to match the mode frame obtained from the mer object.

# Version 0.1.1 

* Created 'adjust_lmList' class to handle fitting separate linear models when a factor is constant across the group.
* 'random_ls_coef' was removed and replaced by the 'coef' method for 'adjust_lmList' objects.
* Added 'rvc' diagnostic.
* Extended 'case_delete' from only handling deletion for fixed effects to also handling deletion for variance components.
* Improved 'dotplot_diag' to handle modified dotplots.
