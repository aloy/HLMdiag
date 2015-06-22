## Test environments
* local OS X install, R 3.1.2
* win-builder (devel and release)

## R CMD check results

On my local OS X install and the release version of windows
there were no ERRORs, WARNINGs, or NOTEs. 

On the development version there were 2 ERRORs, but these
seem to stem from an issue with lme4 (sigma not being exported); 
however, sigma is exported, so there may be some instability
with the development version.
