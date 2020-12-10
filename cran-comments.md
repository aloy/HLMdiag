## Test environments
* local OS X install, R 4.0.3
* local OS X install, R-devel
* Ubuntu 16.04.6 LTS (on travis-ci), R 4.0.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs in R 4.0.3 or R-devel on OS X or win-builder.

There was one NOTE on Ubuntu:   
sub-directories of 1Mb or more:
    libs   5.2Mb

This appears to be caused by the use of C++ code

## Downstream dependencies

I have also run R CMD check on downstream dependencies of HLMdiag.
All packages passed.
