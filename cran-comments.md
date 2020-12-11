## Test environments
* local OS X install, R 4.0.3
* local OS X install, R-devel
* Ubuntu 16.04.6 LTS (on travis-ci), R 4.0.2
* Ubuntu Linux 16.04 LTS, GCC (r-hub), R-release
* Fedora Linux, R-devel, clang, gfortran (r-hub)
* win-builder (devel and 4.0.3)

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs in R 4.0.3 or R-devel on OS X or win-builder.

There was one NOTE on Ubuntu:   
   installed size is  6.1Mb
     sub-directories of 1Mb or more:
       libs   5.2Mb

   This appears to be caused by the use of C++ code

## Downstream dependencies

I have also run R CMD check on downstream dependencies of HLMdiag.
All packages passed.
