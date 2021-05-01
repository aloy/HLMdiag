## Test environments
* local OS X install, R 4.0.5
* Windows Server 2008 R2 SP1, R-devel (R-hub), 32/64 bit
* Fedora Linux, R-devel, clang, gfortran (R-hub)
* Ubuntu 20.04 LTS (on GitHub), R 4.0.5 and devel
* win-builder (devel and 4.0.5)

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs in R 4.0.5 or R-devel on 
  OS X, Windows server on R-hub, Fedora, or win-builder.

There was one NOTE on Ubuntu (relase and devel):   
    installed size is  7.6Mb
    sub-directories of 1Mb or more:
      libs   6.7Mb

   This appears to be caused by the use of templated C++ linear 
   algebra library Armadillo

## Downstream dependencies

I have also run R CMD check on downstream dependencies of HLMdiag.
All packages passed. (https://github.com/aloy/HLMdiag/tree/master/revdep).
