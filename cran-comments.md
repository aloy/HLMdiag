## Test environments
* local OS X install, R 4.0.3
* local OS X install, R-devel
* Ubuntu 20.04.1 LTS (on GitHub), R 4.0.3
* Ubuntu 20.04.1 LTS (on GitHub), R-devel
* win-builder (devel and 4.0.3)

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs in R 4.0.3 or R-devel on 
  OS X or win-builder.

There was one NOTE on Ubuntu (relase and devel):   
   installed size is  7.5Mb
     sub-directories of 1Mb or more:
       libs   6.5Mb

   This appears to be caused by the use of templated C++ linear 
   algebra library Armadillo

## Downstream dependencies

I have also run R CMD check on downstream dependencies of HLMdiag.
All packages passed.
