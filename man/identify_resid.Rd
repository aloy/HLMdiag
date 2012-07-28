\name{identify_resid}
\alias{identify_resid}
\title{Identifying unusual points}
\usage{
identify_resid(eb, ols, identify)
}
\arguments{
  \item{formula}{a formula that can be used with
  \code{lm()}}

  \item{identify}{the percentage of points to identify as
  unusual}
}
\description{
  This function will identify the unusual points found when
  looking at the plot comparing the shrinkage estimates and
  the random coefficients found through OLS.
}
\author{
  Adam Loy \email{aloy@istate.edu}
}
