\name{ppnormBand}
\alias{ppnormBand}
\title{Adding simultaneous confidence bands to normal probability plots}
\usage{
  ppnormBand(x, e)
}
\arguments{
  \item{x}{a vector of probabilities}

  \item{e}{factors for EP simultaneous approx. confidence
  bands (see Nair 1984 or Meeker and Escobar 1998)}
}
\description{
  This function calculates simultaneous confidence bands
  for the estimated cdf using the logistic-transform
  normal-approximation method.
}
\details{
  Information on how to find the approximate factors
  \code{e} can be found in Meeker and Escobar (1998) on
  page 61. They provide a table which summarizes some
  common choices for the factor.
}
\author{
  Adam Loy \email{aloy@iastate.edu}
}
\references{
  Meeker, W. Q. and Escobar, L. A. (1998),
  \emph{Statistical Methods for Reliability Data}, New
  York, NY: Wiley.
}

