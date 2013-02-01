\name{BlockZ}
\alias{BlockZ}
\title{Reorganizing Z matrix}
\usage{
  BlockZ(object)
}
\arguments{
  \item{object}{a fitted model object of class \code{mer}.}
}
\description{
  This function reorganizes the Z matrix obtained from an
  \code{mer} object into the block form discussed by
  Demidenko (2004). Currently, this function assumes there
  is only one level and that units are nested.
}
\author{
  Adam Loy \email{aloy@iastate.edu}
}
\references{
  Demidenko, E. (2004). Mixed Models: Theory and
  Applications.  New York, Wiley.
}

