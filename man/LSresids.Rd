\name{LSresids}
\alias{LSresids}
\title{Calculating case-level (level-1) least squares residuals}
\usage{
  LSresids(object, level, sim = NULL, standardize = NULL)
}
\arguments{
  \item{object}{an object of class \code{mer}.}

  \item{level}{which residuals should be plotted: 1 for
  within-group residuals or the name of the grouping factor
  (as defined in \code{flist} of the \code{mer} object) for
  between-group residuals.}

  \item{sim}{optional argument giving the data frame used
  for LS residuals. This is used mainly for when dealing
  with simulations.}

  \item{standardize}{if \code{TRUE} the standardized
  level-1 residuals will also be returned (if \code{level =
  1}); if \code{"semi"} then the semi-standardized level-1
  residuals will be returned.}
}
\description{
  This function calculates the case-level least squares
  residuals found by fitting separate LS regression models
  to each case.
}
\author{
  Adam Loy \email{aloy@istate.edu}
}
\keyword{models}
\keyword{regression}

