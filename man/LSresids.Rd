\name{LSresids}
\alias{LSresids}
\title{Calculating least squares residuals}
\usage{
  LSresids(object, level, sim = NULL, standardize = FALSE)
}
\arguments{
  \item{object}{an object of class \code{mer}.}

  \item{level}{which residuals should be extracted: 1 for
  case-level residuals or the name of a grouping factor (as
  defined in \code{flist} of the \code{mer} object) for
  between-group residuals.}

  \item{sim}{optional argument giving the data frame used
  for LS residuals. This is used mainly when dealing with
  simulations.}

  \item{standardize}{if \code{TRUE} the standardized
  level-1 residuals will also be returned (if \code{level =
  1}); if \code{"semi"} then the semi-standardized level-1
  residuals will be returned.}
}
\description{
  This function calculates least squares (LS) residuals
  found by fitting separate LS regression models to each
  case. For examples see the documentation for
  \code{HLMresid}.
}
\author{
  Adam Loy \email{aloy@istate.edu}
}
\references{
  Hilden-Minton, J. (1995) Multilevel diagnostics for mixed
  and hierarchical linear models. University of California
  Los Angeles.
}
\seealso{
  \code{\link{HLMresid}}
}
\keyword{models}
\keyword{regression}

