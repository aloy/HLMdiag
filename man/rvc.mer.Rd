\name{rvc.mer}
\alias{rvc}
\alias{rvc.mer}
\title{Relative variance change for mixed/hierarchical linear models}
\usage{
  \method{rvc}{mer} (model, group = NULL, delete = NULL)
}
\arguments{
  \item{model}{fitted model of class \code{mer}}

  \item{group}{variable used to define the group for which
  cases will be deleted.  If \code{group = NULL}, then
  individual cases will be deleted.}

  \item{delete}{index of individual cases to be deleted.
  For higher level units specified in this manner, the
  \code{group} parameter must also be specified.  If
  \code{case = NULL} then all cases are iteratively
  deleted.}
}
\value{
  If \code{delete = NULL} a matrix with columns
  corresponding to the variance components of the model and
  rows corresponding to the deleted observation/group is
  returned.

  If \code{delete} is specified then a named vector is
  returned.

  The residual variance is named \code{sigma2} and the
  other variance componenets are named \code{D**} where the
  trailing digits give the position in the covariance
  matrix of the random effects.
}
\description{
  This function calculates the relative variance change
  (RVC) for mixed/hierarchical linear models fit via
  \code{lmer}.
}
\author{
  Adam Loy \email{aloy@iastate.edu}
}
\references{
  Dillane, D. (2005). ``Deletion Diagnostics for the Linear
  Mixed Model.'' Ph.D. thesis, Trinity College Dublin
}
\seealso{
  \code{\link{leverage.mer}},
  \code{\link{cooks.distance.mer}},
  \code{\link{mdffits.mer}}, \code{\link{covratio.mer}},
  \code{\link{covtrace.mer}}
}
\keyword{models}
\keyword{regression}

