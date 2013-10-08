\name{rvc.default}
\alias{rvc}
\alias{rvc.default}
\alias{rvc.lmerMod}
\alias{rvc.mer}
\title{Relative variance change for HLMs}
\usage{
  \method{rvc}{default} (object, ...)

  \method{rvc}{mer} (object, group = NULL, delete = NULL,
    ...)

  \method{rvc}{lmerMod} (object, group = NULL,
    delete = NULL, ...)
}
\arguments{
  \item{object}{fitted object of class \code{mer} or
  \code{lmerMod}}

  \item{group}{variable used to define the group for which
  cases will be deleted.  If \code{group = NULL}, then
  individual cases will be deleted.}

  \item{delete}{index of individual cases to be deleted. To
  delete specific observations the row number must be
  specified. To delete higher level units the group ID and
  \code{group} parameter must be specified. If \code{delete
  = NULL} then all cases are iteratively deleted.}

  \item{...}{do not use}
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
  (RVC) of hierarchical linear models fit via \code{lmer}.
}
\author{
  Adam Loy \email{loyad01@gmail.com}
}
\references{
  Dillane, D. (2005) Deletion Diagnostics for the Linear
  Mixed Model. Ph.D. thesis, Trinity College Dublin
}
\seealso{
  \code{\link{leverage.mer}},
  \code{\link{cooks.distance.mer}},
  \code{\link{mdffits.mer}}, \code{\link{covratio.mer}},
  \code{\link{covtrace.mer}}
}
\keyword{models}
\keyword{regression}

