\name{mdffits.mer}
\alias{mdffits}
\alias{mdffits.mer}
\title{MDFFITS for mixed/hierarchical linear models}
\usage{
  \method{mdffits}{mer} (object, group = NULL,
    delete = NULL, ...)
}
\arguments{
  \item{object}{fitted object of class \code{mer}}

  \item{group}{variable used to define the group for which
  cases will be deleted.  If \code{group = NULL}, then
  individual cases will be deleted.}

  \item{delete}{index of individual cases to be deleted.
  For higher level units specified in this manner, the
  \code{group} parameter must also be specified.  If
  \code{case = NULL} then all cases are iteratively
  deleted.}

  \item{...}{do not use}
}
\description{
  This function calculates multivariate DFFITS for the
  fixed effects parameters for a mixed/hierarchical model
  fit by \code{lmer}.
}
\author{
  Adam Loy \email{aloy@iastate.edu}
}
\references{
  Schabenberger, O. (2004),``Mixed Model Influence
  Diagnostics,'' in \emph{Proceedings of the Twenty-Ninth
  SAS Users Group International Conference}, SAS Users
  Group International.
}
\seealso{
  \code{\link{leverage.mer}}, \code{\link{mdffits.mer}},
  \code{\link{covratio.mer}}, \code{\link{covtrace.mer}},
  \code{\link{rvc.mer}}
}
\keyword{models}
\keyword{regression}

