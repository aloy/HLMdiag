\name{covratio.mer}
\alias{covratio}
\alias{covratio.mer}
\title{COVRATIO for mixed/hierarchical linear models}
\usage{
  \method{covratio}{mer} (model, group = NULL,
    delete = NULL)
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
  If \code{delete = NULL} then a vector corresponding to
  each deleted observation/group is returned.

  If \code{delete} is specified then a single value is
  returned corresponding to the deleted subset specified.
}
\description{
  This function calculates COVRATIO for the fixed effects
  parameters for a mixed/hierarchical model fit by
  \code{lmer}.
}
\author{
  Adam Loy \email{aloy@iastate.edu}
}
\references{
  Christensen, R., Pearson, L., & Johnson, W. (1992).
  Case-deletion diagnostics for mixed models.
  \emph{Technometrics}, 34(1), 38–45.

  Schabenberger, O. (2004),``Mixed Model Influence
  Diagnostics,'' in \emph{Proceedings of the Twenty-Ninth
  SAS Users Group International Conference}, SAS Users
  Group International.
}
\seealso{
  \code{\link{leverage.mer}},
  \code{\link{cooks.distance.mer}}
  \code{\link{mdffits.mer}}, \code{\link{covtrace.mer}},
  \code{\link{rvc.mer}}
}
\keyword{models}
\keyword{regression}

