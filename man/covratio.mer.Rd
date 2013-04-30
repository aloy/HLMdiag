\name{covratio.default}
\alias{covratio}
\alias{covratio.default}
\alias{covratio.lmerMod}
\alias{covratio.mer}
\alias{covtrace}
\alias{covtrace.default}
\alias{covtrace.lmerMod}
\alias{covtrace.mer}
\title{Influence on precision of fixed effects in HLMs}
\usage{
  \method{covratio}{default} (object, ...)

  \method{covtrace}{default} (object, ...)

  \method{covratio}{mer} (object, group = NULL,
    delete = NULL, ...)

  \method{covratio}{lmerMod} (object, group = NULL,
    delete = NULL, ...)

  \method{covtrace}{mer} (object, group = NULL,
    delete = NULL, ...)

  \method{covtrace}{lmerMod} (object, group = NULL,
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
  If \code{delete = NULL} then a vector corresponding to
  each deleted observation/group is returned.

  If \code{delete} is specified then a single value is
  returned corresponding to the deleted subset specified.
}
\description{
  These functions calculate measures of the change in the
  covariance matrices for the fixed effects based on the
  deletetion of an observation, or group of observations,
  for a hierarchical linear model fit using \code{lmer}.
}
\details{
  Both the covariance ratio (\code{covratio}) and the
  covariance trace (\code{covtrace}) measure the change in
  the covariance matrix of the fixed effects based on the
  deletion of a subset of observations.  The key difference
  is how the variance covariance matrices are compared:
  \code{covratio} compares the ratio of the determinants
  while \code{covtrace} compares the trace of the ratio.
}
\examples{
#' data(sleepstudy, package = 'lme4')
ss <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

# covratio for individual observations
ss.cr1 <- covratio(ss)

# covratio for subject-level deletion
ss.cr2 <- covratio(ss, group = "Subject")

\dontrun{
## A larger example
data(Exam, package = 'mlmRev')
fm <- lmer(normexam ~ standLRT * schavg + (standLRT | school), Exam)

# covratio for individual observations
cr1 <- covratio(fm)

# covratio for school-level deletion
cr2 <- covratio(fm, group = "school")
}
# covtrace for individual observations
ss.ct1 <- covtrace(ss)

# covtrace for subject-level deletion
ss.ct2 <- covtrace(ss, group = "Subject")

\dontrun{
## Returning to the larger example
# covtrace for individual observations
ct1 <- covtrace(fm)

# covtrace for school-level deletion
ct2 <- covtrace(fm, group = "school")
}
}
\author{
  Adam Loy \email{aloy@iastate.edu}
}
\references{
  Christensen, R., Pearson, L., & Johnson, W. (1992)
  Case-deletion diagnostics for mixed models.
  \emph{Technometrics}, \bold{34}(1), 38--45.

  Schabenberger, O. (2004) Mixed Model Influence
  Diagnostics, in \emph{Proceedings of the Twenty-Ninth SAS
  Users Group International Conference}, SAS Users Group
  International.
}
\seealso{
  \code{\link{leverage.mer}},
  \code{\link{cooks.distance.mer}}
  \code{\link{mdffits.mer}}, \code{\link{rvc.mer}}
}
\keyword{models}
\keyword{regression}

