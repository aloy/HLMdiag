\name{mdffits.default}
\alias{cooks.distance}
\alias{cooks.distance.mer}
\alias{mdffits}
\alias{mdffits.default}
\alias{mdffits.mer}
\title{Influence on fixed effects of HLMs}
\usage{
  \method{mdffits}{default} (object, ...)

  \method{cooks.distance}{mer} (model, group = NULL,
    delete = NULL, ...)

  \method{mdffits}{mer} (object, group = NULL,
    delete = NULL, ...)
}
\arguments{
  \item{model}{fitted model of class \code{mer}}

  \item{group}{variable used to define the group for which
  cases will be deleted.  If \code{group = NULL}, then
  individual cases will be deleted.}

  \item{delete}{index of individual cases to be deleted. To
  delete specific observations the row number must be
  specified. To delete higher level units the group ID and
  \code{group} parameter must be specified. If \code{delete
  = NULL} then all cases are iteratively deleted.}

  \item{...}{do not use}

  \item{object}{fitted object of class \code{mer}}
}
\value{
  Both functions return a numeric vector (or single value
  if \code{delete} has been specified) with attribute
  \code{beta_cdd} giving the difference between the full
  and deleted parameter estimates.
}
\description{
  These functions calculate measures of the change in the
  fixed effects estimates based on the deletetion of an
  observation, or group of observations, for a hierarchical
  linear model fit using \code{lmer}.
}
\details{
  Both Cook's distance and MDFFITS measure the change in
  the fixed effects estimates based on the deletion of a
  subset of observations. The key difference between the
  two diagnostics is that Cook's distance uses the
  covariance matrix for the fixed effects from the original
  model while MDFFITS uses the covariance matrix from the
  deleted model.
}
\note{
  Because MDFFITS requires the calculation of the
  covariance matrix for the fixed effects for every model,
  it will be slower.
}
\examples{
data(sleepstudy, package = 'lme4')
ss <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

# Cook's distance for individual observations
ss.cd.lev1 <- cooks.distance(ss)

# Cook's distance for each Subject
ss.cd.subject <- cooks.distance(ss, group = "Subject")

\dontrun{
data(Exam, package = 'mlmRev')
fm <- lmer(normexam ~ standLRT * schavg + (standLRT | school), Exam)

# Cook's distance for individual observations
cd.lev1 <- cooks.distance(fm)

# Cook's distance for each school
cd.school <- cooks.distance(fm, group = "school")

# Cook's distance when school 1 is deleted
cd.school1 <- cooks.distance(fm, group = "school", delete = 1)

}
# MDFFITS  for individual observations
ss.m1 <- mdffits(ss)

# MDFFITS for each Subject
ss.m.subject <- mdffits(ss, group = "Subject")

\dontrun{

# MDFFITS  for individual observations
m1 <- mdffits(fm)

# MDFFITS for each school
m.school <- mdffits(fm, group = "school")
}
}
\author{
  Adam Loy \email{aloy@iastate.edu}
}
\references{
  Christensen, R., Pearson, L., & Johnson, W. (1992)
  Case-deletion diagnostics for mixed models.
  \emph{Technometrics}, \bold{34}, 38--45.

  Schabenberger, O. (2004) Mixed Model Influence
  Diagnostics, in \emph{Proceedings of the Twenty-Ninth SAS
  Users Group International Conference}, SAS Users Group
  International.
}
\seealso{
  \code{\link{leverage.mer}}, \code{\link{covratio.mer}},
  \code{\link{covtrace.mer}}, \code{\link{rvc.mer}}
}
\keyword{models}
\keyword{regression}

