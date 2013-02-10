\name{case_delete}
\alias{case_delete}
\title{Case Deletion for \code{mer} objects}
\usage{
  case_delete(model, group = NULL,
    type = c("both", "fixef", "varcomp"), delete = NULL)
}
\arguments{
  \item{model}{the original hierarchical model fit using
  \code{lmer()}}

  \item{group}{a variable used to define the group for
  which cases will be deleted.  If this is left \code{NULL}
  (default), then the function will delete individual
  observations.}

  \item{type}{the part of the model for which you are
  obtaining deletion diagnostics: the fixed effects
  (\code{"fixef"}), variance components (\code{"varcomp"}),
  or \code{"both"} (default).}

  \item{delete}{index of individual cases to be deleted.
  For higher level units specified in this manner, the
  \code{group} parameter must also be specified.  If
  \code{delete = NULL} then all cases are iteratively
  deleted.}
}
\value{
  a list with the following compontents: \describe{
  \item{\code{fixef.original}}{the original fixed effects
  estimates} \item{\code{ranef.original}}{the original
  predicted random effects} \item{\code{vcov.original}}{the
  original variance-covariance matrix for the fixed
  effects} \item{\code{varcomp.original}}{the original
  estimated variance components}
  \item{\code{fixef.delete}}{a list of the fixed effects
  estimated after case deletion}
  \item{\code{ranef.delete}}{a list of the random effects
  predicted after case deletion}
  \item{\code{vcov.delete}}{a list of the
  variance-covariance matrices for the fixed effects
  obtained after case deletion}
  \item{\code{fitted.delete}}{a list of the fitted values
  obtained after case deletion}
  \item{\code{varcomp.delete}}{a list of the estimated
  variance components obtained after case deletion} }
}
\description{
  This function is used to iteratively delete groups
  corresponding to the levels of a hierarchical linear
  model. It uses \code{lmer()} to fit the models for each
  deleted case (i.e. uses brute force). To investigate
  numerous levels of the model, the function will need to
  be called multiple times, specifying the group (level) of
  interest each time.
}
\examples{
data(sleepstudy, package = 'lme4')
fm <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)

# Deleting every Subject
fmDel <- case_delete(model = fm, group = "Subject", type = "both")

# Deleting only subject 308
del308 <- case_delete(model = fm, group = "Subject", type = "both", delete = 308)

# Deleting a subset of subjects
delSubset <- case_delete(model = fm, group = "Subject", type = "both", delete = 308:310)
}
\author{
  Adam Loy \email{aloy@istate.edu}
}
\references{
  Christensen, R., Pearson, L.M., and Johnson, W. (1992)
  Case-Deletion Diagnostics for Mixed Models,
  \emph{Technometrics}, \bold{34}, 38 -- 45.

  Schabenberger, O. (2004) Mixed Model Influence
  Diagnostics, in \emph{Proceedings of the Twenty-Ninth SAS
  Users Group International Conference}, SAS Users Group
  International.
}
\keyword{models}
\keyword{regression}

