\name{diagnostics}
\alias{cooks.distance.case_delete}
\alias{covratio.case_delete}
\alias{covtrace.case_delete}
\alias{diagnostics}
\alias{mdffits.case_delete}
\alias{rvc.case_delete}
\title{Calculating influence diagnostics for HLMs.}
\usage{
  diagnostics(object)

  \method{cooks.distance}{case_delete} (model, ...)

  \method{mdffits}{case_delete} (object, ...)

  \method{covtrace}{case_delete} (object, ...)

  \method{covratio}{case_delete} (object, ...)

  \method{rvc}{case_delete} (object, ...)
}
\arguments{
  \item{object}{an object containing the output returned by
  \code{case_delete()}}

  \item{model}{an object containing the output returned by
  \code{case_delete()}. This is only named differently to
  agree with the generic.}

  \item{...}{do not use}
}
\description{
  This group of functions is used to compute deletion
  diagnostics for a hierarchical linear model based on the
  building blocks returned by \code{case_delete}.
}
\details{
  The primary function is \code{diagnostics} which returns
  either a list or data frame of influence measures
  depending on whether \code{type = "both"} (\code{list})
  or if only one aspect of the model is selected
  (\code{data.frame}). If \code{type = "both"}, then a list
  with Cook's distance, MDFFITS, COVTRACE, and COVRATIO are
  returned for the fixed effects and relative variance
  change (RVC) is returned for the variance components.

  The methods \code{cooks.distance}, \code{mdffits},
  \code{covtrace}, \code{covratio}, and \code{rvc} can be
  used for direct computation of the corresponding
  diagnostic quantities from an object of class
  \code{case_delete}.
}
\note{
  The results provided by this function will give exact
  values of the diagnostics; however, these are
  computationally very slow. Approximate versions of
  \code{cooks.distance}, \code{mdffits}, \code{covtrace},
  \code{covratio} are implemented in HLMdiag, and can be
  called directly on the \code{mer} object.
}
\examples{
data(sleepstudy, package = 'lme4')
fm <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

# Subject level deletion and diagnostics
subject.del  <- case_delete(model = fm, group = "Subject", type = "both")
subject.diag <- diagnostics(subject.del)
}
\author{
  Adam Loy \email{aloy@iastate.edu}
}
\references{
  Christensen, R., Pearson, L.M., and Johnson, W. (1992)
  ``Case-Deletion Diagnostics for Mixed Models,
  \emph{Technometrics}, \bold{34}, 38 -- 45.

  Dillane, D. (2005), Deletion Diagnostics for the Linear
  Mixed Model,'' Ph.D. thesis, Trinity College Dublin.

  Schabenberger, O. (2004) Mixed Model Influence
  Diagnostics, in \emph{Proceedings of the Twenty-Ninth SAS
  Users Group International Conference}, SAS Users Group
  International.
}
\seealso{
  \code{\link{leverage.mer}},
  \code{\link{cooks.distance.mer}},
  \code{\link{mdffits.mer}}, \code{\link{covratio.mer}},
  \code{\link{covtrace.mer}}
}
\keyword{models}
\keyword{regression}

