\name{diagnostics}
\alias{cooksd_hlm}
\alias{covratio_hlm}
\alias{covtrace_hlm}
\alias{diagnostics}
\alias{mdffits_hlm}
\alias{rvc}
\title{Calculating diagnostics for two-level hierarchical linear models.}
\usage{
  diagnostics(model, delete)

  cooksd_hlm(model, delete)

  mdffits_hlm(model, delete)

  covtrace_hlm(model, delete)

  covratio_hlm(model, delete)

  rvc(delete)
}
\arguments{
  \item{model}{an object contatining the original
  hierarchical model fit using \code{lmer()}}

  \item{delete}{an object containing the output returned by
  \code{case_delete()}}
}
\description{
  This group of functions is used to compute deletion
  diagnostics for a two-level normal hierarchical model at
  both levels of the model.
}
\details{
  The primary function is \code{diagnostics} which returns
  either a list or data frame of influence measures
  depending on whether \code{type = "both"} or if only one
  aspect of the model is selected. If \code{type = "both"},
  then a list with Cook's distance, MDFFITS, COVTRACE, and
  COVRATIO are returned for the fixed effects and relative
  variance change (RVC) is returned for the variance
  components.

  The functions \code{cooksd_hlm}, \code{mdffits_hlm},
  \code{covtrace_hlm}, \code{covratio_hlm}, and \code{rvc}
  can be used for direct computation of the corresponding
  diagnostic quantities.
}
\examples{
data(Oxboys, package = "mlmRev")
fm <- lmer(formula = height ~ age + I(age^2) + (age + I(age^2)| Subject), data = Oxboys)
fmDel <- case_delete(model = fm, group = TRUE, type = "both")
fmDiag <- diagnostics(model = fm, delete = fmDel)

\dontrun{
library(mlmRev)
exm1 <- lmer(normexam ~ standLRT + sex + schgend + (1 | school), data = Exam)
exm1DEL <- case_delete(model = exm1, group = TRUE, type = "both")
exm1DIAG <- diagnostics(model = exm1, delete = exm1DEL)}
}
\author{
  Adam Loy \email{aloy@iastate.edu}
}
\references{
  Christensen, R., Pearson, L.M., and Johnson, W. (1992),
  ``Case-Deletion Diagnostics for Mixed Models,''
  \emph{Technometrics}, 34, 38 -- 45.

  Dillane, D. (2005), ``Deletion Diagnostics for the Linear
  Mixed Model,'' Ph.D. thesis, Trinity College Dublin.

  Schabenberger, O. (2004),``Mixed Model Influence
  Diagnostics,'' in \emph{Proceedings of the Twenty-Ninth
  SAS Users Group International Conference}, SAS Users
  Group International.
}

