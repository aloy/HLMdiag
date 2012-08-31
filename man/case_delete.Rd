\name{case_delete}
\alias{case_delete}
\title{Implementing Case Deletion.}
\usage{
  case_delete(model, group = FALSE, type = c("both",
  "fixef", "varcomp"))
}
\arguments{
  \item{model}{the original hierarchical model fit using
  \code{lmer()}}

  \item{group}{a logical value used whether a group of observations
  or individual observations will be deleted. If this is left
  \code{FALSE}, then the function will delete individual
  observations.}

  \item{type}{the part of the model for which you are
  obtaining deletion diagnostics: the fixed effects
  (\code{fixef}), variance components (\code{varcomp}), or
  \code{both}}
}
\value{
  a list with the following compontents:
  \item{fixef.original}{the original fixed effects}
  \item{ranef.original}{the origingal random effects}
  \item{vcov.original}{the original variance-covariance
  parameters} \item{varcomp.original}{the original variance
  components} \item{fixef.delete}{a list of the fixed
  effects obtained through case deletion}
  \item{ranef.delete}{a list of the random effects obtained
  through case deletion} \item{vcov.delete}{a list of the
  variance-covariance parameters obtained through case
  deletion} \item{fitted.delete}{a list of the fitted
  values obtained through case deletion}
  \item{varcomp.delete}{a list of the variance components
  obtained through case deletion}
}
\description{
  This function is used to iteratively delete groups
  corresponding to the levels of a two-level hierarchical
  linear model. It uses \code{lmer()} to fit the models for
  each deleted case (i.e. uses brute force). To investigate
  numerous levels of the model, the function will need to
  be called multiple times, specifying the group (level) of
  interest each time.
}
\examples{
data(Oxboys, package = 'mlmRev')
fm <- lmer(formula = height ~ age + I(age^2) + (age + I(age^2)| Subject), data = Oxboys)
fmDel <- case_delete(model = fm, group = TRUE, type = "both")

\dontrun{library(mlmRev)
exm1 <- lmer(normexam ~ standLRT + sex + schgend + (1 | school), data = Exam)
exm1DEL <- case_delete(model = exm1, group = TRUE, type = "both")}
}
\author{
  Adam Loy \email{aloy@istate.edu}
}
\references{
  Christensen, R., Pearson, L.M., and Johnson, W. (1992),
  ``Case-Deletion Diagnostics for Mixed Models,''
  \emph{Technometrics}, 34, 38 -- 45.

  Schabenberger, O. (2004),``Mixed Model Influence
  Diagnostics,'' in \emph{Proceedings of the Twenty-Ninth
  SAS Users Group International Conference}, SAS Users
  Group International.
}

