\name{case_delete2}
\alias{case_delete2}
\title{Implementing Case Deletion.}
\usage{
case_delete2(model, group = FALSE, type = c("both", "fixef", "varcomp"))
}
\arguments{
  \item{model}{the original hierarchical model fit using
  \code{lmer()}}

  \item{group}{a variable used to define the group for
  which cases will be deleted.  If this is left
  \code{FALSE}, then the function will delete individual
  observations.}

  \item{type}{the part of the model for which you are
  obtaining deletion diagnostics: the fixed effects
  (\code{fixef}), random effects (\code{ranef}), variance
  components (\code{varcomp}), or \code{all}.}
}
\value{
  a list with the following compontents:

  fixef.original the original fixed effects

  ranef.original the origingal random effects

  vcov.original the original variance-covariance parameters

  varcomp.original the original variance components

  fixef.delete a list of the fixed effects obtained through
  case deletion

  ranef.delete a list of the changes in random effects
  obtained through case deletion

  vcov.delete a list of the variance-covariance parameters
  obtained through case deletion -- currently NULL

  fitted.delete a list of the fitted values obtained
  through case deletion -- Currently NULL

  varcomp.delete a list of the variance components obtained
  through case deletion -- Currently NULL
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
\dontrun{
data(Exam, package = "mlmRev")
fm <- lmer(normexam ~ standLRT + I(standLRT^2) + I(standLRT^3) + schgend + (schgend | school), data = Exam)
fm.del <- case_delete2(fm)
}
}
