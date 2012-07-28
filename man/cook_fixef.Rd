\name{cook_fixef}
\alias{cook_fixef}
\title{Cook's distance for fixed effects}
\usage{
cook_fixef(bb, fixef, sig0)
}
\arguments{
  \item{bb}{list of building blocks for level-1 deletion}

  \item{fixef}{fixed effects from mer object}

  \item{sig0}{}
}
\description{
  Cook's distance for fixed effects
}

\examples{
data(Exam, package = 'mlmRev')
fm <- lmer(normexam ~ standLRT + I(standLRT^2) + I(standLRT^3) + schgend + (schgend | school), data = Exam)
fm.del <- case_delete2(fm)
fm.cd  <- cook_fixef(fm.del$fixef.delete, fixef(fm), sqrt(fm.del$varcomp.original[1]))}
}