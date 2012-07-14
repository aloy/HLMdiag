\name{LSresids}
\alias{LSresids}
\title{Calculating OLS residuals for a two-level hierarchical model}
\usage{
  LSresids(object, level, sim = NULL, semi.standardize = TRUE)
}
\arguments{
  \item{object}{a fitted model object of class \code{mer}. }

  \item{level}{the level from which the residuals should be extracted; either 
\code{1} or \code{2} (both can be specified).}

  \item{sim}{an optional data vector containing a (simulated) response variable.}

  \item{semi.standardize}{if \code{TRUE} the
  semi-standardized residuals will also be returned. This argument is only used when \code{level = 1}.}
}
\value{
  a data frame containing the following:
  \itemize{
    \item When \code{level=1}: the model frame, the LS residuals (\code{LS.resid}),  the fitted values (\code{fitted}), and, if \code{semi.standardize = TRUE}, the diagonal elements of the hat matrix (\code{hat}) and the semi-standardized residuals (\code{semi.std.resid})

    \item When \code{level=2}: the LS residuals corresponding to each random effect.
  }
}
\description{
  This function calculates the least squares
  residuals for a two-level hierarchical model. The
  level-1 residuals are found by fitting separate LS regression models
  to each case. The level-2 residuals are found by treating the random 
  effects as fixed and fitting both separate and global linear models. The 
  difference between the coefficients corresponding to the random effects 
  are the level-2 residuals.
}
\examples{
data(Oxboys, package = 'mlmRev')
fm <- lmer(formula = height ~ age + I(age^2) + (age + I(age^2)| Subject), data = Oxboys)
level1Resids <- LSresids(object = fm, level = c(1,2), semi.standardize = TRUE)

\dontrun{wages.fm1 <- lmer(lnw ~ exper + (exper | id), data = wages)
LSresids(fobject = wages.fm1, level = c(1,2), semi.standardize = TRUE)}
}
\author{
  Adam Loy \email{aloy@istate.edu}
}

