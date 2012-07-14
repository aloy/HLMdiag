\name{HLMresid}
\alias{HLMresid}
\title{Extracting residuals from two-level HLMs}
\usage{
  HLMresid(object, level = c(1, 2), type = c("both", "LS", "EB", "marginal"), sim = NULL, semi.standardize = TRUE)
}
\arguments{
  \item{object}{a fitted model object of class \code{mer}.}

  \item{level}{the level from which the residuals should be extracted; either 
\code{1} or \code{2} (both can be specified).}

  \item{type}{the type of residual to be extracted; \code{"LS"}, \code{"EB"}, \code{"both"} (i.e. LS and EB), or \code{"marginal"}.}

  \item{sim}{an optional data vector containing a (simulated) response variable.}

  \item{semi.standardize}{if \code{TRUE} the
  semi-standardized residuals will also be returned. This argument is only used when \code{level = 1}.}
}
\value{
  If \code{type = "marginal"} is specified, a numeric vector of the marginal residuals is returned. 
  Otherwise, a data frame is returned containing the following:
  \itemize{
    \item When \code{level=1}: the model frame, the residuals (\code{LS.resid} and/or \code{EB.resid}),  the fitted values (\code{fitted}), and, if \code{semi.standardize = TRUE}, the diagonal elements of the hat matrix (\code{hat}), the semi-standardized residuals (\code{semi.std.resid})

    \item When \code{level=2}: the LS and/or EB residuals corresponding to each random effect.

    \item When \code{level = c(1, 2)}: a list with elements \code{level.1} and \code{level.2}. Each element is specified as stated above.
  }
}
\description{
  This is a wrapper function to extract the residuals from a two-level normal hierarchical
  linear model (HLM) fit using \code{lmer}. The EB residuals are extracted using
  using the functions \code{resid} and \code{ranef} for levels 1 and 2, respectively.
  To extract the LS residuals, \code{LSresids} is used.
}
\examples{
data(Oxboys, package = 'mlmRev')
fm <- lmer(formula = height ~ age + I(age^2) + (age + I(age^2)| Subject), data = Oxboys)
level1Resids <- HLMresid(object = fm, level = c(1,2), type = "both", semi.standardize = TRUE)
}
\author{
  Adam Loy \email{aloy@istate.edu}
}