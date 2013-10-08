\name{HLMresid.default}
\alias{HLMresid}
\alias{HLMresid.default}
\alias{HLMresid.lmerMod}
\alias{HLMresid.mer}
\title{Calculating residuals from HLMs}
\usage{
  \method{HLMresid}{default} (object, ...)

  \method{HLMresid}{mer} (object, level, type = "EB",
    sim = NULL, standardize = FALSE, ...)

  \method{HLMresid}{lmerMod} (object, level, type = "EB",
    sim = NULL, standardize = FALSE, ...)
}
\arguments{
  \item{object}{an object of class \code{mer} or
  \code{lmerMod}.}

  \item{level}{which residuals should be extracted: 1 for
  within-group (case-level) residuals, the name of a
  grouping factor (as defined in \code{flist} of the
  \code{mer} object) for between-group residuals, or
  \code{marginal}.}

  \item{type}{how are the residuals predicted: either
  \code{"EB"} or \code{"LS"} (the default is \code{"EB"}).}

  \item{sim}{optional argument giving the data frame used
  for LS residuals. This is used mainly for dealing with
  simulations.}

  \item{standardize}{if \code{standardize = TRUE} the
  standardized residuals will be returned; if
  \code{standardize = "semi"} then the semi-standardized
  level-1 residuals will be returned. Note that for
  higher-level residuals of \code{type = "LS"},
  \code{standardize = TRUE} does not result in standardized
  residuals as they have not been implemented.}

  \item{...}{do not use}
}
\description{
  \code{HLMresid} is a function that extracts residuals
  from a hierarchical linear model fit using \code{lmer}.
  That is, it is a unified framework that
  extracts/calculates residuals from \code{mer} or
  \code{lmerMod} objects.
}
\details{
  This function extracts residuals from the model, and can
  extract residuals estimated using least squares (LS) or
  Empirical Bayes (EB). This unified framework enables the
  analyst to more easily conduct an upward residual
  analysis during model exploration/checking.

  The \code{HLMresid} function provides a wrapper that will
  extract residuals from a fitted \code{mer} or
  \code{lmerMod} object. The function provides access to
  residual quantities already made available by the
  functions \code{resid} and \code{ranef}, but adds
  additional functionality. Below is a list of types of
  residuals that can be extracted. \describe{ \item{raw
  level-1 residuals}{These are equivalent to the residuals
  extracted by \code{resid} if \code{level = 1}, \code{type
  = "EB"}, and \code{standardize = FALSE} is specified. You
  can also specify \code{type = "LS"} for LS residuals that
  are not equivalent to those from \code{resid}.}
  \item{standardized level-1 residuals}{Specify \code{level
  = 1}, and \code{standardize = TRUE}. This works with both
  \code{type = "EB"} or \code{"LS"}.}
  \item{semi-standardized level-1 residuals}{Specify
  \code{level = 1}, \code{type = "LS"} and
  \code{standardize = "semi"}.} \item{raw group level
  residuals}{These are equivalent to extracting the
  predicted random effects for a given group using
  \code{ranef}. Set \code{level} to a grouping factor name
  and \code{type = "EB"}. \code{type = "LS"} can also be
  specified, though this is less common.}
  \item{standardized group level residuals}{Set
  \code{level} to a grouping factor name, \code{type =
  "EB"}, and \code{standardized = TRUE}. This will not
  produce standardized residuals for \code{type = "LS"}.}
  \item{marginal residuals}{The marginal residuals can be
  obtained by setting \code{level = "marginal"}. Only
  \code{type = "EB"} is implemented.} \item{cholesky
  residuals}{These are essentially standardized marginal
  residuals. To obtain cholesky residuals set \code{level =
  "marginal"}, \code{type = "EB"}, and \code{standardize =
  TRUE}.} } Note that \code{standardize = "semi"} is only
  implemented for level-1 LS residuals.
}
\examples{
data(sleepstudy, package = "lme4")
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)

# level-1 residuals
all.equal(HLMresid(object = fm1, level = 1, type = "EB"), resid(fm1)) ## EB
r1LS <- HLMresid(object = fm1, level = 1, type = "LS") ## raw LS resids
head(r1LS)
r1LS.std <- HLMresid(object = fm1, level = 1, type = "LS", standardize = TRUE) ## std. LS resids
head(r1LS.std)

# level-2 residuals
all.equal(r2EB <- HLMresid(object = fm1, level = "Subject", type = "EB"), ranef(fm1)[["Subject"]])
r2EB.std <- HLMresid(object = fm1, level = "Subject", type = "EB", standardize = TRUE)
head(r2EB)
head(r2EB.std)

# marginal residuals
mr <- HLMresid(object = fm1, level = "marginal")
cholr <- HLMresid(object = fm1, level = "marginal", standardize = TRUE) # Cholesky residuals
}
\author{
  Adam Loy \email{loyad01@gmail.com}
}
\references{
  Hilden-Minton, J. (1995) Multilevel diagnostics for mixed
  and hierarchical linear models. University of California
  Los Angeles.

  Houseman, E. A., Ryan, L. M., & Coull, B. A. (2004)
  Cholesky Residuals for Assessing Normal Errors in a
  Linear Model With Correlated Outcomes. \emph{Journal of
  the American Statistical Association}, \bold{99}(466),
  383--394.
}
\seealso{
  \code{\link{LSresids}}, \code{\link{resid}},
  \code{\link{ranef}}
}
\keyword{models}
\keyword{regression}

