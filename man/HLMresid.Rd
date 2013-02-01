\name{HLMresid}
\alias{HLMresid}
\title{Calculating residuals from two-level HLMs}
\usage{
  HLMresid(object, level = c(1, 2),
    type = c("both", "LS", "EB", "marginal"), sim = NULL,
    semi.standardize = TRUE)
}
\arguments{
  \item{object}{an object of class \code{mer}.}

  \item{level}{specification of the residual of interest.}

  \item{how}{the residuals should be estimated; default is
  \code{"all"}.}

  \item{if}{\code{type = "all"} or \code{"LS"}, a linear
  formula used by \code{adjust_lmList} (y ~ x1 + ... + xn |
  g where g is a grouping factor) must be specified.}

  \item{data}{optional argument giving the data frame used
  for LS residuals. This is used mainly for when dealing
  with simulations.}
}
\description{
  \code{HLMresid} is a function that extracts residuals
  from two-level mixed/hierarchical linear models fit using
  \code{lmer}. That is, it is a function that will extract
  residuals from \code{mer} objects in a unified framework.
}
\details{
  This function can extract residuals from either level of
  the model, and can extract residuals estimated using
  least squares (LS), Empirical Bayes (EB), or both. This
  unified framework enables the analyst to more easily
  conduct an upward residual analysis during model
  exploration or checking.
}
\keyword{models}
\keyword{regression}

