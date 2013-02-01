\name{HLMresid}
\alias{HLMresid}
\title{Calculating residuals from two-level HLMs}
\usage{
  HLMresid(object, level, type = "EB", sim = NULL,
    standardize = TRUE)
}
\arguments{
  \item{object}{an object of class \code{mer}.}

  \item{level}{which residuals should be plotted: 1 for
  within-group residuals or the name of the grouping factor
  (as defined in \code{flist} of the \code{mer} object) for
  between-group residuals, or \code{marginal}.}

  \item{type}{how are the residuals predicted: either
  \code{EB} or \code{LS}.}

  \item{sim}{data optional argument giving the data frame
  used for LS residuals. This is used mainly for when
  dealing with simulations.}

  \item{standardize}{if \code{TRUE} the standardized
  level-1 residuals will also be returned (if \code{level =
  1}); if \code{"semi"} then the semi-standardized level-1
  residuals will be returned.}
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
\author{
  Adam Loy \email{aloy@istate.edu}
}
\keyword{models}
\keyword{regression}

