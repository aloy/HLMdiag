\name{leverage.mer}
\alias{leverage}
\alias{leverage.mer}
\title{Leverage for mixed/hierarchical linear models}
\usage{
  \method{leverage}{mer} (model, level = "both")
}
\arguments{
  \item{model}{fitted model of class \code{mer}}

  \item{the}{level at which the leverage should be
  calculated; either 1 or 2 (\code{"both"} can be
  specified)}
}
\description{
  This function calculates the leverage of a
  mixed/hierarchical model fit by \code{lmer}.
}
\author{
  Adam Loy \email{aloy@iastate.edu}
}
\references{
  Nobre, J. S., & Singer, J. M. (2011).  Leverage analysis
  for linear mixed models.  Journal of Applied Statistics,
  38(5), 1063–1072.

  Demidenko, E., & Stukel, T. A. (2005).  Influence
  analysis for linear mixed-effects models.  Statistics in
  Medicine, 24(6), 893–909.
}
\seealso{
  \code{\link{cooks.distance.mer}},
  \code{\link{mdffits.mer}}, \code{\link{covratio.mer}},
  \code{\link{covtrace.mer}}, \code{\link{rvc.mer}}
}

