\name{group_qqnorm}
\alias{group_qqnorm}
\title{Overlaying normal Q-Q plots}
\usage{
  group_qqnorm(x, group, line = NULL, alpha_point = 1,
    alpha_line = 1, ...)
}
\arguments{
  \item{x}{a numeric vector from which quantiles will be
  calculated}

  \item{group}{a vector indicating group membership for
  each value in \code{x}.}

  \item{line}{the method used to fit reference lines. If no
  reference lines are desired, leave the value as
  \code{NULL}. \code{line = "rlm"} will use robust
  regression to fit reference lines. \code{line =
  "quantile"} will fit lines through the first and third
  quartiles.}

  \item{...}{other arguments to be passed to \code{ggplot}}

  \item{alpha_point}{alpha value specified for the points}

  \item{alpha_line}{alpha value specified for the lines}
}
\description{
  This function will overlay multiple normal Q-Q plots on
  the same plot. This will be particulary useful when
  comparing the distribution between groups. In this
  situation, significantly different slopes would indicate
  the normal distributions for the groups do not share a
  common standard deviation.
}
\author{
  Adam Loy \email{aloy@istate.edu}
}
\references{
  Hilden-Minton, J. A. (1995) Mulilevel Diagnostics for
  Mixed and Hierarchical Linear Models, Ph.D. thesis,
  University of California Los Angeles.
}
\keyword{hplot}

