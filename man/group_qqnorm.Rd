\name{group_qqnorm}
\alias{group_qqnorm}
\title{Plotting numerous normal q-q plots on same plot}
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

  \item{...}{other arguments to be passed to ggplot}

  \item{alpha_point}{alpha value specified for the points}

  \item{alpha_line}{alpha value specified for the lines}
}
\description{
  This function will plot multiple normal q-q plots on the
  same plot. This will be particulary useful when comparing
  the distrubtion between groups. As differing slopes would
  indicate the normal distributions for the groups do not
  share a common slope.
}
\author{
  Adam Loy \email{aloy@istate.edu}
}
\references{
  Hilden-Minton, J. A. (1995), ``Mulilevel Diagnostics for
  Mixed and Hierarchical Linear Models,'' Ph.D. thesis,
  University of California Los Angeles.
}

