\name{ggplot_qqnorm}
\alias{ggplot_qqnorm}
\title{Constructing a normal Q-Q plot}
\usage{
  ggplot_qqnorm(x, line = NULL, ...)
}
\arguments{
  \item{x}{a numeric vector}

  \item{line}{the method used to fit a reference line. If
  no reference line is desired, leave the value as
  \code{NULL}. \code{line = "rlm"} will use robust
  regression to fit a reference line. \code{line =
  "quantile"} will fit a line through the first and third
  quartiles. These options are the same as those given for
  the \code{qqPlot} function in the \code{car} package.}

  \item{...}{other arguments to be passed to
  \code{qplot()}}
}
\description{
  This function will construct a normal Q-Q plot within the
  \code{ggplot2} framework. It combines the functionality
  of \code{qqnorm} and \code{qqline}.
}
\author{
  Adam Loy \email{loyad01@gmail.com}
}
\seealso{
  \code{\link{qqnorm}}, \code{\link{qqline}}
}
\keyword{hplot}

