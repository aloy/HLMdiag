\name{ggplot_ppnorm}
\alias{ggplot_ppnorm}
\title{Constructing a normal probability plot}
\usage{
  ggplot_ppnorm(x, line = NULL, ...)
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
  This function will construct a normal probability plot
  within the \code{ggplot} framework.
}
\author{
  Adam Loy \email{aloy@istate.edu}
}

