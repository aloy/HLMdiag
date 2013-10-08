\name{compare_eb_ls}
\alias{compare_eb_ls}
\title{Visually comparing shrinkage and LS estimates}
\usage{
  compare_eb_ls(eb, ols, identify = FALSE, silent = TRUE,
    ...)
}
\arguments{
  \item{eb}{a matrix of random effects}

  \item{ols}{a matrix of the OLS estimates found using
  \code{random_ls_coef}}

  \item{identify}{the percentage of points to identify as
  unusual, \code{FALSE} if you do not want the points
  identified.}

  \item{silent}{logical: should the list of data frames
  used to make the plots be supressed.}

  \item{...}{other arguments to be passed to
  \code{qplot()}}
}
\description{
  This function creates a plot (using \code{qplot()}) where
  the shrinkage estimate appears on the horizontal axis and
  the LS estimate appears on the vertical axis.
}
\examples{
wages.fm1 <- lmer(lnw ~ exper + (exper | id), data = wages)
wages.sepLM <- adjust_lmList(lnw ~ exper | id, data = wages)
rancoef.eb <- coef(wages.fm1)$id
rancoef.ols <- coef(wages.sepLM)
compare_eb_ls(eb = rancoef.eb, ols = rancoef.ols, identify = 0.01)
}
\author{
  Adam Loy \email{loyad01@gmail.com}
}

