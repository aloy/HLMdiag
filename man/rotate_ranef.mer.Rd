\name{rotate_ranef.default}
\alias{rotate_ranef}
\alias{rotate_ranef.default}
\alias{rotate_ranef.lmerMod}
\alias{rotate_ranef.mer}
\title{Calculate s-dimensional rotated random effects}
\usage{
  \method{rotate_ranef}{default} (.mod, ...)

  \method{rotate_ranef}{mer} (.mod, .L, s = NULL,
    .varimax = FALSE, ...)

  \method{rotate_ranef}{lmerMod} (.mod, .L, s = NULL,
    .varimax = FALSE, ...)
}
\arguments{
  \item{.mod}{an object of class \code{mer} or
  \code{lmerMod}.}

  \item{.L}{a matrix defining which combination of random
  effects are of interest.}

  \item{s}{the dimension of the subspace of interest.}

  \item{.varimax}{if \code{.varimax = TRUE} than the raw
  varimax rotation will be applied to the resulting
  rotation.}

  \item{...}{do not use}
}
\description{
  This function calculates reduced dimensional rotated
  random effects. The rotation reduces the influence of the
  residuals from other levels of the model so that
  distributional assessment of the resulting random effects
  is possible.
}
\author{
  Adam Loy \email{loyad01@gmail.com}
}

