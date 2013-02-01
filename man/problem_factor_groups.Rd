\name{problem_factor_groups}
\alias{problem_factor_groups}
\title{Adjusting formula for lmList}
\usage{
  problem_factor_groups(formula, data)
}
\arguments{
  \item{formula}{a linear formula that is used by
  \code{lmList()} e.g. y ~ x1 + x2 + ... + xn | g}

  \item{data}{the model frame from the model fit by
  \code{lmer()}}
}
\value{
  a list containing the adjusted formulas
}
\description{
  When dealing with nested data, sometimes a factor will
  only take on one level within a group. This causes an
  error, since a factor must have 2 or more levels to be
  included as a predictor within a formula for a linear
  model. When a factor takes on only one level, the effect
  of that factor can then be included in the intercept of
  the model. This function takes a formula for use with
  \code{lmList} and will make the necessary adjustments to
  avoid errors while fitting the separate LS models.
}
\author{
  Adam Loy \email{aloy@iastate.edu}
}
\keyword{models}
\keyword{regression}

