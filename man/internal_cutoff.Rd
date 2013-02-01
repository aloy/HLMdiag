\name{internal_cutoff}
\alias{internal_cutoff}
\title{Calculating a cutoff value for diagnostic measures}
\usage{
  internal_cutoff(data, type, name)
}
\arguments{
  \item{data}{an object containing the output from
  \code{diagnostics()}}

  \item{type}{the part of the model the diagnostic
  corresponds to, either \code{"fixef"} or
  \code{"varcomp"}.}

  \item{name}{specification of which diagnostic to plot
  (either COOKSD, MDFFITS, COVTRACE, COVRATIO, or rvc).}
}
\description{
  This function provides cutoff values using internal
  scaling. In other words, a measure of relative standing
  is used (3*IQR) to specify unusual values of the
  diagnostic of interest relative to the vector of
  diagnostics.
}
\author{
  Adam Loy \email{aloy@istate.edu}
}

