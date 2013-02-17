\name{dotplot_diag}
\alias{dotplot_diag}
\title{Dot plots for influence diagnostics}
\usage{
  dotplot_diag(x, index, data, cutoff,
    name = c("cooks.distance", "mdffits", "covratio", "covtrace", "rvc", "leverage"),
    modify = FALSE, ...)
}
\arguments{
  \item{x}{values of the diagnostic of interest}

  \item{index}{index (IDs) of \code{x} values}

  \item{data}{data frame to use (optional)}

  \item{cutoff}{value(s) specifying the boundary for
  unusual values of the diagnostic. The cutoff(s) can
  either be supplied by the user, or automatically
  calculated using measures of internal scaling if
  \code{cutoff = "internal"}}

  \item{name}{what diagnostic is being plotted (one of
  \code{"cooks.distance"}, \code{"mdffits"},
  \code{"covratio"}, \code{"covtrace"}, \code{"rvc"}, or
  \code{"leverage"}). this is used to for the calculation
  of "internal" cutoffs}

  \item{modify}{specifies the \code{geom} to be used to
  produce a space-saving modification: either
  \code{"dotplot"} or \code{"boxplot"}}

  \item{...}{other arguments to be passed to
  \code{qplot()}}
}
\description{
  This is a function that can be used to create (modified)
  dotplots for the diagnostic measures.  The plot allows
  the user to understand the distribution of the diagnostic
  measure and visually identify unusual cases.
}
\note{
  The resulting plot uses \code{coord_flip} to rotate the
  plot, so when adding customized axis labels you will need
  to flip the usage of \code{xlab} and \code{ylab}.
}
\examples{
data(sleepstudy, package = 'lme4')
fm <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

# Subject level deletion and diagnostics
subject.del  <- case_delete(model = fm, group = "Subject", type = "both")
subject.diag <- diagnostics(subject.del)

dotplot_diag(x = COOKSD, index = IDS, data = subject.diag[["fixef_diag"]],
             name = "cooks.distance", modify = FALSE,
             xlab = "Subject", ylab = "Cook's Distance")

dotplot_diag(x = sigma2, , index = IDS, data = subject.diag[["varcomp_diag"]],
             name = "rvc", modify = "dotplot", cutoff = "internal",
             xlab = "Subject", ylab = "Relative Variance Change")

dotplot_diag(x = sigma2, , index = IDS, data = subject.diag[["varcomp_diag"]],
             name = "rvc", modify = "boxplot", cutoff = "internal",
             xlab = "Subject", ylab = "Relative Variance Change")
}
\author{
  Adam Loy \email{aloy@istate.edu}
}
\keyword{hplot}

