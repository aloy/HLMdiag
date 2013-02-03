\name{dotplot_diag}
\alias{dotplot_diag}
\title{Dot plots for influence diagnostics}
\usage{
  dotplot_diag(data, type = c("fixef", "varcomp"), name,
    cutoff = NULL, modify = FALSE, ...)
}
\arguments{
  \item{data}{an object containing the output from
  \code{diagnostics()}.}

  \item{type}{the part of the model the diagnostic
  corresponds to, either \code{"fixef"} or
  \code{"varcomp"}.}

  \item{name}{specification of which diagnostic to plot
  (either COOKSD, MDFFITS, COVTRACE, COVRATIO, or rvc).}

  \item{cutoff}{value specifying unusual values of the
  diagnostic}

  \item{modify}{if \code{TRUE} will produce a space-saving
  modification}

  \item{...}{other arguments to be passed to
  \code{qplot()}}
}
\description{
  This is a function that can be used to create modified
  dotplots for the diagnostic measures.  The plot allows
  the user to understand the distribution of the diagnostic
  measure and visually identify unusual cases.
}
\note{
  The resulting plot is uses \code{coord_flip} to rotate
  the plot, so when adding customized axis labels you will
  need to flip the usage of \code{xlab} and \code{ylab}.
}
\examples{
data(sleepstudy, package = 'lme4')
fm <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

# Subject level deletion and diagnostics
subject.del  <- case_delete(model = fm, group = "Subject", type = "fixef")
subject.diag <- diagnostics(subject.del)

dotplot_diag(data = subject.diag, type = "fixef", name = "COOKSD",
             cutoff = "internal", modify = TRUE, xlab = "Subject",
             ylab = "Cook's Distance")
dotplot_diag(data = subject.diag, type = "fixef", name = "COOKSD",
             cutoff = "internal", modify = FALSE, xlab = "Subject",
             ylab = "Cook's Distance")
}
\author{
  Adam Loy \email{aloy@istate.edu}
}
\keyword{hplot}

