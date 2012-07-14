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
  (either COOKSD, MDFFITS, COVTRACE, COVRATIO, or the name
  of the variance component for rvc).}

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
\examples{
data(Oxboys, package = 'mlmRev')
fm <- lmer(formula = height ~ age + I(age^2) + (age + I(age^2)| Subject), data = Oxboys)
fmDel <- case_delete(model = fm, group = "Subject", type = "both")
fmDiag <- diagnostics(model = fm, delete = fmDel)
dotplot_diag(data = fmDiag, type = "fixef", name = "COOKSD", cutoff = "internal",
 	modify = TRUE, xlab = "Subject", ylab = "Cook's Distance")
dotplot_diag(data = fmDiag, type = "fixef", name = "COOKSD", cutoff = "internal",
 	modify = FALSE, xlab = "Subject", ylab = "Cook's Distance")
}
\author{
  Adam Loy \email{aloy@istate.edu}
}

