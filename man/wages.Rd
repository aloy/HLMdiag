\name{wages}
\docType{data}
\alias{wages}
\title{Wages for male high school dropouts}
\description{
	Data on the labor-market experience of male high school dropouts.
}
\usage{wages}
\format{
	A data frame with 6402 observations on the following 15 variables.
	\describe{
	\item{id}{respondent id - a factor with 888 levels.}
	\item{lnw}{natural log of wages expressed in 1990 dollars.}
	\item{exper}{years of experience in the work force}
	\item{ged}{equals 1 if respondent has obtained a GED as of the time of 
	survey, 0 otherwise}
	\item{postexp}{labor force participation since obtaining a GED (in years) - 
	before a GED is earned postexp = 0, and on the day a GED is earned postexp = 0}
	\item{black}{factor - equals 1 if subject is black, 0 otherwise}
	\item{hispanic}{factor - equals 1 if subject is hispanic, 0 otherwise}
	\item{hgc}{highest grade completed - takes integers 6 through 12}
	\item{hgc.9}{hgc - 9, a centered version of hgc}
	\item{uerate}{ocal area unemployment rate for that year}
	\item{ue.7}{}
	\item{ue.centert1}{}
	\item{ue.mean}{}
	\item{ue.person.cen}{}
	\item{ue1}{}
	}
}
\source{
These data are originally from the 1979 National Longitudinal Survey on Youth (NLSY79) that can be found here \url{http://www.bls.gov/nls/nlsdata.htm}.

Singer and Willett (2003) used these data for examples in chapter (insert info. here) and the data sets used can be found on the UCLA Statistical Computing 
website: \url{http://www.ats.ucla.edu/stat/examples/alda/}

Additionally the data were discussed by Cook and Swayne (2003) and the data can be found on the GGobi website: \url{http://www.ggobi.org/book/}.s
}
\references{
Singer, J. D. and Willett, J. B. (2003), \emph{Applied Longitudinal Data Analysis: Modeling Change and Event Occurence}, New York: Oxford University Press.

Cook, D. and Swayne, D. F. (2007), \emph{Interactive and Dynamic Graphics for Data Analysis with R and GGobi}, Springer.
}
\examples{
str(wages)
summary(wages)
lmer(lnw ~ exper + (exper | id), data = wages)
}