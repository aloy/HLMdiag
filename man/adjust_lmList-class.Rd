\docType{class}
\name{adjust_lmList-class}
\alias{adjust_lmList-class}
\alias{coef,adjust_lmList-method}
\alias{confint,adjust_lmList-method}
\alias{formula,adjust_lmList-method}
\alias{plot,adjust_lmList.confint,ANY-method}
\alias{show,adjust_lmList-method}
\title{Class \code{"adjust_lmList"}}
\description{
  Adjusted List of \code{lm} Objects with a Common Model
}
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("adjust_lmList", ...)} or by
  \code{adjust_lmList}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"}, a list of \code{lm} objects}
    \item{\code{call}:}{Object of class \code{"call"}, the \code{call} used to create the \code{adjust_lmList} object.}
    \item{\code{pool}:}{Object of class \code{"logical"}, a logical expression stating whether the pooled standard deviation should be estimate should be used.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
}
\section{Methods}{
  \describe{
  \item{formula}{\code{signature(x = "adjust_lmList")}: extracts the formula
      	from an adjust_lmList object.}
    \item{coef}{\code{signature(object = "adjust_lmList")}: extracts the model 
		coefficients from the individually fit models.}
    \item{confint}{\code{signature(object = "adjust_lmList")}: computes confidence 
		intervals for the for the parameters of the individually fit models. }
    \item{show}{\code{signature(object = "adjust_lmList")}: extracts the formula 
		used to create the \code{adjust_lmList} object.}
	\item{plot}{\code{signature(object = "adjust_lmList.confint")}: plot the confidence
	 intervals for the parameters of the individually fit models using the \code{ggplot2}
	 framework.} 
	 }
}
\keyword{classes}