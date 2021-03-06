% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DSFM_1D.R
\name{summary.DSFM1D}
\alias{summary.DSFM1D}
\title{Summarizing DSFM Fits for One-Dimensional Data}
\usage{
\method{summary}{DSFM1D}(object, ...)
}
\arguments{
\item{object}{an object of class \code{"DSFM1D"}.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
The function \code{summary.DSFM1D} returns a list of summary
statistics of the fitted DSFM given in \code{object}, using the components
(list elements) \code{"call"} from its arguments, plus:
\item{\code{ZHat}}{a summary of the estimated factor loadings.}
\item{\code{EV}}{the Explained Variance, used to select the approriate
number of factors.}
\item{\code{RMSE}}{the Root Mean Squared Error, used to compare the goodness-
of-fit between models.}
\item{\code{AIC}}{The bandwidth and its selection criteria.}
}
\description{
\code{summary} method for class \code{"DSFM1D"}.
}
\seealso{
\code{\link{DSFM1DData}}, \code{\link{DSFM}}, \code{\link{DSFM1D}},
\code{\link{plot.DSFM1D}}, \code{\link{predict.DSFM1D}}.
}

