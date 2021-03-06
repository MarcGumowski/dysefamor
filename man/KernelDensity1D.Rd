% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DSFM.R
\name{KernelDensity1D}
\alias{KernelDensity1D}
\title{One-Dimensional Kernel Estimation}
\usage{
KernelDensity1D(y, I, J, x1, u, U, h)
}
\arguments{
\item{y}{a numeric matrix of the data.}

\item{I}{the total number of rows in \code{y}.}

\item{J}{the total number of columns in \code{y}.}

\item{x1}{a numeric vector of the covariates.}

\item{u}{a numeric vector of the estimation grid.}

\item{U}{the length of \code{u}.}

\item{h}{a numeric vector of the bandwidth.}
}
\value{
It returns \eqn{\hat{p}_t(u)} and \eqn{\hat{q}_t(u)}.
}
\description{
This function computes the one-dimensional kernel estimation needed for the
\code{\link{DSFM1D}} funtion to work. It is using the
\code{\link{QuarticKernel1D}} function.
}
\references{
Fengler, Matthias R, Wolfgang K Haerdle, and Enno Mammen (2007).
"A Semiparametric Factor Model for Implied Volatility Surface Dynamics".
In: \emph{Journal of Financial Econometrics 5.2}, pp. 189-218.
}
\seealso{
Other kernel.functions: \code{\link{KernelDensity2D}},
  \code{\link{NormalKernel1D}},
  \code{\link{QuarticKernel1D}}
}

