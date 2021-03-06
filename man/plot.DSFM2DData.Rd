% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DSFM_2D.R
\name{plot.DSFM2DData}
\alias{plot.DSFM2DData}
\title{Plot Method for an \code{DSFM2DData} Object}
\usage{
\method{plot}{DSFM2DData}(x, n = 1, pal = "pink", theta = 40,
  border = NA, box = T, shade = 0.2, expand = 0.5,
  ticktype = "simple", ...)
}
\arguments{
\item{x}{an \code{object} of class \code{"DSFM2DData"}.}

\item{n}{number of time indicator to be plotted.}

\item{pal}{the color palette for the fit plot. To choose between
\code{"pink"},\code{"blue"},\code{"light"},\code{"dark"}.}

\item{theta}{angles defining the viewing direction of the fit plot.}

\item{border}{the color of the line drawn around the surface facets of the
fit plot. The default, NULL, corresponds to par("fg"). A value of NA will
disable the drawing of borders: this is sometimes useful when the surface is
shaded.}

\item{box}{should the bounding box for the surface be displayed.
The default is \code{FALSE}.}

\item{shade}{the shade at a surface facet is computed as ((1+d)/2)^shade,
where d is the dot product of a unit vector normal to the facet and a unit
vector in the direction of a light source. Values of shade close to one
yield shading similar to a point light source model and values close to zero
produce no shading. Values in the range 0.5 to 0.75 provide an approximation
to daylight illumination.}

\item{expand}{a expansion factor applied to the z coordinates. Often used
with 0 < expand < 1 to shrink the plotting box in the z direction.}

\item{ticktype}{character: "simple" draws just an arrow parallel to the axis
to indicate direction of increase; "detailed" draws normal ticks as per 2D
plots.}

\item{...}{other parameters to be passed through to plotting
\code{\link{persp}} function.}
}
\description{
Plots the 3D vizualisation of the data set contained in \code{object} of class
\code{DSFM2DData}.
}
\seealso{
\code{\link{DSFM2DData}}, \code{\link{summary.DSFM2DData}}.
}

