% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotMap.R
\name{compareMaps}
\alias{compareMaps}
\title{Plot multiple maps
\code{compareMaps}}
\usage{
compareMaps(r, mfrow, contour = FALSE, ...)
}
\arguments{
\item{r}{List of raster or SpatialData}

\item{mfrow}{mfrow in par()}

\item{contour}{Whether plot the contour.}

\item{...}{More options in par()}
}
\description{
Plot multiple maps
\code{compareMaps}
}
\examples{
library(raster)
data(volcano)
r <- raster(volcano)
extent(r) <- c(0, 610, 0, 870)
r1= sin(r/100)
r2= cos(r/100)
compareMaps(list(r,r1), mfrow=c(1,2))
compareMaps(list(r,r1,r2, r1+r2), mfrow=c(2,2))
compareMaps(list(r,r1,r2, r1+r2), mfrow=c(2,2), contour = TRUE)
}
