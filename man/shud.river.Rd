% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/River.R
\name{shud.river}
\alias{shud.river}
\title{calculate river order, downstream, slope and length
\code{shud.river}}
\usage{
shud.river(sl, dem, rivord = NULL, rivdown = NULL, AREA = NULL)
}
\arguments{
\item{sl}{SpatialLines*}

\item{dem}{Raster of elevation}

\item{rivord}{Order of each river reach}

\item{rivdown}{Downstream Index of each river reach.}

\item{AREA}{Area of the watershed, for estimating the width/depth of river.}
}
\value{
SHUD.RIVER
}
\description{
calculate river order, downstream, slope and length
\code{shud.river}
}
