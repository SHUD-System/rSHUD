% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GIS_Projection.R
\name{crs.long2utmZone}
\alias{crs.long2utmZone}
\title{Get UTM zone frome longitutde
\code{crs.long2utmZone}}
\usage{
crs.long2utmZone(lon)
}
\arguments{
\item{lon}{Longitude in degree}
}
\value{
Zone number
}
\description{
Get UTM zone frome longitutde
\code{crs.long2utmZone}
}
\examples{
long = seq(-180, 179, 6)+1
crs.long2utmZone(long)
}
