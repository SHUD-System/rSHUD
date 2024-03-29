% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WaterBalance.R
\name{wb.lake}
\alias{wb.lake}
\title{Calculate the Change of Storage.
\code{wb.lake}}
\usage{
wb.lake(
  xl = loaddata(varname = paste0("lak", c("ystage", "atop", "qrivin", "qrivout",
    "qsub", "qsurf", "vevap", "vprcp"))),
  lakeid = NULL
)
}
\arguments{
\item{xl}{List of data.}

\item{lakeid}{Index of the lake for waterbalance calculation.}
}
\value{
A list of water balance by each lake; a matrix in each list-item
}
\description{
Calculate the Change of Storage.
\code{wb.lake}
}
