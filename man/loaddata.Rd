% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readout.R
\name{loaddata}
\alias{loaddata}
\title{Read multiply SHUD model output files}
\usage{
loaddata(
  varname = c(paste0("eley", c("surf", "unsat", "gw", "snow")), paste0("elev",
    c("prcp", "infil", "rech")), paste0("elev", c("etp", "eta", "etev", "ettr", "etic")),
    paste0("rivq", c("down", "sub", "surf")), paste0("rivy", "stage")),
  rdsfile = NULL
)
}
\arguments{
\item{varname}{vector of output keywords}

\item{rdsfile}{Save RDS file. NULL=do not save rds file.}
}
\value{
A list of TimeSeries data.
}
\description{
Read multiply SHUD model output files
}
\keyword{output.}
\keyword{read}
