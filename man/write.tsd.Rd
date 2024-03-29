% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/writeInput.R
\name{write.tsd}
\alias{write.tsd}
\title{Write xts data out into file
\code{write.tsd}}
\usage{
write.tsd(x, file, dt_units = "days", append = F, quite = F, header = NULL)
}
\arguments{
\item{x}{xts data}

\item{file}{file name}

\item{dt_units}{Time interval of the 1 unit change in time column.}

\item{append}{whether append}

\item{quite}{TRUE/FALSE, if quiet mode}

\item{header}{Header of the file. Default is the dimension of data.frame.}
}
\description{
Write xts data out into file
\code{write.tsd}
}
