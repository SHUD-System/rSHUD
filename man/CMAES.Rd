% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Calibration_CMAES.R
\name{CMAES}
\alias{CMAES}
\title{Automatically calibrate model, that requires user-defined model and objective function.
\code{CMAES}}
\usage{
CMAES(
  CV,
  cmd,
  objfunc,
  Call_Model,
  NumPop = CV$method$NUMPOP,
  maxstep = CV$method$MAXGEN,
  sigma = CV$method$SIGMA,
  stopfitness = CV$method$STOPFITNESS,
  debug = FALSE,
  ...
)
}
\arguments{
\item{CV}{Input data, list of data.}

\item{cmd}{Command to run the model.}

\item{objfunc}{User-defined objective function which return the objective values.}

\item{Call_Model}{Function that calls model simulation.}

\item{NumPop}{Number of children in each generation}

\item{maxstep}{Maximum generations}

\item{sigma}{Sigma Value to sample (0, 1)}

\item{stopfitness}{The optimal value. When the objective value is smaller than stopfitness, calibration success.}

\item{debug}{Whether debug Model.}

\item{...}{More options passing to objfunc.}
}
\description{
Automatically calibrate model, that requires user-defined model and objective function.
\code{CMAES}
}
