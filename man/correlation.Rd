% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlation.R
\docType{class}
\name{corSAR1-class}
\alias{corSAR1-class}
\alias{corSAR1}
\alias{corAR1-class}
\alias{corAR1}
\alias{corSAR1AR1-class}
\alias{corSAR1AR1}
\title{Correlation Structure}
\usage{
corSAR1(W)

corAR1(nTime)

corSAR1AR1(nTime, W)
}
\arguments{
\item{W}{the row-standardised proximity matrix}

\item{nTime}{(numeric) number of time periods}
}
\description{
Various correlation structures. They can be used inside the \link{rfh}
function to supply an alterantive variance structure to be fitted. For
examples see the documentation of \link{rfh}.
}
\details{
\code{corSAR1} can be used to model a simultanous autoregressive
  process of order one: spatial correlation.

\code{corAR1} can be used to model a autoregressive
  process of order one: temporal correlation.

\code{corSAR1AR1} can be used to model to model spatial and temporal
  correlation
}
\section{Slots}{

\describe{
\item{\code{W}}{the row-standardised proximity matrix}

\item{\code{nTime}}{(numeric) number of time periods}
}}

