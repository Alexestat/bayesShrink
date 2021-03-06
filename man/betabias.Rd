% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/betabias.R
\name{betabias}
\alias{betabias}
\title{Bias of shinkage rule under beta prior.}
\usage{
betabias(theta, alpha, a, b, m, s)
}
\arguments{
\item{theta}{Wavelet coefficients vector.}

\item{alpha}{Weight of the point mass at zero function of the prior.}

\item{a}{Shape parameter of the beta prior.}

\item{b}{Shape parameter of the beta prior.}

\item{m}{Upper value of the beta prior support.}

\item{s}{Standard deviation of the normal random noise.}
}
\value{
Vector of bias of the shrinkage rule.
}
\description{
Provides the bias of the wavelet shrinkage rule under beta prior.
}
\examples{
betabias(c(0,1,2),0.9,2,3,10,1)
}
