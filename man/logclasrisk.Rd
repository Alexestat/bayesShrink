% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logclasrisk.R
\name{logclasrisk}
\alias{logclasrisk}
\title{Classical risk of shinkage rule under logistic prior.}
\usage{
logclasrisk(theta, alpha, t, s)
}
\arguments{
\item{theta}{Wavelet coefficients vector.}

\item{alpha}{Weight of the point mass at zero function of the prior.}

\item{t}{Scale parameter of the logistic prior.}

\item{s}{Standard deviation of the normal random noise.}
}
\value{
Vector of classical risks of the shrinkage rule.
}
\description{
Provides the classical risk of the wavelet shrinkage rule under logistic prior.
}
\examples{
logclasrisk(c(0,1,2),0.9,1,1)
}
