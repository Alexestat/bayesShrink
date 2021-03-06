% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logshrink.R
\name{logshrink}
\alias{logshrink}
\title{Wavelet shrinkage under logistic prior.}
\usage{
logshrink(d, alpha, t, s)
}
\arguments{
\item{d}{The empirical wavelet coefficients vector.}

\item{alpha}{The weight of the point mass at zero function of the prior.}

\item{t}{The scale parameter of the logistic prior.}

\item{s}{The standard deviation of the normal random noise.}
}
\value{
The shrunk wavelet coefficients vector.
}
\description{
Performs bayesian shrinkage under logistic prior on empirical wavelet coefficients.
}
\examples{
logshrink(c(0.5,1,2),0.9,1,1)

}
