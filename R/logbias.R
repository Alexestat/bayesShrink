

#' @title Bias of shinkage rule under logistic prior.
#' @description Provides the bias of the wavelet shrinkage rule under logistic prior.
#' @param theta Wavelet coefficients vector.
#' @param alpha Weight of the point mass at zero function of the prior.
#' @param t Scale parameter of the logistic prior.
#' @param s Standard deviation of the normal random noise.
#'
#' @return Vector of bias of the shrinkage rule.
#' @export
#'
#' @examples logbias(c(0,1,2),0.9,1,1)
logbias = function(theta,alpha,t,s){

  n = length(theta)

  logbias = NA

  for(i in 1:n){

    integrand = function(z){

      logshrink(theta[i]+s*z,alpha,t,s)*dnorm(z)

    }

    expec = integrate(integrand,lower = -10, upper = 10, stop.on.error = FALSE)$value

    logbias[i] = theta[i] - expec

  }

  logbias

}
