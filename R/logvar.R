#' @title Variance of shinkage rule under logistic prior.
#' @description Provides the variance of the wavelet shrinkage rule under logistic prior.
#' @param theta Wavelet coefficients vector.
#' @param alpha Weight of the point mass at zero function of the prior.
#' @param t Scale parameter of the logistic prior.
#' @param s Standard deviation of the normal random noise.
#'
#' @return Vector of variance of the shrinkage rule.
#' @export
#'
#' @examples logvar(c(0,1,2),0.9,1,1)

logvar = function(theta,alpha,t,s){

  n = length(theta)

  logvar = NA

  for(i in 1:n){

    integrand1 = function(z){

      logshrink(theta[i]+s*z,alpha,t,s)*dnorm(z)

    }

    expec = integrate(integrand1,lower = -50, upper = 50)$value

    integrand2 = function(z){

      logshrink(theta[i]+s*z,alpha,t,s)^2*dnorm(z)

    }

    mom2 = integrate(integrand2,lower = -50, upper = 50)$value

    logvar[i] = mom2 - expec^2

  }

  logvar

}
