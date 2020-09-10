#' @title Variance of shinkage rule under beta prior.
#' @description Provides the variance of the wavelet shrinkage rule under beta prior.
#' @param theta Wavelet coefficients vector.
#' @param alpha Weight of the point mass at zero function of the prior.
#' @param a Shape parameter of the beta prior
#' @param b Shape parameter of the beta prior
#' @param m Upper value of the beta prior support.
#' @param s Standard deviation of the normal random noise.
#'
#' @return Vector of variances of the shrinkage rule.
#' @export
#'
#' @examples betavar(c(0,1,2),0.9,2,3,10,1)
betavar = function(theta,alpha,a,b,m,s){

  n = length(theta)

  betavar = NA

  for(i in 1:n){

    integrand1 = function(z){

      betashrink(theta[i]+s*z,alpha,a,b,m,s)*dnorm(z)

    }

    expec = integrate(integrand1,lower = -10, upper = 10)$value

    integrand2 = function(z){

      betashrink(theta[i]+s*z,alpha,a,b,m,s)^2*dnorm(z)

    }

    mom2 = integrate(integrand2,lower = -10, upper = 10)$value

    betavar[i] = mom2 - expec^2

  }

  betavar

}
