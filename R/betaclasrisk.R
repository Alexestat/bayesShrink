
#' @title Classical risk of shinkage rule under beta prior.
#' @description Provides the classical risk of the wavelet shrinkage rule under beta prior.
#' @param theta Wavelet coefficients vector.
#' @param alpha Weight of the point mass at zero function of the prior.
#' @param a Shape parameter of the beta prior.
#' @param b Shape parameter of the beta prior.
#' @param m Upper value of the beta prior support.
#' @param s Standard deviation of the normal random noise.
#'
#' @return Vector of classical risks of the shrinkage rule.
#' @export
#'
#' @examples betaclasrisk(c(0,1,2),0.9,2,3,10,1)
betaclasrisk = function(theta,alpha,a,b,m,s){

  n = length(theta)

  betaclasrisk = NA

  for(i in 1:n){

    integrand = function(z){

      (theta[i] - betashrink(theta[i]+s*z,alpha,a,b,m,s))^2*dnorm(z)
    }


    betaclasrisk[i] = integrate(integrand,lower = -10, upper =10)$value
  }
  betaclasrisk
}
