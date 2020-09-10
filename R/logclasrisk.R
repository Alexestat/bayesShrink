#' @title Classical risk of shinkage rule under logistic prior.
#' @description Provides the classical risk of the wavelet shrinkage rule under logistic prior.
#' @param theta Wavelet coefficients vector.
#' @param alpha Weight of the point mass at zero function of the prior.
#' @param t Scale parameter of the logistic prior.
#' @param s Standard deviation of the normal random noise.
#'
#' @return Vector of classical risks of the shrinkage rule.
#' @export
#'
#' @examples logclasrisk(c(0,1,2),0.9,1,1)

logclasrisk = function(theta,alpha,t,s){

  n = length(theta)

  logclasrisk = NA

  for(i in 1:n){

    integrand = function(z){

      (theta[i] - logshrink(theta[i]+s*z,alpha,t,s))^2*dnorm(z)
    }


    logclasrisk[i] = integrate(integrand,lower = -10, upper =10)$value
  }
  logclasrisk
}
