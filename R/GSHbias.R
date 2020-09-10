#' @title Bias of shinkage rule under GSH prior.
#' @description Provides the bias of the wavelet shrinkage rule under GSH prior.
#' @param theta Wavelet coefficients vector.
#' @param alpha Weight of the point mass at zero function of the prior.
#' @param tau Scale parameter of the GSH prior.
#' @param t Kurtosis parameter of the GSH prior.
#' @param s Standard deviation of the normal random noise.
#'
#' @return Vector of bias of the shrinkage rule.
#' @export
#'
#' @examples GSHbias(c(0,1,2),0.9,1,2,1)

GSHbias = function(theta,alpha,tau,t,s){

  n = length(theta)

  GSHbias = NA

  for(i in 1:n){

    integrand = function(z){

      GSHshrink(theta[i]+s*z,alpha,tau,t,s)*dnorm(z)

    }

    expec = integrate(integrand,lower = -20, upper = 20, stop.on.error = FALSE)$value

    GSHbias[i] = theta[i] - expec

  }

  GSHbias

}
