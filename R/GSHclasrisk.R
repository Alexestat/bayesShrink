#' @title Classical risk of shinkage rule under GSH prior.
#' @description Provides the classical risk of the wavelet shrinkage rule under GSH prior.
#' @param theta Wavelet coefficients vector.
#' @param alpha Weight of the point mass at zero function of the prior.
#' @param tau Scale parameter of the GSH prior.
#' @param t Kurtosis parameter of the GSH prior.
#' @param s Standard deviation of the normal random noise.
#'
#' @return Vector of classical risks of the shrinkage rule.
#' @export
#'
#' @examples GSHclasrisk(c(0,1,2),0.9,1,2,1)

GSHclasrisk = function(theta,alpha,tau,t,s){

  n = length(theta)

  GSHclasrisk = NA

  u = rnorm(1000)

  for(i in 1:n){

    GSHclasrisk[i] = mean((GSHshrink(theta[i]+u,alpha,tau,t,s)-theta[i])^2)

  }
  GSHclasrisk
}
