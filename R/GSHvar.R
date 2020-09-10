#' @title Variance of shinkage rule under GSH prior.
#' @description Provides the variance of the wavelet shrinkage rule under GSH prior.
#' @param theta Wavelet coefficients vector.
#' @param alpha Weight of the point mass at zero function of the prior.
#' @param tau Scale parameter of the GSH prior.
#' @param t Kurtosis parameter of the GSH prior.
#' @param s Standard deviation of the normal random noise.
#'
#' @return Vector of variances of the shrinkage rule.
#' @export
#'
#' @examples GSHvar(c(0,1,2),0.9,1,2,1)

GSHvar = function(theta,alpha,tau,t,s){

  n = length(theta)

  u = rnorm(2000)

  GSHvar = NA

  for(i in 1:n){
    GSHvar[i] = var(GSHshrink(theta[i]+u,alpha,tau,t,s))

  }

  GSHvar

}
