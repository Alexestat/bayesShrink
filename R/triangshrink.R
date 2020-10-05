#' @title Wavelet shrinkage under triangular prior.
#' @description Performs bayesian shrinkage under triangular prior on empirical wavelet coefficients.
#' @param d The empirical wavelet coefficients vector.
#' @param alpha The weight of the point mass at zero function of the prior.
#' @param m The upper value of the beta prior support.
#' @param s The standard deviation of the normal random noise.
#' @return The shrunk wavelet coefficients vector.
#' @export
#' @examples triangshrink(c(0.5,1,2),0.9,10,1)
triangshrink = function(d,alpha,m,s){

  library("triangle")
  n = length(d)

  triangshrink = NA

  for(i in 1:n){

    integrand1 = function(theta){
      theta*dtriangle(theta,-m,m)*dnorm(d[i],theta,s)
    }

    integrand2 = function(theta){
      dtriangle(theta,-m,m)*dnorm(d[i],theta,s)
    }


    num = (1-alpha)*integrate(integrand1,lower = -m, upper = m,stop.on.error = FALSE)$value

    den = alpha*dnorm(d[i],0,s) + (1-alpha)*integrate(integrand2,lower = -m,upper = m,stop.on.error = FALSE)$value

    triangshrink[i] = num/den
    if(abs(triangshrink[i])>m) {triangshrink[i] = d[i]}
  }
  triangshrink

}
