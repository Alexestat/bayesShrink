
#' @title Wavelet shrinkage under logistic prior.
#' @description Performs bayesian shrinkage under logistic prior on empirical wavelet coefficients.
#' @param d The empirical wavelet coefficients vector.
#' @param alpha The weight of the point mass at zero function of the prior.
#' @param t The scale parameter of the logistic prior.
#' @param s The standard deviation of the normal random noise.
#'
#' @return The shrunk wavelet coefficients vector.
#' @export
#'
#' @examples logshrink(c(0.5,1,2),0.9,1,1)
#'
logshrink = function(d,alpha,t,s){

  n = length(d)

  logshrink = NA

  u = rnorm(10000)

  for(i in 1:n){

    x=(s*u+d[i])*(cosh((s*u+d[i])/(2*t)))^(-2)
    int1 = mean(x)

    y=(cosh((s*u+d[i])/(2*t)))^(-2)
    int2 = mean(y)

    num = (1-alpha)*int1

    den = 4*t*alpha*dnorm(d[i],0,s)/s + (1-alpha)*int2

    logshrink[i] = num/den

  }

  logshrink

}
