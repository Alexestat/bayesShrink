#' @title Wavelet shrinkage under GSH prior.
#' @description Performs bayesian shrinkage under GSH prior on empirical wavelet coefficients.
#' @param d Empirical wavelet coefficients vector.
#' @param alpha Weight of the point mass at zero function of the prior.
#' @param tau Scale parameter of the GSH prior.
#' @param t Kurtosis parameter of the GSH prior.
#' @param s Standard deviation of the normal random noise.
#'
#' @return Shrunk wavelet coefficients vector.
#' @export
#'
#' @examples GSHshrink(c(0.5,1,2),0.9,1,2,1)
#'

GSHshrink = function(d,alpha,tau,t,s){

  GSH = function(x,tau,t){
    GSH = NA
    n = length(x)

    if(t>-pi && t<0){
      a = cos(t)
      c2 = sqrt((pi^2-t^2)/3)
      c1 = c2*sin(t)/t

    }

    if(t>0){
      a = cosh(t)
      c2 = sqrt((pi^2+t^2)/3)
      c1 = c2*sinh(t)/t
    }

    for(i in 1:n){
      u = c2*x[i]/tau
      GSH[i] = (c1/tau)*exp(u)/(exp(2*u)+2*a*exp(u)+1)
    }
    GSH
  }


  n = length(d)

  u = rnorm(1000)

  GSHshrink=NA

  for(i in 1:n){

    v = s*u+d[i]

    int1 = mean(v*GSH(v,tau,t))

    int2 = mean(GSH(v,tau,t))

    GSHshrink[i] = (1-alpha)*int1/(alpha*dnorm(d[i],0,s)+(1-alpha)*int2)

  }
  GSHshrink

}
