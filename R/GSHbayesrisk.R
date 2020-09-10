#' @title Bayesian risk of shinkage rule under GSH prior.
#' @description Provides the bayesian risk of the wavelet shrinkage rule under GSH prior.
#' @param alpha Weight of the point mass at zero function of the prior.
#' @param tau Scale parameter of the GSH prior.
#' @param t Kurtosis parameter of the GSH prior.
#' @param s Standard deviation of the normal random noise.
#'
#' @return Bayesian risk value of the shrinkage rule.
#' @export
#'
#' @examples GSHbayesrisk(0.9,1,2,1)

GSHbayesrisk = function(alpha,tau,t,s){

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

  u = runif(100,-100,100)

  int = 200*mean(GSHclasrisk(u,alpha,tau,t,s)*GSH(u,tau,t))

  GSHbayesrisk = (alpha)*GSHclasrisk(0,alpha,tau,t,s) + (1-alpha)*int

  GSHbayesrisk

}
