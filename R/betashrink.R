
#' @title Wavelet shrinkage under beta prior.
#' @description Performs bayesian shrinkage under beta prior on empirical wavelet coefficients.
#' @param d The empirical wavelet coefficients vector.
#' @param alpha The weight of the point mass at zero function of the prior.
#' @param a The shape parameter of the beta prior.
#' @param b The shape parameter of the beta prior.
#' @param m The upper value of the beta prior support.
#' @param s The standard deviation of the normal random noise.
#'
#' @return The shrunk wavelet coefficients vector.
#' @export
#'
#' @examples betashrink(c(0.5,1,2),0.9,2,3,10,1)
betashrink = function(d,alpha,a,b,m,s){

  betadist = function(x,a,b,m){

    betadist = vector()

    for(i in 1:length(x)){

      if(abs(x[i])<=m)
        betadist[i] = (x[i]+m)^(a-1)*(m-x[i])^(b-1)/((2*m)^(a+b-1)*beta(a,b))
      else
        betadist[i]=0

    }

    return(betadist)

  }

  n = length(d)

  betashrink = NA

  for(i in 1:n){

    integrand1 = function(theta){
      theta*betadist(theta,a,b,m)*dnorm(d[i],theta,s)
    }

    integrand2 = function(theta){
      betadist(theta,a,b,m)*dnorm(d[i],theta,s)
    }


    num = (1-alpha)*integrate(integrand1,lower = -m, upper = m,stop.on.error = FALSE)$value

    den = alpha*dnorm(d[i],0,s) + (1-alpha)*integrate(integrand2,lower = -m,upper = m,stop.on.error = FALSE)$value

    betashrink[i] = num/den
    #if(abs(betashrink2[i])>m){betashrink2[i]=x[i]}
  }

  betashrink

}
