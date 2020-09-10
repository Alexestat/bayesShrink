#' @title Bayesian risk of shinkage rule under beta prior.
#' @description Provides the bayesian risk of the wavelet shrinkage rule under beta prior.
#' @param alpha Weight of the point mass at zero function of the prior.
#' @param a Shape parameter of the beta prior.
#' @param b Shape parameter of the beta prior.
#' @param m Upper value of the beta prior support.
#' @param s Standard deviation of the normal random noise.
#'
#' @return The bayesian risk value of the shrinkage rule.
#' @export
#'
#' @examples betabayesrisk(0.9,2,3,10,1)

betabayesrisk = function(alpha,a,b,m,s){

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

  integrand = function(theta){

    betaclasrisk(theta,alpha,a,b,m,s)*betadist(theta,a,b,m)

  }

  int = integrate(integrand,lower = -m, upper = m)$value

  betabayesrisk = (alpha)*betaclasrisk(0,alpha,a,b,m,s) + (1-alpha)*int

  betabayesrisk

}
