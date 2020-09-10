#' @title Bayesian risk of shinkage rule under logistic prior.
#' @description Provides the bayesian risk of the wavelet shrinkage rule under logistic prior.
#' @param alpha Weight of the point mass at zero function of the prior.
#' @param t Scale parameter of the logistic prior.
#' @param s Standard deviation of the normal random noise.
#'
#' @return Bayesian risk value of the shrinkage rule.
#' @export
#'
#' @examples logbayesrisk(0.9,1,1)

logbayesrisk = function(alpha,t,s){

  u = rlogis(100,0,t)

  x = logclasrisk(u,alpha,t,s)

  logbayesrisk = (alpha)*logclasrisk(0,alpha,t,s) + (1-alpha)*mean(x)

  logbayesrisk

}

