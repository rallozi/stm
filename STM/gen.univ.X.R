library(rmutil) #rpareto
library(gld) #rgl

gen.univ.X<-function(n, dist, param) {
  ##########################################
  #INPUT:
  #m=number of spatial locations
  #dist=distribution to simulate from
  #param=list containing parameters for dist
  ##########################################
  
  #Simulate data from a normal distribution.
  #param is a vector of length 2 containing mean, variance
  if(dist=="Normal") { 
    if(!all(c('mean', 'variance') %in% names(param))) {
      stop('Must provide a mean and variance for Normal distribution.')
    }
    mu<-param[['mean']]
    sigma<-sqrt(param[['variance']])
    dat<-rnorm(n, mu, sigma)
  } 
  
  #Simulate data from a lognormal distribution.
  #param is a vector of length 2 containing mean, variance
  else if(dist=="Lognormal") {
    if(!all(c('mean', 'variance') %in% names(param))) {
      stop('Must provide a mean and variance for lognormal distribution.')
    }
    mu<-param[['mean']]
    sigma<-sqrt(param[['variance']])
    dat<-rlnorm(n, meanlog=mu, sdlog=sigma)
  } 
  
  #Simulate data from a Gamma distribution.
  #param is a vector of length 2 containing alpha (shape), beta (scale)
  else if(dist=='Gamma') {
    if(!all(c('shape', 'scale') %in% names(param))) {
      stop('Must provide a shape and scale parameter for Gamma distribution.')
    }
    alpha<-param[['shape']]
    beta<-param[['scale']]
    dat<-rgamma(n, shape=alpha, scale=beta)
  } 
  
  #Simulate data from a Pareto distribution
  #param is a vector of length 2 containing location, shape
  else if(dist=='Pareto') {
    if(!all(c('location', 'shape') %in% names(param))) {
      stop('Must provide a location and shape parameter for Pareto distribution.')
    }
    loc<-param[['location']]
    shape<-param[['shape']]
    dat<-rpareto(n=n,m=loc, s=shape)
  }
  
  #Simulate data from a Beta distribution.
  #param is a vector of length 2 containing shape1, shape2
  else if(dist=="Beta") {
    if(!all(c('shape1', 'shape2') %in% names(param))) {
      stop('Must provide a 2 shape parameters for Beta distribution.')
    }
    shape1<-param[['shape1']]
    shape2<-param[['shape2']]
    dat<-rbeta(n, shape=shape1, scale=shape2)
  } 
  
  else if(dist=="Tukey-Lambda") {
    if(!all(c('location', 'scale', 'shape1', 'shape2') %in% names(param))) {
      stop('Must provide location, scale, and 2 shape parameters for Tukey-Lambda distribution.')
    }
    location<-param[['location']]
    scale<-param[['scale']]
    shape1<-param[['shape1']]
    shape2<-param[['shape2']]
    dat<-rgl(n=n, lambda1=location, lambda2=scale, lambda3=shape1, lambda4=shape2)
  } 
  
  else if(dist=="Bernoulli") {
    if(!all(c('p') %in% names(param))) {
      stop('Must provide a p for Bernoulli variable.')
    }
    p<-param[['p']]
    dat<-rbinom(n, 1, p)
  } 
  
  else if(dist=="Poisson") {
    if(!all(c('lambda') %in% names(param))) {
      stop('Must provide lambda for Poisson variable.')
    }
    lam<-param[['lambda']]
    dat<-rpois(n, lam)
  } 
  
  else if(dist=="Uniform") {
    if(!all(c('a', 'b') %in% names(param))) {
      stop('Must provide min, max for Uniform variable.')
    } 
    a<-param[['a']]
    b<-param[['b']]
    dat<-runif(n,a,b)
  }
  
  else {
    stop('dist must be Normal, Lognormal, Gamma, Pareto, Beta, Tukey-Lambda, Bernoulli, Poisson, or Uniform.')
  }
  return(dat)
}
