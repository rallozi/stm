####################################
####Functions Used in Simulation####
####################################
eval.sim<-function(res, tp) {
  ############################
  #res: simulation results
  #tp:  true parameters
  #############################
  #put tp in vector form
  tp.vec.01<-list.to.vector(tp)
  tp.vec<-tp.vec.01[which(names(tp.vec.01) %in% colnames(res))]
  
  #Simulation mean, standard deviation
  all.mean<-apply(res, 2, mean)
  all.sd<-apply(res, 2, sd)
  
  #Simulation Bias
  all.bias<-all.mean[colnames(res)]-tp.vec[colnames(res)]
  
  #Simulation RMSE
  all.mse.01<-lapply(colnames(res), function(param) {
    ires<-res[,param]
    itp<-tp.vec[param]
    sqrt(sum((ires-itp)^2)/length(ires))
  })
  all.mse<-unlist(all.mse.01)
  names(all.mse)<-colnames(res)
  
  
  #Center 95%
  get.ci<-lapply(1:length(tp.vec), function(ip) {
    par<-names(tp.vec[ip])
    est<-res[,par]
    ci.01<-quantile(est, probs=c(0.025, 0.975))
    ci.range<-format(round(abs(ci.01[2]-ci.01[1]), 4), nsmall=4)
    ci<-paste0("(", format(round(ci.01[1], 4), nsmall=4), ",", format(round(ci.01[2], 4), nsmall=4), ")")
    return(list(ci, ci.range))
  })
  all.ci<-sapply(get.ci, "[[", 1)
  names(all.ci)<-names(tp.vec)
  all.range<-sapply(get.ci, "[[", 2)
  names(all.range)<-names(tp.vec)
  
  all.msd<-paste0(format.num(all.mean[colnames(res)]), " (", format.num(all.sd[colnames(res)]), ")")
  
  all.res<-rbind(tp.vec[colnames(res)], 
                 all.msd, 
                 format.num(all.mse[colnames(res)]), 
                 format.num(all.bias[colnames(res)]),
                 all.ci[colnames(res)],
                 all.range[colnames(res)])
  
  rownames(all.res)<-c("True Value", 'Sim Mean (SE)', 'Sim RMSE', "Sim Bias", "Middle 95%", "Middle 95% Range")
  return(t(all.res))
}


format.num<-function(x) {
  ifelse(abs(x)>1e-4, format(round(x, 4), nsmall=4), formatC(x, format = "e", digits = 2))
}


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
    stop('dist must be Normal, Lognormal, Gamma, Beta, Bernoulli, Poisson, or Uniform.')
  }
  return(dat)
}


getK<-function(Xmat, A, r) {
  ################################################################
  #INPUT
  # X:  a list of length m containing an n(x)p matrix of covariates
  # A:  n(x)n adjacency matrix
  # r:  latent variable dimension
  n<-nrow(A)
  
  SD<-lapply(Xmat, function(Xmati) {
    XiTXi<-t(Xmati) %*% Xmati #X'X
    
    if(!is.positive.definite(XiTXi)) {
      XiTXi<-as.matrix(nearPD(XiTXi,  ensureSymmetry=TRUE)$mat)
    }
    
    Ei<-diag(n)-Xmati %*% solve(XiTXi) %*% t(Xmati) #I-X(X'X)^(-1)X'
    Mi<-Ei %*% A %*% Ei
    
    SD.i<-eigen(Mi)
    Phi.i<-SD.i$vectors
    lambda.i<-SD.i$values
    
    return(list(E.vectors=Phi.i, E.values=lambda.i)) 
  })

  Ki<-lapply(SD, function(SD.i) {
    Phi.i<-SD.i$E.vectors
    Phi.i.r<-Phi.i[,1:r]
    if(r==1) {
      matrix(Phi.i.r, nrow=n, ncol=1)
    } else {
      Phi.i.r
    }
  })
  return(Ki)
}


list.to.vector<-function(dat) {
  out<-lapply(names(dat), function(param) {
    idat<-dat[[param]]
    if(length(idat)>1) {
      idat.vec<-as.vector(idat)
      names(idat.vec)<-paste0(param, ".", 1:length(idat.vec))
    } else {
      idat.vec<-as.vector(idat)
      names(idat.vec)<-param
    }
    return(idat.vec)
  })
  return(unlist(out)) 
}



################################
####Load STM Model Functions####
################################
library(MASS)
library(mvtnorm)

#stm.path is path to folder containg STM codes
setwd(stm.path)
file.sources = list.files(pattern="*.R")
sapply(file.sources,source,.GlobalEnv)


##################################
###Load data used in simulation###
##################################
load(file.path(dat.path, "FD.Rdata")) #functional distances
bessel.moments<-read.csv(file.path(dat.path, "bessel.moments.n=30.csv"))


##########################################
###Simulation to evaluate ECM procedure###
##########################################
p<-4 #intercept, age, FIQ, DX_GROUP
r<-3
n<-30
tp<-60
SigmaEtaForm<-"Diagonal"

N<-10000
#true model parameters
phi <- list(beta=matrix(c(3.65,0.046,0.904, -2),p,1),
            sigma2eps=0.3,
            sigma2omega=0.5,
            sigma2eta=0.25,
            nu=0.6,
            s=0.2,
            G=diag(0.2, r, r)+matrix(rep(0.2, r*r), r, r),
            m0=as.matrix(rep(0, r), r, 1),
            C0=diag(0.7, r, r)+matrix(rep(0.3, r*r), r, r))

m<-100
est1 <- foreach(i=1:N, .packages=c('corpcor', 'Matrix', 'MASS', 'e1071', 'mvtnorm')) %dopar% {
  sim<-Stem.Simulation(n=n,
                       tp=tp,
                       r=r,
                       p=p,
                       m=m,
                       X=list(list(dist='Normal', mean=26.59, variance=32.33), #Age
                              list(dist='Normal', mean=0, variance=1), #FIQ
                              list(dist='Bernoulli', p=0.5)), #DX_GROUP
                       phi=phi,
                       FD=FD[1:n,1:n],
                       SigmaEta=SigmaEtaForm)
  
  phi0<-get.phi0(z=sim$z, X=sim$X, K=sim$K, SigmaEtaForm=SigmaEtaForm)
  
  mod.it <- Stem.Model(z=sim$z,               #a list of length m containing observation matrices of dimension tp*n
                       covariates=sim$X,      #a list of length m containing covariate matrices of dimension n*p
                       funcdist=FD[1:n,1:n],  #n*n matrix of functional distance
                       phi=phi0,              #starting parameters
                       K=sim$K,               #a list of length m containing loading matrices of dimension n*r
                       r=r,                   #r, dimension of the latent process
                       SigmaEta=SigmaEtaForm)
  
  sim.est <- Stem.Estimation(StemModel=mod.it, 
                             precision=1e-10, 
                             max.iter=2000, 
                             SigmaEta=SigmaEtaForm, 
                             learning.rate=0.5)
  
  sim.phi <- list.to.vector(sim.est$estimates$phi.hat)
  sim.loglik<-sim.est$estimates$loglik
  return(list(phi0=phi0, phi=sim.phi, loglik=sim.loglik, est=sim.est$estimates$phi.hat))
}

est1.phi.01<-do.call(rbind, sapply(est1, "[[", 2))
est1.phi.02<-est1.phi.01[complete.cases(est1.phi.01),]
est1.phi.03<-eval.sim(res=est1.phi.02, tp=phi)
est1.phi<-unique(est1.phi.03)

m<-200
est2 <- foreach(i=1:N, .packages=c('corpcor', 'Matrix', 'MASS', 'e1071', 'mvtnorm')) %dopar% {
  sim<-Stem.Simulation(n=n,
                       tp=tp,
                       r=r,
                       p=p,
                       m=m,
                       X=list(list(dist='Normal', mean=26.59, variance=32.33), #Age
                              list(dist='Normal', mean=0, variance=1), #FIQ
                              list(dist='Bernoulli', p=0.5)), #DX_GROUP
                       phi=phi,
                       FD=FD[1:n,1:n],
                       SigmaEta=SigmaEtaForm)
  
  phi0<-get.phi0(z=sim$z, X=sim$X, K=sim$K, SigmaEtaForm=SigmaEtaForm)
  
  mod.it <- Stem.Model(z=sim$z,               #a list of length m containing observation matrices of dimension tp*n
                       covariates=sim$X,      #a list of length m containing covariate matrices of dimension n*p
                       funcdist=FD[1:n,1:n],  #n*n matrix of functional distance
                       phi=phi0,              #starting parameters
                       K=sim$K,               #a list of length m containing loading matrices of dimension n*r
                       r=r,                   #r, dimension of the latent process
                       SigmaEta=SigmaEtaForm)
  
  sim.est <- Stem.Estimation(StemModel=mod.it, 
                             precision=1e-10, 
                             max.iter=2000, 
                             SigmaEta=SigmaEtaForm, 
                             learning.rate=0.5)
  
  sim.phi <- list.to.vector(sim.est$estimates$phi.hat)
  sim.loglik<-sim.est$estimates$loglik
  return(list(phi0=phi0, phi=sim.phi, loglik=sim.loglik, est=sim.est$estimates$phi.hat))
}

est2.phi.01<-do.call(rbind, sapply(est2, "[[", 2))
est2.phi.02<-est2.phi.01[complete.cases(est2.phi.01),]
est2.phi.03<-eval.sim(res=est2.phi.02, tp=phi)
est2.phi<-unique(est2.phi.03)

m<-400
est3 <- foreach(i=1:N, .packages=c('corpcor', 'Matrix', 'MASS', 'e1071', 'mvtnorm')) %dopar% {
  sim<-Stem.Simulation(n=n,
                       tp=tp,
                       r=r,
                       p=p,
                       m=m,
                       X=list(list(dist='Normal', mean=26.59, variance=32.33), #Age
                              list(dist='Normal', mean=0, variance=1), #FIQ
                              list(dist='Bernoulli', p=0.5)), #DX_GROUP
                       phi=phi,
                       FD=FD[1:n,1:n],
                       SigmaEta=SigmaEtaForm)
  
  phi0<-get.phi0(z=sim$z, X=sim$X, K=sim$K, SigmaEtaForm=SigmaEtaForm)
  
  mod.it <- Stem.Model(z=sim$z,               #a list of length m containing observation matrices of dimension tp*n
                       covariates=sim$X,      #a list of length m containing covariate matrices of dimension n*p
                       funcdist=FD[1:n,1:n],  #n*n matrix of functional distance
                       phi=phi0,              #starting parameters
                       K=sim$K,               #a list of length m containing loading matrices of dimension n*r
                       r=r,                   #r, dimension of the latent process
                       SigmaEta=SigmaEtaForm)
  
  sim.est <- Stem.Estimation(StemModel=mod.it, 
                             precision=1e-10, 
                             max.iter=2000, 
                             SigmaEta=SigmaEtaForm, 
                             learning.rate=0.5)
  
  sim.phi <- list.to.vector(sim.est$estimates$phi.hat)
  sim.loglik<-sim.est$estimates$loglik
  return(list(phi0=phi0, phi=sim.phi, loglik=sim.loglik, est=sim.est$estimates$phi.hat))
}

est3.phi.01<-do.call(rbind, sapply(est3, "[[", 2))
est3.phi.02<-est3.phi.01[complete.cases(est3.phi.01),]
est3.phi.03<-eval.sim(res=est3.phi.02, tp=phi)
est3.phi<-unique(est3.phi.03)


###############################
##Table 2: Simulation Results##
###############################
comb1<-data.frame(Parameter=rownames(est1.phi), m=100, est1.phi, stringsAsFactors=FALSE)
comb2<-data.frame(Parameter=rownames(est2.phi), m=200, est2.phi, stringsAsFactors=FALSE)
comb3<-data.frame(Parameter=rownames(est3.phi), m=400, est3.phi, stringsAsFactors=FALSE)

comb<-rbind(comb1, comb2, comb3)
colnames(comb)<-c("Parameter", "True Value", "m", "Mean (SE)", "RMSE", "Bias", "Middle 95%")
comb<-comb[order(comb$Parameter, comb$m),]

