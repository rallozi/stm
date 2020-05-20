source(file.path(source.path, "gen.univ.X.R"))
source(file.path(source.path, "getK.R"))
source(file.path(stem.path, "Gamma.Cov.R"))

library(corpcor)
library(Matrix)

`Stem.Simulation` <-
function (n, tp, r, p, m, X, phi, FD, SigmaEta){

  #n=number of spatial locations
  #tp=number of time points
  #r=dimension of latent variable
  #p=number of fixed-effects parameters (including intercept)
  #m=number of subjects
  #X=a list of length p-1 containing information for covariates
  #phi=true parameters
  #FD=functional distance
  
  phi$b=phi$sigma2eps/phi$sigma2omega #Setting 2

  ##Simulate covariates
  if(length(X)!=(p-1)) {
    stop("Length of X does not match p-1.")
  }
  
  #n by p matrix
  allX<-lapply(1:m, function(i) {
    Xi.01<-lapply(X, function(x) {
      gen.univ.X(n=1, dist=x$dist, param=x[2:length(x)])
    })
    
    Xi.02<-do.call('cbind', Xi.01)
    Xi.03<-cbind(1, Xi.02)
    Xi.04<-matrix(rep(Xi.03, n), nrow=n, ncol=p, byrow=TRUE)
    colnames(Xi.04)<-paste('X', 0:(p-1), sep='')
    
    Xi<-lapply(1:tp, function(tt) {
      Xi.04
    })
    return(Xi)
  })
  names(allX)<-1:m
  
  Xi<-lapply(1:m, function(i) {
    allX[[i]][[1]]
  })
  names(Xi)<-1:m
  
  A<-exp(-0.5*FD)
  allK<-getK(Xmat=Xi, A=A, r=r)
  names(allK)<-1:m
  
  allG<-lapply(1:m, function(i) {
    phi$G
  })
  
  ####################
  ##Variance Matries##
  ####################
  SigmaXi=phi$sigma2omega * Gamma.Cov(n=n, b=phi$b, d=FD, nu=phi$nu, s=phi$s) #Setting 2
  if(SigmaEta=="Unstructured") {
    Wmat=phi$Sigmaeta                           
  } else if(SigmaEta=="AR(1)") {
    if(r==1) {
      Sigmaeta<-1
    } else {
      H <- abs(outer(1:r, 1:r, "-"))
      V <- phi$rho^H
      V[cbind(1:r, 1:r)] <- V[cbind(1:r, 1:r)]
      Sigmaeta=V   
    }
    Wmat=phi$sigma2eta*Sigmaeta #Variance of eta
  } else if(SigmaEta=="Diagonal") {
    Wmat=phi$sigma2eta*diag(r)
  }

  ##########################################
  ##Starting Parameters for Latent Prcoess##
  ##########################################
  m0   = phi$m0                                 
  C0   = phi$C0		                     

  #####################
  allz<-lapply(1:m, function(i) {
    z = matrix(NA,nrow = n,ncol = tp)
    y = matrix(NA,nrow = r,ncol = tp)
    
    K<-allK[[i]]
    G<-allG[[i]]
    X<-allX[[i]]
    
    y0 = matrix(mvrnorm(n=1, mu=m0, Sigma=C0),nrow=r, ncol=1)
    
    #time=1
    y[,1] = G %*% y0 + mvrnorm(n=1, mu=matrix(0,nrow=r,ncol=1), Sigma=Wmat)
    z[,1] = X[[1]] %*% phi$beta + K %*% as.matrix(y[,1], r) + mvrnorm(n=1, mu=matrix(0,nrow=n,ncol=1), Sigma=SigmaXi)
    
    for (tt in 2:tp) {
      y[,tt] = G %*% y[,tt-1] + mvrnorm(n=1, mu=matrix(0,nrow=r,ncol=1), Sigma=Wmat)
      z[,tt] = X[[tt]] %*% phi$beta + K %*% y[,tt] + mvrnorm(n=1, mu=matrix(0,nrow=n,ncol=1),Sigma=SigmaXi)
    }
    
    return(t(z))
  })
  names(allz)<-1:m
  return(list(z=allz, X=Xi, K=allK, G=allG))
}

