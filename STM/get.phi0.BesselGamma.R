################################################
##initial values for bessel spatial covariance##
################################################
library(MASS)
library(e1071)

get.phi0<-function(z, X, K, SigmaEtaForm="Unstructured", bessel.moments) {
  n<-ncol(z[[1]])
  tp<-nrow(z[[1]])
  m<-length(z)
  r<-ncol(K[[1]])
  p<-ncol(X[[1]])
    
  #beta.0
  beta.0.01<-lapply(1:n, function(j){
    dat.n<-lapply(1:m, function(i){
      data.frame(X[[i]], z=z[[i]][j,], id=i)
    })
    dat.all<-do.call(rbind, dat.n)
    myform<-as.formula(paste("z~", paste(colnames(X[[1]])[2:ncol(X[[1]])], collapse="+"), sep=""))
    myfit <- lm(myform, dat.all)
    coef(myfit)
  })
  beta.0.02<-do.call(rbind, beta.0.01)
  beta.0<-apply(beta.0.02, 2, mean)
  
  #estimate latent variables for each individual i across all time points
  #returns a list of length m, with each individual's list containing a list of length tp
  y.i.01<-lapply(1:m, function(i) {
    X.i<-X[[i]] #n(x)p
    K.i<-K[[i]] #n(x)r
    y.itt<-lapply(1:tp, function(tt) {
      z.itt<-z[[i]][tt,] #ith individual, data observed at tt'th time point
      y.itt<-ginv(K.i) %*% (z.itt-X.i %*% beta.0)
    })
    return(y.itt)
  })
  
  y.i<-lapply(1:m, function(i) {
    do.call(cbind, y.i.01[[i]])
  })
  
  #get mean y at each time point
  y.est<-lapply(1:tp, function(tt) {
    y.tt<-lapply(1:m, function(i) {
      y.i[[i]][,tt]
    })
    y.tt<-do.call(rbind, y.tt)
    mu.tt<-apply(y.tt, 2, mean)
    P.tt<-cov(y.tt)
    return(list(mu=mu.tt, P=P.tt))
  })
  
  #m0, C0
  m0.0<-y.est[[1]]$mu
  C0.0<-y.est[[1]]$P
  
  
  #G0
  G.0.01<-lapply(2:tp, function(tt) {
    y.tt<-lapply(1:m, function(i) {
      y.i[[i]][,tt]
    })
    y.tt<-do.call(rbind, y.tt)
    
    y.ttmin1<-lapply(1:m, function(i) {
      y.i[[i]][,tt-1]
    })  
    y.ttmin1<-do.call(rbind,  y.ttmin1)
    
    all.coef<-list()
    for(i in 1:r) {
      y<-y.tt[,i]
      x<-y.ttmin1
      fit<-lm(y~x-1)
      all.coef[[i]]<-coef(fit)
    }
    
    G.est<-do.call(cbind, all.coef)
    return(G.est)
  })
  G.0<-Reduce('+', G.0.01)/(tp-1)
  
  #SigmaEta0
  eta.i<-lapply(2:tp, function(tt) {
    y.tt<-lapply(1:m, function(i) {
      y.i[[i]][,tt]
    })
    y.tt<-do.call(rbind, y.tt)
    
    y.ttmin1<-lapply(1:m, function(i) {
      y.i[[i]][,tt-1]
    })  
    y.ttmin1<-do.call(rbind,  y.ttmin1)
    
    eta.i<-y.tt-t(G.0 %*% t(y.ttmin1))
    return(eta.i)
  })
  
  mu.eta<-sum(unlist(eta.i))/(m*(tp-1))
  
  SigmaEta.0.01<-lapply(1:length(eta.i), function(tt) {
    mymat<-eta.i[[tt]]
    SigmaEta.i<-lapply(1:m, function(i){
      (mymat[i,]-mu.eta) %*% t((mymat[i,]-mu.eta) )
    })
    sumMatrices(SigmaEta.i)
  })
  SigmaEta.0<-sumMatrices(SigmaEta.0.01)/(m*(tp-1))
  
  if(SigmaEtaForm=="AR(1)") {
    sigma2eta.0<-mean(diag(SigmaEta.0))
    
    H <- abs(outer(1:r, 1:r, "-"))
    V<-(abs(SigmaEta.0/sigma2eta.0))^(1/H)
    
    rho.0<-mean(V[upper.tri(V)])
  } else if(SigmaEtaForm=="Diagonal") {
    sigma2eta.0<-mean(diag(SigmaEta.0))
  }

  #sigma2eps, sigma2omega
  mu.xi.t<-lapply(1:tp, function(tt) {
    xi.it<-lapply(1:m, function(i) {
      X.i<-X[[i]]
      K.i<-K[[i]]
      xi.it<-z[[i]][tt,] - X.i %*% beta.0 - K.i %*% y.i[[i]][,tt]
      return(xi.it)
    })
    apply(do.call(cbind, xi.it), 1, mean)
  })
  
  SigmaXi.01<-lapply(1:m, function(i) {
    X.i<-X[[i]]
    K.i<-K[[i]]
    
    sigma2xi.i<-lapply(1:tp, function(tt) {
      xi.it<-z[[i]][tt,] - X.i %*% beta.0 - K.i %*% y.i[[i]][,tt]
      xi.it %*% t(xi.it) - matrix(mu.xi.t[[tt]]) %*% t(matrix(mu.xi.t[[tt]]))
    })  
    out<-sumMatrices(sigma2xi.i)
    return(out)
  })
  SigmaXi.0<-sumMatrices(SigmaXi.01)/(m*tp)
  
  sigma2eps.0<-sigma2omega.0<-mean(diag(SigmaXi.0))/2
  
  #nu, s
  SigmaGamma.0<-SigmaXi.0/sigma2omega.0
  UT.SigmaGamma.0<-SigmaGamma.0[upper.tri(SigmaGamma.0)]
  
  UT.mean<-mean(UT.SigmaGamma.0)
  UT.sd<-sd(UT.SigmaGamma.0)
  UT.sk<-skewness(UT.SigmaGamma.0)
  UT.k<-kurtosis(UT.SigmaGamma.0)
  
  UT.diff<-(UT.mean-bessel.moments$mean)^2+(UT.sd-bessel.moments$sd)^2+(UT.sk-bessel.moments$skewness)^2+(UT.k-bessel.moments$kurtosis)^2
  s.0<-bessel.moments[which(UT.diff==min(UT.diff)),c("scale")]
  nu.0<-bessel.moments[which(UT.diff==min(UT.diff)),c("nu")]
  
  if(SigmaEtaForm=="Unstructured") {
    phi <- list(beta=beta.0,
                sigma2eps=sigma2eps.0,
                sigma2omega=sigma2omega.0,
                nu=nu.0,
                s=s.0,
                G=G.0,
                Sigmaeta=SigmaEta.0,
                m0=matrix(m0.0, r, 1),
                C0=C0.0)
  } else if(SigmaEtaForm=="AR(1)") {
    phi <- list(beta=beta.0,
                sigma2eps=sigma2eps.0,
                sigma2omega=sigma2omega.0,
                nu=nu.0,
                s=s.0,
                G=G.0,
                sigma2eta=sigma2eta.0,
                rho=rho.0,
                m0=matrix(m0.0, r, 1),
                C0=C0.0)
  } else if(SigmaEtaForm=="Diagonal") {
    phi <- list(beta=beta.0,
                sigma2eps=sigma2eps.0,
                sigma2omega=sigma2omega.0,
                nu=nu.0,
                s=s.0,
                G=G.0,
                sigma2eta=sigma2eta.0,
                m0=matrix(m0.0, r, 1),
                C0=C0.0)
  }

  return(phi)
}