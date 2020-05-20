getK<-function(Xmat, A, r) {
  ################################################################
  #INPUT
  # X:  a list of length m containing an n(x)p matrix of covariates
  # A:  n(x)n adjacency matrix
  # r:
  #OUPUT
  #
  ################################################################
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
  
  # #Find optimal number of loadings
  # SD.values<-sapply(SD, "[[", 2)
  # SD.cp<-apply(SD.values, 2, function(x) {
  #   perc.var<-x/sum(x)
  #   cum.var<-cumsum(perc.var)
  #   return(cum.var)
  # })
  # 
  # ri<-apply(SD.cp, 2, function(x) {
  #   min(which(x>=cp))
  # })
  # 
  # r<-max(ri)
  
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



