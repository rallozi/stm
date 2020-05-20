dSigmetadrho<-function(rho, r) {
  H <- abs(outer(1:r, 1:r, "-"))
  V <- H*rho^(H-1)
  V
}

d2Sigmetadrho<-function(rho, r) {
  H <- abs(outer(1:r, 1:r, "-"))
  V <- H*(H-1)*rho^(H-2)
  V
}

d2Qd2rho<-function(m, tp, r, SigmaEtaInv, sigma2eta, rho, E){
  dSigmetadrho_j<-dSigmetadrho(rho=rho, r=r)
  d2Sigmetadrho_j<-d2Sigmetadrho(rho=rho, r=r)
  
  trace.1<-sum(diag(SigmaEtaInv %*% d2Sigmetadrho_j - SigmaEtaInv %*% dSigmetadrho_j %*% SigmaEtaInv %*% dSigmetadrho_j))
  trace.2<-sum(diag(SigmaEtaInv %*% d2Sigmetadrho_j %*% SigmaEtaInv %*% E))
  trace.3<-sum(diag(SigmaEtaInv %*% dSigmetadrho_j %*% SigmaEtaInv %*% dSigmetadrho_j %*% SigmaEtaInv %*% E))
  
  -m*tp*trace.1+(1/sigma2eta)*trace.2-(2/sigma2eta)*trace.3
}


