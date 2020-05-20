dSigmetadrho<-function(rho, r) {
  H <- abs(outer(1:r, 1:r, "-"))
  V <- H*rho^(H-1)
  V
}


dQdrho<-function(m, tp, r, SigmaetaInv, sigma2eta, rho, E) {
  dSigmetadrho_j<-dSigmetadrho(rho=rho, r=r)
  
  trace.1<-sum(diag(SigmaetaInv %*% dSigmetadrho_j))
  trace.2<-sum(diag(SigmaetaInv %*% dSigmetadrho_j %*% SigmaetaInv %*% E))
  
  -m*tp*trace.1+(1/sigma2eta)*trace.2
}