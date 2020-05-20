d2Qdnuds<-function(m, tp, GammaInv, sigma2omega, W, d, s, nu) {
  dGammads_j<-dGammads(nu, s, d)
  dGammadnu_j<-dGammadnu(nu, s, d)
  dGammadnuds_j<-dGammadnuds(nu, s, d)
    
  t1<-sum(diag(GammaInv %*% dGammadnuds_j - GammaInv %*% dGammads_j %*% GammaInv %*% dGammadnu_j))
  t2<-sum(diag(GammaInv %*% dGammadnuds_j %*% GammaInv %*% W))
  t3<-sum(diag(GammaInv %*% dGammadnu_j %*% GammaInv %*% dGammads_j %*% GammaInv %*% W))
 
  -m*tp*t1+(1/sigma2omega)*t2-(2/sigma2omega)*t3   
}