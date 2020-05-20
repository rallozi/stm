d2Qdgammads<-function(m, tp, GammaInv, sigma2omega, W, d, s, nu) {
  dGammads_j<-dGammads(nu, s, d)
  
  t1<-sum(diag(GammaInv %*% dGammads_j %*% GammaInv))
  t2<-sum(diag(GammaInv %*% GammaInv %*% dGammads_j %*% GammaInv %*% W))
  
  m*tp*t1-(2/sigma2omega)*t2
}