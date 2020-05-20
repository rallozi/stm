d2Qdgammadnu<-function(m, tp, GammaInv, sigma2omega, W, d, s, nu) {
  dGammadnu_j<-dGammadnu(nu, s, d)
  
  t1<-sum(diag(GammaInv %*% dGammadnu_j %*% GammaInv))
  t2<-sum(diag(GammaInv %*% GammaInv %*% dGammadnu_j %*% GammaInv %*% W))
  
  m*tp*t1-(2/sigma2omega)*t2
}