dQds<-function(m, tp, GammaInv, sigma2omega, W, d, s, nu) {
  dGammads_j<-dGammads(nu, s, d)
  
  t1<-sum(diag(GammaInv %*% dGammads_j))
  t2<-sum(diag(GammaInv %*% dGammads_j %*% GammaInv %*% W))
    
  -m*tp*t1 + (1/sigma2omega)*t2
}