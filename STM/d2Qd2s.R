d2Qd2s<-function(m, tp, GammaInv, sigma2omega, W, d, s, nu) {
  dGammads_j<-dGammads(nu, s, d)
  d2Gammad2s_j<-d2Gammad2s(nu, s, d)
  
  t1<-sum(diag(GammaInv %*% d2Gammad2s_j - GammaInv %*% dGammads_j %*% GammaInv %*% dGammads_j))
  t2<-sum(diag(GammaInv %*% d2Gammad2s_j %*% GammaInv %*% W))
  t3<-sum(diag(GammaInv %*% dGammads_j %*% GammaInv %*% dGammads_j %*% GammaInv %*% W))
  
  -m*tp*t1+(1/sigma2omega)*t2-(2/sigma2omega)*t3
}
