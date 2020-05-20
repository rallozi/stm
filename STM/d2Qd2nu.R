d2Qd2nu<-function(m, tp, GammaInv, sigma2omega, W, d, s, nu)  {
  dGammadnu_j<-dGammadnu(nu, s, d)
  d2Gammad2nu_j<-d2Gammad2nu(nu, s, d)
  
  t1<-sum(diag(GammaInv %*% d2Gammad2nu_j - GammaInv %*% dGammadnu_j %*% GammaInv %*% dGammadnu_j))
  t2<-sum(diag(GammaInv %*% d2Gammad2nu_j %*% GammaInv %*% W))
  t3<-sum(diag(GammaInv %*% dGammadnu_j %*% GammaInv %*% dGammadnu_j %*% GammaInv %*% W))
  
  -m*tp*t1+(1/sigma2omega)*t2-(2/sigma2omega)*t3
}