dQdgamma<-function(m, tp, sigma2omega, GammaInv, W){
  -m*tp*sum(diag(GammaInv))+sigma2omega^(-1)*sum(diag(GammaInv %*% GammaInv %*% W))
}