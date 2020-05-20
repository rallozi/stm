d2Qd2gamma<-function(m, tp, GammaInv, sigma2omega, W) {
  m*tp*sum(diag(GammaInv %*% GammaInv)) - (2/sigma2omega)*sum(diag(GammaInv %*% GammaInv %*% GammaInv %*% W))
}