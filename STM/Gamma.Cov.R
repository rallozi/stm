Gamma.Cov<-function(n, b, d, nu, s) {
  k<-seq(0,100,1)
  Bessel.cov<-apply(d, 2, function(x) {
    out.01<-lapply(1:length(x), function(i) {
      xx<-x[i]
      if(xx==0) {
        return(0)
      } else {
        J<-sum((((-1)^k)/(factorial(k)*gamma(k+nu+1)))*(xx/(2*s))^(2*k+nu))
        return(2^(nu)*gamma(nu+1)*(xx/s)^(-nu)*J)
      }
    })
    unlist(out.01)
  })
  diag(1+b,n) + Bessel.cov
}