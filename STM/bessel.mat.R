bessel.mat<-function(nu, s, d) {
  k<-seq(0,100,1)
  out<-apply(d, 2, function(x) {
    out.01<-lapply(1:length(x), function(i) {
      xx<-x[i]
      if(xx==0) {
        return(0)
      } else {
        sum((((-1)^k)/(factorial(k)*gamma(k+nu+1)))*(xx/(2*s))^(2*k+nu))
      }
    })
    unlist(out.01)
  })
  return(out)
}