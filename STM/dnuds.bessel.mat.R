dnuds.bessel.mat<-function(nu, s, d) {
  p1<-apply(d, 2, function(x) {
    ifelse(x==0, 0, log(x/(2*s)))
  })
  
  p2<-ds.bessel.mat(nu, s, d)
  p3<-(1/s)*bessel.mat(nu, s, d)
  
  k<-seq(0, 100, 1)
  p4<-apply(d, 2, function(x) {
    out.01<-lapply(1:length(x), function(i) {
      xx<-x[i]
      if(xx==0) {
        return(0)
      } else {
        sum((((-1)^(k))/(factorial(k)*gamma(k+nu+1)))*(digamma(k+nu+1))*(xx/2)^(2*k+nu)*(2*k+nu)*(1/s)^(2*k+nu+1))
      }
    })
    unlist(out.01)
  })
  
  out.01<-lapply(1:ncol(d), function(j) {
    j.p1<-p1[,j]
    j.p2<-p2[,j]
    j.p3<-p3[,j]
    j.p4<-p4[,j]
    j.p1*j.p2-j.p3+j.p4
  })
  out<-do.call(cbind,out.01)
  return(out)
  
}