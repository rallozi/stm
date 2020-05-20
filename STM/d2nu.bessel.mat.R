d2nu.bessel.mat<-function(nu, s, d){
  p1<-apply(d, 2, function(x) {
    ifelse(x==0, 0, log(x/(2*s)))
  })
    
  p2<-dnu.bessel.mat(nu, s, d)

  k<-seq(0,100,1)
  p3<-apply(d, 2, function(x) {
    out.01<-lapply(1:length(x), function(i) {
      xx<-x[i]
      if(xx==0) {
        return(0)
      } else {
        sum( ((((-1)^(k))/(factorial(k)*gamma(k+nu+1)))*(xx/(2*s))^(2*k+nu))*(-digamma(k+nu+1)*log(xx/(2*s))+digamma(k+nu+1)^2-psigamma(k+nu+1, 1)))
      }
    })
    unlist(out.01)
  })
  
  out.01<-lapply(1:ncol(d), function(j) {
    j.p1<-p1[,j]
    j.p2<-p2[,j]
    j.p3<-p3[,j]
    j.p1*j.p2+j.p3
  })
  
  out<-do.call(cbind, out.01)
  return(out)

}