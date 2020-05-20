dGammadnu<-function(nu, s, d) {
  p1<-apply(d, 2, function(x) {
    ifelse(x==0, 0, ((2*s)/(x))^(nu))
  })
  
  p2<-gamma(nu+1) #constant
  p3<-dnu.bessel.mat(nu, s, d)
  p4<-bessel.mat(nu, s, d)
  p5<-digamma(nu+1) #constant
  p6<-apply(d, 2, function(x) {
    ifelse(x==0, 0, log((2*s)/(x)))
  })
    
  out.01<-lapply(1:ncol(d), function(j) {
    j.p1<-p1[,j]
    j.p2<-p2
    j.p3<-p3[,j]
    j.p4<-p4[,j]
    j.p5<-p5
    j.p6<-p6[,j]
    
    j.p1*j.p2*(j.p3+j.p4*(j.p5+j.p6))
  })
  
  out<-do.call(cbind, out.01)
  return(out)
}