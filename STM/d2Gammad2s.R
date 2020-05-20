d2Gammad2s<-function(nu, s, d) {
  p1<-gamma(nu+1) #constant
  p2<-apply(d, 2, function(x) {
    ifelse(x==0, 0, (2/x)^nu)
  })
    
  p3<-(s^(nu))*d2s.bessel.mat(nu, s, d)
  p4<-2*nu*s^(nu-1)*ds.bessel.mat(nu, s, d)
  p5<-nu*(nu-1)*s^(nu-2)*bessel.mat(nu, s, d)
  
  out.01<-lapply(1:ncol(d), function(j) {
    j.p1<-p1
    j.p2<-p2[,j]
    j.p3<-p3[,j]
    j.p4<-p4[,j]
    j.p5<-p5[,j]  
    
    j.p1*j.p2*(j.p3+j.p4+j.p5)
  })
  out<-do.call(cbind,out.01)
  return(out)
}