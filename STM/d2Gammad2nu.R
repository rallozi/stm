d2Gammad2nu<-function(nu, s, d){
  p1<-apply(d, 2, function(x) {
    ifelse(x==0, 0, ((2*s)/x)^(nu))
  })
  
  p2<-gamma(nu+1) #constant
  p3<-d2nu.bessel.mat(nu, s, d)
  p4<-2*dnu.bessel.mat(nu, s, d)
  p5<-digamma(nu+1) #constant
  p6<-apply(d, 2, function(x) {
    ifelse(x==0, 0, log((2*s)/x))
  })
  
  p7<-bessel.mat(nu, s, d)
  p8<-apply(d, 2, function(x) {
    ifelse(x==0, 0, ((2*s/x)^(nu)*gamma(nu+1)*psigamma(nu+1, 1)+(digamma(nu+1)+log(2*s/x))*((2*s/x)^(nu)*gamma(nu+1)*(digamma(nu+1)+log(2*s/x)))) )
  })
    
  out.01<-lapply(1:ncol(d), function(j) {
    j.p1<-p1[,j]
    j.p2<-p2
    j.p3<-p3[,j]
    j.p4<-p4[,j]
    j.p5<-p5
    j.p6<-p6[,j]
    j.p7<-p7[,j]
    j.p8<-p8[,j]
    j.p1*j.p2*(j.p3+j.p4*(j.p5+j.p6))+j.p7*j.p8
  })  
    
  out<-do.call(cbind, out.01)
  return(out)
}