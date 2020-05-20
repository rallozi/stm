d2s.bessel.mat<-function(nu, s, d) {
  k<-seq(0,100,1)
  p1<-apply(d, 2, function(x) {
    out.01<-lapply(1:length(x), function(i) {
      xx<-x[i]
      if(xx==0) {
        return(0)
      } else {
        xxx<-(((-1)^(k))/(factorial(k)*gamma(k+nu+1)))*((2*k+nu)*(2*k+nu+1)*(xx/2)^(2*k+nu-1)*(1/s)^(2*k+nu+2))
        if(any(is.nan(xxx) | xxx==Inf | xxx==-Inf)) {xxx[which(is.nan(xxx) | xxx==Inf | xxx==-Inf)]<-0}
        (xx/2)*sum(xxx)
      }
    })
    unlist(out.01)
  })
  return(p1)
}



