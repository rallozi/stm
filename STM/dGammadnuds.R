bessel<-function(nu, s, d) {
  m<-seq(0,100,1)
  sum((((-1)^m)/(factorial(m)*gamma(m+nu+1)))*(d/(2*s))^(2*m+nu))
}

ds.bessel<-function(nu, s, d) {
  m<-seq(0,100,1)
  -d/(2*s^2)*sum((((-1)^m)/(factorial(m)*gamma(m+nu+1)))*(2*m+nu) *(d/(2*s))^(2*m+nu-1))
}

dn.bessel<-function(nu, s, d) {
  m<-seq(0,100,1)
  bessel(nu, s, d)*log(d/(2*s))-sum((((-1)^m)/(factorial(m)*gamma(m+nu+1)))*digamma(m+nu+1)*(d/(2*s))^(2*m+nu))
}

dnds.bessel<-function(nu, s, d) {
  m<-seq(0,100,1)
  log(d/(2*s))*ds.bessel(nu, s, d)+(1/s)*bessel(nu, s, d)+sum((((-1)^m)/(factorial(m)*gamma(m+nu+1)))*digamma(m+nu+1)*((d/2)^(2*m+nu))*(2*m+nu)*(1/s)^(2*m+nu+1))
}

dsdn.bessel<-function(nu, s, d) {
  m<-seq(0,100,1)
  d*(dn.bessel(nu+1,s,d)-dn.bessel(nu-1,s,d))/(2*s^2)
}

#mine
dGammadnuds<-function(nu, s, d) {
  p1<-gamma(nu+1)

  p2<-apply(d, 2, function(x) {
    out<-rep(NA, length(x))
    for(i in 1:length(x)) {
      xx<-x[i]
      if(xx==0) {
        out[i]<-0
      } else {
        xxx.1<-(2*s/xx)^(nu)*dsdn.bessel(nu, s, d=xx)
        if(any(is.nan(xxx.1) | xxx.1==Inf | xxx.1==-Inf)) {
          xxx.1[which(is.nan(xxx.1) | xxx.1==Inf | xxx.1==-Inf)]<-0
        }
        
        xxx.2<-nu*s^(nu-1)*(2/xx)^(nu)*dn.bessel(nu, s, d=xx)
        if(any(is.nan(xxx.2) | xxx.2==Inf | xxx.2==-Inf)) {
          xxx.2[which(is.nan(xxx.2) | xxx.2==Inf | xxx.2==-Inf)]<-0
        }    
        
        out[i]<-xxx.1+xxx.2
      }
      
    }
    return(out)
  })

  p3<-apply(d, 2, function(x) {
    out<-rep(NA, length(x))
    for(i in 1:length(x)) {
      xx<-x[i]
      if(xx==0) {
        out[i]<-0
      } else {
        xxx.1<-(2*s/xx)^(nu)*ds.bessel(nu, s, d=xx)
        if(any(is.nan(xxx.1) | xxx.1==Inf | xxx.1==-Inf)) {
          xxx.1[which(is.nan(xxx.1) | xxx.1==Inf | xxx.1==-Inf)]<-0
        }
        xxx.2<-nu*s^(nu-1)*(2/xx)^(nu)*bessel(nu, s, d=xx)
        if(any(is.nan(xxx.2) | xxx.2==Inf | xxx.2==-Inf)) {
          xxx.2[which(is.nan(xxx.2) | xxx.2==Inf | xxx.2==-Inf)]<-0
        }    
        out[i]<-xxx.1+xxx.2
      }
    }
    return(out)
  })

  p4<-apply(d, 2, function(x) {
    out<-rep(NA, length(x))
    for(i in 1:length(x)) {
      xx<-x[i]
      
      if(xx==0) {
        out[i]<-0
      } else {
        xxx.1<-(1/s)*(2*s/xx)^(nu)*bessel(nu, s, d=xx)
        if(any(is.nan(xxx.1) | xxx.1==Inf | xxx.1==-Inf)) {
          xxx.1[which(is.nan(xxx.1) | xxx.1==Inf | xxx.1==-Inf)]<-0
        }
        xxx.2<-log(2*s/xx)*((2*s/xx)^(nu)*ds.bessel(nu, s, d=xx)+nu*s^(nu-1)*(2/xx)^(nu)*bessel(nu, s, d=xx))
        if(any(is.nan(xxx.2) | xxx.2==Inf | xxx.2==-Inf)) {
          xxx.2[which(is.nan(xxx.2) | xxx.2==Inf | xxx.2==-Inf)]<-0
        }    
        out[i]<-xxx.1+xxx.2
      }
    }
    return(out)
  })

  out.01<-lapply(1:ncol(d), function(j) {
    j.p1<-p1
    j.p2<-p2[,j]
    j.p3<-p3[,j]
    j.p4<-p4[,j]
    j.p1*(j.p2+digamma(nu+1)*j.p3+j.p4)
  })
  out<-do.call(cbind,out.01)
  return(out)

}


