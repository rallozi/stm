list.to.vector<-function(dat) {
  out<-lapply(names(dat), function(param) {
    idat<-dat[[param]]
    if(length(idat)>1) {
      idat.vec<-as.vector(idat)
      names(idat.vec)<-paste0(param, ".", 1:length(idat.vec))
    } else {
      idat.vec<-as.vector(idat)
      names(idat.vec)<-param
    }
    return(idat.vec)
  })
  return(unlist(out)) 
}