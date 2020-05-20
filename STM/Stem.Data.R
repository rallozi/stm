`Stem.Data` <-
function(...) {
    if (nargs() == 1)
        x <- as.list(...)
    else
        x <- list(...)

# Stem.Data components
  comp <- c("z", "covariates", "funcdist")

#Verify if all the required components are given by user	
  compInd <- match(comp, names(x))
    if (any(is.na(compInd)))
        stop(paste("Component(s)", paste(comp[is.na(compInd)], collapse=", "),
                   "is (are) missing"))
                   
  #Check that values are numeric           
  check.Z<-lapply(x$z, function(i) {
    !apply(i,1,is.numeric)
  })
  if (sum(unlist(check.Z))>0)  stop("Component z must be numeric.")
  
  check.covariates<-lapply(x$covariates, function(i) {
    !apply(i,1,is.numeric)
  })  
  if (sum(unlist(check.covariates))>0)   stop("Component covariates must be numeric.")
  
  check.dist<-!apply(x$funcdist, 1, is.numeric)
  if (sum(unlist(check.dist))>0 || any(check.dist<0))  stop("Component funcdist must be numeric and non-negative.")
  
  
  #Controls over the dimensions of the components  
  get.p<-unlist(lapply(x$covariates, function(i){
    ncol(i)
  }))
  p<-unique(get.p)
  if(length(p)!=1) stop("Every individual must have the same number of covariates.")
  
  get.tp<-unlist(lapply(x$z, function(i){
    nrow(i)
  }))
  tp<-unique(get.tp)
  if(length(tp)!=1) stop("Every individual must have the same number of time points.")
  
  get.n<-unlist(lapply(x$z, function(i){
    ncol(i)
  }))
  n<-unique(get.n)
  if(length(n)!=1) stop("Every individual must have the same number of spatial locations.")
  
  check.cov<-lapply(x$covariates, function(i) {
    nrow(i)
  })
  if(any(unlist(check.cov)!=n)) stop("The number of rows in each covariate matrix must be n") 
  
  #check for missing values
  if(any(unlist(lapply(x$z, is.na))) || any(unlist(lapply(x$covariates, is.na))))
        stop("Missing values are not allowed in components z and covariates")

  out<-list(z=x$z,
            covariates=x$covariates,
            funcdist=x$funcdist,
            p=p,   #dimension of latent process
            tp=tp, #number of timepoints
            n=n)   #number of spatial locations

  class(out) <- "Stem.Data"
  return(out)  
}

