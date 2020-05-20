`Stem.Model` <-
function(...) {
	
     if (nargs() == 1)
        x <- as.list(...)
    else
        x <- list(...)

    #Stem.Model components : skeleton and data
    skeleton <- Stem.Skeleton(x) #(phi, r, K, G)
    data <- Stem.Data(x) #(z, covariates, funcdist, p, tp, n)
    
    if(length(skeleton$phi$beta) != data$p) stop("The length of beta must be equal to the number of columns of covariates")
    
    check.K<-lapply(skeleton$K, function(i) {
      K.nrow<-nrow(i)
      K.ncol<-ncol(i)
      return(list(K.nrow=K.nrow, K.ncol=K.ncol))
    })
    K.nrow<-sapply(check.K, "[[", 1)
    K.ncol<-sapply(check.K, "[[", 2)
    if(!(all(K.nrow == data$n && K.ncol == skeleton$r))) stop("The dimension of every matrix K must be n*r")
    
    x=list(skeleton=skeleton,data=data)
    class(x) <- "Stem.Model"
    return(x)	
}

