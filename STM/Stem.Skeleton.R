`Stem.Skeleton` <-
function(...) {
    if (nargs() == 1)
        x <- as.list(...)
    else
        x <- list(...)

# Stem.Skeleton components
   comp <- c("phi", "r", "K")

#Verify if all the required components are given by user 
	compInd <- match(comp, names(x))
   if (any(is.na(compInd)))
        stop(paste("Component(s)", paste(comp[is.na(compInd)], collapse=", "),
                   "is (are) missing"))
                   
#Controls over the dimensions of the components 
  if (sum(!unlist((lapply(x$phi,is.numeric))))>0)  stop("Component of phi must be numeric and positive")
	if(x$SigmaEta=="Unstructured") {
	  if(!(nrow(x$phi$Sigmaeta) == x$r && ncol(x$phi$Sigmaeta) == x$r)) stop("The dimension of matrix Sigmaeta must be r*r")
	} else if(x$SigmaEta=="AR(1)") {
	  if(!(length(x$phi$sigma2eta) == 1 && (x$phi$sigma2eta) > 0)) stop("Sigma2eta must be scalar and positive")
	  if(r>1) {
	    if(!(length(x$phi$rho) == 1 && (x$phi$rho) <= 1 || (x$phi$rho) >= -1)) stop("rho must be scalar and between -1 and 1")
	  }
	}
  if(!(nrow(x$phi$m0) == x$r && ncol(x$phi$m0) == 1)) stop("The dimension of matrix mu0 must be r*1")
  if(!(nrow(x$phi$C0) == x$r && ncol(x$phi$C0) == x$r)) stop("The dimension of matrix C0 must be r*r")
	if(!(nrow(x$phi$G) == x$r && ncol(x$phi$G) == x$r)) stop("The dimension of matrix G must be r*r")
	if(!(length(x$phi$s) == 1 && (x$phi$s) > 0)) stop("Scale s must be scalar and positive")
  if(!(length(x$phi$sigma2eps) == 1 && (x$phi$sigma2eps) > 0)) stop("Sigma2eps must be scalar and positive")
  if(!(length(x$phi$sigma2omega) == 1 && (x$phi$sigma2omega) > 0)) stop("Sigma2omega must be scalar and positive")
	
#Definition of the class Stem.Skeleton
 out<-list(phi=x$phi,
           r=x$r,
           K=x$K)
 class(out) <- "Stem.Skeleton"
 return(out)
}

