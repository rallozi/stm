`filtering` <-
function(ss) {
	m <- matrix(NA,ss$tp,ss$r) #tp*r matrix, mean of y at each time tp
	C <- vector("list",ss$tp) #covariance matrix of y at each time tp
	
	a <- matrix(NA,ss$tp,ss$r)
	R <- vector("list",ss$tp) #
	
	if(!isSymmetric(ss$Wmat)) {
	  ss$Wmat[lower.tri(ss$Wmat)] = t(ss$Wmat)[lower.tri(ss$Wmat)]
	}
	
	firststep <-
		filterstep(
			z = matrix(ss$z[1,]), #initial observation for all n locations
			Fmat	= ss$Fmat, #r*n matrix
			Gmat = ss$Gmat,  #r*r matrix
			Vt 	= ss$Vmat,   #n*n matrix
			Wt	= ss$Wmat,   #r*r matrix
			mx	= ss$m0, #starting value for mu, mean of y
			Cx	= ss$C0, #starting covariance matrix of y
			XXXcov  = ss$XXX[[1]],
			betacov = ss$beta,
			flag    = ss$flag.cov #indicator for covariates
               )
	m[1,] <- firststep$m #mu_(t=1)
	C[[1]]<- firststep$C #C_(t=1)
	
	a[1,] <- firststep$a #
	R[[1]]<- firststep$C #R_(t=1)
	
	
	loglik<- firststep$loglikterm
	
  ## run the recursion until time tp
	for (tt in 2:ss$tp) {
	nextstep <-
		filterstep(
			z= matrix(ss$z[tt,]),
			Fmat = ss$Fmat,
			Gmat = ss$Gmat,
			Vt = ss$Vmat,
			Wt = ss$Wmat,
			mx 	= matrix(m[tt-1,],nrow=1),
			Cx 	= C[[tt-1]],
			XXXcov  = ss$XXX[[tt]],
			betacov = ss$beta,
			flag    = ss$flag.cov
		)
      m[tt,]  <- nextstep$m
      C[[tt]] <- nextstep$C
      
      a[tt,]  <- nextstep$a
      R[[tt]] <- nextstep$R
      loglik  <- loglik + nextstep$loglikterm #summed over all timepoints
   }
   
	if (is.ts(ss$z))
		ss$m <- ts(m,start(ss$z),end=end(ss$z),frequency=frequency(ss$z))
	else
	ss$m <- m
	ss$C <- C
	
	ss$a <- a
	ss$R <- R
	
	ss$loglik <- loglik
	ss
}

