`smoothing` <-
function(ss) {
  m <- ss$m    #Kalman filter values of y for all time points
  C <- ss$C    #Kalman filter values of variance of y for all time points
  m0<- ss$m0   #initial values for y0
  C0<- ss$C0   #initial values for C0
  n <- ss$n    #number of spatial locations
  tp <- ss$tp  #number of time points

  #To get smoothed values, one runs the following backward recursion for t = T,
  #T − 1, . . . , 1, which is sometimes called the Kalman smoother:
  mu <- matrix(NA, tp, n)
  mu[tp,] <- t(ss$Fmat) %*% m[tp,] #Last time point
  
  for (tt in (tp-1):1)    { 
	if (ss$p == 1){
	  nextstep <-smootherstep.uni(m=m[tt,],
	                              C=C[[tt]],
	                              Gmatx=ss$Gmat,
	                              Wtx=ss$Wmat,
	                              mx=m[tt+1,],
	                              Cx=C[[tt+1]])
	}	else {
	  nextstep <-smootherstep(m=matrix(m[tt,],nrow=1),    #y_(t-1)^(t-1)
	                          C=C[[tt]],                  #P_(t-1)^(t-1)
	                          Gmatx=ss$Gmat,              #G
	                          Wtx=ss$Wmat,                #Sigma Eta
	                          mx=matrix(m[tt+1,],nrow=1), #y_(t)^(T)
	                          Cx=C[[tt+1]]) #P_(t)^(T)
	}

	m[tt,]  <- nextstep$ms #y_(t)^(T)
	C[[tt]] <- nextstep$Cs #P_(t)^(T)
	mu[tt,] <- t(ss$Fmat) %*% m[tt,]
	}
    
	ss$m0 <- m0
	ss$C0 <- C0
	if (is.ts(ss$Z))   ss$m <- ts(m,start(ss$z),end=end(ss$z),frequency=frequency(ss$z))
	else
	ss$m <- m
	ss$C <- C
	ss$mu <- mu
	ss
}

