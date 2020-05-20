##This version does not iterate through NR for each EM iteration
##This version only updates gamma (b) and keeps s, nu constant

`kalman.s` <-
  function (z, dist, r, tp, n, p, phi_j, covariates, SigmaEta, dQds_j, learning.rate) {
    
    zz=lapply(z, function(df) {
      ts(df) #The function ts is used to create time-series objects.
    })
    m<-length(zz)
    
    ####################################
    ###Kalman filtering and smoothing###
    ####################################
    ids<-names(zz)
    
    if(SigmaEta=="Unstructured") {
      Wmat_j=phi_j$Sigmaeta #Sigma Eta
    } else if(SigmaEta=="Diagonal") {
      Wmat_j=phi_j$sigma2eta*diag(r)
    }
    
    kalman.res<-lapply(ids, function(i) {
      SSmodel  = list(z	= zz[[i]], #individual i's time series data
                      Fmat 	= t(phi_j$K[[i]]), #transpose of individual i's loading matrix
                      Gmat 	= phi_j$G[[i]],    #individual i's transition matrix
                      Vmat = phi_j$sigma2omega * Bessel.Cov(n=n, b=phi_j$b, d=dist, nu=phi_j$nu, s=phi_j$s),
                      Wmat 	= Wmat_j,
                      m0   	= t(phi_j$m0),    #mean of y0
                      C0   	= phi_j$C0,       #covariance matrix of y0
                      beta 	= phi_j$beta,     #beta values
                      phi  	= phi_j, 
                      XXX  	= covariates[[i]], #individual i's covariates matrix
                      flag.cov 	= TRUE,
                      n	= n,
                      p	= p,
                      r	= r,
                      tp=tp,
                      m = NA,
                      C	= NA,
                      loglik = NA)
      
      filter.m   	= filtering(SSmodel)
      smoother.m 	= smoothing(filter.m)
      
      ####calculating  B0 and P_0_T using kalman filter, smoother and initial values
      R1    = SSmodel$Gmat %*% SSmodel$C0 %*% t(SSmodel$Gmat) + SSmodel$Wmat  #P_1^0
      B0    = SSmodel$C0 %*% SSmodel$Gmat %*% 	solve(R1)
      P_0_T = SSmodel$C0 + B0 %*% (smoother.m$C[[1]] - R1 ) %*% t(B0)
      
      ########################
      ###Define the elements for B_function (output of mod1$filter)
      ###Use B_function (see functions.R) ---> list of B_t (t=1,...,n)
      ########################
      m    	= filter.m$m
      C    	= filter.m$C
      B    	= list()
      
      for (tt in (tp-1):1) {            
        nextstep = B_function(m     = matrix(m[tt,],nrow=1),
                              C     = C[[tt]],
                              Gmatx = SSmodel$Gmat,
                              Wtx   = SSmodel$Wmat,
                              mx    = matrix(m[tt+1,],nrow=1),
                              Cx    = C[[tt+1]])
        
        m[tt,]  = nextstep$ms #equal to filter.m$m
        C[[tt]] = nextstep$Cs #equal to filter.m$C
        B[[tt]] = nextstep$B  #what I need
      }
      
      #B_n=C_n * t(Gmat) * solve(R_n+1)
      #where  R_n+1=Gmat * C_n * t(Gmat) +Wmat
      
      B[[tp]] = filter.m$C[[tp]] %*% t(SSmodel$Gmat) %*%
        solve(SSmodel$Gmat %*% filter.m$C[[tp]] %*%
                t(SSmodel$Gmat)+SSmodel$Wmat)
      
      ###############################################
      ###Recursion for LAG ONE COVARIANCE SMOOTHER
      ###############################################
      CCC = list() #empty list for lag one covariance values
      
      ###FIRST STEP: Calculate C*_{n,n-1} (start of the iterative procedure --> CCC[[n]])
      ###Pn_n-1=R_n=G_n%*%C_n-1%*%t(G_n)+W_n
      ###A_n=R_n%*%F_n%*%solve(t(F_n)%*%R_n%*%F_n+V_n)
      R_n=SSmodel$Gmat %*% filter.m$C[[tp-1]] %*%t(SSmodel$Gmat)+ SSmodel$Wmat
      A_n=R_n%*%SSmodel$Fmat %*% solve(t(SSmodel$Fmat) %*% R_n %*% SSmodel$Fmat + SSmodel$Vmat)
      CCC[[tp]]=(diag(r)-A_n %*% t(SSmodel$Fmat)) %*% SSmodel$Gmat %*% filter.m$C[[tp-1]]
      
      ################################
      ###loop for calculating the other elements of CCC using cov_lagone function
      ################################
      for (tt in (tp):2) {
        passiCCC =cov_lagone(C_t_minus_1=filter.m$C[[tt-1]],
                             B_t_minus_1=B[[tt-1]],
                             Gmat=SSmodel$Gmat,
                             if (tt>=3)  
                               B_t_minus_2= B[[tt-2]] else B_t_minus_2 = B0,
                             CCCx       		 = CCC[[tt]])
        
        CCC[[tt-1]] = passiCCC$cov
      }
      
      return(list(smoother=smoother.m, filter=filter.m, SSmodel=SSmodel, B0=B0, P_0_T=P_0_T, CCC=CCC, loglik=filter.m$loglik))
    })
    names(kalman.res)<-ids
    
    loglik.01<-lapply(ids, function(i) {
      kalman.res[[i]]$loglik
    })
    loglik_jj<-sum(unlist(loglik.01))
    
    #######################################################
    ################EM Parameters estimates################
    #######################################################
    ################################
    ###PARAMETER 1: m0=y_0_T
    ################################
    y_0_T<-lapply(kalman.res, function(x) {
      t(x$SSmodel$m0) + x$B0 %*% (t(matrix(x$smoother$m[1,],nrow=1))-x$SSmodel$Gmat %*% t(x$SSmodel$m0))
    })
    m0_jj<-apply(do.call(cbind, y_0_T), 1, mean)
    
    
    ################################
    ###PARAMETER 2: G
    ################################
    ###S00 by time for every individual i
    S00.i<-lapply(kalman.res, function(x) {
      S00list = list()
      #first one uses yi0, Pi0
      y0_j<-t(x$SSmodel$m0) + x$B0 %*% (t(matrix(x$smoother$m[1,],nrow=1))-x$SSmodel$Gmat %*% t(x$SSmodel$m0))
      S00list[[1]] = y0_j %*% t(y0_j) + x$P_0_T
      
      for (tt in 2:tp) {
        S00list[[tt]] = matrix(x$smoother$m[(tt-1),],nrow=r) %*% t(matrix(x$smoother$m[(tt-1),],nrow=r)) +  x$smoother$C[[(tt-1)]]
      }
      return(S00list)
    })
    
    S00.m<-lapply(S00.i, function(i) {
      sumMatrices(i)
    })
    
    S00<-sumMatrices(S00.m)
    
    ###S10 by time for every individual i
    S10.i<-lapply(kalman.res, function(x) {
      S10list = list()
      #first uses time 0
      y0_j<-t(x$SSmodel$m0) + x$B0 %*% (t(matrix(x$smoother$m[1,],nrow=1))-x$SSmodel$Gmat %*% t(x$SSmodel$m0))
      S10list[[1]] = matrix(x$smoother$m[1,],nrow=r) %*% t(y0_j) + x$CCC[[1]]
      
      for (tt in 2:tp) {
        S10list[[tt]] = matrix(x$smoother$m[tt,],nrow=r) %*% t(matrix(x$smoother$m[(tt-1),],nrow=r)) +  x$CCC[[tt]]
      }
      return(S10list)
    })
    
    S10.m<-lapply(S10.i, function(i) {
      sumMatrices(i)
    })
    
    S10<-sumMatrices(S10.m)
    
    G_jj<- S10 %*% solve(S00)
    
    
    ################################
    ###PARAMETER 3: Sigmaeta
    ################################
    ###S11 by time for every individual i
    S11.i<-lapply(kalman.res, function(x) {
      S11list = list()
      for (tt in 1:tp) {
        S11list[[tt]] = matrix(x$smoother$m[tt,],nrow=r) %*% t(matrix(x$smoother$m[tt,],nrow=r)) +  x$smoother$C[[tt]]
      }
      return(S11list)
    })
    
    E.01<-lapply(ids, function(i) {
      S00<-sumMatrices(S00.i[[i]])
      S11<-sumMatrices(S11.i[[i]])
      S10<-sumMatrices(S10.i[[i]])
      #G<-kalman.res[[i]]$SSmodel$Gmat
      G<-G_jj
      E.i<-S11 - S10 %*% t(G) - G %*% t(S10) + G %*% S00 %*% t(G)
      return(E.i)
    })
    E<-sumMatrices(E.01)
    
    if(SigmaEta=="Unstructured") {
      Sigmaeta_jj<-E/(m*tp)
    } else if(SigmaEta=="Diagonal") {
      SigmaetaInv<-solve(diag(r))
      sigma2eta_jj<-sum(diag(SigmaetaInv %*% E)) /(m*r*tp)
    }
    
    #############################
    ###PARAMETER 3: sigma2omega
    #############################
    W.i<-lapply(ids, function(i){
      z.it<-lapply(1:tp, function(tt) {
        matrix(as.numeric(z[[i]][tt,]), nrow=n)
      })
      X.i<-covariates[[i]]
      K.i<-phi_j$K[[i]]
      
      y.it<-lapply(1:tp, function(tt) {
        matrix(kalman.res[[i]]$smoother$m[tt,], nrow=r)
      })
      
      P.it<-lapply(1:tp, function(tt) {
        kalman.res[[i]]$smoother$C[[tt]]
      })
      
      W.it<-lapply(1:tp, function(tt) {
        step.01<-z.it[[tt]] - X.i[[tt]] %*% phi_j$beta - K.i  %*% y.it[[tt]]
        step.02<- step.01 %*% t(step.01)
        W.it<-step.02+ K.i %*% P.it[[tt]] %*% t(K.i)
        return(W.it)
      })
      sumMatrices(W.it)
    })
    W<-sumMatrices(W.i)
    
    GammaInv<-solve(Bessel.Cov(n=n, b=phi_j$b, d=dist, nu=phi_j$nu, s=phi_j$s))
    sigma2omega_jj_num<-sum(diag(GammaInv %*% W))
    sigma2omega_jj<-sigma2omega_jj_num/(m*n*tp)
    
    
    ####################################
    ###PARAMETER 4: beta coefficients###
    ####################################
    SigmaXiInv<-solve(sigma2omega_jj * Bessel.Cov(n=n, b=phi_j$b, d=dist, nu=phi_j$nu, s=phi_j$s))
    MM.i<-lapply(ids, function(i) {
      X.i<-covariates[[i]]
      MM.it<-lapply(X.i, function(Xit) {
        t(Xit) %*% SigmaXiInv %*% Xit
      })
      sumMatrices(MM.it)
    })
    MM<-sumMatrices(MM.i)
    # if(m==1) {
    #   MM<-MM.i[[1]]
    # } else {
    #   MM<-sumMatrices(MM.i)
    # }
    
    
    v.i<-lapply(ids, function(i) {
      z.it<-lapply(1:tp, function(tt) {
        matrix(as.numeric(z[[i]][tt,]), nrow=n)
      })
      K.i<-phi_j$K[[i]]
      X.i<-covariates[[i]]
      
      y.it<-lapply(1:tp, function(tt) {
        matrix(kalman.res[[i]]$smoother$m[tt,], nrow=r)
      })
      
      v.it<-lapply(1:tp, function(tt) {
        step.01<-t(X.i[[tt]]) %*% SigmaXiInv
        step.02<-z.it[[tt]] - K.i %*% y.it[[tt]]
        step.01 %*% step.02
      })
      v.i<-sumMatrices(v.it)
      return(v.i)
    })
    v<-sumMatrices(v.i)
    # if(m==1) {
    #   v<-v.i[[1]]
    # } else {
    #   v<-sumMatrices(v.i)
    # }
    beta_jj = solve(MM) %*% v
    
    #if(det(MM) != 0) {beta_jj = solve(MM) %*% v}
    #if(det(MM)  < 10^(-7)) {cat("Error in beta estimation! The matrix can not be inverted!!!!")}
    
    #######################################################
    ################NR Parameters estimates################
    #######################################################
    dQds_jj<-dQds(m=m, tp=tp, GammaInv=GammaInv, sigma2omega=sigma2omega_jj, W=W, d=dist, s=phi_j$s, nu=phi_j$nu)
    d2Qd2s_jj<-d2Qd2s(m=m, tp=tp, GammaInv=GammaInv, sigma2omega=sigma2omega_jj, W=W, d=dist, s=phi_j$s, nu=phi_j$nu)
    
    s_jj<- phi_j$s -learning.rate*(1/(d2Qd2s_jj))*dQds_jj
    NR_jj<- c(phi_j$b, phi_j$nu, s_jj)
    
    new.G<-lapply(1:m, function(i) {
      G_jj
    })
    names(new.G)<-ids
    
    ###################################################
    ###New parameter vector (K, G, C0 don't change)
    ###################################################
    if(SigmaEta=="Unstructured") {
      ##############
      ##SETTING 2###
      ##############
      phi_jj = list(K=phi_j$K, #does not change
                    sigma2omega=sigma2omega_jj,
                    b=NR_jj[1],
                    nu=NR_jj[2],
                    s=NR_jj[3],
                    beta=beta_jj,
                    G=new.G,
                    Sigmaeta=Sigmaeta_jj,
                    m0=m0_jj,
                    C0=phi_j$C0) #does not change
      
    } else if(SigmaEta=="Diagonal") {
      phi_jj = list(K=phi_j$K, #does not change
                    sigma2omega=sigma2omega_jj,
                    b=NR_jj[1],
                    nu=NR_jj[2],
                    s=NR_jj[3],
                    beta=beta_jj,
                    G=new.G,
                    sigma2eta=sigma2eta_jj,
                    m0=m0_jj,
                    C0=phi_j$C0) #does not change
      
    } 
    
    return(list(phi=phi_jj, likelihood=loglik_jj, dQds_j=dQds_jj))
  }

