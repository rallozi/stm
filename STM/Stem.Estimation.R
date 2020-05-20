##############################################################################################################
#Utilizes kalman.b, kalman.nu, kalman.s, and kalman.v2:
#1. Maximizes b, nu, s, and then uses the new values as starting values for the EM algorithm
#2. NR_jj<- est_j - learning.rate * solve(hess_jj) %*% grad_jj
#3. Convergence based on relative absolute change in likelihood
##############################################################################################################

`Stem.Estimation` <-
  function(StemModel, precision, max.iter=1e4, SigmaEta, learning.rate){
    
    z=StemModel$data$z
    m=length(StemModel$data$z)
    covariates=StemModel$data$covariates
    dist=StemModel$data$funcdist
    
    p=StemModel$data$p      #p=number of fixed-effects parameters in model
    n=StemModel$data$n      #n=number of spatial locations
    tp=StemModel$data$tp    #tp=number of time points
    r=StemModel$skeleton$r  #r=dimension of latent process
    
    phi_start=StemModel$skeleton$phi
    
    #Add K and G to phi_start since they will be used in the EM algorithm
    phi_start$K=StemModel$skeleton$K 
    phi_start$G=lapply(1:m, function(i){
      StemModel$skeleton$phi$G
    })
    names(phi_start$G)<-names(z)
    
    #Transformed parameter - Setting 2
    phi_start$b=phi_start$sigma2eps/phi_start$sigma2omega
    
    #Duplicate X for each i tp times
    temp<-covariates
    ids<-names(covariates)
    covariates<-lapply(ids, function(i) {
      covariates.i<-covariates[[i]]
      covariates.it<-lapply(1:tp, function(tt) {
        covariates.i
      })
      return(covariates.it)
    })
    names(covariates)<-ids
    
    #phi0=phi_start
    
    #############################
    ###EM algorithm while loop
    #############################
    converged_EM<-FALSE
    n_iter_EM<-1
    
    allestimates<-list()
    alllikelihood<-list()
    
    my.thresh<-2e3
    
    dQdgamma_j<-1e10
    #maximize b keeping s, nu constant
    while(abs(dQdgamma_j)>my.thresh) {
      EM.j = kalman.b(z=z,
                      p=p,
                      n=n,
                      tp=tp,
                      r=r,
                      dist=dist,
                      phi_j=phi_start,
                      covariates=covariates,
                      SigmaEta=SigmaEta,
                      dQdgamma_j=dQdgamma_j,
                      learning.rate=1)
      
      
      ###Updating the iteration number and the parameter vector
      phi_start 	= EM.j$phi
      dQdgamma_j  = EM.j$dQdgamma_j 
      cat(paste("****************Optimization for b, gradient=",round(dQdgamma_j, 1), sep=""),"\n")
      print(EM.j$phi$b)
    }
    
    phi_start<-EM.j$phi
    dQds_j<-1e10
    while(abs(dQds_j)>my.thresh) {
      EM.j = kalman.s(z=z,
                      p=p,
                      n=n,
                      tp=tp,
                      r=r,
                      dist=dist,
                      phi_j=phi_start,
                      covariates=covariates,
                      SigmaEta=SigmaEta,
                      dQds_j=dQds_j,
                      learning.rate=1)
      
      
      ###Updating the iteration number and the parameter vector
      phi_start 	= EM.j$phi
      dQds_j  = EM.j$dQds_j 
      cat(paste("****************Optimization for s, gradient=",round(dQds_j, 1), sep=""),"\n")
      print(EM.j$phi$s)
    }
    
    phi_start<-EM.j$phi
    dQdnu_j<-1e10
    while(abs(dQdnu_j)>my.thresh) {
      EM.j = kalman.nu(z=z,
                       p=p,
                       n=n,
                       tp=tp,
                       r=r,
                       dist=dist,
                       phi_j=phi_start,
                       covariates=covariates,
                       SigmaEta=SigmaEta,
                       dQdnu_j=dQdnu_j,
                       learning.rate=0.5)
      
      
      ###Updating the iteration number and the parameter vector
      phi_start 	= EM.j$phi
      dQdnu_j  = EM.j$dQdnu_j 
      cat(paste("****************Optimization for nu, gradient=",round(dQdnu_j, 1), sep=""),"\n")
      print(EM.j$phi$nu)
    }
    
    StemModel$skeleton$phi.upd<-EM.j$phi
    phi_start<-StemModel$skeleton$phi.upd
    
    #maximize all with new starting values
    while ((!converged_EM)  && n_iter_EM <= max.iter){
      start.time <- Sys.time()
      cat(paste("****************EM Algorithm - iteration n.",n_iter_EM, sep=""),"\n")
      EM.j = kalman(z=z,
                    p=p,
                    n=n,
                    tp=tp,
                    r=r,
                    dist=dist,
                    phi_j=phi_start,
                    covariates=covariates,
                    SigmaEta=SigmaEta,
                    learning.rate=learning.rate)
      
      print(EM.j$likelihood)
      if(n_iter_EM>1) {
        #Add likelihood
        alllikelihood[[n_iter_EM-1]]<-EM.j$likelihood #first likelihood is for starting values; second likelihood is for iteration 1
      }
      
      if(n_iter_EM>2) { #compare likelihoods after 2 iterations
        converged_EM<-(abs(EM.j$likelihood-alllikelihood[[n_iter_EM-2]])/abs(EM.j$likelihood))<precision
        #converged_EM<-sum((allestimates[[n_iter_EM-1]]-allestimates[[n_iter_EM-2]])^2)<precision
      }
      
      ###Updating the iteration number and the parameter vector
      phi_start 	= EM.j$phi
      n_iter_EM 	= n_iter_EM + 1

      check.list.01<-c("m0", "sigma2omega", "beta", "nu", "s", "G")
      
      ##############
      ##SETTING 2###
      ##############
      if(SigmaEta=="Unstructured") {
        check.list<-c(check.list.01, "Sigmaeta", "b")
      } else if(SigmaEta=="Diagonal") {
        check.list<-c(check.list.01, "sigma2eta", "b")
      } 
      
      EM.j$phi[["G"]]<-EM.j$phi[["G"]][[1]] #all the same
      print(list.to.vector(EM.j$phi[check.list]))
      
      allestimates[[n_iter_EM-1]]<-list.to.vector(EM.j$phi[check.list]) #first iteration is first allestimates 
      
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      print(time.taken)
    } # here the while loop ends
    
    #get updated kalman filter - just need EM.j$all.a from this
    EM.j = kalman.v2(z=z,
                     p=p,
                     n=n,
                     tp=tp,
                     r=r,
                     dist=dist,
                     phi_j=phi_start,
                     covariates=covariates,
                     SigmaEta=SigmaEta,
                     learning.rate=learning.rate)
    StemModel$estimates$yt.tmin1=EM.j$all.a
    StemModel$estimates$Pt.tmin1=EM.j$all.R
    
    #Transform back to sigma2eps
    phi_start$sigma2eps<-phi_start$b*phi_start$sigma2omega
    
    #Remove anything not estimated
    phi_start$G<-phi_start$G[[1]]
    phi.estimated=phi_start[-which(names(phi_start) %in% c("K", "b", "C0"))]
    StemModel$estimates$phi.hat=phi.estimated 
    
    #Add estimates for all iterations
    StemModel$estimates$allestimates=allestimates
    
    #Add likelihoods for all iterations
    StemModel$estimates$loglik=unlist(alllikelihood)
    
    convergence.par=list(loglik.conv=converged_EM,
                         iterEM=n_iter_EM-1)
    StemModel$estimates$convergence.par = convergence.par
    
    return (StemModel)		
  }

