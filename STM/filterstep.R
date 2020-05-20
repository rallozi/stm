`filterstep` <-
function(z,Fmat,Gmat,Vt,Wt,mx,Cx,XXXcov,betacov,flag)
  {
    #F=K'; G=G; V=SigmaXi; W=SigmaEta
    a <- Gmat %*% t(mx)                   #y_(t)=Gy_(t-1)  (**Xu Step 1)
    R <- Gmat %*% Cx %*% t(Gmat) + Wt     #R=P_(t)^(t-1)=var(y_(t))=GC_(t-1)G'+W  (**Xu Step 2)

    #difference between observation z and predicted value 
    f <- t(Fmat) %*% a                    #Ka=Ky_(t)=KGy_(t-1)
    if(flag==FALSE)      e <- z - f       #z-Ky_(t)^(t-1)
    else                 e <- z - (XXXcov %*% betacov) - f #z-XB-Ky_(t)^(t-1)

    Q <- t(Fmat) %*% R %*% Fmat + Vt      #Q=var(z)=KRK'+V
    A <- R %*% Fmat %*% solve(Q)          #(**Xu Step 3)
    m <- a + A%*%e                        #y_(t)^(t)=y_(t)^(t-1)+A_(t)*(z_(t)-X_(t)B-Ky_(t)^(t-1)) #(**Xu Step 4)
    C <- R-A%*%t(Fmat)%*%R                #new C #(**Xu Step 5)

    if(!isSymmetric(Q)) {
      Q[lower.tri(Q)] = t(Q)[lower.tri(Q)]
    }
    
    if(dmvnorm(as.numeric(z),as.numeric(f+(XXXcov %*% betacov)),Q)==0) {
      loglikterm <-log(1e-323)
    } else {
      loglikterm <-log(dmvnorm(as.numeric(z),as.numeric(f+(XXXcov %*% betacov)),Q)) #incomplete data likelihood
    }
    #if (length(z)>1  &&  flag==FALSE)  loglikterm <-log(dmvnorm(as.numeric(z),as.numeric(f),Q))
    #if (length(z)==1 &&  flag==TRUE)  loglikterm <- -0.5*( log(2*pi) + log(Q) + (z - (f-XXXcov %*% betacov)) ^2/Q)
    #if (length(z)==1 &&  flag==FALSE)  loglikterm <- -0.5*( log(2*pi) + log(Q) + (z-f)^2/Q)

   list(m=m,C=C,R=R,a=a,loglikterm=loglikterm)
   #list(m=m,C=C)
  }

