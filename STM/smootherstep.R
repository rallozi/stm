`smootherstep` <-
function(m,C,Gmatx,Wtx,mx,Cx) {
    Rx <- Gmatx %*% C %*% t(Gmatx) + Wtx #P_(t)^(t-1)
    B  <- C %*% t(Gmatx) %*% solve( Rx ) #(**Xu Step 1)
    ms <- t(m) + B%*%(t(mx) - Gmatx %*%t(m)) #(**Xu Step 2)
    Cs <- C + B%*%(Cx-Rx)%*%t(B) #(**Xu Step 3)
    list(ms=ms,Cs=Cs)
  }

