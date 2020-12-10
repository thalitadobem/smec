##################################################################################
##        Semiparametric Mixed-Effects Models with Censored Responses           ##
## ---------------------------------------------------------------------------- ## 
## Article: A semiparametric mixed-effects model for censored longitudinal data ##
## Authors: Thalita B. Mattos, Larissa A. Matos and Victor H. Lachos            ##
## ---------------------------------------------------------------------------- ##
## Last update: 2020-09-15                                                      ##
##################################################################################

# ------------------------------------ #
# Function for estimating a SMEC model #
# ------------------------------------ #

nsmec <- function(y,  # vector of responses (Nx1)
                  cc, # vector of censoring (Nx1) - 0 if non-censored, 1 if censored.
                  x,  # design matrix of fixed effects (Nxp)
                  z,  # design matrix of random effects (N x q)
                  tt, # vector of ordered distinct values of the time points (r x 1)
                  ttc, # vector of the time points (N x 1)
                  nj,  # vector with size of clusters (m x 1)
                  LL,  # vector of lower limits (N x 1)
                  LU,  # vector of upper limits (N x 1)
                  Ns,  # incidence matrix
                  initial, # list with initial values (beta,sigma2, alphas,phi1, phi2, lambda)
                  struc="unc", # correlation structures - "dec" (damped exponential correlation), "ar" (continuous-time autoregressive of order 1), "sym" (compound symmetric), "unc" (uncorrelated)
                  lambda.fixed=TRUE, #logical, if TRUE, smoothing parameter (lambda) will not be estimated.
                  iter.max=200, # maximum number of iterations of the EM algorithm
                  precision=1e-5 # tolerance for the convergence criterion
                  ){

  if(is.null(x)) 
  {
    p <- 1
    x <- matrix(0,nrow = sum(nj),ncol=p)
    beta1 <- 0
  }else{
    p <- dim(x)[2]
    beta1 <- matrix(initial$beta,p,1)
  }
  
  m <- length(nj)[1]
  N <- sum(nj)
  q1 <- dim(z)[2]
  m2 <- m*q1

  # non-parametric
  rr <- length(tt)
  
  t2 <- natural.cubic.spline(tt, nknots=length(tt))
  # non-negative definite matrix that depends only on the knot differences
  K <- t2$K   # penalty Matrix
  
  
  #Initial values
  sigmae <- initial$sigma2
  D1 <- initial$alphas
  iD1 <- solve(D1)
  iD1 <- (iD1 + t(iD1))/2
  qr <- length(D1[lower.tri(D1, diag = T)])
  lambda <- initial$lambda
  ff <- solve(t(Ns)%*%Ns+lambda*sigmae*K)%*%t(Ns)%*%(y-x%*%beta1)
  W <- cbind(x,Ns)
  gamma1 <- as.vector(c(beta1,ff))
  
  ## Initial values - EM-lambda
  
  t.M <- cbind(rep(1,length(tt)),as.matrix(tt))
  if(is.matrix(beta1)) X.est <- cbind(x,Ns%*%t.M) else X.est <- cbind(Ns%*%t.M)
  
  fit1 <- lm( y ~ -1 + X.est)
  beta.est <- as.vector(fit1$coef) 
  
  Q <- t2$Q
  B <- Q%*%solve(t(Q)%*%Q)
  
  ## correlation structures 
  
  if(struc=="dec")
  {
    phi1 <- initial$phi1
    phi2 <- initial$phi2
    if(is.matrix(beta1)) teta <- c(beta1,ff,sigmae,D1[upper.tri(D1, diag = T)],phi1,phi2) else teta <- c(ff,sigmae,D1[upper.tri(D1, diag = T)],phi1,phi2)
  }
  
  if(struc=="ar")
  {
    phi1 <- initial$phi1
    phi2 <- 1
    if(is.matrix(beta1)) teta <- c(beta1,ff,sigmae,D1[upper.tri(D1, diag = T)],phi1) else teta <- c(ff,sigmae,D1[upper.tri(D1, diag = T)],phi1)
  }
  
  if(struc=="sym")
  {
    phi1 <- initial$phi1
    phi2 <- 0
    if(is.matrix(beta1)) teta <- c(beta1,ff,sigmae,D1[upper.tri(D1, diag = T)],phi1) else teta <- c(ff,sigmae,D1[upper.tri(D1, diag = T)],phi1)
  }
  
  if(struc=="unc")
  {
    phi1 <- NULL
    phi2 <- NULL
    if(is.matrix(beta1)) teta <- c(beta1,ff,sigmae,D1[upper.tri(D1, diag = T)]) else teta <- c(ff,sigmae,D1[upper.tri(D1, diag = T)])
  }
  
  criterio <- 1
  count <- 0
  
  loglik <- loglikslmec(y=y,cc=cc,x=x,z=z,Ns=Ns,ttc=ttc,nj=nj,LL=LL,LU=LU,betas=beta1,ff=ff,sigmae=sigmae,D1=D1,phi1=phi1,phi2=phi2,struc=struc)
  loglikp <- as.numeric(loglik - 0.5*lambda*(t(ff)%*%K%*%ff)) 
  
  while(criterio > precision){
    
    count <- count + 1

    #print(count)
    
    soma1 <- matrix(0,q1,q1)
    soma2 <- 0
    soma3 <- matrix(0,p,p)
    soma4 <- matrix(0,p,1)
    soma5 <- matrix(0,rr,rr)
    soma6 <- matrix(0,rr,1)
    
    soma7 <- matrix(0,p,p)
    soma8 <- matrix(0,p,rr)
    soma9 <- matrix(0,rr,rr)
    
    Infbetasff <- matrix(0,p+rr,p+rr)
    
    yest <- matrix(0,N,1)
    ubi <- matrix(0,m2,m)
    ubbi <- matrix(0,m2,m2)
    uybi <- matrix(0,N,m2)
    uyyi <- matrix(0,N,N)
    uyi <- matrix(0,N,m)
    
    #ver <- matrix(0,m,1)
    
    Zdiag <- matrix(0,N,m2)
    OMEGAdiag <- Vy <- matrix(0,N,N)
    Vlambda <- matrix(0,N,N)
    Einv <- matrix(0,N,N)
    Psi.est <- der.psi <- matrix(0,rr-2+m2,rr-2+m2)
    
    for (j in 1:m){

      cc1 <- cc[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      y1 <- y[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      x1 <- matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
      z1 <- matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
      Ns1 <- matrix(Ns[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=rr)
      tt1 <- ttc[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      W1 <- cbind(x1,Ns1)
      
      LL1 <- LL[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      LU1 <- LU[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      
      muii <- W1%*%gamma1
      Gama <- MatDec(tt1,phi1,phi2,struc)
      invGama <- solve(Gama)
      SIGMA <- (sigmae*Gama + (z1)%*%D1%*%t(z1)) 
      SIGMA <-(SIGMA+t(SIGMA))/2
      SIGMAinv <- solve(SIGMA)
      Lambda1 <- solve(iD1 + (t(z1)%*%invGama%*%z1)*(1/sigmae))
      Lambda1 <- (Lambda1 + t(Lambda1))/2
      
      Zdiag[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),((j-1)*q1+1):(j*q1)] <- z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ]
      OMEGAdiag[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),(sum(nj[1:j-1])+1) : (sum(nj[1:j]))] <- sigmae*Gama
      Einv[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),(sum(nj[1:j-1])+1) : (sum(nj[1:j]))] <- invGama
      Vy[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),(sum(nj[1:j-1])+1) : (sum(nj[1:j]))] <- SIGMA
      
      if(sum(cc1)==0)
      {
        uy <- matrix(y1,nj[j],1)
        uyy <- y1%*%t(y1)
        ub <- (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uy-muii)
        ubb <- Lambda1 + (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uyy-uy%*%t(muii)-muii%*%t(uy)+muii%*%t(muii))%*%t(Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)
        ubb <- (ubb + t(ubb))/2
        uyb <- (uyy-uy%*%t(muii))%*%t(Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)
        
        #ver[j,] <- mvtnorm::dmvnorm(x = as.vector(y1),mean = as.vector(muii),sigma = SIGMA)
      }
      
      if(sum(cc1)>=1)
      {
        
        if(sum(cc1)==nj[j])
        {
          muiic <- W1%*%gamma1
          Sc <- SIGMA
          Sc <- (Sc + t(Sc))/2
          
          aux <- MomTrunc::meanvarTMD(lower = LL1,upper = LU1,mu = muiic,Sigma = Sc, dist = "normal")
          uy <- aux$mean
          uyy <- aux$EYY
          
          ub <- (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uy-muii)
          ubb <- Lambda1 + (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uyy-uy%*%t(muii)-muii%*%t(uy)+muii%*%t(muii))%*%t(Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)
          uyb <- (uyy-uy%*%t(muii))%*%t(Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)
          
          #ver[j,] <- mvtnorm::pmvnorm(lower = LL1,upper = LU1,mean = as.vector(muiic),sigma = Sc)[1]
        }
        else{
          muiic <- W1[cc1==1,]%*%gamma1 + SIGMA[cc1==1,cc1==0]%*%solve(SIGMA[cc1==0,cc1==0])%*%(y1[cc1==0]-W1[cc1==0,]%*%gamma1)
          Sc <- SIGMA[cc1==1,cc1==1]-SIGMA[cc1==1,cc1==0]%*%solve(SIGMA[cc1==0,cc1==0])%*%SIGMA[cc1==0,cc1==1]
          Sc <- (Sc+t(Sc))/2
          
          LL1c <- LL1[cc1==1]
          LU1c <- LU1[cc1==1]

          aux <- MomTrunc::meanvarTMD(lower = LL1c,upper = LU1c,mu = as.vector(muiic), Sigma = Sc, dist = "normal")
          w1aux <- aux$mean
          w2aux <- aux$EYY
          uy <- matrix(y1,nj[j],1)
          uy[cc1==1] <- w1aux
          uyy <- y1%*%t(y1)
          uyy[cc1==0,cc1==1] <- y1[cc1==0]%*%t(w1aux) 
          uyy[cc1==1,cc1==0] <- w1aux%*%t(y1[cc1==0])
          uyy[cc1==1,cc1==1] <- w2aux
          uyy <- (uyy + t(uyy))/2
          ub <- (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uy-muii)
          ubb<- Lambda1 +(Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uyy-uy%*%t(muii)-muii%*%t(uy)+muii%*%t(muii))%*%t(Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)
          ubb <- (ubb + t(ubb))/2
          uyb <- (uyy-uy%*%t(muii))%*%t(Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)
          
          #ver[j,] <- dmvnorm(x = as.vector(y1[cc1==0]),mean = as.vector(muii[cc1==0]),sigma = as.matrix(SIGMA[cc1==0,cc1==0]))*as.numeric(mvtnorm::pmvnorm(lower = as.vector(LL1c),upper = as.vector(LU1c),mean =  as.vector(muiic), sigma = as.matrix(Sc))[1])
        }
        
      } #end if sum(cc1)>=1
      
      soma1 <- soma1 + ubb
      soma2 <- soma2 + (sum(diag(uyy%*%invGama)) - t(uy)%*%invGama%*%muii - t(muii)%*%invGama%*%uy - sum(diag(t(uyb)%*%invGama%*%z1)) - sum(diag(uyb%*%t(z1)%*%invGama))
                        + t(muii)%*%invGama%*%z1%*%ub + t(ub)%*%t(z1)%*%invGama%*%muii + t(muii)%*%invGama%*%muii + sum(diag(ubb%*%t(z1)%*%invGama%*%z1)))
      soma3 <- soma3 + (t(x1)%*%invGama%*%x1)
      soma4 <- soma4 + (t(x1)%*%invGama%*%(uy-z1%*%ub-Ns1%*%ff))
      soma5 <- soma5 + (t(Ns1)%*%invGama%*%Ns1)
      soma6 <- soma6 + (t(Ns1)%*%invGama%*%(uy-z1%*%ub-x1%*%beta1))
      
      soma7 <- soma7 + (t(x1)%*%SIGMAinv%*%x1 - t(x1)%*%SIGMAinv%*%(uyy-uy%*%t(uy))%*%SIGMAinv%*%x1)
      soma8 <- soma8 + (t(x1)%*%SIGMAinv%*%Ns1 - t(x1)%*%SIGMAinv%*%(uyy-uy%*%t(uy))%*%SIGMAinv%*%Ns1)
      soma9 <- soma9 + ((t(Ns1)%*%SIGMAinv%*%Ns1) - t(Ns1)%*%SIGMAinv%*%(uyy-uy%*%t(uy))%*%SIGMAinv%*%Ns1)
      
      ubi[(((j-1)*q1)+1) : (j*q1), j] <- ub
      ubbi[(((j-1)*q1)+1) : (j*q1), (((j-1)*q1)+1) : (j*q1)]<- ubb
      uybi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),(((j-1)*q1)+1) : (j*q1)]<- uyb
      uyyi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]<- uyy
      uyi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),j]<- uy
      yest[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))] <- z1%*%ub + muii
      
    }
    
    if(is.matrix(beta1)) beta1 <- solve(soma3)%*%soma4 else beta1 <- 0
    ff <- solve(soma5 + sigmae*lambda*K)%*%soma6
    gamma1 <- as.vector(c(beta1,ff))
    sigmae <- (1/N)*(soma2)
    sigmae <- as.numeric(sigmae)
    D1 <- (1/m)*(soma1)
    iD1 <- solve(D1)
    
    yest11 <- apply(uyi,1,sum)
    yest[cc==1] <- yest11[cc==1]
    
    ######################################
    ## The information matrix
    ######################################
    
    Infbetasff[1:p,1:p] <- soma7
    Infbetasff[1:p,(p+1):(p+rr)] <- soma8
    Infbetasff[(p+1):(p+rr),1:p] <- t(soma8)
    Infbetasff[(p+1):(p+rr),(p+1):(p+rr)] <- soma9 + (lambda^2)*K%*%ff%*%t(ff)%*%K
    
    Infbetasff <- (Infbetasff + t(Infbetasff))/2      
    
    ######################################
    ## Estimation phi1 and phi2 
    ######################################
    
    if(struc=="dec")                                                                     
    {
      phis <- optim(c(phi1,phi2), FCi,lower =c(0.01,0.01), upper=c(0.9,30),method = "L-BFGS-B", hessian=TRUE, beta1=beta1,ff=ff,sigmae=sigmae,ttc=ttc,ubi=ubi,ubbi=ubbi,uybi=uybi,uyyi=uyyi,uyi=uyi,x=x,z=z,Ns=Ns,nj=nj,struc=struc)$par
      phi1 <- phis[1]
      phi2 <- phis[2]
      if(is.matrix(beta1)) teta1 <- c(beta1,ff,sigmae,D1[upper.tri(D1, diag = T)],phi1,phi2) else teta1 <- c(ff,sigmae,D1[upper.tri(D1, diag = T)],phi1,phi2)
    }
    
    
    if(struc=="ar")
    {
      phi2 <- 1                                   
      phi1 <- optimize(f=FCiphi1, lower= 0.0001, upper=0.9, phi2=phi2, beta1=beta1,ff=ff,sigmae=sigmae,ttc=ttc,ubi=ubi,ubbi=ubbi,uybi=uybi,uyyi=uyyi,uyi=uyi,x=x,z=z,Ns=Ns,nj=nj,struc=struc)$minimum
      if(is.matrix(beta1)) teta1 <- c(beta1,ff,sigmae,D1[upper.tri(D1, diag = T)],phi1) else teta1 <- c(ff,sigmae,D1[upper.tri(D1, diag = T)],phi1)
    }
    
    if(struc=="sym")
    {
      phi2 <- 0
      phi1 <- optimize(f=FCiphi1, lower= 0.0001, upper=0.9, phi2=phi2, beta1=beta1,ff=ff,sigmae=sigmae,ttc=ttc,ubi=ubi,ubbi=ubbi,uybi=uybi,uyyi=uyyi,uyi=uyi,x=x,z=z,Ns=Ns,nj=nj,struc=struc)$minimum
      if(is.matrix(beta1)) teta1 <- c(beta1,ff,sigmae,D1[upper.tri(D1, diag = T)],phi1) else teta1 <- c(ff,sigmae,D1[upper.tri(D1, diag = T)],phi1)
    }
    
    if(struc=="unc")
    {
      if(is.matrix(beta1)) teta1 <- c(beta1,ff,sigmae,D1[upper.tri(D1, diag = T)]) else teta1 <- c(ff,sigmae,D1[upper.tri(D1, diag = T)])
    }

    ######################################
    ## Estimation lambda 
    ######################################
    
    if(lambda.fixed==FALSE)
    {
      ymat <- as.matrix(yest11)
      Dalpha.diag <- kronecker(diag(m),D1)
      N <- dim(ymat)[1]
      
      invOMEGAdiag <- solve(OMEGAdiag)
      
      Q <- t2$Q
      B <- Q%*%solve(t(Q)%*%Q)
      Z.est <- cbind(Ns%*%B,Zdiag)
      
      Psi.est[1:(rr-2),1:(rr-2)] <- (sigmae/lambda)*diag(1,rr-2)
      Psi.est[(rr-1):(rr-2+m2),(rr-1):(rr-2+m2)] <- Dalpha.diag
      invPsi.est <- solve(Psi.est)
      
      V.est <- OMEGAdiag + Z.est%*%Psi.est%*%t(Z.est)
      invV <- solve(V.est)
      
      der.psi[1:(rr-2),1:(rr-2)] <- -(sigmae/lambda^2)*diag(1,rr-2)

      b.est <- Psi.est%*%t(Z.est)%*%invV%*%(ymat - X.est%*%beta.est)
      bb.est <- Psi.est%*%t(Z.est)%*%invV%*%(ymat - X.est%*%beta.est)%*%t((ymat-X.est%*%beta.est))%*%invV%*%Z.est%*%Psi.est + Psi.est - Psi.est%*%t(Z.est)%*%invV%*%Z.est%*%Psi.est
      
      beta.est <- solve(t(X.est)%*%invOMEGAdiag%*%X.est)%*%t(X.est)%*%invOMEGAdiag%*%(ymat - Z.est%*%b.est)
      lambda <- (rr-2)/sum(diag(-invPsi.est%*%der.psi%*%invPsi.est%*%bb.est))
      
    } #end lambda
    
    ########
    
    loglik1 <- loglikslmec(y=y,cc=cc,x=x,z=z,Ns=Ns,ttc=ttc,nj=nj,LL=LL,LU=LU,betas=beta1,ff=ff,sigmae=sigmae,D1=D1,phi1=phi1,phi2=phi2,struc=struc)
    loglikp1 <- as.numeric(loglik1 - 0.5*lambda*(t(ff)%*%K%*%ff)) 
    
    #if(count > 1){criterio <- sqrt(((loglikp1/loglikp)-1)%*%((loglikp1/loglikp)-1))}
    if(count > 1){criterio <- sqrt(((teta1/teta)-1)%*%((teta1/teta)-1))}
    if(count==iter.max){criterio <- precision*0.0001}
    
    teta <- teta1
    loglik <- loglik1
    loglikp <- loglikp1
    
  } #end while
  
  dd <- D1[upper.tri(D1, diag = T)]
  
  if(is.matrix(beta1)) Infbetasff1 <- Infbetasff else Infbetasff1 <- Infbetasff[-p,-p]  
  
  ###########################
  ### Criterion selection ###
  ########################### 
  
  npar <- length(c(teta1))
  AICc <- -2*loglikp + 2*npar
  BICc <- -2*loglik + log(N)*npar
  SIC <- -2*loglikp + npar*log(N)
  
  AICp <- AICnormalp(lambda,loglikp,Ns,K,sigmae,Einv,npar-rr)
  
  obj.out <- list(beta1 = beta1, ff=ff, sigmae = sigmae, dd = dd, phi1=phi1, phi2=phi2,lambda=lambda, loglik=loglik, loglikp=loglikp, AIC=AICc, BIC=BICc, SIC=SIC, 
                  AICp=AICp, iter = count, ubi = ubi, ubbi = ubbi, uybi = uybi, uyi = uyi, uyyi = uyyi, Infbetasff=Infbetasff1, yest=yest, Vy=Vy)  
  
  return(obj.out)
}

