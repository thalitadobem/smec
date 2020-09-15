##################################################################################
##        Semiparametric Mixed-Effects Models with Censored Responses           ##
## ---------------------------------------------------------------------------- ## 
## Article: A semiparametric mixed-effects model for censored longitudinal data ##
## Authors: Thalita B. Mattos, Larissa A. Matos and Victor H. Lachos            ##
## ---------------------------------------------------------------------------- ##
## Last update: 2020-09-15                                                      ##
##################################################################################

# ------------------------------------------- #
# packages needed for estimating a SMEC model #
# ------------------------------------------- #

library(mvtnorm)
library(MomTrunc)

# ----------------------------------------------- #
# Auxiliary functions for estimating a SMEC model #
# ----------------------------------------------- #

## -------------------------------------- ##
## Natural Cubic Spline                   ## 
## (adapted function of the ssym package) ##
## -------------------------------------- ##
 
natural.cubic.spline <- function (xx, lambda, nknots, all.knots)
{
  xx <- as.matrix(round(xx, digits = 6))
  difv <- as.matrix(as.numeric(levels(factor(xx))))
  if (length(difv) < 3) 
    stop("Variable in natural cubic spline does not have at least three different values!!", 
         call. = FALSE)
  if (!is.numeric(xx)) 
    stop("Variable in natural cubic spline must be numeric!!", 
         call. = FALSE)
  if (missing(all.knots)) 
    all.knots <- FALSE
  if (all.knots == FALSE) {
    if (missing(nknots)) {
      nknots <- floor(length(xx)^(1/3)) + 3
      xk <- quantile(xx, prob = seq(0, 1, length = nknots))
      difv2 <- as.matrix(as.numeric(levels(factor(xk))))
      while (length(difv2) < nknots) {
        nknots <- nknots - 1
        xk <- quantile(xx, prob = seq(0, 1, length = nknots))
        difv2 <- as.matrix(as.numeric(levels(factor(xk))))
      }
      if (length(difv) < nknots) {
        nknots <- length(difv)
        xk <- difv
      }
    }
    else {
      if (floor(nknots) < 3) 
        stop("Number of knots must be an integer >= 3!!", 
             call. = FALSE)
      nknots <- floor(nknots)
      xk <- quantile(xx, prob = seq(0, 1, length = nknots))
      difv2 <- as.matrix(as.numeric(levels(factor(xk))))
      if (length(difv2) < nknots) 
        stop("Too many knots!!", call. = FALSE)
    }
  }
  else {
    nknots <- length(difv)
    xk <- difv
  }
  if (!missing(lambda)) {
    if (lambda <= 0) 
      stop("Smoothing parameter must be a positive value!!", 
           call. = FALSE)
    else status <- "known"
  }
  else {
    status <- "unknown"
    lambda <- 1
  }
  n <- length(xx)
  m <- nknots
  h <- matrix(0, m - 1, 1)
  Q <- matrix(0, m, m - 2)
  R <- matrix(0, m - 2, m - 2)
  for (i in 1:(m - 1)) {
    h[i] <- xk[i + 1] - xk[i]
  }
  for (j in 2:(m - 1)) {
    Q[j - 1, j - 1] <- 1/h[j - 1]
    Q[j, j - 1] <- -1/h[j - 1] - 1/h[j]
    Q[j + 1, j - 1] <- 1/h[j]
    R[j - 1, j - 1] <- (h[j - 1] + h[j])/3
  }
  for (j in 2:(m - 2)) {
    R[j - 1, j] <- (h[j])/6
    R[j, j - 1] <- (h[j])/6
  }
  K <- Q %*% solve(R) %*% t(Q)
  for (j in 2:(m - 2)) {
    R[j - 1, j] <- (h[j])/6
    R[j, j - 1] <- (h[j])/6
  }
  K <- Q %*% solve(R) %*% t(Q)
  return(list(K=K, Q=Q,R=R))
}

## ----------------------- ##
## construction matrix DEC ##    
## ----------------------- ##

MatDec <- function(tt,phi1,phi2,struc){
  r <- length(tt)
  
  if(struc=="dec" || struc=="ar"){
    if(phi2<=0.0000001){
      W <- matrix(phi1,nrow=r,ncol=r)
      for (i in 1:r){W[i,i]<- 1}
      V <- W }
    else{
      H <- (abs(outer(tt, tt, "-")))^phi2
      V <- (phi1^H)}
  } 
  if(struc=="sym"){
    W <- matrix(phi1,nrow=r,ncol=r)
    diag(W)<-1
    V <- W
  }
  if(struc=="ma"){
    W <- matrix(0,nrow=r,ncol=r)
    for (i in 1:r){
      W[i,i]<- 1
      for(j in 1:r){ 
        dif <- abs(tt[i]-tt[j])
        if(dif==1){W[i,j]= phi1}}}
    V <- W
  }
  if(struc=="unc"){
    W <- diag(1,nrow=r,ncol=r)
    V <- W
  }
  return(V)
}

## --------------------------------- ##
## Function used in the maximization ##
## --------------------------------- ##

## Estimate phi1 and phi2

FCi <- function(phiG,beta1,ff,sigmae,ttc,ubi,ubbi,uybi,uyyi,uyi,x,z,Ns,nj,struc)
{
  if(struc=="dec"){phi1 <- phiG[1]; phi2 <- phiG[2]}
  if(struc=="ar"){phi2 <- 1; phi1 <- phiG}
  if(struc=="sym"){phi2 <- 0; phi1 <- phiG}

  m <- length(nj)[1]
  p <- dim(x)[2]
  q1 <- dim(z)[2]
  r <- dim(Ns)[2]
  gamma1 <- as.vector(c(beta1,ff))
  soma <- 0
  
  for (j in 1:m ){
    
    x1 <- matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
    z1 <- matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
    Ns1 <- matrix(Ns[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=r)
    W1 <- cbind(x1,Ns1)
    muii <- W1%*%gamma1
    tt1 <- ttc[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
    
    ub <- ubi[(((j-1)*q1)+1) : (j*q1), j]
    ubb <- ubbi[(((j-1)*q1)+1) : (j*q1), (((j-1)*q1)+1) : (j*q1)]
    uyb <- uybi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),(((j-1)*q1)+1) : (j*q1)]
    uyy <- uyyi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
    uy <- uyi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),j]
    
    Cii <- MatDec(tt1,phi1,phi2,struc) 
    Cii <- (Cii + t(Cii))/2
    if(det(Cii)<=0){A <- 1}else{A <- det(Cii)}
    invCii <- solve(Cii)
    
    Ai <-  as.vector(sum(diag(uyy%*%invCii)) -t(uy)%*%invCii%*%muii - t(muii)%*%invCii%*%uy - sum(diag(invCii%*%((uyb)%*%t(z1)))) - sum(diag(invCii%*%((uyb)%*%t(z1))))
                     + t(muii)%*%invCii%*%z1%*%ub + t(ub)%*%t(z1)%*%invCii%*%muii + t(muii)%*%invCii%*%muii + sum(diag(ubb%*%t(z1)%*%invCii%*%z1)))
    
    soma <- soma - 0.5*log(A) - (0.5/sigmae)*Ai                
    
  }
  
  return(-soma)
}

## Estimate phi1 

FCiphi1 <- function(phi1,phi2,beta1,ff,sigmae,ttc,ubi,ubbi,uybi,uyyi,uyi,x,z,Ns,nj,struc)
{

  m <- length(nj)[1]
  p <- dim(x)[2]
  q1 <- dim(z)[2]
  r <- dim(Ns)[2]
  gamma1 <- as.vector(c(beta1,ff))
  soma <- 0
  
  for (j in 1:m ){
    
    x1 <- matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
    z1 <- matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
    Ns1 <- matrix(Ns[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=r)
    W1 <- cbind(x1,Ns1)
    muii <- W1%*%gamma1
    tt1 <- ttc[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
    
    ub <- ubi[(((j-1)*q1)+1) : (j*q1), j]
    ubb <- ubbi[(((j-1)*q1)+1) : (j*q1), (((j-1)*q1)+1) : (j*q1)]
    uyb <- uybi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),(((j-1)*q1)+1) : (j*q1)]
    uyy <- uyyi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
    uy <- uyi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),j]
    
    Cii <- MatDec(tt1,phi1,phi2,struc) 
    Cii <- (Cii + t(Cii))/2
    if(det(Cii)<=0){A <- 1}else{A <- det(Cii)}
    invCii <- solve(Cii)
    
    Ai <-  as.vector(sum(diag(uyy%*%invCii)) -t(uy)%*%invCii%*%muii - t(muii)%*%invCii%*%uy - sum(diag(invCii%*%((uyb)%*%t(z1)))) - sum(diag(invCii%*%((uyb)%*%t(z1))))
                     + t(muii)%*%invCii%*%z1%*%ub + t(ub)%*%t(z1)%*%invCii%*%muii + t(muii)%*%invCii%*%muii + sum(diag(ubb%*%t(z1)%*%invCii%*%z1)))
    
    soma <- soma - 0.5*log(A) - (0.5/sigmae)*Ai                
    
  }
  
  return(-soma)
}

## ---- ##
## AICp ##
## ---- ##

AICnormalp <- function(lambda,loglikp,Ns,K,sigmae,Einv,npar)
{
  S <- solve(t(Ns)%*%Einv%*%Ns + lambda*sigmae*K)%*%t(Ns)%*%Einv
  df <- sum(diag(Ns%*%S))
  AICp <- -2*loglikp + 2*(npar + df)
  return(AICp)
}

## ----------------------- ##
## log-likelihood - normal ##
## ----------------------- ##

loglikslmec <-  function(y,cc,x,z,Ns,ttc,nj,LL,LU,betas,ff,sigmae,D1,phi1,phi2,struc)
{
  
  m <- length(nj)[1]
  p <- dim(x)[2]
  q1 <- dim(z)[2]
  m1 <- m*p
  m2 <- m*q1
  rr <- dim(Ns)[2]
  
  gamma1 <- as.vector(c(betas,ff))
  iD1 <- solve(D1)
  iD1 <- (iD1 + t(iD1))/2
  
  ver <- matrix(0,m,1)
  
  for(j in 1:m)
  {
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
    
    if(sum(cc1)==0)
    {
      ver[j,] <- mvtnorm::dmvnorm(x = as.vector(y1),mean = as.vector(muii),sigma = SIGMA)
    }
    if(sum(cc1)>=1)
    {
      
      if(sum(cc1)==nj[j])
      {
        muiic <- W1%*%gamma1
        Sc <- SIGMA
        Sc <- (Sc + t(Sc))/2
        ver[j,] <- mvtnorm::pmvnorm(lower = LL1,upper = LU1,mean = as.vector(muiic),sigma = Sc)[1]
      }
      else{
        
        muiic <- W1[cc1==1,]%*%gamma1 + SIGMA[cc1==1,cc1==0]%*%solve(SIGMA[cc1==0,cc1==0])%*%(y1[cc1==0]-W1[cc1==0,]%*%gamma1)
        Sc <- SIGMA[cc1==1,cc1==1]-SIGMA[cc1==1,cc1==0]%*%solve(SIGMA[cc1==0,cc1==0])%*%SIGMA[cc1==0,cc1==1]
        Sc <- (Sc+t(Sc))/2
        
        LL1c <- LL1[cc1==1]
        LU1c <- LU1[cc1==1]
        
        ver[j,] <- dmvnorm(x = as.vector(y1[cc1==0]),mean = as.vector(muii[cc1==0]),sigma = as.matrix(SIGMA[cc1==0,cc1==0]))*as.numeric(mvtnorm::pmvnorm(lower = as.vector(LL1c),upper = as.vector(LU1c),mean =  as.vector(muiic), sigma = as.matrix(Sc))[1])
      }
      
    } #end if sum(cc1)>=1
  } #end for 
  
  logvero <- sum(log(ver))
  
  return(logvero)
}

