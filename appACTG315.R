##################################################################################
##        Semiparametric Mixed-Effects Models with Censored Responses           ##
## ---------------------------------------------------------------------------- ## 
## Article: A semiparametric mixed-effects model for censored longitudinal data ##
## Authors: Thalita B. Mattos, Larissa A. Matos and Victor H. Lachos            ##
## ---------------------------------------------------------------------------- ##
## Last update: 2020-09-15                                                      ##
##################################################################################

## ====================== ##
## Application - ACTG 315 ##
## ====================== ## 

# loading SMEC functions 

source("smec.R")
source("utils_smec.R")

# loading packages 

require(ggplot2)
require(nlme)      # for initial values

# loading dataset

data1 <- read.csv("Wuhivcompletos.csv")
subjects <- unique(data1$Ind)
cluster <- match(data1$Ind,subjects)
m <- c(length(subjects))
N <- length(cluster)
# vector of responses
y <- c(data1$lgcopy)
# vector of censoring
cc <- (y<2) + 0
y[cc==1] <- 2
# design matrix of fixed effects
cd4cov <- (data1$cd4-mean(data1$cd4))/sd(data1$cd4)
x <- as.matrix(cd4cov)
# design matrix of random effects
z <- cbind(1,data1$day1)
# vector of ordered distinct values of the time points
tt <- as.numeric(levels(as.factor(data1$day1)))
# vector of the time points
ttc <- matrix(data1$day1,ncol=1)
# vector with size of cluster
nj <- matrix(0,m,1)
for (j in 1:m){
  nj[j]= sum(cluster==j)
}
# vector of lower limits
LL <- rep(-Inf,length(y))
# vector of upper limits
LU <- as.vector(y)
# incidence matrix
Ns <- matrix(0,nrow=N,ncol=length(tt)) 
for(j in 1:length(tt))
{
  for(i in 1:length(ttc))
  {
    if(ttc[i] == tt[j]){Ns[i,j] <- 1} else {Ns[i,j] <- 0}
  }
}


### ------- ###
### Figures ###
### ------- ### 

dataplot <- data.frame(y=y,cd4=data1$cd4,tt=data1$day,ind=cluster)

figure1a <- ggplot() + 
  geom_line(data=dataplot, aes(x=tt, y=y, group = ind), colour="black") +
  geom_point(data=dataplot, aes(x=tt, y=y, group = ind), size = 1.5,colour="black") + 
  labs(x = "Days after starting treatment", y="log10 HIV-1 RNA") +
  geom_hline(yintercept = log10(100),color="gray50",linetype="dotted",size=0.7) +
  theme_bw()
figure1a

figure1b <- ggplot(data=dataplot, aes(x=cd4, y=y)) + 
  geom_point() + 
  labs(x="CD4 cell count", y="log10 HIV-1 RNA") + 
  theme_bw()
figure1b

## -------------------------
## [X,NT] full rank 
## -------------------------

t.mat <- cbind(matrix(1,length(tt),1),as.matrix(tt))
dim(cbind(x,Ns%*%t.mat))
qr(cbind(x,Ns%*%t.mat))$rank

## ------------------------
## Initial values
## ------------------------ 

dfinitial <- data.frame(cluster = cluster, ttc = ttc, y = y, x=x)

fitN <- lme(y ~ -1 + x, data= dfinitial, random = list(cluster=(~ ttc)))
betasI <- as.vector(fixed.effects(fitN))
sigma2I <- (sigma(fitN))^2
dd <- as.matrix(getVarCov(fitN, type = "random.effects")[1:2,1:2])
alphasI <- matrix(dd,2,2,dimnames = NULL)

fitNAR <- update(fitN,correlation=corAR1())
phi1I <- abs(as.numeric(coef(fitNAR$modelStruct$corStruct, unconstrained=FALSE)))
phi2I <- 2

initial <- list(betas=betasI,sigma2=sigma2I,alphas=alphasI,phi1=phi1I,phi2=phi2I,lambda=100)

## ------------------------
## Fitting the model 
## ------------------------ 

## UNC - uncorrelated 

set.seed(6987)
fitunc <- nsmec(y=y, cc=cc, x=x, z=z, tt=tt, ttc=ttc, nj=nj, LL=LL, LU=LU, Ns=Ns, initial=initial, struc="unc", lambda.fixed=FALSE, iter.max=200,precision=1e-6)

# Estimate
round(fitunc$beta1,4)
round(fitunc$ff,4)
round(fitunc$sigmae,4)
round(fitunc$dd,4)
round(fitunc$lambda,4)
round(fitunc$loglikp,4)
round(fitunc$AICp,4)

# SE
round(sqrt(diag(solve(fitunc$Infbetasff))),4)


## DEC 

set.seed(6987)
fitdec <- nsmec(y=y, cc=cc, x=x, z=z, tt=tt, ttc=ttc, nj=nj, LL=LL, LU=LU, Ns=Ns, initial=initial, struc="dec", lambda.fixed=FALSE, iter.max=200,precision=1e-6)

# Estimate
round(fitdec$beta1,4)
round(fitdec$ff,4)
round(fitdec$sigmae,4)
round(fitdec$dd,4)
round(fitdec$phi1,4)
round(fitdec$phi2,4)
round(fitdec$lambda,4)
round(fitdec$loglikp,4)
round(fitdec$AICp,4)

# SE
round(sqrt(diag(solve(fitdec$Infbetasff))),4)

# ----------------- 

yest_dec <- fitdec$yest

MAEs <- sum(abs(y - yest_dec))/length(y) 
round(MAEs,4)

MSEs <- sum((y - yest_dec)^2)/length(y)
round(MSEs,4)
 

## AR(1)

set.seed(6987)
fitar <- nsmec(y=y, cc=cc, x=x, z=z, tt=tt, ttc=ttc, nj=nj, LL=LL, LU=LU, Ns=Ns, initial=initial, struc="ar", lambda.fixed=FALSE, iter.max=200,precision=1e-6)

# Estimate
round(fitar$beta1,4)
round(fitar$ff,4)
round(fitar$sigmae,4)
round(fitar$dd,5)
round(fitar$phi1,4)
round(fitar$phi2,4)
round(fitar$lambda,4)
round(fitar$loglikp,4)
round(fitar$AICp,4)

# SE
round(sqrt(diag(solve(fitar$Infbetasff))),4)


## CS

set.seed(6987)
fitcs <- nsmec(y=y, cc=cc, x=x, z=z, tt=tt, ttc=ttc, nj=nj, LL=LL, LU=LU, Ns=Ns, initial=initial, struc="sym", lambda.fixed=FALSE, iter.max=200,precision=1e-6)

# Estimate
round(fitcs$beta1,4)
round(fitcs$ff,4)
round(fitcs$sigmae,4)
round(fitcs$dd,4)
round(fitcs$phi1,4)
round(fitcs$phi2,4)
round(fitcs$lambda,4)
round(fitcs$loglikp,4)
round(fitcs$AICp,4)

# SE
round(sqrt(diag(solve(fitcs$Infbetasff))),4)

## ------------------------
## Figures -  DEC model 
## ------------------------ 

dp.fit <- sqrt(diag(solve(fitdec$Infbetasff)))
dataff <- data.frame(tt=tt,fest=fitdec$ff,lwr=fitdec$ff-1.96*dp.fit[2:11], uppr=fitdec$ff+1.96*dp.fit[2:11] )

figure9a <- ggplot(dataff, aes(tt,fest)) + 
  geom_ribbon(data=dataff, aes(ymin=lwr,ymax=uppr,x = tt), fill = "grey80", alpha=0.8) + 
  geom_line(data = dataff, aes(x=tt,y=fest), colour="black", size=0.7) +
  scale_y_continuous("Fitted nonparametric function", limits = c(2,6.5)) + 
  labs(x = "Days after starting treatment") +
  theme_bw()
figure9a


datafit <- data.frame(cluster=cluster, y=y, yfit= fitdec$yest, tt=ttc)
datafit1 <- datafit[which(datafit$cluster %in% c(2,17,19,24,37,39)),]

figure9b <- ggplot() +  
  geom_line(data=datafit1, aes(x=tt, y=y), colour="black") +
  geom_point(data=datafit1, aes(x=tt, y=y, group = cluster), size = 1.5,colour="black") + 
  geom_line(data=datafit1, aes(x=tt, y=yfit), colour="red", linetype="dotted") +
  geom_point(data=datafit1, aes(x=tt, y=yfit, group = cluster), size = 1.5,colour="red") +
  scale_x_continuous("Days after starting treatment", limits = c(0,200)) +
  facet_wrap(~ cluster) +
  labs(y="log10 HIV-1 RNA") +
  geom_hline(yintercept = log10(100),color="gray50",linetype="dotted",size=0.7) +
  theme_bw() 
figure9b

#######################################################################################



