library(moveHMM)
setwd("C:/Users/timom/anaconda3/envs/madeleine_project")
data <- read.table("imputed_data18004.csv", header = TRUE, sep= ",") #"imputed_data18002.csv"

# Change Column Names
names(data)[names(data) == "long"] <- "x"
names(data)[names(data) == "lat"] <- "y"
names(data)[names(data) == "X"] <- "time"
names(data)[names(data) == "individual.local.identifier"] <- "ID"

#Select AM99
am99 <- data[65218:71065,] #data[1:35083,]

#Change time to POSIXct
am99$time <- as.POSIXct(strptime(am99$time, format='%Y-%m-%d %H:%M:%S'))

# PrepData
am99 <- prepData(am99,coordNames=c("x","y"))

#am105$step <- am105$step*1000

summary(am99)

plot(am99,compact=T)

any(is.na(am99$step))
am99[is.na(am99)] <- 1
any(is.na(am99$step))

## initial parameters for gamma and von Mises distributions
mu0 <- c(1,2) # step mean (two parameters: one for each state)
sigma0 <- c(3,4) # step SD
zeromass0 <- c(0.1,0.05) # step zero-mass
stepPar0 <- c(mu0,sigma0,zeromass0)
angleMean0 <- c(pi,0) # angle mean
kappa0 <- c(1,1) # angle concentration
anglePar0 <- c(angleMean0,kappa0)
## call to fitting function
m <- fitHMM(data=am99,nbStates=2,stepPar0=stepPar0,
            anglePar0=anglePar0,formula=~temp)
m

CI(m)

plot(m, plotCI=TRUE)
##################################################################
#2nd try
##################################################################

# histogram
hist(am99$step,prob=T,breaks=40)

## starting values
theta <- c(0.5,0.5,1,2,3,4)
theta.star <- c(qlogis(theta[1]),qlogis(theta[2]),log(theta[3]),log(theta[4]),log(theta[5]),log(theta[6]))

# ML Evaluation Function
mllk<-function(theta.star,x){
  theta <- c(plogis(theta.star[1]),plogis(theta.star[2]),
             exp(theta.star[3]),exp(theta.star[4]),
             exp(theta.star[5]),exp(theta.star[6]))
  Gamma <- diag(theta[1:2])
  Gamma[1,2] <- 1-Gamma[1,1]
  Gamma[2,1] <- 1-Gamma[2,2] 
  delta <- solve(t(diag(2)-Gamma+1),c(1,1))
  mu <- theta[3:4]
  sigma <- theta[5:6]
  allprobs <- matrix(1,length(x),2)
  ind<-which(!is.na(x))
  allprobs[ind,] <- cbind(
    dgamma(x[ind],shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1]),
    dgamma(x[ind],shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2]))
  foo <- delta%*%diag(allprobs[1,])
  l <- log(sum(foo))
  phi <- foo/sum(foo)
  for (t in 2:length(x)){
    foo <- phi%*%Gamma%*%diag(allprobs[t,])
    l <- l+log(sum(foo))
    phi <- foo/sum(foo)
  }
  return(-l)
}

# Fitting model 
mod <- nlm(mllk,theta.star,x=am99$step,print.level=2)
mod$estimate

#save model parameters (in order to retransform)
theta.star.mle <- mod$estimate

# retransform:
theta.mle <- c(plogis(theta.star.mle[1]),plogis(theta.star.mle[2]),
               exp(theta.star.mle[3]),exp(theta.star.mle[4]),
               exp(theta.star.mle[5]),exp(theta.star.mle[6]))
Gamma.mle <- diag(theta.mle[1:2])
Gamma.mle[1,2] <- 1-Gamma.mle[1,1]
Gamma.mle[2,1] <- 1-Gamma.mle[2,2]
delta.mle <- solve(t(diag(2)-Gamma.mle+1),c(1,1))
mu.mle <- theta.mle[3:4]
sigma.mle <- theta.mle[5:6]

# ml estimates
Gamma.mle
delta.mle
mu.mle
sigma.mle

# histogram
hist(am99$step, probability=TRUE, main='Histogram of step length',
     xlab='step length in metres', breaks=40, col="light grey")

# add the state-dependent distributions
z <- seq(0,10,by=0.01)
lines(z, theta.mle[1]*dgamma(z,shape=theta.mle[3]^2/theta.mle[5]^2,scale=theta.mle[5]^2/theta.mle[3]),col='orange',lwd=2)
lines(z, theta.mle[2]*dgamma(z,shape=theta.mle[4]^2/theta.mle[6]^2,scale=theta.mle[6]^2/theta.mle[4]),col='lightblue',lwd=2)


# add the overall density (sum of the state-dependent densities)
lines(z, theta.mle[1]*dgamma(z,shape=theta.mle[3]^2/theta.mle[5]^2,scale=theta.mle[5]^2/theta.mle[3])
      + theta.mle[2]*dgamma(z,shape=theta.mle[4]^2/theta.mle[6]^2,scale=theta.mle[6]^2/theta.mle[4]),lwd=2, lty=2) 

legend("topright", inset = c(0.05,0.05), legend=c("state 1","state 2", "overall"), lwd=2, col = c("orange","lightblue", "black"), lty = c(1,1,2))
