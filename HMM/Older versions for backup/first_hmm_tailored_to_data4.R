library(moveHMM)
setwd("C:/Users/timom/anaconda3/envs/madeleine_project")
data <- read.table("imputed_data18004.csv", header = TRUE, sep= ",")

# Change Column Names
names(data)[names(data) == "long"] <- "x"
names(data)[names(data) == "lat"] <- "y"
names(data)[names(data) == "X"] <- "time"
names(data)[names(data) == "individual.local.identifier"] <- "ID"

#Set time to POSIXct
data$time <- as.POSIXct(strptime(data$time, format='%Y-%m-%d %H:%M:%S'))

# PrepData
data <- prepData(data,coordNames=c("x","y"))

summary(data)

plot(data,compact=T)

# Checking for NAs
any(is.na(am99$step))
data$step[is.na(data$step)] <- 1 # Set Step length to 1 if NA

# Normalizing temperature
data$temp <-
  (data$temp-mean(data$temp))/sd(data$temp)

## initial parameters for gamma and von Mises distributions
mu0 <- c(1,0.3,5) # step mean
sigma0 <- c(5,2,10) # step SD
zeromass0 <- c(0.1,0.05,0.2) # step zero-mass since some step lengths are 0
stepPar0 <- c(mu0,sigma0,zeromass0)
angleMean0 <- c(0,pi,0) # angle mean
kappa0 <- c(1,1,1) # angle concentration
anglePar0 <- c(angleMean0,kappa0)
## call to fitting function
m <- fitHMM(data=data,nbStates=3,stepPar0=stepPar0,
            anglePar0=anglePar0,formula=~temp)
m

CI(m)

plot(m, plotCI=TRUE)

plotStationary(m, plotCI=TRUE)