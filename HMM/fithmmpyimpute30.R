library(moveHMM)
setwd("C:/Users/timom")
data <- read.table("imputed_dataframe30min.csv", header = TRUE, sep= ",")

indx <- apply(data, 2, function(x) any(is.na(x) | is.infinite(x)))
indx
# Change Column Names
names(data)[names(data) == "location.long"] <- "x"
names(data)[names(data) == "location.lat"] <- "y"
names(data)[names(data) == "Unnamed..0"] <- "time"
names(data)[names(data) == "individual.local.identifier"] <- "ID"
names(data)[names(data) == "external.temperature"] <- "temp"

#Set time to POSIXct
data$time <- as.POSIXct(strptime(data$time, format='%Y-%m-%d %H:%M:%S'))

# PrepData
data <- prepData(data,coordNames=c("x","y"))

summary(data)

plot(data,compact=T)

# Checking for NAs
any(is.na(data$step))
data$step[is.na(data$step)] <- 1 # Set Step length to 1 if NA

# Normalizing temperature
data$temp <-
  (data$temp-mean(data$temp))/sd(data$temp)

am99 <- data[332287:367367,]
am99dum <- data[332287:367367,]
library (sp)
coordinates(am99dum) = ~x +y

proj4string(am99dum) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +units=m")

am99dum <- spTransform(am99dum,CRS="+proj=utm +zone=36J +datum=WGS84 ")

am99dum <- as.data.frame(coordinates(am99dum))

am99['x'] = am99dum['x'] 
am99['y'] = am99dum['y'] 

# PrepData
data <- prepData(am99,type="UTM",coordNames=c("x","y"))

summary(data)

plot(data,compact=T)

## initial parameters for gamma and von Mises distributions
mu0 <- c(50,250,300,750) # step mean 500
sigma0 <- c(100,100,100,200) # step SD
zeromass0 <- c(0.01,0.05,0.2,0.25) # step zero-mass since some step lengths are 0
stepPar0 <- c(mu0,sigma0,zeromass0)
angleMean0 <- c(pi,0,pi,0) # angle mean
kappa0 <- c(1,1,1,1) # angle concentration
anglePar0 <- c(angleMean0,kappa0)
## call to fitting function
stateNames <- c("Resting","Foraging", "Traveling")
m <- fitHMM(data=data,nbStates=4,stepPar0=stepPar0,
            anglePar0=anglePar0,formula=~temp)
m

CI(m)

plot(m, plotCI=TRUE)

plotStationary(m, plotCI=TRUE)
