library(moveHMM)
setwd("C:/Users/timom")
datas <- read.table("imputed_dataframe30min.csv", header = TRUE, sep= ",")

#indx <- apply(data, 2, function(x) any(is.na(x) | is.infinite(x)))
#indx
# Change Column Names
names(datas)[names(datas) == "location.long"] <- "x"
names(datas)[names(datas) == "location.lat"] <- "y"
names(datas)[names(datas) == "Unnamed..0"] <- "time"
names(datas)[names(datas) == "individual.local.identifier"] <- "ID"
names(datas)[names(datas) == "external.temperature"] <- "temp"

#Set time to POSIXct
datas$time <- as.POSIXct(strptime(datas$time, format='%Y-%m-%d %H:%M:%S'))

# PrepData
#data <- prepData(data,coordNames=c("x","y"))

#summary(data)

#plot(data,compact=T)

# Checking for NAs
#any(is.na(data$step))
#data$step[is.na(data$step)] <- 1 # Set Step length to 1 if NA

# Normalizing temperature
datas$temp <-
  (datas$temp-mean(datas$temp))/sd(datas$temp)

am93 <- datas[297205:332285,]
am93dum <- datas[297205:332285,]
library (sp)
coordinates(am93dum) = ~x +y

proj4string(am93dum) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

am93dum <- spTransform(am93dum,CRS="+proj=utm +zone=36J +datum=WGS84 ")

am93dum <- as.data.frame(coordinates(am93dum))

am93['x'] = am93dum['x'] 
am93['y'] = am93dum['y'] 

# PrepData
data <- prepData(am93,type="UTM",coordNames=c("x","y"))

summary(data)

plot(data,compact=T)

data$step[is.na(data$step)] <- 200

# Indices of steps of length zero
whichzero <- which(data$step == 0)
# Proportion of steps of length zero in the data set
length(whichzero)/nrow(data)


## initial parameters for gamma and von Mises distributions
mu0 <- c(21.22671,136.8682,483.87) # step mean 500
sigma0 <- c(17.90707,84.60687,289.3098) # step SD
#zeromass0 <- c(0.000497791,0.00000000000134072,0.000000002520846) # step zero-mass since some step lengths are 0
stepPar0 <- c(mu0,sigma0)
angleMean0 <- c(0,0.1,0) # angle mean
kappa0 <- c(0.355384957,1.53831499,2.59417278) # angle concentration
anglePar0 <- c(angleMean0,kappa0)
## call to fitting function
stateNames <- c("Resting","Foraging", "Traveling")
m <- fitHMM(data=data, nbStates=3,stepPar0=stepPar0,
            anglePar0=anglePar0,formula=~temp)
m

CI(m)

plot(m, plotCI=TRUE)

plotStationary(m, plotCI=TRUE)

###################################
#########################################################################
#parallel
#######
set.seed(34523456)
# Package for parallel computations
library(parallel)
# Create cluster of size ncores
ncores <- detectCores()-1
cl <- makeCluster(getOption("cl.cores", ncores))
# Export objects needed in parallelised function to cluster
clusterExport(cl, list("data", "fitHMM"))
# Number of tries with different starting values
niter <- 6
# Create list of starting values
allPar0 <- lapply(as.list(1:niter), function(x) {
  # Step length mean
  stepMean0 <- runif(3,
                     min = c(10, 150, 350),
                     max = c(30, 350, 700))
  # Step length standard deviation
  stepSD0 <- runif(3,
                   min = c(10, 150, 350),
                   max = c(30, 350, 700))
  
  #  zeromass0 <- c(0.0001,0.0001,0.0001,0.0001)
  # Turning angle mean
  angleMean0 <- c(0, 0, 0)
  # Turning angle concentration
  angleCon0 <- runif(3,
                     min = c(0.1, 1, 3),
                     max = c(1, 2, 4))
  # Return vectors of starting values
  stepPar0 <- c(stepMean0, stepSD0)
  anglePar0 <- c(angleMean0, angleCon0)
  return(list(step = stepPar0, angle = anglePar0))
})

# Fit the niter models in parallel
allm_parallel <- parLapply(cl = cl, X = allPar0, fun = function(par0) {
  m <- fitHMM(data = data, nbStates = 3, stepPar0 = par0$step,
              anglePar0 = par0$angle, formula=~temp)
  return(m)
})
# Then, we can extract the best-fitting model from allm_parallel
# as before

allnllk <- unlist(lapply(allm_parallel, function(m) m$mod$minimum))
allnllk

whichbest <- which.min(allnllk)
whichbest
mbest <- allm_parallel[[whichbest]]
mbest
mbest$mod
CI(mbest)

plot(mbest, plotCI=TRUE)

plotStationary(mbest, plotCI=TRUE)
