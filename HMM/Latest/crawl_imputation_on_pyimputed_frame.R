library(crawl)
library(tidyverse)
library(lubridate)
library(sp)
library(sf)
library(momentuHMM)

setwd("C:/Users/timom")
data <- read.table("imputed_dataframe30min.csv", header = TRUE, sep= ",")

names(data)[names(data) == "location.long"] <- "x"
names(data)[names(data) == "location.lat"] <- "y"
names(data)[names(data) == "Unnamed..0"] <- "time"
names(data)[names(data) == "individual.local.identifier"] <- "ID"
names(data)[names(data) == "external.temperature"] <- "temp"
names(data)[names(data) == "X"] <- "numer"

data <- read.table("crawl_data60.csv", header = TRUE, sep= ",")

data$time <- as.POSIXct(strptime(data$time, format='%Y-%m-%d %H:%M:%S'))
data <- data[ -c(2,3) ]

indx <- apply(data, 2, function(x) any(is.na(x) | is.infinite(x)))
colnames[indx]

data <- na.omit(data)


df <- data[ -c(3, 8,9,10, 12) ]

############################################################################
datas <- read.table("imputed_dataframe30min.csv", header = TRUE, sep= ",")

names(datas)[names(datas) == "location.long"] <- "x"
names(datas)[names(datas) == "location.lat"] <- "y"
names(datas)[names(datas) == "Unnamed..0"] <- "time"
names(datas)[names(datas) == "individual.local.identifier"] <- "ID"
names(datas)[names(datas) == "external.temperature"] <- "temp"
names(datas)[names(datas) == "X"] <- "numer"

datas$time <- as.POSIXct(strptime(datas$time, format='%Y-%m-%d %H:%M:%S'))

indx <- apply(datas, 2, function(x) any(is.na(x) | is.infinite(x)))
colnames[indx]

datas <- na.omit(datas)

indx <- apply(datas, 2, function(x) any(is.na(x) | is.infinite(x)))
colnames[indx]
indx

coordinates(datas) = ~x +y

proj4string(datas) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +units=m")

datas <- spTransform(datas,CRS="+proj=utm +zone=36J +datum=WGS84 ")

datas <- as.data.frame(coordinates(datas))
###############################################################################

df['x'] = datas['x'] #erst untitled 8 laufen lassen
df['y'] = datas['y'] # dort wurde longlat in utm umgewandelt

# not working with kalman_inits (initial.state)
crwOut <- crawlWrap(df, timeStep = "1 hour", theta=c(2, -2), fixPar=c(NA,NA)) #theta=c(6.855, -0.007)

ggplot(df, aes(x, y)) + geom_point() + geom_path() +
  geom_point(aes(mu.x, mu.y), crwOut$crwPredict, col = "firebrick", size = 0.5) +
  geom_path(aes(mu.x, mu.y), crwOut$crwPredict, col = "firebrick", size = 0.3)

data <- prepData(data,type="UTM",coordNames=c("x","y")) #crwOut if imputed in this session

summary(data)

plot(data,compact=T)

indx <- apply(data, 2, function(x) any(is.na(x) | is.infinite(x)))
indx # No NAs in x, y and time column

# klappt
##############################################################################
am99 <- data[134933:152473,] #varies
am99$step[is.na(am99$step)] <- 200
am99$angle[is.na(am99$angle)] <- 0
whichzero <- which(am99$step == 0)
# Proportion of steps of length zero in the data set
length(whichzero)/nrow(am99)

any(is.na(am99$temp))
nbStates <- 3
stepDist <- "gamma" # step distribution
angleDist <- "vm" # turning angle distribution

am99$temp[is.na(am99$temp)] <- 25
## initial parameters for gamma and von Mises distributions
mu0 <- c(21.22671,136.8682,483.87) # step mean 500
sigma0 <- c(17.90707,84.60687,289.3098) # step SD
zeromass0 <- c(0.000497791,0.00000000000134072,0.000000002520846)
stepPar0 <- c(mu0,sigma0,zeromass0)
angleMean0 <- c(0.003592673,0.01941996,0.01499825) # angle mean
kappa0 <- c(0.355384957,1.53831499,2.59417278) # angle concentration
anglePar0 <- c(angleMean0)
## call to fitting function
stateNames <- c("Resting","Foraging", "Traveling")
m <- fitHMM(data=am99, nbStates=nbStates,dist=list(step=stepDist,angle=angleDist), Par0=list(step=stepPar0,angle=anglePar0),formula=~temp, stateNames=stateNames)
m
CIreal(m)

plot(m, plotCI=TRUE, plotStationary=TRUE,lwd=1)

#########################################################
am239 <- data[52625:70164,] #varies
am239$step[is.na(am239$step)] <- 200
am239$angle[is.na(am239$angle)] <- 0
whichzero <- which(am239$step == 0)
# Proportion of steps of length zero in the data set
length(whichzero)/nrow(am239)

any(is.na(am239$temp))
nbStates <- 3
stepDist <- "gamma" # step distribution
angleDist <- "vm" # turning angle distribution

am239$temp[is.na(am239$temp)] <- 25
## initial parameters for gamma and von Mises distributions
mu0 <- c(21.22671,136.8682,483.87) # step mean 500 21.22671,136.8682,483.87
sigma0 <- c(17.90707,84.60687,289.3098) # step SD 17.90707,84.60687,289.3098
zeromass0 <- c(0.0002280502,0.0002280502,0.0002280502) #0.0002280502,0.0002280502,0.0002280502
stepPar0 <- c(mu0,sigma0,zeromass0)
angleMean0 <- c(0.003592673,0.01941996,0.01499825) # angle mean
kappa0 <- c(0.355384957,1.53831499,2.59417278) # angle concentration
anglePar0 <- c(angleMean0)
## call to fitting function
stateNames <- c("Resting","Foraging", "Traveling")
m <- fitHMM(data=am239, nbStates=nbStates,dist=list(step=stepDist,angle=angleDist), Par0=list(step=stepPar0,angle=anglePar0),formula=~temp, stateNames=stateNames)
m
CIreal(m)

plot(m, plotCI=TRUE, plotStationary=TRUE,lwd=1)
##########################################
am253<- data[71517:81805,] #varies, some outliers at the beginning
am253$step[is.na(am253$step)] <- 200
am253$angle[is.na(am253$angle)] <- 0
whichzero <- which(am253$step == 0)
# Proportion of steps of length zero in the data set
length(whichzero)/nrow(am253)

any(is.na(am253$temp))
nbStates <- 3
stepDist <- "gamma" # step distribution
angleDist <- "vm" # turning angle distribution

am253$temp[is.na(am253$temp)] <- 25
## initial parameters for gamma and von Mises distributions
mu0 <- c(21.22671,136.8682,483.87) # step mean 500 21.22671,136.8682,483.87
sigma0 <- c(17.90707,84.60687,289.3098) # step SD 17.90707,84.60687,289.3098
zeromass0 <- c(9.719118e-05,9.719118e-05,9.719118e-05) #0.0002280502,0.0002280502,0.0002280502
stepPar0 <- c(mu0,sigma0,zeromass0)
angleMean0 <- c(0.003592673,0.01941996,0.01499825) # angle mean
kappa0 <- c(0.355384957,1.53831499,2.59417278) # angle concentration
anglePar0 <- c(angleMean0)
## call to fitting function
stateNames <- c("Resting","Foraging", "Traveling")
m <- fitHMM(data=am253, nbStates=nbStates,dist=list(step=stepDist,angle=angleDist), Par0=list(step=stepPar0,angle=anglePar0),formula=~temp, stateNames=stateNames)
m
CIreal(m)

plot(m, breaks=20, plotCI=TRUE, plotStationary=TRUE,lwd=1)
######################################
am255<- data[81808:92271,] #varies, some outliers at the beginning
am255$step[is.na(am255$step)] <- 200
am255$angle[is.na(am255$angle)] <- 0
whichzero <- which(am255$step == 0)
# Proportion of steps of length zero in the data set
length(whichzero)/nrow(am255)
am255 <- am255[-c(3357:3913),] 
any(is.na(am255$temp))
nbStates <- 3
stepDist <- "gamma" # step distribution
angleDist <- "vm" # turning angle distribution

am255$temp[is.na(am255$temp)] <- 25
## initial parameters for gamma and von Mises distributions
mu0 <- c(21.22671,136.8682,483.87) # step mean 500 21.22671,136.8682,483.87
sigma0 <- c(17.90707,84.60687,289.3098) # step SD 17.90707,84.60687,289.3098
stepPar0 <- c(mu0,sigma0)
angleMean0 <- c(0.003592673,0.01941996,0.01499825) # angle mean
kappa0 <- c(0.355384957,1.53831499,2.59417278) # angle concentration
anglePar0 <- c(angleMean0)
## call to fitting function
stateNames <- c("Resting","Foraging", "Traveling")
m <- fitHMM(data=am255, nbStates=nbStates,dist=list(step=stepDist,angle=angleDist), Par0=list(step=stepPar0,angle=anglePar0),formula=~temp, stateNames=stateNames)
m
CIreal(m)

plot(m, plotCI=TRUE, plotStationary=TRUE,lwd=1)
########################
am307<- data[92274:104936,] #varies, some outliers at the beginning
am307$step[is.na(am307$step)] <- 200
am307$angle[is.na(am307$angle)] <- 0
whichzero <- which(am307$step == 0)
# Proportion of steps of length zero in the data set
length(whichzero)/nrow(am307)

any(is.na(am307$temp))
nbStates <- 3
stepDist <- "gamma" # step distribution
angleDist <- "vm" # turning angle distribution

am307$temp[is.na(am307$temp)] <- 25
## initial parameters for gamma and von Mises distributions
mu0 <- c(21.22671,136.8682,483.87) # step mean 500 21.22671,136.8682,483.87
sigma0 <- c(17.90707,84.60687,289.3098) # step SD 17.90707,84.60687,289.3098
zeromass0 <- c(0.0001579405,0.0001579405,0.0001579405) #0.0002280502,0.0002280502,0.0002280502
stepPar0 <- c(mu0,sigma0,zeromass0)
angleMean0 <- c(0.003592673,0.01941996,0.01499825) # angle mean
kappa0 <- c(0.355384957,1.53831499,2.59417278) # angle concentration
anglePar0 <- c(angleMean0)
## call to fitting function
stateNames <- c("Resting","Foraging", "Traveling")
m <- fitHMM(data=am307, nbStates=nbStates,dist=list(step=stepDist,angle=angleDist), Par0=list(step=stepPar0,angle=anglePar0),formula=~temp, stateNames=stateNames)
m
CIreal(m)

plot(m, breaks=20, plotCI=TRUE, plotStationary=TRUE)

write.csv(data,"crawl_data60.csv", row.names = FALSE)
###############################
am308<- data[104939:117390,] #varies, some outliers at the beginning
am308$step[is.na(am308$step)] <- 200
am308$angle[is.na(am308$angle)] <- 0
whichzero <- which(am308$step == 0)
# Proportion of steps of length zero in the data set
length(whichzero)/nrow(am308)

any(is.na(am308$temp))
nbStates <- 3
stepDist <- "gamma" # step distribution
angleDist <- "vm" # turning angle distribution

am308$temp[is.na(am308$temp)] <- 25
## initial parameters for gamma and von Mises distributions
mu0 <- c(21.22671,136.8682,483.87) # step mean 500 21.22671,136.8682,483.87
sigma0 <- c(17.90707,84.60687,289.3098) # step SD 17.90707,84.60687,289.3098
zeromass0 <- c(0.0001606168,0.0001606168,0.0001606168) #0.0002280502,0.0002280502,0.0002280502
stepPar0 <- c(mu0,sigma0,zeromass0)
angleMean0 <- c(0.003592673,0.01941996,0.01499825) # angle mean
kappa0 <- c(0.355384957,1.53831499,2.59417278) # angle concentration
anglePar0 <- c(angleMean0)
## call to fitting function
stateNames <- c("Resting","Foraging", "Traveling")
m <- fitHMM(data=am308, nbStates=nbStates,dist=list(step=stepDist,angle=angleDist), Par0=list(step=stepPar0,angle=anglePar0),formula=~temp, stateNames=stateNames)
m
CIreal(m)

plot(m, plotCI=TRUE, plotStationary=TRUE)
##############################
am93<- data[117393:134931,] #varies, some outliers at the beginning
am93$step[is.na(am93$step)] <- 200
am93$angle[is.na(am93$angle)] <- 0
whichzero <- which(am93$step == 0)
# Proportion of steps of length zero in the data set
length(whichzero)/nrow(am93)

any(is.na(am93$temp))
nbStates <- 3
stepDist <- "gamma" # step distribution
angleDist <- "vm" # turning angle distribution

am93$temp[is.na(am93$temp)] <- 25
## initial parameters for gamma and von Mises distributions
mu0 <- c(30,270,700) # step mean 500 21.22671,136.8682,483.87
sigma0 <- c(30,270,700) # step SD 17.90707,84.60687,289.3098
zeromass0 <- c(5.701579e-05,5.701579e-05,5.701579e-05) #0.0002280502,0.0002280502,0.0002280502
stepPar0 <- c(mu0,sigma0,zeromass0)
angleMean0 <- c(0.1,0.1,0.1) # angle mean
kappa0 <- c(1,2,3) # angle concentration
anglePar0 <- c(angleMean0)
## call to fitting function
stateNames <- c("Resting","Foraging", "Traveling")
m <- fitHMM(data=am93, nbStates=nbStates,dist=list(step=stepDist,angle=angleDist), Par0=list(step=stepPar0,angle=anglePar0),formula=~temp, stateNames=stateNames)
m
CIreal(m)

plot(m, breaks=20,plotCI=TRUE, plotStationary=TRUE)
