library(crawl)
library(tidyverse)
library(lubridate)
library(sp)
library(sf)
library(momentuHMM)

setwd("C:/Users/timom/anaconda3/envs/madeleine_project")
data <- read.table("ThermochronTracking Elephants Kruger 2007.csv", header = TRUE, sep= ",")

data$timestamp <- as.POSIXct(strptime(data$timestamp, format='%Y-%m-%d %H:%M:%S'))

indx <- apply(data, 2, function(x) any(is.na(x) | is.infinite(x)))
colnames[indx]

data <- na.omit(data)

names(data)[names(data) == "timestamp"] <- "time"
df <- data[ -c(1,2, 7:9, 11) ]

############################################################################
datas <- read.table("ThermochronTracking Elephants Kruger 2007.csv", header = TRUE, sep= ",")

datas$timestamp <- as.POSIXct(strptime(datas$timestamp, format='%Y-%m-%d %H:%M:%S'))

indx <- apply(datas, 2, function(x) any(is.na(x) | is.infinite(x)))
colnames[indx]

datas <- na.omit(datas)

indx <- apply(datas, 2, function(x) any(is.na(x) | is.infinite(x)))
colnames[indx]
indx


#datas = subset(datas, select=-c('event.id', 'visible', 'external.temperature',
# 'sensor.type', 'individual.taxon.canonical.name', 'tag.local.identifier', 'individual.local.identifier', 'study.name'))
#
#coordinates(datas) <- c('location.long', "location.lat")

datas$timestamp <- as.POSIXct(strptime(datas$timestamp, format='%Y-%m-%d %H:%M:%S'))

names(datas)[names(datas) == "location.long"] <- "longitude"
names(datas)[names(datas) == "location.lat"] <- "latitude"

coordinates(datas) = ~longitude +latitude

proj4string(datas) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +units=m")

datas <- spTransform(datas,CRS="+proj=utm +zone=36J +datum=WGS84 ")

datas <- as.data.frame(coordinates(datas))
###############################################################################

df['location.long'] = datas['longitude'] #erst untitled 8 laufen lassen
df['location.lat'] = datas['latitude'] # dort wurde longlat in utm umgewandelt

#library(sp)
#library(rgdal)


names(df)[names(df) == "individual.local.identifier"] <- "ID"
names(df)[names(df) == "location.long"] <- "x"
names(df)[names(df) == "location.lat"] <- "y"

#df <- df[ -c(4) ]
#am99 <- df[252240:283688,]

######################
# Parameters of the correlated random walk (for simulation)
#beta <- 0.1
#sigma <- 10

# vorher am99 bei kalman inits
#kalman_inits <- list(a = c(df$x[1], 0, df$y[1], 0),
#                     P = diag(c(0, 100, 0, 100)))

# not working with kalman_inits (initial.state)
crwOut <- crawlWrap(df, timeStep = "1 hour", theta=c(2, -2), fixPar=c(NA,NA)) #theta=c(6.855, -0.007)
# noch gucken ob das passt und long lat koordinaten probieren

ggplot(df, aes(x, y)) + geom_point() + geom_path() +
  geom_point(aes(mu.x, mu.y), crwOut$crwPredict, col = "firebrick", size = 0.5) +
  geom_path(aes(mu.x, mu.y), crwOut$crwPredict, col = "firebrick", size = 0.3)

data <- prepData(crwOut,type="UTM",coordNames=c("x","y"))

summary(data)

plot(data,compact=T)

indx <- apply(data, 2, function(x) any(is.na(x) | is.infinite(x)))
indx # No NAs in x, y and time column

# klappt
##############################################################################
