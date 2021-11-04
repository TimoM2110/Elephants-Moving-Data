library(crawl)
library(tidyverse)
library(lubridate)
library(sp)
library(sf)
library(momentuHMM)

setwd("C:/Users/timom/anaconda3/envs/madeleine_project")
data <- read.table("ThermochronTracking Elephants Kruger 2007.csv", header = TRUE, sep= ",")

names(data)[names(data) == "location.long"] <- "x"
names(data)[names(data) == "location.lat"] <- "y"
names(data)[names(data) == "timestamp"] <- "time"
names(data)[names(data) == "individual.local.identifier"] <- "ID"
names(data)[names(data) == "external.temperature"] <- "temp"

data$time <- as.POSIXct(strptime(data$time, format='%Y-%m-%d %H:%M:%S'))

indx <- apply(data, 2, function(x) any(is.na(x) | is.infinite(x)))
colnames[indx]

data <- na.omit(data)


df <- data[ -c(1, 2,7:9,11) ]

############################################################################
datas <- read.table("ThermochronTracking Elephants Kruger 2007.csv", header = TRUE, sep= ",")

names(datas)[names(datas) == "location.long"] <- "x"
names(datas)[names(datas) == "location.lat"] <- "y"
names(datas)[names(datas) == "Unnamed..0"] <- "time"
names(datas)[names(datas) == "individual.local.identifier"] <- "ID"
names(datas)[names(datas) == "external.temperature"] <- "temp"

datas$time <- as.POSIXct(strptime(datas$time, format='%Y-%m-%d %H:%M:%S'))

indx <- apply(datas, 2, function(x) any(is.na(x) | is.infinite(x)))
indx

datas <- na.omit(datas)

coordinates(datas) = ~x +y

proj4string(datas) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +units=m")

datas <- spTransform(datas,CRS="+proj=utm +zone=36J +datum=WGS84 ")

datas <- as.data.frame(coordinates(datas))
###############################################################################
# project long lat
df['x'] = datas['x']
df['y'] = datas['y'] 

# not working with kalman_inits (initial.state)
crwOut <- crawlWrap(df, timeStep = "1 hour", theta=c(2, -2), fixPar=c(NA,NA)) #theta=c(6.855, -0.007)

ggplot(df, aes(x, y)) + geom_point() + geom_path() +
  geom_point(aes(mu.x, mu.y), crwOut$crwPredict, col = "firebrick", size = 0.5) +
  geom_path(aes(mu.x, mu.y), crwOut$crwPredict, col = "firebrick", size = 0.3)

data <- prepData(data,type="UTM",coordNames=c("x","y")) #crwOut if imputed in this session

summary(data)

plot(data,compact=T)