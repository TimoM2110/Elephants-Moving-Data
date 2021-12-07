

# This is a script for an Ornstein-Uhlenbeck model in order to estimate the
# velocity of the elephants depending on the temperature

# install package dependencies
library(devtools)
install_github("TheoMichelot/smoothSDE")
library(smoothSDE)
library(lubridate)
library(sp)
library(sf)
library(tidyverse)
# set working directory
setwd("C:/Users/timom")

data <- read.table("crawl_data60.csv", header = TRUE, sep= ",")
data <- data[ -c(3) ]

datas <- read.table("crawl_data6002.csv", header = TRUE, sep= ",")
datas <- datas[ -c(3) ]

AM105 <- data[14500:15540,] #
AM107 <- data[17544:20381,] #
AM108 <- datas[40984:42153,] #
AM110 <- data[36084:38623,] #
AM239 <- data[52625:70164,] #
AM253<- data[71517:81805,] #
AM254 <- datas[101480:106018,] #
AM255<- data[85724:91271,] #
AM306 <- read.table("C:/Users/timom/anaconda3/envs/madeleine_project/crawl_data_am306.csv", header = TRUE, sep= ",")
AM306$numer <- 1
AM306$visible <- 'True'
AM306 <- AM306[, c("ID", 'step',"TimeNum","locType", "numer", "time", "visible",
                   "temp", "nu.x", "nu.y", "se.mu.x","se.nu.x","se.mu.y","se.nu.y", "speed", "x", "y")]
AM306 <- AM306[3252:5027,] ###
AM307 <- data[92834:99336,] #
AM308 <- data[110639:114390,] #
AM91 <- read.table("C:/Users/timom/anaconda3/envs/madeleine_project/crawl_data_am91.csv", header = TRUE, sep= ",") ###
AM91$numer <- 1
AM91$visible <- 'True'
AM91 <- AM91[, c("ID", 'step',"TimeNum","locType", "numer", "time", "visible",
                 "temp", "nu.x", "nu.y", "se.mu.x","se.nu.x","se.mu.y","se.nu.y", "speed", "x", "y")]
AM93 <- read.table("C:/Users/timom/anaconda3/envs/madeleine_project/crawl_data_am93.csv", header = TRUE, sep= ",")
AM93$numer <- 1
AM93$visible <- 'True'
AM93 <- AM93[, c("ID", 'step',"TimeNum","locType", "numer", "time", "visible",
                 "temp", "nu.x", "nu.y", "se.mu.x","se.nu.x","se.mu.y","se.nu.y", "speed", "x", "y")]
AM93 <- AM93[3307:11615,] ###
AM99 <- data[134933:152473,] #

data <- rbind(AM105,AM107,AM108,AM110,AM239,AM253,AM254,AM255,AM306,AM307,AM308,AM91,AM93,AM99)

data <- data[-c(3:5,7,9:15)]

data$time <- as.POSIXct(
  strptime(
    data$time,
    format='%Y-%m-%d %H:%M:%S'
  )
)

data <- na.omit(data)

data <- data[data$step >= 1.5, ]  
data <- data[data$step <= 4000, ]

#

whichzero <- which(data$step == 0)
# Proportion of steps of length zero in the data set
length(whichzero)/nrow(data)

data$step[data$step==0] <- 0.001

data$temp[is.na(data$step)] <- 200

data$step[data$step>=5000] <- runif(1, min=4000, max = 4999)

#####################################################
# Fitting a Continuous time correlated random walk ##
#####################################################
names(data)[names(data) == "time"] <- "date"
data$date <- as.POSIXlt(data$date, tz = "GMT")
# First period in order to avoid seasonal effects
#####################################################
# Einschub ##
#####################################################
track105 <- subset(
  data,
  ID == unique(ID)[1] # keep AM99 only
) 
track105$ID <- 1
#track105$time <- as.numeric(track105$date - min(track105$date))/3600

track107 <- subset(
  data,
  ID == unique(ID)[2] # keep AM99 only
) 
track107$ID <- 2
#track107$time <- as.numeric(track107$date - min(track107$date))/3600

track108 <- subset(
  data,
  ID == unique(ID)[3] # keep AM99 only
) 
track108$ID <- 3
#track108$time <- as.numeric(track108$date - min(track108$date))/3600

track110 <- subset(
  data,
  ID == unique(ID)[4] # keep AM99 only
) 
track110$ID <- 4
#track110$time <- as.numeric(track110$date - min(track110$date))/3600

track239 <- subset(
  data,
  ID == unique(ID)[5] # keep AM99 only
) 
track239$ID <- 5
#track239$time <- as.numeric(track239$date - min(track239$date))/3600
track239 <- track239[(8600:9200),]

track253 <- subset(
  data,
  ID == unique(ID)[6] # keep AM99 only
) 
track253$ID <- 6
#track253$time <- as.numeric(track253$date - min(track253$date))/3600
track253 <- track253[(6800:7600),]

track254 <- subset(
  data,
  ID == unique(ID)[7] # keep AM99 only
) 
track254$ID <- 7
#track254$time <- as.numeric(track254$date - min(track254$date))/3600
#track254 <- track254[(3200:4200),]

track255 <- subset(
  data,
  ID == unique(ID)[8] # keep AM99 only
) 
track255$ID <- 8
#track255$time <- as.numeric(track255$date - min(track255$date))/3600
#track255 <- track255[(3500:4000),]

track306 <- subset(
  data,
  ID == unique(ID)[9] # keep AM99 only
) 
track306$ID <- 9
#track306$time <- as.numeric(track306$date - min(track306$date))/3600

track307 <- subset(
  data,
  ID == unique(ID)[10] # keep AM99 only
) 
track307$ID <- 10

track307 <- track307[(1:2687),]
#track307$time <- as.numeric(track307$date - min(track307$date))/3600

track308 <- subset(
  data,
  ID == unique(ID)[11] # keep AM99 only
) 
track308$ID <- 11
#track308$time <- as.numeric(track308$date - min(track308$date))/3600

track91 <- subset(
  data,
  ID == unique(ID)[12] # keep AM99 only
) 
track91$ID <- 12
#track91$time <- as.numeric(track91$date - min(track91$date))/3600

track93 <- subset(
  data,
  ID == unique(ID)[13] # keep AM99 only
) 
track93$ID <- 13
#track93$time <- as.numeric(track93$date - min(track93$date))/3600

track99 <- subset(
  data,
  ID == unique(ID)[14] # keep AM99 only
) 
track99$ID <- 14
#track99$time <- as.numeric(track99$date - min(track99$date))/3600
track99 <- track99[(1:11199),]

data <- rbind(track105,track107,track108,track110,track254,track255,track306,track307,track308,track91,track93,track99)#,track105,track107,track108,track110,track239,track253,track254,track255,track306,track307,track308,track91,track93,track99
#####################################################
#data$time <- as.numeric(data$date - min(data$date))/3600
# Convert projected coordinates (UTM) to km
#library(dplyr)
#data <- data %>% slice(-c(11199))

#data <- data[(1:73128),]


#data$time <- as.numeric(data$date - min(data$date))/3600
#data$time <- as.numeric(data$date)
data$time <- seq(1,62860)
#data <- data[order(data$time), ]

data <- data   # for correct index
rownames(data) <- 1:nrow(data)    
                               
boxplot(data$time,
        ylab = "y"
)


data$x <- data$x/1000
data$y <- data$y/1000

# Data set including ID, time, responses (x, y), and covariate (temp)
head(data)

# Using the Ornstein-Uhlenbeck model (also called Continuous Time
# Correlated Random Walk (CTCRW))
# Model formulas for CTCRW parameters
# 2 parameters, beta for mean reversion and sigma for variance
formulas <- list(
  beta = ~ s(temp, k = 10, bs = "ts"),
  sigma = ~ s(temp, k = 10, bs = "ts")
)

# Type of stochastic differential equation (SDE) model: continuous-time
# correlated random walk
type <- "CTCRW"

# Create SDE model object (inheriting from SDE object)
my_sde <- SDE$new(
  formulas = formulas,
  data = data,
  type = type,
  response = c("x", "y"),
  par0 = c(7,3)
)

# Fit model 
my_sde$fit()

# Plot the model parameters (beta, sigma)
my_sde$plot_par(
  "temp",
  n_post = 100
)

######################
## Plot the results ##
######################

# design matrices for the plot
mats <- my_sde$make_mat_grid(var = "temp")

# Save fitted Parameter estimates
par <- my_sde$par(
  t = "all",
  X_fe = mats$X_fe,
  X_re = mats$X_re
)

# Posterior samples
n_post <- 2000

post <- my_sde$post_par(
  X_fe = mats$X_fe,
  X_re = mats$X_re,
  n_post = n_post
)

# Get nu parameter (mean speed, y-axis of the plot)
nu <- sqrt(pi) * par[,"sigma"] / (2 * sqrt(par[,"beta"]))
post_nu <- sqrt(pi) * post[,"sigma",] / (2 * sqrt(post[,"beta",]))

# Confidence intervals
quants <- apply(
  post_nu,
  1,
  quantile,
  probs = c(0.025, 0.975)
)

new_data <- data.frame(
  temp = seq(
    from = min(data$temp),
    to = max(data$temp),
    length.out = 1000
  )
)

df <- data.frame(
  x = new_data$temp,
  mle = nu,
  low = quants[1,],
  upp = quants[2,]
)

# Plot
p_est <- ggplot(df, aes(x, mle)) +
  xlab("temperature (°C)") + ylab(bquote(nu[t])) + geom_line(size = 1) +
  geom_ribbon(aes(ymin = low, ymax = upp), alpha = 0.3) +
  theme_light() +
  theme(axis.title = element_text(size = 12),
        strip.text = element_text(color = 1, size = 12))

p_est


###################
###################
track105 <- subset(
  data,
  ID == unique(ID)[1] # keep AM99 only
) 
track105$time <- as.numeric(track105$date - min(track105$date))/3600
track105$ID <- 1
keep_rows <- which(
  track105$date > as.POSIXct("2009-04-08 01:30:00") &
    track105$date <= as.POSIXct("2009-05-21 11:00:00")
)
track105 <- track105[keep_rows,]

track107 <- subset(
  data,
  ID == unique(ID)[2] # keep AM99 only
) 
track107$time <- as.numeric(track107$date - min(track107$date))/3600
track107$ID <- 2
keep_rows <- which(
  track107$date > as.POSIXct("2007-08-12 22:30:00") &
    track107$date <= as.POSIXct("2007-12-09 05:00:00")
)
track107 <- track107[keep_rows,]

track108 <- subset(
  data,
  ID == unique(ID)[3] # keep AM99 only
) 
track108$time <- as.numeric(track108$date - min(track108$date))/3600
track108$ID <- 3
keep_rows <- which(
  track108$date > as.POSIXct("2008-04-15 01:00:00") &
    track108$date <= as.POSIXct("2008-06-02 19:00:00")
)
track108 <- track108[keep_rows,]

track110 <- subset(
  data,
  ID == unique(ID)[4] # keep AM99 only
) 
track110$time <- as.numeric(track110$date - min(track110$date))/3600
track110$ID <- 4
keep_rows <- which(
  track110$date > as.POSIXct("2007-09-23 15:30:00") &
    track110$date <= as.POSIXct("2008-01-07 11:30:00")
)
track110 <- track110[keep_rows,]

track239 <- subset(
  data,
  ID == unique(ID)[5] # keep AM99 only
) 
track239$time <- as.numeric(track239$date - min(track239$date))/3600
track239$ID <- 5
keep_rows <- which(
  track239$date > as.POSIXct("2007-08-12 22:30:00") &
    track239$date <= as.POSIXct("2009-08-12 17:30:00")
)
track239 <- track239[keep_rows,]

track253 <- subset(
  data,
  ID == unique(ID)[6] # keep AM99 only
) 
track253$time <- as.numeric(track253$date - min(track253$date))/3600
track253$ID <- 6
keep_rows <- which(
  track253$date > as.POSIXct("2007-10-08 05:30:00") &
    track253$date <= as.POSIXct("2008-12-09 21:30:00")
)
track253 <- track253[keep_rows,]

track254 <- subset(
  data,
  ID == unique(ID)[7] # keep AM99 only
) 
track254$time <- as.numeric(track254$date - min(track254$date))/3600
track254$ID <- 7
keep_rows <- which(
  track254$date > as.POSIXct("2008-04-18 23:30:00") &
    track254$date <= as.POSIXct("2008-10-25 01:30:00")
)
track254 <- track254[keep_rows,]

track255 <- subset(
  data,
  ID == unique(ID)[8] # keep AM99 only
) 
track255$time <- as.numeric(track255$date - min(track255$date))/3600
track255$ID <- 8
keep_rows <- which(
  track255$date > as.POSIXct("2008-01-23 03:30:00") &
    track255$date <= as.POSIXct("2008-09-10 06:30:00")
)
track255 <- track255[keep_rows,]

track306 <- subset(
  data,
  ID == unique(ID)[9] # keep AM99 only
) 
track306$time <- as.numeric(track306$date - min(track306$date))/3600
track306$ID <- 9
keep_rows <- which(
  track306$date > as.POSIXct("2008-07-16 11:01:00") &
    track306$date <= as.POSIXct("2008-09-28 10:01:00")
)
track306 <- track306[keep_rows,]

track307 <- subset(
  data,
  ID == unique(ID)[10] # keep AM99 only
) 
track307$time <- as.numeric(track307$date - min(track307$date))/3600
track307$ID <- 10
keep_rows <- which(
  track307$date > as.POSIXct("2008-03-26 12:00:00") &
    track307$date <= as.POSIXct("2008-12-22 10:00:00")
)
track307 <- track307[keep_rows,]

track308 <- subset(
  data,
  ID == unique(ID)[11] # keep AM99 only
) 
track308$time <- as.numeric(track308$date - min(track308$date))/3600
track308$ID <- 11
keep_rows <- which(
  track308$date > as.POSIXct("2008-11-03 12:00:00") &
    track308$date <= as.POSIXct("2009-04-08 19:00:00")
)
track308 <- track308[keep_rows,]

track91 <- subset(
  data,
  ID == unique(ID)[12] # keep AM99 only
) 
track91$time <- as.numeric(track91$date - min(track91$date))/3600
track91$ID <- 12
keep_rows <- which(
  track91$date > as.POSIXct("2007-08-12 22:33:00") &
    track91$date <= as.POSIXct("2009-08-12 08:33:00")
)
track91 <- track91[keep_rows,]

track93 <- subset(
  data,
  ID == unique(ID)[13] # keep AM99 only
) 
track93$time <- as.numeric(track93$date - min(track93$date))/3600
track93$ID <- 13
keep_rows <- which(
  track93$date > as.POSIXct("2007-12-28 17:30:00") &
    track93$date <= as.POSIXct("2008-12-08 21:30:00")
)
track93 <- track93[keep_rows,]

track99 <- subset(
  data,
  ID == unique(ID)[14] # keep AM99 only
) 
track99$time <- as.numeric(track99$date - min(track99$date))/3600
track99$ID <- 14
keep_rows <- which(
  track99$date > as.POSIXct("2007-08-12 22:30:00") &
    track99$date <= as.POSIXct("2009-08-12 18:30:00")
)
track99 <- track99[keep_rows,]



######
#latest
######
track105 <- subset(
  data,
  ID == unique(ID)[1] # keep AM99 only
) 
track105$time <- as.numeric(track105$date - min(track105$date))/3600
track105$ID <- 1
keep_rows <- which(
  track105$date > as.POSIXct("2008-05-01 00:00:00") &
    track105$date <= as.POSIXct("2008-09-30 23:59:59")
)
track105 <- track105[keep_rows,]

track107 <- subset(
  data,
  ID == unique(ID)[2] # keep AM99 only
) 
track107$time <- as.numeric(track107$date - min(track107$date))/3600
track107$ID <- 2
keep_rows <- which(
  track107$date > as.POSIXct("2009-06-11 00:00:00") &
    track107$date <= as.POSIXct("2009-08-12 23:59:59")
)
track107 <- track107[keep_rows,]

track108 <- subset(
  data,
  ID == unique(ID)[3] # keep AM99 only
) 
track108$time <- as.numeric(track108$date - min(track108$date))/3600
track108$ID <- 3
keep_rows <- which(
  track108$date > as.POSIXct("2008-05-01 00:00:00") &
    track108$date <= as.POSIXct("2008-09-30 23:59:59")
)
track108 <- track108[keep_rows,]

track110 <- subset(
  data,
  ID == unique(ID)[4] # keep AM99 only
) 
track110$time <- as.numeric(track110$date - min(track110$date))/3600
track110$ID <- 4
keep_rows <- which(
  track110$date > as.POSIXct("2008-05-01 00:00:00") &
    track110$date <= as.POSIXct("2008-09-30 23:59:59")
)
track110 <- track110[keep_rows,]

track239 <- subset(
  data,
  ID == unique(ID)[5] # keep AM99 only
) 
track239$time <- as.numeric(track239$date - min(track239$date))/3600
track239$ID <- 5
keep_rows <- which(
  track239$date > as.POSIXct("2007-08-12 22:30:00") &
    track239$date <= as.POSIXct("2009-08-12 17:30:00")
)
track239 <- track239[keep_rows,]

track253 <- subset(
  data,
  ID == unique(ID)[6] # keep AM99 only
) 
track253$time <- as.numeric(track253$date - min(track253$date))/3600
track253$ID <- 6
keep_rows <- which(
  track253$date > as.POSIXct("2007-10-08 05:30:00") &
    track253$date <= as.POSIXct("2008-12-09 21:30:00")
)
track253 <- track253[keep_rows,]

track254 <- subset(
  data,
  ID == unique(ID)[7] # keep AM99 only
) 
track254$time <- as.numeric(track254$date - min(track254$date))/3600
track254$ID <- 7
keep_rows <- which(
  track254$date > as.POSIXct("2008-04-18 23:30:00") &
    track254$date <= as.POSIXct("2008-10-25 01:30:00")
)
track254 <- track254[keep_rows,]

track255 <- subset(
  data,
  ID == unique(ID)[8] # keep AM99 only
) 
track255$time <- as.numeric(track255$date - min(track255$date))/3600
track255$ID <- 8
keep_rows <- which(
  track255$date > as.POSIXct("2008-01-23 03:30:00") &
    track255$date <= as.POSIXct("2008-09-10 06:30:00")
)
track255 <- track255[keep_rows,]

track306 <- subset(
  data,
  ID == unique(ID)[9] # keep AM99 only
) 
track306$time <- as.numeric(track306$date - min(track306$date))/3600
track306$ID <- 9
keep_rows <- which(
  track306$date > as.POSIXct("2008-07-21 11:01:00") &
    track306$date <= as.POSIXct("2008-09-18 10:01:00")
)
track306 <- track306[keep_rows,]

track307 <- subset(
  data,
  ID == unique(ID)[10] # keep AM99 only
) 
track307$time <- as.numeric(track307$date - min(track307$date))/3600
track307$ID <- 10
keep_rows <- which(
  track307$date > as.POSIXct("2008-03-26 12:00:00") &
    track307$date <= as.POSIXct("2008-12-22 10:00:00")
)
track307 <- track307[keep_rows,]

track308 <- subset(
  data,
  ID == unique(ID)[11] # keep AM99 only
) 
track308$time <- as.numeric(track308$date - min(track308$date))/3600
track308$ID <- 11
keep_rows <- which(
  track308$date > as.POSIXct("2008-11-03 12:00:00") &
    track308$date <= as.POSIXct("2009-04-08 19:00:00")
)
track308 <- track308[keep_rows,]

track91 <- subset(
  data,
  ID == unique(ID)[12] # keep AM99 only
) 
track91$time <- as.numeric(track91$date - min(track91$date))/3600
track91$ID <- 12
keep_rows <- which(
  track91$date > as.POSIXct("2007-08-12 22:33:00") &
    track91$date <= as.POSIXct("2009-08-12 08:33:00")
)
track91 <- track91[keep_rows,]

track93 <- subset(
  data,
  ID == unique(ID)[13] # keep AM99 only
) 
track93$time <- as.numeric(track93$date - min(track93$date))/3600
track93$ID <- 13
keep_rows <- which(
  track93$date > as.POSIXct("2007-12-28 17:30:00") &
    track93$date <= as.POSIXct("2008-12-08 21:30:00")
)
track93 <- track93[keep_rows,]

track99 <- subset(
  data,
  ID == unique(ID)[14] # keep AM99 only
) 
track99$time <- as.numeric(track99$date - min(track99$date))/3600
track99$ID <- 14
keep_rows <- which(
  track99$date > as.POSIXct("2007-08-12 22:30:00") &
    track99$date <= as.POSIXct("2009-08-12 18:30:00")
)
track99 <- track99[keep_rows,]

track105_a <- track105[(1:500),]
track105_b <- track105[(500:1000),]
track105_b$ID <- 2

track93_a <- track93[(1:3567),]
track93_b <- track93[(6784:9970),]
track93_b$ID <- 3







# read in the data
data <- read.table(
  "ThermochronTracking Elephants Kruger 2007.csv",
  header = TRUE,
  sep= ","
)

#rename some columns for readability
names(data)[names(data) == "location.long"] <- "x"
names(data)[names(data) == "location.lat"] <- "y"
names(data)[names(data) == "timestamp"] <- "time"
names(data)[names(data) == "individual.local.identifier"] <- "ID"
names(data)[names(data) == "external.temperature"] <- "temp"

# convert time to datetime
data$time <- as.POSIXct(
  strptime(
    data$time,
    format='%Y-%m-%d %H:%M:%S'
  )
)

# delete entries containing NAs
data <- na.omit(data)

# omit irrelevant columns
data <- data[ -c(1, 2,7:9,11) ]

######################################
# Converting the coordinates to UTM ##
######################################
# make a copy of the prepared data for conversion
datas<-data.frame(data)

# project the data
coordinates(datas) = ~x +y
proj4string(datas) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +units=m")
datas <- spTransform(
  datas,
  CRS="+proj=utm +zone=36J +datum=WGS84 "
)
datas <- as.data.frame(coordinates(datas))

# pass converted coordinates to the data frame
data['x'] = datas['x']
data['y'] = datas['y'] 
