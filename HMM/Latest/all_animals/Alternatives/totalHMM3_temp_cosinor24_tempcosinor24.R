library(momentuHMM)
setwd("C:/Users/timom")
data <- read.table("crawl_data60.csv", header = TRUE, sep= ",")

data$time <- as.POSIXct(strptime(data$time, format='%Y-%m-%d %H:%M:%S'))
data <- data[ -c(2,3) ]

data <- prepData(data,type="UTM",coordNames=c("x","y"))

datas <- read.table("crawl_data6002.csv", header = TRUE, sep= ",")
datas$time <- as.POSIXct(strptime(datas$time, format='%Y-%m-%d %H:%M:%S'))
datas <- datas[ -c(2,3) ]
datas <- prepData(datas,type="UTM",coordNames=c("x","y"))

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
AM306 <- AM306[, c("ID", "step", "angle", "TimeNum","locType", "numer", "time", "visible",
                   "temp", "nu.x", "nu.y", "se.mu.x","se.nu.x","se.mu.y","se.nu.y", "speed", "x", "y")]
AM306 <- AM306[3252:5027,] ###
AM307 <- data[92834:99336,] #
AM308 <- data[110639:114390,] #
AM91 <- read.table("C:/Users/timom/anaconda3/envs/madeleine_project/crawl_data_am91.csv", header = TRUE, sep= ",") ###
AM91$numer <- 1
AM91$visible <- 'True'
AM91 <- AM91[, c("ID", "step", "angle", "TimeNum","locType", "numer", "time", "visible",
                 "temp", "nu.x", "nu.y", "se.mu.x","se.nu.x","se.mu.y","se.nu.y", "speed", "x", "y")]
AM93 <- read.table("C:/Users/timom/anaconda3/envs/madeleine_project/crawl_data_am93.csv", header = TRUE, sep= ",")
AM93$numer <- 1
AM93$visible <- 'True'
AM93 <- AM93[, c("ID", "step", "angle", "TimeNum","locType", "numer", "time", "visible",
                 "temp", "nu.x", "nu.y", "se.mu.x","se.nu.x","se.mu.y","se.nu.y", "speed", "x", "y")]
AM93 <- AM93[3307:11615,] ###
AM99 <- data[134933:152473,] #

data <- rbind(AM105,AM107,AM108,AM110,AM239,AM253,AM254,AM255,AM306,AM307,AM308,AM91,AM93,AM99)

library(lubridate)

data$month <- month(as.POSIXlt(data$time, format="%Y/%m/%d %H:%M:%S"))
data$hour <- hour(as.POSIXlt(data$time, format="%Y/%m/%d %H:%M:%S"))

plot(data,compact=T)

# Indices of steps of length zero
whichzero <- which(data$step == 0)
# Proportion of steps of length zero in the data set
length(whichzero)/nrow(data)

nbStates <- 3
stepDist <- "gamma" # step distribution
angleDist <- "vm" # turning angle distribution
data$temp[is.na(data$step)] <- 200
data$temp[is.na(data$temp)] <- 25
data$angle[is.na(data$angle)] <- 0

any(is.na(data$temp))

## initial parameters for gamma and von Mises distributions
mu0 <- c(8.705417e+01,2.387662e+02,6.845123e+02) # 
sigma0 <- c(8.754230e+01,1.357714e+02,4.199985e+02) # 
zeromass0 <- c(2.088501e-04,3.987702e-04,2.754461e-10) #
stepPar0 <- c(mu0,sigma0, zeromass0)
angleMean0 <- c(0.003592673,0.01941996,0.01499825) # angle mean
kappa0 <- c(0.355384957,1.53831499,3.59417278) # angle concentration 0.355384957,1.53831499,3.59417278
anglePar0 <- c(angleMean0)
## call to fitting function
stateNames <- c("Resting","Foraging", "Traveling")
formula = ~temp+cosinor(hour, 24)+ temp*cosinor(hour, 24)
m <- fitHMM(data=data, nbStates=nbStates,dist=list(step=stepDist,angle=angleDist), Par0=list(step=stepPar0,angle=anglePar0),formula=formula, stateNames=stateNames) #
m
CIreal(m)

#m <- readRDS("totalHMM3.rds")
plot(m, plotCI=TRUE, breaks = 40, plotStationary=TRUE,lwd=0.5)

plotPR(m, ncores=7)
AIC(m)

saveRDS(m, "totalHMM3_temp_cosinor24_tempcosinor24.rds")

#####
mod1 <- readRDS("totalHMM2_states.rds")
plot(mod1, plotCI=TRUE, breaks = 40, plotStationary=TRUE,lwd=2)
plotPR(mod1, ncores=7)
mod2 <- readRDS("totalHMM4.rds")
plot(mod2, plotCI=TRUE, breaks = 40, plotStationary=TRUE,lwd=0.5, plotTracks = F)
aic <- numeric(3)
aic[1] <- AIC(mod1)
aic[2] <- AIC(mod2)
aic[3] <- AIC(m)

### BIC
bic <- numeric(3)
T <- dim(data)[1]
bic[1] <- 2*mod1$mod$minimum + log(T)*length(mod1$mod$estimate)
bic[2] <- 2*mod2$mod$minimum + log(T)*length(mod2$mod$estimate)
bic[3] <- 2*m$mod$minimum + log(T)*length(m$mod$estimate)
aic
bic

m <- readRDS("totalHMM3_temp+cosinor24.rds")

quantile(data$temp)

temp_frame = as.data.frame(t(temps));
frame_all = data.frame(t(hours), t(temps))
temps <- temp_frame[, 1]
my_name_vector      = c(rep("hour", 24), "temp")
colnames(frame_all) <- my_name_vector

plotStationary(m)

names <- rep("temp", 131)

setwd("C:/Users/timom/gif_hmm3")

pdf("gif_hmm3.pdf",width=6,height=4);par(mfrow=c(1,1),mar=c(4,4,2,1))
for (zoo in 1:411){
  tt<-5.9+zoo*0.1
  hour <- seq(0,23,1)
  frame_allz = data.frame(t(hour), t(tt))
  my_name_vector      = c(rep("hour", 24), "temp")
  colnames(frame_allz) <- my_name_vector
  plotStationary(m, covs= frame_allz, plotCI = T,lwd=0.5)
}
graphics.off()

######################
## Global Decoding ###
######################
data$locType <- seq.int(nrow(data))
vitstates <- viterbi(mod5)
IDs <- data$ID[!duplicated(data$ID)] 
colors <- c("#E69F00", "#56B4E9", "#009E73")
viterbis <- viterbi(m)
names <- c("AM105","AM107","AM108","AM110","AM239","AM253","AM254","AM255",
           "AM306","AM307","AM308","AM91","AM93","AM99")

sum <- 0
for (i in 1:14){
  a <- sum + 1
  track <- subset(
    data,
    ID == unique(ID)[i] 
  )
  sum <- sum + nrow(track)
  states <- vitstates[a:(sum-1)]
  plot(track$time,track$step,type="h",xlab="time",
       ylab="step length in metres",main=paste(names[i],"decoded states"),
       col=colors[states])
  legend("topright", legend=c("Resting","Foraging", "Traveling"), 
         lwd=2, col = colors, cex = 0.55)
}

