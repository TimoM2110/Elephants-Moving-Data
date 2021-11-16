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

am105 <- data[14500:15540,] #
am107 <- data[17544:20381,] #
am108 <- datas[40984:42153,] #
am110 <- data[36084:38623,] #
am239 <- data[52625:70164,] #
am253<- data[71517:81805,] #
am254 <- datas[101480:106018,] #
am255<- data[85724:91271,] #
am306 <- read.table("C:/Users/timom/anaconda3/envs/madeleine_project/crawl_data_am306.csv", header = TRUE, sep= ",")
am306$numer <- 1
am306$visible <- 'True'
am306 <- am306[, c("ID", "step", "angle", "TimeNum","locType", "numer", "time", "visible",
             "temp", "nu.x", "nu.y", "se.mu.x","se.nu.x","se.mu.y","se.nu.y", "speed", "x", "y")]
am306 <- am306[3252:5027,] ###
am307 <- data[92834:99336,] #
am308 <- data[110639:114390,] #
am91 <- read.table("C:/Users/timom/anaconda3/envs/madeleine_project/crawl_data_am91.csv", header = TRUE, sep= ",") ###
am91$numer <- 1
am91$visible <- 'True'
am91 <- am91[, c("ID", "step", "angle", "TimeNum","locType", "numer", "time", "visible",
                   "temp", "nu.x", "nu.y", "se.mu.x","se.nu.x","se.mu.y","se.nu.y", "speed", "x", "y")]
am93 <- read.table("C:/Users/timom/anaconda3/envs/madeleine_project/crawl_data_am93.csv", header = TRUE, sep= ",")
am93$numer <- 1
am93$visible <- 'True'
am93 <- am93[, c("ID", "step", "angle", "TimeNum","locType", "numer", "time", "visible",
                 "temp", "nu.x", "nu.y", "se.mu.x","se.nu.x","se.mu.y","se.nu.y", "speed", "x", "y")]
am93 <- am93[3307:11615,] ###
am99 <- data[134933:152473,] #

data <- rbind(am105,am107,am108,am110,am239,am253,am254,am255,am306,am307,am308,am91,am93,am99)

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
stepPar0 <- c(mu0,sigma0,zeromass0)
angleMean0 <- c(0.003592673,0.01941996,0.01499825) # angle mean
kappa0 <- c(0.355384957,1.53831499,3.59417278) # angle concentration 0.355384957,1.53831499,3.59417278
anglePar0 <- c(angleMean0)
## call to fitting function
stateNames <- c("Resting","Foraging", "Traveling")
m <- fitHMM(data=data, nbStates=nbStates,dist=list(step=stepDist,angle=angleDist), Par0=list(step=stepPar0,angle=anglePar0),formula=~temp, stateNames=stateNames)
m
CIreal(m)

#m <- readRDS("totalHMM3.rds")
plot(m, plotCI=TRUE, breaks = 40, plotStationary=TRUE,lwd=0.5)

plotPR(m, ncores=7)
AIC(m)

saveRDS(m, "totalHMM3.rds")
#####################
## custom plot ######
#####################

library(ggplot2)

m <- readRDS("totalHMM3.rds")

# calculate frequencies of states
viterbi_ratio <- viterbi(m)
stateFreq <- table(viterbi_ratio) / length(viterbi_ratio)

# default plot.momentuhmm colours
colours.states <- c("#E69F00", "#56B4E9", "#009E73")

# generate sequence for x axis of density functions
x <- seq(0, 6000, length=6000)

# get converged mean and sd for each state 
meanEn <- m$mle$step[1,1]  # 10.87739
sdEn <- m$mle$step[2,1]    # 8.775794

meanEx <- m$mle$step[1,2]    # 60.38591
sdEx <- m$mle$step[2,2]      # 43.38034

meanTr <- m$mle$step[1,3]   # 304.8606 
sdTr <- m$mle$step[2,3]     # 211.0313

# calculate shape and scale of the gamma distributions from mean and sd
sh <- function(mean, sd) { return(mean^2 / sd^2)}
sc <- function(mean, sd) { return(sd^2 / mean)}

# get density functions of the distributions
y_en <- dgamma(x, shape=sh(meanEn,sdEn),  scale=sc(meanEn,sdEn)) * stateFreq[[1]]
y_ex <- dgamma(x, shape=sh(meanEx,sdEx),  scale=sc(meanEx,sdEx)) * stateFreq[[2]]
y_tr <- dgamma(x, shape=sh(meanTr,sdTr),  scale=sc(meanTr,sdTr)) * stateFreq[[3]]
# sum densities to get total
y_tot <- y_en + y_ex + y_tr

# combine densities in a single dataframe for more convenient plotting
df.y_en <- data.frame(dens=y_en, state="Resting", x=x)
df.y_ex <- data.frame(dens=y_ex,  state="Foraging", x=x)
df.y_tr <- data.frame(dens=y_tr,  state="Traveling", x=x)
df.y_tot <- data.frame(dens=y_tot,  state="Total", x=x)

cmb <- rbind(df.y_en,df.y_ex,df.y_tr,df.y_tot)

# reorder factor levels so "total" appears bottom of the legend
cmb$state <- factor(cmb$state, levels=c("Resting","Foraging","Traveling","Total"))

# plot distributions
ggplot() +
  geom_line(data=cmb,aes(x=x,y=dens,colour=state,linetype=state), size=0.9) +
  scale_colour_manual(values=c(colours.states,"#000000")) +
  scale_linetype_manual(values=c("solid","solid","solid","dashed")) +
  scale_y_continuous(limits=c(0,0.0035)) +
  scale_x_continuous(limits=c(0,6000)) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  ylab("density") +
  xlab("step") +
  ggtitle('Fitted HMM - 3 states')

