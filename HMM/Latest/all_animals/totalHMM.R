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

plotStates(m)

#####################################
##  (Global) State Decoding #########
#####################################
# add index
data$locType <- seq.int(nrow(data))

### AM105
vitstates <- viterbi(m)
viterbi105 <- vitstates[1:1041]
track105 <- data[1:1041,]
colors <- c("#E69F00", "#56B4E9", "#009E73")
plot(track105$time,track105$step,type="h",xlab="time",
     ylab="step length in metres",main="AM105 decoded states",
     col=colors[viterbi105])
legend("topright", legend=c("Resting","Foraging", "Traveling"), 
       lwd=2, col = colors, cex = 0.55)

### AM107
viterbi107 <- vitstates[1042:3879]
track107 <- data[1042:3879,]
colors <- c("#E69F00", "#56B4E9", "#009E73")
plot(track107$time,track107$step,type="h",xlab="time",
     ylab="step length in metres",main="AM107 decoded states",
     col=colors[viterbi107])
legend("topright", legend=c("Resting","Foraging", "Traveling"), 
       lwd=2, col = colors, cex = 0.55)

### AM108
viterbi108 <- vitstates[3880:5049]
track108 <- data[3880:5049,]
colors <- c("#E69F00", "#56B4E9", "#009E73")
plot(track108$time,track108$step,type="h",xlab="time",
     ylab="step length in metres",main="AM108 decoded states",
     col=colors[viterbi108])
legend("topright", legend=c("Resting","Foraging", "Traveling"), 
       lwd=2, col = colors, cex = 0.55)

### AM110
viterbi110 <- vitstates[5050:7589]
track110 <- data[5050:7589,]
colors <- c("#E69F00", "#56B4E9", "#009E73")
plot(track110$time,track110$step,type="h",xlab="time",
     ylab="step length in metres",main="AM110 decoded states",
     col=colors[viterbi110])
legend("topright", legend=c("Resting","Foraging", "Traveling"), 
       lwd=2, col = colors, cex = 0.55)

### AM239
viterbi239 <- vitstates[7590:25129]
track239 <- data[7590:25129,]
colors <- c("#E69F00", "#56B4E9", "#009E73")
plot(track239$time,track239$step,type="h",xlab="time",
     ylab="step length in metres",main="AM239 decoded states",
     col=colors[viterbi239])
legend("topright", legend=c("Resting","Foraging", "Traveling"), 
       lwd=2, col = colors, cex = 0.55)

### AM253
viterbi253 <- vitstates[25130:35418]
track253 <- data[25130:35418,]
colors <- c("#E69F00", "#56B4E9", "#009E73")
plot(track253$time,track253$step,type="h",xlab="time",
     ylab="step length in metres",main="AM253 decoded states",
     col=colors[viterbi253])
legend("topright", legend=c("Resting","Foraging", "Traveling"), 
       lwd=2, col = colors, cex = 0.55)

### AM254
viterbi254 <- vitstates[35419:39957]
track254 <- data[35419:39957,]
colors <- c("#E69F00", "#56B4E9", "#009E73")
plot(track254$time,track254$step,type="h",xlab="time",
     ylab="step length in metres",main="AM254 decoded states",
     col=colors[viterbi254])
legend("topright", legend=c("Resting","Foraging", "Traveling"), 
       lwd=2, col = colors, cex = 0.55)

### AM255
viterbi255 <- vitstates[39958:45505]
track255 <- data[39958:45505,]
colors <- c("#E69F00", "#56B4E9", "#009E73")
plot(track255$time,track255$step,type="h",xlab="time",
     ylab="step length in metres",main="AM255 decoded states",
     col=colors[viterbi255])
legend("topright", legend=c("Resting","Foraging", "Traveling"), 
       lwd=2, col = colors, cex = 0.55)

### AM306
viterbi306 <- vitstates[45506:47281]
track306 <- data[45506:47281,]
colors <- c("#E69F00", "#56B4E9", "#009E73")
plot(track306$time,track306$step,type="h",xlab="time",
     ylab="step length in metres",main="AM306 decoded states",
     col=colors[viterbi306])
legend("topright", legend=c("Resting","Foraging", "Traveling"), 
       lwd=2, col = colors, cex = 0.55)

### AM307
viterbi307 <- vitstates[47282:53784]
track307 <- data[47282:53784,]
colors <- c("#E69F00", "#56B4E9", "#009E73")
plot(track307$time,track307$step,type="h",xlab="time",
     ylab="step length in metres",main="AM307 decoded states",
     col=colors[viterbi307])
legend("topright", legend=c("Resting","Foraging", "Traveling"), 
       lwd=2, col = colors, cex = 0.55)

### AM308
viterbi308 <- vitstates[53785:57536]
track308 <- data[53785:57536,]
colors <- c("#E69F00", "#56B4E9", "#009E73")
plot(track308$time,track308$step,type="h",xlab="time",
     ylab="step length in metres",main="AM308 decoded states",
     col=colors[viterbi308])
legend("topright", legend=c("Resting","Foraging", "Traveling"), 
       lwd=2, col = colors, cex = 0.55)

### AM91
viterbi91 <- vitstates[57537:75067]
track91 <- data[57537:75067,]
colors <- c("#E69F00", "#56B4E9", "#009E73")
plot(track91$time,track91$step,type="h",xlab="time",
     ylab="step length in metres",main="AM91 decoded states",
     col=colors[viterbi91])
legend("topright", legend=c("Resting","Foraging", "Traveling"), 
       lwd=2, col = colors, cex = 0.55)

### AM93
viterbi93 <- vitstates[75068:83376]
track93 <- data[75068:83376,]
colors <- c("#E69F00", "#56B4E9", "#009E73")
plot(track93$time,track93$step,type="h",xlab="time",
     ylab="step length in metres",main="AM93 decoded states",
     col=colors[viterbi93])
legend("topright", legend=c("Resting","Foraging", "Traveling"), 
       lwd=2, col = colors, cex = 0.55)

### AM99
viterbi99 <- vitstates[83377:100917]
track99 <- data[83377:100917,]
colors <- c("#E69F00", "#56B4E9", "#009E73")
plot(track99$time,track99$step,type="h",xlab="time",
     ylab="step length in metres",main="AM99 decoded states",
     col=colors[viterbi99])
legend("topright", legend=c("Resting","Foraging", "Traveling"), 
       lwd=2, col = colors, cex = 0.55)

################################
### Global vs Local Decoding ###
################################

sp <- stateProbs(m)

### AM105
sp105 <- sp[1:1041,]

par(mfrow=c(1,1))
plot(viterbi105,pch=16,col=colors[apply(sp105,1,which.max)],ylab="states according to viterbi",main="Decoded time series")
legend("topright", inset = c(0.1,0.1), title="states acc. to state prob.s", legend= c("Resting","Foraging","Traveling"), pch = 16, col=colors)

par(mfrow=c(2,1))
plot(track105$time,track105$step,type="h",xlab="time",ylab="step length in metres",main="Global decoding AM105",col=colors[viterbi105],lwd=1)
plot(track105$time,track105$step,type="h",xlab="time",ylab="step length in metres",main="Local decoding AM105",col=colors[apply(sp105,1,which.max)],lwd=1)

### AM107
sp107 <- sp[1042:3879,]

par(mfrow=c(1,1))
plot(viterbi107,pch=16,col=colors[apply(sp107,1,which.max)],ylab="states according to viterbi",main="Decoded time series")
legend("topright", inset = c(0.1,0.1), title="states acc. to state prob.s", legend= c("Resting","Foraging","Traveling"), pch = 16, col=colors)

par(mfrow=c(2,1))
plot(track107$time,track107$step,type="h",xlab="time",ylab="step length in metres",main="Global decoding AM107",col=colors[viterbi107],lwd=1)
plot(track107$time,track107$step,type="h",xlab="time",ylab="step length in metres",main="Local decoding AM107",col=colors[apply(sp107,1,which.max)],lwd=1)

### AM108
sp108 <- sp[3880:5049,]

par(mfrow=c(1,1))
plot(viterbi108,pch=16,col=colors[apply(sp108,1,which.max)],ylab="states according to viterbi",main="Decoded time series")
legend("topright", inset = c(0.1,0.1), title="states acc. to state prob.s", legend= c("Resting","Foraging","Traveling"), pch = 16, col=colors)

par(mfrow=c(2,1))
plot(track108$time,track108$step,type="h",xlab="time",ylab="step length in metres",main="Global decoding AM108",col=colors[viterbi108],lwd=1)
plot(track108$time,track108$step,type="h",xlab="time",ylab="step length in metres",main="Local decoding AM108",col=colors[apply(sp108,1,which.max)],lwd=1)

### AM110
sp110 <- sp[5050:7589,]

par(mfrow=c(1,1))
plot(viterbi110,pch=16,col=colors[apply(sp110,1,which.max)],ylab="states according to viterbi",main="Decoded time series")
legend("topright", inset = c(0.1,0.1), title="states acc. to state prob.s", legend= c("Resting","Foraging","Traveling"), pch = 16, col=colors)

par(mfrow=c(2,1))
plot(track110$time,track110$step,type="h",xlab="time",ylab="step length in metres",main="Global decoding AM110",col=colors[viterbi110],lwd=1)
plot(track110$time,track110$step,type="h",xlab="time",ylab="step length in metres",main="Local decoding AM110",col=colors[apply(sp110,1,which.max)],lwd=1)

### AM239
sp239 <- sp[7590:25129,]

par(mfrow=c(1,1))
plot(viterbi239,pch=16,col=colors[apply(sp239,1,which.max)],ylab="states according to viterbi",main="Decoded time series")
legend("topright", inset = c(0.1,0.1), title="states acc. to state prob.s", legend= c("Resting","Foraging","Traveling"), pch = 16, col=colors)

par(mfrow=c(2,1))
plot(track239$time,track239$step,type="h",xlab="time",ylab="step length in metres",main="Global decoding AM239",col=colors[viterbi239],lwd=1)
plot(track239$time,track239$step,type="h",xlab="time",ylab="step length in metres",main="Local decoding AM239",col=colors[apply(sp239,1,which.max)],lwd=1)

### AM253
sp253 <- sp[25130:35418,]

par(mfrow=c(1,1))
plot(viterbi253,pch=16,col=colors[apply(sp253,1,which.max)],ylab="states according to viterbi",main="Decoded time series")
legend("topright", inset = c(0.1,0.1), title="states acc. to state prob.s", legend= c("Resting","Foraging","Traveling"), pch = 16, col=colors)

par(mfrow=c(2,1))
plot(track253$time,track253$step,type="h",xlab="time",ylab="step length in metres",main="Global decoding AM253",col=colors[viterbi253],lwd=1)
plot(track253$time,track253$step,type="h",xlab="time",ylab="step length in metres",main="Local decoding AM253",col=colors[apply(sp253,1,which.max)],lwd=1)

### AM254
sp254 <- sp[35419:39957,]

par(mfrow=c(1,1))
plot(viterbi254,pch=16,col=colors[apply(sp254,1,which.max)],ylab="states according to viterbi",main="Decoded time series")
legend("topright", inset = c(0.1,0.1), title="states acc. to state prob.s", legend= c("Resting","Foraging","Traveling"), pch = 16, col=colors)

par(mfrow=c(2,1))
plot(track254$time,track254$step,type="h",xlab="time",ylab="step length in metres",main="Global decoding AM254",col=colors[viterbi254],lwd=1)
plot(track254$time,track254$step,type="h",xlab="time",ylab="step length in metres",main="Local decoding AM254",col=colors[apply(sp254,1,which.max)],lwd=1)

### AM255
sp255 <- sp[39958:45505,]

par(mfrow=c(1,1))
plot(viterbi255,pch=16,col=colors[apply(sp255,1,which.max)],ylab="states according to viterbi",main="Decoded time series")
legend("topright", inset = c(0.1,0.1), title="states acc. to state prob.s", legend= c("Resting","Foraging","Traveling"), pch = 16, col=colors)

par(mfrow=c(2,1))
plot(track255$time,track255$step,type="h",xlab="time",ylab="step length in metres",main="Global decoding AM255",col=colors[viterbi255],lwd=1)
plot(track255$time,track255$step,type="h",xlab="time",ylab="step length in metres",main="Local decoding AM255",col=colors[apply(sp255,1,which.max)],lwd=1)

### AM306
sp306 <- sp[45506:47281,]

par(mfrow=c(1,1))
plot(viterbi306,pch=16,col=colors[apply(sp306,1,which.max)],ylab="states according to viterbi",main="Decoded time series")
legend("topright", inset = c(0.1,0.1), title="states acc. to state prob.s", legend= c("Resting","Foraging","Traveling"), pch = 16, col=colors)

par(mfrow=c(2,1))
plot(track306$time,track306$step,type="h",xlab="time",ylab="step length in metres",main="Global decoding AM306",col=colors[viterbi306],lwd=1)
plot(track306$time,track306$step,type="h",xlab="time",ylab="step length in metres",main="Local decoding AM306",col=colors[apply(sp306,1,which.max)],lwd=1)

### AM307
sp307 <- sp[47282:53784,]

par(mfrow=c(1,1))
plot(viterbi307,pch=16,col=colors[apply(sp307,1,which.max)],ylab="states according to viterbi",main="Decoded time series")
legend("topright", inset = c(0.1,0.1), title="states acc. to state prob.s", legend= c("Resting","Foraging","Traveling"), pch = 16, col=colors)

par(mfrow=c(2,1))
plot(track307$time,track307$step,type="h",xlab="time",ylab="step length in metres",main="Global decoding AM307",col=colors[viterbi307],lwd=1)
plot(track307$time,track307$step,type="h",xlab="time",ylab="step length in metres",main="Local decoding AM307",col=colors[apply(sp307,1,which.max)],lwd=1)

### AM308
sp308 <- sp[53785:57536,]

par(mfrow=c(1,1))
plot(viterbi308,pch=16,col=colors[apply(sp308,1,which.max)],ylab="states according to viterbi",main="Decoded time series",lwd=1)
legend("topright", inset = c(0.1,0.1), title="states acc. to state prob.s", legend= c("Resting","Foraging","Traveling"), pch = 16, col=colors)

par(mfrow=c(2,1))
plot(track308$time,track308$step,type="h",xlab="time",ylab="step length in metres",main="Global decoding AM308",col=colors[viterbi308],lwd=1)
plot(track308$time,track308$step,type="h",xlab="time",ylab="step length in metres",main="Local decoding AM308",col=colors[apply(sp308,1,which.max)],lwd=1)

### AM91
sp91 <- sp[57537:75067,]

par(mfrow=c(1,1))
plot(viterbi91,pch=16,col=colors[apply(sp91,1,which.max)],ylab="states according to viterbi",main="Decoded time series",lwd=1)
legend("topright", inset = c(0.1,0.1), title="states acc. to state prob.s", legend= c("Resting","Foraging","Traveling"), pch = 16, col=colors)

par(mfrow=c(2,1))
plot(track91$time,track91$step,type="h",xlab="time",ylab="step length in metres",main="Global decoding AM91",col=colors[viterbi91],lwd=1)
plot(track91$time,track91$step,type="h",xlab="time",ylab="step length in metres",main="Local decoding AM91",col=colors[apply(sp91,1,which.max)],lwd=1)

### AM93
sp93 <- sp[75068:83376,]

par(mfrow=c(1,1))
plot(viterbi93,pch=16,col=colors[apply(sp93,1,which.max)],ylab="states according to viterbi",main="Decoded time series",lwd=1)
legend("topright", inset = c(0.1,0.1), title="states acc. to state prob.s", legend= c("Resting","Foraging","Traveling"), pch = 16, col=colors)

par(mfrow=c(2,1))
plot(track93$time,track93$step,type="h",xlab="time",ylab="step length in metres",main="Global decoding AM93",col=colors[viterbi93],lwd=1)
plot(track93$time,track93$step,type="h",xlab="time",ylab="step length in metres",main="Local decoding AM93",col=colors[apply(sp93,1,which.max)],lwd=1)

### AM99
sp99 <- sp[83377:100917,]

par(mfrow=c(1,1))
plot(viterbi99,pch=16,col=colors[apply(sp99,1,which.max)],ylab="states according to viterbi",main="Decoded time series",lwd=0.1)
legend("topright", inset = c(0.1,0.1), title="states acc. to state prob.s", legend= c("Resting","Foraging","Traveling"), pch = 16, col=colors)

par(mfrow=c(2,1))
plot(track99$time,track99$step,type="h",xlab="time",ylab="step length in metres",main="Global decoding AM99",col=colors[viterbi99],lwd=1)
plot(track99$time,track99$step,type="h",xlab="time",ylab="step length in metres",main="Local decoding AM99",col=colors[apply(sp99,1,which.max)],lwd=1)
