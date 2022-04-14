library(momentuHMM)
library(sp)
#library(devtools)
#install_github("gavinsimpson/gratia")
library(gratia)
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
#

track <- data[which(data$ID %in% c("AM91")),]

###### Einschub

setwd("C:/Users/timom/anaconda3/envs/madeleine_project")
data <- read.table("crawl_data30_public_best_imp.csv", header = TRUE, sep= ",")

data$time <- as.POSIXct(strptime(data$time, format='%Y-%m-%d %H:%M:%S'))
data <- data[ -c(2,3) ]
data <- prepData(data,type="UTM",coordNames=c("x","y"))

AM105 <- data[29000:31079,] ##
AM107 <- data[35086:40760,] ##
AM108 <- data[81966:84304,] ##
AM110 <- data[98549:104587,] ##
AM239 <- data[132590:167668,] ##
AM253 <- data[170374:190950,] ##
AM254 <- data[202956:212032,] ##
AM255 <- data[233870:244964,] ##
AM306 <- data[253470:257020,] ##
AM307 <- data[272018:285022,] ##
AM308 <- data[307628:315130,] ##
AM93  <- data[327745:344361,] ##
AM99  <- data[356214:391295,] ##
AM91  <- data[391297:426357,] ##

data <- rbind(AM105,AM107,AM108,AM110,AM239,AM253,AM254,AM255,AM306,AM307,AM308,AM91,AM93,AM99)
#######
# Einschub Ende
#########


library(lubridate)

data$month <- month(as.POSIXlt(data$time, format="%Y/%m/%d %H:%M:%S"))
data$hour <- hour(as.POSIXlt(data$time, format="%Y/%m/%d %H:%M:%S"))
data$day <- (as.POSIXlt(data$time, format="%Y/%m/%d %H:%M:%S"))

#plot(data,compact=T)

# load libraries
library(dplyr)
library(purrr)
library(readr)

# stats
library(mgcv)

# spatial
library(sf)
library(raster)

# plotting
library(ggplot2)
library(ggthemes)
library(viridis)

# custom funcs
ci = function(x) 1.96*sd(x, na.rm = T)/sqrt(length(x))

names(data)[names(data) == "speed"] <- "v"
names(data)[names(data) == "ID"] <- "id"

data$season <- data$month
data$season <- ifelse(data$month > 4 & data$month < 11,data$season <- "dry", data$season <- "wet")

# make id a factor
data$id = as.factor(data$id)
data$season = as.factor(data$season)


data$first <- is.character(data$numer)

counts <- data[, .(rowCount = .N), by = id]

data$first[1] <- T #105
data$first[1042] <- T #107
data$first[(1041+2838+1)] <- T #108
data$first[(1041+2838+1170+1)] <- T #110
data$first[(1041+2838+1170+2540+1)] <- T #239
data$first[(7590+17540)] <- T #253
data$first[(25130+10289)] <- T #254 check
data$first[(35419+4539)] <- T #255 check
data$first[(39958+5548)] <- T #306 check
data$first[(45506+1776)] <- T #307 check
data$first[(47282+6503)] <- T #308
data$first[(53785+3752)] <- T #91 check
data$first[(57537+17531)] <- T #93 check
data$first[(75068+8309)] <- T #99 check

data$coshour <- cos(2*pi*data$hour/24)
data$sinhour <- sin(2*pi*data$hour/24)

# run a GAMM using the mgcv package
mod.speed = bam(v ~ s(temp,  k = 40, bs="cp")+
                  season+
                  s(id, bs = "re")+
                  s(hour,k=50, bs= "cp"),
                  rho=0.25,# AR.start = data$first, 
                  data = data, family = "Gamma"(link = "log"),discrete=TRUE, gamma=1)

AIC(mod.speed)
logLik.gam(mod.speed)
AIC(mod.speed)

REML <- rho <- 0.1+0:69/100
for (i in 1:length(rho)) {
  mod.speed = bam(v ~ s(temp,  k = 10, bs="ps")+
                    season+
                    s(id, bs = "re")+
                    s(hour, bs= "cc"),
                  rho=rho[i], AR.start = data$first, 
                  data = data)
  REML[i] <- mod.speed$gcv.ubre
}

min(REML)
#vorher ohne cos1 und cos 2
REML[29]
rho[29]
summary(mod.speed)

draw(mod.speed, rug = F, scales = "fixed")

appraise(mod.speed)

gam.check(mod.speed, breaks = "scott", type = "deviance")

par(mfrow=c(1,2))
hist(residuals.gam(mod.speed), breaks = "scott", main = "Gaussian", ylim = c(0,10000))
hist(residuals.gam(mod.speed2), breaks = "scott", main = "Gamma", ylim = c(0,10000))
par(mfrow=c(1,1))

gaussian <- qq_plot(mod.speed, method = "normal", type = "deviance")
gamma <- qq_plot(mod.speed, method = "normal", type = "deviance")

qq.gam(mod.speed,type=c("deviance"), pch=19)

require(gridExtra)


library(mgcViz)
gaussian.speed <- getViz(mod.speed)
gamma.speed <- getViz(mod.speed)

gamma <- qq(gam.speed, method = "tnormal", type = "deviance", pch=19, discrete=T, CI = "normal", 
   a.qqpoi = list("shape" = 20, "size" = 4))
gaussian <- qq(gam.speed, method = "tnormal", type = "deviance", pch=19, discrete=T, CI = "normal", 
            a.qqpoi = list("shape" = 20, "size" = 4))

print(gaussian,gamma, ncol=2)

gridPrint(qq(gaussian.speed, level=0.95, method = "tnormal", type = "deviance", pch=19, discrete=T, CI = "normal", a.ablin = list(colour = "black"),
                  a.qqpoi = list("shape" = 20, "size" = 3, col = "red"))+ylab("Sample Quantiles"),
          qq(gamma.speed, level=0.95, method = "tnormal", type = "deviance", pch=19, discrete=T, CI = "normal", a.ablin = list(colour = "black"),
                                                                       a.qqpoi = list("shape" = 20, "size" = 3, col = "red"))+ylab("Sample Quantiles"), ncol = 2)


par(mfrow=c(1,2))
hist(residuals(gaussian.speed, type = "deviance"),xlab = "", breaks = 80, main = "Gaussian", ylim = c(0,20000))
hist(residuals(gamma.speed, type = "deviance"),xlab = "", breaks = 80, main = "Gamma", ylim = c(0,20000))
par(mfrow=c(1,1))

BIC(mod.speed)
AIC(mod.speed)
res_gam = residuals.gam(mod.speed)
res_gaus = residuals.gam(mod.speed)
par(mfrow=c(1,2))
acf(res_gaus, ylim=c(0,1), main ="Gaussian") #, main = "bam: AR(1)"
acf(res_gam, ylim=c(0,1),main ="Gamma")




track <- data[which(data$id %in% c("AM105", "AM108")),]

#data <- na.omit(data)

#dif <- diff(data$numer)
#min(dif)

data$index <- seq(from = 1, by = 1, to = 100917)


# prepare data for plotting
ele.speed.temp = 
  data %>%
  mutate(v.pred = predict(mod.speed, 
                          newdata = ., type = "response", 
                          allow.new.levels = T), 
         temp = plyr::round_any(temp,1)) %>%
  ungroup() %>% 
  
  group_by(season, temp) %>%
  
  summarise(v.mean = mean(v), 
            v.sd = sd(v), 
            n.v = length(v), 
            pred.mean = mean(v.pred, na.rm = T), 
            pred.sd=sd(v.pred, na.rm = T), 
            pred.n = length(v.pred)) %>%
  
  mutate(v.ci = qnorm(0.975)*v.sd/sqrt(n.v), 
         ci.pred = qnorm(0.975)*pred.sd/sqrt(pred.n))

ele.speed.temp

#ele.speed.temp$v.mean <- ele.speed.temp$v.mean/1000

mean(ele.speed.temp$v.mean[ele.speed.temp$season == "wet"])
mean(ele.speed.temp$v.mean[ele.speed.temp$season == "dry"])

ele.speed.temp <- na.omit(ele.speed.temp)

mean(ele.speed.temp$v.sd[ele.speed.temp$season == "wet"])
mean(ele.speed.temp$v.sd[ele.speed.temp$season == "dry"])

mean(ele.speed.temp$v.ci[ele.speed.temp$season == "wet"])
mean(ele.speed.temp$v.ci[ele.speed.temp$season == "dry"])

ele.speed.temp2 = 
  data %>%
  mutate(v.pred = predict(mod.speed, 
                          newdata = ., scale = "response", 
                          allow.new.levels = T), 
         temp = plyr::round_any(temp,1)) %>%
  ungroup() %>% 
  
  group_by(season, temp)

ele.speed.temp2$v.pred <- ele.speed.temp2$v.pred/1000

ele.speed.temp2 <- na.omit(ele.speed.temp2)

sd(ele.speed.temp2$v.pred)

sd(ele.speed.temp2$v.pred[ele.speed.temp2$season == "wet"])
sd(ele.speed.temp2$v.pred[ele.speed.temp2$season == "dry"])

# figure for speed and temperature
fig_speed_temp =
  ele.speed.temp %>% 
  filter(temp >= 15 & temp <= 42) %>%
  ggplot()+
  
  geom_rangeframe(data = data_frame(x=c(15,40), y = c(0.15,1.75)),aes(x,y))+
  geom_smooth(aes(x = temp, y = pred.mean*2/1e3, 
                  col = season, fill = season, lty = season), 
              alpha = 0.2, lwd = 0.5)+
  geom_pointrange(aes(x = temp, y = v.mean*2/1e3, 
                      ymin = (v.mean-v.ci)*2/1e3, ymax = (v.mean+v.ci)*2/1e3,
                      col = season, shape = season), 
                  fill = "white", size = 0.4, stroke =0.7, lty = 1, 
                  position = position_dodge(width = 0.3))+
  
  scale_fill_brewer(palette = "Set1")+
#  scale_color_brewer(palette = "Set1")+
  scale_shape_manual(values=c(21,24))+
  scale_linetype_manual(values=c("dashed","solid"))+
  
#  theme_few()+ #panel.border = element_blank(), , title="Speed and temperature / k = 10"
  theme(
        legend.position = "none")+
  coord_cartesian(ylim=c(0.15,1.75))+
  scale_y_continuous(breaks = seq(.15,1.75,.15))+
  scale_x_continuous(breaks = seq(15,40,5))+
  theme_bw()+
  labs(x = "Collar temperature (°C)", y = "Speed (km/h)",
       col = "Season", fill = "Season")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))

fig_speed_temp
require(gridExtra)

grid.arrange(arrangeGrob(fig_speed_temp2, top = 'Gaussian'), arrangeGrob(fig_speed_temp, top = 'Gamma'), ncol=2)


acf(data$v, lag.max = 10,
    plot = TRUE)


# prepare data for plotting
ele.speed.hour = 
  data %>%
  mutate(v.pred = predict(mod.speed, 
                          newdata = ., type = "response", 
                          allow.new.levels = T), 
         hour = plyr::round_any(hour,0.5)) %>%
  ungroup() %>% 
  
  group_by(season, hour) %>%
  
  summarise(v.mean = mean(v), 
            v.sd = sd(v), 
            n.v = length(v), 
            pred.mean = mean(v.pred, na.rm = T), 
            pred.sd=sd(v.pred, na.rm = T), 
            pred.n = length(v.pred)) %>%
  
  mutate(v.ci = qnorm(0.975)*v.sd/sqrt(n.v), 
         ci.pred = qnorm(0.975)*pred.sd/sqrt(pred.n))

fig_speed_hour =
  ele.speed.hour %>% 
#  filter(temp >= 15 & temp <= 42) %>%
  ggplot()+
  
  geom_rangeframe(data = data_frame(x=c(0,24), y = c(0.15,0.75)),aes(x,y))+
  geom_smooth(aes(x = hour, y = pred.mean*1/1e3, 
                  col = season, fill = season, lty = season), 
              alpha = 0.2, lwd = 0.5)+
  geom_pointrange(aes(x = hour, y = v.mean*1/1e3, 
                      ymin = (v.mean-v.ci)*1/1e3, ymax = (v.mean+v.ci)*1/1e3,
                      col = season, shape = season), 
                  fill = "white", size = 0.4, stroke =0.7, lty = 1, 
                  position = position_dodge(width = 0.3))+
  
  scale_fill_brewer(palette = "Set1")+
  #  scale_color_brewer(palette = "Set1")+
  scale_shape_manual(values=c(21,24))+
  scale_linetype_manual(values=c("dashed","solid"))+
  
  #  theme_few()+ #panel.border = element_blank(), , title="Speed and temperature / k = 10"
  theme(
    legend.position = "none")+
  coord_cartesian(ylim=c(0.15,0.75))+
  scale_y_continuous(breaks = seq(.15,0.75,.15))+
  scale_x_continuous(breaks = seq(0,24,1))+
  theme_bw()+
  labs(x = "Time of day (hours)", y = "Speed (km/h)",
       col = "Season", fill = "Season")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))

fig_speed_hour
