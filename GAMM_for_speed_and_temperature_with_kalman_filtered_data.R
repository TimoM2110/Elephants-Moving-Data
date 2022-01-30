library(momentuHMM)
library(sp)
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

# run a GAMM using the mgcv package
mod.speed = bam(v ~ s(temp, k = 4) + 
                  season + 
                  s(id, bs = "re") + 
                  s(hour, bs = "re"), rho=0.25, 
                data = data)


summary(mod.speed)

mod.speed2 = gam(v ~ s(temp, k = 6) + 
                  season + 
                  s(id, bs = "re") + 
                  s(hour, bs = "re"), correlation = corAR1(form = ~ time | id),
                data = data)

summary(mod.speed2)

# prepare data for plotting
ele.speed.temp = 
  data %>%
  mutate(v.pred = predict(mod.speed2, 
                          newdata = ., scale = "response", 
                          allow.new.levels = T), 
         temp = plyr::round_any(temp,2)) %>%
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

# figure for speed and temperature
fig_speed_temp =
  ele.speed.temp %>% 
  filter(temp %in% 15:40) %>%
  ggplot()+
  
  geom_rangeframe(data = data_frame(x=c(15,40), y = c(0.15,1.35)),aes(x,y))+
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
  
#  theme_few()+ #panel.border = element_blank(),
  theme(
        legend.position = "none")+
  coord_cartesian(ylim=c(0.15,1.35))+
  scale_y_continuous(breaks = seq(.15,1.35,.15))+
  labs(x = "Collar temperature (°C)", y = "Speed (km/h)",
       col = "Season", fill = "Season", title="Speed and temperature")

fig_speed_temp


acf(data$v, lag.max = 10,
    plot = TRUE)
