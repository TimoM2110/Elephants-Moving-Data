library(devtools)
install_github("TheoMichelot/smoothSDE")
library(smoothSDE)
setwd("C:/Users/timom/anaconda3/envs/madeleine_project")
data <- read.table("ThermochronTracking Elephants Kruger 2007.csv", header = TRUE, sep= ",")

names(data)[names(data) == "location.long"] <- "x"
names(data)[names(data) == "location.lat"] <- "y"
names(data)[names(data) == "timestamp"] <- "time"
names(data)[names(data) == "individual.local.identifier"] <- "ID"
names(data)[names(data) == "external.temperature"] <- "temp"

data$time <- as.POSIXct(strptime(data$time, format='%Y-%m-%d %H:%M:%S'))

####
#converting to utm
####
library(lubridate)
library(sp)
library(sf)
library(tidyverse)

indx <- apply(data, 2, function(x) any(is.na(x) | is.infinite(x)))
colnames[indx]

data <- na.omit(data)


data <- data[ -c(1, 2,7:9,11) ]

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
data['x'] = datas['x']
data['y'] = datas['y'] 

# Only keep five months to eliminate seasonal effects
track <- subset(data, ID == unique(ID)[12]) # keep AM91 only
names(track)[names(track) == "time"] <- "date"
track$date <- as.POSIXlt(track$date, tz = "GMT")
track$time <- as.numeric(track$date - min(track$date))/3600
keep_rows <- which(track$date > as.POSIXct("2008-05-01 00:00:00") &
                     track$date < as.POSIXct("2008-09-30 23:59:59"))
track_1 <- track[keep_rows,]
# Convert to km
track_1$x <- track_1$x/1000
track_1$y <- track_1$y/1000
# Data set including ID, time, responses (x, y), and covariate (temp)
head(track)

# Using the Ornstein-Uhlenbeck process (also called Continuous Time Correlated Random Walk (CTCRW))
# Model formulas for CTCRW parameters
# 2 parameters, beta for mean reversion and sigma for variance
formulas <- list(beta = ~ s(temp, k = 10, bs = "ts"),
                 sigma = ~ s(temp, k = 10, bs = "ts"))
# Type of SDE model: continuous-time correlated random walk
type <- "CTCRW"

# Create SDE model object
my_sde <- SDE$new(formulas = formulas, data = track_1, type = type,
                  response = c("x", "y"))

# Fit model
my_sde$fit()
# Plot the model
my_sde$plot_par("temp", n_post = 500)

##################
## Plot results ##
##################
# Build design matrices
mats <- my_sde$make_mat_grid(var = "temp")

# Parameter estimates
par <- my_sde$par(t = "all", X_fe = mats$X_fe, X_re = mats$X_re)

# Posterior samples
n_post <- 2000
post <- my_sde$post_par(X_fe = mats$X_fe, X_re = mats$X_re, n_post = n_post)

# Get nu parameter (mean speed)
nu <- sqrt(pi) * par[,"sigma"] / (2 * sqrt(par[,"beta"]))
post_nu <- sqrt(pi) * post[,"sigma",] / (2 * sqrt(post[,"beta",]))

# Confidence bands
quants <- apply(post_nu, 1, quantile, probs = c(0.025, 0.975))
new_data <- data.frame(temp = seq(from = min(track_1$temp), to = max(track_1$temp), length.out = 1000))
df <- data.frame(x = new_data$temp,
                 mle = nu,
                 low = quants[1,],
                 upp = quants[2,])

# Plot
p_est <- ggplot(df, aes(x, mle)) +
  xlab("temperature (°C)") + ylab(bquote(nu[t])) + geom_line(size = 1) +
  geom_ribbon(aes(ymin = low, ymax = upp), alpha = 0.3) +
  theme_light() +
  theme(axis.title = element_text(size = 12),
        strip.text = element_text(color = 1, size = 12))
p_est

################
###2nd period###
################
# 7 months of warmer period
keep_rows_2 <- which(track$date > as.POSIXct("2008-10-01 00:00:00") &
                       track$date < as.POSIXct("2009-04-30 23:59:59"))
track_2 <- track[keep_rows_2,]
# Convert to km
track_2$x <- track_2$x/1000
track_2$y <- track_2$y/1000
# Data set including ID, time, responses (x, y), and covariate (temp)
head(track_2)

# Using the Ornstein-Uhlenbeck process (also called Continuous Time Correlated Random Walk (CTCRW))
# Model formulas for CTCRW parameters
# 2 parameters, beta for mean reversion and sigma for variance
formulas <- list(beta = ~ s(temp, k = 10, bs = "ts"),
                 sigma = ~ s(temp, k = 10, bs = "ts"))
# Type of SDE model: continuous-time correlated random walk
type <- "CTCRW"

# Create SDE model object
my_sde_2 <- SDE$new(formulas = formulas, data = track_2, type = type,
                    response = c("x", "y"))

# Fit model
my_sde_2$fit()
# Plot the model
my_sde_2$plot_par("temp", n_post = 500)

##################
## Plot results ##
##################
# Build design matrices
mats <- my_sde_2$make_mat_grid(var = "temp")

# Parameter estimates
par <- my_sde_2$par(t = "all", X_fe = mats$X_fe, X_re = mats$X_re)

# Posterior samples
n_post <- 2000
post <- my_sde_2$post_par(X_fe = mats$X_fe, X_re = mats$X_re, n_post = n_post)

# Get nu parameter (mean speed)
nu <- sqrt(pi) * par[,"sigma"] / (2 * sqrt(par[,"beta"]))
post_nu <- sqrt(pi) * post[,"sigma",] / (2 * sqrt(post[,"beta",]))

# Confidence bands
quants <- apply(post_nu, 1, quantile, probs = c(0.025, 0.975))
new_data <- data.frame(temp = seq(from = min(track_2$temp), to = max(track_2$temp), length.out = 1000))
df <- data.frame(x = new_data$temp,
                 mle = nu,
                 low = quants[1,],
                 upp = quants[2,])

# Plot
p_est <- ggplot(df, aes(x, mle)) +
  xlab("temperature (°C)") + ylab(bquote(nu[t])) + geom_line(size = 1) +
  geom_ribbon(aes(ymin = low, ymax = upp), alpha = 0.3) +
  theme_light() +
  theme(axis.title = element_text(size = 12),
        strip.text = element_text(color = 1, size = 12))
p_est
