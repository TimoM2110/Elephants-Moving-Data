library(momentuHMM)
library(raster)
setwd("C:/Users/timom/anaconda3/envs/madeleine_project")
databz <- read.table("imputed_dataframe_with_rivers_waterholes_and_distance.csv", header = TRUE, sep= ",")
databz <- databz[ -c(2,3,4) ]
databz$time <- as.POSIXlt(databz$time, tz = "GMT")

setwd("C:/Users/timom")
dataz <- read.table("imputed_data_lonlat_without_waterholes.csv", header = TRUE, sep= ",")
dataz <- dataz[ -c(2,3,4) ]

databz <- databz[which(databz$ID %in% c("AM110")),]
dataz <- dataz[which(dataz$ID %in% c("AM110")),]

#dataz <- subset(
#  dataz,
#  ID == unique(ID)[4] | ID == unique(ID)[4] # keep AM99 only
#) 

databz['x'] = dataz['x']
databz['y'] = dataz['y'] 

databz <- na.omit(databz)
dataz <- na.omit(dataz)
databz <- databz[1:2530,]
dataz <- dataz[1:2530,]


library(sf)
dist2water <- databz[,(c(15,16,14))]
databz <- databz[, -c(14)]
datas <- st_as_sf(databz, coords=c("x", "y"), crs="+proj=longlat +datum=WGS84 +ellps=WGS84") #36J
cords <- st_as_sf(dist2water, coords=c("x", "y"), crs="+proj=longlat +datum=WGS84 +ellps=WGS84") #36J
dist2sabie <- raster(cords, ncol=2530/10) #anpassen
ncell(dist2sabie)
values <- cords$mindw
length(values)
tail(values)
values[17539] <- 100
values[17540] <- 30
values[10290] <- 1000
dist2sabie <- setValues(dist2sabie, values)
plot(dist2sabie)
try <- as.data.frame(dist2sabie, xy=TRUE)
dist2sabie_scaled <- dist2sabie / mean(values(terrain(dist2sabie,
                                                      opt = "slope")),
                                       na.rm = T)

# calculate gradient
D_scaled <- ctmcmove::rast.grad(dist2sabie_scaled)
## W (recharge function covariates)
# near_sabie = indicator for <500m from water
intercept <- raster(dist2sabie)
values(intercept) <- 1
W <- stack(list("intercept" = intercept,
                "near_sabie" = dist2sabie < 0.4e3))
W_names <- names(W)
## orthogonalize W based on locations ----
library(sp)
sp_df <- SpatialPointsDataFrame(databz[,c("x", "y")], databz)
W_ortho <- W
W_path <- extract(x = W, y = matrix(sp_df@coords, ncol = 2))
obstimes <- as.numeric(sp_df$time) / 3600 # numeric hours
W_tilde <- apply(W_path * c(0, diff(obstimes)), 2, cumsum)
W_tilde_svd <- svd(W_tilde)
W_tilde_proj_mat <- W_tilde_svd$v %*% diag(W_tilde_svd$d^(-1))
W_mat <- as.matrix(W)
W_mat_proj <- W_mat %*% W_tilde_proj_mat
for(layer in 1:ncol(W_mat)){
  values(W_ortho[[layer]]) <- W_mat_proj[, layer]
  names(W_ortho[[layer]]) <- paste0("svd", layer)
}

buffaloData <- data.frame(ID = sp_df$ID,
                          time = obstimes,
                          temp = sp_df$temp,
                          x = sp_df@coords[, 1],
                          y = sp_df@coords[, 2]
                          )


#crwOut <- crawlWrap(buffaloData,
#                    theta = c(2,-2),
#                    fixPar = c(NA,NA),
#                    timeStep = "1 hour", # predict at 15 min time steps
#                    attempts = 10)

crwOut <- crawlWrap(buffaloData, timeStep = 1, theta = c(3, -.1),fixPar = c(NA,NA),attempts = 100)
                    
spatialCovs <- list(W_intercept = W_ortho$svd1,
                    W_near_sabie = W_ortho$svd2,
                    dist2sabie = dist2sabie,
                    D.x = D_scaled$rast.grad.x,
                    D.y = D_scaled$rast.grad.y)
# best predicted track data
hmmData <- prepData(crwOut,
                    spatialCovs = spatialCovs,
                    altCoordNames = "mu")
head(hmmData[,c("ID","time", "mu.x", "mu.y",
                "W_intercept","W_near_sabie","dist2sabie",
                "D.x","D.y", "temp")])

nbStates <- 2
stateNames <- c("charged", "discharged")
dist <- list(mu = "rw_mvnorm2") # bivariate normal random walk
# pseudo-design matrix for mu
DM <- list(mu=matrix(c("mu.x_tm1", 0, 0,0,0,0,
                       "mu.x_tm1", 0,"D.x",0,0,0,
                       0,"mu.y_tm1", 0,0,0,0,
                       0,"mu.y_tm1","D.y",0,0,0,
                       0, 0, 0,1,0,0,
                       0, 0, 0,0,1,0,
                       0, 0, 0,0,0,1,
                       0, 0, 0,0,0,1,
                       0, 0, 0,1,0,0,
                       0, 0, 0,0,1,0),
                     5*nbStates,
                     6,byrow=TRUE,
                     dimnames=list(c(paste0("mean.",
                                            rep(c("x_","y_"),
                                                each=nbStates),
                                            1:nbStates),
                                     paste0("sigma.",
                                            rep(c("x_","xy_","y_"),
                                                each=nbStates),
                                            1:nbStates)),
                                   c("x:x_tm1",
                                     "y:y_tm1",
                                     "xy:D",
                                     "sigma_1:(Intercept)",
                                     "sigma_2:(Intercept)",
                                     "sigma_12:(Intercept)"))))
# starting values
Par0=list(mu=c(1, 1, 0, log(85872.66), log(37753.53), 0)) #3. und 6. war 0
g0 <- 0 # recharge function at time 0
theta <- c(0,20,40) # recharge function parameters
## specify recharge formula
# note that theta formula requires an 'intercept' term
formula <- ~ recharge(g0 = ~1,
                      theta = ~W_intercept+W_near_sabie)

## remove Markov property
betaRef <- c(1,1) # make state 1 the reference state
betaCons <- matrix(c(1,2),2,2) # 1 -> 1 = 2 -> 1 and 1 -> 2 = 2 -> 2
## set fixed parameters
fixPar <- list(mu = c(Par0$mu[1:2],NA,NA,NA,Par0$mu[6]),
               beta = matrix(c(0,-1,0,-1),2,2),
               delta = c(0.5,0.5),
               theta = c(0,NA,NA)) # fix extra 'intercept' term to zero
# check recharge model specification
checkPar0(hmmData, nbStates = nbStates, dist = dist,
          formula = formula, Par0 = Par0,
          beta0 = list(beta = fixPar$beta,
                       g0 = g0,
                       theta = theta),
          delta0 = fixPar$delta, fixPar = fixPar,
          DM = DM, betaRef = betaRef, betaCons = betaCons,
          stateNames = stateNames)

buffaloFit <- fitHMM(hmmData, nbStates = nbStates, dist = dist,
                     formula = formula, Par0 = Par0,
                     beta0 = list(g0=g0,
                                  theta=theta),
                     fixPar = fixPar,
                     DM = DM, betaRef = betaRef, betaCons = betaCons,
                     stateNames = stateNames,
                     mvnCoords = "mu",
                     optMethod = "Nelder-Mead",
                     control = list(maxit=1000))

saveRDS(buffaloFit, "recharge_model_AM110.rds")
#buffaloFit <- readRDS("recharge_model_AM253.rds")

plot(buffaloFit,plotCI=T, ask=F)
AIC(buffaloFit)

#plotSpatialCov(buffaloFit,dist2sabie)

bestPar <- getPar(buffaloFit)
buffaloFits <- MIfitHMM(crwOut, nSims=28,
                        spatialCovs = spatialCovs,
                        mvnCoords="mu", altCoordNames = "mu",
                        nbStates=nbStates, dist=dist, formula=formula,
                        Par0=bestPar$Par, beta0=bestPar$beta,
                        fixPar=fixPar, DM=DM,
                        betaRef=betaRef, betaCons=betaCons,
                        stateNames = stateNames,
                        retryFits = 3, retrySD=list(mu=c(0,0,3,0,0,0),
                                                    g0=1,
                                                    theta=c(0,1,1)),
                        optMethod = "Nelder-Mead",
                        control = list(maxit=100000))

saveRDS(buffaloFits, "recharge_model_AM306_2.rds")
#buffaloFits <- readRDS("recharge_model_total.rds")

plot(buffaloFits,plotCI=F, ask=F)
plotSpatialCov(buffaloFits,dist2sabie)

# plot estimates and CIs for Pr(discharged) at each time step

buffaloFits <- readRDS("recharge_model_AM105_to_AM110.rds")

trProbs <- getTrProbs(buffaloFit, getCI=T)

dataz$time <- as.Date(dataz$time, "%Y-%m-%d")
vitstates <- viterbi(buffaloFit)

plot(trProbs$est[1,2,],xaxt="n",type="l", ylim=c(0,1),
     ylab="Pr(discharged)", xlab="Time",
     col=c("#E69F00", "#56B4E9")[buffaloFits$miSum$Par$states],main= paste("ID", as.character(dataz$ID[1]), sep=" "))
arrows(1:dim(trProbs$est)[3],
       trProbs$lower[1,2,],
       1:dim(trProbs$est)[3],
       trProbs$upper[1,2,],
       length=0.025, angle=90, code=3,
       col=c("#E69F00", "#56B4E9")[buffaloFits$miSum$Par$states],
       lwd=1.3)
abline(h=0.5,lty=2)
axis(1, at = seq(1, length(trProbs$est[1,2,]), length.out=5), labels = seq(min(dataz$time), max(dataz$time),length.out=5),cex.axis = 1)

plot(trProbs$est[1,2,],xaxt="n",type="l", ylim=c(0,1),
     ylab="Pr(discharged)", xlab="Time",
     col=c("#E69F00", "#56B4E9")[vitstates],main= paste("ID", as.character(dataz$ID[1]), sep=" "))
arrows(1:dim(trProbs$est)[3],
       trProbs$lower[1,2,],
       1:dim(trProbs$est)[3],
       trProbs$upper[1,2,],
       length=0.025, angle=90, code=3,
       col=c("#E69F00", "#56B4E9")[vitstates],
       lwd=1.3)
abline(h=0.5,lty=2)
axis(1, at = seq(1, length(trProbs$est[1,2,]), length.out=5), labels = seq(min(dataz$time), max(dataz$time),length.out=5),cex.axis = 1)


