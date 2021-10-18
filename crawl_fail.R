library(crawl)
library(tidyverse)
library(lubridate)
library(sp)
library(sf)

setwd("C:/Users/timom/anaconda3/envs/madeleine_project")
data <- read.table("ThermochronTracking Elephants Kruger 2007.csv", header = TRUE, sep= ",")

data$timestamp <- as.POSIXct(strptime(data$timestamp, format='%Y-%m-%d %H:%M:%S'))

indx <- apply(data, 2, function(x) any(is.na(x) | is.infinite(x)))
colnames[indx]

data <- na.omit(data)
datas <- na.omit(datas)

names(data)[names(data) == "timestamp"] <- "date_time"


data['x'] = datas['longitude'] #erst untitled 8 laufen lassen
data['y'] = datas['latitude'] # dort wurde longlat in utm umgewandelt

#names(data)[names(data) == "location.long"] <- "longitude"
#names(data)[names(data) == "location.lat"] <- "latitude"


sf_locs <- sf::st_as_sf(data, coords = c("x","y")) %>%  #tbl locs
  sf::st_set_crs(3857)

library(sp)
library(rgdal)

names(sf_locs)[names(sf_locs) == "individual.local.identifier"] <- "ID"

sf_lines <- sf_locs %>% 
  dplyr::arrange(ID, date_time) %>% 
  sf::st_geometry() %>% 
  sf::st_cast("MULTIPOINT",ids = as.integer(as.factor(sf_locs$ID))) %>% 
  sf::st_cast("MULTILINESTRING") %>% 
  sf::st_sf(deployid = as.factor(unique(sf_locs$ID)))

library(ggplot2)
library(ggspatial)

ggplot() + 
  annotation_map_tile(zoomin = 1,progress = "none") +
  layer_spatial(sf_lines, size = 0.75,aes(color = deployid)) +
  scale_x_continuous(expand = expand_scale(mult = c(.6, .6))) +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "none") +
  ggtitle("Observed Argos Location Paths", 
          subtitle = "harbor seals (n=2), Aleutian Islands, Alaska, USA")

sf_locs <- sf::st_transform(sf_locs, 3857)

sf::st_geometry(sf_locs)

sfc_as_cols <- function(x, names = c("x","y")) {
  stopifnot(inherits(x,"sf") && inherits(sf::st_geometry(x),"sfc_POINT"))
  ret <- sf::st_coordinates(x)
  ret <- tibble::as_tibble(ret)
  stopifnot(length(names) == ncol(ret))
  x <- x[ , !names(x) %in% names]
  ret <- setNames(ret,names)
  dplyr::bind_cols(x,ret)
}

sf_locs <- sf_locs %>% 
  mutate(date_hour = lubridate::floor_date(date_time,'hour')) %>% 
  arrange(ID, date_time) %>% 
  sfc_as_cols()

#initial = list(a=c(coordinates(sf_locs)[1,1],0,
#                   coordinates(sf_locs)[1,2],0),
#               P=diag(c(10000^2,54000^2,10000^2,5400^2)))

initial = list(a = c(sf_locs$x[1], 0,
           sf_locs$y[1], 0),
     P = diag(c(10 ^ 2, 10 ^ 2,
                10 ^ 2, 10 ^ 2)))


############################################################

fit1 <- crwMLE(mov.model=~1,
               data=sf_locs,
               Time.name="date_time",
               initial.state=initial,
               control=list(maxit=30,trace=0),
               initialSANN=list(maxit=200,trace=0))

