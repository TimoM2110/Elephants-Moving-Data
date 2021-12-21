# load libraries
library(dplyr)
library(purrr)
library(lubridate)
library(readr)
library(glue)
library(stringr)

# spatial
library(sf)

# plotting
library(ggplot2)
library(ggthemes)
library(viridis)
setwd("C:/Users/timom/anaconda3/envs/madeleine_project")

# read in the data
data <- read.table(
  "ThermochronTracking Elephants Kruger 2007.csv",
  header = TRUE,
  sep= ","
)
ele$xutm <- 1
ele$yutm <- 2

data_sf = st_as_sf(ele, coords = c("xutm", "yutm")) %>%
  `st_crs<-`(32736)

rivers = st_read("data/rivers/kruger_rivers_cropped.shp")

# load waterholes or point features
wh = st_read("data/waterholes/waterpoints_zambatis_2011.shp")
# filter waterholes by the extent
wh = filter(wh, CURRENT=="Open")

ext = st_read("data/clip/kruger_clip.shp")

data_sf = st_as_sf(ele, coords = c("location.long", "location.lat")) %>%
  `st_crs<-`(32736)

ele = ele %>% 
  select(id = event.id, long = location.long, lat = location.lat, 
         temp = external.temperature, 
         time = timestamp)

ele$time <- as.POSIXct(strptime(ele$time, format='%Y-%m-%d %H:%M:%S'))


library(move)
library(data.table)


data <- getDataRepositoryData("doi:10.5441/001/1.403h24q5")
save(data, file = "data/elephant_data.Rdata")
data_coords <- split(data)

# get data
data_coords <- Map(function(le, tag_id) {
  dt <- data.table(
    cbind(
      coordinates(le),
      timestamps(le)
    ),
    tag_id
  )
  setnames(dt, c("x", "y", "time", "id"))
}, data_coords, names(data_coords))

rm(data)
gc()

geometry <- st_sfc(
  lapply(data_coords, function(x) {
    st_linestring(
      as.matrix(x[, c("x", "y")])
    )
  }),
  crs = 4326
)
geometry <- st_transform(geometry, 32736)
data_sf <- mapply(
  function(df) {
    df[1, c("id")]
  },
  data_coords,
  SIMPLIFY = FALSE
)
# add geometry
data_sf <- rbindlist(data_sf)
data_sf[, geometry := geometry]
# make sf
data_sf <- st_sf(data_sf, crs = 32736)
# save
st_write(data_sf,
         dsn = "data/data_lines_elephants.gpkg",
         append = FALSE
)
dataz <- st_read("data/data_lines_elephants.gpkg")
data_sf = st_as_sf(dataz, coords = c("xutm", "yutm")) %>%
  `st_crs<-`(32736)
# get distance to waterholes and rivers
distwh = as.numeric(st_distance(data_sf, wh))



#################
###############
# prepare bounding box
bbox <- c(
  xmin = 330000,
  xmax = 393000,
  ymin = 7260000,
  ymax = 7298050
)
bbox_sf <- st_bbox(bbox)
bbox_sf <- st_as_sfc(bbox_sf)
st_crs(bbox_sf) <- 32736
##
kruger <- st_read("data/clip/kruger_clip.shp")
kruger <- st_transform(kruger, 32736)
# get inversion
kruger_invert <- st_difference(
  st_as_sfc(st_bbox(kruger)),
  kruger
)
# get kruger point -- this is hardcoded but could also be a centroid
kruger_point <- st_point(c(31.5, -24))
kruger_point <- st_sfc(kruger_point, crs = 4326)
kruger_point <- st_transform(kruger_point, 32736)

# get africa for inset
africa <- st_read("data/africa.gpkg")
africa <- st_transform(africa, 32736)

rivers <- st_transform(rivers[is.na(rivers$seasonal), ], 32736)

river_season = list(distr_seasonal = rivers %>% filter(is.na(seasonal)),
                    distr = (rivers)) %>%
  map(st_union)

distr_seasonal = map(river_season, function(l){
  return(as.numeric(st_distance(data_sf, l)))
})


# add distance to rivers data
data = mutate(data,
              distr = distr_seasonal[["distr"]],
              distr_s = distr_seasonal[["distr_seasonal"]],
              distwh = distwh)

# calculate mindw, the minimum distance to water
data = data %>%
  mutate(mindw = case_when(
    season =="dry"~ifelse(distr_s < distwh, distr_s, distwh),
    season == "wet"~ifelse(distr < distwh, distr, distwh),
    T~as.double(NA)
  ))

drop_geom <- function(x)
{
  if(inherits(x,"sf"))
    ret <- x[,setdiff(names(x),attr(x,'sf_column')),drop=T]
  else
    ret <- x
  
  class(ret) <- 'data.frame'
  return(ret)
}

str(x <- st_sf(a=3:4,b=5:6, geom=st_sfc(st_point(1:2),st_point(3:4))))
str(df <- drop_geom(data_sf))

df <- st_drop_geometry(data_sf)

library(sfheaders)
df <- sf_to_df(data_sf) # to convert sf to data frame

ele$utm_long <- df$x
ele$utm_lat <- df$y

distwh = as.numeric(st_distance(ele, wh))

data_sf = st_as_sf(ele, coords = c("utm_long", "utm_lat")) %>%
  `st_crs<-`(32736)

# get distance to waterholes and rivers
distwh = as.numeric(st_distance(data_sf, wh))

data = mutate(data_sf,
              distwh = distwh)

###################
###riversääääääääää
library(rnaturalearth)
library(osmdata)

# if data does not already exist
if (!file.exists("data/rivers_kruger.gpkg")) {
  # kruger bounding box
  kruger <- st_read("data/clip/kruger_clip.shp")
  q <- opq(bbox = st_bbox(kruger))
  
  # make query
  query_waterways <- add_osm_feature(q,
                                     key = "waterway",
                                     value = c("river", "stream")
  )
  
  # run query
  rivers_kruger <- osmdata_sf(query_waterways)
  
  # get only lines
  rivers_kruger <- rivers_kruger$osm_lines
  
  # assign crs
  st_crs(rivers_kruger) <- 4326
  
  st_write(
    rivers_kruger,
    "data/rivers_kruger.gpkg"
  )
}

river_season = list(distr_seasonal = rivers %>% filter(is.na(seasonal)),
                    distr = (rivers)) %>%
  map(st_union)

distr_seasonal = map(river_season, function(l){
  return(as.numeric(st_distance(data_sf, l)))
})
####
data = mutate(ele,
              distr = distr_seasonal[["distr"]],
              distr_s = distr_seasonal[["distr_seasonal"]],
              distwh = distwh)


index <- st_nearest_feature(x = data_sf, y = wh)
wh <- wh %>% slice(index)
distwhs <- st_distance(x = data_sf, y= wh, by_element = TRUE)
data_sf$distwh <- distwhs
# works!