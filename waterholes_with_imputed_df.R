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

library(move)
library(data.table)

library(tidyverse)
#################################################################################
setwd("C:/Users/timom/")
datazz <- read.table(
  "imputed_data_lonlat_without_waterholes.csv",
  header = TRUE,
  sep= ","
)
df <- st_as_sf(x = datazz,                         
               coords = c("x", "y"),
               crs = "+proj=longlat +datum=WGS84 +ellps=WGS84")
geometry <- st_transform(df, 32736)###################
sfc_df <- as.data.frame(geometry)
library(tidyverse)

separated_coord <- sfc_df %>%
  mutate(x = unlist(map(sfc_df$geometry,1)),
         y = unlist(map(sfc_df$geometry,2)))


write.csv(separated_coord,"imputed_data_UTM_without_waterholes.csv", row.names = FALSE)
#################################################################################


df_sf <- st_as_sf(x = separated_coord,                         
               coords = c("x", "y"),
               crs = 32736)

setwd("C:/Users/timom/anaconda3/envs/madeleine_project")
# Waterholes
wh = st_read("data/waterholes/waterpoints_zambatis_2011.shp")
# filter waterholes by the extent
wh = filter(wh, CURRENT=="Open")

# Distance
index <- st_nearest_feature(x = df_sf, y = wh)
wh <- wh %>% slice(index)
distwhs <- st_distance(x = df_sf, y= wh, by_element = TRUE)
#data_sf$distwh <- distwhs

#############
### boundaries###
#############
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

#################
## Rivers ###
#################
####
# rivers 2
#####
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

rivers <- st_read("data/rivers_kruger.gpkg")

rivers <- st_transform(rivers[is.na(rivers$seasonal), ], 32736)
river_season = list(distr_seasonal = rivers %>% filter(is.na(seasonal)),
                    distr = (rivers)) %>%
  map(st_union)

distr_seasonal = map(river_season, function(l){
  return(as.numeric(st_distance(df_sf, l)))
})


# add distance to rivers data
data = mutate(df_sf,
              distr = distr_seasonal[["distr"]],
              distr_s = distr_seasonal[["distr_seasonal"]],
              distwhs = distwhs)

######
# Seasons
####
data$season <- data$month
data$season <- ifelse(data$month > 4 & data$month < 11,data$season <- "dry", data$season <- "wet")
  
data$distwhs <- as.numeric(data$distwhs)
## minimum distance
# calculate mindw, the minimum distance to water
data$mindw <- data$month
data = data %>%
  mutate(mindw = case_when(
    season =="dry"~ifelse(distr_s < distwhs, distr_s, distwhs),
    season == "wet"~ifelse(distr < distwhs, distr, distwhs),
    T~as.double(NA)
  ))

boxplot(data$mindw)

# get rid of geometry

dataz <- data %>%
  mutate(x = unlist(map(data$geometry,1)),
         y = unlist(map(data$geometry,2)))

dataz <- as.data.frame(dataz)

dataz <- dataz[, -c(19)]

write.csv(dataz,"imputed_dataframe_with_rivers_waterholes_and_distance.csv", row.names = FALSE)
