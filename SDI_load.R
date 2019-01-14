library(sf)
library(sp)

# List of national extents CDL rasters downloaded to machine:
cdl_list <- Sys.glob(paste0("C:/Users/A02256433/Desktop/Data/CDL/", "*", ".img"))

# Shapefile of counties in the contiguous US
counties <- st_read("./data/county.shp") 
counties <- st_transform(counties, 5070)

# Shapefile of states in the contiguous US
states <- st_read("./data/states.shp", stringsAsFactors = F)
states <- st_transform(states, 5070) 
states <- states[states$STATE_FIPS %in% unique(counties$STATEFP),] 

# List of daily PRISM data downloaded to machine:
prism_tmax_list <- Sys.glob(paste0("C:/Users/A02256433/Desktop/Data/PRISM/daily_tmax_data/", "*", ".bil"))
 
# gdd range, from p. 8 in Neil's thesis
gdd_range <- list()
gdd_range[["corn"]] <- c(10,30)
gdd_range[["soy"]] <- c(10, 30)
gdd_range[["wheat"]] <- c(0,35)
gdd_range[["wwheat"]] <- c(0, 30)
gdd_range[["cotton"]] <- c(15.6, 37.8)

# county season stop-start data
stop_start <- readRDS("./out/stopstart.RDS")
stop_start <- spTransform(stop_start, CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

