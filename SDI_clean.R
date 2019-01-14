# Running this script will reconstruct all data necessary to create the landscape-agriculture panels.  
# Computationally-intensive code chunks are commented-out.
# You'll need to first load all data described in SDI_load.R
# For clarity, I've loaded all package::functions required by each function below directly above that function.

source("SDI_func.R")
source("SDI_load.R")

############################################################################################################################
# Split large national CdL rasters into state-sized rasters to facilitate use of velox
############################################################################################################################

library(sf)
library(raster)
library(velox)
library(rgeos) # for buffer
library(parallel)

years <- 2008:2017
for (i in 1:length(years)) {
  cdl <- raster(cdl_list[i])
  library(doParallel)
  num_cores <- detectCores() - 3
  clus <- makeCluster(num_cores)
  registerDoParallel(clus)
  write_states(year = years[i])
  stopCluster(clus)
}

############################################################################################################################
# Extract landscape indices to shapefiles
############################################################################################################################

library(sf)
library(velox)
library(doParallel)

# Load list of state.grd files:
all_states <- Sys.glob(paste0("C:/Users/A02256433/Desktop/Data/CDL/States/", "*", ".grd"))

# Specify a mask (agriculture below)
mask <- c(37, 59, 60, 63:65, 81:83, 87:88, 92, 111:112, 121:124, 131, 141:143, 152, 176, 190, 195)
# mask <- c()

# Select an index_function() that includes masking, e.g. sdi()
ptm <- proc.time()
num_cores <- detectCores() - 4
clus <- makeCluster(num_cores)
registerDoParallel(clus)
out_data <- extract_CDL_values(index_function = sidi)
stopCluster(clus)
colnames(out_data) <- c("YEAR", "SIDI_AG", "GEOID")
saveRDS(out_data, "./out/indices/sidi_agmask2.RDS")
proc.time() - ptm

############################################################################################################################
# Extract climate indices
############################################################################################################################

library(velox)
library(doParallel)
library(sf)

# be sure VOI matches prism_list variable
prism_list <- Sys.glob(paste0("C:/Users/A02256433/Desktop/Data/PRISM/daily_ppt_data/", "*", ".bil"))
counties <- st_transform(counties, 4269)  # project counties to PRISM projection (unprojected)
voi <- "PPT"
# start <- 79 # tmax, I was lazy with regex, so also need to input a "start" index for date values
start <- 77  # ppt

ptm <- proc.time()
num_cores <- detectCores() - 1
clus <- makeCluster(num_cores)
registerDoParallel(clus)
out_data <- extract_PRISM_values()
stopCluster(clus)
colnames(out_data) <- c("GEOID", "Year", "Month", "Day", voi)
saveRDS(out_data, "daily_county_PPT.RDS")
proc.time() - ptm

# checks
tmax <- readRDS("./out/prism/daily_county_TMAX.RDS")
ppt <- readRDS("./out/prism/daily_county_PPT.RDS")
cnt <- ppt %>% group_by(Year, GEOID) %>% summarize(n=n()) #355, 366
nas <- tmax[is.na(tmax$TMAX),] # Falls Church, VA, Manassas Park, VA, Fairfax, VA, Lexington, VA
nasp <- ppt[is.na(ppt$PPT),] # idem
# re-ran extract and didn't work... odd that these polygons don't work.


#############################################################################################################################
# Construct season stop/start 
#############################################################################################################################
 
corn_plant <- raster("./data/corn.plant.start.asc")
corn_harvest <- raster("./data/corn.harvest.end.asc")
soy_plant <- raster("./data/soy.plant.start.asc")
soy_harvest <- raster("./data/soy.harvest.end.asc")
wheat_plant <- raster("./data/wheat.plant.start.asc")
wheat_harvest <- raster("./data/wheat.harvest.end.asc")
wwheat_plant <- raster("./data/wwheat.plant.start.asc")
wwheat_harvest <- raster("./data/wwheat.harvest.end.asc")
stop_start <- brick(corn_plant, corn_harvest,
                    soy_plant, soy_harvest,
                    wheat_plant, wheat_harvest,
                    wwheat_plant, wwheat_harvest)
stop_start@crs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
cty_prj <- spTransform(counties, stop_start@crs)
coi <- extract(stop_start, cty_prj, fun=mean, sp=T, na.rm=T)
coi_data <- coi@data %>% select(GEOID, NAME, corn.plant.start, corn.harvest.end,
                                soy.plant.start, soy.harvest.end,
                                wheat.plant.start, wheat.harvest.end,
                                wwheat.plant.start, wwheat.harvest.end)
coi@data <- coi_data
saveRDS(coi, "./out/stopstart.RDS")

#############################################################################################################################
# Yield
#############################################################################################################################

library(stringr)

corn <- read.csv("./data/yield/corn.csv")
corn_a <- corn %>%
  filter(Data.Item == "CORN, GRAIN - ACRES HARVESTED") %>%
  filter(!is.na(County.ANSI)) %>%   # drops anything not at the county-scale (district, national)
  mutate(STATE = str_pad(as.character(State.ANSI), 2, side = "left", pad = "0"),
         COUNTY = str_pad(as.character(County.ANSI), 3, side = "left", pad = "0"),
         GEOID = as.factor(paste0(STATE,COUNTY)),
         YEAR = Year,
         ACRES = Value) %>%
      select(YEAR, GEOID, ACRES)
corn_y <- corn %>%
  filter(Data.Item == "CORN, GRAIN - YIELD, MEASURED IN BU / ACRE") %>%
  filter(!is.na(County.ANSI)) %>%
  mutate(STATE = str_pad(as.character(State.ANSI), 2, side = "left", pad = "0"),
         COUNTY = str_pad(as.character(County.ANSI), 3, side = "left", pad = "0"),
         GEOID = as.factor(paste0(STATE,COUNTY)),
         YIELD = Value,
         YEAR = Year) %>%
  select(YEAR, GEOID, YIELD)

soy <- read.csv("./data/yield/soy.csv")
soy_a <- soy %>%
  filter(Data.Item == "SOYBEANS - ACRES HARVESTED") %>%
  filter(!is.na(County.ANSI)) %>%
  mutate(STATE = str_pad(as.character(State.ANSI), 2, side = "left", pad = "0"),
         COUNTY = str_pad(as.character(County.ANSI), 3, side = "left", pad = "0"),
         GEOID = as.factor(paste0(STATE,COUNTY)),
         ACRES = Value,
         YEAR = Year) %>%
  select(YEAR, GEOID, ACRES)

soy_y <- soy %>%
  filter(Data.Item == "SOYBEANS - YIELD, MEASURED IN BU / ACRE") %>%
  filter(!is.na(County.ANSI)) %>%
  mutate(STATE = str_pad(as.character(State.ANSI), 2, side = "left", pad = "0"),
         COUNTY = str_pad(as.character(County.ANSI), 3, side = "left", pad = "0"),
         GEOID = as.factor(paste0(STATE,COUNTY)),
         YIELD = Value,
         YEAR = Year) %>%
  select(YEAR, GEOID, YIELD)

wwheat <- read.csv("./data/yield/wwheat.csv")

wwheat_a <- wwheat %>%
  filter(Data.Item == "WHEAT, WINTER - ACRES HARVESTED") %>%
  filter(!is.na(County.ANSI)) %>%
  mutate(STATE = str_pad(as.character(State.ANSI), 2, side = "left", pad = "0"),
         COUNTY = str_pad(as.character(County.ANSI), 3, side = "left", pad = "0"),
         GEOID = as.factor(paste0(STATE,COUNTY)),
         ACRES = Value,
         YEAR = Year) %>%
  select(YEAR, GEOID, ACRES)

wwheat_y <- wwheat %>%
  filter(Data.Item == "WHEAT, WINTER - YIELD, MEASURED IN BU / ACRE") %>%
  filter(!is.na(County.ANSI)) %>%
  mutate(STATE = str_pad(as.character(State.ANSI), 2, side = "left", pad = "0"),
         COUNTY = str_pad(as.character(County.ANSI), 3, side = "left", pad = "0"),
         GEOID = as.factor(paste0(STATE,COUNTY)),
         YIELD = Value,
         YEAR = Year) %>%
  select(YEAR, GEOID, YIELD)

soym <- merge(soy_y, soy_a, by = c("GEOID", "YEAR"), all=T)
cornm <- merge(corn_y, corn_a, by = c("GEOID", "YEAR"), all=T)
wwheatm <- merge(wwheat_y, wwheat_a, by = c("GEOID", "YEAR"), all=T)
saveRDS(soym, "./out/yield/soy.RDS")
saveRDS(cornm, "./out/yield/corn.RDS")
saveRDS(wwheatm, "./out/yield/wwheat.RDS")


#################################################################################################################################
# Extract MIRAD to county, percent agricultural lands irrigated in a county
#################################################################################################################################

library(velox)
library(sf)
library(raster) # total area 

# 2012 as "baseline"
mirad <- velox("C:/Users/A02256433/Desktop/Data/MIRAD/mirad250_12v3_bsq/mirad250_12v3envi.bsq", 
               crs="+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
mirad <- velox("/Users/emily.burchfield/Downloads/mirad250_12v3_bsq/mirad250_12v3envi.bsq")
counties <- st_transform(counties, "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
sumirr <- mirad$extract(sp = counties, fun = function(x) sum(x, na.rm = TRUE), legacy = T)
totalarea <- mirad$extract(sp = counties, fun = function(x) length(x[!is.na(x)]), legacy=T)
irr <- cbind.data.frame(sumirr, totalarea, counties$GEOID)
colnames(irr) <- c("irr", "total", "GEOID")
irr$PERC_IRR <- (irr$irr/irr$total)*100

# total area in acres of counties 
cty <- shapefile("./data/county.RDS")
cty$AREA_SQKM <- area(cty)/1000000  # checked comparing max/min
cty$AREA_ACRES <- cty$AREA_SQKM*247.105

irr <- merge(cty@data, irr, by = "GEOID", all.x=T)
irr <- irr %>% dplyr::select(GEOID, PERC_IRR, AREA_ACRES)  # checked merge
colnames(irr) <- c("GEOID", "PERC_IRR", "TOTAL_AREA_ACRES")

saveRDS(irr, "./out/irrigation.RDS")


#################################################################################################################################
# Build seasonal parameters
#################################################################################################################################

library(lubridate)
library(tidyverse)

ss <- readRDS("./out/stopstart.RDS")
# tmax <- readRDS("./out/prism/daily_county_TMAX.RDS")
# ppt <- readRDS("./out/prism/daily_county_PPT.RDS")
# 
# tmax$DOY <- yday(as.Date(paste(tmax$Year,tmax$Month,tmax$Day,sep="-"))) 
# ppt$DOY <- yday(as.Date(paste(ppt$Year,ppt$Month,ppt$Day,sep="-"))) 
# 
# tmax <- tmax %>% select(GEOID, Year, DOY, TMAX) %>% arrange(Year, DOY, GEOID)
# colnames(tmax) <- c("GEOID", "YEAR", "DOY", "TMAX") # fix lowercase Year
# ppt <- ppt %>% select(GEOID, Year, DOY, PPT) %>% arrange(Year, DOY, GEOID)
# prism <- cbind(tmax, ppt %>% select(PPT))
# prism <- merge(prism, stop_start@data, by.x = "GEOID", by.y = "GEOID", all.x = T)
# saveRDS(prism, "./out/prism/prism.RDS")
prism <- readRDS("./out/prism/prism.RDS")

crop <- "corn"
corn_panel <- seasonal_panels(prism, "corn")
corn_panel$CROP <- "corn"
crop <- "soy"
soy_panel <- seasonal_panels(prism, "soy")
soy_panel$CROP <- "soy"
crop <- "wwheat"
wwheat_panel <- seasonal_panels(prism, "wwheat")
wwheat_panel$CROP <- "wwheat"

final_panel <- rbind(corn_panel, soy_panel, wwheat_panel)
saveRDS(final_panel, "./out/full_panel.RDS")






