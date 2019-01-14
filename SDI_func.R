# This script contains the tools used to construct the landscape-agriculture dataset:

############################################################################################################################
# Extraction from CDL
############################################################################################################################

# Slice national CdL rasters into smaller state rasters for use in velox
write_states <- function(year = 2008) {
  foreach(i = 1:length(unique(states$STATE_FIPS)),
          .export = c("states", "cdl"),
          .packages = c("sf", "rgeos", "raster")) %dopar% {
            stateFIPS <- unique(states$STATE_FIPS)
            state <- states[states$STATE_FIPS == stateFIPS[i],]
            state <- sf::st_as_sf(rgeos::gBuffer(as(state, "Spatial"), width = 50000, byid=T))  # add 50K buffer around the state 
            raster::crop(cdl, state, filename = paste0("C:/Users/A02256433/Desktop/Data/CDL/States/", year, "_", state$STATE_NAME, ".gri"))
          }
}

# Extract indices from CDL using index_function to county shapefiles using velox
extract_CDL_values <- function(index_function) {
  foreach(i = 1:length(all_states),
          .combine = rbind,
          .export = c("all_states", "counties", "states", "mask"),
          .packages = c("velox")) %dopar% {
            year <- as.numeric(substr(all_states[i], 44, 47)) 
            state_name <- gsub(".grd", "", sub(paste0(".*", year,"_"), "", all_states[i]))
            cty <- counties[counties$STATEFP == states$STATE_FIPS[states$STATE_NAME == state_name],]
            ras <- velox(all_states[i])
            ev <- ras$extract(sp = cty, fun = index_function, legacy=T) 
            remove(ras)
            cbind.data.frame(rep(year, nrow(ev)), ev, cty$GEOID)
          }
}

# Extract indices from CDL using index_function to grid (for visualization purposes)
extract_CDL_to_grd <- function(index_function) {
  foreach(i = 1:length(all_states),
          .combine = rbind,
          .export = c("all_states", "counties", "states", "mask", "grd"),
          .packages = c("velox")) %dopar% {
            year <- as.numeric(substr(all_states[i], 44, 47)) 
            state_name <- gsub(".grd", "", sub(paste0(".*", year,"_"), "", all_states[i]))
            cty <- counties[counties$STATEFP == states$STATE_FIPS[states$STATE_NAME == state_name],]
            ras <- velox(all_states[i])
            ev <- ras$extract(sp = grd, fun = index_function, legacy=T) 
            remove(ras, grds)
            cbind.data.frame(ev, rep(year, nrow(ev)), grd$ID)
          }
}

############################################################################################################################
# Landscape indices (include mask option for extraction)
############################################################################################################################

# Shannon Diversity Index
sdi <- function(ext) {
  ls <- 0
  ext[ext == 0] <- NA  # In CDL, zero is "background"/water in coastal areas
  ext[ext %in% mask] <- NA # Be sure to add the mask to future indices
  ext <- ext[complete.cases(ext)] # remove NAs, focus only on land in ROI
  
  for (z in unique(ext)) {
    p_i <- length(ext[ext == z])/length(ext) 
    y <- -(p_i * log(p_i))
    ls <- ls + y
  }
  return(ls)
  remove(ext, y, p_i)
}

# Simpson Diversity Index

#ext <- c(rep(1,2), rep(2, 8), 3, 4, rep(5, 3))
sidi  <- function(ext) {
  # tested with this: https://geographyfieldwork.com/SimpsonsDiversityIndex.htm
  nsum <- 0
  nsum2 <- 0
  ext[ext == 0] <- NA
  ext[ext %in% mask] <- NA
  for (z in unique(ext)) {
    n <- length(ext[ext==z])
    nsum <- nsum+n
    n2 <- n*(n-1)
    nsum2 <- nsum2 + n2
  }
  out <- 1 - (nsum2/(nsum*(nsum-1)))
  return(out)
}

# Richness

rich  <- function(ext) {
  ext[ext == 0] <- NA
  ext[ext %in% mask] <- NA
  out <- length(unique(ext))
  return(out)
}


#################################################################################################################################
# Extraction from PRISM
#################################################################################################################################

extract_PRISM_values <- function() {
  foreach(i = 1:length(prism_list),
          .combine = rbind,
          .export = c("counties", "start", "prism_list"),
          .packages = c("velox")) %dopar%  {
            ev <- velox(prism_list[i])$extract(sp = counties, fun = function(x) mean(x, na.rm = TRUE), legacy=T) 
            year <- as.numeric(substr(prism_list[i], start, start+3))
            month <- as.numeric(substr(prism_list[i], start+4, start+5))
            day <- as.numeric(substr(prism_list[i], start+6, start+7))
            cbind.data.frame(counties$GEOID, rep(year, nrow(ev)), rep(month, nrow(ev)), rep(day, nrow(ev)), ev)
          }
}

#################################################################################################################################
# Climate Indices
#################################################################################################################################

# Subset extracted daily PRISM data for a crop-specific growing season
time_subset <- function(df, crop)  {
  df <- df %>% filter(YEAR %in% 2008:2017)  # filter years CDL data available
  sta_vn <- paste0(crop, ".plant.start")
  sto_vn <- paste0(crop, ".harvest.end")
  
  if (crop == "wwheat") {
     df[,"PPT"][df$DOY >= df[,sto_vn] & df$DOY <= df[,sta_vn]] <- NA 
     df[,"TMAX"][df$DOY >= df[,sto_vn] & df$DOY <= df[,sta_vn]] <- NA 
   } else {
     df[,"PPT"][df$DOY <= df[,sta_vn]] <- NA
     df[,"PPT"][df$DOY >= df[,sto_vn]] <- NA
     df[,"TMAX"][df$DOY <= df[,sta_vn]] <- NA
     df[,"TMAX"][df$DOY >= df[,sto_vn]] <- NA
   }
  
  return(df)
}

# Compute growing degree days (GDD) for daily climate data masked for a specific crop
gdd <- function(df, crop, gdd_range) {
  tmax <- df$TMAX
  gdd <- ifelse(tmax <= gdd_range[[crop]][1], 0, tmax)
  gdd <- ifelse(tmax >= gdd_range[[crop]][1] & tmax <= gdd_range[[crop]][2], tmax - gdd_range[[crop]][1], gdd)
  gdd <- ifelse(tmax >= gdd_range[[crop]][2], gdd_range[[crop]][2] - gdd_range[[crop]][1], gdd)
  return(gdd)
}

# Compute stress degree days (SDD) for daily climate data masked for a specific crop
sdd <- function(df, crop, gdd_range) {
  tmax <- df$TMAX
  sdd <- ifelse(tmax >= gdd_range[[crop]][2], tmax - gdd_range[[crop]][2], 0)
  return(sdd)
}

# Compute effective precipitation (EfP) for daily climate data masked for a specific crop
efp <- function(df, crop, th) {
  rain <- df$PPT
  efp <- ifelse(rain >= th, th, rain)
  return(efp)
}

wheat_subset <- function(df) {
  season <- df %>% mutate(YEAR2 = YEAR, YEAR = NA) 
  season$YEAR[season$DOY >= season$wwheat.plant.start] <- season$YEAR2[season$DOY >= season$wwheat.plant.start]
  season$YEAR[season$DOY <= season$wwheat.harvest.end] <- season$YEAR2[season$DOY <= season$wwheat.harvest.end] - 1
  season <- season %>% filter(!is.na(season$YEAR))  # drops DOYs between stop/start where YEAR remains NA
  
  season <- season %>% group_by(YEAR, GEOID) %>%
    summarize(GDD = sum(gdd, na.rm=T), SDD = sum(sdd, na.rm=T), 
              EfP = sum(efp, na.rm = T), TP = sum(PPT, na.rm=T)) %>%
    mutate(ExP = TP - EfP) 
  comment(season) <- crop
  season <- season %>% filter(YEAR > 2007)
  # only have the season available for these years
  season$GDD[season$YEAR == 2017] <- NA
  season$SDD[season$YEAR == 2017] <- NA
  season$EfP[season$YEAR == 2017] <- NA
  season$TP[season$YEAR == 2017] <- NA
  season$ExP[season$YEAR == 2017] <- NA
  
  return(season)
}

yield_merge <- function(df, crop) {
  y <- readRDS(paste0("./out/yield/", crop, ".RDS"))
  dfm <- merge(df, y, by = c("GEOID", "YEAR"), all = T)
}

seasonal_panels <- function(prism, crop) {
  
  x <- time_subset(prism, crop)
  x$gdd <- gdd(x, crop, gdd_range[crop]) 
  x$sdd <- sdd(x, crop, gdd_range[crop])
  x$efp <- efp(x, crop, th = 20)
  
  if (crop == "wwheat") {
    season <- wheat_subset(x)
  } else {
    season <- x %>% group_by(YEAR, GEOID) %>%
      summarize(GDD = sum(gdd, na.rm=T), SDD = sum(sdd, na.rm=T), 
                EfP = sum(efp, na.rm = T), TP = sum(PPT, na.rm=T)) %>%
      mutate(ExP = TP - EfP)
  }
  
  season <- yield_merge(season, crop)
  comment(season) <- crop  
  
  sdi <- readRDS("./out/indices/sdi_agmask.RDS")
  colnames(sdi) <- c("YEAR", "SDI_AG", "GEOID")
  sdina <- readRDS("./out/indices/sdi_nomask.RDS")
  colnames(sdina) <- c("YEAR", "SDI_ALL", "GEOID")
  sidi <- readRDS("./out/indices/sidi_agmask.RDS")
  colnames(sidi) <- c("YEAR", "SIDI_AG", "GEOID")
  sidina <- readRDS("./out/indices/sidi_nomask.RDS")
  colnames(sidina) <- c("YEAR", "SIDI_ALL", "GEOID")
  rich <- readRDS("./out/indices/rich_agmask.RDS")
  colnames(rich) <- c("YEAR", "RICH_AG", "GEOID")
  richna <- readRDS("./out/indices/rich_nomask.RDS")
  colnames(richna) <- c("YEAR", "RICH_ALL", "GEOID")
  ir <- readRDS("./out/irrigation.RDS")
  irrig <- do.call("rbind", replicate(length(unique(sdi$YEAR)), ir, simplify = FALSE))
  irrig$YEAR <- rep(unique(sdi$YEAR), each = nrow(ir))
  season <- merge(season, sdina, by = c("GEOID", "YEAR"), all.x=T)
  season <- merge(season, sdi, by = c("GEOID", "YEAR"), all.x=T)
  season <- merge(season, sidina, by = c("GEOID", "YEAR"), all.x=T)
  season <- merge(season, sidi, by = c("GEOID", "YEAR"), all.x=T)
  season <- merge(season, irrig, by = c("GEOID", "YEAR"), all.x = T)
  season <- merge(season, rich, by = c("GEOID", "YEAR"), all.x=T)
  season <- merge(season, richna, by = c("GEOID", "YEAR"), all.x = T)
  return(season)
}


