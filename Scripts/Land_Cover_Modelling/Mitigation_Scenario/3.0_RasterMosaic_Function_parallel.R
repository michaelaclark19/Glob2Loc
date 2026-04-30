#!/usr/bin/env Rscript

# Read me ----
# Mosaicing the rasters created in "ConversionScript"
# This is for Business-as-usual, change the scenario on L29, 134, 142

# Packages ---
library(raster)
library(data.table)
library(doMC)
library(plyr)
library(dplyr)
rasterOptions(maxmemory = 1e+11)

# Setting working directory
setwd('/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity')

# Getting 2010 crop ----
crop_2010 <- raster(paste0(getwd(),"/ESA LandCov Maps/Crop2010_Corrected_GlobCov.tif"))
pasture_2010 <- raster(paste0(getwd(),"/ESA LandCov Maps/Pasture2010_Corrected_GlobCov.tif"))
crop_2010_vector <- getValues(crop_2010)

# Getting row index vectors [DON'T CHANGE THESE FOR EACH SCENARIO] -----
# Rasters for the row indices of individual countries
row_index_rasters <- list.files(path = paste0(getwd(),"/Other Data Inputs/Row_Index_New"),
                                pattern = "Index",
                                full.names = TRUE) %>%
  .[grep('.RData',.)]

# Get a global raster of row indices

# Setting up data table
world_dt <- data.table(
  cell = 1:289285431, # Number of cells in global raster
  crop_orig = rep(-1, 289285431))

world_dt$crop_new <- world_dt$crop_orig

setkey(world_dt, cell)

# Country raster
country_id_raster_global <- raster(paste0(getwd(),'/Global Mollweide Maps/MollweideCountryID_1.5km.tif'))


# Run across multiple cores  -----
# Creating loop around the parallel function
# This loop just subsets rasters to feed into the function
land.covs <- c("pasture.*mean",'pasture.*sd',
               'crop.*mean','crop.*sd')

# List of years
years <- paste(c(seq(2010, 2045, by = 5)),
               c(seq(2015, 2050, by = 5)),
               sep = "_")

registerDoMC(cores = length(years))

# list of scenarios
directories <- list.files(path = paste0(getwd(),"/Land Forecast Outputs/"),pattern = 'EAT')

# Function to moasic the reasters together
registerDoMC(cores = length(years))

# Function to moasic the reasters together

# Function to run -----
mosaic_raster_fun <- function(year){
  
  
  crop_rasters <- raster_list[grep(year, raster_list)]
  mean_sd <- if(grepl("mean", crop_rasters[1])) "mean" else "sd"
  # Create a copy of the data.table, so we don't have to re-make it
  dt_copy <- copy(world_dt)
  
  # Cycle through the crop rasters...
  for(i in 1:length(crop_rasters)){
    print(i) # Print progress
    
    # Getting row index raster
    # This loads a raster named tmp.country
    # Need an exception for the countries that span world regions (e.g. OAO, etc)
    if(length(grep(substring(gsub("^.*/tifs/", "",
                                  crop_rasters[[i]]),
                             first = 1, last =3),
                   row_index_rasters)) > 1) {
      # List of all rasters that meet the country ISO3
      tmp.load <- row_index_rasters[grep(substring(gsub("^.*/tifs/", "",
                                                        crop_rasters[[i]]),
                                                   first = 1, last =3),
                                         row_index_rasters)]
      # Getting region of our rasters
      tmp.region <- (substring(gsub("^.*/tifs/", "",
                                    crop_rasters[[i]]),
                               first = 4, last =8))
      # And subsetting to the raster in the correct region
      tmp.load <- tmp.load[grep(tmp.region,tmp.load)]
      # Loading
      load(tmp.load)
    } else {
      load(row_index_rasters[grep(substring(gsub("^.*/tifs/", "",
                                                 crop_rasters[[i]]),
                                            first = 1, last =3),
                                  row_index_rasters)])
    }
    
    # Create a data.table with the proportion of cropland and the row indexes.
    tmp_dt <- data.table(crop = getValues(raster(crop_rasters[[i]])),
                         row_index = getValues(tmp.country))
    # Hopefully not needed, but better safe than sorry
    tmp_dt <- tmp_dt[order(tmp_dt$row_index)]
    
    dt_copy[tmp_dt[!is.na(tmp_dt$crop), row_index],
            # Returns the row_index of the cells in the temporary data.table 
            #     where the cropland cell is not NA
            crop_new := tmp_dt[!is.na(tmp_dt$crop), crop]]
    # the ":=" syntax is some madcap binary search thing which is lightning
  }
  
  rm(tmp.country)
  
  # Making raster
  tmp_new <- raster(matrix(dt_copy$crop_new,
                           nrow = nrow(crop_2010),
                           byrow = TRUE), crs = crs(crop_2010))
  # Changing extent
  extent(tmp_new) <- extent(crop_2010)
  # Getting rid of oceans...
  tmp_new[tmp_new %in% -1] <- NA
  #tmp_new[is.na(tmp_new) & !is.na(crop_2010)] <- crop_2010[is.na(tmp_new) & !is.na(crop_2010)]


  # if else statement to make sure we still have all areas covered
  if(grepl('crop',land)) {
	  if(grepl('sd',land)) {tmp_new[is.na(tmp_new) & !is.na(crop_2010)]<-0}
	  if(!grepl('sd',land)) {tmp_new[is.na(tmp_new) & !is.na(crop_2010)]<-crop_2010[is.na(tmp_new) & !is.na(crop_2010)]}
  
  } else {

	  if(grepl('sd',land)) {tmp_new[is.na(tmp_new) & !is.na(pasture_2010)]<-0}
          if(!grepl('sd',land)) {tmp_new[is.na(tmp_new) & !is.na(pasture_2010)]<-pasture_2010[is.na(tmp_new) & !is.na(pasture_2010)]}
  
  } # End if else statement to ensure all areas are covered
  
  # return(tmp_new)
  # And saving the raster
  writeRaster(tmp_new,
              paste0(getwd(),"/Land Forecast Outputs/",directory,"/Global tifs/",gsub("\\.\\*","_",land),year,".tif"),
              overwrite = TRUE)
}


# And running the script
for(directory in directories) {
  # Creating directory
  dir.create(paste0(getwd(),"/Land Forecast Outputs/",directory,"/Global tifs"))
  
  # Getting full list of files
  raster.files <- list.files(paste0(getwd(),"/Land Forecast Outputs/",directory,"/tifs"), full.names = TRUE)
  
  # Looping through land cover types
  for(land in land.covs) {
    # Subsetting based on land cover
    raster_list = raster.files[grepl(land,raster.files)]
    
    # Doing this in a loop across years
    # For whatever reason this isn't working in parallel with either mclapply or foreach
    # This isn't ideal, but the process doesn't take too long
    # So probably not a huge issue
    for(y in 1:length(years)) {
      mosaic_raster_fun(year = years[y])
    }
    
    # And now doing this in parallel
    # Note that this is processor intensive, so be a bit careful on how much you parallelize this
    # trial.2 <- mclapply(years[1:2], FUN = raster_fun)
    # 
    # trial.3 = foreach(j = 1:length(years)) %dopar% raster_fun2(year = years)
  }
}
