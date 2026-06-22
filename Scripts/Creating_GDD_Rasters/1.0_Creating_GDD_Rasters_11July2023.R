#!/usr/bin/env Rscript

#######
# Interpolating GDD rasters between time periods for RCP scenarios
# GDD affects the potential location of crop production
# GDDs rasters are made by climate scenario (RCP, etc) by year (2010, 2015, etc)
#######

# loading libraries
library(raster)
library(rgdal)
library(plyr)
library(dplyr)
library(stringr)

# setting working directory
setwd('/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity')

# getting list of ssps
ssp.list <- list.files(paste0(getwd(),'/CMIP6_Climate_Data'), pattern = 'SSP', full.names = TRUE)

# getting years
years <- 
  list.files(ssp.list[1]) %>%
  .[!grepl('19',.)] %>%
  .[!grepl('Interpolated',.)]

# getting crop raster for gdd weighting
crop_2010 <- raster(paste0(getwd(),'/ESA LandCov Maps/Crop2010_Corrected_GlobCov.tif'))

# Looping through GDDs and SSPs and years
for(gdd in c('5C','10C')) {
  
  # looping through ssps
  for(ssp in ssp.list) {
    
    # Looping through years
    for(y in years) {
      cat(gdd,ssp,y,'\n')
 	     
      # need to load historic if y is 2021
      if(y %in% years[1]) {
        t0_raster <- 
          list.files(paste0(getwd(),'/CMIP6_Climate_Data/Historic/1995-2014/Managed_Rasters'),
                     pattern = 'GDD_', full.names = TRUE) %>%
          .[grepl('_Mollweide.tif',.)] %>% # Need rasters project in mollweide
          .[grepl(gdd,.)] %>% # need raster for appropriate gdd scenario
          raster() # loading as raster
        
        # time of t0
        mean_t0 <- 2005
      } else {
        t0_raster <- 
          list.files(paste0(ssp,'/',years[which(years %in% y)-1],'/Managed_Rasters'), full.names = TRUE) %>%
          .[grepl(gdd,.)] %>% # limiting to appropriate gdd scenario
          .[grepl('Mollweide.tif',.)] %>%
	  raster() # loading as raster
        
        # time of t0
	mean_t0 <- round(mean(as.numeric(unlist(str_extract_all(years[which(years %in% y)-1], '[0-9]{4,4}')))))
      } # End if statement to import t0_raster
      

      # importing t1_raster
      t1_raster <-
        list.files(paste0(ssp,'/',y,'/Managed_Rasters'), full.names = TRUE) %>%
        .[grepl(gdd,.)] %>% # appropriate gdd scenario
        .[grepl('Mollweide.tif',.)] %>%
        raster() # importing as raster
      cat('End Import t0 and t1 Rasters')

      # getting time in t1
      mean_t1 <- round(mean(as.numeric(unlist(str_extract_all(y, '[0-9]{4,4}')))))
      
      # getting years to interpolate
      years_interpolate <-
        seq(from = mean_t0,
            to = mean_t1,
            by = 5)
      
      # weightings for the different rasters
      weighting_t0 <-
        seq(from = 1,
            to = 0,
            length = length(years_interpolate))
      
      weighting_t1 <-
        seq(from = 0,
            to = 1,
            length = length(years_interpolate))
      
      cat('Looping through yy years\n')
      # looping through these years
      #for(yy in years_interpolate[!(years_interpolate %in% c(mean_t0, mean_t1))]) {
      for(yy in years_interpolate[!(years_interpolate %in% mean_t0)]) {
	      cat(gdd,ssp,y,yy,'\n')

        # weighting raster
        gdd_out <-
          t0_raster * weighting_t0[which(years_interpolate %in% yy)] +
          t1_raster * weighting_t1[which(years_interpolate %in% yy)]
        
        # getting gdd cutoffs if 2010
        # else applying gdd cutoffs to other years
        if(yy %in% 2010) {
          # overlap between gdd and crop
          # creating a data frame to do this
          gdd_crop_df <-
            data.frame(gdd = getValues(gdd_out),
                       crop = getValues(crop_2010)) %>%
            filter(!is.na(gdd)) %>% # removing cells with no gdd value
            filter(!is.na(crop)) %>% # removing cells with no crop value (e.g. these are oceans and lakes)
            filter(crop > 0) %>% # removing cells with no crop area
            filter(gdd > 0) %>%
            mutate(gdd_crop = gdd * crop) %>% # getting weighted average of these
            arrange(gdd) %>%
            mutate(sum_gdd_crop = cumsum(gdd_crop)) %>%
            mutate(prop_sum = sum_gdd_crop / max(sum_gdd_crop))
          
          # getting 2.5th and 97.5th percentile cutoffs
          gdd_low = 
            gdd_crop_df %>% 
            filter(prop_sum < .02505 & prop_sum > .02495)
          gdd_low = mean(gdd_low$gdd)
            
          gdd_high = 
            gdd_crop_df %>% 
            filter(prop_sum < .97505 & prop_sum > .97495)
          gdd_high = mean(gdd_high$gdd)
          
          # creating gdd binary raster
          gdd_binary <- gdd_out
          gdd_binary[gdd_binary < gdd_low] <- NA
          gdd_binary[gdd_binary > gdd_high] <- NA
          gdd_binary[!is.na(gdd_binary)] <- 1

          # we do not do this for pasture
          # given that e.g. goats and camels and many ruminants can live just about anywhere
        } else { # for other time periods
          # creating gdd binary raster
          gdd_binary <- gdd_out
          gdd_binary[gdd_binary < gdd_low] <- NA
          gdd_binary[gdd_binary > gdd_high] <- NA
          gdd_binary[!is.na(gdd_binary)] <- 1
        }

	# converting non-land areas to NAs
	gdd_binary[is.na(crop_2010)] <- NA
	gdd_out[is.na(crop_2010)] <- NA
        
        # saving raster
        writeRaster(gdd_out,
                    paste0(ssp,'/Interpolated_Rasters/GDD_',gdd,'_',yy,'.tif'),
                    overwrite = TRUE)
        
        
        # saving gdd binary raster
        writeRaster(gdd_binary,
                    paste0(ssp,'/Interpolated_Rasters/GDD_',gdd,'_binary_',yy,'.tif'),
                    overwrite = TRUE)
      } # End years loop
    } # End loop years
  } # End loop SSPs
} # End loop


