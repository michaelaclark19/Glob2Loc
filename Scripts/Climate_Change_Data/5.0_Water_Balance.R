#!/usr/bin/env Rscript

# Getting water balance
# Needed for SDMs

# Libraries
library(plyr)
library(dplyr)
library(stringr)
library(raster)
library(parallel)

# Number of cores to parallelise
n_cores = 5

# Setting work directory
setwd('/data/pubh-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity')
# setwd('/Users/macuser/Desktop/')

# List of directories
dir.list <- 
  list.dirs(paste0(getwd(),'/CMIP6_Climate_Data')) %>%
  .[grepl('[0-9]{4,4}-[0-9]{4,4}$',.)] %>%
  .[grepl('Historic.*1995|SSP.*2021|SSP.*2041',.)]

# Checking if any already exist
already.have <- unlist(lapply(dir.list, function(i){return(list.files(paste0(i,'/Managed_Rasters'),full.names = TRUE, pattern = 'Wat.*Bal.*Mollweide.tif'))}))
already.have <- gsub("/Managed.*","",already.have) %>% unique()
# Which do we need?
if(length(already.have)>0) {dir.list <- dir.list[!(dir.list %in% already.have)]}

# Template raster to repoject maps
template.raster <- raster(paste0(getwd(),'/Global Mollweide Maps/IUCN.Region.Numerc Mollweide 1.5km.tif'))

# Making function
wat.balance.function <-
  function(dir) {
	  cat(dir)
    # List of months
    month.list <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
    # Getting data frame of precipitation
    precip.files <- 
      list.files(paste0(dir,'/Total_precip'),full.names=TRUE) %>%
      .[!grepl('.aux',.)] %>%
      .[!grepl('Mollweide',.)]
    
    # Logic check - don't run if no files
    if(length(precip.files) %in% 0 |
       grepl('Historic',dir) & grepl('20[0-9]{2,2}-2[0-9]{3,3}',dir) |
       grepl('SSP',dir) & grepl('18[0-9]{2,2}|19[0-9]{2,2}',dir)) {
      
      # Do nothing
      return(NULL)
      
    } else {
      # Template raster for number of rows
      tmp.raster <- raster(precip.files[1])
      
      # Making data frame
      precip.df <- data.frame(matrix(ncol = 12, nrow = (tmp.raster@nrows*tmp.raster@ncols)))
      names(precip.df) <- month.list
      
      # Adding to the data frame
      for(m in month.list) {
        precip.df[,m] <- getValues(raster(precip.files[grepl(m,precip.files)]))
      }

	# Multiplying by days in each month
	precip.df <- precip.df * c(31,28.25,31,30,31,30,31,31,30,31,30,31)

      # Getting data frame of temperature
      # Getting data frame of precipitation
      temp.files <- 
        list.files(paste0(dir,'/Mean_temp'),full.names=TRUE) %>% 
        .[!grepl('.aux',.)] %>%
        .[!grepl('Mollweide',.)]
      
      # Making data frame
      temp.df <- data.frame(matrix(ncol = 12, nrow = (tmp.raster@nrows*tmp.raster@ncols)))
      names(temp.df) <- month.list
      
      # Adding to the data frame
      for(m in month.list) {
        temp.df[,m] <- getValues(raster(temp.files[grepl(m,temp.files)]))
      }
      
      # Now calculating water balance
      # This is from Skov and Skenning, 2004, Ecography
      wat_balance <- precip.df - temp.df * 58.93 / 12 # This is precipitation - potential evapotranspiration, where PET = temp * 58.93/12
      
      # Annual water balance
      tot_wat_balance <- 
        rowSums(wat_balance, na.rm = TRUE)
      
      # And making raster
      wat_balance_raster = raster(matrix(tot_wat_balance, 
                                         nrow = tmp.raster@nrows,
                                         ncol = tmp.raster@ncols,
                                         byrow = TRUE),
                                  crs = crs(tmp.raster))
      extent(wat_balance_raster) = extent(tmp.raster)
      
      # For checking
      # return(wat_balance_raster)
      
      # # And saving
      writeRaster(wat_balance_raster,
                  paste0(dir,'/Managed_Rasters/Water_Balance.tif'),
                  overwrite = TRUE)
      
      # And now reprojecting
      tmp.raster.out <-
        projectRaster(wat_balance_raster,
                      template.raster,
                      method = 'bilinear')
      # And saving
      writeRaster(tmp.raster.out,
                  paste0(dir,'/Managed_Rasters/Water_Balance_Mollweide.tif'),
                  overwrite = TRUE)
      # Returning list to check if it works
      return(list(wat_balance_raster,tmp.raster.out))
    }
  }
 
# For checking
# for(dir in dir.list[1:9]) {
#   wat.balance.function(dir)
# }

# And running in parallel
mclapply(dir.list,
         wat.balance.function,
         mc.cores=n_cores)
# END




