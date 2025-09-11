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
  .[grepl('[0-9]{4,4}-[0-9]{4,4}$',.)]
# Which do we have?
dirs.have <-
	lapply(paste0(dir.list,'/Managed_Rasters'), list.files, full.names = TRUE) %>%
	.[(grepl('Annual_P',.))] %>%
	unique() %>%
	unlist() %>%
	gsub('/Managed_Rasters.*','',.) %>%
	unique()
# Ones we need
dir.list <- dir.list[!(dir.list %in% dirs.have)]
# Which have managed rasters already
dir.list.managed <- 
	unlist(lapply(dir.list,list.files,full.names=TRUE)) %>%
	.[grepl('Managed_',.)] %>%
	gsub('/Managed.*','',.) %>%
	unique()
# And filtering again
dir.list <- dir.list[dir.list %in% dir.list.managed]


# Template raster to repoject maps
template.raster <- raster(paste0(getwd(),'/Global Mollweide Maps/IUCN.Region.Numerc Mollweide 1.5km.tif'))

# Making function
precip.function <-
  function(dir) {
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

      # Multiplying by days in the month
      precip.df <-
	      precip.df *
	      c(31,28.25,31,30,31,30,31,31,30,31,30,31)
      # Getting total precip
      tot_precip <-
        rowSums(precip.df, na.rm = TRUE)

      # And making raster
      precip_raster = raster(matrix(tot_precip, 
                                         nrow = tmp.raster@nrows,
                                         ncol = tmp.raster@ncols,
                                         byrow = TRUE),
                                  crs = crs(tmp.raster))
      extent(precip_raster) = extent(tmp.raster)
      
      # For checking
      # return(wat_balance_raster)
      
      # # And saving
      writeRaster(precip_raster,
                  paste0(dir,'/Managed_Rasters/Annual_Precipitation.tif'),
                  overwrite = TRUE)
      
      # And now reprojecting
      tmp.raster.out <-
        projectRaster(precip_raster,
                      template.raster,
                      method = 'bilinear')
      # And saving
      writeRaster(tmp.raster.out,
                  paste0(dir,'/Managed_Rasters/Annual_Precipitation_Mollweide.tif'),
                  overwrite = TRUE)
      # Returning list to check if it works
      return(list(precip_raster,tmp.raster.out))
    }
  }
 
# For checking
#   wat.balance.function(dir)
# }

# And running in parallel
mclapply(dir.list,
         precip.function,
         mc.cores=n_cores)
# END




