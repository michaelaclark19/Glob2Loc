#!/usr/bin/env Rscript

### 
# Managing global climate maps
# Getting GDD @ 5C and 10C
# Getting variation in precipitation
# Getting variation in temperature

library(plyr)
library(dplyr)
library(stringr)
library(raster)
library(parallel)
library(envirem)
library(matrixStats)

# Setting working directory
setwd('/data/pubh-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity')
# setwd('/Users/macuser/Desktop/')

# Number of cores to parallelise
n_cores = 10

# List of directories - only need mean temp and total precip
tmp <- list.dirs(paste0(getwd(),"/CMIP6_Climate_Data")) %>% .[grepl('[0-9]{4,4}$',.)]
tmp <- tmp %>% .[!(grepl('Historic.*20[0-9]{2,2}-2',.))] %>% .[!grepl('SSP.*19',.)]

# Function to stack rasters
raster.stack.fun <-
  function(raster.files) {
    # Order to import files
    month.order <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
    # Creating stack
    raster.stack <- raster::stack()
    for(m in month.order) {raster.stack <- raster::stack(raster.stack,raster(raster.files[grepl(m,raster.files)]))}
    names(raster.stack) <- month.order
    # Creating data frame
    raster.df <- data.frame(matrix(ncol = 12, nrow = raster.stack[[1]]@nrows*raster.stack[[1]]@ncols))
    # names(raster.df) <- names(raster.stack)
    for(k in 1:12) {raster.df[,k] <- getValues(raster.stack[[k]])}
    names(raster.df) <- names(raster.stack)
    
    # item to return
    raster.stack <- raster::stack(raster.stack)
    return.list <- list(raster.stack,raster.df)
    names(return.list) <- c('Rasters','Data.Frame')
    # And returning the list
    return(return.list)
  }

# Main function to get GDD, var temp, and var precip
mclapply.function <-
  function(i) {
    # List of precip and temp files
    precip_files <-list.files(paste0(i,'/Total_precip'), full.names = TRUE)
    temp_files <- list.files(paste0(i,'/Mean_temp'), full.names = TRUE)
    
    
    
    # Dropping unwanted rasters
    if(sum(grepl('Mollweide',precip_files))>0) {precip_files <- precip_files[!grepl('Mollweide',precip_files)]}
    if(sum(grepl('Mollweide',temp_files))>0) {temp_files <- temp_files[!grepl('Mollweide',temp_files)]}
    
    if(sum(grepl('.aux',precip_files))>0) {precip_files <- precip_files[!grepl('.aux',precip_files)]}
    if(sum(grepl('.aux',temp_files))>0) {temp_files <- temp_files[!grepl('.aux',temp_files)]}
    
    
    
    # Getting stacks of the rasters
    precip.info <- raster.stack.fun(precip_files)
    # precip.info <- raster.stack.fun(list.files(paste0(i,'/Total_precip'), full.names = TRUE))
    temp.info <- raster.stack.fun(temp_files)
    # temp.info <- raster.stack.fun(list.files(paste0(i,'/Mean_temp'), full.names = TRUE))
    
    # Getting GDD
    gdd_5c = growingDegDays(temp.info[['Rasters']],5,tempScale=1) # tempScale indicates in degrees C
    gdd_10c = growingDegDays(temp.info[['Rasters']],10,tempScale=1) # tempScale indicates in degrees C
    
    # And saving
    writeRaster(gdd_5c,paste0(i,'/Managed_Rasters/GDD_5C.tif'),overwrite = TRUE)
    writeRaster(gdd_10c,paste0(i,'/Managed_Rasters/GDD_10C.tif'),overwrite = TRUE)
    
    # Getting var temp
    temp_var <- rowSds(as.matrix(temp.info[['Data.Frame']]))
    temp_var <- temp_var * 100
    
    # Getting var precip
    precip_sd <- rowSds(as.matrix(precip.info[['Data.Frame']]))
    precip_mean <- rowMeans(as.matrix(precip.info[['Data.Frame']]))
    precip_var <- (precip_sd / precip_mean) * 100
    
    
    # Converting temp_var and precip_var back into rasters
    temp_var_raster <- raster(matrix(temp_var,ncol = precip.info[['Rasters']]@ncols,byrow=TRUE),crs = crs(precip.info[['Rasters']]))
    extent(temp_var_raster) <- extent(precip.info[['Rasters']])
    
    precip_var_raster <- raster(matrix(precip_var,ncol = precip.info[['Rasters']]@ncols,byrow=TRUE),crs = crs(precip.info[['Rasters']]))
    extent(precip_var_raster) <- extent(precip.info[['Rasters']])
    
    # Saving rasters
    writeRaster(temp_var_raster,paste0(i,'/Managed_Rasters/Variance_of_temperature.tif'),overwrite = TRUE)
    writeRaster(precip_var_raster,paste0(i,'/Managed_Rasters/Variance_of_precipitation.tif'),overwrite = TRUE)
    
    # And getting mean temperature
    # Need to multiply by days in each month
    month.order <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
    month.length <- c(31,28.25,31,30,31,30,31,31,30,31,30,31)
    
    # Multiplying
    raster.df <- temp.info[['Data.Frame']]
    for(k in 1:12) {raster.df[,k] <- raster.df[,k] * month.length[k]}
    
    # Row summing
    temp_sum <- rowSums(raster.df[,month.order])
    # Averaging by days in the year
    temp_avg <- temp_sum / sum(month.length)
    
    # Making raster
    temp_avg_raster <- raster(matrix(temp_avg, nrow = temp.info[['Rasters']]@nrows, byrow = TRUE), crs = crs(temp.info[['Rasters']]))
    extent(temp_avg_raster) = extent(temp.info[['Rasters']])
    
    # And saving raster
    writeRaster(temp_avg_raster,
                paste0(i,'/Managed_Rasters/Mean_annual_temp.tif'),
                overwrite = TRUE)
  }

# for(i in tmp[1:3]) {mclapply.function(i) 
#   cat(i)}

mclapply(tmp, mclapply.function, mc.cores = n_cores)
