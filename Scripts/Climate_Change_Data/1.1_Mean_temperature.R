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

# Number of cores to parallelise
n_cores = 20

# List of directories - only need mean temp and total precip
tmp <- list.dirs(paste0(getwd(),
		 '/CMIP6_Climate_Data')) %>% 
  .[grepl('[0-9]{4,4}$',.)] %>%
  .[grepl('Historic.*1995|SSP.*2021|SSP.*2041',.)]

# Making function to loop through
temp.fun <-
	function(i) {
	# Order of months
	month.order <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

	# List of rasters to import
	raster.files <- list.files(paste0(i,'/Mean_temp'),full.names=TRUE) %>% .[!grepl('Mollweide',.)] %>% .[!grepl('.aux',.)] #%>% .[!grepl('.tiff',.)]

	# Making empty stack to append the raster for each month
	raster.stack <- raster::stack()
    for(m in month.order) {raster.stack <- raster::stack(raster.stack,raster(raster.files[grepl(m,raster.files)]))}
    names(raster.stack) <- month.order
    # Creating data frame
    raster.df <- data.frame(matrix(ncol = 12, nrow = raster.stack[[1]]@nrows*raster.stack[[1]]@ncols))
    # names(raster.df) <- names(raster.stack)
    for(k in 1:12) {raster.df[,k] <- getValues(raster.stack[[k]])}
    names(raster.df) <- names(raster.stack)

    # Need to multiply by days in each month
    month.length <- c(31,28.25,31,30,31,30,31,31,30,31,30,31)

    # Multiplying
    raster.df[,month.order] <- raster.df[,month.order] * month.length

    # Row summing
    temp_sum <- rowSums(raster.df[,month.order])
    # Averaging by days in the year
    temp_avg <- temp_sum / sum(month.length)

    # Making raster
    temp_avg_raster <- raster(matrix(temp_avg, nrow = raster.stack@nrows, byrow = TRUE), crs = crs(raster.stack))
    extent(temp_avg_raster) = extent(raster.stack)

    # And saving raster
    writeRaster(temp_avg_raster,
		paste0(i,'/Managed_Rasters/Mean_annual_temp.tif'),
		overwrite = TRUE)
	}

# And running
mclapply(tmp,temp.fun,mc.cores = n_cores)
