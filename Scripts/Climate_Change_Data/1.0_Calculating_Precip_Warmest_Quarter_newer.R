### 
# Managing climate data
# Getting precip in warmest quarter
# variation in precip
# variation in temperature
# and GDD

# Getting precip in warmest quarter before reprojecting
# Getting others after reprojecting

library(plyr)
library(dplyr)
library(stringr)
library(raster)
library(parallel)

# number of cores
n_cores = 10

# working directory
setwd('/data/pubh-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity')

# Getting a template raster for reprojecting, etc
raster.files <- list.files(paste0(getwd(),'/CMIP6_Climate_Data'))
# Template raster
template.raster <- raster(raster.files[1])

# Precip in warmest quarter
raster.stack.fun <-
  function(raster.files) {
    # Order to import files
    month.order <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
	month.length <- c(31,28.25,31,30,31,30,31,31,30,31,30,31)
    # Creating stack
    raster.stack <- stack()
    raster.stack <- lapply(month.order, function(i) {raster.stack <- stack(raster.stack,raster(raster.files[grepl(i,raster.files)]))})
    names(raster.stack) <- month.order
    # Creating data frame
    raster.df <- data.frame(matrix(ncol = 12, nrow = raster.stack[[1]]@nrows*raster.stack[[1]]@ncols))
    names(raster.df) <- names(raster.stack)
    for(i in 1:12) {raster.df[,i] <- getValues(raster.stack[[i]])}
    
    # And returning the df
    return(raster.df)
  }

tmp.df <- raster.stack.fun(raster.files)

# Temperature in warmest quarter
temp.warmest.quarter <-
  function(df) {
    # Making a new data frame
    temp.df <- data.frame(matrix(nrow = nrow(df),ncol = ncol(df)))
    names(temp.df) <- names(df)
    # Rbinding Jan and Feb to end of data frame - makes this easier
    df <- cbind(df,df[,1:2])
    # Adding precip over 3 month time period
    for(i in 1:12) {temp.df[,i] <- rowSums(df[,seq(from = i, to = i+2, by = 1)])}
    # Which of these has max precipitation
    df$max_temp <- apply(df,MARGIN = 1, FUN = which.max)
    # Returning
    return(df)
  }

# Now getting precip in warmest quarter
precip.warmest.fun <- 
  function(df) {
    # # Making a new data frame
    # precip.df <- data.frame(matrix(nrow = nrow(df),ncol = ncol(df)))
    # names(temp.df) <- names(df)
    # Rbinding Jan and Feb to end of data frame - makes this easier
    df <- cbind(df,df[,1:2])
    # Getting warmest quarter
    df$warmest_quarter <- temp.df.warmest$max_temp
    # Getting new data frame
    df$warmest_precip <- NA
    # looping to sum
    for(i in 1:nrow(df)) {df$warmest_precip[i] <- rowSums(df[i,seq(from=df$warmest_quarter[i],to=df$warmest_quarter[i]+2,by=1)])}
    
    # And returning data frame
    return(df)
    
  } # End function

# Now creating wrapper to do this across time periods
mclapply.function <-
  function(i) {
    
    # Logical check to see if we nee dto make the raster
    if(length(list.files(paste0(i,'/Total_precip'),full.names=TRUE)) != 12) {
      # Do nothing
    } else { # Do everything!
      # Stacking and extracting precip and temp data
      # Precip
      precip.df <- raster.stack.fun(list.files(paste0(i,'/Total_precip'),full.names=TRUE))
    precip.df <- precip.df * c(31,28.25,31,30,31,30,31,31,30,31,30,31)
      # Temperature
      temp.df <- raster.stack.fun(list.files(paste0(i,'/Mean_temp'),full.names=TRUE))
      
      # Warmest quarter
      temp.df.warmest <- temp.warmest.quarter(temp.df)
      # Precipitation in warmest quarter
      precip.df.warmest <- precip.warmest.fun(precip.df)
      
      # And creating raster
      save.raster <- raster(matrix(precip.df.warmest$warmest_precip,nrow = template.raster@nrows,byrow=TRUE), crs = template.raster@crs)
      extent(save.raster) <- extent(template.raster)
      # Creating directory
      dir.create(paste0(i,'/Managed_Rasters'))
      # Then saving
      writeRaster(save.raster,
                  paste0(i,'/Managed_Rasters/Precip_in_warmest_quarter.tiff'),
                  overwrite = TRUE)
    }
  }

# List of directories to loop over
tmp <- list.dirs(paste0(getwd(),"/CMIP6_Climate_Data")) %>% .[grepl('[0-9]{4,4}$',.)]

# And running
mclapply(tmp,mclapply.function,mc.cores=n_cores)

# END
   

