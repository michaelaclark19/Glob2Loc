#!/usr/bin/env Rscript

# Reprojecting minimum temperature maps
# Needed for SDMs

# Libraries
library(plyr)
library(dplyr)
library(stringr)
library(raster)
library(parallel)
# Number of cores to parallelise
n_cores = 10

# New cmip maps
# new.maps <- list.files('/Users/macuser/Desktop/Min_Temp_CMIP6', full.names = TRUE)
new.maps <- list.files('/data/pubh-glob2loc/pubh0329/Min_Temp_CMIP6', full.names = TRUE)
# Directories for old maps
old.dirs <- 
  list.dirs('/data/pubh-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/CMIP6_Climate_Data') %>%
  # list.dirs('/Users/macuser/Desktop/CMIP6_Climate_Data') %>%
  .[grepl('[0-9]{4,4}-[0-9]{4,4}$',.)]

# Setting working directory
setwd('/data/pubh-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity')

# Looping through dirs
for(dir in old.dirs) {
  # getting ssp of directory
  ssp = 
    gsub('.*CMIP6_Climate_Data/','',dir) %>%
    gsub('/[0-9]{4,4}-[0-9]{4,4}','',.)
  # Year of directory
  year = str_extract(dir, '[0-9]{4,4}-[0-9]{4,4}')
  # If ssp = historic
  if(ssp %in% 'Historic') {
    move.files <- new.maps[!grepl('SSP',new.maps)]
    move.files <- move.files[grepl(year,move.files)]
  } else {
    move.files <- new.maps[grepl(ssp,new.maps)]
    move.files <- move.files[grepl(year,move.files)]
  }
  
  # Making new directory
  dir.create(paste0(dir,'/Minimum_Temperature'))
  
  # And moving files
  file.copy(from = move.files,
            to = gsub('/data/pubh-glob2loc/pubh0329/Min_Temp_CMIP6',paste0(dir,'/Minimum_Temperature'),move.files))
}

# And reprojecting files

# Function to reproject
# Reproject function
reproject.fun.min_temp <-
  function(i) {
    #if(str_detect(i,'warmest_quarter')) {
    #  dest_file <- gsub('.tiff','_Mollweide.tiff',i)
    #} else {
    #  dest_file <- gsub(" -.*-",'',i) %>% gsub('\\([0-9]{1,2} models\\)','Mollweide',.)
    #}
    
    # Time period
    year = str_extract(i,'[0-9]{4,4}-[0-9]{4,4}')
    dest_file <- gsub('Minimum_Temperature/.*','',i)
    dest_file <- paste0(dest_file,'Managed_Rasters/Minimum_Temperature_Mollweide.tif')
    
    
    # dest_file <- gsub("CMIP.*Rasters/",'',i)
    # dest_file <- gsub('.tiff','.tif',dest_file)
    
    # And reprojecting Using gdal
    # gdalwarp(i,
    #          dest_file,
    #          # t_srs = 'EPSG:6350',
    #          t_srs = crs(template.raster),
    #          r = 'near',
    #          tr = res(template.raster),
    #          te = extent(template.raster),
    #          # compress = 'DEFLATE',
    #          pred = 2,
    #          zlevel = 3)
    
    # And deflating after rgdal - this results in very large files - need to do this by importing and then saving the file
    # tmp.raster.out <- raster(dest_file)
    
    # And alternative if gdal is not available...
    tmp.raster.out <-
      projectRaster(raster(i),
                    template.raster,
                    method = 'bilinear')
    
    # And writing raster
    writeRaster(tmp.raster.out, dest_file, overwrite = TRUE)
    
    # checking progress
    cat(i)
  }

# get wd
# setwd('/Users/macuser/Desktop')
# setwd('/data/pubh-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity')
# template raster
template.raster <- raster(paste0(getwd(),'/Global Mollweide Maps/IUCN.Region.Numerc Mollweide 1.5km.tif'))
# files to reproject
files.reproject <- list.dirs(paste0(getwd(),'/CMIP6_Climate_Data')) %>% .[grepl('Minimum_Temperature',.)]
files.reproject <- lapply(files.reproject, list.files, full.names = TRUE) %>% unlist()


# And reprojecting
mclapply(files.reproject,
         reproject.fun.min_temp, 
         mc.cores = n_cores)

# END
