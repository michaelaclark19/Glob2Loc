#!/usr/bin/env Rscript

### 
# Managing climate maps
# Reprojecting into global mollweide

library(plyr)
library(dplyr)
library(stringr)
library(raster)
library(parallel)
library(rgdal)
library(gdalUtils)
# Number of cores to parallelise
n_cores = 4
setwd('/data/pubh-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity')
# Setting working directory
setwd('/data/pubh-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity')

# List of directories
tmp <- list.dirs(paste0(getwd(),'/CMIP6_Climate_Data')) %>% .[grepl('Managed',.)]
tmp.files <- lapply(tmp,list.files, full.names = TRUE) %>% unlist() %>% .[grepl('.tif',.)] %>% .[!grepl('Mollweide',.)]

# Template raster
template.raster <- raster(paste0(getwd(),'/Global Mollweide Maps/IUCN.Region.Numerc Mollweide 1.5km.tif'))

### For testing
# template.raster <- raster::aggregate(template.raster,15)
# template.raster.aggregate <- raster::aggregate(template.raster,3)
# template.raster <- template.raster.aggregate

# Reproject function
reproject.fun <-
  function(i) {
    #if(str_detect(i,'warmest_quarter')) {
    #  dest_file <- gsub('.tiff','_Mollweide.tiff',i)
    #} else {
    #  dest_file <- gsub(" -.*-",'',i) %>% gsub('\\([0-9]{1,2} models\\)','Mollweide',.)
    #}
	  dest_file <- gsub('.tif.*','_Mollweide.tif',i)
    
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

split.files <- split(tmp.files,             # Applying split() function
      ceiling(seq_along(tmp.files) / n_cores))
# Doing this in a loop, because doing this in parallel makes the cluster unhappy...
for(l in rev(split.files[[1]])) {
  cat('\n',l,'\n')
  reproject.fun(l)
}

#split.files <- split(tmp.files,             # Applying split() function
#      ceiling(seq_along(tmp.files) / n_cores))

# Doing in a loop to check progress
#for(l in 1:length(split.files)) {
#	cat(split.files[[l]])
#	mclapply(split.files[[l]], reproject.fun, mc.cores = n_cores)
#}

# END
