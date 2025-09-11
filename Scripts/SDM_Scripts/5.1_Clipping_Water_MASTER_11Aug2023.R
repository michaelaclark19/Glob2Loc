#!/usr/bin/env Rscript

###
# Script to take 5.0 outputs
# which contain maps of potentially suitable habitat
# that include ocean and water
# to remove ocean and water from these maps


# libraries
library(raster)
library(parallel)
library(plyr)
library(dplyr)

# working directory
setwd('/data/pubh-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/')

# creating folders
dirs <-
  paste0(getwd(),'/ESH_RCPs/SSP2-4.5/') %>%
  paste0(.,'Climate_Suitable_') %>%
  paste0(.,c('Amphibians','Birds','Mammals','Reptiles'))


lapply(dirs, dir.create)

# template raster
template <- raster(paste0(getwd(),'/Ecoregions_Feb2023/biomes_realms_raster.tif'))
bio <- raster(paste0(getwd(),'/Ecoregions_Feb2023/biomes_raster.tif'))

# converting oceans (1) and africa (48)to NAs
template[template %in% 1] <- NA
template[is.na(bio)] <- NA


# list of files from 5.0 script that contain water
file.list <-
  do.call(c,lapply(paste0(getwd(),'/ESH_RCPs/SSP2-4.5/') %>%
           paste0(.,c('Amphibians','Birds','Mammals','Reptiles')),
         list.files,
         full.names = TRUE)) %>%
  .[!grepl('migrate',.)] %>%
  gsub('.*SSP2-4.5/','',.)

# list of files already converted
files.have <-
  do.call(c,lapply(paste0(getwd(),'/ESH_RCPs/SSP2-4.5/') %>%
           paste0(.,'Climate_Suitable_') %>%
           paste0(.,c('Amphibians','Birds','Mammals','Reptiles')),
         list.files,
         full.names = TRUE)) %>%
  gsub('.*Climate_Suitable_','',.)

# list of files that need to be converted
files.need <- file.list[!(file.list %in% files.have)]




# function
remove_water_function <-
  function(i) {
    cat(i,'\n')
    # importing raster
    tmp <- raster(paste0(getwd(),'/ESH_RCPs/SSP2-4.5/',i))

    # cropping template
    tmp_template <- raster::crop(template,tmp)

    # removing water
    tmp[is.na(tmp_template)] <- NA

    # saving raster
    writeRaster(tmp,
                paste0(getwd(),'/ESH_RCPs/SSP2-4.5/Climate_Suitable_',i))

    # removing rasters
    rm(tmp,tmp_template)
  }

# testing
# lapply(files.need[1:100], remove_water_function)

# chunking to run in parallel
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
files.need <- chunk2(files.need,5)

# running in parallel
#for(i in files.need[[1]]) {
#	cat(i,'\n')
#	remove_water_function(i)
#}



mclapply(files.need[[1]], remove_water_function, mc.cores = 20)
