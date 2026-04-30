#!/usr/bin/env Rscript

###
# Projecting Future SDM maps
# Taking clipped 2010 map
# Buffering around 2010 map based on migration and dispersal distances
# Clipping future maps based on buffer around 2010 map and potentially suitable habitat from SDMs
###

# Getting python environment
# Sys.setenv(RETICULATE_PYTHON = "/apps/system/easybuild/software/Anaconda3/2022.05/bin/python")

# Setting working directory
setwd("/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity")

# Loading functions
source(paste0(getwd(),'/Scripts/SDM Scripts/0.0_SDM_Functions_2024.07.13.R'))

# Setting chunk_n - variable to indicate which chunk of species to loop over
chunk_n <- 1

# Loading libraries
library(reticulate)
library(raster)
library(plyr)
library(dplyr)
library(stringr)
library(parallel)
library(future.apply)

# Loading python libraries
py_run_string("import numpy as np") # Load numpy  package
py_run_string("import pandas as pd") # Load pandas package
py_run_string("from skimage import io, morphology, measure")
py_run_string("from scipy.ndimage import convolve")
py_run_string("from PIL import Image")
py_run_string("import copy")

# Loading realm+biome raster and extents
biome_realm <- raster(paste0(getwd(),'/Global Mollweide Maps/TNC Ecoregions Mollweide 1.5km.tif')) # Importing raster for ecoregions - used to clip 2005 maps

# Creating buffer for dispersal capacity
# 5 km buffer for amphibians
# circle_amps = matrix(1,nrow = 12, ncol = 12)
# py_run_string("xx, yy = np.mgrid[:r.circle_amps.shape[0],:r.circle_amps.shape[1]]")
# py_run_string("circle_amps = (xx - 6) ** 2 + (yy - 6) ** 2")
# py_run_string("circle_amps[circle_amps<=16]=1")
# py_run_string("circle_amps[circle_amps>16]=0")

circle_amps = matrix(1,nrow = 9, ncol = 9)
py_run_string("xx, yy = np.mgrid[:r.circle_amps.shape[0],:r.circle_amps.shape[1]]")
py_run_string("circle_amps = (xx - 4) ** 2 + (yy - 4) ** 2")
py_run_string("circle_amps[circle_amps<=9]=1")
py_run_string("circle_amps[circle_amps>9]=0")

# 15 km buffer for birds
# circle_mards = matrix(1,nrow = 33, ncol = 33)
# py_run_string("xx, yy = np.mgrid[:r.circle_mards.shape[0],:r.circle_mards.shape[1]]")
# py_run_string("circle_mards = (xx - 16) ** 2 + (yy - 16) ** 2")
# py_run_string("circle_mards[circle_mards<=225]=1")
# py_run_string("circle_mards[circle_mards>225]=0")

circle_mards = matrix(1,nrow = 21, ncol = 21)
py_run_string("xx, yy = np.mgrid[:r.circle_mards.shape[0],:r.circle_mards.shape[1]]")
py_run_string("circle_mards = (xx - 10) ** 2 + (yy - 10) ** 2")
py_run_string("circle_mards[circle_mards<=100]=1")
py_run_string("circle_mards[circle_mards>100]=0")

# List of all species
species.list <-
  do.call(c,lapply(paste0(getwd(),'/ESH_RCPs/SSP1-2.6/',c('Amphibians','Mammals','Birds','Reptiles')),
                   list.files,
                   full.names = TRUE)) %>%
  .[grepl('csv',.)] %>%
  .[grepl('completed',.)] %>%
  .[!grepl('_not_',.)] %>%
  .[!grepl('_not\\b',.)] %>%
  gsub('.*SSP1-2.6/','',.) %>%
  gsub('_completed.*','',.) %>%
  gsub('_climate_suitable','',.)

# List of species we have
species.have <-
  do.call(c,lapply(paste0(getwd(),'/ESH_RCPs/SSP1-2.6/',c('Amphibians','Mammals','Birds','Reptiles')),
                   list.files,
                   full.names = TRUE)) %>%
  .[grepl('csv',.)] %>%
  .[grepl('migrate',.)] %>%
  gsub('.*SSP1-2.6/','',.) %>%
  gsub('_migrate.*','',.)

# If we have species, then don't need to do anything with them
if(length(species.have) > 0) {
  species.need <- species.list[!(species.list %in% species.have)]
} else {
  species.need <- species.list
}

# dropping specices that have not been completed
species.need <- species.need[!grepl('_not\\b|_error\\b',species.need)]


# Getting time periods to loop
years.loop <-
  list.files(paste0(getwd(),'/ESH_RCPs/SSP1-2.6/Amphibians')) %>%
  str_extract('[0-9]{4,4}') %>% unique() %>% as.numeric() %>% .[!is.na(.)] %>%
  .[. <= 2050]

# List of SSPs
ssps <-
  list.files('/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/ESH_RCPs', pattern = 'SSP1')

# List of threshold types
thresholds <- c('prevalence','specsens')

# Splitting into sublists
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

# Chunking
# n_cores = 40
# species.need.list <- chunk2(species.need, n_cores)

cat(ssps,'\n')
cat(years.loop,'\n')

cat('Amphibians:', sum(grepl('Amphibian',species.have)), 'Species Completed \n')
cat('Birds:', sum(grepl('Bird',species.have)), 'Species Completed \n')
cat('Mammals:', sum(grepl('Mammal',species.have)), 'Species Completed \n')
cat('Reptiles:', sum(grepl('Reptile',species.have)), 'Species Completed \n')

cat('Amphibians:', sum(grepl('Amphibian',species.need)), 'Species Needed \n')
cat('Birds:', sum(grepl('Bird',species.need)), 'Species Needed \n')
cat('Mammals:', sum(grepl('Mammal',species.need)), 'Species Needed \n')
cat('Reptiles:', sum(grepl('Reptile',species.need)), 'Species Needed \n')

# Randomising order
species.need <- sample(species.need, length(species.need), replace = FALSE)

for(i in species.need) { # For testing; chunk_n is a variable that indicates which chunk to loop over - doing this so it is easier to change the script as needed
  # Checking if we have the species again...doing this for parallelisation, when order species loop through is randomised
  species.have <-
    do.call(c,lapply(paste0(getwd(),'/ESH_RCPs/SSP1-2.6/',c('Amphibians','Mammals','Birds','Reptiles')),
                     list.files,
                     full.names = TRUE)) %>%
    .[grepl('csv',.)] %>%
    .[grepl('migrate',.)] %>%
    gsub('.*SSP1-2.6/','',.) %>%
    gsub('_migrate.*','',.)
  
  if(i %in% species.have) {
     # Do Nothing
  } else { # Do everything
    
    # species_wrap_fun <- function(i) { # For running
    # Getting taxon
    taxa <- str_extract(i, 'Amphibians|Birds|Mammals|Reptiles')
    
    # Getting name
    tmp.name <-
      gsub(paste0('.*',taxa),'',i) %>%
      gsub('\\/','',.) %>%
      gsub('2010.*','',.)
    
    # Getting buffer
    if(taxa %in% 'Amphibians' | taxa %in% 'Reptiles') {
      species_buffer <- py$circle_amps
    } else {
      species_buffer <- py$circle_mards
    }
    
    # If else statement - some species throw an error, so only running on species that did not throw an error!
    if(grepl('_error\\b',tmp.name)) {
      cat('Skipping: ', taxa, tmp.name, '\n')
      # do nothing
      
    } else { # do everything!
      
      # Getting esh
      clipped.esh <-
        raster(list.files(paste0(getwd(),'/ESH_Tifs_12Oct/',taxa,'/'), pattern = tmp.name, full.names = TRUE))
      
      # Getting overlap with ecoregion map
      species.ecoregion <- raster::crop(biome_realm, clipped.esh)
      
      # Which ecoregions is a species in? Used to clip 2005 map
      species.ecoregion[!(clipped.esh %in% 1)] <- NA
      eco.keep <- 
        getValues(species.ecoregion) %>% 
        unique() %>% 
        .[!is.na(.)] %>% #
        .[.>1] # this last is a double check to remove oceans
      
      
      
      # Looping through threshold types
      for(t in thresholds) {
        # Looping through SSPs
        for(s in ssps) {
          
          # Looping through years
          for(y in years.loop) {
            # Tracking progress
            cat(taxa,tmp.name,t,s,y,'\n')
            
            # Importing sdm output
            sdm.year <-
              list.files(paste0(getwd(),'/ESH_RCPs/',s,'/',taxa), pattern = tmp.name, full.names = TRUE) %>%
              .[grepl(y,.)] %>%
              .[grepl(t,.)] %>%
              .[grepl('[0-9]{4,4}.tif',.)] %>%
              raster(.)
            
            # Converting to 0s and 1s for buffering
            sdm.year[is.na(sdm.year)] <- 0
            
            if(y %in% '2005') { # Only need to do this once
              # Clipping species sdm to be within current ecoregions
              
              # Getting ecoregion map to overlap with the SDMs
              eco.sdm <- raster::crop(biome_realm, sdm.year)
              # Removing cells outside of the ecoregion in which a species currrently exists
              sdm.year[!(eco.sdm %in% eco.keep)] <- 0
              
              # Saving raster in 2005
              writeRaster(sdm.year,
                          paste0(getwd(),'/ESH_RCPs/',s,'/',taxa,'/',tmp.name,'_',s,'_',y,'_',t,'_migrate.tif'),
                          overwrite = TRUE)
              
              # Converting to matrix to buffer
              # This is used as input to buffer future maps
              sdm_matrix <- matrix(sdm.year, nrow = sdm.year@nrows, byrow = TRUE)
              
            }  else { # End statement checking for 2005
              
              # Buffering SDM matrix from previous year
              py_run_string("species_buffer = copy.copy(r.species_buffer)")
              py_run_string("input_array = copy.copy(r.sdm_matrix)")
              py_run_string("input_array[np.isnan(input_array)] = 0")
              py_run_string("esh_buffer = convolve(input_array, species_buffer, mode='constant')")
              
              # Converting back to 1s and 0s (presence and absence)
              py_run_string("esh_buffer[esh_buffer>0]=1")
              
              # Converting buffered SDM back to a raster
              sdm.buffer <- raster(py$esh_buffer, crs = crs(sdm.year))
              extent(sdm.buffer) <- extent(eco.sdm)
              
              # Clipping SDM based on the buffer (migration distance)
              sdm.year[!(sdm.buffer %in% 1)] <- 0
              
              # And clipping to remove oceans
              tmp_clip_biome_realm <- raster::crop(biome_realm,sdm.year)
              tmp_clip_biome_realm[!(tmp_clip_biome_realm %in% eco.keep)] <- NA
              sdm.year[is.na(tmp_clip_biome_realm)] <- NA
              
              
              cat('Writing Raster \n')
              # writing raster
              writeRaster(sdm.year,
                          paste0(getwd(),'/ESH_RCPs/',s,'/',taxa,'/',tmp.name,'_',s,'_',y,'_',t,'_migrate.tif'),
                          overwrite = TRUE)
              
              # And updating matrix
              # This is used as the input to the next time period
              sdm_matrix <- matrix(sdm.year, nrow = sdm.year@nrows, byrow = TRUE)
              
              gc()
              
            } # End statement for other years
          } # End loop through years
        } # End loop through SSPs
      } # End loop through thresholds
      
      cat('Writing CSV \n')
      # Saving csv to keep tabs on species
      write.csv(data.frame(species = tmp.name, taxa = taxa, completed = 'yes'),
                paste0(getwd(),'/ESH_RCPs/SSP1-2.6/',taxa,'/',tmp.name,'_migrate.csv'))
      
      
      
    } # End loop through species
  } # End if statement seeing if we need to loop through species
} 
 

    
