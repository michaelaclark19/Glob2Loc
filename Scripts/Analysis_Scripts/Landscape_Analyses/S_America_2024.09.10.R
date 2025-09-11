#!/usr/bin/env Rscript

#####
# Getting landscape level outcomes, for species in the Amazon, for different stressors
# Doing this in three steps:
# (1) Outcomes for each stressor scenario
# (2) By taxa
# (3) Which will then be cropped to the amazonian basin

#####
# Libraries
library(raster)
library(plyr)
library(dplyr)
library(parallel)

#####
### General data management
# Setting wd
setwd('/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/')

# How many cores to use
mc_cores = 20

# Getting list of scenarios
scen_list_full <-
  c('current.extent.esh.tif',
    'future.extent.esh.exp.int.urb.tif',
    'future.extent.esh.exp.tif',
    'future.extent.esh.int.tif',
    'future.extent.esh.urb.tif')

scen_list_2020 <- 'future.extent.esh.exp.int.urb.tif'

# Getting list of taxa
taxa_list <- c('Amphibians','Birds','Mammals','Reptiles')

# Getting theshold list
thresh_list <- c('prevalence','specsens')

# Getting outcome types
outcome_list <- c('hab_availability','pop_density')

# Getting list of output files
# full_file_list <-
#   do.call(c,
#           lapply(paste0(getwd(),'/Outputs/Raster_Outputs/BAU/',taxa_list),
#                  list.files,
#                  full.names=TRUE))

# making template raster
# s_am_raster <- raster('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Ecoregions_Feb2023/biomes_realms_raster.tif')
s_am_raster <- raster(paste0(getwd(),'/Ecoregions_Feb2023/biomes_realms_raster.tif'))
s_am_raster <-
  raster::crop(s_am_raster,
               c(-8.5e6,-3e6,
                 -7e6,3e6))
s_am_raster[s_am_raster %in% c(1)] <- NA # 1 = ocean
s_am_raster[!is.na(s_am_raster)] <- 0



#####
###
# Making function
# This will be parallised later
species_funcion <-
  function(species_list) {
    
    # Getting template rasters
    out_raster <- s_am_raster
    out_raster_richness <- s_am_raster
    
    # Looping across list of species
    for(ss in species_list) {
      
      ### First for habitat availability
      # Getting raster
      species_raster <- raster(ss)
      
      # Checking whether we need to loop over the species
      if(tryCatch(!is.null(raster::crop(species_raster,extent(s_am_raster))), error=function(e) return(FALSE))) { # This returns a true/false statement, TRUE if rasters intersect, FALSE if they do not
        
        # Cropping raster
        species_raster <- raster::crop(species_raster,extent(s_am_raster))
        
        # Getting species richness
        species_richness <- !is.na(species_raster)
        
        # Removing NAs
        species_raster[is.na(species_raster)] <- 0
        species_raster[species_raster<0] <- 0
        
        # Extending both rasters to extent of south america
        species_raster <- raster::extend(species_raster, s_am_raster, value = 0)
        species_richness <- raster::extend(species_richness, s_am_raster, value = 0)
        
        # Adding species raster to the out rasters
        out_raster <- out_raster + species_raster
        out_raster_richness <- out_raster_richness + species_richness
        
        } # End of try catch if statement
     } # End of species loop
    
    # Making list to return object
    out_list <- 
      list(out_raster,
           out_raster_richness)
    
    names(out_list) <- c('outcome','richness')
    
    ### Returning files
    return(out_list)
  }

#####
###
# Now big loop to run this

# Loop through threshold types
# Species types
# Years (only 2020 2050)
# Outcome type (hab area, pop abundance)

for(taxa in taxa_list) { # Looping through taxa
  # Limiting file list...
  file_list_taxa <-
    list.files(paste0(getwd(),'/Outputs/Raster_Outputs/BAU/',taxa),
               full.names = TRUE)
  
  for(thresh in thresh_list) { # Looping through taxa
    # Limiting file list...
    file_list_thresh <-
      file_list_taxa %>%
      .[grepl(thresh,.)]
    
    for(year in c(2020,2050)) { # Looping through years...
      # Limiting file list...
      file_list_year <-
        file_list_thresh %>%
        .[grepl(year,.)]
      
      # Updating scen list
      # Don't need to do all stress scenarios in 2020...
      scen_list <- if(year %in% 2020) scen_list_2020 else scen_list_full
      
      for(scen in scen_list) { # Looping through stress scenarios
        # Limiting file list...
        file_list_scenario <-
          file_list_year %>%
          .[grepl(scen,.)]
        
        for(outcome in outcome_list) { # Looping through type of otucomes
          # Progress tracking...
          cat('Starting:', taxa, thresh, year, scen, outcome,'\n')
          
          # Limiting file list...
          file_list_outcome <-
            file_list_scenario %>%
            .[grepl(outcome,.)]
          
          # Chunking file list into n parts for parallelisation
          chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
          file_list_chunked <- chunk2(file_list_outcome, n = mc_cores)
          
          # Running function
          # This returns a list of rasters
          out_raster_list <- mclapply(file_list_chunked, species_funcion, mc.cores = mc_cores)
          # out_raster_list <- lapply(file_list_chunked[1], species_funcion) # For testing, outside of parallel
          
          # Getting list of rasters we need...
          out_raster_stack <- s_am_raster
          out_raster_richness <- s_am_raster
          
          # Loop through the list and add the 'outcome' and 'richness' layers
          for(ii in 1:length(out_raster_list)) {out_raster_stack <- out_raster_stack + out_raster_list[[ii]][['outcome']]}
          for(ii in 1:length(out_raster_list)) {out_raster_richness <- out_raster_richness + out_raster_list[[ii]][['richness']]}
          
          # Saving rasters
          # writeRaster(out_raster_stack,'/data/ouce-glob2loc/pubh0329/outcome_check.tif')
          # writeRaster(out_raster_richness,'/data/ouce-glob2loc/pubh0329/richness_check.tif')
          writeRaster(out_raster_stack,
                      paste0(getwd(),'/Analyses/Landscape_Analyses/S_America_',
                            thresh,'_',taxa,'_',year,'_',scen,'_',outcome,'_',
                            Sys.Date(),'.tif'),
                      overwrite = TRUE)

          writeRaster(out_raster_richness,
                      paste0(getwd(),'/Analyses/Landscape_Analyses/S_America_',
                            thresh,'_',taxa,'_',year,'_',scen,'_',outcome,'_richness_',
                            Sys.Date(),'.tif'),
                      overwrite = TRUE)
          
          # Clearing files
          rm(out_raster_richness,out_raster_stack,file_list_chunked,file_list_outcome)
          
        } # End loop for outcome type - e.g. habitat area or population abundance
      } # End loop through stress scenarios
    } # End loop through years
  } # End loop through thresholds
} # End big loop through taxa