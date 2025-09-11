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
mc_cores = 10

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
s_am_raster <- rast(paste0(getwd(),'/Ecoregions_Feb2023/biomes_realms_raster.tif'))
s_am_raster <-
  terra::crop(s_am_raster,
               c(-8.5e6,-3e6,
                 -7e6,3e6))
s_am_raster[s_am_raster %in% c(1)] <- NA # 1 = ocean
s_am_raster[!is.na(s_am_raster)] <- 0



#####
###
# Making function
# This will be parallised later
species_funcion <-
  function(index_number) {
    
    # Getting template rasters
    # out_raster <- s_am_raster
    # out_raster_richness <- s_am_raster
    
    # getting species list
    species_list <- file_list_chunked[[index_number]]
    
    # Looping across list of species
    for(ss in species_list) {
      cat(ss,'\n')
      species <-
        gsub(paste0('_',thresh,'.*'),'',ss) %>%
        gsub(paste0('.*',outcome,'_'),'',.)
      
      ### First for habitat availability
      # Getting raster
      species_raster_2020 <- 
        file_list_thresh %>%
        .[grepl(species,.)] %>%
        .[grepl(2020,.)] %>%
        .[grepl('.exp.int.urb',.)] %>%
        .[grepl(outcome,.)] %>%
        raster()
      
      species_raster_2050 <- raster(ss)
      
      intersected_extent <- raster::intersect(extent(species_raster_2020),extent(species_raster_2050))
      
      # Checking whether we need to loop over the species
      if(tryCatch(!is.null(raster::intersect(intersected_extent,extent(s_am_raster))), error=function(e) return(FALSE))) { # This returns a true/false statement, TRUE if rasters intersect, FALSE if they do not
      # if(tryCatch(!is.null(raster::crop(species_raster_2050,extent(s_am_raster))), error=function(e) return(FALSE))) { # This returns a true/false statement, TRUE if rasters intersect, FALSE if they do not
        
        # Cropping raster
        species_raster_2020 <- raster::crop(species_raster_2020,extent(s_am_raster))
        species_raster_2050 <- raster::crop(species_raster_2050,extent(s_am_raster))
        
        # Negative numbers should be 0
        species_raster_2020[species_raster_2020<0] <- NA
        species_raster_2050[species_raster_2050<0] <- NA
        
        # Combined extent of species rasters
        combined_extent <- extent(raster::union(extent(species_raster_2020), extent(species_raster_2050)))
        
        # Extending both
        species_raster_2020 <- raster::extend(species_raster_2020, combined_extent)
        species_raster_2050 <- raster::extend(species_raster_2050, combined_extent)
        
        # Getting species richness - present in either 2020 or 2050
        species_richness_2020 <- !is.na(species_raster_2020)
        species_richness_2050 <- !is.na(species_raster_2050)
        species_richness <- species_richness_2020 + species_richness_2050
        species_richness <- species_richness > 0
        species_richness[species_richness %in% 0] <- NA
        
        
        # Removing NAs
        species_raster_2020[is.na(species_raster_2020) & species_richness %in% 1] <- 0
        # species_raster_2020[species_raster_2020<0] <- 0
        
        species_raster_2050[is.na(species_raster_2050) & species_richness %in% 1] <- 0
        # species_raster_2050[species_raster_2050<0] <- 0
        
        # Adding small number to both - so that we do not get errors/NAs when taking logs
        species_raster_2020 <- species_raster_2020 + 1e-10
        species_raster_2050 <- species_raster_2050 + 1e-10
        
        # Getting log
        species_raster <- log2(species_raster_2050 / species_raster_2020)
        
        # Extending both rasters to extent of south america
        species_raster <- raster::extend(species_raster, s_am_raster)
        species_richness <- raster::extend(species_richness, s_am_raster)
        
        # And writing - because saving as raster grid does not work
        writeRaster(species_raster,
                    paste0(getwd(),'/Analyses/Landscape_Analyses/S_Am_Species_Rasters/',taxa,'_',species,'_',
                           thresh,'_',year,'_',scen,'_',outcome,'_geometric_',
                           Sys.Date(),'.tif'),
                    overwrite = TRUE)
        
        # writeRaster(species_richness,
        #             paste0(getwd(),'/Analyses/Landscape_Analyses/S_Am_Species_Rasters/',taxa,'_',species,'_',
        #                    thresh,'_',year,'_',scen,'_',outcome,'_geometric_richness_',
        #                    Sys.Date(),'.tif'),
        #             overwrite = TRUE)
        } # End of try catch if statement
     } # End of species loop
  } # End of function

#####
###
# Now big loop to run this

# Loop through threshold types
# Species types
# Years (only 2020 2050)
# Outcome type (hab area, pop abundance)

for(taxa in taxa_list[2]) { # Looping through taxa
  # Limiting file list...
  file_list_taxa <-
    list.files(paste0(getwd(),'/Outputs/Raster_Outputs/BAU/',taxa),
               full.names = TRUE)
  
  # For testing
  # file_list_taxa <-
  #   file_list_taxa %>%
  #   .[grepl('exp.int.urb',.)] %>%
  #   .[grepl('Amphibians',.)] %>%
  #   # .[grepl('pop_density',.)] %>%
  #   # .[!grepl('richness',.)] %>%
  #   .[grepl('Brachycephalus|Ischnoch|Adelophryne|Bokermannohyla|Corythomantis|Crossodactylus|Dendropsophus|Hylodes|Scinax|Proceratophrys|Aplastodiscus',.)]
    
  
  for(thresh in thresh_list) { # Looping through taxa
    # Limiting file list...
    file_list_thresh <-
      file_list_taxa %>%
      .[grepl(thresh,.)]
    
    for(year in c(2050)) { # Looping through years...
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
          file_list_chunked <- chunk2(file_list_outcome[130:1000], n = mc_cores)
          
          # Running function
          # This returns a list of rasters
          out_raster_list <- mclapply(1:mc_cores, species_funcion, mc.cores = mc_cores)
          out_raster_list <- lapply(1, species_funcion) # For testing, outside of parallel
          
          # # Getting list of rasters we need...
          # out_raster_stack <- s_am_raster
          # out_raster_richness <- s_am_raster
          # 
          # # Loop through the list and add the 'outcome' and 'richness' layers
          # for(ii in 1:length(out_raster_list)) {out_raster_stack <- out_raster_stack + out_raster_list[[ii]][['outcome']]}
          # for(ii in 1:length(out_raster_list)) {out_raster_richness <- out_raster_richness + out_raster_list[[ii]][['richness']]}
          # 
          # # Saving rasters
          # # writeRaster(out_raster_stack,'/data/ouce-glob2loc/pubh0329/outcome_check.tif')
          # # writeRaster(out_raster_richness,'/data/ouce-glob2loc/pubh0329/richness_check.tif')
          # writeRaster(out_raster_stack,
          #             paste0(getwd(),'/Analyses/Landscape_Analyses/S_America_',
          #                   thresh,'_',taxa,'_',year,'_',scen,'_',outcome,'_',
          #                   Sys.Date(),'.tif'),
          #             overwrite = TRUE)
          # 
          # writeRaster(out_raster_richness,
          #             paste0(getwd(),'/Analyses/Landscape_Analyses/S_America_',
          #                   thresh,'_',taxa,'_',year,'_',scen,'_',outcome,'_richness_',
          #                   Sys.Date(),'.tif'),
          #             overwrite = TRUE)
          
          # Clearing files
          # rm(out_raster_richness,out_raster_stack,file_list_chunked,file_list_outcome)
          
        } # End loop for outcome type - e.g. habitat area or population abundance
      } # End loop through stress scenarios
    } # End loop through years
  } # End loop through thresholds
} # End big loop through taxa