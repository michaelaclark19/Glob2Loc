#!/usr/bin/env Rscript


#####
###
# Looking at connectivity + dispersal between patches
# Taking saved rasters for each species
# Pairing with dispersal distances
# And doing:
# No dispersal
# Half max dispersal (don't need to do if dispersal < 3000)
# Max dispersal (don't need to do if dispersal < 1500)
###
#####


###
# Libraries
library(reticulate)
library(raster)
library(terra)
library(viridis)
library(plyr)
library(dplyr)
library(sf)
library(RColorBrewer)
library(stringr)
library(readr)

###
# python management
# Loading python packages
py_run_string("import numpy as np")
py_run_string("import pandas as pd")
py_run_string("from skimage import io, morphology, measure")
py_run_string("from scipy.ndimage import convolve")
py_run_string("from PIL import Image")
py_run_string("import copy")

###
# setting working directory
setwd('/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity')

###
# functions
source(paste0(getwd(),"/Scripts/Calculating_Biodiversity_Scripts/0.0_ESH_Loss_Functions_2024.08.22.R"))

###
# List of taxa in analysis
taxa.list <- c('Amphibians','Birds','Mammals','Reptiles')

###
# List of species to loop through
species_all <-
  do.call(c,lapply(paste0(getwd(),'/Outputs/CSV_File_Outputs/BAU/',taxa.list),list.files,full.names = TRUE)) %>%
  .[grepl('_completed',.)] %>%
  .[!grepl('Climate_',.)] %>%
  .[grepl('2020|2050',.)] %>%
  gsub('_2020.*|_2050.*','',.) %>%
  unique()

threshold.list <- c('prevalence','specsens')

###
# Getting dispersal distance frames
mam_dispersal_distance_frame = read.csv(paste0(getwd(),"/Other Data Inputs/DispersalDistances_Estimated.csv"))
amp_dispersal_distance_frame = read.csv(paste0(getwd(),"/Other Data Inputs/Amphibian Migration Distance 29Nov2019.csv"))

###
# Creating directories
dir.create(paste0(getwd(),'/Outputs/CSV_File_Outputs_Migration'))
lapply(paste0(getwd(),'/Outputs/CSV_File_Outputs_Migration/',taxa.list),dir.create)

###
# Looping through taxa
#for(tt in sample(taxa.list,length(taxa.list))) {
for(tt in 'Reptiles') {  
  ### 
  # Looping through thresholds
  for(thresh in threshold.list[2]) {
    
    ###
    # Looping through years
    for(yy in c(2020,2050)) {
      
      # Which species are completed
      species_have <-
        list.files(paste0(getwd(),'/Outputs/CSV_File_Outputs_Migration/',tt)) %>%
        # .[!grepl('Climate_',.)] %>%
        .[grepl(thresh,.)] %>%
        .[grepl(yy,.)] %>%
        gsub(paste0('_',thresh,'.*'),'',.) %>%
        unique()
      
      # Hab avail rasters
      hab_avail_rasters <-
        list.files(paste0(getwd(),'/Outputs/Raster_Outputs/BAU/',tt),
                   full.names = TRUE) %>%
        .[!grepl('climate_',.,ignore.case = TRUE)] %>%
        .[grepl(yy,.)] %>%
        .[grepl(thresh,.)] %>%
        .[grepl('.int.urb',.)] %>%
        .[!grepl('pop_dens',.)]
      
      # Pop abund rasters
      pop_dens_rasters <-
        list.files(paste0(getwd(),'/Outputs/Raster_Outputs/BAU/',tt),
                   full.names = TRUE) %>%
        .[!grepl('climate_',.,ignore.case=TRUE)] %>%
        .[grepl(yy,.)] %>%
        .[grepl(thresh,.)] %>%
        .[grepl('.int.urb',.)] %>%
        .[grepl('pop_dens',.)]
      
      
      # Limiting to species we need
      species.list <- 
        species_all[grepl(tt,species_all)] %>%
	gsub(paste0('.*',tt,'/'),'',.) %>%
        .[!(. %in% species_have)] %>%
        unique()

cat('Species Completed: ',length(species_have),'\n')
cat('Species Remaining: ',length(species.list),'\n')
      
cat('Starting: ',tt,' ',thresh,' ',yy,'\n')
      ###
      # Looping through species
      for(ss in sample(species.list,length(species.list))) {
        
        cat('Starting: ',tt,ss,'\n')
        
        species_have <-
          list.files(paste0(getwd(),'/Outputs/CSV_File_Outputs_Migration/',tt)) %>%
          # .[!grepl('Climate_',.)] %>%
          .[grepl(thresh,.)] %>%
          .[grepl(yy,.)] %>%
          gsub(paste0('_',thresh,'.*'),'',.) %>%
          unique()
        
        # And only running if we don't have the specices
        if(!(ss %in% species_have)) {
          # Creating species frame to get migration data, etc
          species.frame <-
            data.frame(Species = ss,
                       Raster = ss,
                       taxon = tt,
                       Have_SDM = 1) %>%
            mutate(tetra_taxon = ifelse(taxon %in% 'Birds','Aves', # Adding taxon identifier for tetra density species
                                        ifelse(taxon %in% 'Mammals','Mammalia',
                                               ifelse(taxon %in% 'Amphibians','Amphibia','Reptilia'))))
          
          # getting body mass (used for patch dispersal)
          species.frame <-
            pop.dens.coef.function(species.frame)
          
          # getting dispersal distance
          k = 1 # from holdover of previous function
          dispersal_distance <-
            dispersal_distance_function(species = species.frame$Species[1], species.frame)
          
          # max limit to dispersal distance - otherwise function takes a crazy amount of time to run
          dispersal_distance <- min(dispersal_distance,25000)
          
          
          # Getting rasters for species
          hab_avail <-
            raster(hab_avail_rasters %>%
                     .[grepl(ss,.)])
          hab_avail[hab_avail < 0] <- 0
          
          pop_dens <-
            raster(pop_dens_rasters %>%
                     .[grepl(ss,.)])
          pop_dens[pop_dens < 0] <- 0
          
          ###
          # Background data management to label patches
          vals_pops_2020 <- getValues(hab_avail)
          trial_ncols = hab_avail@ncols
          
          ###
          # Labelling patches
          py_run_string("d = r.dispersal_distance") # Dispersal distance
          py_run_string("dp = d / 1000 / 1.5") # Numbers of cells
          py_run_string("im = r.vals_pops_2020") # Creating rasters with ones and 0s
          py_run_string("im = np.reshape(im,(-1,r.trial_ncols))") # Creating raster
          py_run_string("im_migrate = copy.copy(im)") # Copying
          py_run_string("im_migrate[np.isnan(im_migrate)] = 0") # NAs to 0s so convolve works
          py_run_string("im_migrate[im_migrate>=.2] = 1") # Adjusting for migration ability
          py_run_string("im_migrate[im_migrate<.2] = 0") # Assumes cell needs >=20% suitable habitat for species for a species to migrate through it
          
          # Assuming there is no migration between habitat patches
          py_run_string("nim = im_migrate != im_migrate.min()")
          py_run_string("bim_no_migrate = morphology.binary_dilation(nim, morphology.disk(0, dtype=bool))")
          py_run_string("labim_no_migrate = measure.label(bim_no_migrate)")
          py_run_string("labim_no_migrate[~nim] = 0")
          
          patches_raster = raster(py$labim_no_migrate, crs = crs(hab_avail))
          extent(patches_raster) = extent(hab_avail)
          # removing non habitat areas
          patches_raster[is.na(hab_avail)] <- NA
          
          ###
          # Labelling populations (w/ half dispersal)
          if(dispersal_distance <= 3000) {
            pops_half_dispersal_raster <- patches_raster
          } else {
            py_run_string("d = r.dispersal_distance/2") # Dispersal distance
            py_run_string("dp = d / 1000 / 1.5") # Numbers of cells
            
            py_run_string("bim = morphology.binary_dilation(nim, morphology.disk(dp, dtype=bool))")
            py_run_string("labim = measure.label(bim)")
            py_run_string("labim[~nim] = 0")
            
            pops_half_dispersal_raster = raster(py$labim, crs = crs(hab_avail))
            extent(pops_half_dispersal_raster) = extent(hab_avail)
            # removing non habitat areas
            pops_half_dispersal_raster[is.na(hab_avail)] <- NA
          }
          
          ###
          # Labelling populations (w/ full dispersal)
          if(dispersal_distance <= 1500) {
            pops_full_dispersal_raster <- patches_raster
          } else {
            py_run_string("d = r.dispersal_distance") # Dispersal distance
            py_run_string("dp = d / 1000 / 1.5") # Numbers of cells
            
            py_run_string("bim = morphology.binary_dilation(nim, morphology.disk(dp, dtype=bool))")
            py_run_string("labim = measure.label(bim)")
            py_run_string("labim[~nim] = 0")
            
            pops_full_dispersal_raster = raster(py$labim, crs = crs(hab_avail))
            extent(pops_full_dispersal_raster) = extent(hab_avail)
            # removing non habitat areas
            pops_full_dispersal_raster[is.na(hab_avail)] <- NA
          }
          
          ### cell area
          cell_area = res(patches_raster)[1] * res(patches_raster)[2] / 1e6
          
          ###
          # Getting summary stats by patch, half population, and population
          out_df <-
            data.frame(patch_id = getValues(patches_raster),
                       half_pop_id = getValues(pops_half_dispersal_raster),
                       pop_id = getValues(pops_full_dispersal_raster),
                       hab_avail = getValues(hab_avail),
                       pop_dens = getValues(pop_dens)) %>%
            mutate(pop_dens = pop_dens * hab_avail) %>%
            dplyr::group_by(patch_id, half_pop_id, pop_id) %>%
            dplyr::summarise(hab_avail = sum(hab_avail, na.rm = TRUE),
                             pop_abund = sum(pop_dens, na.rm = TRUE)) %>%
            mutate(hab_avail = hab_avail * cell_area,
                   pop_abund = pop_abund * cell_area)
          
          # And saving file
          write.csv(out_df,
                    paste0(getwd(),'/Outputs/CSV_File_Outputs_Migration/',tt,'/',ss,'_',thresh,'_',yy,'.csv'),
                    row.names = FALSE)
        }
      }
    }
  }
}
  
