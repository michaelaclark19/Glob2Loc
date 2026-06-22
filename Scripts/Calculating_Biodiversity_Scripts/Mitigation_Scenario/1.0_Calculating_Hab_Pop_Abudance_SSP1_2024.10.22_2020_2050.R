#!/usr/bin/env Rscript
### ESH Script

#-------------------------------------------------------------------------------
# Name:        ESH and agric
# Purpose:      Calculate ESH for species based on range polygons, altitudes
# and land cover associatsions.
#
# Author:      Mike Clark
#
# Created:     21/Dec/2023
# Copyright:   (c) MC 2023
#
#-------------------------------------------------------------------------------

# importing libraries
library(raster)
library(gdalUtils)
library(rgdal)
library(scales)
library(SDMTools)
library(plyr)
library(dplyr)
library(stringr)
library(countrycode)
library(parallel)
library(future)
library(future.apply)
#####
# Python management
# Telling R to get the correct python directory
# Need to do this at the beginning of the script beacuse otherwise R might not recognise python properly
# Sys.setenv(RETICULATE_PYTHON = "/apps/system/easybuild/software/Anaconda3/2022.05/bin/python")
#library(reticulate)

# Loading python packages
#py_run_string("import numpy as np")
#py_run_string("import pandas as pd")
#py_run_string("from skimage import io, morphology, measure")
#py_run_string("from scipy.ndimage import convolve")
#py_run_string("from PIL import Image")
#py_run_string("import copy")


#####
# General management
# Working directory, functions, list of taxa, species, and years

# Setting working directory
setwd("/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity")

# Loading scripts that contain functions to run the analyses
source(paste0(getwd(),"/Scripts/Calculating_Biodiversity_Scripts/0.0_ESH_Loss_Functions_2024.08.22.R"))
source(paste0(getwd(),"/Scripts/Calculating_Biodiversity_Scripts/0.0 Patch Connectivity Functions.R"))

# Directory to write raster outputs - e.g. where the analyses are saved
# Only need to create these if the directory doesn't already exist
if('Raster Outputs' %in% paste0(getwd(),'/Outputs/')) {
  # Do nothing
} else {
  # Create directories
  raster.write.dir <- paste0(getwd(),'/Outputs/Raster_Outputs')
  patch.estimate.dir <- paste0(getwd(),'/Outputs/Patch_Outputs')
  csv.estimate.dir <- paste0(getwd(),'/Outputs/CSV_File_Outputs')
}

# Which scenarios to save?
# Two options here
# 'All' = save all potential combinations
# 'Main' = save single stressors, and all stressors together
# Change these to save amount of time spent forecasting outcomes
scens.save <- 'All'

# Which assumption on climate dispersal from current habitat to future suitable habitat
# Should be one of:
# 'no_limts'
# 'limits'
migration <- climate_migration_function('limits')

# List of SSPs
ssp.list <- list.files(paste0(getwd(),'/ESH_RCPs'), pattern = 'SSP')

# List of taxa
taxa.list <- 
  list.files(paste0(getwd(),'/ESH_RCPs/',ssp.list[1])) %>% 
  .[!grepl('.tif',.)] %>%
  .[!grepl('_',.)]

# List of years to calculate
years.list <- 
  list.files(paste0(getwd(),'/ESH_RCPs/',ssp.list[1],'/',taxa.list[1])) %>%
  str_extract(.,'[0-9]{4,4}') %>%
  unique() %>%
  as.numeric() %>%
  .[.<=2050] %>%
  .[!is.na(.)] %>%
  .[.>=2020]

# List of years in ssp climate forecast
ssp.years <-
  list.files(paste0(getwd(),'/CMIP6_Climate_Data/SSP1-2.6')) %>%
  str_extract_all(.,'[0-9]{4,4}', simplify = TRUE) %>%
  as.data.frame() %>%
  mutate(year_1 = as.numeric(V1),
         year_2 = as.numeric(V2)) %>%
  mutate(mean_year = round_any((year_1+year_2)/2,5)) %>%
  dplyr::select(year_1, year_2, mean_year) %>%
  filter(!is.na(mean_year))

# List of threshold types
thresh.list <- 
  list.files(paste0(getwd(),'/ESH_RCPs/',ssp.list[1],'/',taxa.list[1]), pattern = 'migrate') %>% 
  .[grepl('.tif\\b',.)] %>%
  gsub('_migrate.*','',.) %>%
  gsub('.*[0-9]{4,4}_','',.) %>%
  unique()

#####
# Data input management
# Location of projections of stressors
# And world map - e.g. if world map == NA, then all other layers == NA

# Ag expansion
ag_exp_tifs <- paste0(getwd(),"/Land Forecast Outputs/EAT_Lancet/Global tifs/")

# Ag intensification
ag_int_tifs <- paste0(getwd(),"/Ag Intensity Outputs/Cropland Intensity Forecasts/EAT_Lancet")

# Urbanisation
urb_tifs <- paste0(getwd(),"/Urbanization_Forecasts")

# Climate change
clim_tifs <- paste0(getwd(),"/CMIP6_Climate_Data/", ssp.list)

# Species richness files
species.richness_tifs <- paste0(getwd(),'/Species_Richness_Forecasts/', ssp.list)

# And species richness
if('tot_richness.tif' %in% list.files(paste0(getwd(),"/Other Data Inputs/Species Richness Rasters"))) {
  # Importing
  richness <- raster(paste0(getwd(),"/Other Data Inputs/Species Richness Rasters/tot_richness.tif"))
} else {
  # Creating
  richness <- 
    raster(paste0(getwd(),"/Other Data Inputs/Species Richness Rasters/Amphibians_Richness.tif")) +
    raster(paste0(getwd(),"/Other Data Inputs/Species Richness Rasters/Birds_Richness.tif")) +
    raster(paste0(getwd(),"/Other Data Inputs/Species Richness Rasters/Mammals_Richness.tif")) +
    raster(paste0(getwd(),"/Other Data Inputs/Species Richness Rasters/Reptile_Richness.tif"))
  
  writeRaster(richness,
              paste0(getwd(),"/Other Data Inputs/Species Richness Rasters/tot_richness.tif"))
}

# NPP tif
# This remains constant through climate
# Because we do no thave a good way of estimating how this might change in the future
npp_raster = raster(paste0(getwd(),"/Global Mollweide Maps/NPP_Current.tif"))

# Snap Raster - world map, used to indicate land vs not land
Snap.Raster <- raster(paste0(getwd(),"/ESA LandCov Maps/Crop2010_Corrected_GlobCov.tif"))

# Richness for the time being...
richness[is.na(Snap.Raster)] <- NA # And converting oceans to NAs
richness[richness < 0] <- NA # And values less than 0 to NAs

# And removing
rm(Snap.Raster)


#####
# Files for SDMs/ESH maps
# Getting list of species with SDMs
# And species without SDMs
# Saving these for later
# These are baseline maps of species' habitat area

# These are used to identify where to import species habitat range maps (SDMs or AOH)
# In 2010, and for species with SDMs, both the AOH and SDM maps will be used in calculations
# For a random selection of xxxx species, the non-weighted SDMs will be used for sensitivity analyses
# Getting list of species with SDMs
sdm.species <-
  do.call(c,
          lapply(taxa.list,
                 function(i) {
                   list.files(paste0(getwd(),'/ESH_RCPs/SSP1-2.6/',i),
                              pattern = 'csv') %>%
                     .[grepl('migrate',.)] %>% 
                     paste0(i,'/',.) %>%
                     gsub(migration['sdm.species.search'],'',.)
                 }))

# Getting list of species with AOH maps (e.g. all species)
aoh.species <-
  do.call(c,
          lapply(taxa.list, 
                 function(i) {
                   list.files(paste0(getwd(),'/ESH_Tifs_12Oct/',i)) %>%
                     gsub('_HabEl.*','',.) %>% # Removing nomenclature from the AOH habitats
                     gsub('_P[0-9].*','',.) %>% # Removing nomenclature from the AOH habitats
                     gsub('.tif','',.) %>% # Removing file type
                     paste0(i,'/',.) %>% # Adding taxa back on
                     unique() # Unique observations only
                 })) %>%
  .[!grepl('.aux',.)]


#####
# Data frame management for outputs
# Creating data frames that will be used to store outputs
esh.table <- data.frame(species = NA,
                        subspecies = NA,
                        year = NA,
                        esh_area_orig = NA,
                        esh_loss = NA)


# Data set up for output data frame ----
species.frame <- 
  data.frame(Species = aoh.species,
             Raster = aoh.species) %>%
  mutate(taxon = gsub("/.*",'',Species)) %>% # Identifying taxon
  left_join(.,
            data.frame(Species = sdm.species, Have_SDM = 1)) %>%
  unique() %>%
  mutate(tetra_taxon = ifelse(taxon %in% 'Birds','Aves', # Adding taxon identifier for tetra density species
                              ifelse(taxon %in% 'Mammals','Mammalia',
                                     ifelse(taxon %in% 'Amphibians','Amphibia','Reptilia'))))

# Which species do we already have ----
# Identifying these, and removing from the analsyis
species.have <-
  do.call(c,
          lapply(taxa.list,
                 function(i) {
                   list.files(paste0(getwd(),'/Outputs/CSV_File_Outputs/EAT_Lancet/',i),
                              pattern = 'completed') %>%
                     paste0(i,'/',.)}))

# Flagging these species
species.frame <- 
  species.frame %>%
  mutate(completed = ifelse(Species %in% species.have, 1, 0)) #%>%
#filter(completed %in% 0)


# Checking which species haven't merged...
# which(sdm.species %in% species.frame$Species[!is.na(species.frame$Have_SDM)])

#####
# Loading species specific information 

# First up is:
# Climate coefficient and body mass information

# Climate coefficients
# These are outputs from regression models
# Assessing how the population density 
# Of different types of species responds to climate change

# This is a wrapper function to avoid clutter
# Loads regression model coefficients
# And body mass information
species.frame <-
  pop.dens.coef.function(species.frame)

# Second up are habitat preferences
# These are based on the IUCN habitat classifications
# This is a wrapper function to avoid clutter
aoh.habs <- hab.prefs.function('yay') 

# Third is species traits
# Species traits are based on 'commonality' and 'rarity'
# These are used to assess how species respond to agricultural intensification
species.traits <- species.traits.function('yay') %>% unique() # This is a wrapper function to avoid clutter - see function file for actual code

# Fourth is climate coefficients
coef.climate <-
  rbind(read.csv(paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_fixed_Mammals.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Mammals'),
        read.csv(paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_fixed_Amphibians.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Amphibians'),
        read.csv(paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_fixed_Birds.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Birds'),
        read.csv(paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_fixed_Reptiles.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Reptiles'))


#####
# Importing population density information from TetraDensity
# This is used to set upper and lower estimates of population density
tetra_density = read.csv(paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/TetraDENSITY copy2.csv")) %>%
  mutate(binomial = paste(.$Genus,"_",.$Species, sep = ""))

# Dispersal distance
mam_dispersal_distance_frame = read.csv(paste0(getwd(),"/Other Data Inputs/DispersalDistances_Estimated.csv"))
amp_dispersal_distance_frame = read.csv(paste0(getwd(),"/Other Data Inputs/Amphibian Migration Distance 29Nov2019.csv"))

# Limits species frame to species identified in grepl function...
# species.frame <-
#   species.frame %>%
#   filter(grepl('Panthera_pardus|Acinonyx_jubatus|...',Species))

# Making climate maps for intermediate years

# If statement to check if the climate rasters exist
if(sum(grepl('pvar_raster_2045',list.files(paste0(getwd(),'/CMIP6_Climate_Data/',ssp.list[1],'/Interpolated_Rasters')))) %in% 0) {
  for(y in years.list) {
    cat(y,'\n')
    
    # If already have ssp year - then don't need to do anything
    if(y %in% ssp.years$mean_year) { 
      # Do nothing
    } else { # Do something
      # Years to import
      year_1 <-
        paste0(max(ssp.years$year_1[ssp.years$mean_year < y]),
               '-',
               max(ssp.years$year_2[ssp.years$mean_year < y]))
      
      year_2 <-
        paste0(min(ssp.years$year_1[ssp.years$mean_year > y]),
               '-',
               min(ssp.years$year_2[ssp.years$mean_year > y]))
      
      # Historic or SSP folder?
      if(y < 2030) {ssp_folder_1 <- 'Historic'} else {ssp_folder_1 <- 'SSP1-2.6'}
      
      # Import year 1 rasters
      pwarmest_1 <- 
        list.files(paste0(getwd(),'/CMIP6_Climate_Data/',ssp_folder_1),
                   pattern = year_1,
                   full.names = TRUE) %>%
        paste0(.,'/Managed_Rasters') %>%
        list.files(path = ., pattern = 'Precip.*warm.*Mollweide', full.names = TRUE) %>%
        .[!grepl('squared',.)] %>%
        raster()
      
      pvar_1 <- 
        list.files(paste0(getwd(),'/CMIP6_Climate_Data/',ssp_folder_1),
                   pattern = year_1,
                   full.names = TRUE) %>%
        paste0(.,'/Managed_Rasters') %>%
        list.files(path = ., pattern = 'Var.*precip.*Moll', full.names = TRUE) %>%
        .[!grepl('squared',.)] %>%
        raster()
      
      
      # Import year 2 rasters
      pwarmest_2 <- 
        list.files(paste0(getwd(),'/CMIP6_Climate_Data/SSP1-2.6'),
                   pattern = year_2,
                   full.names = TRUE) %>%
        paste0(.,'/Managed_Rasters') %>%
        list.files(path = ., pattern = 'Precip.*warm.*Mollweide', full.names = TRUE) %>%
        .[!grepl('squared',.)] %>%
        raster()
      
      pvar_2 <- 
        list.files(paste0(getwd(),'/CMIP6_Climate_Data/SSP2-4.5'),
                   pattern = year_2,
                   full.names = TRUE) %>%
        paste0(.,'/Managed_Rasters') %>%
        list.files(path = ., pattern = 'Var.*precip.*Moll', full.names = TRUE) %>%
        .[!grepl('squared',.)] %>%
        raster()
      
      
      
      # Sequence to Interpolate
      int_seq <- 
        seq(from = 1, 
            to = 0, 
            length = (min(ssp.years$mean_year[ssp.years$mean_year > y]) -
                        max(ssp.years$mean_year[ssp.years$mean_year < y])) / 5 +
              1)
      
      year_seq <-
        seq(from = max(ssp.years$mean_year[ssp.years$mean_year < y]),
            to = min(ssp.years$mean_year[ssp.years$mean_year > y]),
            by = 5)
      
      # Interpolating
      pwarmest_raster <-
        pwarmest_1 * int_seq[which(year_seq %in% y)] +
        pwarmest_2 * rev(int_seq)[which(year_seq %in% y)]
      
      pvar_raster <-
        pvar_1 * int_seq[which(year_seq %in% y)] +
        pvar_2 * rev(int_seq)[which(year_seq %in% y)]
      
      # And saving interpolated rasters
      # pwarmest raster
      writeRaster(pwarmest_raster,
                  paste0(getwd(),
                         '/CMIP6_Climate_Data/',
                         ssp.list[1],
                         '/Interpolated_Rasters/pwarmest_raster_',
                         y,
                         '.tif'),
                  overwrite = TRUE)
      
      # pvar raster
      writeRaster(pvar_raster,paste0(getwd(),
                                     '/CMIP6_Climate_Data/',
                                     ssp.list[1],
                                     '/Interpolated_Rasters/pvar_raster_',
                                     y,
                                     '.tif'),
                  overwrite = TRUE)
      
      # Clearing memory
      removeTmpFiles(h = .05)
      gc()
    } # End loop to do something
  } # End loop to interpolate climate rasters
} # End if statement to check whether climate rasters already exist


#####
# Creating directories to save files
# Some files are saved after each species (rasters)
# Some are saved after they reach a certain size (data frames)
dir.create(paste0(getwd(),'/Outputs/'))
dir.create(raster.write.dir)
dir.create(patch.estimate.dir)
dir.create(csv.estimate.dir)

# Creating directories for updated map files
ssps_tmp <- paste0(getwd(),'/ESH_RCPs/',ssp.list,'/')
taxas_tmp <- paste0(taxa.list,'_Updated')

lapply(apply(expand.grid(ssps_tmp, taxas_tmp), 1, paste, collapse=""),
       dir.create)

# Creating folders for interpolated weather rasters
lapply(paste0(getwd(),'/CMIP6_Climate_Data/',ssp.list,'/Interpolated_Rasters'), dir.create)

# Combination of stressors to assess spatially
# Doing this to limit processing time
stressors <-
  c('esh.*abs', # Contribution from individual stressors, ignoring any contribution from other stressors
    'esh.current.extent.esh$', # climate only
    'future.*exp$', # ag expansion only
    'future.*esh.int$', # intensification only
    'future.*esh.urb$', # urban only
    'future.*esh.exp.int.urb') # All

# List of species, scenarios, years, and thresholds that have finished!
# Doing this to avoid repeating analyses that have already been completed
species_completed <-
  do.call(c,
          lapply(
            list.files(paste0(getwd(),'/Outputs/CSV_File_Outputs/EAT_Lancet'),
                       full.names = TRUE) %>%
              .[!grepl('Old',.)],
            list.files, full.names = TRUE, pattern = 'completed')) %>%
  .[grep('csv',.)] %>%
  .[grep('completed',.)] %>%
  .[grepl(migration['species_completed.search'],.)] %>%
  .[grepl('SSP',.)]

# And limiting to only species with an SDM to start...
# And likewise to mammals
species.frame.full <- 
  species.frame %>%
  filter(Have_SDM %in% 1) # %>%
# filter(taxon %in% 'Mammals') 

# And updating binomial - some of these are NAs
species.frame.full <-
  species.frame.full %>%
  mutate(binomial = ifelse(is.na(binomial),
                           gsub('.*/','',Species),
                           binomial))

# Updating body mass info - for species that don't have it
species.frame.full <- 
  body_mass_update(species.frame.full) %>%
  #filter(taxon %in% 'Birds') %>%
  mutate(sort_col = rnorm(nrow(.))) %>%
  arrange(sort_col) %>%
  dplyr::select(-sort_col)

# getting coo data
# Use this to limit which species to incorporate in the analysis
coo_data <-
  read.csv(paste0(getwd(),'/COO_Data/NewCountryOfOccurrenceData.csv')) %>%
  mutate(region = countrycode(ISO3,origin = 'iso3c',destination = 'region')) %>%
  filter(taxon %in% 'Birds') %>%
  # filter(grepl('Vulpes',species)) %>%
  #filter(grepl('America',region,ignore.case=TRUE)) %>%
  #filter(region %in% 'Sub-Saharan Africa') %>%
  dplyr::select(species, region) %>%
  unique()

cat(nrow(coo_data),'\n')

#####
# And creating a loop that goes through species
# This loop calculates biodiversity outcomes for each species, under each potential combination of stressors
# This includes:
# Habitat Fragmentation
# Patch Connectivity
# Habitat area (total, by patch, by cell)
# Population abundance (total, and by patch)

# These are tracked in each time period
# For each potential combination of stressors

# This loops across:
# Years
# Species
# Thresholds (for the SDM models)
# Stressor combinations (climate, ag intensity, ag expansion, urban, fragmentation)

# I've also added a loop for scenarios
# Which will be used in a later project

# This is not in a function
# R stalls when python is called in a function
# Python is used to identify patches
# So need to do this in a loop, rather than a big function

scenarios.list <- 'EAT_Lancet'

###
# And gap filling body mass data where needed
# This has resulted from updated species taxons
cat('Number Species Completed:', length(species_completed %>% .[grepl('specsens',.)]))

###
# Making function for hab availability and pop abundance estimates
species_loop_function <-
  function(k) {
# for(k in 1:nrow(species.frame)) {
  # for(k in cheetah) { # For testing  
  # Progress tracking 
  cat('Starting:',species.frame$binomial[k],"\n")
    
    # doing if else check if species has finished
    species_tmp_check <-
      do.call(c,
              lapply(
                list.files(paste0(csv.estimate.dir,'/EAT_Lancet'),
                           full.names = TRUE) %>%
                  .[!grepl('Old',.)],
                list.files, full.names = TRUE, pattern = 'completed')) %>%
      .[grep('csv',.)] %>%
      .[grep('completed',.)] %>%
      .[grepl(migration['species_completed.search'],.)] %>%
      .[grepl('SSP',.)]
    
    species_tmp_check <-
      species_tmp_check %>%
      gsub('.*CSV_File_Outputs/EAT_Lancet/','',.) %>%
      # gsub(migration['species_completed_year.search'],'',.) %>%
      .[grep(y,.)] %>%
      .[grep(ssp,.)] %>%
      .[grep('specsens',.)] %>%
      gsub('_[0-9]{4,4}','',.) %>%
      gsub('_specsens.*','',.)

    # If else check to save time
    if(species.frame$Species[k] %in% species_tmp_check) {
      # Do nothing
    } else { # Do everything
      
      # cat('Does Species Have SDM: ',species.frame$Have_SDM[k],'\n')
      # Clearing memorgy
      gc()
      
      # Getting location of habitat maps to import
      # AOH map
      aoh_map_file <-
        list.files(paste0(getwd(),'/ESH_Tifs_12Oct/',species.frame$taxon[k]),
                   pattern = gsub(paste0(species.frame$taxon[k],'/'),'',species.frame$Species[k]),
                   full.names = TRUE)
      
      # cat(aoh_map_file)
      
      # SDM map
      # Loading SDM map
      if(species.frame$Have_SDM[k] %in% 1) {
        sdm_map_file <-
          list.files(paste0(getwd(),'/ESH_RCPs/',ssp,'/',species.frame$taxon[k]), 
                     # pattern = gsub(paste0(species.frame$taxon[k],'/'),'',species.frame$Species[k]),
                     # pattern = species.frame$Species[k],
                     full.names = TRUE) %>%
          .[grepl(species.frame$Species[k],.)] %>%
          # .[grepl(paste0(ssp,'_',y,'_'),.)] %>%
          .[grepl(substr(ssp,1,4),.)] %>%
          .[grepl(y,.)] %>%
          .[grepl(migration['sdm_map_file.search'],.)]
      } # End if statement to get file path for SDM maps
      
      # Loop through thresholds for SDMs
      # Using prevalence threshold for BAU
      # And spec sens threshold as a sensitivity analysis
      for(t in thresh.list) {
        
        # Removing temporary files
        # To clear memory space
        # removeTmpFiles(h = .01)
        gc()
        
        
        # Loading aoh map
        aoh_map <- raster(aoh_map_file)
        # If year is 2010, need to load both SDM and ESH raster
        # Assuming both SDM and ESH raster are present
        if(y %in% 2020) {
          
          # Loading SDM map
          if(species.frame$Have_SDM[k] %in% 1) {
            # cat('\n','Importing SDM Map:', sdm_map_file %>%
            #       .[grepl(t,.)],'\n')
            
            # Import SDM Map
            sdm_map <-
              sdm_map_file %>%
              .[grepl(t,.)] %>%
              raster()
            
            # Extend AOH map
            aoh_map <-
              raster::extend(aoh_map,
                             sdm_map,
                             values = NA)
          } else {
            sdm_map <- raster(aoh_map_file)
          }
          
          # And getting 'current map' which is same as 2020 map
          sdm_map_current <- sdm_map
          
        } else { # If not 2020, only need to load either AOH raster or SDM raster
          
          # Loading aoh if its present
          if(species.frame$Have_SDM[k] %in% 1) {
            # Importing sdm map
            sdm_map <-
              sdm_map_file %>%
              .[grepl(t,.)] %>%
              raster()
            
            # Current 2020 map
            sdm_map_file_current <- 
              list.files(paste0(getwd(),'/ESH_RCPs/',ssp,'/',species.frame$taxon[k]), 
                         # pattern = gsub(paste0(species.frame$taxon[k],'/'),'',species.frame$Species[k]),
                         # pattern = species.frame$Species[k],
                         full.names = TRUE) %>%
              .[grepl(species.frame$Species[k],.)] %>%
              # .[grepl(paste0(ssp,'_',y,'_'),.)] %>%
              .[grepl(substr(ssp,1,4),.)] %>%
              .[grepl(2020,.)] %>%
              .[grepl(migration['sdm_map_file.search'],.)] %>%
              .[grepl(t,.)]
            
            # And importing file
            sdm_map_current <- raster(sdm_map_file_current)
            
          } else {
            
            # Importing AOH map
            sdm_map <- raster(aoh_map_file)
            sdm_map_current <- sdm_map
            
          } # End if statement loading habitat maps
        } # End loop checking for 2020 to import species raster
        
        # Doing checks on realm/biome combinations in sdm projections
        if(species.frame$Have_SDM[k] %in% 1) { # Only need to do this if using sdm maps
          # Running a check - this checks to see if realm biome combinations in projected sdm are the same as the realm biome in current aoh
          sdm_map <- 
            realm_biome_check(aoh_file = aoh_map_file,
                              sdm_map_raster = sdm_map,
                              sdm_raster_map_file = sdm_map_file %>% .[grepl(t,.)],
                              species_taxa = species.frame$taxon[k])
          
          # If 2020, only need to do this once
          if(y %in% 2020) {
            sdm_map_current <- sdm_map
          } else { # Need to do the check
            sdm_map_current <-
              realm_biome_check(aoh_file = aoh_map_file,
                                sdm_map_raster = sdm_map_current,
                                sdm_raster_map_file = sdm_map_file_current %>% .[grepl(t,.)],
                                species_taxa = species.frame$taxon[k])
          } # End 2020 check
        } # End realm/biome in sdm projections check
        
        
        # Cropping SDM based on presence/absence
        # Identifying min and max rows in which species exists
        # Doing this to limit size of the map on which calculations are performed - e.g. no reason to do maths on rows and columns where the species doesn't exist!
        # helps speed up performance by quite a bit!
        
        # for future sdm
        which_colsums <- which(colSums(raster::as.matrix(sdm_map),na.rm=TRUE) > 0) # Columns species has habitat
        which_rowsums <- which(rowSums(raster::as.matrix(sdm_map),na.rm=TRUE) > 0) # Rows species has habitat
        
        # for current sdm
        which_colsums_current <- which(colSums(raster::as.matrix(sdm_map_current),na.rm=TRUE) > 0) # Columns species has habitat
        which_rowsums_current <- which(rowSums(raster::as.matrix(sdm_map_current),na.rm=TRUE) > 0) # Rows species has habitat
        
        
        # Adding exception - if no habitat left, just get a 10 x 10 raster for the rest of the script to run
        if(max(which_colsums,which_colsums_current) < 0) {
          # cropping to 10 x 10 raster
          sdm_map <-
            raster::crop(sdm_map,
                         extent(sdm_map,
                                1,10,1,10))
        } else { # Actually cropping sdm map
          # Cropping sdm map
          sdm_map <- 
            raster::crop(sdm_map,
                         extent(sdm_map,
                                min(which_rowsums,which_rowsums_current),
                                max(which_rowsums,which_rowsums_current),
                                min(which_colsums,which_colsums_current),
                                max(which_colsums,which_colsums_current)))
        } # End if statement to make exception for cropping raster
        
        
        # And cropping aoh and current sdm maps
        #aoh_map <- raster::crop(aoh_map, sdm_map)
        sdm_map_current <- raster::crop(sdm_map_current, sdm_map)
        
        # Cropping land cover raster maps
        tmp_crop_raster <- raster::crop(crop_raster, sdm_map)
        tmp_pasture_raster <- raster::crop(pasture_raster,sdm_map)
        tmp_urban_raster <- raster::crop(urban_raster, sdm_map)
        
        # Extending if necessary
        if(extent(tmp_crop_raster)[1] > extent(sdm_map)[1] |
           extent(tmp_crop_raster)[3] > extent(sdm_map)[3] |
           extent(tmp_crop_raster)[2] < extent(sdm_map)[2] |
           extent(tmp_crop_raster)[4] < extent(sdm_map)[4]) {
          # Extending
          tmp_crop_raster <- raster::extend(tmp_crop_raster, sdm_map, value = NA)
          tmp_pasture_raster <- raster::extend(tmp_pasture_raster, sdm_map, value = NA)
          tmp_urban_raster <- raster::extend(tmp_urban_raster, sdm_map, value = NA)
        } # End if statement to extend raster
        
        # Getting land cover maps in 2020
        tmp_crop_2020 <- raster::crop(crop_raster_2020, sdm_map)
        tmp_past_2020 <- raster::crop(pasture_raster_2020, sdm_map)
        tmp_urb_2020 <- raster::crop(urban_raster_2020, sdm_map)
        
        # Getting values of rasters
        # Habitat maps
        sdm_values <- getValues(sdm_map)
        
        # Land cover forecasts
        crop_values <- getValues(tmp_crop_raster)
        pasture_values <- getValues(tmp_pasture_raster)
        pasture_intens_values <- pasture_values
        urb_values <- getValues(tmp_urban_raster)
        urb_intens_values <- urb_values
        
        # Checking if the intens raster  exists
        # Ran into issues with this before, so making an exception
        
        # cat('Intens Raster\n')
        if(!exists('intens_raster')) {
          intens_raster <- raster(list.files(ag_int_tifs, full.names = TRUE) %>% .[grep(paste0(y,'.tif'),.)])
        }
        
        crop_intens_values <- getValues(raster::crop(intens_raster, tmp_crop_raster))
        crop_intens_values[is.na(crop_values)] <- NA
        crop_intens_values[crop_values %in% 0] <- 0
        
        # Land cover 2020
        crop_2020_values <- getValues(tmp_crop_2020)
        past_2020_values <- getValues(tmp_past_2020)
        past_intens_2020_values <- past_2020_values
        urb_2020_values <- getValues(tmp_urb_2020)
        urb_intens_2020_values <- urb_2020_values
        
        
        # cat('Intens Raster 2020\n')
        
        # Checking if the intens raster 2010 exists
        # Ran into issues with this, so making an exception
        if(!exists('intens_raster_2020')) {
          intens_raster_2020 <- raster(list.files(ag_int_tifs, full.names = TRUE) %>% .[grep(paste0('2020.tif'),.)])
        }
        
        crop_intens_2020_values <- getValues(raster::crop(intens_raster_2020, tmp_crop_raster))
        crop_intens_2020_values[is.na(crop_2020_values)] <- NA
        crop_intens_2020_values[crop_2020_values %in% 0] <- 0
        
        
        
        # cat('End Intens Raster Value Getting\n')
        # Getting habitat preferences for each species
        # 1 = suitable habitat
        # .5 = marginal habitat
        # Other = not suitable or marginal
        tmp_crop_pref <- aoh.habs$Crop_Tolerance[aoh.habs$species == species.frame$Species[k]]
        tmp_past_pref <- aoh.habs$Pasture_Tolerance[aoh.habs$species == species.frame$Species[k]]
        tmp_urban_pref <- aoh.habs$Urban_Tolerance[aoh.habs$species %in% species.frame$Species[k]]
        
        # Setting values to 0 if not found
        # max(NULL) returns -Inf
        # So checking on max values, and updating values to 0
        if(max(tmp_crop_pref) < 0) {
          tmp_crop_pref <- 0
        } else {tmp_crop_pref <- max(tmp_crop_pref)}
        if(max(tmp_past_pref) < 0) {
          tmp_past_pref <- 0
        } else {tmp_past_pref <- max(tmp_past_pref)}
        if(max(tmp_urban_pref) < 0) {
          tmp_urban_pref <- 0
        } else {tmp_urban_pref <- max(tmp_urban_pref)}
        
        # And updating to 0s and 1s
        # For this, we only care about whether species can (marginal or suitable) or cannot (other) exist in a habitat
        # These correspond with 1s and 0s
        # So taking ceiling of returned values, as both 1 and .5 will return a value of 1
        # These values will be used later, to identify habitat remaining for each species
        tmp_urban_pref <- ceiling(tmp_urban_pref)
        tmp_crop_pref <- ceiling(tmp_crop_pref)
        tmp_past_pref <- ceiling(tmp_past_pref)
        
        ###
        # Now Updating values for species responses to different land covers
        # These are based on regressions recreated from Sykes et al
        # That estimate how species with different traits
        # Respond to different types of human modified habitats
        
        # If species cannot exist in crop or pasture
        # Then don't need to worry about crop or pasture intensity
        
        # Cropland first
        # cat('Before Cropland Intens Value Getting\n')
        if(tmp_crop_pref %in% 0) {
          crop_intens_2020_values <- rep(1, length(crop_intens_2020_values))
          crop_intens_values <- rep(1, length(crop_intens_values))
        } else { # Need to update values for species to exist in croplands of different intensities
          # High intensity cropland - counterfactual in 2010
          crop_intens_2020_values[crop_intens_2020_values %in% c(2,3)] <-
            species.traits$Mod_High_Cropland[paste0(species.traits$Class,'/',species.traits$binomial) %in% 
                                               paste0(species.frame$tetra_taxon[k], '/', species.frame$binomial[k]) &
                                               species.traits$estimate_type %in% 'Mean']
          # Projection in other years
          # species.traits[paste0(species.traits$Class,'/',species.traits$binomial) %in% 
          #                                    paste0(species.frame$tetra_taxon[k], '/', species.frame$binomial[k]) &
          #                                    species.traits$estimate_type %in% 'Mean',]
          
          # High intensity cropland - projection
          crop_intens_values[crop_intens_values %in% c(2,3)] <-
            species.traits$Mod_High_Cropland[paste0(species.traits$Class,'/',species.traits$binomial) %in% 
                                               paste0(species.frame$tetra_taxon[k], '/', species.frame$binomial[k]) &
                                               species.traits$estimate_type %in% 'Mean']
          
          # Low intensity cropland - counterfactual in 2010
          crop_intens_2020_values[crop_intens_2020_values %in% c(1)] <-
            species.traits$Low_Cropland[paste0(species.traits$Class,'/',species.traits$binomial) %in% 
                                          paste0(species.frame$tetra_taxon[k], '/', species.frame$binomial[k]) &
                                          species.traits$estimate_type %in% 'Mean']
          
          # Low intensity cropland - projection
          crop_intens_values[crop_intens_values %in% c(1)] <-
            species.traits$Low_Cropland[paste0(species.traits$Class,'/',species.traits$binomial) %in% 
                                          paste0(species.frame$tetra_taxon[k], '/', species.frame$binomial[k]) &
                                          species.traits$estimate_type %in% 'Mean']
          
          # And removing NAs
          crop_intens_2020_values[is.na(crop_intens_2020_values)] <- 0
          crop_intens_values[is.na(crop_intens_values)] <- 0
          
          #crop_intens_2020_values[crop_intens_2020_values < 0] <- 0
          #crop_intens_values[crop_intens_values < 0] <- 0
          
        } # End if statement for cropland intensity values
        
        
        # Pastureland second
        # cat('Before Pastureland Intens Responses\n')
        if(tmp_past_pref %in% 0) {
          past_intens_2020_values <- rep(1, length(past_intens_2020_values))
          pasture_intens_values <- rep(1, length(pasture_intens_values))
        } else { # Need to update value for species' capacity to exist in pastureland
          
          # Pasture intensity - counterfactual in 2010
          past_intens_2020_values[past_2020_values > 0] <- 
            species.traits$Pasture[paste0(species.traits$Class,'/',species.traits$binomial) %in% 
                                     paste0(species.frame$tetra_taxon[k], '/', species.frame$binomial[k]) &
                                     species.traits$estimate_type %in% 'Mean']
          
          # Pasture intensity - projection
          pasture_intens_values[pasture_values > 0] <- 
            species.traits$Pasture[paste0(species.traits$Class,'/',species.traits$binomial) %in% 
                                     paste0(species.frame$tetra_taxon[k], '/', species.frame$binomial[k]) &
                                     species.traits$estimate_type %in% 'Mean']
          # past.intens.tmp[past.intens.2010%in%0] <- 0
          # And removing NAs
          past_intens_2020_values[is.na(past_intens_2020_values)] <- 0
          pasture_intens_values[is.na(pasture_intens_values)] <- 0
          
          #past_intens_2020_values[past_intens_2020_values < 0] <- 0
          #pasture_intens_values[pasture_intens_values < 0] <- 0
        } # End if statement for pastureland intensity values
        
        # Urban land cover third
        # cat('Before Urban Land Cover Responses\n')
        if(tmp_urban_pref %in% 0) { 
          urb_intens_2020_values <- rep(1, length(urb_2020_values))
          urb_intens_values <- rep(1, length(urb_values))
        } else { # Need to update value for species' capacity to exist in Urban
          
          # Urban intensity - counterfactual
          urb_intens_2020_values[urb_2020_values > 0] <- 
            species.traits$Urban[paste0(species.traits$Class,'/',species.traits$binomial) %in% 
                                   paste0(species.frame$tetra_taxon[k], '/', species.frame$binomial[k]) &
                                   species.traits$estimate_type %in% 'Mean']
          
          # Urban intensity - projection
          urb_intens_values[urb_values > 0] <- 
            species.traits$Urban[paste0(species.traits$Class,'/',species.traits$binomial) %in% 
                                   paste0(species.frame$tetra_taxon[k], '/', species.frame$binomial[k]) &
                                   species.traits$estimate_type %in% 'Mean']
          
          # And removing NAs
          urb_intens_2020_values[is.na(urb_2020_values)] <- 0
          urb_intens_values[is.na(urb_values)] <- 0
          
          #urb_intens_2020_values[urb_intens_2020_values < 0] <- 0
          #urb_intens_values[urb_intens_values < 0] <- 0
        } # End if statement for pastureland intensity values
        
        ### Getting cell by cell estimate of (1) remaining habitat, and (2) population abundance modifiers 
        # Converting it all into a big data frame for maths
        # Easier to do it this way
        # Than to do the math on a series of vectors
        # Converting all of these into a data frame
        # Then will calculate change in ESH area, pop size, and suitable habitats from here
        
        
        ### Checking sdm maps at this point
        # cat('Checking SDM Maps')
        # cat(sum(getValues(sdm_map),na.rm=TRUE),'\n')
        # cat(sum(getValues(sdm_map_current),na.rm=TRUE),'\n')
        
        
        
        # Remaining Habitat Area
        # For Climate Change
        species_df_climate <-
          data.frame(esh.values = sdm_values,
                     crop.2010 = crop_2020_values,
                     past.2010 = past_2020_values,
                     urban.2010.values = urb_2020_values,
                     crop.values = crop_values,
                     pasture.values = pasture_values,
                     urban.values = urb_values,
                     crop.pref = tmp_crop_pref,
                     past.pref = tmp_past_pref,
                     crop.intens.2010 = rep(1-tmp_crop_pref, length(crop_intens_2020_values)),
                     past.intens.2010 = rep(1-tmp_past_pref, length(past_intens_2020_values)),
                     urb.intens.2010 = rep(1-tmp_urban_pref, length(urb_intens_2020_values)),
                     crop.intens.tmp = rep(1-tmp_crop_pref, length(crop_intens_values)),
                     past.intens.tmp = rep(1-tmp_past_pref, length(pasture_intens_values)),
                     urb.intens.tmp = rep(1-tmp_urban_pref, length(urb_intens_values)))
        
        # Converting 0s and 1s to NAs if outside the species' habitat range
        # This ensures pressures outside a species' habitat range are not added into the calcs
        species_df_climate[(species_df_climate$esh.values %in% 0),
                           !grepl('esh.values',names(species_df_climate))] <- NA
        
        # Repeating for no climate change
        species_df_noclimate <- 
          data.frame(esh.values = getValues(sdm_map_current),
                     crop.2010 = crop_2020_values,
                     past.2010 = past_2020_values,
                     urban.2010.values = urb_2020_values,
                     crop.values = crop_values,
                     pasture.values = pasture_values,
                     urban.values = urb_values,
                     crop.pref = tmp_crop_pref,
                     past.pref = tmp_past_pref,
                     crop.intens.2010 = rep(1-tmp_crop_pref, length(crop_intens_2020_values)),
                     past.intens.2010 = rep(1-tmp_past_pref, length(past_intens_2020_values)),
                     urb.intens.2010 = rep(1-tmp_urban_pref, length(urb_intens_2020_values)),
                     crop.intens.tmp = rep(1-tmp_crop_pref, length(crop_intens_values)),
                     past.intens.tmp = rep(1-tmp_past_pref, length(pasture_intens_values)),
                     urb.intens.tmp = rep(1-tmp_urban_pref, length(urb_intens_values)))
        
        # Converting 0s and 1s to NAs if outside the species' habitat range
        # This ensures pressures outside a species' habitat range are not added into the calcs
        species_df_noclimate[(species_df_noclimate$esh.values %in% 0),
                             !grepl('esh.values',names(species_df_noclimate))] <- NA
        
        # cat('Before Non Spatial AOH Calcs\n')	
        # Getting non-spatial habitat loss estimates
        # This is a high-level summary
        # Across the species' entire habitat range
        species_df_climate <-
          esh.loss.function(current = 'yes', # Estimate current ESH area?
                            stresses = c('intensification','expansion','urbanization'), # What stresses to estimate?
                            species.df = species_df_climate, # input data frame
                            tmp.crop.pref = tmp_crop_pref, # tolerance to cropland (as definied by suitability in IUCN)
                            tmp.past.pref = tmp_past_pref, # tolerance to pastureland (as definied by suitability in IUCN)
                            tmp.urban.pref = tmp_urban_pref)  # tolerance to urban (as defined by IUCN)
        
        # Repeating for no climate change
        species_df_noclimate <-
          esh.loss.function(current = 'yes', # Estimate current ESH area?
                            stresses = c('intensification','expansion','urbanization'), # What stresses to estimate?
                            species.df = species_df_noclimate, # input data frame
                            tmp.crop.pref = tmp_crop_pref, # tolerance to cropland (as definied by suitability in IUCN)
                            tmp.past.pref = tmp_past_pref, # tolerance to pastureland (as definied by suitability in IUCN)
                            tmp.urban.pref = tmp_urban_pref)  # tolerance to urban (as defined by IUCN)
        
        # And summarising remaining habitat area for entire species
        # Under each combination of stressors
        
        # Sum by column
        out_sum_climate <- colSums(species_df_climate,na.rm=TRUE) 
        # Transposing
        out_sum_climate <- t(out_sum_climate)
        # Converting to data frame
        out_sum_climate <- 
          out_sum_climate %>% 
          as.data.frame(.) %>% 
          mutate(species = species.frame$Raster[k], 
                 taxon = species.frame$taxon[k])
        
        out_sum_climate <- out_sum_climate %>% mutate(year = y)
        
        # Repeating for no climate
        # Sum by column
        out_sum_noclimate <- colSums(species_df_noclimate,na.rm=TRUE) 
        # Transposing
        out_sum_noclimate <- t(out_sum_noclimate)
        # Converting to data frame
        out_sum_noclimate <- 
          out_sum_noclimate %>% 
          as.data.frame(.) %>% 
          mutate(species = species.frame$Raster[k], 
                 taxon = species.frame$taxon[k])
        
        out_sum_noclimate <- out_sum_noclimate %>% mutate(year = y)
        
        # Changing names to work with rest of code
        final.dat.df <- 
          rbind(out_sum_climate %>% mutate(climate_change = 'yes'),
                out_sum_noclimate %>% mutate(climate_change = 'no'))
        
        ### Getting non-spatial population density estimates
        # This is a high-level summary
        # Across the species' entire habitat range
        
        # Logic check - if species cannot exist in any of the habitats, then non-spatial abundance will be same as non-spatial habitat
        if(tmp_crop_pref %in% 0 & tmp_past_pref %in% 0 & tmp_urban_pref %in% 0) {
          # Will be exactly the same number of cells as habitat area
          species_df_population_climate <- species_df_climate
          species_df_population_noclimate <- species_df_noclimate
          
          # And rbinding - will be exactly same number of cells
          final.dat.df.pop <-
            rbind(out_sum_climate %>% mutate(climate_change = 'yes'),
                  out_sum_noclimate %>% mutate(climate_change = 'no'))
          
          
        } else { # Creating population abundance modifiers
          
          # Population Abundance Modifiers - with climate dispersal
          species_df_population_climate <-
            data.frame(esh.values = sdm_values,
                       crop.2010 = crop_2020_values,
                       past.2010 = past_2020_values,
                       urban.2010.values = urb_2020_values,
                       crop.values = crop_values,
                       pasture.values = pasture_values,
                       urban.values = urb_values,
                       crop.pref = tmp_crop_pref,
                       past.pref = tmp_past_pref,
                       crop.intens.2010 = crop_intens_2020_values,
                       past.intens.2010 = past_intens_2020_values,
                       urb.intens.2010 = urb_intens_2020_values,
                       crop.intens.tmp = crop_intens_values,
                       past.intens.tmp = pasture_intens_values,
                       urb.intens.tmp = urb_intens_values)
          
          # Converting to NAs
          species_df_population_climate[(species_df_population_climate$esh.values %in% 0),
                                        !grepl('esh.values',names(species_df_population_climate))] <- NA
          
          # Population Abundance Modifiers - without climate dispersal
          species_df_population_noclimate <-
            data.frame(esh.values = getValues(sdm_map_current),
                       crop.2010 = crop_2020_values,
                       past.2010 = past_2020_values,
                       urban.2010.values = urb_2020_values,
                       crop.values = crop_values,
                       pasture.values = pasture_values,
                       urban.values = urb_values,
                       crop.pref = tmp_crop_pref,
                       past.pref = tmp_past_pref,
                       crop.intens.2010 = crop_intens_2020_values,
                       past.intens.2010 = past_intens_2020_values,
                       urb.intens.2010 = urb_intens_2020_values,
                       crop.intens.tmp = crop_intens_values,
                       past.intens.tmp = pasture_intens_values,
                       urb.intens.tmp = urb_intens_values)
          
          # Converting to NAs
          species_df_population_noclimate[(species_df_population_noclimate$esh.values %in% 0),
                                          !grepl('esh.values',names(species_df_population_noclimate))] <- NA
          
          
          
          # And getting summary - with climate change
          species_df_population_climate <-
            esh.loss.function(current = 'yes', # Estimate current ESH area?
                              stresses = c('intensification','expansion','urbanization'), # What stresses to estimate?
                              species.df = species_df_population_climate, # input data frame
                              tmp.crop.pref = tmp_crop_pref, # tolerance to cropland (as definied by suitability in IUCN)
                              tmp.past.pref = tmp_past_pref, # tolerance to pastureland (as definied by suitability in IUCN)
                              tmp.urban.pref = tmp_urban_pref)  # tolerance to urban (as defined by IUCN)
          
          # And getting summary - without climate change
          species_df_population_noclimate <-
            esh.loss.function(current = 'yes', # Estimate current ESH area?
                              stresses = c('intensification','expansion','urbanization'), # What stresses to estimate?
                              species.df = species_df_population_noclimate, # input data frame
                              tmp.crop.pref = tmp_crop_pref, # tolerance to cropland (as definied by suitability in IUCN)
                              tmp.past.pref = tmp_past_pref, # tolerance to pastureland (as definied by suitability in IUCN)
                              tmp.urban.pref = tmp_urban_pref)  # tolerance to urban (as defined by IUCN)
          
          # Sum by column - with climate
          out_sum_pop_climate <- colSums(species_df_population_climate,na.rm=TRUE) 
          # Transposing
          out_sum_pop_climate <- t(out_sum_pop_climate)
          # Converting to data frame
          out_sum_pop_climate <- 
            out_sum_pop_climate %>% 
            as.data.frame(.) %>% 
            mutate(species = species.frame$Raster[k], 
                   taxon = species.frame$taxon[k])
          
          out_sum_pop_climate <- out_sum_pop_climate %>% mutate(year = y)
          
          # Sum by column - without climate
          out_sum_pop_noclimate <- colSums(species_df_population_noclimate,na.rm=TRUE) 
          # Transposing
          out_sum_pop_noclimate <- t(out_sum_pop_noclimate)
          # Converting to data frame
          out_sum_pop_noclimate <- 
            out_sum_pop_noclimate %>% 
            as.data.frame(.) %>% 
            mutate(species = species.frame$Raster[k], 
                   taxon = species.frame$taxon[k])
          
          out_sum_pop_noclimate <- out_sum_pop_noclimate %>% mutate(year = y)
          
          # if(is.na(species.frame$Have_SDM[k])) {
          #   # out_sum[,names(final.dat.df)[!(names(final.dat.df) %in% names(out_sum))]] <- NA
          #   out_sum <- out_sum[,names(final.dat.df)]
          # }
          
          final.dat.df.pop <-
            rbind(out_sum_pop_climate %>% mutate(climate_change = 'yes'),
                  out_sum_pop_noclimate %>% mutate(climate_change = 'no'))
        } # End statement to make non-spatial modified population abundance estimate
        
        ### Saving csv files
        #write.csv(species_df_climate,'/data/ouce-glob2loc/pubh0329/glob2loc_checks/Axis_porcinus_species_df.csv')
        #write.csv(species_df_noclimate,'/data/ouce-glob2loc/pubh0329/glob2loc_checks/Axis_porcinus_species_df_noclimate.csv')
        #write.csv(species_df_population_noclimate,'/data/ouce-glob2loc/pubh0329/glob2loc_checks/Axis_porcinus_species_df_population_noclimate.csv')
        #write.csv(species_df_population_climate,'/data/ouce-glob2loc/pubh0329/glob2loc_checks/Axis_porcinus_species_df_population_climate.csv')
        
        
        
        
        # Getting dispersal distance for species
        #dispersal_distance <-
        #  dispersal_distance_function(species = species.frame$Species[k], species.frame)
        
        # Getting population density estimates for species
        # First part is data management
        # - Need to limit estimate pop density to within bounds (+- a bit) of observed densities
        
        # Min and max placeholders
        hab_min = 0
        hab_max = 0
        
        # Min will always be 0 individuals
        # Identifying max
        # This comes from the tetra density database
        
        # Checking if there is data from the species
        hab_max = max(tetra_density$Density[tetra_density$binomial %in% species.frame$Species[k]], na.rm = TRUE) * 1.2
        
        # If no data from the species, then going to the genus
        # updating for genus if not binomial
        if(!is.finite(hab_max) | hab_min == hab_max) {
          hab_max = max(tetra_density$Density[tetra_density$Genus %in% species.frame$genus[k]], na.rm = TRUE) * 1.2
        }
        
        # If no data from the genus, then going to the family
        if(!is.finite(hab_max) | hab_min == hab_max) {
          hab_max = max(tetra_density$Density[tetra_density$Family %in% str_to_title(species.frame$family[k])], na.rm = TRUE) * 1.2
        } 
        
        # And if no data from the family, then using the order
        if(!is.finite(hab_max) | hab_min == hab_max) {
          hab_max = max(tetra_density$Density[tetra_density$Order %in% str_to_title(species.frame$order[k])], na.rm = TRUE) * 1.2
        } 
        
        # And repeating for class if not order
        if(!is.finite(hab_max) | hab_min == hab_max) {
          # hab_min = min(tetra_density$Density[tetra_density$Family %in% str_to_title(species.frame$family[k])], na.rm = TRUE)
          hab_max = max(tetra_density$Density[tetra_density$Class %in% str_to_title(species.frame$tetra_taxon[k])], na.rm = TRUE) * 1.2
          # density.frame = density.frame %>% mutate(Species = species.frame$Species[k], Density_Min_Max = 'Family') %>%
          #   mutate(Density = squish(Density, range = c(hab_min, hab_max)))
        } 
        
        
        # Estimating population densities
        # This is done by recreating the framework from Santini et al 
        # The output from this is a raster map
        # Containing estimates of population density per sq km for the species
        # Doing this for both climate and no climate
        # cat('Before Pop Density Raster\n')
        
        density_raster_climate <- 
          pop_density_function(coefs = coef.climate %>% filter(taxon %in% species.frame$taxon[k]),
                               coef_order = species.frame$order_intercept_adj_weighted[k],
                               coef_family = species.frame$family_intercept_adj_weighted[k],
                               coef_binomial = species.frame$binomial_intercept_adj_weighted[k],
                               body_mass_grams = species.frame$est_mass_kg[k]*1000,
                               npp_input = crop(npp_raster, sdm_map),
                               pcv_input = crop(pvar_raster, sdm_map),
                               pwarmest_input = crop(pwarmest_raster, sdm_map),
                               richness_input = crop(richness, sdm_map),
                               hab_min = hab_min,
                               hab_max = hab_max,
                               nrows = sdm_map@nrows)
        # Updating crs
        crs(density_raster_climate) <- crs(sdm_map)
        extent(density_raster_climate) <- extent(sdm_map)
        
        # Now with current climate
        density_raster_noclimate <- 
          pop_density_function(coefs = coef.climate %>% filter(taxon %in% species.frame$taxon[k]),
                               coef_order = species.frame$order_intercept_adj_weighted[k],
                               coef_family = species.frame$family_intercept_adj_weighted[k],
                               coef_binomial = species.frame$binomial_intercept_adj_weighted[k],
                               body_mass_grams = species.frame$est_mass_kg[k]*1000,
                               npp_input = crop(npp_raster, sdm_map_current),
                               pcv_input = crop(pvar_raster_2020, sdm_map_current),
                               pwarmest_input = crop(pwarmest_raster_2020, sdm_map_current),
                               richness_input = crop(richness, sdm_map_current),
                               hab_min = hab_min,
                               hab_max = hab_max,
                               nrows = sdm_map_current@nrows)
        # Updating crs
        crs(density_raster_noclimate) <- crs(sdm_map_current)
        extent(density_raster_noclimate) <- extent(sdm_map_current)
        
        # cat('After Pop Density Raster\n')
        
        
        # saving density rasters
        #writeRaster(density_raster_climate,paste0('/data/ouce-glob2loc/pubh0329/glob2loc_checks/Axis_porcinus_Density_Climate_',y,'.tif'),overwrite=TRUE)
        #writeRaster(density_raster_noclimate,paste0('/data/ouce-glob2loc/pubh0329/glob2loc_checks/Axis_porcinus_Density_pClimate_',y,'.tif'),overwrite=TRUE)
        # writeRaster(sdm_map,paste0('/data/ouce-glob2loc/pubh0329/glob2loc_checks/',species.frame$binomial[k],'_2050.tif'),overwrite=TRUE)
        # writeRaster(sdm_map_current,paste0('/data/ouce-glob2loc/pubh0329/glob2loc_checks/',species.frame$binomial[k],'_2020.tif'),overwrite=TRUE)
        
        # Getting non-spatial density estimates
        # This is used for post-hoc analysis
        # And to save approximate population density estimates
        # Across a species' entire habitat range
        density.frame.climate <- 
          density.function.new(esh.list = stack(sdm_map),
                               density.list = stack(density_raster_climate)) %>%
          mutate(species = species.frame$Species[k], scenario = i, year = y, climate_change = 'yes')
        
        density.frame.noclimate <- 
          density.function.new(esh.list = stack(sdm_map_current),
                               density.list = stack(density_raster_noclimate)) %>%
          mutate(species = species.frame$Species[k], scenario = i, year = y, climate_change = 'no')
        
        # Can't properly do density estimates if no information on a species' body mass
        if(is.na(species.frame$est_mass_kg[k])) {
          density.frame.climate[,which(names(density.frame.climate) %in% 'mean_density') : which(names(density.frame.climate) %in% 'density_95th')] <- 'No body mass data'
          density.frame.noclimate[,which(names(density.frame.noclimate) %in% 'mean_density') : which(names(density.frame.noclimate) %in% 'density_95th')] <- 'No body mass data'
        }
        
        # cat('Binding Density Data Frames \n')
        # And rbinding to save outcomes
        final.density.frame <- rbind(final.density.frame, density.frame.climate, density.frame.noclimate)
        
        # List of scenarios to loop through
        # These correspond with columns of species_df
        scen.list <- grep('esh.*abs|esh.*future', names(species_df_climate))
        
        # Number of columns in the esh raster
        # Need this to make the python script work
        trial_ncols = sdm_map@ncols
        
        # Creating empty raster stacks
        # These are saved after the analysis for each species is conducted
        tmp.out.stack <- stack() # For population estimates
        tmp.loss.stack <- stack() # For remaining habitat estimates
        
        # And looping through scenarios to get spatial estimates of biodiversity outcomes
        # This loops gets, for each habitat patch:
        # Habitat fragmentation
        # Geographic location and area
        # Population abundance
        
        # Only need scenarios for individual stressors and all stressors
        scen.list <- 
          grep(paste(stressors,collapse = '|'),names(species_df_climate))# %>%
        #.[grepl('current|crop|int.urb',.)]
        
        scen.list <- scen.list[c(5,6,9,10,11)]
        #print(names(species_df_climate)[scen.list])
        # cat(scen.list,'\n')
        
        # Which frame to use for the spatial analyses?
        #if(dispersal_distance > 200000) {
        # Writing file
        #write.csv(data.frame(taxa = species.frame$taxon[k],
        #                     binomial = species.frame$binomial[k],
        #                     used_dispersal_distance = 200000,
        #                     actual_dispersal_distance = dispersal_distance),
        #          paste0(csv.estimate.dir,"/",i,'/',species.frame$taxon[k],"/New_migration_distance_",
        #                 species.frame$binomial[k],
        #                 ".csv"))
        
        # Updating dispersal distance
        #  dispersal_distance <- 200000
        #}
        
        for(scen in scen.list) { # for testing
          # for(scen in scen.list) {
          
          # Removing temporary files
          # To clear memory space
          # removeTmpFiles(h = .01)
          gc()
          
          # Progress tracking for sanity
          scen_r = names(species_df_climate)[scen]
          cat(paste(scen_r,Sys.time(),"\n"))
          
          # Getting correct data frame containing population abundance and habitat density
          if(grepl('current|esh.exp.int.urb',scen_r)) {
            species_df <- species_df_climate
            species_df_pop <- species_df_population_climate
            density_raster_scen <- density_raster_climate
          } else {
            species_df <- species_df_noclimate
            species_df_pop <- species_df_population_noclimate
            density_raster_scen <- density_raster_noclimate
          }
          
          # Making area of remaining habitat availability
          species.raster.loss = 
            raster(matrix(species_df[,scen_r],
                          nrow = sdm_map@nrows,
                          byrow = TRUE))
          crs(species.raster.loss) = crs(sdm_map)
          extent(species.raster.loss) = extent(sdm_map)
          
          # Identifying habitat patches
          # Input to this is a binary raster, where 0 indicates not in a patch; and 1 indicates in a patch
          # stats for remaining population outside of patches is calculated later
          species.raster.loss.tmp <- species.raster.loss > .2
          new.patch_no_migrate = raster::clump(species.raster.loss.tmp)
          extent(new.patch_no_migrate) = extent(sdm_map)
          
          # Habitat patch, with migration raster
          # new.patch = raster(py$labim, crs = crs(sdm_map))
          # extent(new.patch) = extent(sdm_map)
          
          # Updating density raster, to be specific to the scenario, to account for the abundance responses to human modified habitats
          # Getting values
          scen_density_raster <-
            species_df_pop[,scen_r] *
            getValues(density_raster_scen)
          

    #cat('CHECKING VALUES: Population Multipliers\n')
    #print(range(species_df_pop[,scen_r],na.rm=TRUE))

          # Checking for negative values before converting into raster
          # Can be negative because some of these use land covers from different time periods
          # E.g. pasture in 2010, with crop in 2050, when crop expanded into existing pasture
          
          # if(min(scen_density_raster,na.rm=TRUE) < 0) {
          #   scen_density_raster[scen_density_raster < 0] <- 0
          # }
          
          # Converting back to raster
          scen_density_raster <- raster(matrix(scen_density_raster, nrow =  sdm_map@nrows, byrow = TRUE))
          crs(scen_density_raster) <- crs(sdm_map)
          extent(scen_density_raster) <- extent(sdm_map)
          
          # Saving raster outputs
          #writeRaster(scen_density_raster,paste0(getwd(),'/Outputs/Raster_Outputs/EAT_Lancet/',species.frame$taxon[k],'/pop_density_',species.frame$binomial[k],'_',t,'_',y,'_',scen_r,'.tif'),overwrite=TRUE)
          #writeRaster(species.raster.loss,paste0(getwd(),'/Outputs/Raster_Outputs/EAT_Lancet/',species.frame$taxon[k],'/hab_availability_',species.frame$binomial[k],'_',t,'_',y,'_',scen_r,'.tif'),overwrite=TRUE)
          
          # Getting stats by patch
          # This includes habitat area and population abundance
          # This does not yet account for dispersal between nearby patches - this occurs in a later script
          # This calculates habitat area/etc for locations outside of habitat patches
          if(!is.na(maxValue(new.patch_no_migrate))) {
            save.df.spatial.pop <-
              pop_estimate_spatial(density_raster = scen_density_raster,
                                   patch_no_migrate = new.patch_no_migrate,
                                   #patch_migrate = new.patch_no_migrate,
                                   habitat_availability = species.raster.loss) %>%
              mutate(scenario = scen_r) %>%
              mutate(binomial = species.frame$Species[k]) %>%
              mutate(have_sdm = species.frame$Have_SDM[k], taxon = species.frame$taxon[k], year = y, land_scen = i)
          } else {
            # creating empty df
            save.df.spatial.pop <- data.frame()
          }
          
          
          # And rbinding stats by patch
          if(nrow(save.df.spatial.pop) > 0) {
            if(is.na(species.frame$est_mass_kg[k])) { # Can't estimate population sizes without body mass estimates
              save.df.spatial.pop$estimated_population <- 'No body mass data'
            }
            save.spatial.pop.final <-
              rbind(save.spatial.pop.final,
                    save.df.spatial.pop)
          } else {
            save.spatial.pop.final <-
              rbind(save.spatial.pop.final,
                    data.frame(patch_no_migrate = 'no remaining habitat',
                               #patch_migrate = 'no remaining habitat',
                               estimated_population = 'no remaining habitat',
                               habitat_availability_sqkm = 'no remaining habitat',
                               scenario = scen_r, binomial = species.frame$Species[k], have_sdm = species.frame$Have_SDM[k], taxon = species.frame$taxon[k],
                               year = y, land_scen = i))
          }
          
          # Removing habitat availability and pop density raster
          rm(species.raster.loss, scen_density_raster)
          
        } # End of loop for scenario combinations of different stressors within a single species
        
        
        # Removing items to clear memory space
        # Python objects first
        #py_run_string("input_array, sdm_buffer, esh_buffer, labim, bim, labim_no_migrate, bim_no_migrate, nim, im, im_migrate, buffer_raster = 'NONE', 'NONE', 'NONE', 'NONE', 'NONE', 'NONE', 'NONE', 'NONE', 'NONE', 'NONE', 'NONE'")
        
        # R objects second
        rm(sdm_map)
        rm(tmp_crop_raster, tmp_pasture_raster, tmp_urban_raster,
           tmp_crop_2020, tmp_past_2020, tmp_urb_2020)
        rm(crop_values, pasture_values, urb_values,
           crop_2020_values, past_2020_values, urb_2020_values,
           crop_intens_values, pasture_intens_values,
           crop_intens_2020_values, past_intens_2020_values,
           sdm_values)
        
        # Saving data frames for each species if they exceed a certain number of rows
        # Doing this because appending additional rows to the bottom
        # Of a really big data frame becomes increasingly slow
        # As the size of the data frame increases
        # if(nrow(save.spatial.pop.final) > 1e6) {
        # Spatial population estimates by patch
        
        write.csv(save.spatial.pop.final,
                  paste0(csv.estimate.dir,"/",i,'/',species.frame$taxon[k],"/Spatial_Population_Estimates_csv_",
                         species.frame$binomial[k],
                         '_',y,'_',t,'_',ssp,
                         migration['climate_file_save'],
                         ".csv"))
        # Non spatial population estimates
        write.csv(final.density.frame,
                  paste0(csv.estimate.dir,"/",i,'/',species.frame$taxon[k],"/Non_Spatial_Density_Estimates_",
                         species.frame$binomial[k],
                         '_',y,'_',t,'_',ssp,
                         migration['climate_file_save'],
                         ".csv"))
        # Habitat size estimates
        write.csv(final.dat.df,
                  paste0(csv.estimate.dir,"/",i,'/',species.frame$taxon[k],"/Habitat_Size_Estimates_",
                         species.frame$binomial[k],
                         '_',y,'_',t,'_',ssp,
                         migration['climate_file_save'],
                         ".csv"))
        # Modified population abundance estimates
        write.csv(final.dat.df.pop,
                  paste0(csv.estimate.dir,"/",i,'/',species.frame$taxon[k],"/Pop_Abundance_Estimates_",
                         species.frame$binomial[k],
                         '_',y,'_',t,'_',ssp,
                         migration['climate_file_save'],
                         ".csv"))
        
        # And recreating the data frames
        save.spatial.pop.final <- data.frame()
        final.patch.df <- data.frame()
        final.density.frame <- data.frame()
        final.dat.df <- data.frame()
        final.dat.df.pop <- data.frame()
        
        # Pinging a csv to indicate the species is complete
        write.csv(data.frame(species = species.frame$Species, completed = 'yes'),
                  paste0(csv.estimate.dir,"/",i,'/',species.frame$taxon[k],'/',species.frame$binomial[k],
                         '_',y,'_',t,'_',ssp,
                         migration['climate_file_save'],
                         "_completed.csv"))
      } # End loop for thresholds
      
      rm(aoh_map_file, sdm_map_file)
    } # End loop for if else statement that checks whether species need to be looped over
} # End function that loops over species


# Reimporting SNAP Raster
# Used to limit crop and pasture extent, to be safe
Snap_Raster <- raster(paste0(getwd(),'/Global Mollweide Maps/MollweideCountryID_1.5km.tif'))

# Getting list of species endemic to biodiversity hotspots
species_by_biodiv_hotspot <- read.csv(paste0(getwd(),'/Analyses/species_by_biodiv_hotspot_esh_maps.csv'))

  # endemic to hotspots
  species_by_biodiv_hotspot <-
    species_by_biodiv_hotspot %>%
    distinct() %>%
    filter(!(ecoregion_id %in% -99)) %>%
    #dplyr::group_by(species,taxa) %>%
    #dplyr::summarise(prop_in_hotspot = sum(prop_cells)) %>%
    #distinct() %>%
    #filter(prop_in_hotspot >=1) %>%
    # filter(ecoregion_id %in% c(1,8,9)) %>%
    dplyr::group_by(species,taxa) %>%
    dplyr::summarise(prop_cells = sum(prop_cells,na.rm=TRUE)) %>%
    filter(prop_cells %in% 1)



# Start loop
for(i in scenarios.list) { # Loop through scenarios
  
  # Checking progress
  # cat(paste("Processing",i,'\n'))
  
  # Creating folder for the scenario
  dir.create(paste0(raster.write.dir,'/',i))
  dir.create(paste0(patch.estimate.dir,'/',i))
  dir.create(paste0(csv.estimate.dir,'/',i))
  
  # List of directories for taxa
  raster.write.scen.taxa.dir <- paste(raster.write.dir,i, unique(species.frame$taxon), sep = '/')
  patch.estimate.scen.taxa.dir <- paste(patch.estimate.dir,i, unique(species.frame$taxon), sep = '/')
  csv.estimate.scen.taxa.dir <- paste(csv.estimate.dir,i, unique(species.frame$taxon), sep = '/')
  
  # Creating directories
  lapply(raster.write.scen.taxa.dir, dir.create)
  lapply(patch.estimate.scen.taxa.dir, dir.create)
  lapply(csv.estimate.scen.taxa.dir, dir.create)
  
  ###
  # Creating empty data frame containing MVP patches, size of the patches ----
  # These will be appended to later
  final.patch.df <- data.frame() # Total habitat area by patch
  final.dat.df <- data.frame()  # Total habitat area across all patches
  final.density.frame <- data.frame() # Estimates of pop density
  save.spatial.pop.final <- data.frame() # Spatial estimates of pop density
  
  # Looping through ssps
  # Using SSP2 for BAU
  # With SSP 3 for sensitivity for a subset of species
  for(ssp in ssp.list[1]) {
    
    # Loop through years
    # for(y in years.list[years.list > 2005 & years.list < 2055]) { # Only 2010 to 2050. Have climate maps to 2090
    for(y in c(2020,2050)) { # For testing
      # Checking progress
      cat(paste("Processing",i,y,'\n'))
      
      # Importing raster maps in time period of forecast
      if(y %in% '2020') { # If 2020
        crop_raster <- raster(list.files(ag_exp_tifs, full.names = TRUE, pattern = 'crop.*mean') %>% .[grep(paste0(2020,'.tif'),.)])
        pasture_raster <- raster(list.files(ag_exp_tifs, full.names = TRUE, pattern = 'past.*mean') %>% .[grep(paste0(2020,'.tif'),.)])
        urban_raster = raster(list.files(urb_tifs, full.names = TRUE) %>% .[grep(2020,.)])
        intens_raster <- raster(list.files(ag_int_tifs, full.names = TRUE) %>% .[grep(paste0(2020,'.tif'),.)])
       

	crop_raster[is.na(Snap_Raster)] <- NA
	pasture_raster[is.na(Snap_Raster)] <- NA

        # No intensity if no cropland
        intens_raster[crop_raster %in% 0] <- 0
        intens_raster[is.na(crop_raster)] <- NA
        
      } else { # If not 2020
        crop_raster <- raster(list.files(ag_exp_tifs, full.names = TRUE, pattern = 'crop.*mean') %>% .[grep(paste0(y,'.tif'),.)])
        pasture_raster <- raster(list.files(ag_exp_tifs, full.names = TRUE, pattern = 'past.*mean') %>% .[grep(paste0(y,'.tif'),.)])
        urban_raster = raster(list.files(urb_tifs, full.names = TRUE) %>% .[grep(y,.)])
        intens_raster <- raster(list.files(ag_int_tifs, full.names = TRUE) %>% .[grep(paste0(y,'.tif'),.)])
        

	crop_raster[is.na(Snap_Raster)] <- NA
        pasture_raster[is.na(Snap_Raster)] <- NA

        # No crop intensity if no crop
        intens_raster[crop_raster %in% 0] <- 0
        intens_raster[is.na(crop_raster)] <- NA
        
      } # End import raster maps
      
      # Importing 2020 maps
      crop_raster_2020 <- raster(list.files(ag_exp_tifs, full.names = TRUE, pattern = 'crop.*mean') %>% .[grep(paste0(2020,'.tif'),.)])
      pasture_raster_2020 <- raster(list.files(ag_exp_tifs, full.names = TRUE, pattern = 'past.*mean') %>% .[grep(paste0(2020,'.tif'),.)])
      urban_raster_2020 = raster(list.files(urb_tifs, full.names = TRUE) %>% .[grep(2020,.)])
      intens_raster_2020 <- raster(list.files(ag_int_tifs, full.names = TRUE) %>% .[grep(paste0(2020,'.tif'),.)])
      
      # No intensity if no cropland
      crop_raster_2020[is.na(Snap_Raster)] <- NA
      pasture_raster_2020[is.na(Snap_Raster)] <- NA
      
      intens_raster_2020[crop_raster_2020 %in% 0] <- 0
      intens_raster_2020[is.na(crop_raster_2020)] <- NA
      
      # Importing climate rasters
      # These are used to estimate the population density of different species
      # Need two variables:
      # Precipitation in warmest quarter
      # And variation of precipitation throughout the year
      
      # If year is provided in the SSPs
      if(y %in% ssp.years$mean_year) {
        pwarmest_raster <- 
          list.files(paste0(getwd() ,'/CMIP6_Climate_Data/',ssp),
                     pattern = paste0(ssp.years$year_1[ssp.years$mean_year %in% y],
                                      '.*',
                                      ssp.years$year_2[ssp.years$mean_year %in% y]),
                     full.names = TRUE) %>%
          paste0(.,'/Managed_Rasters') %>%
          list.files(path = ., pattern = 'Precip.*warm.*Mollweide', full.names = TRUE) %>%
          raster()
        
        pvar_raster <- 
          list.files(paste0(getwd() ,'/CMIP6_Climate_Data/',ssp),
                     pattern = paste0(ssp.years$year_1[ssp.years$mean_year %in% y],
                                      '.*',
                                      ssp.years$year_2[ssp.years$mean_year %in% y]),
                     full.names = TRUE) %>%
          paste0(.,'/Managed_Rasters') %>%
          list.files(path = ., pattern = 'Var.*precip.*Mollweide', full.names = TRUE) %>%
          raster()
        
        # And for 2020 data
        pwarmest_raster_2020 <-
          raster(paste0(getwd(),
                        '/CMIP6_Climate_Data/',
                        ssp,
                        '/Interpolated_Rasters/pwarmest_raster_2020',
                        '.tif'))
        
        # Import pvar raster
        pvar_raster_2020 <- 
          raster(paste0(getwd() ,
                        '/CMIP6_Climate_Data/',
                        ssp,
                        '/Interpolated_Rasters/pvar_raster_2020',
                        '.tif'))
        
      } else { # If year not provided by the SSPs
        # If yes, import pwarmest raster
        pwarmest_raster <-
          raster(paste0(getwd(),
                        '/CMIP6_Climate_Data/',
                        ssp,
                        '/Interpolated_Rasters/pwarmest_raster_',
                        y,
                        '.tif'))
        
        # Import pvar raster
        pvar_raster <- 
          raster(paste0(getwd() ,
                        '/CMIP6_Climate_Data/',
                        ssp,
                        '/Interpolated_Rasters/pvar_raster_',
                        y,
                        '.tif'))
        
        # And for 2020 data
        pwarmest_raster_2020 <-
          raster(paste0(getwd(),
                        '/CMIP6_Climate_Data/',
                        ssp,
                        '/Interpolated_Rasters/pwarmest_raster_2020',
                        '.tif'))
        
        # Import pvar raster
        pvar_raster_2020 <- 
          raster(paste0(getwd() ,
                        '/CMIP6_Climate_Data/',
                        ssp,
                        '/Interpolated_Rasters/pvar_raster_2020',
                        '.tif'))
      } # End climate raster imports
      
      
      # Which species have been completed in this year
      species_completed_year <-
        species_completed %>%
        gsub('.*CSV_File_Outputs/EAT_Lancet/','',.) %>%
        #gsub(migration['species_completed_year.search'],'',.) %>%
        .[grep(y,.)] %>%
        .[grep(ssp,.)] %>%
        .[grep('specsens',.)] %>%
        gsub('_[0-9]{4,4}','',.) %>%
        gsub('_specsens.*','',.)
      
      # Filtering data frame based on species that have not been completed for the prevalence threshold - specsens is used for sensitivity
      species.frame <-
        species.frame.full %>%
        # filter(binomial %in% coo_data$species) %>% # Uncomment this if you are running for limited set of species
	# filter(!(binomial %in% ssa_species_have))
        filter(!(Species %in% species_completed_year)) %>%
	 filter(Species %in% paste0(species_by_biodiv_hotspot$taxa,'/',species_by_biodiv_hotspot$species)) 
        #filter(grepl('Brachycephalus',Species))
        # filter(grepl('Leopardus_wiedii',Species))

species.frame <-
        species.frame %>%
        sample_n(nrow(.))

cat('Number Species Completed in',y,':',length(species_completed_year),'\n')
      cat('Number Species Remaining in',y,':',nrow(species.frame),'\n')
      
# For testing
      # species.frame <- species.frame.full[c(556,14000,15151,25000),]
      
      # Loop through species
      if(nrow(species.frame) > 0) { # If statement seeing if any species need to be analysed
        #for(i in 1:nrow(species.frame)) {species_loop_function(i)}
	      #plan(multicore, workers = 5)
# set.seed(0xBEEF)
#future_lapply(1:nrow(species.frame), species_loop_function)
#plan(sequential)
	      trycatch_fun <-
                    function(ss) {
                    tryCatch(species_loop_function(ss),
                    error = function(e) {cat('Error:',ss,'\n')})
                    }

            mclapply(1:nrow(species.frame),trycatch_fun,mc.cores=1)
	      #mclapply(1:nrow(species.frame),species_loop_function,mc.cores = 1) # End species list
      } # End if statement checking to see if any species need to be analysed...
    } # End years loop
  } # End loop through SSPs
} # End scenario loop








