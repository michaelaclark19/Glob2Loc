#!/usr/bin/env Rscript

# Read me ----
# SDM modeling. Written by MC

# Layout of data frame
# col 1 = longitude
# Col 2 = latitude
# Need presence/absence data
# ESH = presence
# Lack of ESH = absence
# Next columns = climate variables
# Conversely, can also use a raster stack for bioclim, domain, and maxent

###
#

###
# Model methods
# Bioclim, domain, and Maxent work on raster stacks
# GLM and GBMrandom forest need data frames

###
# Bioclim uses only presence values
# pred_nf <- dropLayer(predictors,'biome')
# group <- kfold(bradypus, 5)
# pres_train <- bradypus[group != 1, ]
# pres_test <- bradypus[group == 1, ]
# bc <- bioclim(pred_nf, pres_train

###
# Java memory allocation
# options(java.parameters = "-Xmx30000m")

# Setting working directory
#setwd("/Users/maclark/Desktop/Multiple Stresses of Biodiversity")
setwd("/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity")
###
# Libraries
library(rlang)
library(dismo)
# library(rJava)
library(randomForest)
library(plyr)
library(dplyr)
library(purrr)
library(rgdal)
library(raster)
library(gdalUtils)
library(doMC)
#library(SpaDES)
library(SpaDES.tools)
library(countrycode)
library(readr)
library(maxnet)
library(R.utils)
library(stringr)
library(future)
library(future.apply)


###
# Functions
source(paste0(getwd(),"/Scripts/SDM Scripts/0.0_SDM_Functions_2024.07.13.R"))
# source(paste0(getwd(),"/Scripts/SDM Scripts/0.0_SDM_Functions_3Dec2021.R"))
n_cores <- 35
registerDoMC(cores=35)

# Creating climate maps for SDM files -----
# Only need to do this once
# So making an if exception so that this loop only runs the first time the script is run
if('GDD_5C_Mollweide_squared.tif' %in% list.files(paste0(getwd(),'/CMIP6_Climate_Data/SSP5-8.5/2081-2100/Managed_Rasters'))) { # Checking to see if a managed climate raster exists - if it does, do not need to do the loop
  # Do nothing
} else {
  # Make the climate maps
  
  # List of directories
  dirs.loop <-
    list.dirs(paste0(getwd(),'/CMIP6_Climate_Data')) %>% # All directories
    .[grepl('Managed',.)] # Only need directories containing managed data
  
  # Function to loop through files that need to be managed
  # This allows everything to be run in parallel
  raster.manage.function <- 
    function(d) {
      # Getting list of files - just need to square all of these rasters
      files.loop <- list.files(d, pattern = 'Mollweide',full.names = TRUE)
      # Looping through these files
      for(f in files.loop) {
        # Importing raster
        tmp.raster <- raster(f)
        # Squaring raster
        tmp.raster <- tmp.raster^2
        # Writing raster
        writeRaster(tmp.raster, gsub('.tif','_squared.tif',f))
        # Clearing memory
        rm(tmp.raster)
      } # End loop for files
    } # End loop for function
  
  # Running in parallel
  mclapply(rev(dirs.loop), raster.manage.function, mc.cores = n_cores)
  # Checking to see if it works
  lapply(dirs.loop[1],raster.manage.function)
} # End loop to make managed climate data

###
# Making directories to save the ESH_RCPs
if('ESH_RCPs' %in% list.files(getwd())) {
  # Nothing
} else {
  # List of directories
  dirs.loop <-
    list.dirs(paste0(getwd(),'/CMIP6_Climate_Data')) %>% # All directories
    .[grepl('Managed',.)]
  # Make directory
  dir.create(paste0(getwd(),'/ESH_RCPs'))
  dir.create(paste0(getwd(),'/ESH_RCPs/Mod_Accuracy'))
  
  # Getting list of directories by RCP
  dirs.to.create <- 
    gsub('.*Climate_Data/','',dirs.loop) %>% 
    gsub('/[0-9]{4,4}.*','',.) %>%
    unique()
  
  # And creating directories for each taxon by RCP
  for(i in dirs.to.create) {
    dir.create(paste0(getwd(),'/ESH_RCPs/',i)) # Directory for SDMs
    if(i != 'Mod_Accuracy') {
      sub.dirs <- paste0(i, c('/Amphibians','/Mammals','/Birds','/Reptiles')) # Taxons in the analysis
      for(j in sub.dirs) {
        dir.create(paste0(getwd(),'/ESH_RCPs/',j))
      }
    }
  }
}

###
# Making directories to save raw SDM models
dir.create(paste0(getwd(),'/ESH_RCPs/Raw/'))
dir.create(paste0(getwd(),'/ESH_RCPs/Raw/Amphibians'))
dir.create(paste0(getwd(),'/ESH_RCPs/Raw/Mammals'))
dir.create(paste0(getwd(),'/ESH_RCPs/Raw/Birds'))

# Update on progress of script
cat('Done Managing Climate Maps\n Loading Biome Maps')

### Importing Other data sets ---
# Ecoregions and biomes
# List of ESH files
# Country of occurrence information

# Using this to limit where species can exist in the future
biome.map <- raster(paste0(getwd(),"/Ecoregions_Feb2023/biomes_raster.tif"))
realm.map <- raster(paste0(getwd(),"/Ecoregions_Feb2023/realms_raster.tif"))
tnc.eco <- raster(paste0(getwd(),"/Ecoregions_Feb2023/biomes_realms_raster.tif"))

# Ecoregions map - used to randomly select points near species' habitat
ecoregions.map <- raster(paste0(getwd(),'/Global Mollweide Maps/TNC Ecoregions Mollweide 1.5km.tif'))
###
# Getting intersection between biome and realm
# This is used to limit SDM maps




# Loading ESH maps -----
esh.reptiles <-
  list.files(path = paste0(getwd(),"/ESH_Tifs_12Oct/Reptiles"),
             full.names = TRUE)

esh.all <- esh.reptiles
taxon <- rep("Reptiles", length(esh.all))
# esh.all <- esh.amps
# taxon <- c(rep("Amphibians", length(esh.amps)))

# Country of occurrence data 
# Used to drop species from analysis if we don't have iucn information on them
coo.dat <- read.csv(paste0(getwd(),"/COO_Data/NewCountryOfOccurrenceData.csv"),stringsAsFactors = FALSE) %>%
  mutate(continent = countrycode(ISO3, origin = 'iso3c',destination = 'continent')) %>%
  mutate(region = countrycode(ISO3,origin = 'iso3c', destination = 'region'))

# Next check
cat('\n Managing Species Rasters')

# Managing species raster information
# Need to make sure we're only importing rasters, and not web files
esh.tot <-
  data.frame(ESH_tif = esh.all,
             taxa = taxon, 
             stringsAsFactors = FALSE) %>%
  filter(!grepl(".xml", ESH_tif)) %>%
  # Some of the file names are a bit messed
  mutate(species_name = gsub("\\.tif$", "", ESH_tif)) %>%
  mutate(species_name = gsub("(^.*/)(.*)(_HabElev.*$)", "\\2", species_name)) %>%
  mutate(species_name = gsub("_P.*$", "", species_name)) %>%
  mutate(species_name = gsub(".*/", "", species_name))

# Which species do we already have?
# Dropping species we already have from the analysis
# Don't need to rerun the SDM models if we already have them
already.have <-
  c(list.files(path = paste0(getwd(),"/ESH_RCPs/Mod_Accuracy/"),full.names = FALSE)) %>%
  gsub(".*Accuracy_","",.) %>%
  gsub(".csv","",.) %>%
  gsub("_Reptiles",'',.)

# And filtering out species we already have
if(length(already.have) > 0) {
  esh.tot <- 
    esh.tot %>%
    # filter(species_name %in% coo.dat$species) %>%
    filter(!(species_name %in% already.have))
}

### ---
# Drop species with fewer than 100 cells
# Can't fit SDMs for these species because of lack of data
cat('\n Counting Number of Cells in Habitat Range Maps')

# Function to count number of cells in each raster
cell_count <- function(x){
  ncell(raster(x))
}

# Running the function takes a while 
if('Number_of_cells_by_species_reptiles.csv' %in% list.files(paste0(getwd(),'/Other Data Inputs/Cell_Count'))) {
  esh.tot <-
    read.csv(paste0(getwd(),'/Other Data Inputs/Cell_Count/Number_of_cells_by_species_reptiles.csv'))  %>%
    filter(!(species_name %in% already.have))	       
  
} else {
  # Counting number of cells and saving file
  esh.tot <- #esh.tot %>%
    esh.tot %>% # For testing and troubleshooting
    mutate(no_cells = map_dbl(ESH_tif, cell_count))# %>%
  # filter(no_cells >= 100 & no_cells <= 50000000) %>%
  #mutate(new_extent = map2(ESH_tif, taxa, extent_fun))
  
  dir.create(paste0(getwd(),'/Other Data Inputs/Cell_Count'))
  write.csv(esh.tot,paste0(getwd(),'/Other Data Inputs/Cell_Count/Number_of_cells_by_species_reptiles.csv'),row.names=FALSE)
}

### ---
# Managing ecoregion raster - this is really biome by realm combinations - calling ecoregion here to integrate with rest of script
# Converting raster to a matrix
# Doing this to prevent future migration beyond the extent of the current ecoregion the species are found in
cat('\n Managing Ecoregion Matrix')
tnc.eco <- raster(paste0(getwd(),"/Ecoregions_Feb2023/biomes_realms_raster.tif"))
# tnc.eco.matrix <- matrix(tnc.eco,nrow = tnc.eco@nrows,byrow = TRUE)

# extents for biome_realm intersections
realm.biome.extents <- 
	read.csv(paste0(getwd(),'/Ecoregions_Feb2023/biomes_realms_extents.csv')) %>%
	dplyr::rename(realm_biome_id = eco.keep)


### ---
# Creating temporary files used to save individual SDM tiles
# Need to do this because otherwise memory space runs out and the script stops...
# These are used to save intermediate raster outputs
dir.create(paste0(getwd(),'/TMP_SDM_FILE_1'))
dir.create(paste0(getwd(),'/TMP_SDM_FILE_1_MOSAIC'))

# Saving path names so I don't need to copy and paste a bunch later
path.sdms.write.tmp <- paste0(getwd(),'/TMP_SDM_FILE_1') %>% gsub('ouce-glob2loc','pubh-glob2loc',.)
path.sdms.write.mosaiced.tmp <- paste0(getwd(),'/TMP_SDM_FILE_1_MOSAIC') %>% gsub('ouce-glob2loc','pubh-glob2loc',.)

# List of ssps
ssp.list <- list.files(paste0(getwd(),'/CMIP6_Climate_Data'))


### ---
# Now making a wrapper to run the SDM forecast
# The function uses a species name as an input
# Then for each species, it projects SDMs in the different RCPs/SSPs 
# And within each RCP/SSP, it projects SDMs in different time periods
# It then interpolates SDMs beetween the time periods (RCPs/SSPs provide snapshots at certain periods, we're deriving SDMs at 5-year interval for each species)
cat('\n Building SDM Wrap Function')
# And here's the function
sdm_wrap <- function(k) {
  
	rm(auc.df) 
      
	# Clearing memory space
  gc()
  removeTmpFiles()
  # Setting seed for reproducability
  set.seed(19)
  
  # Progress update - sometimes useful to see these outputs
  cat(df$species_name[k])
  # print(df$species_name[k])
  
  # Importing esh raster for the species
  # This is the map of where the species is found
  esh.raster <- raster(df$ESH_tif[k])
  
  # Empty object - used for logic check
  test.mods <- c()

  # Importing realm_biome raster and extent

 
 cat('Afer test.mods') 
  # Changing 0s to NAs
  # Need to do this for some ESH rasters
  # NAs indicate outside of habitat range, 1s indicate in habitat range
  if(minValue(esh.raster) %in% 0) {
    esh.raster[esh.raster %in% 0] <- NA
  }
  
  # Getting total number of cells in which the species exists
  # Using a minimum threshold of 10 cells for the modelling
  # Using such a small amount of data risks overfitting
  # But there are other checks for overfitting throughout the script
  # And There will be sensitivity analyses later on to check this assumptions
  esh.check <- sum(getValues(esh.raster),na.rm=TRUE)
  
  cat(df$species_name[k],' Checking number of cells in habitat\n')
  if(esh.check >= 25) { # If species' habitat contains more than 50 cells - run the models
    
    # Getting new habitat extent
    # This is based on overlaps between the esh raster and eco regions
    # Species' extent is limited to the ecoregions in which a species currently exists
    # Possibly revisit assumption for long-term dynamics - could be better to use realm + biome combinations here, and then restrict species to ecoregions later
    tmp.eco <- raster::crop(tnc.eco,esh.raster) # Cropping ecoregion map based on species' range 
    eco.keep <- unique(getValues(tmp.eco)[getValues(esh.raster) %in% 1]) %>% .[!is.na(.)] # Getting list of eco regions in which species exists
    
    # The eco region map only contains terrestrial ecoregions
    # If length of eco.keep is 0, it means the species is marine/aquatic
    # So don't forecast SDMs of these species
    if(length(eco.keep) > 0) { # If species is at least somewhat terrestrial...
      # Filtering df to only include ecoregions in the species' habitat range
      new.extent.df <-
        realm.biome.extents %>% 
	filter(realm_biome_id %in% eco.keep) %>%
	filter(realm_biome_id != 1)
      # Creating new extent
      new.extent <- c(min(new.extent.df$x_min),
                      max(new.extent.df$x_max),
                      min(new.extent.df$y_min),
                      max(new.extent.df$y_max))
      
      # Extending esh raster to allow for migration beyond current habitat range
      esh.raster <- raster::extend(esh.raster, 
                                   extent(tnc.eco, new.extent[3],new.extent[4],new.extent[1],new.extent[2]), 
                                   value = NA)
      
      
      
      # Cropping the ecoregion map based on the new extent
      tmp.eco <- raster::crop(tnc.eco,esh.raster)
      tmp.realm <- raster::crop(realm.map, esh.raster)
      tmp.biome <- raster::crop(biome.map, esh.raster)
      tmp.ecoregions <- raster::crop(ecoregions.map, esh.raster)
      # esh.raster <- raster::extend(x = esh.raster, 
      #                              y = extent(tmp.eco) + c(-22500,22500,-22500,22500), # With an additional buffer for sanity
      #                              value = NA)

      
    } # End check on whether species has some habitat that overlaps with land
    
    # 
    if(esh.raster@ncols * esh.raster@nrows < 25 | length(eco.keep) %in% 0) {
      # Saving accuracy as NAs to keep track of which species have been asssessed
      # Otherwiese, do no modelling
      auc.df <-
        data.frame(species = df$species_name[k],
                   taxa = df$taxa[k],
                   bioclim = 'total raster includes <25 cells',
                   maxent = 'total raster includes <25 cells',
                   rf = 'total raster includes <25 cells',
                   glm.gaussian = 'total raster includes <25 cells')
      
      write.csv(auc.df,
                paste0(getwd(),"/ESH_RCPs/Mod_Accuracy/Model_Accuracy_",
                       df$species_name[k], "_", df$taxa[k],
                       '.csv'),
                row.names = FALSE)
      
      test.thresholds <-
        data.frame(kappa = 'total raster includes <25 cells',
                   spec_sens = 'total raster includes <25 cells',
                   no_omission = 'total raster includes <25 cells',
                   prevalence = 'total raster includes <25 cells',
                   equal_sens_spec = 'total raster includes <25 cells',
                   sensitivity = 'total raster includes <25 cells',
                   model_name = c('bioclim','maxent','random.forest','glm.gaussian'))
      
      write.csv(test.thresholds,
                paste0(getwd(),"/ESH_RCPs/Mod_Accuracy/Model_Thresholds_",
                       df$species_name[k], "_", df$taxa[k],
                       '.csv'),
                row.names = FALSE)
      cat(df$species_name[k]) 
    } else { # Do a bunch of SDM modelling
      cat(df$species_name[k], ' Importing 2005 Climate Data \n')
      
      # Structure of this is as follows
      # Loop through RCPs/SSPs
      # For each RCP/SSP, predict SDM at t0 and at t1
      # Interpoalte between t0 and t1 to get SDMs at 5-year intervals
      # Save rasters
      # Remove to
      # Set t1 to t0, then start the next iteration of the loop
      
      # Doing historic maps first
      # Last time period for historic maps is 2010 (or so)
      # This output is saved (not deleted)
      # Because it is used as the first t0 for the SDM projections in the future SSPs/RCPs
      
      # There are various quality checks throughout
      # E.g. if models fit poorly, don't forecast
      # Limit species' habitat to remain within same ecoregions
      
      # Fitting models using 2010 data
      # If models fit well, then run predictions
      # If not, then do not run predictions
      # Developing the model
      # Getting climate data
      # Order to add climate data
      order.import <- c('GDD_10C_Mollweide','Minimum_Temperature_Mollweide','Annual_Precipitation_Mollweide','Water_Balance_Mollweide',
                        'GDD_10C_Moll.*squared','Minimum_Temperature_Moll.*squared','Annual_Precipitation_Moll.*squared','Water_Balance_Moll.*squared')
      # Making empty raster stack
      tmp.stack.bc.xm <- raster::stack()
      # Looping to add climate information to raster stack
      for(o in order.import) {
        tmp.stack.bc.xm <-
          raster::stack(tmp.stack.bc.xm,
                        raster(list.files(paste0(getwd(),'/CMIP6_Climate_Data/Historic/1995-2014/Managed_Rasters'),
                                          full.names = TRUE,
                                          pattern = paste0(o,'.tif'))) %>%
                          raster::crop(.,esh.raster))
      } # End loop to import climate data
      
      # Updating names of raster stack
      # Doing this to work with rest of the script
      names(tmp.stack.bc.xm) <- 
        c("GDD_2010",'Min_Temp_2010','Precip_2010','Water_2010',
          'GDD_2010_sq','Min_Temp_2010_sq','Precip_2010_sq','Water_2010_sq')
      
      # Creating data frame - this is used to create the sdm models
      tmp.df.bc.xm <-
        data.frame(cell_num = 1:(esh.raster@ncols * esh.raster@nrows),
                   esh = getValues(esh.raster),
                   GDD = getValues(tmp.stack.bc.xm$GDD_2010),
                   Precip = getValues(tmp.stack.bc.xm$Precip_2010),
                   Water_Bal = getValues(tmp.stack.bc.xm$Water_2010),
                   Min_temp = getValues(tmp.stack.bc.xm$Min_Temp_2010),
                   GDD_sq = getValues(tmp.stack.bc.xm$GDD_2010_sq),
                   Precip_sq = getValues(tmp.stack.bc.xm$Precip_2010_sq),
                   Water_Bal_sq = getValues(tmp.stack.bc.xm$Water_2010_sq),
                   Min_temp_sq = getValues(tmp.stack.bc.xm$Min_Temp_2010_sq))
      
      # Now running the models
      # This might return a warning message
      # On fitting models with only one/two outputs - this is ok, as we're dealing with presence/absence data
      
      cat(df$species_name[k],'Predicting 2005 SDMs \n')
      
      # Setting seed for reproducability
      set.seed(19)
      # T1 = Sys.time()
      # Checking for unique climate observations
      # Min threshold is 25 values
      # In the species ESH map
      # But only making this check if the species has fewer than 500 cells in its entire habitat range
      if(esh.check <= 500) {
	      if(nrow(tmp.df.bc.xm %>% 
		      dplyr::select(esh, Precip, Min_temp) %>%
		      mutate(Precip = round(Precip/10,1),
			     Min_temp = round(Min_temp,1)) %>%
		      filter(esh %in% 1) %>%
		      unique()) < 5) {
		      # Do nothing
		      test.mods <- 'limited variation in climate data'
		      test.acc <- 0

		      cat('Limited Variation in Climate Data\n')
	      } else {
	      # Tracking progress...
cat('Predicting 2005 SDMs\n')
# Creating temporary holder - this is used in logical check below to run test mods
test.mods <- 'tmp holder'
counter = 1
while(counter < 5 & length(test.mods) <= 1 & !is.null(test.mods)) {
  test.mods <- 
    tryCatch(sdm.mod.fun.eco(#raster.stack = tmp.stack,
      raster.stack.bc.xm = tmp.stack.bc.xm,
      #df.dat = tmp.df,
      df.dat.bc.xm = tmp.df.bc.xm,
      esh.raster = esh.raster,
      realm = tmp.realm,
      biome = tmp.biome,
      eco = tmp.eco),
      error = function(e) 'Try again')
  
  counter = counter + 1
}

if(!is.null(test.mods)) {if(test.mods[1] %in% 'Try again') {test.acc <- 0}}
# Tracking progress...
cat('Finished Predicting Models')
	      }

      } else {
   # Do the modelling
# Tracking progress...
cat('Predicting 2005 SDMs\n')
# Creating temporary holder - this is used in logical check below to run test mods
test.mods <- 'tmp holder'
counter = 1
while(counter < 5 & length(test.mods) <= 1 & !is.null(test.mods)) {
  test.mods <- 
    tryCatch(sdm.mod.fun.eco(#raster.stack = tmp.stack,
      raster.stack.bc.xm = tmp.stack.bc.xm,
      #df.dat = tmp.df,
      df.dat.bc.xm = tmp.df.bc.xm,
      esh.raster = esh.raster,
      realm = tmp.realm,
      biome = tmp.biome,
      eco = tmp.eco),
      error = function(e) 'Try again')
  
  counter = counter + 1
}

if(!is.null(test.mods)) {if(is.character(test.mods[1])) {test.acc <- 0}}
# Tracking progress...
cat('Finished Predicting Models')

      }
       

     cat('L600\n') 
      # Logic check
      # This is for species that have <50 cells that overlap with land (land as per realm/biome maps)
      # Or for species with fewer than 25 unique observations of climate data

      if(is.null(test.mods)) {
        test.mods <- 'yay' # Placeholder signifying how to save the model accuracy info
        test.acc <- 0
      } else if(length(test.mods) %in% 1) { 
              # Do nothing
      } else {
        cat('L632')
	      # Getting accuracy
        test.acc <- 
          c(test.mods$AUC[['bioclim']]@auc,
            test.mods$AUC[['maxent']]@auc,
            test.mods$AUC[['random.forest']]@auc,
            test.mods$AUC[['glm.gaussian']]@auc)
      cat('\n Creating thresholds, L592')
        
        test.thresholds <-
          do.call(rbind,
                  lapply(1:length(test.mods$Thresholds),
                         function(i) {
                           test.mods$Thresholds[[i]] <-
                             test.mods$Thresholds[[i]] %>%
                             mutate(model_name = names(test.mods$Thresholds)[i])
                         }))
      }

     cat('After Test Accuracy, L630\n')
     #test.acc
      # Were there models that fit well?
      # If yes, loop through SSPs
      # If not, then don't loop
      if(!length(test.mods) > 1 | max(test.acc[test.acc < 1]) < 0.75) { # Poorly fitting models
        cat('Checking Poorly Fitting Models, L635\n')
	      
	      # And writing test preds accuracy
        if(test.mods[1] %in% 'yay') {
		cat('Saving Poorly Fitting Models, L640\n')
          # This is for species with fewer than 10 cells that overlap with land
          auc.df <-
            data.frame(species = df$species_name[k],
                       taxa = df$taxa[k],
                       bioclim = 'not enough cells on land to fit SDMs',
                       maxent = 'not enough cells on land to fit SDMs',
                       rf = 'not enough cells on land to fit SDMs',
                       glm.gaussian = 'not enough cells on land to fit SDMs')
          
          test.thresholds <-
            data.frame(kappa = 'not enough cells on land to fit SDMs',
                       spec_sens = 'not enough cells on land to fit SDMs',
                       no_omission = 'not enough cells on land to fit SDMs',
                       prevalence = 'not enough cells on land to fit SDMs',
                       equal_sens_spec = 'not enough cells on land to fit SDMs',
                       sensitivity = 'not enough cells on land to fit SDMs',
                       model_name = c('bioclim','maxent','random.forest','glm.gaussian'))
          
          write.csv(auc.df,
                    paste0(getwd(),"/ESH_RCPs/Mod_Accuracy/Model_Accuracy_",
                           df$species_name[k], "_", df$taxa[k],
                           '.csv'),
                    row.names = FALSE)
          
          write.csv(test.thresholds,
                    paste0(getwd(),"/ESH_RCPs/Mod_Accuracy/Model_Thresholds_",
                           df$species_name[k], "_", df$taxa[k],
                           '.csv'),
                    row.names = FALSE)
          
          
          
          rm(auc.df)
        } else if(test.mods[1] %in% 'limited variation in climate data') {
	cat('Saving Poorly Fitting Models, L670\n')
	# This is for species with limited variation in climate data
          auc.df <-
            data.frame(species = df$species_name[k],
                       taxa = df$taxa[k],
                       bioclim = 'limited variation in climate data',
                       maxent = 'limited variation in climate data',
                       rf = 'limited variation in climate data',
                       glm.gaussian = 'limited variation in climate data')
          
          test.thresholds <-
            data.frame(kappa = 'limited variation in climate data',
                       spec_sens = 'limited variation in climate data',
                       no_omission = 'limited variation in climate data',
                       prevalence = 'limited variation in climate data',
                       equal_sens_spec = 'limited variation in climate data',
                       sensitivity = 'limited variation in climate data',
                       model_name = c('bioclim','maxent','random.forest','glm.gaussian'))

          write.csv(auc.df,
                    paste0(getwd(),"/ESH_RCPs/Mod_Accuracy/Model_Accuracy_",
                           df$species_name[k], "_", df$taxa[k],
                           '.csv'),
                    row.names = FALSE)
          
          write.csv(test.thresholds,
                    paste0(getwd(),"/ESH_RCPs/Mod_Accuracy/Model_Thresholds_",
                           df$species_name[k], "_", df$taxa[k],
                           '.csv'),
                    row.names = FALSE)

          rm(auc.df)
        } else if(test.mods[1] %in% 'Try again') {
        # This is for species with limited variation in climate data
          auc.df <-
            data.frame(species = df$species_name[k],
                       taxa = df$taxa[k],
                       bioclim = 'Maxent AUC error',
                       maxent = 'Maxent AUC error',
                       rf = 'Maxent AUC error',
                       glm.gaussian = 'Maxent AUC error')

          test.thresholds <-
            data.frame(kappa = 'Maxent AUC error',
                       spec_sens = 'Maxent AUC error',
                       no_omission = 'Maxent AUC error',
                       prevalence = 'Maxent AUC error',
                       equal_sens_spec = 'Maxent AUC error',
                       sensitivity = 'Maxent AUC error',
                       model_name = c('bioclim','maxent','random.forest','glm.gaussian'))

          write.csv(auc.df,
                    paste0(getwd(),"/ESH_RCPs/Mod_Accuracy/Model_Accuracy_",
                           df$species_name[k], "_", df$taxa[k],
                           '.csv'),
                    row.names = FALSE)

          write.csv(test.thresholds,
                    paste0(getwd(),"/ESH_RCPs/Mod_Accuracy/Model_Thresholds_",
                           df$species_name[k], "_", df$taxa[k],
                           '.csv'),
                    row.names = FALSE)

          rm(auc.df)
        } else { # This is for species with poorly fitting models
          cat('Saving Poorly Fitting Models, L705\n')
		
		# Saving data
          auc.df <-
            data.frame(species = df$species_name[k],
                       taxa = df$taxa[k],
                       bioclim = test.mods$AUC[['bioclim']]@auc,
                       maxent = test.mods$AUC[['maxent']]@auc,
                       rf = test.mods$AUC[['random.forest']]@auc,
                       glm.gaussian = test.mods$AUC[['glm.gaussian']]@auc)
          cat('\n Creating Thresholds L684')
          test.thresholds <-
            do.call(rbind,
                    lapply(1:length(test.mods$Thresholds),
                           function(i) {
                             test.mods$Thresholds[[i]] <-
                               test.mods$Thresholds[[i]] %>%
                               mutate(model_name = names(test.mods$Thresholds)[i])
                           }))
          
          write.csv(auc.df,
                    paste0(getwd(),"/ESH_RCPs/Mod_Accuracy/Model_Accuracy_",
                           df$species_name[k], "_", df$taxa[k],
                           '.csv'),
                    row.names = FALSE)
          
          write.csv(test.thresholds,
                    paste0(getwd(),"/ESH_RCPs/Mod_Accuracy/Model_Thresholds_",
                           df$species_name[k], "_", df$taxa[k],
                           '.csv'),
                    row.names = FALSE)
          
          rm(auc.df)
        } # End else statement to save species with poorly fitting models
        
      } else if (sum(test.acc %in% 1) == 4) { # Saving models if they are excessively accurate
      cat('Saving Poorly Fitting Models, L740\n')
	     
	      # This is for species with poorly fitting models
          auc.df <-
            data.frame(species = df$species_name[k],
                       taxa = df$taxa[k],
                       bioclim = test.mods$AUC[['bioclim']]@auc,
                       maxent = test.mods$AUC[['maxent']]@auc,
                       rf = test.mods$AUC[['random.forest']]@auc,
                       glm.gaussian = test.mods$AUC[['glm.gaussian']]@auc)
         cat('\n Creating Thresholds L718') 
          test.thresholds <-
            do.call(rbind,
                    lapply(1:length(test.mods$Thresholds),
                           function(i) {
                             test.mods$Thresholds[[i]] <-
                               test.mods$Thresholds[[i]] %>%
                               mutate(model_name = names(test.mods$Thresholds)[i])
                           }))

          write.csv(auc.df,
                    paste0(getwd(),"/ESH_RCPs/Mod_Accuracy/Model_Accuracy_",
                           df$species_name[k], "_", df$taxa[k],
                           '.csv'),
                    row.names = FALSE)
          
          write.csv(test.thresholds,
                    paste0(getwd(),"/ESH_RCPs/Mod_Accuracy/Model_Thresholds_",
                           df$species_name[k], "_", df$taxa[k],
                           '.csv'),
                    row.names = FALSE)

          rm(auc.df)
      } else { 
	      cat('Beginning Model Forecasts, L755\n')
	      # Model SSPs
        # Beginning loop through each 
        # Loop through SSPs
        test.preds.historic <- c()
        for(ssp in ssp.list) {
          gc()
		# Getting years for the SSP scenario
          year.list <- 
            list.dirs(paste0(getwd(),'/CMIP6_Climate_Data/',ssp)) %>%
            .[grepl('Managed',.)]
          # List of years
          year.list <- str_extract(year.list,'[0-9]{4,4}-[0-9]{4,4}')
          
          # Going backward through historic years 
          # Loop below is set up to go from 2010 to 2100, 
          # using 2010 as the current habitat map
          # So need to go from 2010 to 1850 in historic climate data
          # which means reversing years
          if(ssp %in% 'Historic') {
            year.list <- rev(year.list) # Reversing years
            year.list <- year.list[year.list != '1981-2010'] # Dropping year - this results in duplicate estimates to 1986-2005 (both have average year of 1995!)
          }
          
          # Looping through years
          for(year in year.list) {
            cat(ssp,year, '\n')
            
            # Importing climate data to predict models
            order.import <- c('GDD_10C_Mollweide','Minimum_Temperature_Mollweide','Annual_Precipitation_Mollweide','Water_Balance_Mollweide',
                              'GDD_10C_Mollweide_squared','Minimum_Temperature_Mollweide_squared','Annual_Precipitation_Mollweide_squared','Water_Balance_Mollweide_squared')
            # Making empty raster stack
            tmp.stack.bc.xm <- raster::stack()
            # Looping to add climate information to raster stack
            for(o in order.import) {
              tmp.stack.bc.xm <-
                raster::stack(tmp.stack.bc.xm,
                              raster(list.files(paste0(getwd(),'/CMIP6_Climate_Data/',
                                                       ssp,'/',
                                                       year,'/Managed_Rasters'),
                                                
                                                full.names = TRUE,
                                                pattern = paste0(o,'.tif'))) %>%
                                raster::crop(.,esh.raster))
            } # End loop to import climate data
            
            # Updating names of raster stack
            # Doing this to work with rest of the script
            names(tmp.stack.bc.xm) <- 
              c("GDD_2010",'Min_Temp_2010','Precip_2010','Water_2010',
                'GDD_2010_sq','Min_Temp_2010_sq','Precip_2010_sq','Water_2010_sq')
            
            # Creating data frame - this is used to create the sdm models
            tmp.df.bc.xm <-
              data.frame(cell_num = 1:(esh.raster@ncols * esh.raster@nrows),
                         esh = getValues(esh.raster),
                         GDD = getValues(tmp.stack.bc.xm$GDD_2010),
                         Precip = getValues(tmp.stack.bc.xm$Precip_2010),
                         Water_Bal = getValues(tmp.stack.bc.xm$Water_2010),
                         Min_temp = getValues(tmp.stack.bc.xm$Min_Temp_2010),
                         GDD_sq = getValues(tmp.stack.bc.xm$GDD_2010_sq),
                         Precip_sq = getValues(tmp.stack.bc.xm$Precip_2010_sq),
                         Water_Bal_sq = getValues(tmp.stack.bc.xm$Water_2010_sq),
                         Min_temp_sq = getValues(tmp.stack.bc.xm$Min_Temp_2010_sq))
            
            
            # Now predicting habitat extent
            # If statement here avoids running code if there are no models to predict from
            # e.g. if models are performing poorly - this will save some time
            test.preds <- c()
            
            if(length(test.mods) > 1) {
             cat('\n Predicting',df$species_name[k],ssp, year,'Habitat Range') 
              # raster.stack.pred = tmp.stack
              raster.stack.bc.xm.pred = tmp.stack.bc.xm
              
              # Splitting habitat range into a raster stack
              # Doing this to reduce memory space needed to run the models
              # Memory space isn't a huge issue for the glm models, but is for maxent
              
	      if(ssp %in% ssp.list[1] &
		 year %in% year.list[1]) {
	      split.current.extents <-
                split.raster.function(raster.stack = raster.stack.bc.xm.pred,
                                      esh.raster = tmp.eco)

              # Limiting split current extents to only the ecoregions in which the species exists
              # Hopefully this saves a bit of time
              tmp.eco.fun <- function(i) {if(sum(unique(getValues(crop(tmp.eco,i))) %in% eco.keep) >= 1) {return(i)}}
              split.current.extents <- lapply(split.current.extents, tmp.eco.fun) # Dropping tiles w/o eco region
              split.current.extents <- plyr::compact(split.current.extents) # And removing nulls then condensing the list
	      
	      }
              
              # Saving species names
              # used to save raster files and model outputs later on
              species_name_save = df$species_name[k]
              mod.list = test.mods
              
              # And predicting range
              test.preds <-
                sdm.pred.fun(#raster.stack.pred = raster.stack.pred,
                  raster.stack.bc.xm.pred = raster.stack.bc.xm.pred,
                  #df.dat = tmp.df,
                  df.dat.bc.xm = tmp.df.bc.xm,
                  mod.list = test.mods,
                  path.sdms.write.tmp = path.sdms.write.tmp,
                  path.sdms.write.mosaiced.tmp = path.sdms.write.mosaiced.tmp,
                  split.current.extents = split.current.extents,
                  species_name_save = df$species_name[k])
              
	      cat('\n Finished Predicted Models', ssp, year)


              # Saving outputs for interpolating
              test.preds.historic <- test.preds
              # test.preds.t0 <- test.preds
              save.raster.stack.names <- names(tmp.stack.bc.xm)
              
              
              # Saving raw outputs
              raw.preds.stack <- test.preds.historic
              
              if(!is.null(raw.preds.stack)) {
                # Updating names
               
		 names(raw.preds.stack) <- paste0(ssp,'_',year,'_',names(raw.preds.stack))
               
	      cat('Saving Raster Stack')
	      raw.preds.stack

	      # Checking names
	      paste0(getwd(),'ESH_RCPs/Raw/',df$taxa[k],'/',df$species_name[k],'_',names(raw.preds.stack),'.tif')
	      # And saving raster

	      for(ii in 1:length(names(raw.preds.stack))) {
  writeRaster(raw.preds.stack[[ii]],
              paste0(getwd(),"/ESH_RCPs/Raw/",
                     df$taxa[k], "/",
                     df$species_name[k],
                     "_",
                     names(raw.preds.stack)[ii],
                     '.tif') %>%
                gsub('ouce-glob2loc','pubh-glob2loc',.),
              overwrite = TRUE)
}


              } # End if statement to save models
              
              # And saving accuracy - only doing this for the last scenario for book keeping
              if(ssp %in% 'SSP5-8.5' &
                 year %in% '2081-2100') {


		      cat('SSP5 2100 Check')
                # Some more book keeping - making sure to update maxent if it ran into an error
		      if(test.mods$AUC[['maxent']]@auc %in% 1.5) {
			      auc.df <-
                  data.frame(species = df$species_name[k],
                             taxa = df$taxa[k],
                             bioclim = test.mods$AUC[['bioclim']]@auc,
                             maxent = 'Could not fit model - error running maxnet function',
                             rf = test.mods$AUC[['random.forest']]@auc,
                             glm.gaussian = test.mods$AUC[['glm.gaussian']]@auc)


		      test.mods$Thresholds[['maxent']][,1:ncol(test.mods$Thresholds[['maxent']])] <- 'Could not fit model - error running maxnet function'
		      } else {
		              # And writing test preds accuracy
                auc.df <-
                  data.frame(species = df$species_name[k],
                             taxa = df$taxa[k],
                             bioclim = test.mods$AUC[['bioclim']]@auc,
                             maxent = test.mods$AUC[['maxent']]@auc,
                             rf = test.mods$AUC[['random.forest']]@auc,
                             glm.gaussian = test.mods$AUC[['glm.gaussian']]@auc)
		      }
		      
                write.csv(auc.df,
                          paste0(getwd(),"/ESH_RCPs/Mod_Accuracy/Model_Accuracy_",
                                 df$species_name[k], "_", df$taxa[k],
                                 '.csv'),
                          row.names = FALSE)
                
                test.thresholds <-
                  do.call(rbind,
                          lapply(1:length(test.mods$Thresholds),
                                 function(i) {
                                   test.mods$Thresholds[[i]] <-
                                     test.mods$Thresholds[[i]] %>%
                                     mutate(model_name = names(test.mods$Thresholds)[i])
                                 }))
                
                write.csv(test.thresholds,
                          paste0(getwd(),"/ESH_RCPs/Mod_Accuracy/Model_Thresholds_",
                                 df$species_name[k], "_", df$taxa[k],
                                 '.csv'),
                          row.names = FALSE)
                
                rm(auc.df)
                rm(test.thresholds)
              }
             
              # And clearing memory space
              rm(test.preds)
              rm(tmp.df.bc.xm)
              rm(tmp.stack.bc.xm)
              rm(raw.preds.stack)
            } # End if statement to predict models
          } # End loop around years
        } # End loop around SSPs
      } # End else statement to predict SDMs (e.g. check if SDM predictions worked)
    } # End else statement to model SDMs (e.g. check if adequate cells in the species' habitat range)
  } else if(esh.check < 25 | !exists("eco.keep")) { # Saving information for species with poorly fitting models
    auc.df <-
      data.frame(species = df$species_name[k],
                 taxa = df$taxa[k],
                 bioclim = 'habitat includes < 25 cells',
                 maxent = 'habitat includes < 25 cells',
                 rf = 'habitat includes < 25 cells',
                 glm.gaussian = 'habitat includes < 25 cells')
    
    write.csv(auc.df,
              paste0(getwd(),"/ESH_RCPs/Mod_Accuracy/Model_Accuracy_",
                     df$species_name[k], "_", df$taxa[k],
                     '.csv'),
              row.names = FALSE)
    
    test.thresholds <-
      data.frame(kappa = 'habitat includes < 25 cells',
                 spec_sens = 'habitat includes < 25 cells',
                 no_omission = 'habitat includes < 25 cells',
                 prevalence = 'habitat includes < 25 cells',
                 equal_sens_spec = 'habitat includes < 25 cells',
                 sensitivity = 'habitat includes < 25 cells',
                 model_name = c('bioclim','maxent','random.forest','glm.gaussian'))
    
    write.csv(test.thresholds,
              paste0(getwd(),"/ESH_RCPs/Mod_Accuracy/Model_Thresholds_",
                     df$species_name[k], "_", df$taxa[k],
                     '.csv'),
              row.names = FALSE)
    cat(df$species_name[k])
  } # End check on extent of species habitat range
} # End sdm wrap function
   

cat('\n Starting Predictions')
set.seed(19)
df <- 
	esh.tot %>% 
	arrange(no_cells)  %>%
	# filter(no_cells >= 1000000) %>%
	# filter(no_cells < 10000000) %>%
	filter(grepl('Reptile',taxa)) #%>%
	# filter(!grepl('Phyllastrephus_albigula',species_name))



#sdm_wrap(1)
# df <- df %>% filter(grepl('Ischnoc|Brachycephalus',species_name))

nrow(df)
table(df$taxa)

chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
chunk_list <- chunk2(1:nrow(df),5)


tmp.fun <- function(i) {Sys.sleep(i)}
plan(multisession, workers = 5)
t1 = Sys.time()
future.apply::future_lapply(rep(5,5),tmp.fun)
Sys.time() - t1
future.apply::future_lapply(chunk_list[[1]], sdm_wrap)
plan(sequential)

# df$partition <- kfold(df,15)
#df <- df %>% filter(partition %in% 1)
#nrow(df)
#table(df$taxa)
#head(df$species_name)

#for(k in 1:nrow(df)) {
#  cat(k, df$species_name[k])
#	tryCatch(sdm_wrap(k), error = function(e) {'Whoopsies!'})  
#}

#lapply(1:nrow(df),sdm_wrap)
#cat('\n Starting Predictions')
# lapply(1:2, sdm_wrap)
#tmp.fun <- function(i) {Sys.sleep(i)}
#t1 = Sys.time()
#mclapply(rep(5,15),tmp.fun, mc.cores = 15)
#Sys.time() - t1
#mclapply(1:nrow(df),sdm_wrap, mc.cores = 15)

#tmp.k <- grep('Pristimantis_mars',df$species_name)
#sdm_wrap(tmp.k)

#plan(multisession, workers = n_cores)
#t1 = Sys.time()
#future.apply::future_lapply(rep(5,15),tmp.fun)
#Sys.time() - t1
#future.apply::future_lapply(1:nrow(df), sdm_wrap)
#plan(sequential)

#t1 = Sys.time()
#foreach(i = rep(5,15)) %dopar% tmp.fun(i)
#Sys.time() - t1
#foreach(i = 1:nrow(df)) %dopar% sdm_wrap(i)



