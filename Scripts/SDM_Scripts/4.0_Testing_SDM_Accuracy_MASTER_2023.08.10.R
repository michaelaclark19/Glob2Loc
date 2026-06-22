#!/usr/bin/env Rscript

#######
# Testing accuracy based on 2010 SDM outputs
# Doing this based on model thresholds of the different threshold types
# Is there one that is consistently best?
# This only does 2010 habitat maps, and restricts 2010 maps to remain within the same ecoregion in which species currently exist
# The outputs from this are not needed to run the 5.0 or 6.0 scripts
#######

###
# Libraries
library(plyr)
library(dplyr)
library(raster)
library(parallel)

# Installing functions
setwd("/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity")
source(paste0(getwd(),'/Scripts/SDM Scripts/0.0_SDM_Functions_2024.07.13.R'))

# Have multiple working directories because of memory issues
# Needed to get creative!
wd.list <- getwd()

# Function to import list of species
species.list.function <-
  function(i) {
    return(list.files(i, full.names = TRUE) %>%
             .[grepl('2014',.)] %>%
             gsub('__Historic.*','',.) %>%
             gsub('_Historic.*','',.) %>%
             unique())
  }


species.list <-
  do.call(c, 
          lapply(paste0(wd.list,'/ESH_RCPs/Raw/') %>% paste0(.,sort(rep(c('Amphibians','Mammals','Birds','Reptiles'),2))),
                 species.list.function))

#
cat('Species list completed\n')

# Checking for overlaps
# If in /data/ouce, drop from /data/pubh
data_ouce <- species.list[grep(getwd(),species.list)] %>% gsub('.*Mammals|.*Amphibians|.*Birds|.*Reptiles','',.)
data_pubh <- species.list[grep(getwd(),species.list)] %>% gsub('.*Mammals|.*Amphibians|.*Birds|.*Reptiles','',.)

if(sum(data_ouce %in% data_pubh) > 0 | sum(data_pubh %in% data_ouce) > 0) {
  # If overlaps, then filter, giving priority to ouce-glob2loc
  data_pubh_drop <- data_pubh[(data_pubh %in% data_ouce)]
  data_pubh_keep <- data_pubh[!(data_pubh %in% data_ouce)]
  
  # pubh species
  data_pubh_out <-
    do.call(c,
            lapply(data_pubh_drop,
                   function(i) {
                     return(species.list[grepl(getwd(),species.list)] %>%
                              .[grepl(i,.)])
                   }))
  
  # And updating species list
  species.list <- species.list[!(species.list %in% data_pubh_out)]
}


cat('Getting list of outputs in 2014\n')
# List of sdm outputs in 2014
sdms.list <-
  do.call(c,
          lapply(paste0(wd.list,'/ESH_RCPs/Raw/') %>% paste0(.,sort(rep(c('Amphibians','Mammals','Birds','Reptiles'),2))),
                 list.files,
                 full.names = TRUE)) %>%
  .[grepl('2014',.)]

# Making sure these are also in species list
species.list.sdms <-
  gsub('.*Mammals|.*Amphibians|.*Birds|.*Reptiles','',sdms.list) %>%
  gsub('__.*','',.) %>%
  gsub('_Historic.*','',.) %>%
  unique()

# Then checking for overlaps between pubh and ouce directories
# If in /data/ouce, drop from /data/pubh
# Checking for overlaps
# If in /data/ouce, drop from /data/pubh
data_ouce <- sdms.list[grep(getwd(),sdms.list)] %>% gsub('.*Mammals|.*Amphibians|.*Birds|.*Reptiles','',.)
data_pubh <- sdms.list[grep(getwd(),sdms.list)] %>% gsub('.*Mammals|.*Amphibians|.*Birds|.*Reptiles','',.)


cat('Checking overlaps from ouce and pubh directories\n')
if(sum(data_ouce %in% data_pubh) > 1 | sum(data_pubh %in% data_ouce) > 1) {
  # If overlaps, then filter, giving priority to ouce-glob2loc
  data_pubh_drop <- data_pubh[(data_pubh %in% data_ouce)]
  
  # pubh species
  data_pubh_out <-
    do.call(c,
            lapply(data_pubh_drop,
                   function(i) {
                     return(sdms.list[grepl(getwd(),sdms.list)] %>%
                              .[grepl(i,.)])
                   }))
  
  # And updating species list
  sdms.list <- sdms.list[!(sdms.list %in% data_pubh_out)]
}

# list of accuracy and threshold files
auc.list <-
  list.files(paste0(getwd(),'/ESH_RCPs/Mod_Accuracy'), full.names = TRUE)

# list of AOH rasters
aoh.list <-
  list.dirs(paste0(getwd(),'/ESH_Tifs_12Oct')) %>%
  list.files(.,full.names=TRUE)

# Ecoregion map
tnc.eco <- raster(paste0(getwd(),'/Global Mollweide Maps/TNC Ecoregions Mollweide 1.5km.tif'))
realm_biome <- raster(paste0(getwd(),'/Ecoregions_Feb2023/biomes_realms_raster.tif'))


# Creating new folders to save model threshold accuracy
if(!('Weighted_Threshold_Accuracy' %in% list.files(paste0(getwd(),'/ESH_RCPs')))) {
  # For model accuracy
  dir.create(paste0(getwd(),'/ESH_RCPs/Weighted_Threshold_Accuracy'))
  dir.create(paste0(getwd(),'/ESH_RCPs/Weighted_Threshold_Accuracy/Amphibians'))
  dir.create(paste0(getwd(),'/ESH_RCPs/Weighted_Threshold_Accuracy/Birds'))
  dir.create(paste0(getwd(),'/ESH_RCPs/Weighted_Threshold_Accuracy/Mammals'))
  dir.create(paste0(getwd(),'/ESH_RCPs/Weighted_Threshold_Accuracy/Reptiles'))
}

###
# Removing species we already have...
cat('List of species we have\n')
species.have <- 
  c(list.files(paste0(getwd(),'/ESH_RCPs/Weighted_Threshold_Accuracy/Amphibians'), full.names = TRUE),
    list.files(paste0(getwd(),'/ESH_RCPs/Weighted_Threshold_Accuracy/Birds'), full.names = TRUE),
    list.files(paste0(getwd(),'/ESH_RCPs/Weighted_Threshold_Accuracy/Mammals'), full.names = TRUE),
    list.files(paste0(getwd(),'/ESH_RCPs/Weighted_Threshold_Accuracy/Reptiles'), full.names = TRUE)) %>%
  gsub('.csv','',.) %>%
  gsub('.*Weighted_Threshold_Accuracy','',.) %>%
  unique() 

# Overlap between species list and species already have
species.need <- which(!(gsub('.*Raw','',species.list) %in% species.have))

# And removing species we already have
species.list <- species.list[species.need]
head(species.list)

cat('Species Have: ',length(species.have),'\n')
cat('Species Need: ',length(species.list),'\n')

cat('Amphibians Need: ', sum(grepl('Amphibian',species.list)),'\n')
cat('Birds Need: ', sum(grepl('Bird',species.list)),'\n')
cat('Mammals Need: ', sum(grepl('Mammal',species.list)),'\n')
cat('Reptiles Need: ', sum(grepl('Reptile',species.list)),'\n')



# Function that loops through species
species_wrap_fun <-
  function(i) {
    
    # Species taxa
    tmp.taxa <-
      gsub('.*/ESH_RCPs/Raw/','',i) %>%
      gsub('\\/.*','',.)
    
    # Species name
    tmp.name <-
      gsub(paste0('.*',tmp.taxa),'',i) %>%
      gsub('\\/','',.) %>%
      gsub('__His.*','',.)
    
    # Getting models
    model.files <-
      sdms.list[grepl(i,sdms.list)]
    
    if(sum(grepl('__Historic',model.files)) > 1 & sum(grepl('[a-z]_Historic',model.files)) > 1) {
      model.files <- model.files[grepl('[a-z]_Historic',model.files)]
    }
    
    
    
    
    # Stacking if needed
    models <- stack()
    for(k in model.files) {
      models <- raster::stack(models, raster(k))
    }
    
    # Names of models
    names(models) <-
      gsub('.*2014_','',names(models))
    
    # Checking
    if(length(grep(tmp.name,auc.list)) %in% 0) {
      # Do nothing
    } else {
      
      # Getting thresholds and accuracy
      # Accuracy
      auc.csv <-
        read.csv(auc.list[grepl(tmp.name,auc.list)] %>%
                   .[grepl('_Accuracy_',.)])
      
      
      # And another check
      if(length(grep(paste0('Threshold.*',tmp.name),auc.list)) > 0) {
        
        # Thresholds
        thresh.csv <-
          read.csv(auc.list[grepl(tmp.name,auc.list)] %>%
                     .[grepl('_Thresholds_',.)])
        
        # ESH raster
        esh.raster <-
          raster(aoh.list[grepl(tmp.taxa,aoh.list)] %>%
                   .[grepl(tmp.name,.)])
        
        # Getting ecoregions in which the species exists
        tmp.eco <- raster::crop(tnc.eco,esh.raster) # Cropping ecoregion map based on species' range 
        eco.keep <- unique(getValues(tmp.eco)[getValues(esh.raster) %in% 1]) %>% .[!is.na(.)] # Getting list of eco regions in which species exists
        
        # Recropping ecoregions to the extent of the SDM models
        # The SDM models already account for realm/biome combinations
        # So they don't need to be cropped again
        # This will be used on the weighted models to clip where species can/can't be found in 2005
        tmp.eco <- raster::crop(tnc.eco, models)
        
        # Repeating for realms and biomes
        realm_keep <- unique(getValues(raster::crop(realm_biome, esh.raster))[getValues(esh.raster) %in% 1]) %>% .[!is.na(.)] %>% .[.!=1]
        tmp.realm <- raster::crop(realm_biome, models)
        
        # Accuracy for weighting
        names(auc.csv)[names(auc.csv) %in% 'rf'] <- 'random.forest'
        
        cat('Before auc.weight\n')
        cat(names(auc.csv))
        cat('\n')
        cat(names(models))
        cat('\n')
        auc.weight <- as.numeric(auc.csv[,names(models)])
        
        cat('After auc.weight\n')
        ###
        # Next getting model average based on auc values
        names(auc.weight) <- names(models)
        # Getting weighted model means
        # List of model aucs
        auc <- auc.weight
        # Weighting based on accuracy
        # AUC = 0.5 means model is as good as random
        # Weighting based on square of the distance fomr auc = 0.5
        # So more weighting on more accurate models
        w <- (auc-0.5)^2
        
        
        # Weighting the models
        # mods.weighted <- weighted.mean(models, w)
        # This is slow, the below returns the same estimates, with the exception of rounding errors on the order of 1e-16
        w <- w/sum(w)
        
        # Speeding up if only one SDM fit well
        if(length(names(models)) %in% 1) {
          mods.weighted <- models
        } else {
          mods.weighted <- weight.fun.raster(w, models)
        }
        
        # Now Clipping by ecoregion
        # NAs outside the ecoregion
        # numeric values inside the ecoregion
        # mods.weighted[!(tmp.eco %in% eco.keep)] <- NA
        
        # extending esh raster
        esh.raster <- raster::extend(esh.raster, mods.weighted, value = NA)
        
        # Loop to clip by thresholds
        # Recording accuracy for book keeping
        test.acc <- data.frame()
        for(k in c('prevalence','equal_sens_spec')) { # Start loop through thresholds
          # Thresholds for weighting
          tr.weight <- as.numeric(thresh.csv[thresh.csv$model_name %in% names(models),k])
          
          # Getting weighted thresholds
          tr.weight <-
            weighted.mean(tr.weight,
                          w)
          
          
          # Testing overlap between models and ESH rasters
          # How much overlap is there?
          # Or much overlap is there not?
          test.df <-
            data.frame(esh = getValues(esh.raster),
                       sdm = getValues(mods.weighted),
                       eco = getValues(tmp.eco),
                       realm = getValues(tmp.realm)) %>%
            mutate(eco = ifelse(eco %in% eco.keep,1,0)) %>% # Getting count of ecoregions
            mutate(realm = ifelse(realm %in% realm_keep,1,0)) %>%
            mutate(esh_eco = esh,
                   esh_realm = esh) %>%
            mutate(esh_eco = ifelse(is.na(esh_eco),0,esh_eco)) %>% # Adjusting ESH to 1s (present) 0s (absent, but in ecoregion) and NAs (outside of ecoregion)
            mutate(esh_eco = ifelse(eco %in% 0, NA, esh_eco)) %>%
            mutate(esh_realm = ifelse(is.na(esh_realm),0,esh_realm)) %>% # Adjusting ESH to 1s (present) 0s (absent, but in realm_biome) and NAs (outside of realm_biome)
            mutate(esh_realm = ifelse(realm %in% 0, NA, esh_realm))
          
          
          cat(names(test.df))
          cat('\n')
          
          
          names(test.df)[2] <- 'sdm'
          
          # Converting estimates to 1s (present) 0s (absent, but in ecoregion) and NAs (outside of ecoregion)
          test.df <-
            test.df %>%
            mutate(sdm = ifelse(sdm < tr.weight, 0, 1)) %>%
            mutate(sdm_eco = sdm,
                   sdm_realm = sdm) %>%
            mutate(sdm_eco = ifelse(is.na(esh_eco), NA, sdm_eco),
                   sdm_realm = ifelse(is.na(esh_realm),NA,sdm_realm))
          
          # And testing
          test.acc <-
            rbind(test.acc,
                  data.frame(in_both = sum(test.df$esh_eco %in% 1 & test.df$sdm_eco %in% 1), # How many are in both SDM and ESH?
                             esh_not_sdm = sum(test.df$esh_eco %in% 1 & test.df$sdm_eco %in% 0), # How many are in ESH but not SDM?
                             sdm_not_esh = sum(test.df$esh_eco %in% 0 & test.df$sdm_eco %in% 1), # How many are in SDM but not ESH?
                             out_both = sum(test.df$esh_eco %in% 0 & test.df$sdm_eco %in% 0), # How many are in neither?
                             threshold = k, # Which threshold
                             region = 'ecoregion'), # which geographic region
                  data.frame(in_both = sum(test.df$esh_realm %in% 1 & test.df$sdm_realm %in% 1), # How many are in both SDM and ESH?
                             esh_not_sdm = sum(test.df$esh_realm %in% 1 & test.df$sdm_realm %in% 0), # How many are in ESH but not SDM?
                             sdm_not_esh = sum(test.df$esh_realm %in% 0 & test.df$sdm_realm %in% 1), # How many are in SDM but not ESH?
                             out_both = sum(test.df$esh_realm %in% 0 & test.df$sdm_realm %in% 0), # How many are in neither?
                             threshold = k, # Which threshold
                             region = 'realm_biome')) # which geographic region
          
          # Making weighted raster
          sdm.weighted.save <- raster(matrix(test.df$sdm_eco, byrow = TRUE, nrow = mods.weighted@nrows))
          crs(sdm.weighted.save) <- crs(mods.weighted)
          extent(sdm.weighted.save) <- extent(mods.weighted)
          
          # And saving
          #          writeRaster(sdm.weighted.save,
          #                      paste0(getwd(),'/ESH_RCPs/Historic/',tmp.taxa,'/',tmp.name,'2010_threshold_',k,'.tif') %>%
          #			      gsub('ouce-glob2loc','pubh-glob2loc',.),
          #                      overwrite = TRUE)
          
        } # End loop through types of thresholds
        
        # Saving data frame with model accuracy
        test.acc <-
          test.acc %>%
          mutate(taxa = tmp.taxa,
                 species = tmp.name) %>%
          mutate(tot_accuracy = (in_both + out_both) / (in_both + out_both + esh_not_sdm + sdm_not_esh))
        
        # Saving test.acc
        write.csv(test.acc,
                  paste0(getwd(),'/ESH_RCPs/Weighted_Threshold_Accuracy/',tmp.taxa,'/',tmp.name,'.csv'),
                  row.names = FALSE)
        
      } # End check for thresholds
    } # End else statement checking whether accuracy data for species SDM exists
  } # End function

length(species.list)

#species.list <- species.list[grepl('Reptile',species.list)]

cat('Running loop\n')
#for(i in species.list) {
#	cat(i,'\n')
#	species_wrap_fun(i)}

# chunking species into multiple lists
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
chunk_list <- chunk2(species.list,3)

# And building trycatch around the accuracy checking function
try_catch_function <-
  function(k) {
    tryCatch(species_wrap_fun(k),
             error = function(e) { # Try catch error function
               tmp_species <- gsub('.*ESH_RCPs/Raw','',k) 
               
               # which species did it not work for?
               cat('Whoopsies:',tmp_species, '\n')
               cat(tmp_species,'\n')
               cat('Saving file:', paste0(getwd(),"/ESH_RCPs/Weighted_Threshold_Accuracy/",tmp_species,"_error.csv"), '\n')
               # tracking species it did not work for
               
               write.csv(data.frame(taxa_species = tmp_species, error = 'Error'),
                         paste0(getwd(),"/ESH_RCPs/Weighted_Threshold_Accuracy/",tmp_species,"_error.csv"),
                         row.names = FALSE)
               
             }) # End try catch 
  }


length(chunk_list)
head(chunk_list)
#for(chunk in chunk_list) {
#cat('Species: ',chunk,'\n')
#try_catch_function(chunk)
#}


# try_catch_function(species.list[1])
out.df <- mclapply(chunk_list[[1]], try_catch_function, mc.cores=10)


cat('\n End Script \n')
