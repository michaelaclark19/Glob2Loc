#####
# SCRIPT FOR SENSITIVITY ANALYSIS ACROSS SDM MODELS
# TO SEE HOW MANY DO (OR DON'T) HAVE ESTIMATED CLIMATE SUITABLE HABITAT FOR EACH SPECIES
#####

###
# Libraries
library(plyr)
library(dplyr)
library(parallel)
library(raster)
library(stringr)
library(future.apply)



###
# List of taxa
taxa_list <- c('Amphibians','Birds','Mammals','Reptiles')
taxa_list <- c('Reptiles')

###
# Getting folder for SDM files
sdm_files <-
  do.call(c,
          lapply(paste0('/migration/pubh-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/ESH_RCPs/Raw/',taxa_list),
                 list.files,
                 full.names = TRUE))
  

###
# Getting files for SDM thresholds
thresh_files <-
  list.files(path = '/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/ESH_RCPs/Mod_Accuracy',
             full.names = TRUE) %>%
  .[grepl('_Thresh',.)]

###
# Getting files for mod accuracy
acc_files <-
  list.files(path = '/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/ESH_RCPs/Mod_Accuracy',
             full.names = TRUE) %>%
  .[!grepl('_Thresh',.)]

###
# List of species
species_list <-
  acc_files %>%
  gsub('.*_Accuracy_','',.) %>%
  gsub('.csv','',.)

###
# List of ESH maps
esh_maps <- 
  do.call(c,
          lapply(paste0('/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/ESH_Tifs_12Oct/',taxa_list),
                 list.files,
                 full.names = TRUE))

###
# Getting realm biome combinations
realm_biomes <- raster('/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/Ecoregions_Feb2023/biomes_realms_raster.tif')



###
# Making function to use for lapply/mclapply
species_fun <-
  function(ss) {
    ###
    # Tracking progress
    cat(ss,'\n')
    ###
    # Making empty data frame
    out_df <- data.frame()
    
    ###
    # Name of specices
    taxa_name <- str_extract(ss,'Amphibians|Birds|Mammals|Reptiles')
    species_name <- str_remove(ss,'_Amphibians|_Birds|_Mammals|_Reptiles')
    
    ###
    # SDM files for the species
    species_sdms <- 
      sdm_files %>%
      .[grepl(species_name,.)] %>%
      .[grepl(taxa_name,.)]
    
    # Getting files in 2005 and 2050
    species_sdms <-
      species_sdms %>%
      .[grepl('2014|2041',.)] %>%
      .[grepl('Historic|SSP2',.)]
    
    # Checking whether we need to loop
    if(length(species_sdms) >= 2) {
      
      ###
      # Accuracy for species
      species_acc <-
        acc_files %>%
        .[grepl(species_name,.)] %>%
        .[grepl(taxa_name,.)] %>%
        read.csv(.)
      
      # changing to numeric
      species_acc[,c("bioclim", "maxent", "rf", "glm.gaussian")] <- as.numeric(species_acc[,c("bioclim", "maxent", "rf", "glm.gaussian")])
      
      # Checking whether any models accurate enough for further analyiss
      which_mods <- 
        unique(str_extract(species_sdms,
                           'bioclim|maxent|random.forest|glm.gaussian'))
      ###
      # Thresholds for species
      species_thresh <-
        thresh_files %>%
        .[grepl(species_name,.)] %>%
        .[grepl(taxa_name,.)] %>%
        read.csv()
      
      
      # Looping across mod types
      for(mm in which_mods) {
        cat(mm,'\n')
        
        ###
        # Updating rf to work with threshold data
        mod_name = ifelse(mm %in% 'random.forest','rf',mm)
        
        # 2010 model...
        tmp_mod_2010 <- 
          species_sdms %>%
          # .[grepl(tt,.)] %>%
          .[grepl('2014',.)] %>%
          .[grepl(mm,.)] %>%
          .[1] %>%
          raster(.)
        
        # 2050 model...
        tmp_mod_2050 <- 
          species_sdms %>%
          # .[grepl(tt,.)] %>%
          .[grepl('2041',.)] %>%
          .[grepl(mm,.)] %>%
          .[1] %>%
          raster(.)
        
        ###
        # For both models, removing realm biome combos not in the species current habitat
        
        # esh map
        tmp_esh <-
          esh_maps %>%
          .[grepl(species_name,.)] %>%
          .[grepl(taxa_name,.)] %>%
          raster(.)
        
        tmp_biomes <- raster::crop(realm_biomes, tmp_esh)
        tmp_biomes[is.na(tmp_esh)] <- NA
        tmp_biome_values <- unique(getValues(tmp_biomes)) %>% .[.>1] %>% .[!is.na(.)]
        
        ### 
        # And clipping habitat maps
        tmp_biomes <- raster::crop(realm_biomes, tmp_mod_2010)
        tmp_biomes[!(tmp_biomes %in% tmp_biome_values)] <- NA
        tmp_mod_2010[is.na(tmp_biomes)] <- NA
        
        tmp_biomes <- raster::crop(realm_biomes, tmp_mod_2050)
        tmp_biomes[!(tmp_biomes %in% tmp_biome_values)] <- NA
        tmp_mod_2050[is.na(tmp_biomes)] <- NA
        
        ### And getting values
        tmp_vals_2010 <- getValues(tmp_mod_2010)
        tmp_vals_2050 <- getValues(tmp_mod_2050)
        
        # Looping across threshold types
        for(tt in c('prevalence','equal_sens_spec')) {
          
          # cat(tt,'\n')
          ###
          # Getting thresholds for the model
          thresh <-
            species_thresh %>%
            filter(model_name %in% mm) %>%
            .[,tt]
          
          ### 
          # Getting number of values below threshold
          # And saving as a dataframe
          tmp_df <-
            data.frame(taxa = taxa_name,
                       species_name = species_name,
                       mod = mm,
                       thresh = tt,
                       year = c(2010,2050),
                       values = c(sum(tmp_vals_2010 >= thresh, na.rm = TRUE), 
                                  sum(tmp_vals_2050 >= thresh, na.rm = TRUE)))
          
          ###
          # Rbinding outcomes
          out_df <-
            rbind(out_df,
                  tmp_df)
          
        } # End loop through thresholds
      } # End loop through model type
    } # End if statement to see whether we need to do analysis on species
    
    ###
    # Returning file
    return(out_df)
  } # End function


###
# Testing
species_list_amps <- species_list %>% .[grepl('Amphibians',.)] %>% .[grepl('Brachycephalus|Ischnoc',.)]

###
# Testing in parallel with future_lapply
plan(multicore, workers = 10)
t1 = Sys.time()
save_df <-
  do.call(rbind,
          future_lapply(species_list_amps,
                        species_fun,
                        future.seed = NULL))
Sys.time() - t1
plan(sequential)


for(ss in species_list) {species_fun(ss)}
###
# Running in parallel
plan(multicore, workers = 20)
# t1 = Sys.time()
save_df <-
  do.call(rbind,
          future_lapply(species_list,
                        species_fun,
                        future.seed = NULL))
Sys.time() - t1
# plan(sequential)

write.csv(save_df,
          '/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/Analyses/SDM_Mod_Uncertainty.csv',
          row.names = FALSE)

