#!/usr/bin/env Rscript

###
# LIBRARIES
library(raster)
library(plyr)
library(dplyr)
library(parallel)
library(stringr)
###
# SETTING WD
setwd('/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity')

###
# ECOREGION MAP
# ecoregion_map <- raster('/Users/michael/Desktop/TNC_Ecoregions_Map.tif')
ecoregion_map <- raster(paste0(getwd(),'/Ecoregions_Feb2023/Biodiversity_Hotspots_Raster.tif'))
template_map <- raster(paste0(getwd(),'/Global Mollweide Maps/MollweideCountryID_1.5km.tif'))
ecoregion_map[is.na(ecoregion_map) & !is.na(template_map)] <- -99
###
# LIST OF SDM FILES AND SPECIES

# List of taxa in the analysis
taxa_list <- 
  list.files(paste0(getwd(),'/ESH_RCPs/SSP2-4.5'),
             full.names = TRUE) %>%
  .[!grepl('Updated',.)] %>%
  .[!grepl('Climate',.)]

# List of sdm files in 2020 and 2050
sdm_files <-
  do.call(c,mclapply(taxa_list, list.files, mc.cores = 4, full.names = TRUE)) %>%
  .[grepl('2020|2050',.)]

# Getting taxa list
taxa_list <-
  gsub('.*SSP2-4.5/','',taxa_list) 

###
# MAKING FUNCTION TO LOOP ACROSS SPECIES
species_loop_fun <-
  function(ss) {
    cat(which(species_list %in% ss),':',ss,'\n')

dir.create <- paste0(getwd(),'/Analyses//Temp_Folder_Species_By_Biodiv_Hotspot/')

already_have_tmp <-
        list.files(paste0(getwd(),'/Analyses//Temp_Folder_Species_By_Biodiv_Hotspot/'),
                   full.names = TRUE) %>%
.[grepl(mm,.)] %>%
.[grepl(thresh,.)] %>%
.[grepl(tt,.)] %>%
.[!grepl('esh_map',.)] %>%
gsub('.*Temp_Folder_Species_By_Biodiv_Hotspot//','',.) %>%
gsub(paste0('_',mm,'.*'),'',.)

    if(!(ss %in% already_have_tmp)) {
    # List of map files for species
    species_files <- 
      taxa_files %>%
      .[grepl(ss,.)]
    
    # Creating empty df
    out_df <- data.frame()

    years <- str_extract(species_files,'2020|2050')
      
    
    # Looping through species files
    for(yy in years) {
      # Getting file for year
      year_file <- 
        species_files %>%
        .[grepl(yy,.)]
      
      # Importing raster
      year_raster <- raster(year_file)
      
      # If statement to minimise time requirements
      if(maxValue(year_raster) %in% 1) {
        # Cropping ecoregions
        tmp_ecoregion <- raster::crop(ecoregion_map, year_raster)
        
        # Removing NAs
        tmp_ecoregion[!(year_raster %in% 1)] <- NA
        
        # Getting unique values
        eco_vals <- sort(unique(getValues(tmp_ecoregion)))
        
        # If statement if no overlap with land
        if(length(eco_vals) >= 1) {
          # rbinding to dataframe
          out_df <-
            rbind(out_df,
		  data.frame(esh_values = getValues(year_raster),
           ecoregion_id = getValues(tmp_ecoregion)) %>%
  filter(!is.na(esh_values),
         !is.na(ecoregion_id)) %>%
  dplyr::group_by(ecoregion_id) %>%
  dplyr::summarise(prop_cells_in_ecoregion = n()) %>%
  dplyr::ungroup() %>%
  mutate(prop_cells = prop_cells_in_ecoregion / sum(prop_cells_in_ecoregion)) %>%
  mutate(climate_migration = mm,
           year = yy,
           taxa = tt,
           thresh = thresh,
           species = ss))
        }
      } # End if statement to minimise time requirements
    } # End loop across years
    
    # Returning data frame

    write.csv(out_df,
                paste0(getwd(),'/Analyses//Temp_Folder_Species_By_Biodiv_Hotspot/',ss,'_',mm,'_',thresh,'_',tt,'_',Sys.Date(),'.csv'),
                row.names = FALSE)
  } # End if statement for checking whether we have the species
  } # End function for species loop
###
# ASSUMPTIONS ON CLIMATE MIGRATION, THRESHOLDS, ETC
clim_migration <- 
  c('climate migration limits')

threshold_list <-
  c('specsens','prevalence')

###
# MAKING LOOP TO GO ACROSS MIGRATION, TAXA, AND THRESHOLD
for(mm in clim_migration) {
  # migration search term
  mig_search <-
    ifelse(mm %in% 'no climate migration limits',
           '__',
           'migrate')
  # limiting to assumptions on migration
  migration_files <-
    sdm_files %>%
    .[grepl(mig_search,.)]
  
  ### 
  # Loop across thresholds
  for(thresh in rev(threshold_list)) {
    # getting files for thresholds
    thresh_files <-
      migration_files %>%
      .[grepl(thresh,.)]
    
    ###
    # Loop across taxa
    for(tt in rev(taxa_list)) {
      cat(mm, thresh, tt, '\n')
      # Limiting files to taxa
      taxa_files <-
        thresh_files %>%
        .[grepl(tt,.)]
      
      # List of species
      species_list <- 
        taxa_files %>%
        gsub(paste0('.*',tt,'/'),'',.) %>%
        gsub('_SSP.*','',.) %>%
        unique()

already_have <-
	list.files(paste0(getwd(),'/Analyses//Temp_Folder_Species_By_Biodiv_Hotspot/'),
		   full.names = TRUE) %>%
.[grepl(mm,.)] %>%
.[grepl(thresh,.)] %>%
.[grepl(tt,.)] %>%
.[!grepl('esh_map',.)] %>%
gsub('.*Temp_Folder_Species_By_Biodiv_Hotspot//','',.) %>%
gsub(paste0('_',mm,'.*'),'',.)

species_list <- species_list[!(species_list %in% already_have)] 
species_list <- sample(species_list, length(species_list),replace = FALSE)

cat('Species Completed:',length(already_have),'\n')
cat('Species Need:',length(species_list),'\n')

if(length(species_list) > 0) {
trycatch_function <-
  tryCatch(species_loop_fun,
           error = function(e) cat('ERROR WHOOPS:',ss,'\n'))
      
      # Running function in parallel
      t1 = Sys.time()
                #mclapply(species_list,
                #         trycatch_function,
                #         mc.cores = 5)
for(ss in species_list){trycatch_function(ss)}

      #print(head(out_df_save))
      Sys.time() - t1
} # End check to loop across species - only need to do this if there are species left to complete!
      
      # for(ss in species_list[158:1000]) {species_loop_fun(ss)}
      
      # Writing file
      #write.csv(out_df_save,
      #          paste0(getwd(),'/Analyses/Temp_Folder_Species_By_Ecoregion/',ss,'_',mm,'_',thresh,'_',tt,'_',Sys.Date(),'.csv'),
      #          row.names = FALSE)
      
    } # End loop for taxa
  } # End loop for threshold types
} # End loop for climate migration
