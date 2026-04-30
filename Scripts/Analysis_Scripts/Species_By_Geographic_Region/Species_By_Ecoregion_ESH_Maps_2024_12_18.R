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
ecoregion_map <- raster(paste0(getwd(),'/Ecoregions_Feb2023/TNC_Ecoregions_Map.tif'))
template_map <- raster(paste0(getwd(),'/Global Mollweide Maps/MollweideCountryID_1.5km.tif'))
ecoregion_map[is.na(ecoregion_map) & !is.na(template_map)] <- -99
###
# LIST OF SDM FILES AND SPECIES

# List of taxa in the analysis
taxa_list <-
  list.files(paste0(getwd(),'/ESH_Tifs_12Oct'),
             full.names = TRUE) 

# Randomising order of taxa list
taxa_list <- sample(taxa_list, 4)


# List of sdm files in 2020 and 2050
sdm_files <-
  do.call(c,mclapply(taxa_list, list.files, mc.cores = 4, full.names = TRUE)) %>%
  .[grepl('.tif$',.)]

# Getting taxa list
taxa_list <-
  gsub('.*ESH_Tifs_12Oct/','',taxa_list)


###
# MAKING FUNCTION TO LOOP ACROSS SPECIES
species_loop_fun <-
  function(ss) {
    cat(which(species_list %in% ss),':',ss,'\n')
    
    # dir.create <- paste0(getwd(),'/Analyses//Temp_Folder_Species_By_Biodiv_Hotspot/')
    
    already_have_tmp <-
      list.files(paste0(getwd(),'/Analyses/Temp_Folder_Species_By_Ecoregion_ESH/'),
                 full.names = TRUE) %>%
      gsub('.*Temp_Folder_Species_By_Ecoregion_ESH//','',.) %>%
      gsub(paste0('_esh_map','.*'),'',.)
    if(!(ss %in% already_have_tmp)) {
      # List of map files for species
      species_files <-
        taxa_files %>%
        .[grepl(paste0(ss,'.tif'),.)]
      
      # Creating empty df
      out_df <- data.frame()
        
        # Importing raster
        year_raster <- raster(species_files)
        raster_values <- unique(getValues(year_raster))
        
        # If statement to minimise time requirements - 2 if statements
        if(!is.na(maxValue(year_raster))) {
          if(maxValue(year_raster) > 0) {
            # Cropping ecoregions
            tmp_ecoregion <- raster::crop(ecoregion_map, year_raster)
            
            # Removing NAs
            if(length(raster_values) >= 3) {
              tmp_ecoregion[is.na(year_raster)] <- NA
            } else {
              tmp_ecoregion[!(year_raster %in% 1)] <- NA
            }
            
            
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
                        mutate(esh_map= 'esh_map',
                               species = ss,
                               taxa = tt))
            }
          } # End if statement to minimise time requirements 
        }
      
      # Returning data frame
      write.csv(out_df,
                paste0(getwd(),'/Analyses/Temp_Folder_Species_By_Ecoregion_ESH/',ss,'_',tt,'_esh_map_',Sys.Date(),'.csv'),
                row.names = FALSE)
    } # End if statement for checking whether we have the species
  } # End function for species loop

###
# MAKING LOOP TO GO ACROSS MIGRATION, TAXA, AND THRESHOLD

    ###
    # Loop across taxa
    for(tt in (taxa_list)) {
      cat(tt, '\n')
      # Limiting files to taxa
      taxa_files <-
        sdm_files %>%
        .[grepl(tt,.)]
      
      # List of species
      species_list <-
        taxa_files %>%
        gsub(paste0('.*',tt,'/'),'',.) %>%
        gsub('_SSP.*','',.) %>%
	gsub('.tif$','',.) %>%
        unique()
      
      already_have <-
        list.files(paste0(getwd(),'/Analyses/Temp_Folder_Species_By_Ecoregion_ESH/'),
                   full.names = TRUE) %>%
        .[grepl(tt,.)] %>%
        .[grepl('esh_map',.)] %>%
        gsub('.*Temp_Folder_Species_By_Ecoregion_ESH//','',.) %>%
        gsub(paste0('_esh_map_','.*'),'',.) %>%
	gsub(paste0('_',tt,'.*'),'',.)
      
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

