#!/usr/bin/env Rscript


#####
# script to aggregate glob2loc patch level results
# in a way that doesn't use too much space
#####


###
# libraries
library(plyr)
library(dplyr)
library(parallel)
library(stringr)

###
# setting working directory
files_stored_wd <- '/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity'
files_write_wd <- '/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity'
# setwd('/Users/michael/Desktop/Research/Glob2Loc/tmp_cheetah_files')

###
# list of taxas
taxa_list <- c('Mammals','Birds','Amphibians','Reptiles')

###
# list of thresholds
threshold_list <- c('specsens','prevalence')

### 
# list of scenarios we want
scen_list <-
  c('esh.future.extent.esh.exp.int.urb',
    'esh.future.extent.esh.exp',
    'esh.future.extent.esh.cropexp',
    'esh.future.extent.esh.pastexp',
    'esh.future.extent.esh.int',
    'esh.future.extent.esh.urb',
    'esh.current.extent.esh')

###
# species_files <- 
#   list.files(paste0(getwd(),'/Outputs/CSV_File_Outputs/BAU/',),
#              pattern = 'Spatial_Population') %>%
#   .[grepl('prevalence',.)] %>%
#   .[!grepl('nomigration',.)]

# creating directories for file storage
dir.create('/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/Outputs')
dir.create('/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/Outputs/CSV_File_Outputs_Migration')
dir.create('/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/Outputs/Managed_CSV_Files_Migration')
dir.create('/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/Outputs/Managed_CSV_Files_Migration/Abs_Impacts/')
dir.create('/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/Outputs/Managed_CSV_Files_Migration/Prop_Impacts/')
lapply(taxa_list, function(tt) {dir.create(paste0(files_write_wd,'/Outputs/Managed_CSV_Files_Migration/Prop_Impacts/',tt))})
lapply(taxa_list, function(tt) {dir.create(paste0(files_write_wd,'/Outputs/Managed_CSV_Files_Migration/Abs_Impacts/',tt))})

###
# function to get data summaries for each species in each year
species_function <-
  function(df) {

    # summarising
    tot_dat_species <-
      df %>%
      filter(scenario %in% scen_list) %>%
      mutate(patches_area_morethan_1km = ifelse(habitat_availability_sqkm >= 1, 1, 0),
             patches_area_morethan_5km = ifelse(habitat_availability_sqkm >= 5, 1, 0),
             patches_area_morethan_25km = ifelse(habitat_availability_sqkm >= 25, 1, 0),
	     patches_area_morethan_100km = ifelse(habitat_availability_sqkm >= 100, 1, 0)) %>%
      mutate(patches_populations_morethan_50 = ifelse(estimated_population >= 50, 1, 0),
             patches_populations_morethan_250 = ifelse(estimated_population >= 250, 1, 0),
	     patches_populations_morethan_500 = ifelse(estimated_population >= 500, 1, 0),
             patches_populations_morethan_1000 = ifelse(estimated_population >= 1000, 1, 0),
	     patches_populations_morethan_5000 = ifelse(estimated_population >= 5000, 1, 0)) %>%
      dplyr::group_by(scenario, year, binomial) %>%
      dplyr::summarise(num_patches = max(patches,na.rm=TRUE),
                       tot_pop = sum(estimated_population),
                       tot_area = sum(habitat_availability_sqkm),
                       patches_area_morethan_1km = sum(patches_area_morethan_1km),
                       patches_area_morethan_5km = sum(patches_area_morethan_5km),
                       patches_area_morethan_25km = sum(patches_area_morethan_25km),
		       patches_area_morethan_100km = sum(patches_area_morethan_100km),
                       patches_populations_morethan_50 = sum(patches_populations_morethan_50),
                       patches_populations_morethan_250 = sum(patches_populations_morethan_250),
		       patches_populations_morethan_500 = sum(patches_populations_morethan_500),
                       patches_populations_morethan_1000 = sum(patches_populations_morethan_1000),
		       patches_populations_morethan_5000 = sum(patches_populations_morethan_5000))
    
    # patch with largest pop and area
    largest_patch <-
      left_join(df %>%
                  filter(scenario %in% scen_list) %>%
                  dplyr::group_by(scenario) %>%
                  arrange(desc(estimated_population)) %>%
                  slice(1) %>%
                  dplyr::select(scenario, pop_largest_patch = estimated_population),
                df %>%
                  filter(scenario %in% scen_list) %>%
                  dplyr::group_by(scenario) %>%
                  arrange(desc(habitat_availability_sqkm)) %>%
                  slice(1) %>%
                  dplyr::select(scenario, area_largest_patch = habitat_availability_sqkm),
                          by = join_by(scenario))
    
    # patch with largest prop of population
    prop_dat_species <-
      left_join(df,tot_dat_species, by = join_by(scenario, binomial, year)) %>%
      filter(scenario %in% scen_list) %>%
      mutate(prop_pop = estimated_population / tot_pop,
             prop_area = habitat_availability_sqkm / tot_area)
    
    # num patches needed to account for > 50, 75, 90, and 95% of population
    prop_dat_species_area <-
      prop_dat_species %>%
      dplyr::group_by(scenario) %>%
      arrange(-prop_area) %>%
      dplyr::summarise(cumsum_area = cumsum(prop_area)) %>%
      mutate(area_prop_50 = ifelse(cumsum_area <= .5, 1, 0),
             area_prop_75 = ifelse(cumsum_area <= .75, 1, 0),
             area_prop_90 = ifelse(cumsum_area <= .90, 1, 0),
             area_prop_95 = ifelse(cumsum_area <= .95, 1, 0)) %>%
      mutate(area_prop_50 = ifelse(is.na(area_prop_50),-1,area_prop_50),
             area_prop_75 = ifelse(is.na(area_prop_75),-1,area_prop_75),
             area_prop_90 = ifelse(is.na(area_prop_90),-1,area_prop_90),
             area_prop_95 = ifelse(is.na(area_prop_95),-1,area_prop_95)) %>%
      dplyr::group_by(scenario) %>%
      dplyr::summarise(num_patches_for_50percent_area = sum(area_prop_50) + 1,
                       num_patches_for_75percent_area = sum(area_prop_75) + 1,
                       num_patches_for_90percent_area = sum(area_prop_90) + 1,
                       num_patches_for_95percent_area = sum(area_prop_95) + 1)
    
    # and likewise for area
    prop_dat_species_pop <-
      prop_dat_species %>%
      dplyr::group_by(scenario) %>%
      arrange(-prop_pop) %>%
      dplyr::summarise(cumsum_pop = cumsum(prop_pop)) %>%
      mutate(pop_prop_50 = ifelse(cumsum_pop <= .5, 1, 0),
             pop_prop_75 = ifelse(cumsum_pop <= .75, 1, 0),
             pop_prop_90 = ifelse(cumsum_pop <= .90, 1, 0),
             pop_prop_95 = ifelse(cumsum_pop <= .95, 1, 0)) %>%
      mutate(pop_prop_50 = ifelse(is.na(pop_prop_50),-1,pop_prop_50),
             pop_prop_75 = ifelse(is.na(pop_prop_75),-1,pop_prop_75),
             pop_prop_90 = ifelse(is.na(pop_prop_90),-1,pop_prop_90),
             pop_prop_95 = ifelse(is.na(pop_prop_95),-1,pop_prop_95)) %>%
      dplyr::group_by(scenario) %>%
      dplyr::summarise(num_patches_for_50percent_pop = sum(pop_prop_50) + 1,
                       num_patches_for_75percent_pop = sum(pop_prop_75) + 1,
                       num_patches_for_90percent_pop = sum(pop_prop_90) + 1,
                       num_patches_for_95percent_pop = sum(pop_prop_95) + 1)
    
    # merging back in
    tot_dat_species <-
      left_join(tot_dat_species,
                prop_dat_species_area, by = join_by(scenario)) %>%
      left_join(.,
                prop_dat_species_pop, by = join_by(scenario)) %>%
      left_join(.,
                largest_patch, by = join_by(scenario))

    # reutnring data frame
    return(tot_dat_species)
  } # end function for species




# function to parallel lapply across species
species_loop_function <-
  function(ss) {
    # Tracking progress
    cat('Beginning:',tt,ss,'\n')
    
    # getting list of files
    species_files <-
      taxa_files[grepl(ss,taxa_files)]
    
    # list of years for the species
    years_list <- str_extract(species_files,'20[0-9]{2,2}')
    
    # check to see if pop is numeric
    pop_update <- 'NO'

    years_list <- c(2020, years_list[!(years_list %in% 2020)])
    
    # and now looping through years
    for(yy in years_list) {
      # making data frame for species
      if(yy %in% 2020) {
        out_df_nomigrate <- data.frame()
        #out_df_migrate <- data.frame()
        out_df_nomigrate_prop <- data.frame()
        #out_df_migrate_prop <- data.frame()
      } # end if statement making data frame
      
      # getting 2020 file as baseline
      # aggregating into single row by scenario
      dat_species <-
        read.csv(species_files[grepl(yy,species_files)]) %>%
	dplyr::rename(habitat_availability_sqkm = hab_avail,
		      tot_pop = pop_abund) %>%
	dplyr::group_by(pop_id) %>%
	dplyr::summarise(habitat_availability_sqkm = sum(habitat_availability_sqkm),
			 estimated_population = sum(tot_pop)) %>%
	dplyr::rename(patch_no_migrate = pop_id) %>%
	mutate(scenario = 'esh.future.extent.esh.exp.int.urb') %>%
	mutate(year = yy,
	       binomial = paste0(tt,'/',ss))
			 
      
      dat_species <-
        dat_species %>%
        mutate(estimated_population = ifelse(grepl('no remaining habitat',habitat_availability_sqkm),0,estimated_population)) %>%
        mutate(patch_no_migrate = ifelse(grepl('no remaining habitat',habitat_availability_sqkm),0,patch_no_migrate)) %>%
        mutate(habitat_availability_sqkm = ifelse(grepl('no remaining habitat',habitat_availability_sqkm),0,habitat_availability_sqkm)) %>%
        mutate(estimated_population = as.numeric(estimated_population),
               habitat_availability_sqkm = as.numeric(habitat_availability_sqkm),
               patch_no_migrate = as.numeric(patch_no_migrate))
      
      
      # for no patch migration
      dat_no_migrate <- species_function(dat_species %>%
                                           dplyr::rename(patches = patch_no_migrate)) %>%
        mutate(prop_pop_largest_patch = pop_largest_patch / tot_pop,
               prop_area_largest_patch = area_largest_patch / tot_area)
      
      # for patch migration
      #dat_migrate <- 
      #  species_function(dat_species %>%
      #                     dplyr::rename(patches = patch_migrate) %>%
      #                     dplyr::group_by(scenario, patches, binomial, year) %>%
      #                     dplyr::summarise(estimated_population = sum(estimated_population),
      #                                      habitat_availability_sqkm = sum(habitat_availability_sqkm))) %>%
      #  mutate(prop_pop_largest_patch = pop_largest_patch / tot_pop,
      #         prop_area_largest_patch = area_largest_patch / tot_area)
      
      # rbinding absolute impacts
      out_df_nomigrate <- rbind(out_df_nomigrate, dat_no_migrate)
      #out_df_migrate <- rbind(out_df_migrate, dat_migrate)
      
      # getting proportional impacts relative to 2020
      dat_nomigrate_prop <- dat_no_migrate
      #dat_migrate_prop <- dat_migrate
      
      # looping to divide
      for(ii in which(!(names(dat_no_migrate) %in% c('scenario','year','binomial','pop_largest_patch','area_largest_patch')))) {
        # for no migration
        dat_nomigrate_prop[,ii] <-
          dat_nomigrate_prop[,ii] /
          as.numeric(out_df_nomigrate[out_df_nomigrate$scenario %in% 'esh.future.extent.esh.exp.int.urb' &
                                        out_df_nomigrate$year %in% 2020,
                                      ii])
        
        # for migration
        # dat_migrate_prop[,ii] <-
        #   dat_migrate_prop[,ii] /
        #   as.numeric(out_df_migrate[out_df_migrate$scenario %in% 'esh.future.extent.esh.exp.int.urb' & 
        #                               out_df_migrate$year %in% 2020,
        #                             ii])
      } # end statement for prop impacts
      
      # and rbinding proportional values
      out_df_nomigrate_prop <- rbind(out_df_nomigrate_prop, dat_nomigrate_prop)
      #out_df_migrate_prop <- rbind(out_df_migrate_prop, dat_migrate_prop)
      
    } # end year loop
    
    # updating to pop to NAs if needed
   out_df_nomigrate <-
	  out_df_nomigrate %>%
	 mutate(pop_dispersal = 'full_patch_migration') 

 out_df_nomigrate_prop <-
          out_df_nomigrate_prop %>%
         mutate(pop_dispersal = 'full_patch_migration')
    
    # writing csv files
    write.csv(out_df_nomigrate,
              paste0(files_write_wd,'/Outputs/Managed_CSV_Files_Migration/Abs_Impacts/',tt,'/',ss,'_full_patch_migration_',thresh,'.csv'),
              row.names = FALSE)
    
    #write.csv(out_df_migrate,
    #          paste0(files_write_wd,'/Outputs/Managed_CSV_Files/Abs_Impacts/',tt,'/',ss,'_patchmigration_',thresh,'.csv'),
    #         row.names = FALSE)
    
    write.csv(out_df_nomigrate_prop,
              paste0(files_write_wd,'/Outputs/Managed_CSV_Files_Migration/Prop_Impacts/',tt,'/',ss,'_full_patch_migration_',thresh,'.csv'),
              row.names = FALSE)
    
    #write.csv(out_df_migrate_prop,
    #          paste0(files_write_wd,'/Outputs/Managed_CSV_Files/Prop_Impacts/',tt,'/',ss,'_patchmigration_',thresh,'.csv'),
    #          row.names = FALSE)
    
    # removing files
    rm(out_df_nomigrate,
       out_df_nomigrate_prop)
    
    # Tracking progress
    cat('Finished:',tt,ss,'\n')
  } # End function to lapply across species



###
# Setting option to avoid dplyr giving indications on summarise and joining
options(dplyr.summarise.inform = FALSE) # get rid of info on summarising data frames
options(dplyr.left_join.inform = FALSE)

# Try catch function to catch any errors, where needed
try_catch_function <-
  function(ss) {
    tryCatch(species_loop_function(ss),
             error = function(e) cat('Error: ',ss,'\n'))
  }


# trying in lapply
###
# looping through thresholds
for(thresh in threshold_list[1]) {
  # looping through taxa
  for(tt in taxa_list[4]) {
	  cat(thresh, tt, '\n')
    # getting list of species files
    taxa_files <-
      list.files(paste0(files_stored_wd,'/Outputs/CSV_File_Outputs_Migration/',tt),
                 full.names = TRUE) %>%
      .[grepl(thresh,.)] %>%
      .[!grepl('Climate_Suitable',.)]

    
    # Full list of species
    species_list_2020 <-
      taxa_files %>%
      gsub('Spatial_Population_Estimates_csv_','',.) %>%
      .[grepl('2020',.)] %>%
      gsub('_[0-9]{4,4}.*','',.) %>%
      gsub(paste0('.*',tt,'/'),'',.) %>%
      gsub(paste0('_',thresh),'',.) %>%
      unique()

    
    species_list_2050 <-
      taxa_files %>%
      gsub('Spatial_Population_Estimates_csv_','',.) %>%
      .[grepl('2050',.)] %>%
      gsub('_[0-9]{4,4}.*','',.) %>%
      gsub(paste0('.*',tt,'/'),'',.) %>%
       gsub(paste0('_',thresh),'',.) %>%
      unique()

    
    species_list <-
      species_list_2020[species_list_2020 %in% species_list_2050]
   

    # Species we have
    species_have <-
      list.files(paste0(files_write_wd,'/Outputs/Managed_CSV_Files_Migration/Prop_Impacts/',tt,'/')) %>%
      .[grepl(thresh,.)] %>%
      .[grepl('_full_patch',.)] %>%
      .[!grepl('Climate_Suitable',.)] %>%
      gsub('_full_patch.*','',.) %>%
      gsub('_patchmigration.*','',.) %>%
      unique()
    
    # Now limiting to only species we need
    species_list <- 
      species_list[!(species_list %in% species_have)]

# List of species where biodiversity analysis is completed
# adding this to make sure we don't accidentally include species with old results
#species_completed <- 
#	list.files(paste0(files_stored_wd,'/Outputs/CSV_File_Outputs_Migration/',tt),
#                 full.names = TRUE) %>%
#      .[grepl('_completed',.)] %>%
#      .[grepl(thresh,.)] %>%
#      .[!grepl('Climate_Suitable',.)] %>%
#      .[grepl('2020|2050',.)] %>%
#      gsub(paste0('.*',tt,'/'),'',.) %>%
#      gsub('_20.*','',.)

#species_completed <-
#	data.frame(species = species_completed) %>%
#	dplyr::group_by(species) %>%
#	dplyr::summarise(count = n()) %>%
#	filter(count %in% 2) %>%
#	dplyr::select(species)

#species_list <- species_list[species_list %in% species_completed$species]
    


# lapply across species
    mclapply(species_list, try_catch_function, mc.cores = 10)
  } # end taxa loop
}

