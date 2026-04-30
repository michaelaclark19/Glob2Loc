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
library(countrycode)

#####
### General data management
# Setting wd
setwd('/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/')

# How many cores to use
mc_cores = 10

# name of country for this script
country <- 'AUS'
country_iso <- countrycode(country, origin = 'iso3c', destination = 'iso3n')

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
thresh_list <- c('prevalence')#,'specsens')

# Getting outcome types
outcome_list <- c('hab_availability','pop_density')

# Getting list of output files
# full_file_list <-
#   do.call(c,
#           lapply(paste0(getwd(),'/Outputs/Raster_Outputs/BAU/',taxa_list),
#                  list.files,
#                  full.names=TRUE))


# making template raster
country_rast <- raster(paste0(getwd(),"/Global Mollweide Maps/MollweideCountryID_1.5km.tif"))


# Getting global map
country_rast <- raster(paste0(getwd(),"/Global Mollweide Maps/MollweideCountryID_1.5km.tif"))
country_rast[!(country_rast %in% country_iso)] <- NA

# Converting to matrix to find rows and cols with NAs
country_mat <- as.matrix(country_rast)

rows_all_na <- apply(country_mat, 1, function(x) all(is.na(x)))
cols_all_na <- apply(country_mat, 2, function(x) all(is.na(x)))

# Getting max and mins for extent
min_y <- yFromRow(country_rast, max(which(rows_all_na %in% FALSE)) + 20)
max_y <- yFromRow(country_rast, min(which(rows_all_na %in% FALSE)) - 20)
min_x <- xFromCol(country_rast, min(which(cols_all_na %in% FALSE)) - 20)
max_x <- xFromCol(country_rast, max(which(cols_all_na %in% FALSE)) + 20)

# And cropping
country_rast_cropped <-
  crop(country_rast,
       extent(min_x,max_x,
              min_y,max_y)) # xmin, xmax, ymin, ymax

country_rast_cropped[!is.na(country_rast_cropped)] <- 0

###
# Getting body mass data
# bird/amp/mammal body mass
body_mass <-
  read.csv(paste0(getwd() %>% gsub('pubh-glob2loc','ouce-glob2loc',.),"/Other Data Inputs/Pop Density Inputs/Body Mass Estimates 20February2020 Updated Taxonomy.csv"), stringsAsFactors = FALSE) %>% # Importing body mass file
  dplyr::select(order = Order, family = Family, genus = Genus, species_merge = binomial, est_mass_kg, Family_Mass_kg,
Genus_Mass_kg) %>%  # Only keeping necessary columns
  mutate(binomial = species_merge)

# reptile body mass
rep_body_mass <-
  read.csv('/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/Other Data Inputs/Pop Density Inputs/reptile_body_masses.csv',
           stringsAsFactors = FALSE) %>%
  dplyr::rename(order = Order,
                family = Family) %>%
  mutate(est_mass_kg = `mass..g.` / 1000) %>%
  mutate(Family_Mass_kg = est_mass_kg,
         Genus_Mass_kg = est_mass_kg) %>%
  mutate(binomial = gsub(' ','_', binomial)) %>%
  #mutate(binomial = gsub(' ','_', `X...binomial`)) %>%
  mutate(genus = gsub('_.*','',binomial),
         species_merge = binomial)

  # rbinding the body mass estimates for reps and non reps
  body_mass_merge <-
  rbind(body_mass,
        rep_body_mass[,names(body_mass)])



#####
###
# Making function
# This will be parallised later
species_funcion <-
  function(species_list) {
    
    # Getting template rasters
    out_raster <- country_rast_cropped
    out_raster_richness <- country_rast_cropped
    if(outcome %in% 'pop_density'){out_raster_biomass <- country_rast_cropped}
    # Looping across list of species
    for(ss in species_list) {
      
      ### First for habitat availability
      # Getting raster
      species_raster <- raster(ss)
    tmp_species_name <-
	    ss %>%
	    gsub(paste0('.*',outcome,'_'),'',.) %>%
	    gsub(paste0('_',thresh,'.*'),'',.)
     #cat(tmp_species_name,'\n') 
      # Checking whether we need to loop over the species
      if(tryCatch(!is.null(raster::crop(species_raster,extent(country_rast_cropped))), error=function(e) return(FALSE))) { # This returns a true/false statement, TRUE if rasters intersect, FALSE if they do not
        
        # Cropping raster
        species_raster <- raster::crop(species_raster,extent(country_rast_cropped))
        
        # Getting species richness
        species_richness <- !is.na(species_raster)
        
        # Removing NAs
        species_raster[is.na(species_raster)] <- 0
        species_raster[species_raster<0] <- 0
        
        # Extending both rasters to extent of south america
        species_raster <- raster::extend(species_raster, country_rast_cropped, value = 0)
        species_richness <- raster::extend(species_richness, country_rast_cropped, value = 0)
        
        # Adding species raster to the out rasters
        out_raster <- out_raster + species_raster
        out_raster_richness <- out_raster_richness + species_richness

	if(outcome %in% 'pop_density') {
		mass_species <- body_mass_merge$est_mass_kg[body_mass_merge$species_merge %in% tmp_species_name]
		if(length(mass_species) >= 1) {
			if(is.na(mass_species)) {mass_species <- 0}
			if(is.null(mass_species)) {mass_species <- 0}
			   biomass_raster <- species_raster * mass_species
			   biomass_raster <- raster::extend(biomass_raster,country_rast_cropped,value=0)
			   out_raster_biomass <- out_raster_biomass + biomass_raster
		} # End if statement for no body mass data
	} # End if statement for density
        
        } # End of try catch if statement
     } # End of species loop
    
    # Making list to return object
    if(outcome %in% 'pop_density') {
	    out_list <- 
		    list(out_raster,
			 out_raster_richness,
			 out_raster_biomass)
    
    names(out_list) <- c('outcome','richness','biomass')
    } else {
	    out_list <-
                    list(out_raster,
                         out_raster_richness)

    names(out_list) <- c('outcome','richness')
    
    }
    
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
  #file_list_taxa <-
  #  list.files(paste0(getwd(),'/Outputs/Raster_Outputs/BAU/',taxa),
  #             full.names = TRUE)

	cmd <-
		paste0("find /data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/Outputs/Raster_Outputs/BAU/",taxa,"/ -type f")
 
      file_list_taxa <-  
	system(
  cmd,
  intern = TRUE
)



  for(thresh in thresh_list[1]) { # Looping through taxa
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
            .[grepl(outcome,.)] %>%
	    .[!grepl('climate_suitable',.,ignore.case=TRUE)]
          
          # Chunking file list into n parts for parallelisation
          chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
          file_list_chunked <- chunk2(file_list_outcome, n = mc_cores)
          # Running function
          # This returns a list of rasters
          out_raster_list <- mclapply(file_list_chunked, species_funcion, mc.cores = mc_cores)
          # out_raster_list <- lapply(file_list_chunked[1], species_funcion) # For testing, outside of parallel
          
          # Getting list of rasters we need...
          out_raster_stack <- country_rast_cropped
          out_raster_richness <- country_rast_cropped
          out_raster_biomass <- country_rast_cropped
          # Loop through the list and add the 'outcome' and 'richness' layers
          for(ii in 1:length(out_raster_list)) {out_raster_stack <- out_raster_stack + out_raster_list[[ii]][['outcome']]}
          for(ii in 1:length(out_raster_list)) {out_raster_richness <- out_raster_richness + out_raster_list[[ii]][['richness']]}
          if(outcome %in% 'pop_density') {
		  for(ii in 1:length(out_raster_list)) {out_raster_biomass <- out_raster_biomass + out_raster_list[[ii]][['biomass']]}

	  
	  }
          # Saving rasters
          # writeRaster(out_raster_stack,'/data/ouce-glob2loc/pubh0329/outcome_check.tif')
          # writeRaster(out_raster_richness,'/data/ouce-glob2loc/pubh0329/richness_check.tif')
          writeRaster(out_raster_stack,
                      paste0(getwd(),'/Analyses/Landscape_Analyses/',country,'_',
                            thresh,'_',taxa,'_',year,'_',scen,'_',outcome,'_',
                            Sys.Date(),'.tif'),
                      overwrite = TRUE)

          writeRaster(out_raster_richness,
                      paste0(getwd(),'/Analyses/Landscape_Analyses/',country,'_',
                            thresh,'_',taxa,'_',year,'_',scen,'_',outcome,'_richness_',
                            Sys.Date(),'.tif'),
                      overwrite = TRUE)

          if(outcome %in% 'pop_density') {
		  cat('Min Value:',minValue(out_raster_biomass),'\n')
		  cat('Max Value:',maxValue(out_raster_biomass),'\n')
		   writeRaster(out_raster_biomass,
                      paste0(getwd(),'/Analyses/Landscape_Analyses/',country,'_',
                            thresh,'_',taxa,'_',year,'_',scen,'_biomass_',
                            Sys.Date(),'.tif'),
                      overwrite = TRUE)
	  rm(out_raster_biomass)
	  }


          # Clearing files
          rm(out_raster_richness,out_raster_stack,file_list_chunked,file_list_outcome)
          
        } # End loop for outcome type - e.g. habitat area or population abundance
      } # End loop through stress scenarios
    } # End loop through years
  } # End loop through thresholds
} # End big loop through taxa
