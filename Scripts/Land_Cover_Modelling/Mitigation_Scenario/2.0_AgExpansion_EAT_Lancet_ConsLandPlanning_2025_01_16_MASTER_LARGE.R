#!/usr/bin/env Rscript

# Read me ----
# Runs the conversion of cropland and pastureland for Business-as-usual
# To run for alternative scenarios, alter the relevant code on L367 and for saving the 
# rasters (L645 onwards)

# Packages ----
library(dplyr)
library(scales)
library(raster)
library(countrycode)
library(data.table)
library(matrixStats)
library(doMC)
library(future.apply)
library(stringr)
library(terra)

registerDoMC(cores=9)
rasterOptions(maxmemory = 1e+11)
options(future.globals.maxSize= 32*1024^3)

# Clean up ----
rm(list=ls())
# Increasing size of exported data needed to run this in parallel via the futures pacckage
# Default max size is 500MB, this sets it to 3GB
# options(future.globals.maxSize = 3000 * 1024^2)
# Load models, functions and data -----
setwd("/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity")
# setwd("/Users/maclark/Desktop/Multiple Stresses of Biodiversity")
source("/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/Scripts/Land Cover Modelling/SSP1_Land_Cover_Modelling/0.0_LandUseChangeFunctions_ConsLandPlanning_2025_01_14.R")

# Is conservation land use planning implemented?
conservation_land_planning = 'yes'

# List for GBM outputs
file.list.mods_crop <- list.files(pattern = "Crop", path = paste0(getwd(),"/Land Forecasting Coefs and Plots/Model Coefficients"), full.names = TRUE) %>% .[grep('July2021',.)]
file.list.mods_pasture <- list.files(pattern = "Pasture", path = paste0(getwd(),"/Land Forecasting Coefs and Plots/Model Coefficients"), full.names = TRUE) %>% .[grep('July2021',.)]

# List for Regional data
file.list.data <- list.files(path = paste0(getwd(),"/Forecasting_Data"), pattern = '29July2021')
targets <- read.csv(paste0(getwd(),"/Land Forecasting Targets/Land_Targets_2025_01_22.csv"),
                    header = TRUE,
                    stringsAsFactors = FALSE) 
region_lookup <- read.csv(paste0(getwd(),"/Other Data Inputs/RegionLookUp_IMPACT_v2.csv"),
                          header = TRUE,
                          stringsAsFactors = FALSE)
iso3_lookup <- read.csv(paste0(getwd(),"/Other Data Inputs/ISO_N_to_Impact_N_translation_v2.csv"),
                        header = TRUE,
                        stringsAsFactors = FALSE)

# Getting region list to assign names to the model object lists -----
region.names <- 
  gsub(" Crop Model.*","",file.list.mods_crop) %>%
  gsub(".*Coefficients/","",.) %>%
  trimws('both') %>%
  unique(.)

# Creating directory for outputs
if('EAT_Lancet' %in% list.files(paste0(getwd(),'/Land Forecast Outputs'))) {
  # Nothing
} else {
  # Create directory
  dir.create(paste0(getwd(),'/Land Forecast Outputs/EAT_Lancet'))
  dir.create(paste0(getwd(),'/Land Forecast Outputs/EAT_Lancet/tifs'))
}

# FSwitch this depending on climate scenario you're using to forecast
# This can either be '45' for RCP 45 or '85' for RCP85
gdd.cutoff <- '45'
# Getting only small countries
# To run on 25 cores
row.index.rasters <- list.files(path = paste0(getwd(),"/Other Data Inputs/Row_Index_New"),
                                pattern = ".RData", full.names = TRUE)
row.index.bmk <- load(row.index.rasters[grep('JPN',row.index.rasters)])
# The above line loads an object named tmp.country
row.index.bmk <- tmp.country
# List of countries to keep
# These are relatively small countries
countries.keep <- c()
for(i in row.index.rasters) {
  load(i)
  tot_cells <- tmp.country@nrows * tmp.country @ncols
  if(tot_cells <= 1000000) {
    # Nothing, country too big
  } else {
    countries.keep <- c(countries.keep, gsub('.*/Row_Index_New/',"",i) %>% substr(.,1,3))
  }
}

# Getting rasters for conservation land use planning
if(conservation_land_planning %in% 'yes') {
  contraction_areas <-
    c(paste0(getwd(),'/Scenario_Analysis_Inputs/Land_Cover_Layers/Conservation_Areas_SSP1/Biodiversity_Areas.tif'),
      paste0(getwd(),'/Scenario_Analysis_Inputs/Land_Cover_Layers/Conservation_Areas_SSP1/Wilderness_Areas.tif'))
  
  no_expansion_areas <-
    c(paste0(getwd(),'/Scenario_Analysis_Inputs/Land_Cover_Layers/Conservation_Areas_SSP1/Climate_Areas.tif'),
      paste0(getwd(),'/Scenario_Analysis_Inputs/Land_Cover_Layers/Conservation_Areas_SSP1/Climate_Corridors.tif'))
  
  cons_areas <-
    cons_area_function(contraction_areas = contraction_areas,
                       no_expansion_areas = no_expansion_areas)
}

# Starting loop -----
t1 <- Sys.time()
# for(k in c(1:5,9:12,8,7,6)) { # Odd ordering is to complete the small regions first to check the script
for(k in sample(1:12,12)) { # Odd ordering is to complete the small regions first to check the script
#for(k in which(grepl('South America',region.names))) {
  #for(k in 1) {
  cat(region.names[k])	
  # Check if we need to run the region -----
  # Creating iso3c list to import specific row index rasters for each nation -----
  iso3c <- unique(region_lookup[region_lookup$IUCN_region %in% region.names[k], "ISO3_IMPACT"])
  print(paste(iso3c, "in region"))
  # Limiting iso3c to only countries for which we have forecasts, 
  # excluding those that we've already run and converting back into a list
  completed_countries <- list.files(path = paste0(getwd(),"/Land Forecast Outputs/EAT_Lancet/tifs"))
  completed_countries <- c(unique(substring(completed_countries, 1, 3)))
  if(length(iso3c[iso3c %in% completed_countries]) >0) {
    cat(paste(iso3c[iso3c %in% completed_countries], "completed"))
  } else {
    cat("none completed")
  }
  
  targets <- read.csv(paste0(getwd(),"/Land Forecasting Targets/Land_Targets_2025_01_22.csv"),
                      header = TRUE,
                      stringsAsFactors = FALSE) 
  
  targets <- targets %>% 
    filter(!is.na(prop_change_2010_2015))
  iso3c <- iso3c[iso3c %in% targets$ISO3]
  iso3c <- iso3c[!(iso3c %in% completed_countries)]
  iso3c <- as.character(iso3c)
  
  # Limiting country selection to only those that have cropland change in at least one time interval -----
  targets.tmp = targets
  # Getting targets for scenario we're forecasting
  targets.tmp = targets.tmp[!is.na(targets.tmp$ISO3) & targets.tmp$scenario == "FLX",]
  
  targets.tmp$Count <-  0
  
  # Getting only small countries to run
  iso3c <- iso3c[iso3c %in% countries.keep]
  # iso3c <- iso3c[!(iso3c %in% c('LBY','LSO','OAO','SLB','OIO','BDI','ISL','VUT','SLB','FJI','DJI','COD'))]
  
  iso3c <- iso3c[!(iso3c %in% c('OAO','OIO','LSO'))]
  # Getting countries with no change in crop extent in any time period
  # Doing this to not loop through countries that don't have any change in cropland
  # Looping through time period
  
  if(length(iso3c) == 0){
    cat(paste("no countries needed for ", region.names[k]))
  } else {
    cat(paste0("Beginning ", region.names[k]))
    # Loading specifc regional data set -----
    # print("Loading regional data")
    
    region.data <- fread(paste0(getwd(),"/Forecasting_Data/",
                                file.list.data[k]),
                         stringsAsFactors = FALSE) %>%
      mutate(country_id = country_id_impact)
    region.data <- as.data.frame(region.data)
    
    
    # Load targets -----
    targets <- read.csv(paste0(getwd(),"/Land Forecasting Targets/Land_Targets_2025_01_22.csv"),
                        header = TRUE,
                        stringsAsFactors = FALSE) 
    # Loading specific model output
    # First for crops
    load(file.list.mods_crop[k])
    
    coefs.crop.increase <- coef.increase.mean.se
    coefs.crop.decrease <- coef.decrease.mean.se
    coefs.crop.prob <- coef.multinom.mean
    
    # And repeating for pastures
    load(file.list.mods_pasture[k])
    
    coefs.pasture.increase <- coef.increase.mean.se
    coefs.pasture.decrease <- coef.decrease.mean.se
    coefs.pasture.prob <- coef.multinom.mean
    
    # Manipulating regional data to update it ----
    cat("Manipulating regional data")
    # #eliminating row.name column by name, in case order changes
    region.data <-
      region.data %>%
      dplyr::select(-Pasture2000) %>%
      dplyr::select(-Crop1995) %>%
      dplyr::select(-Crop2000) %>%
      dplyr::select(-delta_crop_1995_2000) %>%
      dplyr::select(-delta_crop_2000_2005) %>%
      dplyr::select(-delta_pasture_1995_2000) %>%
      dplyr::select(-delta_pasture_2000_2005) %>%
      dplyr::select(-urban_2010)
    
    # Emulating layout of data_w_preds_crop
    region.data$prop_to_predict_from_cell_delta_crop <- region.data$delta_crop_2005_2010
    region.data$prop_to_predict_from_cell_delta_pasture <- region.data$delta_pasture_2005_2010
    
    # Updating value of IUCN.Region.Name
    region.data$IUCN.Region.Name <- region.names[k]
    # Dropping cells for which Iucn.Region is NA
    region.data$IUCN.Region.Name[is.na(region.data$IUCN.Region)] <- NA
    
    # Adding in ISO3 codes to emulate layout of data_w_preds_crop
    # (data set used to build models for land cover change)
    region.data$ISO3 <- countrycode(region.data$country_id, "iso3n","iso3c")
    iso3_lookup_tmp <- iso3_lookup[iso3_lookup$region == region.names[k],]
    region.data[is.na(region.data$ISO3), "ISO3"] <- 
      iso3_lookup_tmp$IMPACT.Code[match(region.data[is.na(region.data$ISO3), "country_id"],
                                        iso3_lookup_tmp$ISO_numeric_new)]
    
    
    # Renaming columns to match col names in script below
    region.data <- dplyr::rename(region.data,
                                 ag_suitability = suitability,
                                 dist_50k = distance_cities, # Travel time to nearest city despite name. Name is holdover from older analyses
                                 PA_binary = wdpa,
                                 prop_to_predict_from_crop = Crop2010,
                                 prop_to_predict_from_pasture = Pasture2010,
                                 prop_adj_to_predict_from_crop = adj_Crop2010,
                                 prop_adj_to_predict_from_pasture = adj_Pasture2010,
                                 IUCN.Region.Number = iucn_region,
                                 Prop.NonArable.2010 = non_rable_2010)
    
    
    
    # Eliminating rounding errors in MODIS data
    # Some cells have total land cover > 1
    # For these cells, dividing proportion of each land cover type by sum of land cover
    region.data$tmp <- region.data$prop_to_predict_from_crop +
      region.data$prop_to_predict_from_pasture + 
      region.data$Prop.NonArable.2010 +
      region.data$UrbanExtent2010
    # Getting index of cells that have land cover > 1
    tmp_index <- which(round(region.data$tmp, digits = 2) > 1 & !is.na(region.data$tmp))
    # Counter for number of digits to round to
    tmp_counter = 8
    # putting in a while loop
    while(length(tmp_index) > 1 & tmp_counter >= 2) {
      # And dividing land cover of each type by sum of all land cover types for these cells
      region.data$prop_to_predict_from_crop[tmp_index] <- round(region.data$prop_to_predict_from_crop[tmp_index] / region.data$tmp[tmp_index], digits = tmp_counter)
      region.data$prop_to_predict_from_pasture[tmp_index] <- round(region.data$prop_to_predict_from_pasture[tmp_index] / region.data$tmp[tmp_index], digits = tmp_counter)
      region.data$Prop.NonArable.2010[tmp_index] <- round(region.data$Prop.NonArable.2010[tmp_index] / region.data$tmp[tmp_index], digits = tmp_counter)
      region.data$Prop.Urban.2010[tmp_index] <- round(region.data$Prop.Urban.2010[tmp_index] / region.data$tmp[tmp_index], digits = tmp_counter)
      # And updating region.data$tmp
      region.data$tmp <- region.data$prop_to_predict_from_crop +
        region.data$prop_to_predict_from_pasture + 
        region.data$Prop.NonArable.2010 +
        region.data$Prop.Urban.2010
      # and updating tmp_index
      tmp_index <- which(region.data$tmp > 1 & !is.na(region.data$tmp))
      # Updating counter
      tmp_counter = tmp_counter - 1
    }
    
    # Removing tmp_counter and tmp_index
    rm(tmp_counter, tmp_index)
    
    # And removing column tmp from region.data
    region.data <- dplyr::select(region.data, -tmp) %>% mutate()
    
    # And recreating values for prop.crop.2010 and prop.pasture.2010
    region.data$Prop.Pasture.2010 <- region.data$prop_to_predict_from_pasture
    region.data$Prop.Crop.2010 <- region.data$prop_to_predict_from_crop
    
    for(i in 3:10) {
      # Creating vector of target year
      counter.vector = targets.tmp[,i]
      # Getting rows that don't change
      rows.nochange = which(counter.vector == 0)
      # Adding one to counter
      targets.tmp$Count[rows.nochange] <- targets.tmp$Count[rows.nochange] + 1
    }
    
    # Limiting to only countries that change cropland extent in any time period
    targets.tmp <- targets.tmp[targets.tmp$Count != 10,]
    
    # Drop countries we don't need
    iso3c <- iso3c[iso3c %in% targets.tmp$ISO3]
    
    # And dropping countries that do not run simultaneously on 25 cores
    # Doing this to save wall time when using ARC
    # These will be run in a later script
    # iso3c <- iso3c[!(iso3c %in% c('LBY','OAO','SLB','OIO','BDI','ISL','VUT','SLB','FJI','DJI','COD'))]
    iso3c <- iso3c[!(iso3c %in% c('OAO','OIO','LSO'))]
    
    # Looping with specific countries -----
    for(c in sample(1:length(iso3c),length(iso3c))) {
    #for(c in which(iso3c %in% 'ECU')) {
      cat(paste("Beginning", iso3c[c]))
      # Targets 
      targets <- read.csv(paste0(getwd(),"/Land Forecasting Targets/Land_Targets_2025_01_22.csv"),
                          header = TRUE,
                          stringsAsFactors = FALSE)
      # Getting name of raster to import
      row.index.name <- list.files(path = paste0(getwd(),"/Other Data Inputs/Row_Index_New/"),
                                   pattern = ".RData")
      row.index.name <- row.index.name[grep(iso3c[c],row.index.name)]
      # Importing country raster with row.index values
      # row.index.raster <- raster(paste0("/Users/maclark/Desktop/Row_Index_New//",
      #                                  row.index.name))
      
      # Loading row index raster
      load(paste0(getwd(),"/Other Data Inputs/Row_Index_New/",
                  row.index.name))
      row.index.raster <- tmp.country
      
      # Creating row.index.frame that contains only values of row.index associated with the country-specific raster
      row.index.frame <- data.frame(Row.Index = getValues(row.index.raster))
      
      # Merging with row.index.frame
      # Limits data frame to only country-level cells
      model_data_all <- left_join(row.index.frame,
                                  region.data %>% dplyr::rename(Row.Index = row_index))
      
      # All non PA areas should have value of one. 
      # This apparently isn't the case for all countries, so inserting manually here
      model_data_all$PA_binary[-which(model_data_all$PA_binary == 1)] <- 0
      # Likewise for GDD data
      # model_data_all$gdd_cutoff <- 
      #   model_data_all[,paste0('gdd_cutoff_',gdd.cutoff)]
      # model_data_all$gdd_cutoff[-which(model_data_all$gdd_cutoff >= 2010 & !is.na(model_data_all$gdd_cutoff_45))] <- 0
      
      # Dropping all data not from the country
      row.nums <- which(model_data_all$ISO3 %in% iso3c[c])
      
      model_data_all[-row.nums,] <- NA
      
      # Removing row.index.frame
      rm(row.index.frame, row.nums)
      
      # 1) Put data into one dataframe and rename if necessary ------
      # Changes factor variables into character strings
      model_data_all <- model_data_all %>%
        mutate_if(is.factor, as.character)
      
      # Rename some columns so that they work with the script (easier to change
      #     names than the script)
      model_data_all <- model_data_all %>%
        dplyr::rename(IUCN_region = IUCN.Region.Name,
                      row_index = Row.Index, 
                      ISO_numeric = country_id,
                      prop_after_change_crop = Prop.Crop.2010,
                      prop_after_change_pasture = Prop.Pasture.2010) %>%
        dplyr::select(row_index, ISO_numeric, IUCN.Region.Number:prop_after_change_crop) %>%
        as.data.frame()
      
      # Order by row_index and rename the rows - this is for checking things later
      model_data_all <- model_data_all[order(model_data_all$row_index),]
      row.names(model_data_all) <- 1:nrow(model_data_all)
      
      
      # Data prep: Set constants: -----
      # Predictors
      predictor_vars_crop <-
        coefs.crop.increase$variable %>%
        .[-grep('Intercept',.)] %>%
        .[-grep('ISO',.)] %>%
        as.character(.) %>%
        c(.,'ISO_numeric')
      
      
      predictor_vars_pasture <- 
        coefs.pasture.increase$variable %>%
        .[-grep('Intercept',.)] %>%
        .[-grep('ISO',.)] %>%
        as.character(.) %>%
        c(.,'ISO_numeric')
      
      # Set scenarios to investigate
      unique(targets$scenario)
      targets <- subset(targets, scenario == "FLX")
      
      # Dimensions are in metres, so convert to square km
      cell_area <- res(row.index.raster)[1]/1000 * res(row.index.raster)[2]/1000
      rows <- nrow(row.index.raster)
      cols <- ncol(row.index.raster)
      
      # Convolution filter  
      adj_cell_mat <- matrix(nrow=3, ncol = 3, data = 1)
      
      # Data prep: Convert data to NAs from countries that we are not forecasting -----
      # We are forecasting country by country, so this removes any data associated with any country that is not the target country
      # Getting list of row numbers associated with the specific target country
      row.nums <- grep(iso3c[c], model_data_all$ISO3)
      # Converting all other rows to NAs
      model_data_all[-row.nums,-1] <- NA
      
      # Data prep: Add area of each cell ----
      model_data_all$cell_area <- cell_area
      # Put NAs in, so it throws an error if necessary
      model_data_all[is.na(model_data_all$prop_to_predict_from_crop), "cell_area"] <- NA
      rm(row.nums)
      
      ### merging in GDD binary information
      # file list
      gdd_files <- list.files(paste0(getwd() %>% gsub('pubh-glob2loc','ouce-glob2loc',.),'/CMIP6_Climate_Data/SSP1-2.6/Interpolated_Rasters'), pattern = 'GDD_5C_binary', full.names = TRUE)
      
      # years
      years <- sort(str_extract(gdd_files,'binary_[0-9]{4,4}')) %>% gsub('binary_','',.)
      
      # making data frame
      for(gdd_y in years[years <= 2050]) {
        
        # getting binary raster
        tmp_raster <- raster(gdd_files[grepl(gdd_y,gdd_files)]) %>% crop(.,row.index.raster)
        
        if(gdd_y %in% '2010') { # for 2010
          # making data frame
          gdd_frame <- 
            data.frame(row_index = getValues(row.index.raster),
                       gdd_binary = getValues(tmp_raster))
          
          # updating name
          names(gdd_frame)[ncol(gdd_frame)] <- paste0('GDD_Binary_',gdd_y)
          
        } else { # for other years
          gdd_frame[,paste0('GDD_Binary_',gdd_y)] <- getValues(tmp_raster)
        } # end else statement
        
        # and removing raster
        rm(tmp_raster)
      } # end making gdd frame
      
      # merging gdd frame in
      model_data_all <-
        left_join(model_data_all,
                  gdd_frame)
      
      # Merging in conservation areas, if needed
      if(conservation_land_planning %in% 'yes') {
        cons_areas_df <-
          data.frame(row_index = getValues(row.index.raster),
                     no_expansion_areas = getValues(raster::crop(cons_areas$no_expansion_areas_raster,row.index.raster)),
                     contraction_areas = getValues(raster::crop(cons_areas$contraction_areas_raster,row.index.raster))) %>%
          mutate(cons_planning_areas = no_expansion_areas + contraction_areas)
        
        # Merging in
        model_data_all <-
          left_join(model_data_all,
                    cons_areas_df)
        
        # Clearing to remove area
        rm(cons_areas_df)
      }
      
      # converting to GDD raster to NAs
      
      # Data prep: Update the proportion of ag land with new values -----
      # Update "prop_to_predict_from" to current values
      # Moves next time step data into current time step data to get predictions
      model_data_current <- model_data_all
      model_data_current$prop_to_predict_from_crop <- model_data_current$prop_after_change_crop
      model_data_current$prop_to_predict_from_pasture <- model_data_current$prop_after_change_pasture
      # Add in a reference column for proportion of ag at the start of the loop
      model_data_current$prop_start_loop_crop <- model_data_current$prop_to_predict_from_crop
      model_data_current$prop_start_loop_pasture <- model_data_current$prop_to_predict_from_pasture
      
      # Updating columns to reflect code in lines 408 - 412
      model_data_current$prop_change_cell_crop <- model_data_current$prop_to_predict_from_cell_delta_crop
      model_data_current$prop_change_cell_pasture <- model_data_current$prop_to_predict_from_cell_delta_pasture
      
      # Some variables go to NAs, this messes things up. Sets everything to NA if a lake or something equivalent
      # Cut out lakes and other areas where not all predictors are present 
      model_data_current[is.na(model_data_current$prop_to_predict_from_crop),
                         c('dist_50k', 'ag_suitability', 'PA_binary',
                           'prop_adj_to_predict_from_crop', 'prop_adj_to_predict_from_pasture',
                           'prop_to_predict_from_crop', 'prop_to_predict_from_pasture',
                           'prop_to_predict_from_cell_delta_crop', 'prop_to_predict_from_cell_delta_pasture',
                           'ISO_numeric',
                           "prop_change_cell_crop", #"prop_change_cell_pasture",
                           "prop_after_change_crop", "prop_after_change_pasture",
                           "cell_area")] <- NA
      
      # Need to run convolution again to update area of cropland in each cell, but have not adjusted adjacent area
      # Run convolution to get new surrounding values for cropland...
      # Should not need to do this
      model_data_current <- adj.ag.f(m = model_data_current, 
                                     r = rows, 
                                     c = cols, 
                                     crop = "crop")
      
      # ...and pasture
      model_data_current <- adj.ag.f(m = model_data_current, 
                                     r = rows, 
                                     c = cols, 
                                     crop = "pasture")
      
      # Data prep: Calculate targets based on MODIS data ----
      targets <- model_data_current %>%
        dplyr::select(ISO_numeric, ISO3, prop_start_loop_crop, cell_area) %>%
        mutate(area_cropland_km = prop_start_loop_crop * cell_area) %>%
        dplyr::group_by(ISO_numeric) %>%
        plyr::summarise(area_cropland_km_2010 = sum(area_cropland_km[!is.na(area_cropland_km)])) %>%
        # Add ISO3, so we can bind to targets
        mutate(ISO3 = as.character(iso3c[c])) %>%
        # Bind to proportional increases in cropland
        left_join(., targets %>%
                    dplyr::select(-grep("area", names(.)))) %>%
        as.data.frame()
      
      # Calculate area increases
      # Do first area separately and then loop
      # Multiplies proportional increase by current area in cropland
      targets$target_2010_2015 <- targets$area_cropland_km_2010 * targets$prop_change_2010_2015
      targets$area_cropland_km_2015 <- targets$area_cropland_km_2010 + targets$target_2010_2015
      
      # Loop calculates future targets of cropland by 5 year interval
      for(i in 2:length(grep("prop_change", names(targets)))){
        # Get the interval and variable names
        interval <- names(targets)[grep("prop_change", names(targets))][i]
        n_target <- gsub("prop_change", "target", interval)
        n_area <- paste("area_cropland_km", sub("^.*_", "", n_target), sep = "_")
        n_previous <- paste("area_cropland_km", 
                            gsub("(^.+_)(.*)(_.+$)", "\\2", interval), 
                            sep = "_")
        targets[[n_target]] <- targets[,names(targets)[grep("prop_change", names(targets))][i]] * 
          targets[[n_previous]]
        
        targets[[n_area]] <- targets[[n_target]] + targets[[n_previous]]
      }

      print(targets)
      
      # Data prep: Add in proportional land demand ----- 
      # Adds in targets to the big data frame
      model_data_current <- left_join(model_data_current, 
                                      targets[,-grep("area_cropland", names(targets))])
      
      # At this point we have a data frame with all of projection variables updated to 2010
      # Also contains amount of urban land
      # Also contains forecasts of future crop demand for 5 year intervals
      
      # Get new probabilities for a) change b) amount of change -----
      # Cut regions without predictions
      na.data <- model_data_current[is.na(model_data_current$ISO_numeric),]
      model_data_current <- model_data_current[!is.na(model_data_current$ISO_numeric),]
      
      # Add in new demand
      # Updates predictor variable for model forecasts
      # Changed the next line because Dave told me to
      model_data_current$change_demand_crop <- model_data_current$prop_change_2010_2015 + 1
      
      # Creating column named urban_prop to mesh with the functions script
      # Might be able to move this further up when renaming columns
      model_data_current$urban_prop <- model_data_current$UrbanExtent2010
      
      # Flip this if forecasting in climate rcp45 or rcp85
      model_data_current$gdd_cutoff <- 
        model_data_current[,paste0('gdd_cutoff_',gdd.cutoff)]
      
      # And updating GDD for 2010
      model_data_current$GDD_binary <- model_data_current$GDD_Binary_2010
      
      # a) Update cropland
      # Calls from pred_fun in 0_Functions
      model_data_current <- model_data_current %>%
        do(pred_fun(p = ., crop = "crop")) %>% 
        as.data.frame()
      
      # model_data_current <- model_data_current %>%
      #   do(pred_fun(p = ., crop = "pasture")) %>% 
      #   as.data.frame()
      
      # Add columns for pasture loss - currently not doing this -----
      # What this will do is take model data, and calculate how much space is available for pasture
      model_data_current$prop_for_pasture <- round(1 - rowSums(model_data_current[,c('urban_prop',
                                                                                     'prop_to_predict_from_crop',
                                                                                     'Prop.NonArable.2010')],
                                                               na.rm = TRUE),
                                                   digits = 7)
      # Not dealing with pasture outside of country borders...
      model_data_current$prop_for_pasture[is.na(model_data_current$ISO_numeric)] <- NA
      # And really making sure it isn't negative
      # Although the rounding above should take care of that
      model_data_current$prop_for_pasture[model_data_current$prop_for_pasture < 0] <- 0
      
      # Check that nothing mad is being predicted -----
      model_data_current$tmp <- 
        model_data_current$prop_to_predict_from_crop + 
        model_data_current$urban_prop + 
        # model_data_current$prop_to_predict_from_pasture + 
        model_data_current$pred_change_cell_prop_crop_increase
      # range(model_data_current$tmp, na.rm = TRUE)
      
      model_data_current$tmp <- 
        model_data_current$prop_to_predict_from_crop + 
        # model_data_current$prop_to_predict_from_pasture + 
        model_data_current$pred_change_cell_prop_crop_decrease
      # range(model_data_current$tmp, na.rm = TRUE)
      
      model_data_current <- dplyr::select(model_data_current, -tmp)
      
      # c) Update areas without data
      # Adding columns to the NA data set to mesh with the model_data_current data
      na.data$pred_no_change_crop <- NA
      na.data$pred_to_change_crop <- NA
      na.data$area_converted_crop <- NA
      na.data$pred_no_change_pasture <- NA
      na.data$pred_to_change_pasture <- NA
      na.data$area_converted_pasture <- NA
      
      # na.data$ISO_numeric <- as.numeric(na.data$ISO_numeric)
      
      # d) Bind data back together
      model_data_current <- bind_rows(model_data_current, na.data)
      # model_data_current2 <- data.table::rbindlist(model_data_current,na.data)
      
      model_data_current <- model_data_current[order(model_data_current$row_index),]
      rm(na.data)
      
      # Limiting columns that not needed for the forecast script
      model_data_current <- dplyr::select(model_data_current,
                                          row_index,
                                          ISO3,
                                          ISO_numeric,
                                          prop_to_predict_from_crop,
                                          prop_to_predict_from_pasture,
                                          Prop.NonArable.2010,
                                          urban_prop,
                                          prop_adj_to_predict_from_crop,
                                          prop_adj_to_predict_from_pasture,
                                          change_demand_crop,
                                          dist_50k,
                                          ag_suitability,
                                          PA_binary,
                                          GDD_binary,
                                          prop_to_predict_from_cell_delta_crop,
                                          prop_to_predict_from_cell_delta_pasture,
                                          pred_to_increase_crop,
                                          pred_to_decrease_crop,
                                          pred_change_cell_prop_crop_decrease,
                                          pred_change_cell_prop_crop_increase,
                                          area_converted_crop_decrease,
                                          area_converted_crop_increase,
                                          prop_for_pasture,
                                          prop_change_2010_2015,
                                          prop_change_2015_2020,
                                          prop_change_2020_2025,
                                          prop_change_2025_2030,
                                          prop_change_2030_2035,
                                          prop_change_2035_2040,
                                          prop_change_2040_2045,
                                          prop_change_2045_2050,
                                          target_2010_2015,
                                          target_2015_2020,
                                          target_2020_2025,
                                          target_2025_2030,
                                          target_2030_2035,
                                          target_2035_2040,
                                          target_2040_2045,
                                          target_2045_2050,
                                          gdd_cutoff,
                                          UrbanExtent2015,
                                          UrbanExtent2020,
                                          UrbanExtent2025,
                                          UrbanExtent2030,
                                          UrbanExtent2035,
                                          UrbanExtent2040,
                                          UrbanExtent2045,
                                          UrbanExtent2050,
                                          GDD_Binary_2010,
                                          GDD_Binary_2015,
                                          GDD_Binary_2020,
                                          GDD_Binary_2025,
                                          GDD_Binary_2030,
                                          GDD_Binary_2035,
                                          GDD_Binary_2040,
                                          GDD_Binary_2045,
                                          GDD_Binary_2050,
                                          other_2010,
                                          no_expansion_areas,
                                          contraction_areas,
                                          cons_planning_areas)
      
      # Run the future conversion ----
      # Run the conversion
      rm(model_data_all)
      # Run this if using parallel
      # test.list <- c()
      # # Making exceptions for geographically large countries
      # # Need to do this becasue some countries will take all of the RAM available on the computer if doing 5 iterations
      # # And this won't work...as in the computer will turn off. 
      # # Obviously adjust this depending on the memory available. Because R is so inefficient at raster data, this requires
      # #    a LOT of memory (>200GB of working memory), so be warned!
      # if(iso3c[c] == 'IND' | iso3c[c] == 'BRA' | iso3c[c] == "FJI") {
      #   # First set of forecasts
      #   test.list.tmp <- foreach(i = 1:13) %dopar% year_fun(data.frame(input = "model_runs"))
      #   test.list <- c(test.list, test.list.tmp)
      #   rm(test.list.tmp)
      #   # Second set of forecasts
      #   test.list.tmp <- foreach(i = 1:12) %dopar% year_fun(data.frame(input = "model_runs"))
      #   test.list <- c(test.list, test.list.tmp)
      #   rm(test.list.tmp)
      # } else if(iso3c[c] == 'AUS' | iso3c[c] == 'CHM' | iso3c[c] == 'CAN') {
      #   # Looping through first nine
      #   # For total of twenty four model runs
      #   test.list.tmp <- foreach(i = 1:9) %dopar% year_fun(data.frame(input = "model_runs"))
      #   test.list <- c(test.list, test.list.tmp)
      #   rm(test.list.tmp)
      #   # Second set of forecasts
      #   # Looping through second nine
      #   # For total of 18 model runs
      #   test.list.tmp <- foreach(i = 1:9) %dopar% year_fun(data.frame(input = "model_runs"))
      #   test.list <- c(test.list, test.list.tmp)
      #   rm(test.list.tmp)
      #   # Looping through final seven
      #   # For total of twenty five model runs
      #   test.list.tmp <- foreach(i = 1:7) %dopar% year_fun(data.frame(input = "model_runs"))
      #   test.list <- c(test.list, test.list.tmp)
      #   rm(test.list.tmp)
      #   
      # } else if(iso3c[c] == 'RUS' | iso3c[c] == 'USA')  {
      #   print(iso3c[c])
      #   # Looping through three lots of seven and one four
      #   for(n.iter in 1:3) {
      #     print(n.iter)
      #     test.list.tmp <-  foreach(i = 1:7) %dopar% year_fun(data.frame(input = "model_runs"))
      #     # And appending forecasts to a list
      #     test.list <- c(test.list, test.list.tmp)
      #     rm(test.list.tmp)
      #   }
      #   test.list.tmp <- foreach(i = 1:4) %dopar% year_fun(data.frame(input = "model_runs"))
      #   test.list <- c(test.list, test.list.tmp)
      #   rm(test.list.tmp)
      # } else {
      #   test.list.tmp <-  foreach(i = 1:5) %dopar% year_fun(data.frame(input = "model_runs"))
      #   # And appending forecasts to a list
      #   test.list <- c(test.list, test.list.tmp)
      #   rm(test.list.tmp)
      # }    
      
      test.list <- c()
      
      # making function to return rasters contain prop to increase estimates
      
      
      # Making exceptions for geographically large countries
      # Need to do this becasue some countries will take all of the RAM available on the computer if doing 5 iterations
      # And this won't work...as in the computer will turn off
      for(n.iter in 1:25) {
      #for(n.iter in 1:1) { # For Testing...
        print(c(iso3c[c], n.iter))
        t1 = Sys.time()	      
        #plan(multisession, workers = 5)
        #test.list.tmp <- future.apply::future_lapply(paste0('model_runs_',1:5), year_fun)
        test.list.tmp <- lapply(paste0('model_runs_',1),year_fun)          
        #plan(sequential)
        t1 = Sys.time() - t1
        
        test.list <- c(test.list, test.list.tmp)
        rm(test.list.tmp)
      }
      print(length(test.list))
      #tmp.length <- lapply(test.list, length)
      #print(tmp.length)
      
      # Need to put this into a matrix  
      # Creating matrix that contains row index, mean cell value, and sd cell value
      # Doing this through a function that is in the 0_Functions script
      # First getting names of columns we need to take mean and sd of
      # Have to do independently for crop and pasture
      column.names.crop = names(test.list[[1]]$LandUse)[grep("prop_crop_20", names(test.list[[1]]$LandUse))]
      # column.names.crop = column.names.crop[-grep("pred", column.names.crop)]
      column.names.pasture = names(test.list[[1]]$LandUse)[grep("prop_pasture_20", names(test.list[[1]]$LandUse))]
      # column.names.pasture = column.names.pasture[-grep("pred", column.names.pasture)]
      # Running functions that creates rasters
      raster.list.crop <- mean.sd.fun(column.names = column.names.crop)
      raster.list.pasture <- mean.sd.fun(column.names = column.names.pasture)
      
      # And writing rasters
      for(i in 1:length(raster.list.crop)) {
        # Crops. Un-comment if you want to produce images as well
        # png(filename = paste("Outputs/BAU/pngs/",
        #                      iso3c[c],
        #                      region.names[k],
        #                      names(raster.list.crop)[i], 
        #                      "_BAU_Crop.png",
        #                      sep = ""),
        #     height = 20, width = 40, units = "cm", res = 72)
        # plot(raster.list.crop[[i]], zlim = c(-.1,1.1))
        # dev.off()
        # Writing rasters
        writeRaster(raster.list.crop[[i]],
                    paste0(getwd(),"/Land Forecast Outputs/EAT_Lancet/tifs/",
                           iso3c[c],
                           region.names[k],
                           names(raster.list.crop)[i],
                           "_EAT_Lancet_Crop.tif"),
                    overwrite = TRUE)
        # Pasture
        # png(filename = paste("Outputs/BAU_EATLancet/pngs/",
        #                      iso3c[c],
        #                      region.names[k],
        #                      names(raster.list.pasture)[i], 
        #                      "_BAU_EATLancet_Pasture.png",
        #                      sep = ""),
        #     height = 20, width = 40, units = "cm", res = 72)
        # plot(raster.list.pasture[[i]], zlim = c(-.1,1.1))
        # dev.off()
        # Writing rasters
        writeRaster(raster.list.pasture[[i]],
                    paste0(getwd(),"/Land Forecast Outputs/EAT_Lancet/tifs/",
                           iso3c[c],
                           region.names[k],
                           names(raster.list.pasture)[i],
                           "_EAT_Lancet_Pasture.tif"),
                    overwrite = TRUE)
      }
      # Removing items to clear space
      rm(raster.list.crop, raster.list.pasture, model_data_current, test.list, row.index.raster)
      
      # removing tmp raster files periodically
      rasterTmpFile(prefix = 'r_tmp_')
      # showTmpFiles()
      removeTmpFiles(h = 0)
      
    }
  }
}


Sys.time() - t1 # Don't do this, it's just depressing

# Quitting r to reset memory
# quit(save = "no", status = 0, runLast = TRUE)
# 
# par(mfrow = c(3,3))
# 
# for(i in seq(1,15,by = 2)) {
#   plot(raster.list.crop[[i]],
#        main = names(raster.list.crop)[[i]],
#        zlim = c(0,1))
# }
# 
# dev.off()


