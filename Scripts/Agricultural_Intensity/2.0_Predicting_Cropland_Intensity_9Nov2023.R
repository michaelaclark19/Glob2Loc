#!/usr/bin/env Rscript

#####
# Script for predicting future cropland intensity
# This takes models that have been built
# And applies them to forecasts of agricultural land use
##### 

# Importing the right version of python
# Sys.setenv(RETICULATE_PYTHON = "/apps/system/easybuild/software/Anaconda3/2022.05/bin/python")

# libraries
library(raster)
library(plyr)
library(dplyr)
library(reticulate)
library(countrycode)
library(randomForest)

# Python
library(reticulate)

# Loading python packages
py_run_string("import numpy as np")
py_run_string("import pandas as pd")
py_run_string("from skimage import io, morphology, measure")
py_run_string("from scipy.ndimage import convolve")
py_run_string("from PIL import Image")
py_run_string("import copy")

# Setting working directory
setwd("/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity")

# Scenario info
luc_scenarios = c('BAU','Diets','Combined')
yields_scenario = c('BAU', 'BAU','Closed')

# Importing functions
source(paste0(getwd(),"/Scripts/Agricultural Intensity/0.0_Functions_for_Predicting_Future_Cropland_Intensity.R"))

# Creating directory for outputs
dir.create(paste0(getwd(),"/Ag Intensity Outputs/Cropland Intensity Forecasts"))

for(i in luc_scenarios) {
  dir.create(paste0(getwd(),"/Ag Intensity Outputs/Cropland Intensity Forecasts/",i))
}


# Python management ----
 # Used to colve buffer around crop extent
# Importing packages
py_run_string("import numpy as np") # Load numpy  package
py_run_string("import pandas as pd") # Load pandas package
py_run_string("from skimage import io, morphology, measure")
py_run_string("from scipy.ndimage import convolve")
py_run_string("from PIL import Image")
py_run_string("import copy")

# Creating buffer for conolve
circle = matrix(1, nrow = 3, ncol = 3)
py_run_string("circle = copy.copy(r.circle)")


# Data management ----
# Importing cropland forecasts
crop.forecasts <- 
  list.files(path = paste0(getwd(),'/Land Forecast Outputs/BAU/Global tifs'),
	     full.names = TRUE)

# Loading models
load(paste0(getwd(),'/Ag Intensity Outputs/Classification_Model_Outputs/Random Forest Model List 28March2023.rda'))
load(paste0(getwd(),'/Ag Intensity Outputs/Classification_Model_Outputs/Random Forest Model Accuracy 28March2023.rda'))

# Getting current crop intensity
intens <- raster(paste0(getwd(),"/Global Mollweide Maps/Current Cropland Intensity.tif"))

# Country raster
country.iso <- 
  raster(paste0(getwd(),"/Global Mollweide Maps/MollweideCountryID_1.5km.tif"))
# ANd reprojecting to intens
country.iso <- projectRaster(country.iso, intens, method = 'ngb')

# Updating intens raster to account for regions without current crop production
# Doing this by assuming that with no current intensity have an intensity equivalent ot he median value of cells within the same country
# Need to do this to make sure that
# if crop expands, there is a baseline intensity that wil increase/decrease
update.intens.function <-
  function(intens_raster, country_raster) {
    # List of iso3s
    tmp.isos <- sort(unique(getValues(country_raster)))
    # Looping through these
    for(iso in tmp.isos) {
      update.value <-
        median(intens_raster[country_raster %in% iso], na.rm = TRUE)
      intens_raster[is.na(intens_raster) & country_raster %in% iso] <- update.value
    }
    return(intens_raster)
  }

intens <- update.intens.function(intens, country.iso)

# Looping through each year
for(i in luc_scenarios) {
  # List of forecasts
  # Limiting to crop expansion
  crop.forecasts <-
    list.files(path = paste0(getwd(),"/Land Forecast Outputs/",i,"/Global tifs"), full.names = TRUE, pattern = 'mean') %>%
    .[grepl('crop',.)]
  
  # Years to loop over
  years <- 
    gsub(".*crop_mean","",crop.forecasts) %>%
    gsub(".tif","",.) %>% unique(.)

# adding in 2010
  years <- c('2010',years)
  # Importing crop yield forecasts
  # Which scneario is this
  yields <-
    read.csv(list.files(path = paste0(getwd(),"/Other Data Inputs/Yield Forecasts/"),
                        pattern = yields_scenario[which(luc_scenarios %in% i)], full.names = TRUE))
  
  # And projecting cropland intnesity in lapply
  # This an be converted to mclapply to parallelize
  # But this takes lots of memory space
  # And risks returning "NULL" values in the list
  # Running this list over lapply isn't ideal
  # But it avoids the potential issue of returning null values after running the script
  # out.projections <-
  #  lapply(years, project.intensity.function)
  
  for(y in years) {
	  # exception for 2010 raster
	  if(y %in% 2010) {
		  crop.forecasts.years <- 
			  paste0(getwd(),'/ESA LandCov Maps/Crop2010_Corrected_GlobCov.tif')
	  } else {
		   # Getting raster for target year
    crop.forecasts.years <-
      crop.forecasts[grepl(y, crop.forecasts)]
	  }

    # Importing raster
    tmp.raster.keep <- raster(crop.forecasts.years)

    # Cropping and reprojecting
    tmp.raster <- projectRaster(tmp.raster.keep, intens, method = 'bilinear')

    # Getting adjacent cropland extent
    # Converting to a matrix for convolution
    input_array <- matrix(getValues(tmp.raster), nrow = tmp.raster@nrows, byrow = TRUE)
    # Crop convolved
    py_run_string("input_array = copy.copy(r.input_array)")
    py_run_string("input_array[np.isnan(input_array)] = 0")
    # py_run_string("input_array[input_array > 0] = 1")
    py_run_string("crop_convolved = convolve(input_array, circle, mode='constant')")
    # Recreating raster
    crop_convolved <- raster(matrix(py$crop_convolved,byrow = FALSE, nrow = tmp.raster@nrows),
                             crs = crs(tmp.raster))
    extent(crop_convolved) <- extent(tmp.raster)

     # Creating data frame and predicting
    cover.df <-
      data.frame(intens = raster::getValues(intens),
                 crop = raster::getValues(tmp.raster),
                 adj_crop = raster::getValues(crop_convolved),
                 iso = as.numeric(getValues(country.iso))) %>%
      mutate(adj_crop_by_intens = adj_crop * intens) %>%
      left_join(., yields %>% mutate(iso = countrycode(ISO3, origin = 'iso3c', destination = 'iso3n')) %>% mutate(iso = ifelse(ISO3 %in% 'ZZZ',728,iso)) %>% unique(.)) %>%
      mutate(yields_increase = .[,which(names(.) %in% paste0('X',gsub(".*_","",y)))]) %>%
      mutate(yields_increase = ifelse(iso %in% 729, unique(.$yields_increase[.$iso %in% 728]),yields_increase)) %>% # updating yields increase for sudan and south sudan
      mutate(yields_increase = ifelse(is.na(yields_increase),1,yields_increase)) %>% # Converting countires with no projections to a yields increase of 0%. This is conservative, but not a terrible issue beause almost all of these countries are island nations
      mutate(intens = intens * yields_increase)

    # Predicting intensity Random forest
    cover.df <-
      weight.fun(mod.list = mod.list.rf,
                 weight.list = model.acc.rf1,
                 df = cover.df,
                 mod.type = 'RandomForest')

    # Classifying into low medium and high
    # ANd making sure intensity predictions only exist in cells where there is cropland
    cover.df <-
      cover.df %>%
      mutate(predicted_intens = ifelse(Predicted_Classification_RandomForest %in% 'Minimal',1,
                                       ifelse(Predicted_Classification_RandomForest %in% 'Moderate',2,
                                              ifelse(Predicted_Classification_RandomForest %in% 'Intense',3,0)))) %>%
      mutate(predicted_intens = ifelse(crop %in% 0 | is.na(crop), 0, predicted_intens))

table(cover.df$Predicted_Classification_RandomForest)

 # Making intensity raster
    intens.cov.raster <-
      raster(matrix(cover.df$predicted_intens, byrow = TRUE, nrow = intens@nrows), crs = crs(intens))
    extent(intens.cov.raster) <- extent(intens)
    intens.cov.raster[intens.cov.raster %in% 0] <- NA

    # and reprojecting
    intens.cov.raster <- projectRaster(intens.cov.raster, tmp.raster.keep, method = 'ngb', na.rm = TRUE)

    # And saving
    writeRaster(intens.cov.raster,
                paste0(getwd(),"/Ag Intensity Outputs/Cropland Intensity Forecasts/",i,"/Cropland_Intensity_",y,".tif"),
                overwrite = TRUE)


  }



    # # And saving
    # writeRaster(intens.cov.raster,
    #             paste0(getwd(),))
  }

###
# And projecting current intensity
# And writing it into each of the scenario directories
crop.forecasts <- paste0(getwd(),"/ESA LandCov Maps/Crop2010_Corrected_GlobCov.tif")
years <- c(2010)
# Setting directory to run on BAU file
i = 'BAU'
raster_2010 <- project.intensity.function(years)

# and saving raster
writeRaster(raster_2010,
                paste0(getwd(),"/Ag Intensity Outputs/Cropland Intensity Forecasts/",i,"/Cropland_Intensity_",y,".tif"),
                overwrite = TRUE)
