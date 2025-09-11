#######
# Developing framework to predict cropland intensity in different ag scenarios
# raw data input is land cover map from kehoe et al 2018, crop yields from dietrich et al 2012, and cropland extent in 2010 as reported by MODIS
# Land cover from Kehoe: https://box.hu-berlin.de/d/053f45f377/
# Dietrich data from: https://www4.unfccc.int/sites/NAPC/Country%20Documents/General/1-s2.0-S0304380012001093-main.pdf

# Approach uses regression classification models to estimate land cover in Kehoe maps
# Based on crop yields from dietrich and crop extent from MODIS
# Tested multinomial model and random forest model
# Further subdivided data into fifths, then trained and tested data on each fifth
# Accuracy of both models is ~70% (compared to 33% by random chance)
# When weighting the models, it is clear that random forest is much much better
# Confusion matrices for these are plotted for use in SI and methods
#######

# Importing the right version of python
Sys.setenv(RETICULATE_PYTHON = "/apps/system/easybuild/software/Anaconda3/2022.05/bin/python")

# Libraries ----
# Way more than needed...
library(reticulate)
library(raster)
library(plyr)
library(dplyr)
library(dismo)
library(ggplot2)
library(randomForest)
library(nnet)

# Setting working directory
setwd("/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity")

# Importing functions
source(paste0(getwd(),"/Scripts/Agricultural_Intensity/0.0_Functions_for_Predicting_Future_Cropland_Intensity.R"))

# Python management ----
# Importing packages
py_run_string("import numpy as np") # Load numpy  package
py_run_string("import pandas as pd") # Load pandas package
py_run_string("from skimage import io, morphology, measure")
py_run_string("from scipy.ndimage import convolve")
py_run_string("from PIL import Image")
py_run_string("import copy")

# Creating buffer for conolve ----
circle = matrix(1, nrow = 3, ncol = 3)
py_run_string("circle = copy.copy(r.circle)")


# Data management - Importing and reprojecting rasters----
# Doing this to get everything to the same crs/resolution/extent/etc
# Before building the classification models

### Crop extent and convolution
# Importing crop extent in 2010
crop2010 <- raster(paste0(getwd(),"/ESA LandCov Maps/Crop2010_Corrected_GlobCov.tif"))
# Converting to a matrix for convolution
input_array <- matrix(getValues(crop2010), nrow = crop2010@nrows, byrow = TRUE)
# Crop convolved
py_run_string("input_array = copy.copy(r.input_array)")
py_run_string("input_array[np.isnan(input_array)] = 0")
# py_run_string("input_array[input_array > 0] = 1")
py_run_string("crop_convolved = convolve(input_array, circle, mode='constant')")
# Recreating raster
crop_convolved <- raster(matrix(py$crop_convolved,byrow = FALSE, nrow = crop2010@nrows),
                         crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
extent(crop_convolved) <- extent(crop2010)

### Land cover raster
# Importing land cover intensity raster from Kehoe et al 2017
cover.raster <- raster(paste0(getwd(),"/Global Mollweide Maps/Land_Cover_Classification_GLS_v02.bj.tif"))
# Crop intensity raster from Dietrich 2012
intens <- raster(paste0(getwd(),"/Global Mollweide Maps/Current Cropland Intensity.tif"))
# Reprojecting
cover.raster <- projectRaster(cover.raster, intens, method = 'ngb') 

### Reprojecting crop rasters
crs(crop2010) = crs(crop_convolved) # Setting crs so that crop 2010 can be reprojected
crop2010 <- projectRaster(crop2010, intens, method = 'bilinear') # Reprojecting crop 2010
crop_convolved <- projectRaster(crop_convolved, intens, method = 'bilinear') # reprojecting crop convolved crop 2010

# Data management - Creating data frame ----
# This is what is used as inputs into the classification models

# setting seed for replicability
set.seed(19)

# And creating df
cover.df <-
  data.frame(land_cover_numeric = raster::getValues(cover.raster), # Extracting values from the rasters
             intens = raster::getValues(intens), 
             crop = raster::getValues(crop2010),
             adj_crop = raster::getValues(crop_convolved)) %>% 
  mutate(land_cover_classification = ifelse(land_cover_numeric %in% c(1:3,12,16),'Minimal', # Creating classification based on the SI from Kehoe et al 2018
                                            ifelse(land_cover_numeric %in% c(4:6,10:11,13,15,17),'Moderate', # This is cropland intensity
                                                   ifelse(land_cover_numeric %in% c(7:9,14,18),'Intense','Other')))) %>% # Or any other land cover type
  filter(land_cover_classification != 'Other') %>% # Limiting to cropland, don't care about other land covers
  mutate(adj_crop_by_intens = adj_crop * intens) %>% # Interaction between crop extent and yields from Dietrich
  transform(land_cover_classification = factor(land_cover_classification, levels = c('Minimal','Moderate','Intense'), ordered = TRUE)) %>% # Ordering land cover intensity for modeling - these are the levels in PREDICTs
  filter(!is.na(adj_crop_by_intens)) %>% # Limiting to complete cases
  filter(!is.na(land_cover_classification)) %>% # Getting rid of oceans
  mutate(Group = kfold(., 5, by = land_cover_classification)) # And splitting df into five groups for the regression modelling

# Making a few quick plots for sense checks ----
# These make sense
# intensity (from dietrich) = has some correlation with intensity from kehoe
ggplot(dat = cover.df, aes(x = land_cover_classification, y = intens)) + geom_boxplot()
# Extent of surrounding classification does not have a strong correlation
ggplot(dat = cover.df, aes(x = land_cover_classification, y = adj_crop)) + geom_boxplot()
# Extent of in-cell crop does not have strong correlation
ggplot(dat = cover.df, aes(x = land_cover_classification, y = crop)) + geom_boxplot()
# Intensity * adjacent has strong correlation
ggplot(dat = cover.df, aes(x = land_cover_classification, y = adj_crop_by_intens)) + geom_boxplot()

# Building classification models ----
# Doing a random forest
# Anda  multinomial classification model

# Creating empty lists to save models and model accuracy 
model.acc = list()
mod.list = list()

model.acc.rf = list()
mod.list.rf = list()

# Defining the model for random forest
model = land_cover_classification ~ intens + adj_crop + adj_crop_by_intens

# Running and testing through each fifth of the data
# This first round does not have equal sampling from the land cover types
for(i in 1:5) {
  # Creating dfs for train and test groups
  test.df <- cover.df[cover.df$Group %in% i,]
  train.df <- cover.df[cover.df$Group != i,]
  
  ###
  # Multinomial model
  # Building regression model
  test.mod.multi = 
    multinom(land_cover_classification ~
               intens +
               adj_crop +
               adj_crop_by_intens,
             dat = train.df,
             model = FALSE,
             maxit = 250)
  
  # Appending to list
  mod.list[[i]] = test.mod.multi
  
  # Checking accuracy against test group
  # Where accuracy is defined as observed = predicted
  test.df$predicted = predict(test.mod.multi, test.df)
  test.df1 <- test.df %>% filter(!is.na(predicted))
  
  # Appending accuracy to list to weight models later
  model.acc[[i]] = sum(as.character(test.df1$predicted) == as.character(test.df1$land_cover_classification)) / nrow(test.df1)
  
  ###
  # Repeating for random forest
  # Building regression model
  test.rf <- randomForest(model,
                          data = train.df[,c('land_cover_classification','adj_crop','intens','adj_crop_by_intens')],
                          ntree = 500)
  
  # Predicted to check accuracy
  rf.pred <- predict(test.rf, test.df)
  test.df$predicted_rf = predict(test.rf, test.df)
  test.df1 <- test.df %>% filter(!is.na(predicted))
  
  # Appending accuracy to list to get model weight later
  model.acc.rf[[i]] = sum(as.character(test.df1$predicted_rf) == as.character(test.df1$land_cover_classification)) / nrow(test.df1)
  mod.list.rf[[i]] = test.rf
}

# Looking at accuracy
# Accuracy looks decent on both (~70%), multinomial model is nominally higher (by ~2%)
# Digging a bit deeper
model.acc
model.acc.rf

# Weighting models, where weight is accuracy[model i] ^ 2 / sum (accuracy ^ 2)
model.acc1 <- c()
model.acc.rf1 <- c()
for(i in 1:length(model.acc)) {
  model.acc1 <- c(model.acc1,model.acc[[i]][1])
  model.acc.rf1 <- c(model.acc.rf1, model.acc.rf[[i]][1])
}

model.acc1 <- model.acc1^2 / sum(model.acc1^2)
model.acc.rf1 <- model.acc.rf1^2 / sum(model.acc.rf1^2)

# Predicting outcomes based on relative model accuracy... ----
# Multinomial model
# This function adds columns to the data frame
cover.df <-
  weight.fun(mod.list = mod.list,
             weight.list = model.acc1,
             df = cover.df,
             mod.type = 'Multinomial')

# Random forest
cover.df <-
  weight.fun(mod.list = mod.list.rf,
             weight.list = model.acc.rf1,
             df = cover.df,
             mod.type = 'RandomForest')

# Building level plots/confusion plots to visualise accuracy of the model ----
rf.plot <- conf.plot.function(df = cover.df, model_type = 'RandomForest')
multi.plot <- conf.plot.function(df = cover.df, model_type = 'Multinomial')

# And saving the model outputs ----
dir.create(paste0(getwd(),"/Ag Intensity Outputs/Classification_Model_Outputs"))

save(mod.list.rf, file = paste0(getwd(),'/Ag Intensity Outputs/Classification_Model_Outputs/Random Forest Model List 28March2023.rda'))
save(model.acc.rf1, file = paste0(getwd(),'/Ag Intensity Outputs/Classification_Model_Outputs/Random Forest Model Accuracy 28March2023.rda'))

dir.create(paste0(getwd(),'/Ag Intensity Outputs/Level Plots'))

rf.plot
ggsave(paste0(getwd(),'/Ag Intensity Outputs/Level Plots/Random Forest Ag Intensity Accuracy 28March2023.pdf'),width = 8, height = 6)

multi.plot
ggsave(paste0(getwd(),'/Ag Intensity Outputs/Level Plots/Multinomial Ag Intensity Accuracy 28March2023.pdf'),width = 8, height = 6)

