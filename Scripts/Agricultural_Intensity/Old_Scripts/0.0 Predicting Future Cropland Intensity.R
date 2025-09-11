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
# Confusion matrices for these are plotted
# As well as a spiffed up nice-looking version usable in the SI of a paper
#######

Sys.setenv(RETICULATE_PYTHON = "/usr/local/bin/python3")

# Libraries ----
# Way more than needed...
library(raster)
library(plyr)
library(dplyr)
library(reticulate)
library(OpenImageR)
library(parallel)
library(MASS)
library(doMC)
library(rlang)
library(dismo)
library(rJava)
library(randomForest)
library(plyr)
library(dplyr)
library(purrr)
library(rgdal)
library(raster)
library(gdalUtils)
library(doMC)
library(SpaDES)
library(SpaDES.tools)
library(XLConnect)
library(R.utils)
library(parallel)
library(MASS)
library(snow)
library(doParallel)
library(gbm)
library(nnet)
library(caret)
library(ggspatial)

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


# Importing and reprojecting rasters----
### 
# Crop extent and convolution
crop2010 <- raster("/Users/maclark/Desktop/Mike's Files/Global Mollweide 1.5km/Crop2010.Mollweide.1.5km.tif")
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

### 
# Land cover raster
cover.raster <- raster("/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Land_Cover_Classification/GLS_v02.bj.tif")
# Crop intensity raster from Dietrich 2012
intens <- raster("/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Cropland_Intensity/Current Cropland Intensity.tif")
cover.raster <- projectRaster(cover.raster, intens, method = 'ngb') 

###
# Reprojecting rasters
crs(crop2010) = crs(crop_convolved)
crop2010 <- projectRaster(crop2010, intens, method = 'bilinear')
crop_convolved <- projectRaster(crop_convolved, intens, method = 'bilinear')

# Creating data frame ----
cover.df <-
  data.frame(land_cover_numeric = raster::getValues(cover.raster), # Making data frame
             intens = raster::getValues(intens), # And extracting data from rasters
             crop = raster::getValues(crop2010),
             adj_crop = raster::getValues(crop_convolved)) %>% 
  mutate(land_cover_classification = ifelse(land_cover_numeric %in% c(1:3,12,16),'Minimal', # Creating classification based on the SI from Kehoe et al 2018
                                            ifelse(land_cover_numeric %in% c(4:6,10:11,13,15,17),'Moderate',
                                                   ifelse(land_cover_numeric %in% c(7:9,14,18),'Intense','Other')))) %>%
  filter(land_cover_classification != 'Other') %>% # Limiting to cropland
  mutate(adj_crop_by_intens = adj_crop * intens) %>% # Interaction between crop extent and yields from Dietrich
  transform(land_cover_classification = factor(land_cover_classification, levels = c('Minimal','Moderate','Intense'), ordered = TRUE)) %>% # Ordering land cover intensity for modeling - these are the levels in PREDICTs
  filter(!is.na(adj_crop_by_intens)) %>% # Limiting to complete cases
  filter(!is.na(land_cover_classification)) # Getting rid of oceans


# Checking to make sure plots make sense...
# i.e. intens is typically higher in more 'intense' croplands
# And extent is typically higher in more 'intense' croplands
# Looks good
ggplot(dat = cover.df, aes(x = land_cover_classification, y = intens)) +
  geom_boxplot()
ggplot(dat = cover.df, aes(x = land_cover_classification, y = adj_crop)) +
  geom_boxplot()
ggplot(dat = cover.df, aes(x = land_cover_classification, y = crop)) +
  geom_boxplot()
ggplot(dat = cover.df, aes(x = land_cover_classification, y = adj_crop_by_intens)) +
  geom_boxplot()
ggplot(dat = cover.df, aes(x = land_cover_numeric, y = intens * adj_crop * crop, group = land_cover_numeric)) +
  geom_point() +
  geom_boxplot()

# Building regression models ----
# Purpose is to classify cropland intensity as defined in PREDICTs based on yields from dietrich and cropland extent
# Random forest and multinomial models
# Training and testing on fifths of the data

# Setting seed for replication
set.seed(5)
cover.df$Group <-
  kfold(cover.df, 5, by = cover.df$land_cover_numeric)

# Creating lists to save models and model accuracy 
model.acc = list()
mod.list = list()

model.acc.rf = list()
mod.list.rf = list()

# Defining the model for random forest
model = land_cover_classification ~ intens + adj_crop + adj_crop_by_intens

# Running and testing through each fifth of the data
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
# Accuracy looks decent on both (~70%), multinomial model is consistently a litter higher
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

# Creating function to weight models based on probability of each outcome (e.g. minimal, moderate, and high cropland intensity)
weight.fun <-
  function(mod.list, # list of models
           weight.list, # list of model weights
           df, # data frame used to predict model
           mod.type) { # type of model
    for(i in 1:length(mod.list)) { # Looping through model
      probs.tmp = predict(mod.list[[i]], df, 'prob') %>% # predicted probabilities
        as.data.frame(.) # Converting to data frame
      names(probs.tmp) = c('Minimal','Moderate','Intense') # Adding names
      
      probs.keep <- # Data frame with weighted probability for model i
        probs.tmp %>%
        mutate(Minimal = Minimal * weight.list[[i]], # Weighting outcome by the model weight
               Moderate = Moderate * weight.list[[i]],
               Intense = Intense * weight.list[[i]])
      
      if(i == 1) {
        probs.out = probs.keep # And data frame with weighted probability across all models
      } else {
        probs.out = probs.out + probs.keep
      }
    }
    
    probs.out <-
      probs.out %>%
      mutate(Predicted_Classification = # Predicted classification is an ordered factor
               ifelse(Minimal > Moderate & Minimal > Intense, 'Minimal',
             ifelse(Moderate > Intense, 'Moderate', 'Intense')))
    
    # Changing names to be reflective of the model type
    names(probs.out) <- 
      paste(names(probs.out),
            "_",
            mod.type,
            sep = "")
    
    # Adding columns to data frame
    df <- cbind(df, probs.out) 
    
    return(df)
  }

###
# Weighting models and testing accuracy across all 5 models
# Multinomial model
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

# Predicted output is an ordered factor
cover.df <-
  cover.df %>%
  transform(Predicted_Classification_RandomForest = 
              factor(Predicted_Classification_RandomForest, levels = c('Minimal','Moderate','Intense'), ordered = TRUE)) %>%
  transform(Predicted_Classification_Multinomial = 
              factor(Predicted_Classification_Multinomial, levels = c('Minimal','Moderate','Intense'), ordered = TRUE))

# Creating confusion matrixes to look at model fit
conf.matrix.rf = caret::confusionMatrix(cover.df$land_cover_classification, cover.df$Predicted_Classification_RandomForest)
conf.matrix.multi = caret::confusionMatrix(cover.df$land_cover_classification, cover.df$Predicted_Classification_Multinomial)
# Extracting the table
conf.plot.rf = conf.matrix.rf$table
conf.plot.multi = conf.matrix.multi$table

# Converting to relative values
conf.plot.rf[1,] <- conf.plot.rf[1,]/sum(conf.plot.rf[1,])
conf.plot.rf[2,] <- conf.plot.rf[2,]/sum(conf.plot.rf[2,])
conf.plot.rf[3,] <- conf.plot.rf[3,]/sum(conf.plot.rf[3,])

conf.plot.multi[1,] <- conf.plot.multi[1,]/sum(conf.plot.multi[1,])
conf.plot.multi[2,] <- conf.plot.multi[2,]/sum(conf.plot.multi[2,])
conf.plot.multi[3,] <- conf.plot.multi[3,]/sum(conf.plot.multi[3,])

# And plotting to see if either of these are good
# Built in function
# Ideally want values of 1 (or close to it) on a diagonal from lower left to upper right
# And values of 0 (or close to 0) elsewhere
levelplot(conf.plot.rf)
levelplot(conf.plot.multi)

# Oddly enough, the weighted random forest model is predicting much much better across the board
# Even though estimated model accuracy is worse than for the multinomial model
# Multinomial model is having issues predicting intense cropland cover (bins most into either moderate or minimal intensity)

# Looking into random forest matrix in more depth to make sure it's actually working well
# Ideally mcnemars test has a p-value < 0
# This indicates that marginal probabilities of predicting different outcomes are different
# In other words model is working as intended
conf.matrix.rf

# Now making this look good
# Making confusion matrix in a way that meshes with ggplot
multi.sum <-
  cover.df %>%
  group_by(land_cover_classification, Predicted_Classification_Multinomial) %>%
  dplyr::summarise(Count = n()) %>%
  left_join(.,
            cover.df %>% group_by(land_cover_classification) %>%dplyr::summarise(Count_Actual = n())) %>%
  mutate(plot.x = ifelse(land_cover_classification %in% 'Minimal',1,
                         ifelse(land_cover_classification %in% 'Moderate',2,3))) %>%
  mutate(plot.y = ifelse(Predicted_Classification_Multinomial %in% 'Minimal',1,
                         ifelse(Predicted_Classification_Multinomial %in% 'Moderate',2,3))) %>%
  transform(Predicted_Classification_Multinomial = factor(Predicted_Classification_Multinomial, levels = c('Minimal','Moderate','Intense'), ordered = TRUE)) %>%
  mutate(percent_predicted = Count/Count_Actual) %>%
  as.data.frame(.) %>%
  rbind(.,
        data.frame(land_cover_classification = 'Minimal',
                   Predicted_Classification_Multinomial = 'Intense',
                   Count = 0, Count_Actual = sum(cover.df$land_cover_classification %in% 'Minimal'),
                   plot.x = 1, plot.y = 3,
                   percent_predicted = 0)) %>%
  mutate(label = paste('Predicted: ',Count,'\n',
                       'Actual: ', Count_Actual,'\n',
                       'Percent Predicted: ',round(percent_predicted, digits = 2)*100,"%",
                       sep = ""))

# Color scale
my.palette = c('#E6E1B9','#5BC2AE','#0A6164') %>%
  colorRampPalette(.) 
palette.table = 
  data.frame(value = 0:100,
             color = my.palette(101))

# Summarizing values
multi.sum1 <-
  multi.sum %>%
  mutate(value = round(percent_predicted, digits = 2) * 100) %>%
  left_join(.,
            palette.table)
# Manging data to make the plot
multi.sum.plot <- 
  multi.sum1 %>%
  dplyr::select(land_cover_classification, Predicted_Classification_Multinomial,plot.x,plot.y,percent_predicted,color, label) %>%
  mutate(plot.y = as.numeric(land_cover_classification)) %>%
  mutate(plot.x = as.numeric(Predicted_Classification_Multinomial)) %>%
  mutate(color = ifelse(land_cover_classification == Predicted_Classification_Multinomial, 'white','black'))

# And plotting
multi.conf.plot = ggplot(multi.sum.plot) +
  geom_tile(aes(x = plot.x, y = plot.y, fill = percent_predicted * 100)) +
  geom_text(aes(x = plot.x, y = plot.y, label = label, colour = color), size = 2.5, show.legend = FALSE) +
  scale_fill_gradientn(colours=my.palette(100), limits = c(0,100), expand = c(0,0)) +
  scale_colour_manual(values = c('black','white')) +
  scale_x_continuous(breaks = c(1,2,3),labels = (c('Minimal','Moderate','Intense'))) +
  scale_y_reverse(breaks = c(1,2,3),labels = rev(c('Intense','Moderate','Minimal'))) +
  labs(x = 'Predicted Crop Intensity', y = 'Observed Crop Intensity', fill = 'Relative\nPercent') +
  theme(axis.text = element_text(size = 7.5)) +
  theme(axis.title = element_text(size = 10)) +
  ggtitle("Weighted Multinomial Model")

rf.sum <-
  cover.df %>%
  group_by(land_cover_classification, Predicted_Classification_RandomForest) %>%
  dplyr::summarise(Count = n()) %>%
  left_join(.,
            cover.df %>% group_by(land_cover_classification) %>% dplyr::summarise(Count_Actual = n())) %>%
  mutate(plot.x = ifelse(land_cover_classification %in% 'Minimal',1,
                         ifelse(land_cover_classification %in% 'Moderate',2,3))) %>%
  mutate(plot.y = ifelse(Predicted_Classification_RandomForest %in% 'Minimal',1,
                         ifelse(Predicted_Classification_RandomForest %in% 'Moderate',2,3))) %>%
  transform(Predicted_Classification_RandomForest = factor(Predicted_Classification_RandomForest, levels = c('Minimal','Moderate','Intense'), ordered = TRUE)) %>%
  mutate(percent_predicted = Count/Count_Actual) %>%
  as.data.frame(.) %>%
  rbind(.,
        data.frame(land_cover_classification = 'Minimal',
                   Predicted_Classification_RandomForest = 'Intense',
                   Count = 0, Count_Actual = sum(cover.df$land_cover_classification %in% 'Minimal'),
                   plot.x = 1, plot.y = 3,
                   percent_predicted = 0)) %>%
  mutate(label = paste('Predicted: ',Count,'\n',
                       'Actual: ', Count_Actual,'\n',
                       'Percent Predicted: ',round(percent_predicted, digits = 2)*100,"%",
                       sep = ""))

# Color scale
my.palette = c('#E6E1B9','#5BC2AE','#0A6164') %>%
  colorRampPalette(.) 
palette.table = 
  data.frame(value = 0:100,
             color = my.palette(101))

# Summarizing values
rf.sum1 <-
  rf.sum %>%
  mutate(value = round(percent_predicted, digits = 2) * 100) %>%
  left_join(.,
            palette.table)
# Manging data to make the plot
rf.sum.plot <- 
  rf.sum1 %>%
  dplyr::select(land_cover_classification, Predicted_Classification_RandomForest,plot.x,plot.y,percent_predicted,color, label) %>%
  mutate(plot.y = as.numeric(land_cover_classification)) %>%
  mutate(plot.x = as.numeric(Predicted_Classification_RandomForest)) %>%
  mutate(color = ifelse(land_cover_classification == Predicted_Classification_RandomForest, 'white','black'))

# And plotting
rf.conf.plot = ggplot(rf.sum.plot) +
  geom_tile(aes(x = plot.x, y = plot.y, fill = percent_predicted * 100)) +
  geom_text(aes(x = plot.x, y = plot.y, label = label, colour = color), size = 2.5, show.legend = FALSE) +
  scale_fill_gradientn(colours=my.palette(100), limits = c(0,100), expand = c(0,0)) +
  scale_colour_manual(values = c('black','white')) +
  scale_x_continuous(breaks = c(1,2,3),labels = (c('Minimal','Moderate','Intense'))) +
  scale_y_reverse(breaks = c(1,2,3),labels = rev(c('Intense','Moderate','Minimal'))) +
  labs(x = 'Predicted Crop Intensity', y = 'Observed Crop Intensity', fill = 'Relative\nPercent') +
  theme(axis.text = element_text(size = 7.5)) +
  theme(axis.title = element_text(size = 10)) +
  ggtitle('Weighted Random Forest Model')

conf.plots = plot_grid(rf.conf.plot, multi.conf.plot,nrow = 1, align = 'hv')

ggsave("/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Land_Cover_Classification/OutputPlots/ConfusionMatrixPlots_Crop_27January2020.pdf",
       width = 14, height = 6, units = 'in')

# Now using models to forecast 2050 cropland intensity in the different food system scenarios ----
# Getting projected crop extent
projected.crop.rasters <- 
  list.files("/Users/maclark/Desktop/WorldAgProjections_EATLancet",full.names=TRUE) %>%
  paste(.,'/Crop_2045_2050.tif', sep = '') %>%
  .[-grep('Store',.)]

# Getting list of crop scenarios
scenarios <-
  list.files("/Users/maclark/Desktop/WorldAgProjections_EATLancet") %>%
  .[-grep('Store',.)]

# Adding 2010
projected.crop.rasters <-
  c(projected.crop.rasters,"/Users/maclark/Desktop/Mike's Files/Global Mollweide 1.5km/Crop2010.Mollweide.1.5km.tif")
scenarios <- 
  c(scenarios, 'Crop2010')

# List of intensity rasters
intensity.rasters <-
  c("/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Cropland_Intensity/Current Cropland Intensity.tif",
    "/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Cropland_Intensity/Projecte Intensity 2050 BAU.tif",
    "/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Cropland_Intensity/Projecte Intensity 2050 Closed.tif")
intens.year = c(2010,2050,2050)

# Making loop to estimate cropland intensity in the different scenarios
# This takes a while
# Roughly 40 minutes for each cropland scenario
for(i in 1:length(projected.crop.rasters)) {
  # Importing raster
  tmp.raster = raster(projected.crop.rasters[[i]])
  
  # Convolving raster
  # Creating buffer for conolve
  circle = matrix(1, nrow = 3, ncol = 3)
  py_run_string("circle = copy.copy(r.circle)")
  
  # Converting raster to matrix
  input_array <- matrix(raster::getValues(tmp.raster), nrow = tmp.raster@nrows, byrow = TRUE)
  # Crop convolved
  py_run_string("input_array = copy.copy(r.input_array)")
  py_run_string("input_array[np.isnan(input_array)] = 0")
  # py_run_string("input_array[input_array > 0] = 1")
  py_run_string("crop_convolved = convolve(input_array, circle, mode='constant')")
  
  # Converting back into a raster
  crop_convolved <- raster(matrix(py$crop_convolved,byrow = FALSE, nrow = tmp.raster@nrows),
                           crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
  extent(crop_convolved) <- extent(tmp.raster)
  
  # Getting appropriate intensity raster intensity raster
  if(scenarios[[i]] %in% c('Yields','Combined')) {
    intens.raster <- raster(intensity.rasters[[3]])
  } else if(scenarios[[i]] %in% 'Crop2010') {
    intens.raster <- raster(intensity.rasters[[1]])
  } else {
    intens.raster = raster(intensity.rasters[[2]])
  }
  
  # Reprojecting intensity raster
  # If statement to save time/memory
  if(extent(intens.raster) == extent(crop_convolved) & 
     res(intens.raster)[1] == res(crop_convolved)[1] &
     res(intens.raster)[2] == res(crop_convolved)[2]) {
    
  } else {
    intens.raster <-
      projectRaster(intens.raster, crop_convolved, method = 'bilinear')  
  }
  
  # Clearing memory space
  removeTmpFiles()
  gc()
  
  # Aggregating rasters for memory space
  # if your computer is large enough, you can skip this step
  # Aggregating first one to get template
  intens.raster1 <-
    aggregate(intens.raster, fact = 10, fun = 'mean', na.rm = TRUE)
  
  # Reprojecting others so that I use consistent methodology across the rasters
  intens.raster1 <-
    projectRaster(intens.raster, intens.raster1, method = 'bilinear')
  crop_convolved1 <-
    projectRaster(crop_convolved, intens.raster1, method = 'bilinear')
  
  # Need to do the same with countries...
  country.raster <- raster("/Users/maclark/Desktop/Mike's Files/Global Mollweide 1.5km/MollweideCountryID_1.5km.tif")
  
  country.raster1 <-
    projectRaster(country.raster, intens.raster1, method = 'ngb')
  
  # Assigning intensity to cropland cells with NA values for the intens method from Dietrich
  # Assuming these cells are assigned the average intensity within the country
  
  # Creating data frame
  # Contains data for country, intensity, and crop extent
  tmp.df <- 
    data.frame(country = raster::getValues(country.raster1),
               intens = raster::getValues(intens.raster1),
               cropextent = raster::getValues(crop_convolved1)) %>%
    mutate(intens = intens * cropextent) %>%
    filter(!is.na(intens)) %>%
    mutate(count = 1)
  
  # Summary table
  # Used to get average intensity by country
  tmp.sum <-
    rowsum(tmp.df, group = tmp.df$country, na.rm = TRUE) %>%
    mutate(intens = intens / cropextent) %>%
    mutate(country = country/count)
  
  # Replacing values in the raster
  # With average intensity by country in cells where this dataa is not available
  tmp.df <-
    data.frame(country = raster::getValues(country.raster1),
               intens = raster::getValues(intens.raster1),
               crop = raster::getValues(crop_convolved1),
               cell_index = 1:(intens.raster1@nrows * intens.raster1@ncols))
    # filter(!is.na(country)) %>%
    # filter(crop > 0) %>%
    # filter(!is.na(crop))
  
  # Looping through countries to update values
  for(k in 1:nrow(tmp.sum)) {
    tmp.df$intens[tmp.df$country %in% tmp.sum$country[k] &
                    is.na(tmp.df$intens)] <-
      tmp.sum$intens[k]
  }
  
  # Recreating the intens raster
  intens.raster1 <-
    raster(matrix(tmp.df$intens,byrow=TRUE,nrow=intens.raster1@nrows),
           crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
  
  extent(intens.raster1) = extent(crop_convolved1)
  
  # Cells without cropland have no cropland intensity
  
  # Clearing memory space
  removeTmpFiles()
  gc()
  
  # And reprojecting to original extent/etc
  intens.raster <-
    projectRaster(intens.raster1, intens.raster, method = 'bilinear')
  
  # Clearing memory space
  removeTmpFiles()
  gc()
  
  # Creating a data frame with cropland values
  # This data frame is used to predict classification of intensity
  predict.df <-
    data.frame(intens = raster::getValues(intens.raster),
               adj_crop = raster::getValues(crop_convolved),
               crop_values = raster::getValues(tmp.raster),
               cell_index = 1:(intens.raster@ncols * intens.raster@nrows)) %>% # Index used to order cells later
    filter(!is.na(crop_values)) %>% # Filtering to save time and memory requirements
    filter(!is.na(intens)) %>%
    filter(adj_crop > 0) %>%
    mutate(adj_crop_by_intens = intens * adj_crop) %>%
    mutate(split_index = 1:nrow(.)) %>%
    mutate(split_group = ceiling(split_index/1000000)) # Creating group to split into individual data frames
  
  # Getting missing indices
  # This is a list of the cell values that do not contain cropland
  missing.indices <-
    1:(intens.raster@ncols * intens.raster@nrows) %>%
    .[!(. %in% predict.df$cell_index)]
  
  # Recreating raster to make sure this works properly
  # predict.df <-
  #   predict.df[order(predict.df$cell_index),]
  # 
  # check.raster <- raster(matrix(predict.df$adj_crop_by_intens, byrow = TRUE, nrow = crop_convolved@nrows),
  #                        crs = crs(intens.raster))
  # extent(check.raster) = extent(intens.raster)
  
  # Splitting df into list of data frames
  # Where each data from contains 1 million cells
  split.dfs <-
    split(predict.df, predict.df$split_group)
  
  # Clearing memory space
  rm(predict.df)
  gc()
  
  # Parallelizing
  # First making sure this works
  # split.dfs1 = mclapply(1:10,function(x) {
  #   trial = weight.fun(mod.list = mod.list.rf,
  #                      weight.list = model.acc.rf1,
  #                      df = split.dfs[[x]],
  #                      mod.type = 'RandomForest') %>%
  #     dplyr::select(cell_index,Predicted_Classification_RandomForest)
  #   
  #   # return(split.dfs)
  # }, mc.cores = 10)
  # 
  
  # Forecasting cropland intensity for each cell
  split.dfs1 = mclapply(1:length(split.dfs),function(x) {
    trial = weight.fun(mod.list = mod.list.rf,
                       weight.list = model.acc.rf1,
                       df = split.dfs[[x]],
                       mod.type = 'RandomForest') %>%
      dplyr::select(cell_index,Predicted_Classification_RandomForest)

    # return(split.dfs)
  }, mc.cores = 10)
  
  # Combining the different dfs into one big one
  for(j in 1:length(split.dfs1)) {
    if(j %in% 1) {
      out.df <- split.dfs1[[j]]
    } else {
      out.df <- 
        rbind(out.df,
              split.dfs1[[j]])
    }
  }

  # Adding in cells dropped from original data set
  out.df1 <-
    out.df %>%
    dplyr::select(cell_index,Predicted_Classification_RandomForest) %>%
    rbind(.,
          data.frame(cell_index = missing.indices,
                     Predicted_Classification_RandomForest = 'None'))
  
  # Ordering to make raster
  out.df1 <-
    out.df1[order(out.df1$cell_index),]
  
  # Converting to numeric values
  # Classification is minimal/moderate/intense, need values of 1/2/3 to convert into a raster
  out.df1 <-
    out.df1 %>%
    mutate(intens_values = ifelse(Predicted_Classification_RandomForest %in% 'None',0,
                                  ifelse(is.na(Predicted_Classification_RandomForest),0,
                                  ifelse(Predicted_Classification_RandomForest %in% 'Minimal',1,
                                         ifelse(Predicted_Classification_RandomForest %in% 'Moderate',2,3)))))
  
  
  # Checking to make sure classification worked
  # unique(out.df1$Predicted_Classification_RandomForest[out.df1$intens_values %in% 0])
  # unique(out.df1$Predicted_Classification_RandomForest[out.df1$intens_values %in% 1])
  # unique(out.df1$Predicted_Classification_RandomForest[out.df1$intens_values %in% 2])
  # unique(out.df1$Predicted_Classification_RandomForest[out.df1$intens_values %in% 3])
  
  

  # Making and saving raster
  out.raster <-
    raster(matrix(out.df1$intens_values, byrow = TRUE, nrow = intens.raster@nrows),
           crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
  extent(out.raster) = extent(tmp.raster)
  
  # Converting ocean cells to NAs
  out.raster[is.na(tmp.raster)] <- NA
  out.raster[tmp.raster == 0] <- 0
  
  # Plotting
  # plot(out.raster)
  
  # And saving
  writeRaster(out.raster,
              paste("/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Cropland_Intensity/Predicted_Intensity/",
                    'CropIntensity_',
                    scenarios[i],
                    '27January.tif',
                    sep = ''))
  
  # And clearing memory space
  removeTmpFiles()
  gc()
}

# And plotting to make sure these make sense
# Looks good!
file.list = list.files("/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Cropland_Intensity/Predicted_Intensity",
                       full.names = TRUE)

tmp.stack = stack()
for(i in 1:length(file.list)) {
  tmp.raster = raster(file.list[i])
  tmp.stack = stack(tmp.stack,tmp.raster)
}

names(tmp.stack) <- 
  gsub("CropIntensity_","",names(tmp.stack))

plot(tmp.stack)

