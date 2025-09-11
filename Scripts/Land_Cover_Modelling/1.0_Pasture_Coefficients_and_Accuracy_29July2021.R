#!/usr/bin/env Rscript

# Read me ----
# Modified version of 
# "GBM Model IUCN Region prob expansion outcome Two-step 24.August.R" plus
# "GBM 22.Sept Crop_MC.R" and "GBM 22.Sept Pasture_MC.R"

# I have changed some of the variable names to make them slightly more 
#     intelligible when it comes to the land-clearance loop
# Some of this comes from "ProbabilitiesOfLandClearance_Current.R"
#
# This version (16 October) fits the pasture GBM after the cropland one
#     AND fits a trimmed-down version of the GBM:
# - Does not include ISO Numeric but does include overall land demand
# - Fits a global model (not separate models to countries that have had cropland expansion, 
#   compared to contraction)

# Packages -----
library(plyr)
library(tidyverse)
library(reshape2)
library(caret)
library(pROC)
library(gbm)
library(e1071)
library(glmnet)
library(raster)
library(rgdal)
library(gridExtra)
library(doMC)
library(dplyr)
library(plyr)
library(nnet)
library(DHARMa)
library(SuppDists)
library(readr)
library(ggplot2)
library(dismo)
registerDoMC(cores = 4)

# Clean up -----
rm(list = ls())
# change cutoff for where we denote if a cell has experienced a change in crop extent
cutoff <- .025

# setting wd
setwd("/Users/maclark/Desktop/Multiple Stresses of Biodiversity")

# Loading multinomial estimate function
source(paste0(getwd(),'/Scripts/Land Cover Modelling/1.0 Multinomial Estimate Function.R'))
source("Scripts/Land Cover Modelling/LandUseChangeFunctions.R")

# getting list of GBM file names
file.list <-
  list.files(path = paste0(getwd(),'/Forecasting_Data'),
             full.names = TRUE, pattern = '29July')


### Prep: creating lists used to assign model outputs
# creating list used to assign values to model output if cells experienced a change in cropland extent
region_list <-
  gsub(".*29July2021","",file.list) %>%
  gsub(".csv",'',.)


region_list_model_change <- paste(gsub(" ","",region_list),"model_change",sep = "_") 
# list used to assign values to model output for model that estimates the change in cropland extent by cell
region_list_model_extentchange_increase <- paste0(region_list,"_model_extent_change_increase")
region_list_model_extentchange_decrease <- paste0(region_list,"_model_extent_change_decrease") 

tmp <- read_csv(file.list[1])

# Data prep: Set predictor names ----
# I am selecting the columns by name not index to account for possible changes
# This is currently only used to test for scaled predictors
predictor_vars_crop <- c('dist_50k_log' ,
                         'dist_50k_log_sq' ,
                         'ag_suitability' ,
                         'PA_binary' ,
                         'prop_to_predict_from_crop' ,
                         'prop_crop_sq' ,
                         'prop_to_predict_from_pasture' ,
                         'prop_pasture_sq' ,
                         'prop_to_predict_from_cell_delta_pasture' ,
                         'prop_cell_delta_pasture_sq' ,
                         'prop_adj_to_predict_from_pasture' ,
                         'prop_adj_pasture_sq' ,
                         'ISO_numeric')
# And plot variables to test residuals
predictor_vars_crop <- c('dist_50k_log' ,
                         'dist_50k_log_sq' ,
                         'ag_suitability' ,
                         'PA_binary' ,
                         'prop_to_predict_from_crop' ,
                         'prop_crop_sq' ,
                         'prop_to_predict_from_pasture' ,
                         'prop_pasture_sq' ,
                         'prop_to_predict_from_cell_delta_pasture' ,
                         'prop_cell_delta_pasture_sq' ,
                         'prop_adj_to_predict_from_pasture' ,
                         'prop_adj_pasture_sq' ,
                         'ISO_numeric')

# creating lists to savemodels
change_list_crop <- list() # first model step
proportion_list_crop_increase <- list() # second model step
proportion_list_crop_decrease <- list() # second model step

# Creating list of names used to save files
file.names.save <- paste0(region_list, "_29July2021_Cutoff_025_Pasture.rda")

# List of files
mods.prob <- list()
mods.increase <- list()
mods.decrease <- list()

# List of coefficients
coef.prob <- list()
coef.increase <- list()
coef.decrease <- list()

# Timing process
t1 = Sys.time()

### big loop to cycle through each file in file.list
for (k in c(1:12)) {
  # Setting seed for replicability
  set.seed(19)
  # Update on progress
  cat(paste0(region_list[k],"\n"))
  
  # importing files and changing names to fit with coefficient names
  # read_csv imports some of the numeric data as a character
  # so using read.csv, even if slower
  dat_orig <- 
    read.csv(file.list[k]) %>%
    dplyr::rename(dist_50k = distance_cities,
                  ag_suitability = suitability,
                  PA_binary = wdpa,
                  prop_to_predict_from_cell_delta_pasture = delta_pasture_2000_2005) %>%
    mutate(country_id = country_id_impact)
  
  # Updating NA values in suitability
  # For some reason, when this data is written, some of the very small numbers (i.e. x e-5) are written as NAs rather than numbers
  # So going back and reconverting these into numbers
  # read.csv imports these cells as NAs
  # read_csv imports these cells as character values, e.g. "5e-10"
  # So using values in read_csv to update values from read.csv
  dat_update <- read_csv(file.list[k]) %>% mutate(country_id = country_id_impact)
  update.index <- (grep("e",dat_update$suitability))
  update.values <-
    as.numeric(gsub('e.*','',dat_update$suitability[update.index])) *
    10^(as.numeric(gsub('.*e','',dat_update$suitability[update.index])))
  
  # And updating
  dat_orig$ag_suitability[update.index] <-
    update.values
  
  # All non PA areas have value of 0
  dat_orig$PA_binary[which(is.na(dat_orig$PA_binary))] <- 0
  
  # Creating increase and decrease columns
  dat_orig[,'has_changed_pasture'] <- NA
  dat_orig[which(dat_orig[,'delta_pasture_2005_2010'] >= cutoff  & !is.na(dat_orig[,'delta_pasture_2005_2010'])),'has_changed_pasture'] <- 'Increase'
  dat_orig[which(dat_orig[,'delta_pasture_2005_2010'] <= -cutoff  & !is.na(dat_orig[,'delta_pasture_2005_2010'])),'has_changed_pasture'] <- 'Decrease'
  dat_orig[which(is.na(dat_orig[,'has_changed_pasture']) & !is.na(dat_orig[,'prop_to_predict_from_cell_delta_pasture'])),'has_changed_pasture'] <- 'aNoChange'
  
  # Data prep: Remove unnecessary columns, rename predictor_vars ----
  # I am renaming a) for consistency b) so that we can more easily understand
  dat_orig <- dplyr::select(dat_orig,
                            row_index,
                            ISO_numeric = country_id,
                            IUCN_region = iucn_region,
                            dist_50k = dist_50k,
                            ag_suitability = ag_suitability,
                            PA_binary = PA_binary,
                            # change_demand_crop = CountryIncrease.Crop.2002.2007,
                            prop_adj_to_predict_from_pasture = adj_Pasture2005,
                            prop_to_predict_from_crop = Crop2005,
                            prop_to_predict_from_cell_delta_pasture = prop_to_predict_from_cell_delta_pasture,
                            pasture_extent_change_2010.2005 = delta_pasture_2005_2010,
                            prop_to_predict_from_pasture = Pasture2005,
                            prop_crop_cell_2000 = Crop2000,
                            prop_crop_cell_2010 = Crop2010,
                            prop_change_cell_pasture_2005.2010 = delta_pasture_2005_2010,
                            has_changed_pasture = has_changed_pasture)
  
  # Setting factor variables
  dat_orig <-
    dat_orig %>%
    mutate(ISO_numeric = factor(ISO_numeric, levels = sort(unique(dat_orig$ISO_numeric)), ordered = FALSE))
  
  # Data prep: Getting model data
  # Limiting to specific variables and complete cases
  dat_mod <- dat_orig[,c('dist_50k','ag_suitability','PA_binary','prop_adj_to_predict_from_pasture','prop_to_predict_from_crop','prop_to_predict_from_pasture',
                         'prop_to_predict_from_cell_delta_pasture',"has_changed_pasture", "prop_change_cell_pasture_2005.2010","ISO_numeric")]
  dat_mod <- dat_mod[complete.cases(dat_mod),]
  
  # Removing original dat
  rm(dat_orig)
  
  # Adding columns for square transformations
  dat_mod$prop_crop_sq <- dat_mod$prop_to_predict_from_crop^2
  dat_mod$prop_pasture_sq <- dat_mod$prop_to_predict_from_pasture^2
  dat_mod$prop_cell_delta_pasture_sq <- dat_mod$prop_to_predict_from_cell_delta_pasture^2
  dat_mod$prop_adj_pasture_sq <- dat_mod$prop_adj_to_predict_from_pasture^2
  # Adding columns for log transformations
  dat_mod$dist_50k_log <- log(dat_mod$dist_50k + .001)
  dat_mod$dist_50k_log_sq <- dat_mod$dist_50k_log^2
  # dat_mod$ISO_numeric <- as.character(dat_mod$ISO_numeric)
  
  # Creating index for testing
  dat_mod$tmp_index <- 1:nrow(dat_mod)
  
  # Getting multinomial model
  outcomeName <- "has_changed_pasture" # outcome factor variable
  # Attaching data set
  
  # Getting count by iso3n
  dat.ison <- 
    table(dat_mod$ISO_numeric) %>%
    .[. >= 5] %>% names(.)
  
  # Filtering countries with fewer than 5 data points
  dat_mod <-
    dat_mod %>%
    filter(ISO_numeric %in% dat.ison)
  
  # Limiting to equal number of each of no change, increase, or decrease in ag land cover
  sample.size = min(table(dat_mod$has_changed_pasture))
  
  dat_mod_new <-
    rbind(sample_n(dat_mod %>% filter(has_changed_pasture %in% 'aNoChange'),sample.size),
          sample_n(dat_mod %>% filter(has_changed_pasture %in% 'Decrease'),sample.size),
          sample_n(dat_mod %>% filter(has_changed_pasture %in% 'Increase'),sample.size))
  
  
  # Creating kfold group
  dat_mod_new$group <- kfold(dat_mod_new, k = 5, by = dat_mod_new$ISO_numeric)
  
  coef.multinom.list <- list()
  coef.se.list <- list()
  coef.acc.list <- list()
  
  # And looping through five iterations to weight based on relative accuracy
  for(l in 1:5) {
    # Training data
    dat_train <- dat_mod_new %>% filter(!(group %in% l))
    # Testing data
    # (all other regional data that isn't in the subset)
    dat_test <- rbind(dat_mod_new %>% filter(group %in% l) %>% dplyr::select(-group),
                      dat_mod %>% filter(!(tmp_index %in% dat_mod_new$tmp_index)))
    
    # Model
    # Getting multinomial regression model
    # Getting multinomial regression model
    objModel <- multinom(has_changed_pasture ~ dist_50k_log +
                           dist_50k_log_sq +
                           # ag_suitability +
                           PA_binary +
                           prop_to_predict_from_crop +
                           prop_crop_sq +
                           prop_to_predict_from_pasture +
                           prop_pasture_sq +
                           prop_adj_to_predict_from_pasture +
                           prop_adj_pasture_sq +
                           prop_to_predict_from_cell_delta_pasture +
                           prop_cell_delta_pasture_sq +
                           factor(ISO_numeric, ordered = FALSE),
                         dat = dat_train,
                         model = FALSE,
                         maxit = 250)
    
    # Saving to list
    coef.multinom.mean <- as.data.frame(summary(objModel)$coefficients)
    coef.multinom.se <- as.data.frame(summary(objModel)$standard.errors)
    
    coef.multinom.list[[l]] <- coef.multinom.mean
    coef.se.list[[l]] <- coef.multinom.se
    
    
    
    # Getting accuracy and plotting ----
    iso3s <-
      sort(unique(dat_mod$ISO_numeric))
    
    # Gtting variable list
    vars = names(coef.multinom.mean)[2:(min(grep('factor',names(coef.multinom.mean))) - 1)]
    
    # Looping through iso3s
    # Cannot use R's predict function because not every country is in the data frame used to predict
    # Output data frame
    dat.out <- data.frame()
    for(i in iso3s) {
      pred_tmp <- multinom.predict.2b.trial(dat = dat_test %>% filter(ISO_numeric %in% i),
                                            coefs = coef.multinom.mean,
                                            vars = vars,
                                            factor_vars = 'ISO_numeric')
      dat.out <- 
        rbind(dat.out,
              dat_test %>% filter(ISO_numeric %in% i) %>% cbind(pred_tmp))
    }
    
    # And checking this works really quick
    # Not that this only works when k == 1
    # preds.check <-
    #   predict(objModel,
    #           dat_test %>% filter(ISO_numeric %in% 192),
    #           'prob')
    # 
    # preds.check <-
    #   cbind(dat.out %>% filter(ISO_numeric %in% 332),
    #         preds.check)
    
    
    dat.out <-
      dat.out %>%
      transform(has_changed_pasture = factor(has_changed_pasture, levels = c('aNoChange','Decrease','Increase'))) %>%
      mutate(predicted_values = ifelse(aNoChange > Decrease & aNoChange > Increase,'aNoChange',
                                       ifelse(Decrease > Increase, 'Decrease','Increase'))) %>%
      transform(predicted_values = factor(predicted_values, levels = c('aNoChange','Decrease','Increase')))
    
    # Accuracy
    accuracy <- caret::confusionMatrix(data = dat.out$predicted_values,
                                       reference = dat.out$has_changed_pasture)
    
    coef.acc.list[[l]] <- accuracy
  }
  
  # Weighting models based on RMSE
  out.weights.increase <- 
    data.frame(variable = unique(c(colnames(coef.multinom.list[[1]]),
                                   colnames(coef.multinom.list[[2]]),
                                   colnames(coef.multinom.list[[3]]),
                                   colnames(coef.multinom.list[[4]]),
                                   colnames(coef.multinom.list[[5]]))),
               mod_1 = NA, mod_2 = NA, mod_3 = NA, mod_4 = NA, mod_5 = NA)
  
  # Looping through to calculate weighted coefficients
  for(l in 1:5) {
    # Making data frame saying which variables are in which model
    out.weights.increase[out.weights.increase$variable %in% colnames(coef.multinom.list[[l]]),paste0('mod_',l)] <- coef.acc.list[[l]]$overall[1]
    
  }
  # Weighting
  out.weights.increase[is.na(out.weights.increase)] <- 0
  out.weights.increase[,which(names(out.weights.increase) %in% 'mod_1') : which(names(out.weights.increase) %in% 'mod_5')] <-
    out.weights.increase[,which(names(out.weights.increase) %in% 'mod_1') : which(names(out.weights.increase) %in% 'mod_5')]^2 /
    rowSums(out.weights.increase[,which(names(out.weights.increase) %in% 'mod_1') : which(names(out.weights.increase) %in% 'mod_5')]^2)
  
  out.weights.increase$Estimate <- 0
  out.weights.increase$`Std. Error` <- 0
  
  out.weights.decrease <- out.weights.increase
  
  # And looping again
  for(l in 1:5) {
    for(j in colnames(coef.multinom.list[[l]])) {
      # Increasing coefficients
      # Estimate
      out.weights.increase$Estimate[out.weights.increase$variable %in% j] <-
        out.weights.increase$Estimate[out.weights.increase$variable %in% j] +
        coef.multinom.list[[l]][which(rownames(coef.multinom.list[[l]]) %in% 'Increase'),colnames(coef.multinom.list[[l]]) %in% j] * 
        out.weights.increase[out.weights.increase$variable %in% j,paste0('mod_',l)]
      # Standard error
      out.weights.increase$`Std. Error`[out.weights.increase$variable %in% j] <-
        out.weights.increase$`Std. Error`[out.weights.increase$variable %in% j] +
        coef.se.list[[l]][which(rownames(coef.se.list[[l]]) %in% 'Increase'),colnames(coef.se.list[[l]]) %in% j] * 
        out.weights.increase[out.weights.increase$variable %in% j,paste0('mod_',l)]
      
      # Decreasing models
      # Estimate
      out.weights.decrease$Estimate[out.weights.decrease$variable %in% j] <-
        out.weights.decrease$Estimate[out.weights.decrease$variable %in% j] +
        coef.multinom.list[[l]][which(rownames(coef.multinom.list[[l]]) %in% 'Decrease'),colnames(coef.multinom.list[[l]]) %in% j] * 
        out.weights.decrease[out.weights.decrease$variable %in% j,paste0('mod_',l)]
      # Standard error
      out.weights.decrease$`Std. Error`[out.weights.decrease$variable %in% j] <-
        out.weights.decrease$`Std. Error`[out.weights.decrease$variable %in% j] +
        coef.se.list[[l]][which(rownames(coef.se.list[[l]]) %in% 'Decrease'),colnames(coef.se.list[[l]]) %in% j] * 
        out.weights.decrease[out.weights.decrease$variable %in% j,paste0('mod_',l)]
      
      
    }
  }
  
  # Checking for NAs in standard errors...
  # Increase first
  if(sum(is.na(out.weights.increase$`Std. Error`)) != 0) {
    # Which variables?
    variable.list <- as.character(out.weights.increase$variable[is.na(out.weights.increase$`Std. Error`)])
    
    # Looping through these to first reweight
    for(j in variable.list) {
      mods.check <- which(out.weights.increase[out.weights.increase$variable %in% j,which(names(out.weights.increase) %in% 'mod_1') : which(names(out.weights.increase) %in% 'mod_5')] > 0)
      # Looping through models
      for(mod in mods.check) {
        # Updating weighting for this variable
        
        # If st error not NaN
        if(!is.na(coef.se.list[[mod]][which(rownames(coef.se.list[[mod]]) %in% 'Increase'), colnames(coef.se.list[[l]]) %in% j])) {
          # Do nothing
        } else { # Change to 0 if NaN
          out.weights.increase[out.weights.increase$variable %in% j, paste0('mod_',mod)] <- 0
        }
      }
      # Reweighting
      out.weights.increase[out.weights.increase$variable %in% j,which(names(out.weights.increase) %in% 'mod_1') : which(names(out.weights.increase) %in% 'mod_5')] <-
        out.weights.increase[out.weights.increase$variable %in% j,which(names(out.weights.increase) %in% 'mod_1') : which(names(out.weights.increase) %in% 'mod_5')] /
        sum(out.weights.increase[out.weights.increase$variable %in% j,which(names(out.weights.increase) %in% 'mod_1') : which(names(out.weights.increase) %in% 'mod_5')])
      
      # And recalculating standard error
      out.weights.increase$`Std. Error`[out.weights.increase$variable %in% j] <- 0
      for(mod in mods.check) {
        
        if(!is.na(coef.se.list[[mod]][which(rownames(coef.se.list[[mod]]) %in% 'Increase'),colnames(coef.se.list[[mod]]) %in% j])) {
          # If not NaN
          out.weights.increase$`Std. Error`[out.weights.increase$variable %in% j] <-
            out.weights.increase$`Std. Error`[out.weights.increase$variable %in% j] +
            coef.se.list[[mod]][which(rownames(coef.se.list[[mod]]) %in% 'Increase'),colnames(coef.se.list[[mod]]) %in% j] * 
            out.weights.increase[out.weights.increase$variable %in% j,paste0('mod_',l)]
        } else {
          # Nothing
        }
      }
    }
  }
  
  # Repeating for decrease
  if(sum(is.na(out.weights.decrease$`Std. Error`)) != 0) {
    # Which variables?
    variable.list <- as.character(out.weights.decrease$variable[is.na(out.weights.decrease$`Std. Error`)])
    
    # Looping through these to first reweight
    for(j in variable.list) {
      mods.check <- which(out.weights.decrease[out.weights.decrease$variable %in% j,which(names(out.weights.decrease) %in% 'mod_1') : which(names(out.weights.decrease) %in% 'mod_5')] > 0)
      # Looping through models
      for(mod in mods.check) {
        # Updating weighting for this variable
        
        # If st error not NaN
        if(!is.na(coef.se.list[[mod]][which(rownames(coef.se.list[[mod]]) %in% 'Decrease'), colnames(coef.se.list[[l]]) %in% j])) {
          # Do nothing
        } else { # Change to 0 if NaN
          out.weights.decrease[out.weights.decrease$variable %in% j, paste0('mod_',mod)] <- 0
        }
      }
      # Reweighting
      out.weights.decrease[out.weights.decrease$variable %in% j,which(names(out.weights.decrease) %in% 'mod_1') : which(names(out.weights.decrease) %in% 'mod_5')] <-
        out.weights.decrease[out.weights.decrease$variable %in% j,which(names(out.weights.decrease) %in% 'mod_1') : which(names(out.weights.decrease) %in% 'mod_5')] /
        sum(out.weights.decrease[out.weights.decrease$variable %in% j,which(names(out.weights.decrease) %in% 'mod_1') : which(names(out.weights.decrease) %in% 'mod_5')])
      
      # And recalculating standard error
      out.weights.decrease$`Std. Error`[out.weights.decrease$variable %in% j] <- 0
      for(mod in mods.check) {
        
        if(!is.na(coef.se.list[[mod]][which(rownames(coef.se.list[[mod]]) %in% 'Decrease'),colnames(coef.se.list[[mod]]) %in% j])) {
          # If not NaN
          out.weights.decrease$`Std. Error`[out.weights.decrease$variable %in% j] <-
            out.weights.decrease$`Std. Error`[out.weights.decrease$variable %in% j] +
            coef.se.list[[mod]][which(rownames(coef.se.list[[mod]]) %in% 'Decrease'),colnames(coef.se.list[[mod]]) %in% j] * 
            out.weights.decrease[out.weights.decrease$variable %in% j,paste0('mod_',l)]
        } else {
          # Nothing
        }
      }
    }
  }
  
  
  # Getting accuracy again
  # Now for the weighted model
  # Getting accuracy and plotting ----
  coef.multinom.mean <- matrix(rep(0,nrow(out.weights.increase)*2), nrow = 2) %>% as.data.frame(.)
  names(coef.multinom.mean) <- out.weights.increase$variable
  rownames(coef.multinom.mean) <- rownames(coef.multinom.list[[1]])
  
  # Updating layout of the coefficients
  for(j in colnames(coef.multinom.mean)) {
    # Updating increase coefficient
    coef.multinom.mean[rownames(coef.multinom.mean) %in% 'Increase',colnames(coef.multinom.mean) %in% j] <-
      out.weights.increase$Estimate[out.weights.increase$variable %in% j]
    
    # Updating decrease coefficient
    coef.multinom.mean[rownames(coef.multinom.mean) %in% 'Decrease',colnames(coef.multinom.mean) %in% j] <-
      out.weights.decrease$Estimate[out.weights.decrease$variable %in% j]
  }
  
  # List of iso3s for accuracy
  iso3s <-
    sort(unique(dat_mod$ISO_numeric))
  
  # Gtting variable list
  vars = names(coef.multinom.mean)[2:(min(grep('factor',names(coef.multinom.mean))) - 1)]
  
  # Looping through iso3s
  # Cannot use R's predict function because not every country is in the data frame used to predict
  # Output data frame
  dat.out <- data.frame()
  for(i in iso3s) {
    pred_tmp <- multinom.predict.2b.trial(dat = dat_mod %>% filter(ISO_numeric %in% i),
                                          coefs = coef.multinom.mean,
                                          vars = vars,
                                          factor_vars = 'ISO_numeric')
    dat.out <- 
      rbind(dat.out,
            dat_mod %>% filter(ISO_numeric %in% i) %>% cbind(pred_tmp))
  }
  
  # And checking this works really quick
  # Not that this only works when k == 1
  # preds.check <-
  #   predict(objModel,
  #           dat_mod %>% filter(ISO_numeric %in% 332),
  #           'prob')
  # 
  # preds.check <-
  #   cbind(dat.out %>% filter(ISO_numeric %in% 332),
  #         preds.check)
  
  
  dat.out <-
    dat.out %>%
    transform(has_changed_pasture = factor(has_changed_pasture, levels = c('aNoChange','Decrease','Increase'))) %>%
    mutate(predicted_values = ifelse(aNoChange > Decrease & aNoChange > Increase,'aNoChange',
                                     ifelse(Decrease > Increase, 'Decrease','Increase'))) %>%
    transform(predicted_values = factor(predicted_values, levels = c('aNoChange','Decrease','Increase')))
  
  # Accuracy
  accuracy <- caret::confusionMatrix(data = dat.out$predicted_values,
                                     reference = dat.out$has_changed_pasture)
  
  # Making the plot look good
  # Colour palette
  my.palette = c('#ffffff','#007dff') %>%
    colorRampPalette(.)
  
  acc.sum <-
    dat.out %>%
    group_by(predicted_values, has_changed_pasture) %>%
    dplyr::summarise(Count = n()) %>%
    left_join(.,
              dat_mod %>% group_by(has_changed_pasture) %>% dplyr::summarise(Count_Actual = n())) %>%
    mutate(plot.x = ifelse(has_changed_pasture %in% 'aNoChange',2,
                           ifelse(has_changed_pasture %in% 'Decrease',1,3))) %>%
    mutate(plot.y = ifelse(predicted_values %in% 'aNoChange',2,
                           ifelse(predicted_values %in% 'Decrease',1,3))) %>%
    transform(predicted_values = factor(predicted_values, levels = c('aNoChange','Decrease','Increase'), ordered = TRUE)) %>%
    mutate(percent_predicted = Count/Count_Actual) %>%
    as.data.frame(.) %>%
    mutate(label = paste('Predicted: ',Count,'\n',
                         'Actual: ', Count_Actual,'\n',
                         'Percent Predicted: ',round(percent_predicted, digits = 2)*100,"%",
                         sep = ""))
  
  ggplot(acc.sum) +
    theme_bw() +
    geom_tile(aes(x = plot.x, y = plot.y, fill = round(percent_predicted * 100, digits = 2))) +
    geom_text(aes(x = plot.x, y = plot.y, label = label), size = 4.5) +
    scale_fill_gradientn(colours=my.palette(100), limits = c(0,100), expand = c(0,0)) +
    scale_colour_manual(values = c('black','white')) +
    scale_x_continuous(breaks = c(1,2,3),labels = (c('Decrease','No Change','Increase'))) +
    scale_y_continuous(breaks = c(1,2,3),labels = (c('Decrease','No Change','Increase'))) +
    coord_flip() +
    labs(x = 'Actual Pastureland Change', y = 'Predicted Pastureland Change', fill = 'Percent Accuracy') +
    theme(axis.text = element_text(size = 12)) +
    theme(axis.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12))
  
  
  # Creating directory for plots
  if('Land Forecasting Coefs and Plots' %in% list.files(getwd())) {
    # Nothing
  } else {
    dir.create(paste0(getwd(),'/Land Forecasting Coefs and Plots'))
    dir.create(paste0(getwd(),'/Land Forecasting Coefs and Plots/Model Coefficients'))
    dir.create(paste0(getwd(),'/Land Forecasting Coefs and Plots/Model Plots'))
    dir.create(paste0(getwd(),'/Land Forecasting Coefs and Plots/Multinomial Accuracy'))
  }
  Sys.sleep(5)
  
  ggsave(paste0(getwd(),'/Land Forecasting Coefs and Plots/Model Plots/Pasture Multinomial Accuracy 29July2021 ',
                region_list[k],'.pdf'), width = 12, height = 10)
  
  # %>%
  #   mutate(label = paste('Predicted: ',Count,'\n',
  #                        'Actual: ', Count_Actual,'\n',
  #                        'Percent Predicted: ',round(percent_predicted, digits = 2)*100,"%",
  #                        sep = ""))
  
  
  
  
  # 
  # coef.prob[[region_list[k]]] <- as.data.frame(exp(coef(objModel)))
  # mods.prob[[region_list[k]]] <- objModel
  # 
  # # Changing outcome name for extent models
  # outcomeName <- "prop_change_cell_crop_2007.2012"
  # Processing: Creating data frames that only contains cells which have increased or decreasd in extent
  dat_region_increase <- dat_mod[dat_mod$has_changed_pasture %in% 'Increase',]
  dat_region_decrease <- dat_mod[dat_mod$has_changed_pasture %in% 'Decrease',]
  
  # removing data
  rm(dat_mod)
  
  # Dropping rows with values <0 for increaes and > 0 for decrease
  dat_region_increase <- dat_region_increase[dat_region_increase$prop_change_cell_pasture_2005.2010 > 0,]
  dat_region_decrease <- dat_region_decrease[dat_region_decrease$prop_change_cell_pasture_2005.2010 < 0,]
  
  # arcsining response variable
  # dat_region_increase$prop_change_cell_crop_2007.2012 <- asin(sqrt(dat_region_increase$prop_change_cell_crop_2007.2012))
  # dat_region_decrease$prop_change_cell_crop_2007.2012 <- asin(sqrt(abs(dat_region_decrease$prop_change_cell_crop_2007.2012)))
  
  # Creating groups for testing and training
  # Filtering to countries with more than 5 observations
  # Getting count by iso3n
  dat.ison <- 
    table(dat_region_increase$ISO_numeric) %>%
    .[. >= 5] %>% names(.)
  
  # Filtering countries with fewer than 5 data points
  dat_region_increase <-
    dat_region_increase %>%
    filter(ISO_numeric %in% dat.ison)
  
  dat_region_increase$group <-
    kfold(dat_region_increase, k = 5, by = dat_region_increase$ISO_numeric)
  
  # list of RMSE to weight models
  rmse.list <- list()
  coef.list <- list()
  
  # Looping through these to weight coefficients
  for (l in 1:5) {
    # Training dat
    train.dat <-
      dat_region_increase %>%
      filter(!(group %in% l))
    # Testing dat
    test.dat <-
      dat_region_increase %>%
      filter(group %in% l)
    # Getting model for extent increase
    objModel_increase <- glm(prop_change_cell_pasture_2005.2010 ~ 
                               dist_50k_log +
                               dist_50k_log_sq +
                               # ag_suitability +
                               PA_binary +
                               prop_to_predict_from_crop +
                               prop_crop_sq +
                               prop_to_predict_from_pasture +
                               prop_pasture_sq +
                               prop_adj_to_predict_from_pasture +
                               prop_adj_pasture_sq +
                               prop_to_predict_from_cell_delta_pasture +
                               prop_cell_delta_pasture_sq +
                               factor(ISO_numeric, ordered = FALSE),
                             family = Gamma(link = 'log'),
                             dat = train.dat,
                             model = FALSE,
                             control = list(maxit = 250))
    
    coef.list[[l]] <- as.data.frame(summary(objModel_increase)$coefficients)
    
    # Checking to see if there are any missing factor variables
    # ANd if so, setting these factor variables to 0
    isos.have <- sort(unique(train.dat$ISO_numeric))
    isos.need <- sort(unique(test.dat$ISO_numeric))
    
    # Checking to see if we need to update
    # isos.need <- isos.need[!(isos.need %in% isos.have)]
    # if(length(isos.needs > 0)) {
    #   # Coefficients to add
    #   tmp.add <- rep(0,length(isos.need))
    #   # Giving names
    #   names(tmp.add) <- paste0('factor(ISO_numeric, ordered = FALSE)',isos.need)
    #   objModel_increase$coefficients <-
    #     c(objModel_increase$coefficients,
    #       tmp.add)
    # } else {
    #   # Do nothing
    # }
    
    # Getting rmse
    rmse.list[[l]] <- RMSE(obs = test.dat$prop_change_cell_pasture_2005.2010[test.dat$ISO_numeric %in% train.dat$ISO_numeric], # Observed values
                           pred = exp(predict(objModel_increase, test.dat %>% filter(ISO_numeric %in% train.dat$ISO_numeric)))) # Predicted values
    
    
    # # Testing predictions to maake sure function works
    # predictor_vars_pasture <- 
    #   rownames(coef.list[[l]]) %>%
    #   .[-grep('Intercept',.)] %>%
    #   .[-grep('ISO',.)] %>%
    #   as.character(.) %>%
    #   c(.,'ISO_numeric')
    # # 
    # tmp.coefs = coef.list[[l]] %>% mutate(variable = rownames(.))
    # dat.out = data.frame()
    # 
    # iso3s = levels(test.dat$ISO_numeric)
    # for(i in iso3s) {
    #   dat.out <- rbind(dat.out,
    #                    test.dat %>% filter(ISO_numeric %in% i) %>%
    #                      mutate(pred_change_cell_prop_crop_decrease = 
    #                               glm.predict.2b(dat = test.dat %>% filter(ISO_numeric %in% i),
    #                                               coefs = tmp.coefs,
    #                                               vars = predictor_vars_pasture[1:(length(predictor_vars_pasture) - 1)],
    #                                               factor_vars = 'ISO_numeric')))
    #   
    #   
    # }
    # 
    # # Predicting normal values
    # test.dat$preds <- exp(predict(objModel_increase, test.dat))
    
    # These work fine
    
    
  }
  
  # Weighting models based on RMSE
  out.weights <- 
    data.frame(variable = unique(c(rownames(coef.list[[1]]),
                                   rownames(coef.list[[2]]),
                                   rownames(coef.list[[3]]),
                                   rownames(coef.list[[4]]),
                                   rownames(coef.list[[5]]))),
               mod_1 = NA, mod_2 = NA, mod_3 = NA, mod_4 = NA, mod_5 = NA)
  
  # Looping through to calculate weighted coefficients
  for(l in 1:5) {
    # Making data frame saying which variables are in which model
    out.weights[out.weights$variable %in% rownames(coef.list[[l]]),paste0('mod_',l)] <- rmse.list[[l]]
  }
  # Weighting
  
  out.weights[,which(names(out.weights) %in% 'mod_1') : which(names(out.weights) %in% 'mod_5')] <-
    (1 - out.weights[,which(names(out.weights) %in% 'mod_1') : which(names(out.weights) %in% 'mod_5')]^2) /
    rowSums((1 - out.weights[,which(names(out.weights) %in% 'mod_1') : which(names(out.weights) %in% 'mod_5')]^2), na.rm = TRUE)
  
  out.weights[is.na(out.weights)] <- 0
  out.weights$Estimate <- 0
  out.weights$`Std. Error` <- 0
  
  # And looping again
  for(l in 1:5) {
    for(j in rownames(coef.list[[l]])) {
      out.weights$Estimate[out.weights$variable %in% j] <-
        out.weights$Estimate[out.weights$variable %in% j] +
        coef.list[[l]]$Estimate[rownames(coef.list[[l]]) %in% j] * 
        out.weights[out.weights$variable %in% j,paste0('mod_',l)]
      
      out.weights$`Std. Error`[out.weights$variable %in% j] <-
        out.weights$`Std. Error`[out.weights$variable %in% j] +
        coef.list[[l]]$`Std. Error`[rownames(coef.list[[l]]) %in% j] * 
        out.weights[out.weights$variable %in% j,paste0('mod_',l)]
    }
  }
  
  
  # testing something
  # RMSE(obs = dat_region_increase$prop_change_cell_pasture_2005.2010,
  #      pred = exp(predict(objModel_increase, dat_region_increase)))
  # 
  objModel_increase$coefficients[1:nrow(out.weights)] <- out.weights$Estimate
  # 
  # RMSE(obs = dat_region_increase$prop_change_cell_pasture_2005.2010,
  #      pred = exp(predict(objModel_increase, dat_region_increase)))
  
  coef.increase.mean.se <- out.weights
  coef.list.increase <- coef.list
  rmse.list.increase <- rmse.list
  
  
  # Repeating for decrease
  # Creating groups for testing and training
  # Filtering to countries with more than 5 observations
  # Getting count by iso3n
  dat.ison <- 
    table(dat_region_decrease$ISO_numeric) %>%
    .[. >= 5] %>% names(.)
  
  # Filtering countries with fewer than 5 data points
  dat_region_decrease <-
    dat_region_decrease %>%
    filter(ISO_numeric %in% dat.ison)
  # And creating groups
  dat_region_decrease$group <-
    kfold(dat_region_decrease, k = 5, by = dat_region_decrease$ISO_numeric)
  
  # Getting absolute of decrease in crop extent
  dat_region_decrease$prop_change_cell_pasture_2005.2010 <- abs(dat_region_decrease$prop_change_cell_pasture_2005.2010)
  
  # list of RMSE to weight models
  rmse.list <- list()
  coef.list <- list()
  
  # Looping through these to weight coefficients
  for (l in 1:5) {
    # Training dat
    train.dat <-
      dat_region_decrease %>%
      filter(!(group %in% l))
    # Testing dat
    test.dat <-
      dat_region_decrease %>%
      filter(group %in% l)
    # Getting model for extent increase
    objModel_decrease <- glm(prop_change_cell_pasture_2005.2010 ~ 
                               dist_50k_log +
                               dist_50k_log_sq +
                               # ag_suitability +
                               PA_binary +
                               prop_to_predict_from_crop +
                               prop_crop_sq +
                               prop_to_predict_from_pasture +
                               prop_pasture_sq +
                               prop_adj_to_predict_from_pasture +
                               prop_adj_pasture_sq +
                               prop_to_predict_from_cell_delta_pasture +
                               prop_cell_delta_pasture_sq +
                               factor(ISO_numeric, ordered = FALSE),
                             family = Gamma(link = 'log'),
                             dat = train.dat,
                             model = FALSE,
                             control = list(maxit = 250))
    
    # Saving in a list
    coef.list[[l]] <- as.data.frame(summary(objModel_decrease)$coefficients)
    
    # Getting rmse
    rmse.list[[l]] <- RMSE(obs = test.dat$prop_change_cell_pasture_2005.2010[test.dat$ISO_numeric %in% train.dat$ISO_numeric], # Observed values
                           pred = exp(predict(objModel_decrease, test.dat %>% filter(ISO_numeric %in% train.dat$ISO_numeric))))
    
    
  }
  
  # Weighting models based on RMSE
  out.weights <- 
    data.frame(variable = unique(c(rownames(coef.list[[1]]),
                                   rownames(coef.list[[2]]),
                                   rownames(coef.list[[3]]),
                                   rownames(coef.list[[4]]),
                                   rownames(coef.list[[5]]))),
               mod_1 = NA, mod_2 = NA, mod_3 = NA, mod_4 = NA, mod_5 = NA)
  
  # Looping through to calculate weighted coefficients
  for(l in 1:5) {
    # Making data frame saying which variables are in which model
    out.weights[out.weights$variable %in% rownames(coef.list[[l]]),paste0('mod_',l)] <- rmse.list[[l]]
  }
  # Weighting
  
  out.weights[,which(names(out.weights) %in% 'mod_1') : which(names(out.weights) %in% 'mod_5')] <-
    (1 - out.weights[,which(names(out.weights) %in% 'mod_1') : which(names(out.weights) %in% 'mod_5')]^2) /
    rowSums((1 - out.weights[,which(names(out.weights) %in% 'mod_1') : which(names(out.weights) %in% 'mod_5')]^2), na.rm = TRUE)
  
  out.weights[is.na(out.weights)] <- 0
  out.weights$Estimate <- 0
  out.weights$`Std. Error` <- 0
  
  # And looping again
  for(l in 1:5) {
    for(j in rownames(coef.list[[l]])) {
      out.weights$Estimate[out.weights$variable %in% j] <-
        out.weights$Estimate[out.weights$variable %in% j] +
        coef.list[[l]]$Estimate[rownames(coef.list[[l]]) %in% j] * 
        out.weights[out.weights$variable %in% j,paste0('mod_',l)]
      
      out.weights$`Std. Error`[out.weights$variable %in% j] <-
        out.weights$`Std. Error`[out.weights$variable %in% j] +
        coef.list[[l]]$`Std. Error`[rownames(coef.list[[l]]) %in% j] * 
        out.weights[out.weights$variable %in% j,paste0('mod_',l)]
    }
  }
  
  
  # Checking weighted RMSE
  # Predicting
  # testing something
  # RMSE(obs = dat_region_decrease$prop_change_cell_pasture_2005.2010,
  #      pred = exp(predict(objModel_decrease, dat_region_decrease)))
  # 
  # objModel_decrease$coefficients[1:nrow(out.weights)] <- out.weights$Estimate
  # 
  # RMSE(obs = dat_region_decrease$prop_change_cell_pasture_2005.2010,
  #      pred = exp(predict(objModel_decrease, dat_region_decrease)))
  
  coef.decrease.mean.se <- out.weights
  coef.list.decrease <- coef.list
  rmse.list.decrease <- rmse.list
  
  
  # Saving coefficients
  # write.csv(coef.multinom.mean,paste0(getwd(),'/Land Forecasting Coefs and Plots/Model Coefficients/',region_list[k],' Pasture Multinomial Mean Coefficients.csv'), row.names = FALSE)
  # write.csv(coef.multinom.se,paste0(getwd(),'/Land Forecasting Coefs and Plots/Model Coefficients/',region_list[k],' Pasture Multinomial se Coefficients.csv'), row.names = FALSE)
  # write.csv(coef.increase.mean.se,paste0(getwd(),'/Land Forecasting Coefs and Plots/Model Coefficients/',region_list[k],' Pasture GLM Increase Coefficients.csv'), row.names = FALSE)
  # write.csv(coef.decrease.mean.se,paste0(getwd(),'/Land Forecasting Coefs and Plots/Model Coefficients/',region_list[k],' Pasture GLM Decrease Coefficients.csv'), row.names = FALSE)
  
  # And saving multinomial accuracy
  save(file = paste0(getwd(),'/Land Forecasting Coefs and Plots/Multinomial Accuracy/',region_list[k],' Pasture Multinomial Accuracy 29July2021.RData'),
       accuracy)
  
  save(file = paste0(getwd(),'/Land Forecasting Coefs and Plots/Model Coefficients/',region_list[k],' Pasture Model Coefficients 29July2021.RData'),
       coef.multinom.mean, coef.acc.list, coef.increase.mean.se, coef.decrease.mean.se, coef.multinom.se, coef.list.increase, coef.list.decrease, rmse.list.increase, rmse.list.decrease)
  
  
}

Sys.time() - t1


