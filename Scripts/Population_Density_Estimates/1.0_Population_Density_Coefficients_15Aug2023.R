#!/usr/bin/env Rscript

### Refitting a single model to the most important variables in the Santini dataset
# This estimates population density of terrestrial vertebrates based on 
# life history traits
# and climate conditions

# Installing packages
library(raster)
library(blmeco)
library(data.table)
library(lme4)
library(MuMIn)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(Hmisc)
library(stringr)


# System options... 
options(tibble.width = Inf)
rm(list=ls())

# Setting working directory
setwd('/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity')

# Only need to do this first part once
# Putting an if statement here to save time
if('SantiniModellingData_12May2023.csv' %in% list.files(paste0(getwd(),'/Other Data Inputs/Pop Density Inputs'))) { # if it exists
  # Do all the below
  # Load density data from original data frame ----
  # tetra_data <- read_csv(paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/TetraDENSITY Copy2.csv")) %>%
  tetra_data <- read_csv(paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/TetraDENSITY_v.1.csv")) %>%
  # tetra_data <- read_csv("/Users/macuser/Desktop/Pop_Density_Inputs/TetraDENSITY_v.1.csv") %>%
    filter(Class %in% c("Mammalia", 'Aves', 'Amphibia','Reptilia')) %>%
    mutate(density_corrected =
             ifelse(Density_unit %in% 'ind/ha', Density * 100,
                    ifelse(Density_unit %in% 'males/ha', Density * 2 * 100,
                           ifelse(Density_unit %in% 'pairs/km2', Density * 2, Density)))) %>%
    mutate(Density = density_corrected) %>%
    dplyr::select(.,Class:Species, Longitude, Latitude, Density) %>%
    filter(complete.cases(.))

  names(tetra_data) <- tolower(names(tetra_data))
  
  # Get coordinates for spatial data
  coords <- tetra_data %>%
    dplyr::select(.,x = longitude,
                  y = latitude)
  
  # Load explanatory variables and bind to density data ----
  body_mass_data <- read.csv(paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/Body Mass Estimates 14November2019.csv"), stringsAsFactors = FALSE)
  # body_mass_data <- read.csv("/Users/macuser/Desktop/Pop_Density_Inputs/Body Mass Estimates 14November2019.csv", stringsAsFactors = FALSE)
  name_conversion <- read.csv(paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/SantiniNameConversion.csv"), stringsAsFactors = FALSE)
  # name_conversion <- read.csv("/Users/macuser/Desktop/Pop_Density_Inputs/SantiniNameConversion.csv", stringsAsFactors = FALSE)
  
  
  # Importing data for current climate ----
  temp_mean <- raster(paste0(getwd(),'/CMIP6_Climate_Data/Historic/1995-2014/Managed_Rasters/Mean_annual_temp.tif'))
  precip_var <- raster(paste0(getwd(),'/CMIP6_Climate_Data/Historic/1995-2014/Managed_Rasters/Variance_of_precipitation.tif'))
  precip_warmest <- raster(paste0(getwd(),'/CMIP6_Climate_Data/Historic/1995-2014/Managed_Rasters/Precip_in_warmest_quarter.tiff'))
  temp_var <- raster(paste0(getwd(),'/CMIP6_Climate_Data/Historic/1995-2014/Managed_Rasters/Variance_of_temperature.tif'))
  npp <- raster(paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/npp_santini.tif"))
  
  #temp_mean <- raster(paste0('/Users/macuser/Desktop/Historic_Climate_Data/Mean_annual_temp.tif'))
  #precip_var <- raster(paste0('/Users/macuser/Desktop/Historic_Climate_Data/Variance_of_precipitation.tif'))
  #precip_warmest <- raster('/Users/macuser/Desktop/Historic_Climate_Data/Precip_in_warmest_quarter.tiff')
  #temp_var <- raster('/Users/macuser/Desktop/Historic_Climate_Data/Variance_of_temperature.tif')
  #npp <- raster("/Users/macuser/Desktop/Pop_Density_Inputs/npp_santini.tif")
  
  # Richness rasters ----
  # These are based off AOH maps, but can be recreated with other data inputs as needed
  mam.richness <- raster(paste0(getwd(),"/Other Data Inputs/Species Richness Rasters/Mammals_Richness.tif"))
  bird.richness <- raster(paste0(getwd(),"/Other Data Inputs/Species Richness Rasters/Birds_Richness.tif"))
  amp.richness <- raster(paste0(getwd(),"/Other Data Inputs/Species Richness Rasters/Amphibians_Richness.tif"))
 rep.richness <- raster(paste0(getwd(),"/Other Data Inputs/Species Richness Rasters/Reptile_Richness.tif")) 
  # For local computer...
  #mam.richness <- raster("/Users/macuser/Desktop/Species_Richness_Rasters/Mammals_Richness.tif")
  #bird.richness <- raster("/Users/macuser/Desktop/Species_Richness_Rasters/Birds_Richness.tif")
  #amp.richness <- raster("/Users/macuser/Desktop/Species_Richness_Rasters/Amphibians_Richness.tif")
  
  mam.richness <- raster::aggregate(mam.richness,15)
  bird.richness <- raster::aggregate(bird.richness,15)
  amp.richness <- raster::aggregate(amp.richness,15)
  rep.richness <- raster::aggregate(rep.richness,15)
  
  richness = mam.richness + bird.richness + amp.richness + rep.richness
  
  # Resampling rasters to appropriate resolution and CRS ----
  mam.richness <- round(projectRaster(mam.richness, npp, method = 'bilinear', na.rm = TRUE), digits = 0)
  amp.richness <- round(projectRaster(amp.richness, npp, method = 'bilinear', na.rm = TRUE), digits = 0)
  bird.richness <- round(projectRaster(bird.richness, npp, method = 'bilinear', na.rm = TRUE), digits = 0)
  rep.richness <- round(projectRaster(rep.richness, npp, method = 'bilinear', na.rm = TRUE), digits = 0)
  richness <- round(projectRaster(richness, npp, method = 'bilinear', na.rm = TRUE), digits = 0)
  
  temp_mean <- projectRaster(temp_mean, npp, method = 'bilinear', na.rm = TRUE)
  temp_var <- projectRaster(temp_var, npp, method = 'bilinear', na.rm = TRUE)
  precip_var <- projectRaster(precip_var, npp, method = 'bilinear', na.rm = TRUE)
  precip_warmest <- projectRaster(precip_warmest, npp, method = 'bilinear', na.rm = TRUE)
  
  # Manipulating data ----
  tetra_data <- tetra_data %>%
    mutate_at(.vars = c("class", "order", "family"), toupper) %>%
    mutate(santini_binomial = paste(genus, species, sep = "_")) %>%
    left_join(., name_conversion) %>%
    mutate(binomial = ifelse(is.na(iucn_binomial),
                             yes = santini_binomial,
                             no = iucn_binomial)) %>%
    left_join(., body_mass_data %>%
                dplyr::select(binomial, est_mass_kg)) %>%
    # Need to manually add BM estimates for a couple of species
    # mutate(combined_body_mass_g = case_when(binomial == "Camelus_dromedarius" ~ 488000/1000,
    #                                binomial == "Alouatta_seniculus" ~ 6398.31/1000,
    #                                binomial == "Cebus_capucinus" ~ 3010/1000,
    #                                binomial == "Pithecia_monachus" ~ 2110/1000,
    #                                TRUE ~ est_mass_kg)) %>%
    mutate(combined_body_mass_g = est_mass_kg * 1000) %>%
    mutate(log10_density = log10(density),
           log10_body_mass_g = log10(combined_body_mass_g),
           npp = raster::extract(npp, coords),
           pcv = raster::extract(precip_var, coords),
           pwarmest = raster::extract(precip_warmest, coords),
           tvar = raster::extract(temp_var, coords),
           tmean = raster::extract(temp_mean, coords),
           richness_mammals = raster::extract(mam.richness, coords),
           richness_birds = raster::extract(bird.richness, coords),
           richness_amps = raster::extract(amp.richness, coords),
	   richness_reps = raster::extract(rep.richness, coords),
           richness = raster::extract(richness, coords)) %>%
    dplyr::select(.,class, order, family, genus, species, binomial,
                  density, log10_density,
                  combined_body_mass_g, log10_body_mass_g,
                  longitude, latitude, npp, pcv, pwarmest, tvar, tmean,
                  richness_mammals, richness_birds, richness_amps, richness_reps, richness) %>%
    filter(complete.cases(.)) %>%
    mutate(richness_mammals = ceiling(richness_mammals),
           richness_birds = ceiling(richness_birds),
           richness_amps = ceiling(richness_amps),
	   richness_reps = ceiling(richness_reps),
           richness = ceiling(richness))
  
  # Updating some estimates manually
  tetra_data$combined_body_mass_g[tetra_data$binomial == "Camelus_dromedarius"] <- 488000/1000
  tetra_data$combined_body_mass_g[tetra_data$binomial == "Alouatta_seniculus"] <- 6398.31/1000
  tetra_data$combined_body_mass_g[tetra_data$binomial == "Cebus_capucinus"] <- 3010/1000
  tetra_data$combined_body_mass_g[tetra_data$binomial == "Pithecia_monachus"] <- 2110/1000
  
  # mutate(combined_body_mass_g = case_when(binomial == "Camelus_dromedarius" ~ 488000/1000,
  #                                binomial == "Alouatta_seniculus" ~ 6398.31/1000,
  #                                binomial == "Cebus_capucinus" ~ 3010/1000,
  #                                binomial == "Pithecia_monachus" ~ 2110/1000,
  #                                TRUE ~ est_mass_kg)) %>%
  
  tetra_data$richness_amps[tetra_data$richness_amps %in% -1] <- 0
  tetra_data$richness_birds[tetra_data$richness_birds %in% -1] <- 0
  tetra_data$richness_mammals[tetra_data$richness_mammals %in% -1] <- 0
  tetra_data$richness_reps[tetra_data$richness_reps %in% -1] <- 0
  
  tetra_data$tot_richness <-
    tetra_data$richness
  
  tetra_data <-
    tetra_data %>% filter(tot_richness >= 1)
  
  tetra_data$richness_amps[tetra_data$richness_amps == 0] <- NA
  tetra_data$richness_mammals[tetra_data$richness_mammals == 0] <- NA
  tetra_data$richness_birds[tetra_data$richness_birds == 0] <- NA
  tetra_data$richness_reps[tetra_data$richness_reps == 0] <- NA

  tetra_data <-
    tetra_data %>%
    dplyr::rename(richness_amphibia = richness_amps,
                  richness_mammalia = richness_mammals,
                  richness_aves = richness_birds,
		  richness_reptilia = richness_reps)
  
  write_csv(tetra_data,
             paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/SantiniModellingData_12May2023.csv"))
  
  # write_csv(tetra_data, "/Volumes/Citadel/Oxford/Research Projects/Multiple Drivers of Biodiversity/Pop Density Data Sets/SantiniModellingData_12May2023.csv")
  
  
} else { # if it doesn't exist
  # Do nothing
}


# Other functions
vif.lme <- function (fit) {
  ## adapted from rms::vif
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)] }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v 
}

# Function to extract random effects
ran_ef_fun <- function(l){
  for(i in 1:length(l)){
    name <- gsub(":.*$", "", names(l)[i])
    l[[i]][[name]]<- gsub(":.*$", "", row.names(l[[i]]))
    row.names(l[[i]]) <- NULL
    names(l[[i]])[1] <- paste(name, "intercept_adj", sep = "_")
  }
  return(l)
}

# Function to weight across models for random effects
weight_fun <- function(df){
  # df$w_tmp <- df[,names(df)[grepl('intercept_adj',names(df))]] * df$weight

  df$w_tmp <- df[,2] * df$weight

  weight_col <- names(df)[!grepl('model_number|intercept_adj|weight|w_tmp',names(df))]

  # converting to data frame to aggregate - throwing fits with data table =\
  df <- as.data.frame(df)

  df_new <-
	  rowsum(df$w_tmp, 
		 #group = df[,3])
		# df$binomial)
		  group = df[,weight_col])
  

  # changing back into a data frame...
  df_new <- 
	  as.data.frame(df_new) %>%
	  mutate(weight_col = row.names(.))

  # updating names
  names(df_new)[names(df_new) %in% 'V1'] <- paste0(weight_col, '_', 'intercept_adj_weighted')
  names(df_new)[names(df_new) %in% 'weight_col'] <- paste0(weight_col)

  # and reordering columns
  df_new <- df_new[,2:1]

  # returning
  return(df_new)
}


###
# Making the function that runs the models
# Output is a weighted model ensemble, with coefficients thereof
# Getting data
tetra_data <- read_csv(paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/SantiniModellingData_12May2023.csv"))
# tetra_data <- read_csv("/Volumes/Citadel/Oxford/Research Projects/Multiple Drivers of Biodiversity/Pop Density Data Sets/SantiniModellingData_12May2023.csv")

pop_density_function <-
  function(taxon) {
    # Filtering by taxon
    tmp_data <- tetra_data %>% filter(class %in% toupper(taxon))
    
    # Dropping columns and changing name of richness for the taxon
    # Dropping non-needed richness columns
    richness.drop <- names(tmp_data)[grepl('richness',names(tmp_data), ignore.case=TRUE)] %>% .[!grepl(taxon,.,ignore.case=TRUE)]
    tmp_data <- tmp_data[,!(names(tmp_data) %in% richness.drop)]
    tmp_data$richness <- unlist(tmp_data[,grep('richness',names(tmp_data),ignore.case=TRUE)])
    
    # And mutating columns
    tmp_data <- 
      tmp_data %>%
      filter(complete.cases(.)) %>% 
      mutate(class = tolower(class), family = tolower(family)) %>% 
      mutate(class = Hmisc::capitalize(class), family = Hmisc::capitalize(family))
    
    # fitting model
    # Exception for taxon - different coefficients
    if(taxon %in% c('aves','mammalia')) {
      m_max <-
        lmer(data = tmp_data, 
             formula = log10_density ~ 
               log10_body_mass_g + I(log10_body_mass_g^2) + I(log10_body_mass_g^3) +
               npp + I(npp^2) + 
               pcv + I(pcv^2) +
               pwarmest + I(pwarmest^2) +
               richness  +
               (1|order/family/binomial), 
             na.action = na.fail,
             REML = FALSE)
      
      # Fit models, ensuring that quadratics can only appear with linear effects
      m_subset <- expression(dc(log10_body_mass_g, `I(log10_body_mass_g^2)`, `I(log10_body_mass_g^3)`) &
                               dc(`npp`, `I(npp^2)`) &
                               dc(`pcv`, `I(pcv^2)`) &
                               dc(`pwarmest`, `I(pwarmest^2`) &
                               dc(`richness`))
    } else { # if amphibians
      m_max <- lmer(data = tmp_data, 
                    formula = log10_density ~ 
                      log10_body_mass_g + I(log10_body_mass_g^2) + I(log10_body_mass_g^3) +
                      npp + I(npp^2) + 
                      pcv + I(pcv^2) +
                      pwarmest + I(pwarmest^2) +
                      richness  +
                      (1|family), 
                    na.action = na.fail,
                    REML = FALSE)
      
      # Fit models, ensuring that quadratics can only appear with linear effects
      m_subset <- expression(dc(log10_body_mass_g, `I(log10_body_mass_g^2)`, `I(log10_body_mass_g^3)`) &
                               dc(`npp`, `I(npp^2)`) &
                               dc(`pcv`, `I(pcv^2)`) &
                               dc(`pwarmest`, `I(pwarmest^2`) &
                               dc(`richness`))
    }
    
    # Extracting coefficients
    dredge_table <- dredge(global.model = m_max, 
                           subset = m_subset,
                           beta = 'none',
                           evaluate = TRUE)
    
    # Extract models and coefficients and weight them -----
    averaged_models <- model.avg(dredge_table, subset = weight > .01)
    models <- get.models(dredge_table, subset = weight > 0.01)
    
    # Getting summary statistics across models
    # R2
    map(models, r.squaredGLMM) 
    # Pretty damn similar: ~0.3 for marginal; 0.77 for conditional
    map(models, dispersion_glmer)
    # Coefficients by model...
    map(models, vif.lme) 
    
    # Getting model weights
    weights <- data.frame(subset(dredge_table, weight > 0.01)) %>%
      dplyr::select(weight) %>%
      mutate(model_number = row.names(.))
    
    # Fixed effects ------
    fixed_effects <- model.avg(dredge_table, subset = weight > 0.001)
    fixed_effects <- data.frame(fixed_effects$coefficients)
    fixed_effects <- fixed_effects[row.names(fixed_effects) %in% 'full',]
    
    # Updating to make sure names are not getting mismatched
    names(fixed_effects)[names(fixed_effects) == 'X.Intercept.'] = 'intercept_weighted_coef'
    names(fixed_effects)[names(fixed_effects) == 'log10_body_mass_g'] = 'log10_body_mass_g_weighted_coef'
    names(fixed_effects)[names(fixed_effects) == 'I.log10_body_mass_g.2.'] = 'log10_body_mass_g2_weighted_coef'
    names(fixed_effects)[names(fixed_effects) == 'I.log10_body_mass_g.3.'] = 'log10_body_mass_g3_weighted_coef'
    names(fixed_effects)[names(fixed_effects) == 'npp'] = 'npp_weighted_coef'
    names(fixed_effects)[names(fixed_effects) == 'I.npp.2.'] = 'npp2_weighted_coef'
    names(fixed_effects)[names(fixed_effects) == 'pcv'] = 'pcv_weighted_coef'
    names(fixed_effects)[names(fixed_effects) == 'I.pcv.2.'] = 'pcv2_weighted_coef'
    names(fixed_effects)[names(fixed_effects) == 'pwarmest'] = 'pwarmest_weighted_coef'
    names(fixed_effects)[names(fixed_effects) == 'I.pwarmest.2.'] = 'pwarmest2_weighted_coef'
    names(fixed_effects)[names(fixed_effects) == 'richness'] = 'richness_weighted_coef'
    
    # Managing data frame
    row.names(fixed_effects) <- NULL
    fixed_effects <- fixed_effects %>%
      mutate(class = toupper(taxon))
    
    # Random effects -----
    random_effects <- map(models, ranef)
    random_effects <- map(random_effects, ran_ef_fun)
    
    # random_effects_weighted <- random_effects[[1]]
    random_effects_weighted <- list()
    for(i in 1:length(random_effects[[1]])){
      random_effects_weighted[[names(random_effects[[1]][i])]] <- rbindlist(map(random_effects, `[[`, names(random_effects[[1]][i])),
                                                                            idcol = "model_number")
    }
    # Getting random effects
    random_effects_weighted <- map(random_effects_weighted, left_join, weights)
    random_effects_weighted <- map(random_effects_weighted, weight_fun) 
    
    # Add coefficients to data to check predictions -----
    # Exception for taxon
    if(taxon %in% c('aves','mammalia')) { # birds and mammals
      tetra_check <- tmp_data %>%
        left_join(., fixed_effects %>% mutate(class = str_to_sentence(class))) %>%
        left_join(., random_effects_weighted[[1]]) %>%
        left_join(., random_effects_weighted[[2]]) %>%
        left_join(., random_effects_weighted[[3]]) %>%
        mutate(richness = richness, 
               richness_weighted_coef = richness_weighted_coef) %>%
        # dplyr::rename(binomial_intercept_adj_weighted = binomial_intercept_adj,
        #               family_intercept_adj_weighted = family_intercept_adj,
        #               order_intercept_adj_weighted = order_intercept_adj) %>%
        mutate(predicted_model = predict(model.avg(dredge_table, subset = weight > 0.001, fit = TRUE)),
               # mutate(predicted_model = predict(dredge_table %>% filter(weight > 0.1), fit = TRUE),
               predicted_coefficients = intercept_weighted_coef + 
                 binomial_intercept_adj_weighted + family_intercept_adj_weighted + order_intercept_adj_weighted + 
                 log10_body_mass_g * log10_body_mass_g_weighted_coef +
                 log10_body_mass_g^2 * log10_body_mass_g2_weighted_coef + 
                 log10_body_mass_g^3 * log10_body_mass_g3_weighted_coef + 
                 npp * npp_weighted_coef + 
                 npp^2 * npp2_weighted_coef + 
                 pcv * pcv_weighted_coef + 
                 pcv^2 * pcv2_weighted_coef + 
                 pwarmest * pwarmest_weighted_coef +
                 pwarmest^2 * pwarmest2_weighted_coef + 
                 richness * richness_weighted_coef,
               # coef_predict_dens = 10^predict_coefficients,
               check = predicted_coefficients - predicted_model)

	order_effects <- random_effects_weighted[[3]] %>%
      as.data.frame()
    } else { # amphibians
      tetra_check <- tmp_data %>%
        left_join(., fixed_effects %>% mutate(class = str_to_sentence(class))) %>%
        left_join(., random_effects_weighted[[1]]) %>%
        mutate(richness = richness, 
               richness_weighted_coef = richness_weighted_coef) %>%
        # dplyr::rename(binomial_intercept_adj_weighted = binomial_intercept_adj,
        #               family_intercept_adj_weighted = family_intercept_adj,
        #               order_intercept_adj_weighted = order_intercept_adj) %>%
        mutate(predicted_model = predict(model.avg(dredge_table, subset = weight > 0.001, fit = TRUE)),
               # mutate(predicted_model = predict(dredge_table %>% filter(weight > 0.1), fit = TRUE),
               predicted_coefficients = intercept_weighted_coef + 
                 family_intercept_adj_weighted +
                 log10_body_mass_g * log10_body_mass_g_weighted_coef +
                 log10_body_mass_g^2 * log10_body_mass_g2_weighted_coef + 
                 log10_body_mass_g^3 * log10_body_mass_g3_weighted_coef + 
                 npp * npp_weighted_coef + 
                 npp^2 * npp2_weighted_coef + 
                 pcv * pcv_weighted_coef + 
                 pcv^2 * pcv2_weighted_coef + 
                 pwarmest * pwarmest_weighted_coef +
                 pwarmest^2 * pwarmest2_weighted_coef + 
                 richness * richness_weighted_coef,
               # coef_predict_dens = 10^predict_coefficients,
               check = predicted_coefficients - predicted_model)
    }
    
    
    cat(taxon,':',range(round(tetra_check$check, digits = 10), na.rm = TRUE),'\n')
    # OK! All sorted.
    
    # Some orders don't have any data in the Tetra database, so model these as something else. -----
    # AFROSORICIDA: Model as EULIPOTYPHLA?
    # CHIROPTERA: Cut?
    # DERMOPTERA: Model as PRIMATES?
    # MICROBIOTHERIA: Model as DIDELPHIMORPHIA
    # MONOTREMATA: Model as DIPROTODONTIA
    # NOTORYCTEMORPHIA: Model as DASYUROMORPHIA
    # PAUCITUBERCULATA: Model as DIDELPHIMORPHIA
    # PHOLIDOTA: Model as CINGULATA
    # TUBULIDENTATA: Model as CINGULATA
   
    # For mammals - need to add some taxon manually
    if(taxon %in% 'mammalia') {
      order_tmp <- data.frame(order = c("AFROSORICIDA", #"CHIROPTERA", 
                                        "DERMOPTERA", "MICROBIOTHERIA", 
                                        "MONOTREMATA", "NOTORYCTEMORPHIA",
                                        "PAUCITUBERCULATA", "PHOLIDOTA", 
                                        "TUBULIDENTATA"),
                              order_intercept_adj_weighted = c(order_effects[order_effects$order == "EULIPOTYPHLA", "order_intercept_adj_weighted"], 
                                                               order_effects[order_effects$order == "PRIMATES", "order_intercept_adj_weighted"], 
                                                               order_effects[order_effects$order == "DIDELPHIMORPHIA", "order_intercept_adj_weighted"], 
                                                               order_effects[order_effects$order == "DIPROTODONTIA", "order_intercept_adj_weighted"], 
                                                               order_effects[order_effects$order == "DASYUROMORPHIA", "order_intercept_adj_weighted"], 
                                                               order_effects[order_effects$order == "DIDELPHIMORPHIA", "order_intercept_adj_weighted"], 
                                                               order_effects[order_effects$order == "CINGULATA", "order_intercept_adj_weighted"], 
                                                               order_effects[order_effects$order == "CINGULATA", "order_intercept_adj_weighted"]))
      
      order_effects <- bind_rows(order_effects, order_tmp) 
    } else { 
      # Do nothing
    }
    
    # Save the model-averaged coefficients -----
    # Creating directory
    dir.create(paste0(getwd(),'/Population Density Estimates'))
   
   cat('Saving files') 
    # Saving files
    # Exception by taxon
    if(taxon %in% c('aves','mammalia')) {
      write_csv(fixed_effects,
                path = paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_fixed_",taxon,".csv"))
      write_csv(random_effects_weighted[[1]],
                path = paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_binomial_",taxon,".csv"))
      write_csv(random_effects_weighted[[2]],
                path = paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_family_",taxon,".csv"))
      write_csv(order_effects,
                path = paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_order_",taxon,".csv"))
      
      # I'm also saving the total data to check things
      write_csv(tetra_check,
                path = paste0(getwd(),"/Population Density Estimates/TetraDatabaseForChecking_",taxon,".csv"))
    } else { # for amphibians
      write_csv(fixed_effects,
                path = paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_fixed_",taxon,".csv"))
      write_csv(random_effects_weighted[[1]],
                path = paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_family_",taxon,".csv"))
      # write_csv(random_effects_weighted[[2]],
      #           path = "/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Pop_Density_Data/RefittedSantiniModelCoefficients_random_family_Amphibians.csv")
      # write_csv(order_effects,
      #           path = "/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Pop_Density_Data/RefittedSantiniModelCoefficients_random_order_Birds.csv")
      
      # I'm also saving the total data to check things
      write_csv(tetra_check,
                path = paste0(getwd(),"/Population Density Estimates/TetraDatabaseForChecking_",taxon,".csv"))

    }
  } # End model fitting function

# Now running function
for(t in c('amphibia','mammalia','aves')) {
  pop_density_function(t)
}

# END




