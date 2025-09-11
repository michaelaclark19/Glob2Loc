# Refitting a single model to the most important variables in the Santini dataset

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


# System options... 
options(tibble.width = Inf)
rm(list=ls())

# Setting working directory
setwd("/Users/maclark/Desktop/Multiple Stresses of Biodiversity")


# Only need to run the first two chunks the first time
# Load density data ----
tetra_data <- read_csv(paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/TetraDENSITY Copy2.csv")) %>%
  filter(Class %in% c("Mammalia", 'Aves', 'Amphibia')) %>%
  dplyr::select(.,Class:Species, Longitude, Latitude, Density) %>%
  filter(complete.cases(.))

names(tetra_data) <- tolower(names(tetra_data))

# Get coordinates for spatial data
coords <- tetra_data %>%
  dplyr::select(.,x = longitude,
                  y = latitude)

# Load explanatory variables and bind to density data ----
body_mass_data <- read.csv(paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/Body Mass Estimates 14November2019.csv"), stringsAsFactors = FALSE)

name_conversion <- read.csv(paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/SantiniNameConversion.csv"), stringsAsFactors = FALSE)

# Importing data for current climate ----
temp_mean <- raster("/Users/maclark/Desktop/Multiple Stresses of Biodiversity/Other Data Inputs/Pop Density Inputs/mean_temp_santini.tif")
precip_var <- raster("/Users/maclark/Desktop/Multiple Stresses of Biodiversity/Other Data Inputs/Pop Density Inputs/precip_var_santini.tif")
precip_warmest <- raster("/Users/maclark/Desktop/Multiple Stresses of Biodiversity/Other Data Inputs/Pop Density Inputs/precip_warmest_santini.tif")
temp_var <- raster("/Users/maclark/Desktop/Multiple Stresses of Biodiversity/Other Data Inputs/Pop Density Inputs/temp_var_santini.tif")
npp <- raster(paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/npp_santini.tif"))

# Richness rasters ----
# These are based off AOH maps, but can be recreated with other data inputs as needed
mam.richness <- raster(paste0(getwd(),"/Other Data Inputs/Species Richness Rasters/Mammals_Richness.tif"))
bird.richness <- raster(paste0(getwd(),"/Other Data Inputs/Species Richness Rasters/Birds_Richness.tif"))
amp.richness <- raster(paste0(getwd(),"/Other Data Inputs/Species Richness Rasters/Amphibians_Richness.tif"))

richness = mam.richness + bird.richness + amp.richness

# Resampling rasters to appropriate resolution and CRS ----
mam.richness <- round(projectRaster(mam.richness, npp, method = 'bilinear', na.rm = TRUE), digits = 0)
amp.richness <- round(projectRaster(amp.richness, npp, method = 'bilinear', na.rm = TRUE), digits = 0)
bird.richness <- round(projectRaster(bird.richness, npp, method = 'bilinear', na.rm = TRUE), digits = 0)
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
         richness = raster::extract(richness, coords)) %>%
  dplyr::select(.,class, order, family, genus, species, binomial,
         density, log10_density,
         combined_body_mass_g, log10_body_mass_g,
         longitude, latitude, npp, pcv, pwarmest, tvar, tmean,
         richness_mammals, richness_birds, richness_amps, richness) %>%
  filter(complete.cases(.)) %>%
  mutate(richness_mammals = ceiling(richness_mammals),
         richness_birds = ceiling(richness_birds),
         richness_amps = ceiling(richness_amps),
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

tetra_data$tot_richness <-
  tetra_data$richness

tetra_data <-
  tetra_data %>% filter(tot_richness >= 1)

tetra_data$richness_amps[tetra_data$richness_amps == 0] <- NA
tetra_data$richness_mammals[tetra_data$richness_mammals == 0] <- NA
tetra_data$richness_birds[tetra_data$richness_birds == 0] <- NA

write_csv(tetra_data,
          paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/SantiniModellingData SantiniCoefficients_SantiniClimateMaps.csv"))
# Define the maximum model ----
rm(list=ls())
tetra_data <- read_csv(paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/SantiniModellingData SantiniCoefficients_SantiniClimateMaps.csv"))

tetra_mam <- tetra_data %>% filter(class %in% 'MAMMALIA') %>% dplyr::select(.,-richness_birds) %>% dplyr::select(.,-richness_amps) %>%  filter(complete.cases(.)) %>% mutate(class = tolower(class), family = tolower(family)) %>% mutate(class = Hmisc::capitalize(class), family = Hmisc::capitalize(family))
m_max <- lmer(data = tetra_mam, 
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
# Extracting coefficients
dredge_table <- dredge(global.model = m_max, 
                       subset = m_subset,
                       beta = 'none',
                       evaluate = TRUE)


# m_max1 <- lmer(data = tetra_mam, 
#               formula = log10_density ~ 
#                 log10_body_mass_g + I(log10_body_mass_g^2) + I(log10_body_mass_g^3) +
#                 npp + I(npp^2) + 
#                 pcv + I(pcv^2) +
#                 # pwarmest + I(pwarmest^2) +
#                 richness_mammals  +
#                 (1|order/family/binomial), 
#               na.action = na.fail,
#               REML = FALSE)
# 
# # Fit models, ensuring that quadratics can only appear with linear effects
# m_subset1 <- expression(dc(log10_body_mass_g, `I(log10_body_mass_g^2)`, `I(log10_body_mass_g^3)`) &
#                          dc(`npp`, `I(npp^2)`) &
#                          dc(`pcv`, `I(pcv^2)`) &
#                          # dc(`pwarmest`, `I(pwarmest^2`) &
#                          dc(`richness_mammals`))
# # Extracting coefficients
# dredge_table1 <- dredge(global.model = m_max1, 
#                        subset = m_subset1,
#                        beta = 'none',
#                        evaluate = TRUE)

# # Do I want to rescale the variables?! - Not second time around ----
# mean_sd_data <- data.frame(var = c("log10_density", "log10_body_mass_g", "npp", "pcv", "richness_mammals"),
#                            mean = NA,
#                            sd = NA,
#                            stringsAsFactors = FALSE)
# 
# for(i in 1:nrow(mean_sd_data)){
#   mean_sd_data$mean[i] <- mean(tetra_mam[[mean_sd_data[i,"var"]]], na.rm = TRUE)
#   mean_sd_data$sd[i] <- sd(tetra_mam[[mean_sd_data[i,"var"]]], na.rm = TRUE)
# }
# 
# for(i in 1:nrow(mean_sd_data)){
#   n <- paste(mean_sd_data[i,"var"], "rescaled", sep = "_")
#   tetra_mam[[n]] <- (tetra_mam[[mean_sd_data[i,"var"]]] - mean_sd_data$mean[i]) / mean_sd_data$sd[i]
# }
# 
# # Refit
# m_max_rescaled <- lmer(data = tetra_mam,
#                        formula = log10_density_rescaled ~ log10_body_mass_g_rescaled + I(log10_body_mass_g_rescaled^2) + I(log10_body_mass_g_rescaled^3) +
#                 npp_rescaled + I(npp_rescaled^2) +
#                 pcv_rescaled + I(pcv_rescaled^2) +
#                 richness_mammals_rescaled +
#                 (1|order/family/binomial),
#               na.action = na.fail,
#               REML = FALSE)
# 
# m_subset_rescaled <- expression(dc(log10_body_mass_g_rescaled, `I(log10_body_mass_g_rescaled^2)`, `I(log10_body_mass_g_rescaled^3)`) &
#                          dc(`npp_rescaled`, `I(npp_rescaled^2)`) &
#                          dc(`pcv_rescaled`, `I(pcv_rescaled^2)`))
# 
# dredge_table_rescaled <- dredge(global.model = m_max_rescaled, subset = m_subset_rescaled)

# Which are the best models? -----
subset(dredge_table, cumsum(weight) < 0.999)
# subset(dredge_table_rescaled, cumsum(weight) < 0.999)
subset(dredge_table, weight > 0.01)
# cor(tetra_mam[, c("log10_density_rescaled", "npp_rescaled", "pcv_rescaled", "richness_mammals_rescaled")])
# Moderate correlation (0.66) between NPP and richness—may want to think about this? 

# rm(m_max_rescaled, dredge_table_rescaled)
# Same two models both times with identical deltas and weights.
# So use these models and the coefficients from the non-rescaled models. Hurrah!
# NB. can use the weights too as all other models are so poor in comparison. Hurrah again!

# Extract models and coefficients and weight them -----
averaged_models <- model.avg(dredge_table, subset = weight > .01)
averaged_models <- dredge_table %>% filter(weight > .01)


models <- get.models(dredge_table, subset = weight > 0.01)
map(models, r.squaredGLMM) 
# Pretty damn similar: ~0.3 for marginal; 0.77 for conditional
map(models, dispersion_glmer) # 0.62, OK
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
map(models, vif.lme) # I can't remember what this means...

weights <- data.frame(subset(dredge_table, weight > 0.01)) %>%
  dplyr::select(weight) %>%
  mutate(model_number = row.names(.))

# Fixed effects ------
fixed_effects <- model.avg(dredge_table, subset = weight > 0.001)
fixed_effects <- data.frame(fixed_effects$coefficients)
fixed_effects <- fixed_effects[row.names(fixed_effects) %in% 'full',]

# fixed_effects <- dredge_table %>% filter(weight > .01)
# fixed_effects <- data.frame(fixed_effects)


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


row.names(fixed_effects) <- NULL
fixed_effects <- fixed_effects %>%
  mutate(class = "MAMMALIA")

# Random effects -----
random_effects <- map(models, ranef)

ran_ef_fun <- function(l){
  for(i in 1:length(l)){
    name <- gsub(":.*$", "", names(l)[i])
    l[[i]][[name]]<- gsub(":.*$", "", row.names(l[[i]]))
    row.names(l[[i]]) <- NULL
    names(l[[i]])[1] <- paste(name, "intercept_adj", sep = "_")
  }
  return(l)
}

random_effects <- map(random_effects, ran_ef_fun)
# random_effects_weighted <- random_effects[[1]]
random_effects_weighted <- list()
for(i in 1:length(random_effects[[1]])){
  random_effects_weighted[[names(random_effects[[1]][i])]] <- rbindlist(map(random_effects, `[[`, names(random_effects[[1]][i])),
                                                                        idcol = "model_number")
}
random_effects_weighted <- map(random_effects_weighted, left_join, weights)

weight_fun <- function(df){
  df$w_tmp <- df[,2] * df$weight
  df_new <- df %>%
    group_by_(names(df)[3]) %>%
    dplyr::summarise(intercept_adj_weighted = sum(w_tmp)) %>%
    ungroup() 
  names(df_new)[2] <- paste(names(df_new)[1],
                            names(df_new)[2], 
                            sep = "_")
  return(df_new)
}

random_effects_weighted <- map(random_effects_weighted, weight_fun)

# Add coefficients to data to check predictions -----
tetra_mam1 <- tetra_mam %>%
  left_join(., fixed_effects %>% mutate(class = 'Mammalia')) %>%
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
         check = predicted_coefficients - predicted_model) %>%
  mutate(predicted_model_density = 10^predicted_model,
         predicted_density = 10^predicted_coefficients)

View(tetra_mam1 %>% dplyr::select(binomial, density, predicted_model_density, predicted_density, log10_density, predicted_model, predicted_coefficients) %>% filter(binomial %in% 'Acinonyx_jubatus'))

range(round(tetra_mam1$check, digits = 10), na.rm = TRUE)
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

order_effects <- random_effects_weighted[[3]] %>%
  as.data.frame()
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

# Save the model-averaged coefficients -----
# Creating directory
dir.create(paste0(getwd(),'/Population Density Estimates'))
write_csv(fixed_effects,
          path = paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_fixed_Mammals_SantiniClimateMaps.csv"))
write_csv(random_effects_weighted[[1]],
          path = paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_binomial_Mammals_SantiniClimateMaps.csv"))
write_csv(random_effects_weighted[[2]],
          path = paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_family_Mammals_SantiniClimateMaps.csv"))
write_csv(order_effects,
          path = paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_order_Mammals_SantiniClimateMaps.csv"))

# I'm also saving the total data to check things
write_csv(tetra_mam1,
          path = paste0(getwd(),"/Population Density Estimates/TetraDatabaseForChecking_Mammals_SantiniClimateMaps.csv"))

# Repeating for birds ----
# Define the maximum model ----
rm(list=ls())
tetra_data <- read_csv(paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/SantiniModellingData.csv"))

tetra_birds <- tetra_data %>% filter(class %in% 'AVES') %>% dplyr::select(.,-richness_mammals) %>% dplyr::select(.,-richness_amps) %>% dplyr::select(.,-tot_richness) %>% filter(complete.cases(.)) %>% mutate(class = tolower(class), family = tolower(family)) %>% mutate(class = Hmisc::capitalize(class), family = Hmisc::capitalize(family))
m_max <- lmer(data = tetra_birds, 
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
# Extracting coefficients
dredge_table <- dredge(global.model = m_max, 
                       subset = m_subset,
                       beta = 'none',
                       evaluate = TRUE)

dredge_table1 = dredge_table
# 
# m_max1 <- lmer(data = tetra_birds, 
#               formula = log10_density ~ 
#                 log10_body_mass_g + I(log10_body_mass_g^2) + I(log10_body_mass_g^3) +
#                 npp + I(npp^2) + 
#                 pcv + I(pcv^2) +
#                 # pwarmest + I(pwarmest^2) +
#                 richness_birds  +
#                 (1|order/family/binomial), 
#               na.action = na.fail,
#               REML = FALSE)
# 
# # Fit models, ensuring that quadratics can only appear with linear effects
# m_subset1 <- expression(dc(log10_body_mass_g, `I(log10_body_mass_g^2)`, `I(log10_body_mass_g^3)`) &
#                          dc(`npp`, `I(npp^2)`) &
#                          dc(`pcv`, `I(pcv^2)`) &
#                          # dc(`pwarmest`, `I(pwarmest^2`) &
#                          dc(`richness_birds`))
# # Extracting coefficients
# dredge_table1 <- dredge(global.model = m_max1, 
#                        subset = m_subset1,
#                        beta = 'none',
#                        evaluate = TRUE)

# Because for some reason this gets 

# # Do I want to rescale the variables?! - Not second time around ----
# mean_sd_data <- data.frame(var = c("log10_density", "log10_body_mass_g", "npp", "pcv", "richness_mammals"),
#                            mean = NA,
#                            sd = NA,
#                            stringsAsFactors = FALSE)
# 
# for(i in 1:nrow(mean_sd_data)){
#   mean_sd_data$mean[i] <- mean(tetra_mam[[mean_sd_data[i,"var"]]], na.rm = TRUE)
#   mean_sd_data$sd[i] <- sd(tetra_mam[[mean_sd_data[i,"var"]]], na.rm = TRUE)
# }
# 
# for(i in 1:nrow(mean_sd_data)){
#   n <- paste(mean_sd_data[i,"var"], "rescaled", sep = "_")
#   tetra_mam[[n]] <- (tetra_mam[[mean_sd_data[i,"var"]]] - mean_sd_data$mean[i]) / mean_sd_data$sd[i]
# }
# 
# # Refit
# m_max_rescaled <- lmer(data = tetra_mam,
#                        formula = log10_density_rescaled ~ log10_body_mass_g_rescaled + I(log10_body_mass_g_rescaled^2) + I(log10_body_mass_g_rescaled^3) +
#                 npp_rescaled + I(npp_rescaled^2) +
#                 pcv_rescaled + I(pcv_rescaled^2) +
#                 richness_mammals_rescaled +
#                 (1|order/family/binomial),
#               na.action = na.fail,
#               REML = FALSE)
# 
# m_subset_rescaled <- expression(dc(log10_body_mass_g_rescaled, `I(log10_body_mass_g_rescaled^2)`, `I(log10_body_mass_g_rescaled^3)`) &
#                          dc(`npp_rescaled`, `I(npp_rescaled^2)`) &
#                          dc(`pcv_rescaled`, `I(pcv_rescaled^2)`))
# 
# dredge_table_rescaled <- dredge(global.model = m_max_rescaled, subset = m_subset_rescaled)

# Which are the best models? -----
subset(dredge_table, cumsum(weight) < 0.999)
# subset(dredge_table_rescaled, cumsum(weight) < 0.999)
subset(dredge_table, weight > 0.01)
# cor(tetra_mam[, c("log10_density_rescaled", "npp_rescaled", "pcv_rescaled", "richness_mammals_rescaled")])
# Moderate correlation (0.66) between NPP and richness—may want to think about this? 

# rm(m_max_rescaled, dredge_table_rescaled)
# Same two models both times with identical deltas and weights.
# So use these models and the coefficients from the non-rescaled models. Hurrah!
# NB. can use the weights too as all other models are so poor in comparison. Hurrah again!

# Extract models and coefficients and weight them -----
averaged_models <- model.avg(dredge_table, subset = weight > .01)

models <- get.models(dredge_table, subset = weight > 0.01)
map(models, r.squaredGLMM) 
# Pretty damn similar: ~0.3 for marginal; 0.77 for conditional
map(models, dispersion_glmer) # 0.62, OK
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
map(models, vif.lme) # I can't remember what this means...

weights <- data.frame(subset(dredge_table, weight > 0.01)) %>%
  dplyr::select(weight) %>%
  mutate(model_number = row.names(.))

# Fixed effects ------
fixed_effects <- model.avg(dredge_table, subset = weight > 0.01)
fixed_effects <- data.frame(fixed_effects$coefficients)[1,]

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


row.names(fixed_effects) <- NULL
fixed_effects <- fixed_effects %>%
  mutate(class = "AVES")

# Random effects -----
random_effects <- map(models, ranef)

ran_ef_fun <- function(l){
  for(i in 1:length(l)){
    name <- gsub(":.*$", "", names(l)[i])
    l[[i]][[name]]<- gsub(":.*$", "", row.names(l[[i]]))
    row.names(l[[i]]) <- NULL
    names(l[[i]])[1] <- paste(name, "intercept_adj", sep = "_")
  }
  return(l)
}

random_effects <- map(random_effects, ran_ef_fun)
# random_effects_weighted <- random_effects[[1]]
random_effects_weighted <- list()
for(i in 1:length(random_effects[[1]])){
  random_effects_weighted[[names(random_effects[[1]][i])]] <- rbindlist(map(random_effects, `[[`, names(random_effects[[1]][i])),
                                                                        idcol = "model_number")
}
random_effects_weighted <- map(random_effects_weighted, left_join, weights)

weight_fun <- function(df){
  df$w_tmp <- df[,2] * df$weight
  df_new <- df %>%
    group_by_(names(df)[3]) %>%
    dplyr::summarise(intercept_adj_weighted = sum(w_tmp)) %>%
    ungroup() 
  names(df_new)[2] <- paste(names(df_new)[1],
                            names(df_new)[2], 
                            sep = "_")
  return(df_new)
}

random_effects_weighted <- map(random_effects_weighted, weight_fun)

# Add coefficients to data to check predictions -----
tetra_birds1 <- tetra_birds %>%
  left_join(., fixed_effects %>% mutate(class = 'Aves')) %>%
  left_join(., random_effects_weighted[[1]]) %>%
  left_join(., random_effects_weighted[[2]]) %>%
  left_join(., random_effects_weighted[[3]]) %>%
  mutate(richness = richness, 
         richness_weighted_coef = richness_weighted_coef) %>%
  mutate(predicted_model = predict(model.avg(dredge_table, subset = weight > 0.01, fit = TRUE)),
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
         coef_predict_dens = 10^predicted_coefficients,
         check = predicted_coefficients - predicted_model)
range(round(tetra_birds1$check, digits = 10), na.rm = TRUE)
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
# 
order_effects <- random_effects_weighted[[3]] %>%
  as.data.frame()
# order_tmp <- data.frame(order = c("AFROSORICIDA", #"CHIROPTERA", 
#                                   "DERMOPTERA", "MICROBIOTHERIA", 
#                                   "MONOTREMATA", "NOTORYCTEMORPHIA",
#                                   "PAUCITUBERCULATA", "PHOLIDOTA", 
#                                   "TUBULIDENTATA"),
#                         order_intercept_adj_weighted = c(order_effects[order_effects$order == "EULIPOTYPHLA", "order_intercept_adj_weighted"], 
#                                                          order_effects[order_effects$order == "PRIMATES", "order_intercept_adj_weighted"], 
#                                                          order_effects[order_effects$order == "DIDELPHIMORPHIA", "order_intercept_adj_weighted"], 
#                                                          order_effects[order_effects$order == "DIPROTODONTIA", "order_intercept_adj_weighted"], 
#                                                          order_effects[order_effects$order == "DASYUROMORPHIA", "order_intercept_adj_weighted"], 
#                                                          order_effects[order_effects$order == "DIDELPHIMORPHIA", "order_intercept_adj_weighted"], 
#                                                          order_effects[order_effects$order == "CINGULATA", "order_intercept_adj_weighted"], 
#                                                          order_effects[order_effects$order == "CINGULATA", "order_intercept_adj_weighted"]))

# order_effects <- bind_rows(order_effects, order_tmp)

# Save the model-averaged coefficients -----
write_csv(fixed_effects,
          path = paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_fixed_Birds_SantiniClimateMaps.csv"))
write_csv(random_effects_weighted[[1]],
          path = paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_binomial_Birds_SantiniClimateMaps.csv"))
write_csv(random_effects_weighted[[2]],
          path = paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_family_Birds_SantiniClimateMaps.csv"))
write_csv(order_effects,
          path = paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_order_Birds_SantiniClimateMaps.csv"))

# I'm also saving the total data to check things
write_csv(tetra_birds1,
          path = paste0(getwd(),"/Population Density Estimates/TetraDatabaseForChecking_Birds_SantiniClimateMaps.csv"))







# Repeating for mammals ----
# Define the maximum model ----
rm(list=ls())
tetra_data <- read_csv(paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/SantiniModellingData.csv"))

tetra_amps <- tetra_data %>% filter(class %in% 'AMPHIBIA') %>% dplyr::select(.,-richness_mammals) %>% dplyr::select(.,-richness_birds) %>% dplyr::select(.,-tot_richness) %>% filter(complete.cases(.)) %>% mutate(class = tolower(class), family = tolower(family)) %>% mutate(class = Hmisc::capitalize(class), family = Hmisc::capitalize(family))
m_max <- lmer(data = tetra_amps, 
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
# Extracting coefficients
dredge_table <- dredge(global.model = m_max, 
                       subset = m_subset,
                       beta = 'none',
                       evaluate = TRUE)

# Because for some reason this gets 

# # Do I want to rescale the variables?! - Not second time around ----
# mean_sd_data <- data.frame(var = c("log10_density", "log10_body_mass_g", "npp", "pcv", "richness_mammals"),
#                            mean = NA,
#                            sd = NA,
#                            stringsAsFactors = FALSE)
# 
# for(i in 1:nrow(mean_sd_data)){
#   mean_sd_data$mean[i] <- mean(tetra_mam[[mean_sd_data[i,"var"]]], na.rm = TRUE)
#   mean_sd_data$sd[i] <- sd(tetra_mam[[mean_sd_data[i,"var"]]], na.rm = TRUE)
# }
# 
# for(i in 1:nrow(mean_sd_data)){
#   n <- paste(mean_sd_data[i,"var"], "rescaled", sep = "_")
#   tetra_mam[[n]] <- (tetra_mam[[mean_sd_data[i,"var"]]] - mean_sd_data$mean[i]) / mean_sd_data$sd[i]
# }
# 
# # Refit
# m_max_rescaled <- lmer(data = tetra_mam,
#                        formula = log10_density_rescaled ~ log10_body_mass_g_rescaled + I(log10_body_mass_g_rescaled^2) + I(log10_body_mass_g_rescaled^3) +
#                 npp_rescaled + I(npp_rescaled^2) +
#                 pcv_rescaled + I(pcv_rescaled^2) +
#                 richness_mammals_rescaled +
#                 (1|order/family/binomial),
#               na.action = na.fail,
#               REML = FALSE)
# 
# m_subset_rescaled <- expression(dc(log10_body_mass_g_rescaled, `I(log10_body_mass_g_rescaled^2)`, `I(log10_body_mass_g_rescaled^3)`) &
#                          dc(`npp_rescaled`, `I(npp_rescaled^2)`) &
#                          dc(`pcv_rescaled`, `I(pcv_rescaled^2)`))
# 
# dredge_table_rescaled <- dredge(global.model = m_max_rescaled, subset = m_subset_rescaled)

# Which are the best models? -----
subset(dredge_table, cumsum(weight) < 0.999)
# subset(dredge_table_rescaled, cumsum(weight) < 0.999)
subset(dredge_table, weight > 0.01)
# cor(tetra_mam[, c("log10_density_rescaled", "npp_rescaled", "pcv_rescaled", "richness_mammals_rescaled")])
# Moderate correlation (0.66) between NPP and richness—may want to think about this? 

# rm(m_max_rescaled, dredge_table_rescaled)
# Same two models both times with identical deltas and weights.
# So use these models and the coefficients from the non-rescaled models. Hurrah!
# NB. can use the weights too as all other models are so poor in comparison. Hurrah again!

# Extract models and coefficients and weight them -----
averaged_models <- model.avg(dredge_table, subset = weight > .01)

models <- get.models(dredge_table, subset = weight > 0.01)
map(models, r.squaredGLMM) 
# Pretty damn similar: ~0.3 for marginal; 0.77 for conditional
map(models, dispersion_glmer) # 0.62, OK
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
map(models, vif.lme) # I can't remember what this means...

weights <- data.frame(subset(dredge_table, weight > 0.01)) %>%
  dplyr::select(weight) %>%
  mutate(model_number = row.names(.))

# Fixed effects ------
fixed_effects <- model.avg(dredge_table, subset = weight > 0.01)
fixed_effects <- data.frame(fixed_effects$coefficients)[1,]
# Does not have pwarmest^2, pcv^2, or npp^2 in the model
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


row.names(fixed_effects) <- NULL
fixed_effects <- fixed_effects %>%
  mutate(class = "AMPHIBIA")

# Random effects -----
random_effects <- map(models, ranef)

ran_ef_fun <- function(l){
  for(i in 1:length(l)){
    name <- gsub(":.*$", "", names(l)[i])
    l[[i]][[name]]<- gsub(":.*$", "", row.names(l[[i]]))
    row.names(l[[i]]) <- NULL
    names(l[[i]])[1] <- paste(name, "intercept_adj", sep = "_")
  }
  return(l)
}

random_effects <- map(random_effects, ran_ef_fun)
# random_effects_weighted <- random_effects[[1]]
random_effects_weighted <- list()
for(i in 1:length(random_effects[[1]])){
  random_effects_weighted[[names(random_effects[[1]][i])]] <- rbindlist(map(random_effects, `[[`, names(random_effects[[1]][i])),
                                                                        idcol = "model_number")
}
random_effects_weighted <- map(random_effects_weighted, left_join, weights)

weight_fun <- function(df){
  df$w_tmp <- df[,2] * df$weight
  df_new <- df %>%
    group_by_(names(df)[3]) %>%
    dplyr::summarise(intercept_adj_weighted = sum(w_tmp)) %>%
    ungroup() 
  names(df_new)[2] <- paste(names(df_new)[1],
                            names(df_new)[2], 
                            sep = "_")
  return(df_new)
}

random_effects_weighted <- map(random_effects_weighted, weight_fun)

# Add coefficients to data to check predictions -----
tetra_amps1 <- tetra_amps %>%
  left_join(., fixed_effects %>% mutate(class = 'Amphibia')) %>%
  left_join(., random_effects_weighted[[1]]) %>%
  # left_join(., random_effects_weighted[[2]]) %>%
  # left_join(., random_effects_weighted[[3]]) %>%
  mutate(richness = richness, 
         richness_weighted_coef = richness_weighted_coef) %>%
  mutate(predicted_model = predict(model.avg(dredge_table, subset = weight > 0.01, fit = TRUE)),
         predicted_coefficients = intercept_weighted_coef + 
           # binomial_intercept_adj_weighted +
           family_intercept_adj_weighted + 
           # order_intercept_adj_weighted + 
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
         coef_predict_dens = 10^predicted_coefficients,
         check = predicted_coefficients - predicted_model)
range(round(tetra_amps1$check, digits = 10), na.rm = TRUE)
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
# 
# order_effects <- random_effects_weighted[[3]] %>%
#   as.data.frame()
# order_tmp <- data.frame(order = c("AFROSORICIDA", #"CHIROPTERA", 
#                                   "DERMOPTERA", "MICROBIOTHERIA", 
#                                   "MONOTREMATA", "NOTORYCTEMORPHIA",
#                                   "PAUCITUBERCULATA", "PHOLIDOTA", 
#                                   "TUBULIDENTATA"),
#                         order_intercept_adj_weighted = c(order_effects[order_effects$order == "EULIPOTYPHLA", "order_intercept_adj_weighted"], 
#                                                          order_effects[order_effects$order == "PRIMATES", "order_intercept_adj_weighted"], 
#                                                          order_effects[order_effects$order == "DIDELPHIMORPHIA", "order_intercept_adj_weighted"], 
#                                                          order_effects[order_effects$order == "DIPROTODONTIA", "order_intercept_adj_weighted"], 
#                                                          order_effects[order_effects$order == "DASYUROMORPHIA", "order_intercept_adj_weighted"], 
#                                                          order_effects[order_effects$order == "DIDELPHIMORPHIA", "order_intercept_adj_weighted"], 
#                                                          order_effects[order_effects$order == "CINGULATA", "order_intercept_adj_weighted"], 
#                                                          order_effects[order_effects$order == "CINGULATA", "order_intercept_adj_weighted"]))

# order_effects <- bind_rows(order_effects, order_tmp)

# Save the model-averaged coefficients -----
write_csv(fixed_effects,
          path = paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_fixed_Amphibians_SantiniClimateMaps.csv"))
write_csv(random_effects_weighted[[1]],
          path = paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_family_Amphibians_SantiniClimateMaps.csv"))
# write_csv(random_effects_weighted[[2]],
#           path = "/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Pop_Density_Data/RefittedSantiniModelCoefficients_random_family_Amphibians.csv")
# write_csv(order_effects,
#           path = "/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Pop_Density_Data/RefittedSantiniModelCoefficients_random_order_Birds.csv")

# I'm also saving the total data to check things
write_csv(tetra_amps1,
          path = paste0(getwd(),"/Population Density Estimates/TetraDatabaseForChecking_Amphibians_SantiniClimateMaps.csv"))
