#!/usr/bin/env Rscript


###
# libraries

library(raster)
library(plyr)
library(dplyr)
library(parallel)
library(stringr)
library(scales)

###
# setting wd
setwd('/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity')

# importing functions
# for pop density estimates
source(paste0(getwd(),"/Scripts/Calculating_Biodiversity_Scripts/0.0_ESH_Loss_Functions_17Aug2023.R"))

### (1) habitat area estimates
# function to get files
import_function <-
  function(taxon) {
    return(list.files(paste0(getwd(),'/ESH_Tifs_12Oct/',taxon),full.names = TRUE))
  }

# getting list of file paths
aoh_maps <-
  do.call(c,
          mclapply(c('Amphibians','Reptiles','Mammals','Birds'),
                   import_function,
                   mc.cores=4))

### (3) pop density estimates
#
# importing climate maps
pwarmest_raster = raster(paste0(getwd(),'/CMIP6_Climate_Data/Historic/1995-2014/Managed_Rasters/Precip_in_warmest_quarter_Mollweide.tif'))
pvar_raster = raster(paste0(getwd(),'/CMIP6_Climate_Data/Historic/1995-2014/Managed_Rasters/Variance_of_precipitation_Mollweide.tif'))
npp_raster = raster(paste0(getwd(),"/Global Mollweide Maps/NPP_Current.tif"))
temp_var = raster(paste0(getwd(),'/CMIP6_Climate_Data/Historic/1995-2014/Managed_Rasters/Variance_of_temperature_Mollweide.tif'))

# These are based off AOH maps, but can be recreated with other data inputs as needed
richness.mams <- raster(paste0(getwd(),"/Other Data Inputs/Species Richness Rasters/Mammals_Richness.tif"))
richness.bird <- raster(paste0(getwd(),"/Other Data Inputs/Species Richness Rasters/Birds_Richness.tif"))
richness.amp <- raster(paste0(getwd(),"/Other Data Inputs/Species Richness Rasters/Amphibians_Richness.tif"))
richness.rep <- raster(paste0(getwd(),"/Other Data Inputs/Species Richness Rasters/Reptile_Richness.tif"))

# Getting total richness
richness <- 
  richness.mams +
  richness.bird +
  richness.amp +
  richness.rep

# importing density data - used to set upper bound for population density predictions
tetra_density <- 
  read.csv(paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/TetraDENSITY_v.1.csv")) %>%
  mutate(binomial = paste(.$Genus,"_",.$Species, sep = "")) %>%
  mutate(density_corrected =
           ifelse(Density_unit %in% 'ind/ha', Density * 100,
                  ifelse(Density_unit %in% 'males/ha', Density * 2 * 100,
                         ifelse(Density_unit %in% 'pairs/km2', Density * 2, Density)))) %>%
  mutate(Density = density_corrected)

# list of species with AOH
# Getting list of species with AOH maps (e.g. all species)
taxa.list <- list.files(paste0(getwd() %>% gsub('ouce-glob2loc','pubh-glob2loc',.),'/ESH_RCPs/SSP2-4.5')) %>% .[!grepl('.tif',.)] %>% .[!grepl('Climate',.)]

aoh.species <-
  do.call(c,
          lapply(taxa.list, 
                 function(i) {
                   list.files(paste0(getwd(),'/ESH_Tifs_12Oct/',i)) %>%
                     gsub('_HabEl.*','',.) %>% # Removing nomenclature from the AOH habitats
                     gsub('_P[0-9].*','',.) %>% # Removing nomenclature from the AOH habitats
                     gsub('.tif','',.) %>% # Removing file type
                     paste0(i,'/',.) %>% # Adding taxa back on
                     unique() # Unique observations only
                 })) %>%
  .[!grepl('.aux',.)]

# Data set up for output data frame ----
species.frame <- 
  data.frame(Species = aoh.species,
             Raster = aoh.species) %>%
  mutate(taxon = gsub("/.*",'',Species)) %>% # Identifying taxon
  unique() %>%
  mutate(tetra_taxon = ifelse(taxon %in% 'Birds','Aves', # Adding taxon identifier for tetra density species
                              ifelse(taxon %in% 'Mammals','Mammalia',
                                     ifelse(taxon %in% 'Amphibians','Amphibia','Reptilia'))))

# getting coefficients for density estimates - these are the random effects for family/order/binomial
species.frame <-
  pop.dens.coef.function(species.frame)

# Getting climate coefficients - these are fixed effects for e.g. body mass, climate variables, etc
coef.climate <-
  rbind(read.csv(paste0(getwd() %>% gsub('pubh-glob2loc','ouce-glob2loc',.),"/Population Density Estimates/RefittedSantiniModelCoefficients_fixed_Mammals.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Mammals'),
        read.csv(paste0(getwd() %>% gsub('pubh-glob2loc','ouce-glob2loc',.),"/Population Density Estimates/RefittedSantiniModelCoefficients_fixed_Amphibians.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Amphibians'),
        read.csv(paste0(getwd() %>% gsub('pubh-glob2loc','ouce-glob2loc',.),"/Population Density Estimates/RefittedSantiniModelCoefficients_fixed_Birds.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Birds'),
        read.csv(paste0(getwd() %>% gsub('pubh-glob2loc','ouce-glob2loc',.),"/Population Density Estimates/RefittedSantiniModelCoefficients_fixed_Reptiles.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Reptiles'))

# function to make this work
density_function <- 
  function(k) {
    
    cat(k,'\n')
    
    # importing sdm map
    sdm_map <- raster(aoh_maps[grepl(species.frame$Raster[k],aoh_maps)])
    
    # getting habitat density range
    hab_min = 0
    hab_max = 0
    
    # Min will always be 0 individuals
    # Identifying max
    # This comes from the tetra density database
    
    # Checking if there is data from the species
    hab_max = max(tetra_density$Density[tetra_density$binomial == # is there a match
                                          species.frame$Species[k] %>% gsub('.*/','',.)], # removing taxon from species frame
                  na.rm = TRUE) * 2 # and adding flexibility to upper limit, assuming sampling does not fully capture actual range of population densities
    
    # If no data from the species, then going to the genus
    # updating for genus if not binomial
    if(!is.finite(hab_max) | hab_min == hab_max) {
      hab_max = max(tetra_density$Density[tetra_density$Genus %in% species.frame$genus[k]], na.rm = TRUE) * 2
    }
    
    # If no data from the genus, then going to the family
    if(!is.finite(hab_max) | hab_min == hab_max) {
      hab_max = max(tetra_density$Density[tetra_density$Family %in% str_to_title(species.frame$family[k])], na.rm = TRUE) * 2
    } 
    
    # And if no data from the family, then using the order
    if(!is.finite(hab_max) | hab_min == hab_max) {
      hab_max = max(tetra_density$Density[tetra_density$Order %in% str_to_title(species.frame$order[k])], na.rm = TRUE) * 2
    } 
    
    # And repeating for class if not order
    if(!is.finite(hab_max) | hab_min == hab_max) {
      # hab_min = min(tetra_density$Density[tetra_density$Family %in% str_to_title(species.frame$family[k])], na.rm = TRUE)
      hab_max = max(tetra_density$Density[tetra_density$Class %in% str_to_title(species.frame$tetra_taxon[k])], na.rm = TRUE) * 2
      # density.frame = density.frame %>% mutate(Species = species.frame$Species[k], Density_Min_Max = 'Family') %>%
      #   mutate(Density = squish(Density, range = c(hab_min, hab_max)))
    } 
    
    
    # running density estimates
    density_raster <- 
      pop_density_function(coefs = coef.climate %>% filter(taxon %in% species.frame$taxon[k]),
                           coef_order = species.frame$order_intercept_adj_weighted[k],
                           coef_family = species.frame$family_intercept_adj_weighted[k],
                           coef_binomial = species.frame$binomial_intercept_adj_weighted[k],
                           body_mass_grams = species.frame$est_mass_kg[k]*1000,
                           npp_input = crop(npp_raster, sdm_map),
                           pcv_input = crop(pvar_raster, sdm_map),
                           pwarmest_input = crop(pwarmest_raster, sdm_map),
                           richness_input = crop(richness, sdm_map),
                           hab_min = hab_min,
                           hab_max = hab_max,
                           # hab_max = 5000,
                           nrows = sdm_map@nrows)
    
    # getting mean value within species' habitat
    density.frame <- 
      density.function.new(esh.list = stack(sdm_map),
                           density.list = stack(density_raster)) %>%
      mutate(species = species.frame$Species[k])
    
    # writing density frame
    write.csv(density.frame,
              paste0(getwd(),
                     '/tmp_density_estimates/',
                     gsub('/.*','',density.frame$species),
                     '_',
                     gsub('.*/','',density.frame$species),
                     '.csv'),
              row.names = FALSE)
    
    return(density.frame)
  } # end density function



# now making wrapper function
try_catch_function <-
  function(k) {
    tryCatch(density_function(k),
             error = function(e) { # catching errors
               tmp_df <- data.frame(species = species.frame$Raster[k])
               write.csv(tmp_df,
                         paste0(getwd(),
                                '/tmp_density_estimates/ERROR_OH_NO_',
                                gsub('/.*','',species.frame$species),
                                '_',
                                gsub('.*/','',species.frame$species),
                                '.csv'),
                         row.names = FALSE)
             }) # end error catch
  } # end trycatch

# now running in parallel
mclapply(nrow(species.frame):1, try_catch_function, mc.cores=20)
