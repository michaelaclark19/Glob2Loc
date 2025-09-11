#####
# Script defines species as common or rare under three metrics
# Based on methodoogy from Sykes et al 2019
# https://conbio.onlinelibrary.wiley.com/doi/pdf/10.1111/cobi.13419
# (1) Pop density
# (2) ecological breadth (how many IUCN habitats are suitable)
# (3) Range size
# Difference from Sykes is three-fold
# (1) Doing this by taxa (birds, mammals, and amphibians), rather than not separating by taxa
# (2) Across all ~20,000 bird + mammal + amphibian species reported in IUCN, rather than only using species reported in PREDICTs
# (3) using estimated population density data rather than reported data from TETRADensity in order to incorporate all species
#####

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
library(stringr)

# Setting working directory
setwd("/Users/maclark/Desktop/Multiple Stresses of Biodiversity")

# Species rarity data management ----
# Three parts
# 1) Number of habitat preferences
# 2) Pop density estimates
# 3) Habitat range size

# Number habitat preferences
amp.pref <- read.csv(paste0(getwd(),"/Other Data Inputs/Habitat_Preferences/AmphibianHabitats2018.csv"), stringsAsFactors = FALSE) %>% # Data
  mutate(Count = ifelse(Suitability %in% c('Suitable'),1,ifelse(Suitability %in% 'Marginal',.3,0))) %>% # Classifying habitat tolerance
  group_by(species, Species_ID) %>% # Grouping by species
  dplyr::summarise(Count = sum(Count)) %>% # Summing number of habitats
  mutate(taxon = 'Amphibia') %>% # Changing taxon
  mutate(Rare = ifelse(Count > median(.$Count),'No','Yes')) # Creating column for rarity

# Repeating for birds
bird.pref <- read.csv(paste0(getwd(),"/Other Data Inputs/Habitat_Preferences/BirdHabitats2017.csv"), stringsAsFactors = FALSE) %>%
  mutate(Hab_importance = paste(Suitability, Major.importance.)) %>% 
  mutate(Count = ifelse(Hab_importance %in% c('Suitable ', 'Suitable Yes'),1,ifelse(Hab_importance %in% c('Suitable No'),.5,ifelse(Hab_importance %in% c('Marginal ','Unknown ',' '),.3,0)))) %>%
  mutate(species = Scientific.name, Species_ID = SIS.ID) %>%
  group_by(species, Species_ID) %>%
  dplyr::summarise(Count = sum(Count)) %>%
  mutate(taxon = 'Aves') %>%
  mutate(Rare = ifelse(Count > median(.$Count),'No','Yes'))

# Repeating for mammals
mam.pref <- read.csv(paste0(getwd(),"/Other Data Inputs/Habitat_Preferences/MammalHabitats2017.csv"), stringsAsFactors = FALSE) %>%
  mutate(Hab_importance = paste(Suitability, Major.importance.)) %>% 
  mutate(Count = ifelse(Hab_importance %in% c('Suitable ', 'Suitable Yes'),1,ifelse(Hab_importance %in% c('Suitable No'),.5,ifelse(Hab_importance %in% c('Marginal ','Unknown ',' '),.3,0)))) %>%
  mutate(species = Species, Species_ID = assessment.id) %>%
  group_by(species, Species_ID) %>%
  dplyr::summarise(Count = sum(Count)) %>%
  mutate(taxon = 'Mammalia') %>%
  mutate(Rare = ifelse(Count > median(.$Count),'No','Yes'))

# Stacking
hab.pref <-
  rbind(amp.pref,
        bird.pref,
        mam.pref) 

# 2) Habitat extent, as estimate from ESH rasters
# Getting list of esh rasters or habitat range maps
# These aren't provided in the SI - for access email:
# Birds: Graeme Buchanan
# Mammals: Carlo Rondinini
# Amphibians: Francesco Ficetola

# Alternatively - use the range maps available from IUCN
esh.df <-
  rbind(
    data.frame(taxon = 'Aves',
               Raster = list.files("/Users/maclark/Desktop/ESH_Tifs_12Oct/Birds",
                                   pattern = "\\.tif$", full.name = TRUE)),
    data.frame(taxon = 'Amphibia',
               Raster = list.files("/Users/maclark/Desktop/ESH_Tifs_12Oct/Amphibians",
                                   pattern = "\\.tif$", full.name = TRUE)),
    data.frame(taxon = 'Mammalia',
               Raster = list.files("/Users/maclark/Desktop/ESH_Tifs_12Oct/Mammals",
                                   pattern = "\\.tif$", full.name = TRUE)))

# Managing data, mostly to get the correct species name
esh.df <- esh.df %>%
  mutate(binomial = gsub("_Hab.*","",Raster)) %>%
  mutate(binomial = gsub("_P.*","",binomial)) %>%
  mutate(binomial = gsub("\\.tif$","",binomial)) %>%
  mutate(binomial = gsub("/Users/maclark/Desktop/ESH_Tifs_12Oct/Mammals/","",binomial)) %>%
  mutate(binomial = gsub("/Users/maclark/Desktop/ESH_Tifs_12Oct/Amphibians/","",binomial)) %>%
  mutate(binomial = gsub("/Users/maclark/Desktop/ESH_Tifs_12Oct/Birds/","",binomial)) %>%
  mutate(esh_area = NA) %>%
  mutate(Raster = as.character(Raster))

# Function to get number of cells of viable habitat
# Doing it in parallel
# This takes ~40 mins
t1 = Sys.time()
r2 <- mclapply(esh.df$Raster, function(x) {
  cellStats(raster(x), stat = 'sum', na.rm = TRUE)  
}, mc.cores = 10)  
Sys.time() - t1

# Creating empty vector to fill cell values in
r3 <- rep(NA, length(r2))

# ANd converting vector to list
for(i in 1:length(r3)) {
  r3[i] <- r2[[i]][1]
}

# Updating esh area
esh.df$esh_area <- r3

# 3) Pop density estimates
# Unfortunately need to do this using refitted equations from Santini et al 2018...
# Importing climate rasters
# First climate maps under current climate
pwarmest_current = raster(paste0(getwd(),"/Global Climate Maps/Precip_in_WarmestCurrent.tif"))
pcv_current = raster(paste0(getwd(),"/Global Climate Maps/CoefVar_PrecipCurrent.tif"))
npp_raster = raster(paste0(getwd(),"/Global Climate Maps/NPP_Current.tif"))
# These are based off AOH maps, but can be recreated with other data inputs as needed
richness.mams <- raster(paste0(getwd(),"/Other Data Inputs/Species Richness Rasters/Mammals_Richness.tif"))
richness.bird <- raster(paste0(getwd(),"/Other Data Inputs/Species Richness Rasters/Birds_Richness.tif"))
richness.amp <- raster(paste0(getwd(),"/Other Data Inputs/Species Richness Rasters/Amphibians_Richness.tif"))

# Total richness
richness <- 
  richness.mams +
  richness.bird + 
  richness.amp

# Creating raster stack
tmp.stack <- stack(pwarmest_current,
                   pcv_current,
                   npp_raster,
                   richness)

# Coefficients for models used to estimate population density
# Now getting coefficients estimates from regression models
# These models are rebuilt off Santini et al 2018
# Individual models for birds, mammals, and amphibians
coef.order <- 
  rbind(read.csv(paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_order_Mammals.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Mammalia'),
        read.csv(paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_order_Birds.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Aves')) %>%
  mutate(order = toupper(order))

coef.family <- 
  rbind(read.csv(paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_family_Mammals.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Mammalia'),
        read.csv(paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_family_Birds.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Aves'),
        read.csv(paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_family_Amphibians.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Amphibia')) %>%
  mutate(family = toupper(family))

coef.binomial <- 
  rbind(read.csv(paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_binomial_Mammals.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Mammalia'),
        read.csv(paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_binomial_Birds.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Aves'))

coef.climate <- 
  rbind(read.csv(paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_fixed_Mammals.csv"), stringsAsFactors = FALSE) %>% dplyr::rename(.,richness_weighted_coef = richness_weighted_coef) %>% dplyr::select(.,intercept_weighted_coef,log10_body_mass_g_weighted_coef,log10_body_mass_g2_weighted_coef,log10_body_mass_g3_weighted_coef,npp_weighted_coef,npp2_weighted_coef,pcv_weighted_coef,pcv2_weighted_coef,pwarmest_weighted_coef,pwarmest2_weighted_coef,richness_weighted_coef = richness_weighted_coef, class) %>% mutate(taxon = 'Mammalia'),
        read.csv(paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_fixed_Birds.csv"), stringsAsFactors = FALSE) %>% dplyr::select(.,intercept_weighted_coef,log10_body_mass_g_weighted_coef,log10_body_mass_g2_weighted_coef,log10_body_mass_g3_weighted_coef,npp_weighted_coef,npp2_weighted_coef,pcv_weighted_coef,pcv2_weighted_coef,pwarmest_weighted_coef,pwarmest2_weighted_coef,richness_weighted_coef = richness_weighted_coef, class) %>% mutate(taxon = 'Aves'),
        read.csv(paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_fixed_Amphibians.csv"), stringsAsFactors = FALSE)  %>% dplyr::select(.,intercept_weighted_coef,log10_body_mass_g_weighted_coef,log10_body_mass_g2_weighted_coef,log10_body_mass_g3_weighted_coef,npp_weighted_coef,npp2_weighted_coef,pcv_weighted_coef,pcv2_weighted_coef,pwarmest_weighted_coef,pwarmest2_weighted_coef,richness_weighted_coef = richness_weighted_coef, class) %>% mutate(taxon = 'Amphibia'))

# Saving
write.csv(esh.df,
          paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/Estimated Rarity Values 20February2020.csv"))

esh.df = read.csv(paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/Estimated Rarity Values 20February2020.csv"),stringsAsFactors = FALSE)

# Now getting estimated body mass
body_mass_data <- read.csv(paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/Body Mass Estimates 20February2020 Updated Taxonomy.csv"), stringsAsFactors = FALSE)

# Getting updates for mismatched taxonomies
mismatched.birds <-
  rbind(read.csv(paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/Matched From Birdlife.csv"),stringsAsFactors = FALSE),
        read.csv(paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/Mismatched After Birdlife 20Feb2020.csv"),stringsAsFactors = FALSE) %>%
          dplyr::select(AOH_name = binomial, Synonym = Body.mass.species)) %>%
          mutate(Synonym = gsub("_"," ",Synonym)) %>%
          mutate(AOH_name = gsub("_"," ", AOH_name))

esh.df <- 
  left_join(esh.df,
            mismatched.birds %>% mutate(binomial = gsub(" ","_",AOH_name))) %>%
  mutate(Synonym = ifelse(Synonym == '',NA,Synonym)) %>%
  mutate(body.mass.species = ifelse(is.na(Synonym), binomial, Synonym)) %>%
  mutate(body.mass.species = gsub(" ","_", body.mass.species)) %>%
  mutate(body.mass.species = ifelse(body.mass.species == '',NA,body.mass.species))

# Getting order, family, etc data for each species and merging into the data frame
esh.df1 <-
  left_join(esh.df,
            dplyr::select(body_mass_data, order = Order, family = Family, genus = Genus, body.mass.species = binomial, est_mass_kg, Family_Mass_kg, Genus_Mass_kg)) %>%
  mutate(order = toupper(order), family = toupper(family)) 

esh.df1 <-
  esh.df1 %>%
  mutate(merge_genus = gsub("_.*","",body.mass.species)) %>%
  left_join(.,
            dplyr::select(body_mass_data, order1 = Order, family1 = Family, merge_genus = Genus, Genus_Mass_kg1 = Genus_Mass_kg, Family_Mass_kg1 = Family_Mass_kg) %>% unique(.)) %>%
  mutate(family = ifelse(is.na(family), family1, family)) %>%
  mutate(order = ifelse(is.na(order), order1, order)) %>%
  mutate(Genus_Mass_kg = ifelse(is.na(Genus_Mass_kg), Genus_Mass_kg1,Genus_Mass_kg)) %>%
  mutate(Family_Mass_kg = ifelse(is.na(Family_Mass_kg), Family_Mass_kg1, Family_Mass_kg)) %>%
  mutate(est_mass_kg = ifelse(is.na(est_mass_kg), Genus_Mass_kg, est_mass_kg)) %>%
  mutate(est_mass_kg = ifelse(is.na(est_mass_kg), Family_Mass_kg, est_mass_kg)) %>%
  dplyr::select(taxon, Raster, binomial, esh_area, AOH_name, Synonym, body.mass.species, order, family, genus,
                est_mass_kg, Family_Mass_kg, Genus_Mass_kg) %>%
  group_by(taxon, Raster, binomial, esh_area, AOH_name, Synonym, body.mass.species, order, family, genus) %>%
  dplyr::summarise(est_mass_kg = mean(est_mass_kg,na.rm=TRUE),
            Genus_Mass_kg = mean(Genus_Mass_kg, na.rm=TRUE),
            Family_Mass_kg = mean(Family_Mass_kg,na.rm=TRUE))


esh.df1 <-
  esh.df1%>%
  left_join(.,coef.order) %>%
  left_join(.,coef.family) %>%
  left_join(.,coef.binomial) %>%
  left_join(.,coef.climate) %>%
  mutate(log10_body_mass = log10(est_mass_kg*1000))

# Converting NAs to 0s
esh.df1$order_intercept_adj_weighted[is.na(esh.df1$order_intercept_adj_weighted)] <- 0
esh.df1$family_intercept_adj_weighted[is.na(esh.df1$family_intercept_adj_weighted)] <- 0
esh.df1$binomial_intercept_adj_weighted[is.na(esh.df1$binomial_intercept_adj_weighted)] <- 0

# Getting non-species specific coefficients to save space and time later
esh.df1 <-
  esh.df1 %>%
  mutate(body_mass_taxa_coef = 
           log10_body_mass * log10_body_mass_g_weighted_coef +
           log10_body_mass^2 * log10_body_mass_g2_weighted_coef +
           log10_body_mass^3 * log10_body_mass_g3_weighted_coef +
           order_intercept_adj_weighted + family_intercept_adj_weighted + binomial_intercept_adj_weighted +
           intercept_weighted_coef)

# Filtering species without estimated body masses or habitat ranges
# Can't predict densities for these
esh.df2 <-
  esh.df1 %>%
  filter(!is.na(est_mass_kg)) %>%
  filter(esh_area > 0)

# Creating function that only takes a raster as an input
# Writing so that it can be parallelized
density.fun <-
  function(x) {
    # Raster
    tmp.raster = raster(esh.df2$Raster[x])
    # Cropping to habitat extent
    tmp.stack.cropped = crop(tmp.stack, tmp.raster)
    # Creating data frame
    tmp.df <-
      data.frame(pwarm = getValues(tmp.stack.cropped[[1]]),
                 pcv = getValues(tmp.stack.cropped[[2]]),
                 npp = getValues(tmp.stack.cropped[[3]]),
                 richness = getValues(tmp.stack.cropped[[4]])) %>%
      mutate(check = rowSums(.[,1:4])) %>%
      filter(!is.na(check)) %>%
      filter(richness > 0, pwarm > 0, npp > 0, pcv > 0) 
    # Manging data and names of data table
    reg.df <- 
      as.data.frame(summary(tmp.df))
    reg.df <-
      reg.df[c(grep('Median',reg.df$Freq), grep('Mean', reg.df$Freq)),] %>%
      mutate(Var1 = gsub(":.*","",Freq)) %>%
      mutate(Freq = gsub(".*: ","",Freq)) %>%
      mutate(Freq = gsub(".*:","",Freq)) %>%
      mutate(Var2 = trimws(Var2, 'both')) %>%
      mutate(Freq = trimws(Freq, 'both')) %>%
      mutate(Freq = as.numeric(Freq))
    
    reg.df.mean <- reg.df[grep('Mean',reg.df$Var1),]
    reg.df.median <- reg.df[grep('Median',reg.df$Var1),]
    # Estimating population density
    est.dens = 
      data.frame(taxon = esh.df2$taxon[x],
                 Species = esh.df2$binomial[x],
                 Mean = 10 ^ (esh.df2$body_mass_taxa_coef[x] +
                                reg.df.mean$Freq[reg.df.mean$Var2 %in% 'npp'] * esh.df2$npp_weighted_coef[x] +
                                (reg.df.mean$Freq[reg.df.mean$Var2 %in% 'npp'])^2 * esh.df2$npp2_weighted_coef[x] +
                                reg.df.mean$Freq[reg.df.mean$Var2 %in% 'pcv'] * esh.df2$pcv_weighted_coef[x] +
                                (reg.df.mean$Freq[reg.df.mean$Var2 %in% 'pcv'])^2 * esh.df2$pcv2_weighted_coef[x] +
                                reg.df.mean$Freq[reg.df.mean$Var2 %in% 'pwarm'] * esh.df2$pwarmest_weighted_coef[x] +
                                (reg.df.mean$Freq[reg.df.mean$Var2 %in% 'pwarm'])^2 * esh.df2$pwarmest2_weighted_coef[x] +
                                reg.df.mean$Freq[reg.df.mean$Var2 %in% 'richness'] * esh.df2$richness_weighted_coef[x]),
                 Median = 10 ^ (esh.df2$body_mass_taxa_coef[x] +
                                  reg.df.median$Freq[reg.df.median$Var2 %in% 'npp'] * esh.df2$npp_weighted_coef[x] +
                                  (reg.df.median$Freq[reg.df.median$Var2 %in% 'npp'])^2 * esh.df2$npp2_weighted_coef[x] +
                                  reg.df.median$Freq[reg.df.median$Var2 %in% 'pcv'] * esh.df2$pcv_weighted_coef[x] +
                                  (reg.df.median$Freq[reg.df.median$Var2 %in% 'pcv'])^2 * esh.df2$pcv2_weighted_coef[x] +
                                  reg.df.median$Freq[reg.df.median$Var2 %in% 'pwarm'] * esh.df2$pwarmest_weighted_coef[x] +
                                  (reg.df.median$Freq[reg.df.median$Var2 %in% 'pwarm'])^2 * esh.df2$pwarmest2_weighted_coef[x] +
                                  reg.df.median$Freq[reg.df.median$Var2 %in% 'richness'] * esh.df2$richness_weighted_coef[x]))
    
    return(est.dens)
  }


# Trial to make sure the function works, not in parallel
# for(i in 1:5) {
#   trial = density.fun(x = i)
# }

# Making sure function works in parallel
# Alright this works
trial = mclapply(1:5,density.fun,mc.cores=5)

# And now running the whole thing
# This takes a while
trial = mclapply(1:nrow(esh.df2),density.fun,mc.cores = 10)

# Making new data frame in case I accidentally overwrite the previous one
trial2 = trial
esh.df3 <- 
  esh.df1 %>%
  mutate(median_density = NA,
         mean_density = NA)

# And converting list into a vector, and appending data into the data frame containing species-specific info
for(i in 1:length(trial2)) {
  esh.df3$mean_density[esh.df3$binomial == trial2[[i]]$Species] <-
   trial2[[i]]$Mean
  esh.df3$median_density[esh.df3$binomial == trial2[[i]]$Species] <-
    trial2[[i]]$Median
}

# 
# for(i in 1:length(trial2)) {
#   esh.df3$mean_density[esh.df3$binomial == trial2[[i]]$Species] <-
#     trial2[[i]]$Mean
#   esh.df3$median_density[esh.df3$binomial == trial2[[i]]$Species] <-
#     trial2[[i]]$Median
# }

# Creating directory
dir.create(paste0(getwd(),"/Ag Intensity Outputs/Response to Habitat Intensification"))

# Saving file
write.csv(esh.df3,
          paste0(getwd(),"/Ag Intensity Outputs/Response to Habitat Intensification/Estimated Values for Rarity 1May2021.csv"), row.names = FALSE)


esh.df3 <- read.csv(paste0(getwd(),"/Ag Intensity Outputs/Response to Habitat Intensification/Estimated Values for Rarity 1May2021.csv"), stringsAsFactors = FALSE)
# Getting rare/common by niche breadth, range size, and estimated pop density ----
# Species are rare if below median niche breadth, range size, and estimated pop density
esh.df3 <-
  left_join(esh.df3,
            dplyr::rename(hab.pref, binomial = species) %>% as.data.frame() %>% mutate(binomial = gsub(" ","_",binomial)))

# Everything is automatically rare, and updating to common later
esh.df3$Density_Common <- 'Common'
esh.df3$Niche_Breadth_Common <- 'Common'
esh.df3$Range_Size_Common <- 'Common'
# Estimated density
taxon.list <- c('Amphibia','Aves','Mammalia')

tetra_density <- 
  read.csv(paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/TetraDENSITY copy2.csv"),stringsAsFactors = FALSE) %>%
  mutate(binomial = paste(Genus,"_",Species,sep = ""))

# Add min and max habitat densities reported in tetra density
esh.df3$min_density = 0
esh.df3$max_density = 0
for(i in 1:nrow(esh.df3)) {
  tmp.min = min(tetra_density$Density[tetra_density$binomial %in% esh.df3$binomial[i]])
  tmp.max = max(tetra_density$Density[tetra_density$binomial %in% esh.df3$binomial[i]])
  
  if(!is.finite(tmp.min) | tmp.min == tmp.max) {
    tmp.min = min(tetra_density$Density[tetra_density$Genus %in% esh.df3$genus[i]])
    tmp.max = max(tetra_density$Density[tetra_density$Genus %in% esh.df3$genus[i]])
  }
  
  if(!is.finite(tmp.min) | tmp.min == tmp.max) {
    tmp.min = min(tetra_density$Density[tetra_density$Family %in% esh.df3$family[i]])
    tmp.max = max(tetra_density$Density[tetra_density$Family %in% esh.df3$family[i]])
  }

  if(!is.finite(tmp.min) | tmp.min == tmp.max) {
    tmp.min = min(tetra_density$Density[tetra_density$Order %in% str_to_title(esh.df3$order[i])])
    tmp.max = max(tetra_density$Density[tetra_density$Order %in% str_to_title(esh.df3$order[i])])
  }
  
  if(!is.finite(tmp.min) | tmp.min == tmp.max) {
    tmp.min = min(tetra_density$Density[tetra_density$Class %in% str_to_title(esh.df3$class[i])])
    tmp.max = max(tetra_density$Density[tetra_density$Class %in% str_to_title(esh.df3$class[i])])
  }
  
  esh.df3$min_density[i] <- tmp.min
  esh.df3$max_density[i] <- tmp.max
}

# Allowing for +- 20% uncertainty
esh.df3$median_density[esh.df3$median_density < esh.df3$min_density*.8 & !is.na(esh.df3$median_density)] <- 
  esh.df3$min_density[esh.df3$median_density < esh.df3$min_density*.8 & !is.na(esh.df3$median_density)]

esh.df3$median_density[esh.df3$median_density > esh.df3$max_density*1.2 & !is.na(esh.df3$median_density)] <- 
  esh.df3$max_density[esh.df3$median_density > esh.df3$max_density*1.2 & !is.na(esh.df3$median_density)]

esh.df3$mean_density[esh.df3$mean_density < esh.df3$min_density*.8 & !is.na(esh.df3$mean_density)] <- 
  esh.df3$min_density[esh.df3$mean_density < esh.df3$min_density*.8 & !is.na(esh.df3$mean_density)]

esh.df3$mean_density[esh.df3$mean_density > esh.df3$max_density*1.2 & !is.na(esh.df3$mean_density)] <- 
  esh.df3$max_density[esh.df3$mean_density > esh.df3$max_density*1.2 & !is.na(esh.df3$mean_density)]

# Looping through taxon to get rare and common species
for(i in taxon.list) {
  tmp.df <- esh.df3[esh.df3$taxon %in% i,]
  
  tmp.df$Density_Common[tmp.df$mean_density < median(tmp.df$mean_density,na.rm=TRUE)] <- "Rare"
  tmp.df$Range_Size_Common[tmp.df$esh_area < median(tmp.df$esh_area[tmp.df$esh_area>0], na.rm=TRUE)] <- 'Rare'
  tmp.df$Range_Size_Common[tmp.df$esh_area == 0] <- 'Common'
  tmp.df$Niche_Breadth_Common[tmp.df$Count < median(tmp.df$Count, na.rm = TRUE)] <- 'Rare'
  
  if(i == 'Amphibia') {
    out.frame = tmp.df
  } else {
    out.frame = rbind(out.frame,tmp.df)
  }
}

# And writing file
write.csv(out.frame,
          paste0(getwd(),"/Ag Intensity Outputs/Response to Habitat Intensification/Estimated Values for Rarity Ag Intensity Estimates 1May2021.csv"), row.names = FALSE)


# 
# write.csv(esh.df3[is.na(esh.df3$order),],
#           '/Users/maclark/Desktop/Mismatched_Species_Taxonomy.csv',
#           row.names=FALSE)
# 
# write.csv(hab.pref,
#           '/Users/maclark/Desktop/Hab_Prefs_20February2020.csv',
#           row.names=FALSE)
