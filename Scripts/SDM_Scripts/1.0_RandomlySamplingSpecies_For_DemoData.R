###
# Script for randomly selecting 25 species per taxon that are found only in meso america
# Selecting species only in meso america for two reasons
# (a) memory space required to run the SDMs (i.e. sdms may not run on laptops for species with larger geographic ranges)
# (b) for consistency with land use forecasting

###
# Setting seed for randomisation
set.seed(19)

###
# Libraries
library(readr)
library(plyr)
library(dplyr)


### 
# COO data
tmp_coo <-
  rbind(read_csv(paste0(getwd(),'/COO_Data/NewCountryOfOccurrenceData.csv')),
        read_csv(paste0(getwd(),'/COO_Data/coo_dat_reptiles.csv')))

###
# Which species are found outside of meso america
tmp_coo <-
  tmp_coo %>%
  # Tagging which species found outside meso america
  mutate(not_meso = !grepl('MEX|CRI|GTM|BLZ|HND|SLV|NIC|PAN',ISO3)) %>%
  # Summing number of countries species are found in, where the country is not in meso america
  dplyr::group_by(taxon, species) %>%
  dplyr::summarise(num_countries_not_meso = sum(not_meso)) %>%
  # Removing species that are found outside of mesoamerica
  filter(num_countries_not_meso %in% 0) %>%
  # ANd randomly sampling 25 species per taxon
  dplyr::ungroup() %>%
  dplyr::group_by(taxon) %>%
  sample_n(25)

###
# And moving files
# For mammals, reptiles, and amphibians
esh_tif_dir <- '/Volumes/Citadel/Multiple_Stresses_of_Biodiversity/Multiple_Stresses_of_Biodiversity/ESH_Tifs_12Oct/'

###
# Adding taxon to species name - needed to copy and paste files
tmp_coo_move <-
  tmp_coo %>%
  mutate(binomial_move = 
           case_when(taxon == 'Reptiles' ~ species,
                     .default = paste0(taxon,'/',species)))

###
# Getting files to move
files_from <- c()
for(ii in 1:nrow(tmp_coo_move)) {
  files_from <- 
    c(files_from, 
      list.files(paste0('/Volumes/Citadel/Multiple_Stresses_of_Biodiversity/Multiple_Stresses_of_Biodiversity/ESH_Tifs_12Oct/',tmp_coo_move$taxon[ii]),
                 full.names = TRUE) %>%
        .[grepl(tmp_coo_move$binomial_move[ii],.)])
}

# not reptiles
files_to <- gsub('/Volumes/Citadel/Multiple_Stresses_of_Biodiversity/Multiple_Stresses_of_Biodiversity/',
                 '/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Files_For_Katia/',
                 files_from)

# Moving
file.copy(files_from, files_to)

###
# Now for reptiles
files_from <- paste0('/Volumes/Castle/Desktop/Research/Glob2Loc/Reptile_Rasters/', tmp_coo_move$binomial_move[tmp_coo_move$taxon == 'Reptiles'] %>% gsub('Reptiles/','',.),'.tif')
files_to <- gsub('/Volumes/Castle/Desktop/Research/Glob2Loc/Reptile_Rasters/',
                 '/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Files_For_Katia/ESH_Tifs_12Oct/Reptiles/',
                 files_from)

# Moving
file.copy(files_from, files_to)
        