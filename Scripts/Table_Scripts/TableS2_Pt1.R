#!/usr/bin/env Rscript



# importing libraries
library(terra)
library(plyr)
library(dplyr)
library(stringr)
library(parallel)
library(future)
library(future.apply)

# Setting working directory
setwd("/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity")

# Creating directory to save csv outputs
dir.create(paste0(getwd(),'/Analyses/SDM_Analyses/Species_Results'),recursive=TRUE)

# List of taxa
taxa.list <-
  list.files(paste0(getwd(),'/ESH_RCPs/','SSP2-4.5')) %>%
  .[!grepl('.tif',.)] %>%
  .[!grepl('_',.)]

# List of years to calculate
years.list <-
  c(2020,2050)

# List of threshold types
thresh.list <-
  list.files(paste0(getwd(),'/ESH_RCPs/','SSP2-4.5','/','Amphibians'), pattern = 'migrate') %>%
  gsub('_migrate.*','',.) %>%
  gsub('.*[0-9]{4,4}_','',.) %>%
  unique()

# List of species with sdms
sdm.species <-
  do.call(c,
          lapply(taxa.list,
                 function(i) {
                   list.files(paste0(getwd(),'/ESH_RCPs/SSP5-8.5/',i),
                              pattern = 'csv') %>%
                     .[grepl('migrate',.)] %>%
                     paste0(i,'/',.) %>%
                     gsub('_migrate','',.)
                 }))



# Getting list of species with AOH maps (e.g. all species)
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
  left_join(.,
            data.frame(Species = sdm.species %>% gsub('.csv','',.), Have_SDM = 1)) %>%
  unique() %>%
  mutate(tetra_taxon = ifelse(taxon %in% 'Birds','Aves', # Adding taxon identifier for tetra density species
                              ifelse(taxon %in% 'Mammals','Mammalia',
                                     ifelse(taxon %in% 'Amphibians','Amphibia','Reptilia')))) %>%
  filter(!is.na(Have_SDM)) %>%
  mutate(Species = str_remove(Species,'Amphibians/|Mammals/|Birds/|Reptiles/'))

# realms and biomes
realm_biome_map <-
  rast(paste0(getwd(),'/Ecoregions_Feb2023/biomes_realms_raster.tif'))

# list of all sdm files
sdm_map_files_all <-
  do.call(c,
          lapply(taxa.list,
                 FUN = function(ii) {
                   list.files(paste0(getwd(),'/ESH_RCPs/','SSP2-4.5','/',ii),
                              full.names = TRUE)
                 }))

###
# Creating function that loops across species
sdm_area_function <-
  function(ss) {

    cat('Starting ',ss,' : ',species.frame$Species[ss],'\n')

    # getting aoh map for species
    # doing this to find the realm/biome combinations species is currently in
    aoh_map_file <-
      list.files(paste0(getwd(),'/ESH_Tifs_12Oct/',species.frame$taxon[ss]),
                 pattern = gsub(paste0(species.frame$taxon[ss],'/'),'',species.frame$Species[ss]),
                 full.names = TRUE)

    # Getting sdm map for species
    sdm_map_files <-
      # list.files(paste0(getwd() %>% gsub('ouce-glob2loc','pubh-glob2loc',.),'/ESH_RCPs/','SSP2-4.5','/',species.frame$taxon[ss]),
      #            # pattern = gsub(paste0(species.frame$taxon[k],'/'),'',species.frame$Species[k]),
      #            # pattern = species.frame$Species[k],
      #            full.names = TRUE) %>%
      #fs::dir_ls(paste0(getwd() %>% gsub('ouce-glob2loc','pubh-glob2loc',.),'/ESH_RCPs/','SSP2-4.5','/',species.frame$taxon[ss])) %>%
      #.[grepl(species.frame$Species[ss],.)] %>%
      # .[grepl(paste0(ssp,'_',y,'_'),.)] %>%
      #.[grepl('SSP2',.)] %>%
      # .[grepl(ss,.)] %>%
      # .[!grepl('migrate',.)] %>%
      #.[grepl('2020|2050',.)]
	    sdm_map_files_all %>%
	    .[grepl('2020|2050',.)] %>%
	    .[grepl(species.frame$taxon[ss],.)] %>%
	    .[grepl(species.frame$Species[ss],.)]

    # importing
    rast_stack <-
      rast(sdm_map_files)

    # extending
    aoh_map <-
      extend(rast(aoh_map_file),
             rast_stack)

    # getting realm biome species is currently found in
    tmp_realm_biome <-
      crop(realm_biome_map,
           aoh_map)

    # and getting values
    tmp_eco_aoh_values <-
      values(tmp_realm_biome)[values(aoh_map) %in% 1] %>%
      unique() %>%
      .[. >1] %>% # oceans
      .[.!=48] %>% # african greatlakes
      .[!is.na(.)]

    # removing other realms and biomes
    tmp_realm_biome[!(tmp_realm_biome %in% tmp_eco_aoh_values)] <- NA

    # clipping these from the raster stack
    rast_stack[is.na(tmp_realm_biome)] <- NA

    # And getting total number of cells available for species
    num_cells <-
      do.call(rbind,
              lapply(rast_stack,
             FUN = function(ii) {
               return(data.frame(scenario = names(ii),
                                 num_cells = sum(values(ii) %in% 1)))
             })
      ) %>%
      mutate(year = str_extract(scenario,'[0-9]{4,4}'),
             threshold = str_extract(scenario,'prevalence|specsens'),
             climate_migration_limits =
               case_when(grepl('migrate',scenario) ~ 'Yes',
                         .default = 'No')) %>%
      mutate(species = gsub('_SSP2.*','',scenario),
             taxon = species.frame$taxon[ss]) %>%
      dplyr::select(taxon, species, year, threshold, climate_migration_limits, num_cells)

    write.csv(num_cells,
              paste0(getwd(),
                     '/Analyses/SDM_Analyses/Species_Results/',
                     species.frame$taxon[ss],'_',
                     species.frame$Species[ss],'_SDM_Estimated_Habitat_Area.csv'),
              row.names = FALSE)
    # return(num_cells)
  }

###
# Getting files saved
species_have <-
  list.files(paste0(getwd(),
                    '/Analyses/SDM_Analyses/Species_Results/')) %>%
  gsub('_SDM.*','',.)

species_have <-
  paste0(str_extract(species_have,'Amphibians_|Birds_|Mammals_|Reptiles_') %>%
           gsub('_','/',.),
         str_remove(species_have,'Amphibians_|Birds_|Mammals_|Reptiles_'))


# removing secies we have
species.frame <-
  species.frame %>%
  filter(!(Raster %in% species_have)) %>%
  sample_n(nrow(.))


t1 = Sys.time()
mclapply(1:nrow(species.frame),sdm_area_function,mc.cores=10)
Sys.time() - t1
