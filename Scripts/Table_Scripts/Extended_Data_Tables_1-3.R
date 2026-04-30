#####
###
# Building GLM model to regress outcomes against realm / biome / etc

###
# Libraries
library(fmsb)
library(plyr)
library(dplyr)
library(terra)
library(sf)
library(RColorBrewer)
library(lme4)

###
# Setting working directory
setwd("/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity")

#####
###
# Background data management

###
# Getting list of realms and biomes
tmp_sf <- read_sf('/Users/michael/Downloads/Terrestrial_Ecoregions/Terrestrial_Ecoregions.shp')

# Converting into a data frame
realms_biomes <-
  data.frame(ecoregion_id = tmp_sf$ECO_ID_U,
             ecoregion_name = tmp_sf$ECO_NAME,
             realm = tmp_sf$WWF_REALM2,
             biome = tmp_sf$WWF_MHTNAM)

###
# Getting species-level outcomes
file_list <-
  list.files(paste0(getwd(),'/Aggregated_CSV_Files/'),
             pattern = 'Prop',
             full.names = TRUE) %>%
  .[grepl('prevalence',.)] %>%
  .[grepl('nopatch',.)] %>%
  .[!grepl('Climate_Suitable',.)]

# importing files
out_df <-
  do.call(rbind, lapply(file_list, read_csv)) %>%
  mutate(taxon = gsub('/.*','',binomial)) %>%
  mutate(species = gsub('.*/','',binomial))

###
# which species have estimated area in 2020
species_list <-
  out_df %>%
  filter(year %in% 2020) %>%
  dplyr::select(taxon, species, tot_pop, tot_area) %>%
  unique() %>%
  mutate(tot_pop = as.numeric(tot_pop),
         tot_area = as.numeric(tot_area)) %>%
  # filter(!is.na(tot_pop) &
  #        !is.na(tot_area)) %>%
  filter(!is.na(tot_area)) %>%
  filter(tot_area > 0)

###
# Realm / biome combinations with >10% species habitat
realm_biome_species <- 
  read_csv('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/species_by_ecoregion_esh_maps.csv') %>%
  mutate(species = gsub('_[A-Z].*|_[0-9].*','',species)) %>%
  mutate(binomial = paste0(taxa,'/',species)) %>%
  dplyr::select(-c(species,taxa)) %>%
  left_join(.,realms_biomes) %>%
  # dplyr::group_by(binomial, realm, biome) %>%
  # dplyr::summarise(prop_cells = sum(prop_cells)) %>%
  # filter(!is.na(realm),
  #        !is.na(biome)) %>%
  filter(prop_cells > .25) %>%
  dplyr::select(ecoregion_id, ecoregion_name, binomial) %>% 
  distinct()

###
# Outcomes in 2050
out_df_2050_all_species <-
  out_df %>%
  filter(paste0(taxon,species) %in% paste0(species_list$taxon,species_list$species)) %>%
  filter(year %in% 2050) %>%
  left_join(., 
            realm_biome_species %>% distinct()) %>%# Merging in ecoregions species are present in
  #           read_csv('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/species_by_ecoregion_esh_maps.csv') %>%
  #             mutate(species = gsub('_[A-Z].*|_[0-9].*','',species)) %>%
  #             mutate(binomial = paste0(taxa,'/',species)) %>%
  #             dplyr::select(-c(species,taxa))) %>%
  # # filter(grepl('.int.urb',scenario)) %>% # Getting only total outcomes
  # left_join(., # Merging in biomes by ecoregion
  #           realms_biomes) %>%
  # filter(ecoregion_id > 0) %>%
  # filter(!is.na(ecoregion_id),
  #        !is.na(realm),
  #        !is.na(biome)) %>%
  # dplyr::select(-ecoregion_id) %>%
  # dplyr::select(-c(prop_cells_in_ecoregion,prop_cells,esh_map)) %>%
  distinct() %>%
  dplyr::group_by(scenario, ecoregion_id, ecoregion_name) %>%
  dplyr::summarise(median_pop = round(median(tot_pop,na.rm=TRUE),digits = 2),
                   pop_25 = round(quantile(tot_pop,.25,na.rm=TRUE),digits = 2),
                   pop_75 = round(quantile(tot_pop,.75,na.rm=TRUE),digits = 2),
                   sd_pop = round(sd(tot_pop,na.rm=TRUE),digits = 2),
                   median_area = round(median(tot_area,na.rm=TRUE),digits = 2),
                   area_25 = round(quantile(tot_area,.25,na.rm=TRUE),digits=2),
                   area_75 = round(quantile(tot_area,.75,na.rm=TRUE),digits=2),
                   sd_area = round(sd(tot_area,na.rm=TRUE),digits=2),
                   num_species = n()) %>%
  filter(!is.na(ecoregion_name)) %>%
  mutate(change_in_population_abundance = paste0(median_pop,' (',pop_25,' - ',pop_75,'); sd = ',sd_pop),
         change_in_habitat_area = paste0(median_area,' (',area_25,' - ',area_75,'); sd = ',sd_area)) %>%
  mutate(scenario = 
           case_when(grepl('current',scenario) ~ 'Climate Change',
                     grepl('int.urb',scenario) ~ 'All Stressors',
                     grepl('int',scenario) ~ 'Agricultural Intensification',
                     grepl('exp',scenario) ~ 'Agricultural Expansion',
                     grepl('urb',scenario) ~ 'Urban Expansion')) %>%
  mutate(taxon = 'All Species') %>%
  dplyr::select(scenario, ecoregion_id, ecoregion_name, taxon, change_in_population_abundance, change_in_habitat_area, number_of_species = num_species)

out_df_2050_by_taxon <-
  out_df %>%
  filter(paste0(taxon,species) %in% paste0(species_list$taxon,species_list$species)) %>%
  filter(year %in% 2050) %>%
  left_join(., 
            realm_biome_species %>% distinct()) %>%# Merging in ecoregions species are present in
  #           read_csv('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/species_by_ecoregion_esh_maps.csv') %>%
  #             mutate(species = gsub('_[A-Z].*|_[0-9].*','',species)) %>%
  #             mutate(binomial = paste0(taxa,'/',species)) %>%
  #             dplyr::select(-c(species,taxa))) %>%
  # # filter(grepl('.int.urb',scenario)) %>% # Getting only total outcomes
  # left_join(., # Merging in biomes by ecoregion
  #           realms_biomes) %>%
  # filter(ecoregion_id > 0) %>%
  # filter(!is.na(ecoregion_id),
  #        !is.na(realm),
  #        !is.na(biome)) %>%
  # dplyr::select(-ecoregion_id) %>%
  # dplyr::select(-c(prop_cells_in_ecoregion,prop_cells,esh_map)) %>%
  distinct() %>%
  dplyr::group_by(scenario, taxon, ecoregion_id, ecoregion_name) %>%
  dplyr::summarise(median_pop = round(median(tot_pop,na.rm=TRUE),digits = 2),
                   pop_25 = round(quantile(tot_pop,.25,na.rm=TRUE),digits = 2),
                   pop_75 = round(quantile(tot_pop,.75,na.rm=TRUE),digits = 2),
                   sd_pop = round(sd(tot_pop,na.rm=TRUE),digits = 2),
                   median_area = round(median(tot_area,na.rm=TRUE),digits = 2),
                   area_25 = round(quantile(tot_area,.25,na.rm=TRUE),digits=2),
                   area_75 = round(quantile(tot_area,.75,na.rm=TRUE),digits=2),
                   sd_area = round(sd(tot_area,na.rm=TRUE),digits=2),
                   num_species = n()) %>%
  filter(!is.na(ecoregion_name)) %>%
  mutate(change_in_population_abundance = paste0(median_pop,' (',pop_25,' - ',pop_75,'); sd = ',sd_pop),
         change_in_habitat_area = paste0(median_area,' (',area_25,' - ',area_75,'); sd = ',sd_area)) %>%
  mutate(scenario = 
           case_when(grepl('current',scenario) ~ 'Climate Change',
                     grepl('int.urb',scenario) ~ 'All Stressors',
                     grepl('int',scenario) ~ 'Agricultural Intensification',
                     grepl('exp',scenario) ~ 'Agricultural Expansion',
                     grepl('urb',scenario) ~ 'Urban Expansion')) %>%
  dplyr::select(scenario, ecoregion_id, ecoregion_name, taxon, change_in_population_abundance, change_in_habitat_area, number_of_species = num_species)



table_function <-
  function(sdm_threshold,
           climate_migration_limits) {
    file_list <-
      list.files(paste0(getwd(),'/Aggregated_CSV_Files/'),
                 pattern = 'Prop',
                 full.names = TRUE) %>%
      .[grepl(sdm_threshold,.)] %>%
      .[grepl('nopatch',.)] 
    
    if(climate_migration_limits %in% 'Yes') {
      file_list <-
        file_list %>%
        .[!grepl('Climate_Suitable',.)]
    } else if(climate_migration_limits %in% 'No') {
      file_list <-
        file_list %>%
        .[grepl('Climate_Suitable',.)]
    }
    
    # importing files
    out_df <-
      do.call(rbind, lapply(file_list, read_csv)) %>%
      mutate(taxon = gsub('/.*','',binomial)) %>%
      mutate(species = gsub('.*/','',binomial))
   
    out_df_2050_all_species <-
      out_df %>%
      filter(paste0(taxon,species) %in% paste0(species_list$taxon,species_list$species)) %>%
      filter(year %in% 2050) %>%
      left_join(., 
                realm_biome_species %>% distinct()) %>%# Merging in ecoregions species are present in
      #           read_csv('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/species_by_ecoregion_esh_maps.csv') %>%
      #             mutate(species = gsub('_[A-Z].*|_[0-9].*','',species)) %>%
      #             mutate(binomial = paste0(taxa,'/',species)) %>%
      #             dplyr::select(-c(species,taxa))) %>%
      # # filter(grepl('.int.urb',scenario)) %>% # Getting only total outcomes
      # left_join(., # Merging in biomes by ecoregion
      #           realms_biomes) %>%
      # filter(ecoregion_id > 0) %>%
      # filter(!is.na(ecoregion_id),
      #        !is.na(realm),
      #        !is.na(biome)) %>%
      # dplyr::select(-ecoregion_id) %>%
      # dplyr::select(-c(prop_cells_in_ecoregion,prop_cells,esh_map)) %>%
      distinct() %>%
      dplyr::group_by(scenario, ecoregion_id, ecoregion_name) %>%
      dplyr::summarise(median_pop = round(median(tot_pop,na.rm=TRUE),digits = 2),
                       pop_25 = round(quantile(tot_pop,.25,na.rm=TRUE),digits = 2),
                       pop_75 = round(quantile(tot_pop,.75,na.rm=TRUE),digits = 2),
                       sd_pop = round(sd(tot_pop,na.rm=TRUE),digits = 2),
                       median_area = round(median(tot_area,na.rm=TRUE),digits = 2),
                       area_25 = round(quantile(tot_area,.25,na.rm=TRUE),digits=2),
                       area_75 = round(quantile(tot_area,.75,na.rm=TRUE),digits=2),
                       sd_area = round(sd(tot_area,na.rm=TRUE),digits=2),
                       num_species = n()) %>%
      filter(!is.na(ecoregion_name)) %>%
      mutate(change_in_population_abundance = paste0(median_pop,' (',pop_25,' - ',pop_75,'); sd = ',sd_pop),
             change_in_habitat_area = paste0(median_area,' (',area_25,' - ',area_75,'); sd = ',sd_area)) %>%
      mutate(stressor = 
               case_when(grepl('current',scenario) ~ 'Climate Change',
                         grepl('int.urb',scenario) ~ 'All Stressors',
                         grepl('int',scenario) ~ 'Agricultural Intensification',
                         grepl('exp',scenario) ~ 'Agricultural Expansion',
                         grepl('urb',scenario) ~ 'Urban Expansion')) %>%
      mutate(taxon = 'All Species') %>%
      mutate(sdm_threshold = sdm_threshold,
             limits_to_climate_dispersal = climate_migration_limits) %>%
      dplyr::select(sdm_threshold, limits_to_climate_dispersal, stressor, ecoregion_id, ecoregion_name, taxon, change_in_population_abundance, change_in_habitat_area, number_of_species = num_species)
    
    out_df_2050_by_taxon <-
      out_df %>%
      filter(paste0(taxon,species) %in% paste0(species_list$taxon,species_list$species)) %>%
      filter(year %in% 2050) %>%
      left_join(., 
                realm_biome_species %>% distinct()) %>%# Merging in ecoregions species are present in
      #           read_csv('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/species_by_ecoregion_esh_maps.csv') %>%
      #             mutate(species = gsub('_[A-Z].*|_[0-9].*','',species)) %>%
      #             mutate(binomial = paste0(taxa,'/',species)) %>%
      #             dplyr::select(-c(species,taxa))) %>%
      # # filter(grepl('.int.urb',scenario)) %>% # Getting only total outcomes
      # left_join(., # Merging in biomes by ecoregion
      #           realms_biomes) %>%
      # filter(ecoregion_id > 0) %>%
      # filter(!is.na(ecoregion_id),
      #        !is.na(realm),
      #        !is.na(biome)) %>%
      # dplyr::select(-ecoregion_id) %>%
      # dplyr::select(-c(prop_cells_in_ecoregion,prop_cells,esh_map)) %>%
      distinct() %>%
      dplyr::group_by(scenario, taxon, ecoregion_id, ecoregion_name) %>%
      dplyr::summarise(median_pop = round(median(tot_pop,na.rm=TRUE),digits = 2),
                       pop_25 = round(quantile(tot_pop,.25,na.rm=TRUE),digits = 2),
                       pop_75 = round(quantile(tot_pop,.75,na.rm=TRUE),digits = 2),
                       sd_pop = round(sd(tot_pop,na.rm=TRUE),digits = 2),
                       median_area = round(median(tot_area,na.rm=TRUE),digits = 2),
                       area_25 = round(quantile(tot_area,.25,na.rm=TRUE),digits=2),
                       area_75 = round(quantile(tot_area,.75,na.rm=TRUE),digits=2),
                       sd_area = round(sd(tot_area,na.rm=TRUE),digits=2),
                       num_species = n()) %>%
      filter(!is.na(ecoregion_name)) %>%
      mutate(change_in_population_abundance = paste0(median_pop,' (',pop_25,' - ',pop_75,'); sd = ',sd_pop),
             change_in_habitat_area = paste0(median_area,' (',area_25,' - ',area_75,'); sd = ',sd_area)) %>%
      mutate(stressor = 
               case_when(grepl('current',scenario) ~ 'Climate Change',
                         grepl('int.urb',scenario) ~ 'All Stressors',
                         grepl('int',scenario) ~ 'Agricultural Intensification',
                         grepl('exp',scenario) ~ 'Agricultural Expansion',
                         grepl('urb',scenario) ~ 'Urban Expansion')) %>%
      mutate(sdm_threshold = sdm_threshold,
             limits_to_climate_dispersal = climate_migration_limits) %>%
      dplyr::select(sdm_threshold, limits_to_climate_dispersal,stressor, ecoregion_id, ecoregion_name, taxon, change_in_population_abundance, change_in_habitat_area, number_of_species = num_species)
    
    # Making long data frame
    out_df_out <-
      rbind(out_df_2050_all_species,
            out_df_2050_by_taxon) %>%
      dplyr::ungroup() %>%
      dplyr::select(-scenario) %>%
      transform(stressor = factor(stressor,
                                  levels = c('All Stressors', 'Agricultural Expansion','Agricultural Intensification','Climate Change','Urban Expansion'))) %>%
      arrange(ecoregion_id, stressor, taxon)
      
    
    # and returning
    return(out_df_out)
  }


specsens_limits <-
  table_function(sdm_threshold = 'specsens',
                 climate_migration_limits = 'Yes')

specsens_nolimits <-
  table_function(sdm_threshold = 'specsens',
                 climate_migration_limits = 'No')

prevalence_limits <-
  table_function(sdm_threshold = 'prevalence',
                 climate_migration_limits = 'Yes')

prevalence_nolimits <-
  table_function(sdm_threshold = 'prevalence',
                 climate_migration_limits = 'No')

out_df_ecoregion <-
  rbind(prevalence_limits,
        prevalence_nolimits,
        specsens_limits,
        specsens_nolimits)



###
# Repeating with biodiversity hotspots
###
# hotspots species
hotspots_species <-
  read_csv('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/species_by_biodiv_hotspot_esh_maps.csv')

# Vectors
hotspots_vect <- vect('/Users/michael/Downloads/hotspots_2016_1/hotspots_2016_1.shp')
ecoregions_vect <- vect('/Users/michael/Downloads/Terrestrial_Ecoregions/Terrestrial_Ecoregions.shp')
hotspots_map <- raster('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Ecoregions_Feb2023/Biodiversity_Hotspots_Raster.tif')
realms_map <- raster('/Users/michael/Desktop/Research/Glob2Loc/Ecoregions_Feb2023/realms_raster.tif')

hotspots_realms_df <-
  data.frame(hotspot = getValues(hotspots_map),
             realms = getValues(realms_map)) %>%
  filter(!is.na(hotspot)) %>%
  filter(!is.na(realms)) %>%
  distinct() %>%
  left_join(.,
            data.frame(realm_name = ecoregions_vect$WWF_REALM2) %>%
              distinct() %>%
              mutate(realms = 
                       case_when(grepl('Neotropic',realm_name) ~ 1,
                                 grepl('Palearctic',realm_name) ~ 2,
                                 grepl('Nearctic',realm_name) ~ 3,
                                 grepl('Indo-Malay',realm_name) ~ 5,
                                 grepl('Afrotropic',realm_name) ~ 6,
                                 grepl('Oceania',realm_name) ~ 7,
                                 grepl('Austral',realm_name) ~ 8,
                                 grepl('Antarctic',realm_name) ~ 9))) %>%
  mutate(hotspots_name = 
           case_when(hotspot %in% 1 ~ 'Atlantic Forest',
                     hotspot %in% 2 ~ 'California Floristic Province',
                     hotspot %in% 3 ~ 'California Floristic Province (Outer Limits)',
                     hotspot %in% 4 ~ 'Cape Floristic Region',
                     hotspot %in% c(5) ~ 'Caribbean Islands',
                     hotspot %in% c(6) ~ 'Caribbean Islands (Outer Limits',
                     hotspot %in% 7 ~ 'Caucasus',
                     hotspot %in% 8 ~ 'Cerrado',
                     hotspot %in% 9 ~ 'Chilean Winter Rainfall and Valdivian Forests',
                     hotspot %in% 10 ~ 'Chilean Winter Rainfall and Valdivian Forests (Outer Limits)',
                     hotspot %in% 11 ~ 'Coastal Forests of Eastern Africa',
                     hotspot %in% 12 ~ 'East Melanesian Islands',
                     hotspot %in% 13 ~ 'East Melanesian Islands (Outer Limits)',
                     hotspot %in% 14 ~ 'Eastern Afromontane',
                     hotspot %in% 15 ~ 'Guinean Forests of West Africa',
                     hotspot %in% 16 ~ 'Himalaya',
                     hotspot %in% 17 ~ 'Horn of Africa',
                     hotspot %in% 18 ~ 'Horn of Africa (Outer Limits)',
                     hotspot %in% 19 ~ 'Indo-Burma',
                     hotspot %in% 20 ~ 'Indo-Burma (Outer Limits)',
                     hotspot %in% 21 ~ 'Irano-Anatolian',
                     hotspot %in% 22 ~ 'Japan',
                     hotspot %in% 23 ~ 'Japan (Outer Limits)',
                     hotspot %in% 24 ~ 'Madagascar +',
                     hotspot %in% 25 ~ 'Madagascar +',
                     hotspot %in% 26 ~ 'Madrean Pine-Oak Woodlands',
                     hotspot %in% 27 ~ 'Maputaland-Pondoland-Albany',
                     hotspot %in% 28 ~ 'Mediterranean Basin',
                     hotspot %in% 29 ~ 'Mediterranean Basin (Outer Limits)',
                     hotspot %in% 30 ~ 'Mesoamerica',
                     hotspot %in% 31 ~ 'Mesoamerica (Outer Limits)',
                     hotspot %in% 32 ~ 'Mountains of Central Asia',
                     hotspot %in% 33 ~ 'Mountains of Southwest China',
                     hotspot %in% 34 ~ 'New Caledonia',
                     hotspot %in% 35 ~ 'New Caledonia (Outer Limits)',
                     hotspot %in% 36 ~ 'New Zealand',
                     hotspot %in% 37 ~ 'New Zealand (Outer Limits)',
                     hotspot %in% 38 ~ 'Philippines',
                     hotspot %in% 39 ~ 'Philippines (Outer Limit)',
                     hotspot %in% 40 ~ 'Polynesia-Micronesia',
                     hotspot %in% 41 ~ 'Polynesia-Micronesia (Outer Limit)',
                     hotspot %in% 42 ~ 'Southwest Australia',
                     hotspot %in% 43 ~ 'Succulent Karoo',
                     hotspot %in% 44 ~ 'Sundaland',
                     hotspot %in% 45 ~ 'Sundaland (Outer Limits)',
                     hotspot %in% 46 ~ 'Tropical Andes',
                     hotspot %in% 47 ~ 'Tumbes-Choco-Magdalena',
                     hotspot %in% 48 ~ 'Tumbes-Choco-Magdalena (Outer Limits)',
                     hotspot %in% 49 ~ 'Wallacea',
                     hotspot %in% 50 ~ 'Wallacea (Outer Limits)',
                     hotspot %in% 51 ~ 'Western Ghats and Sri Lanka',
                     hotspot %in% 52 ~ 'Forests of East Australia',
                     hotspot %in% 53 ~ 'North American Coastal Plain'))

# formatting data frame to work with function
hotspots_species <-
  left_join(hotspots_species,
            hotspots_realms_df %>%
              dplyr::select(ecoregion_id = hotspot,ecoregion_name = hotspots_name)) %>%
  distinct() %>%
  mutate(binomial = paste0(taxa,'/',species)) %>%
  dplyr::select(-c(species,taxa)) %>%
  filter(prop_cells >= .25) %>%
  dplyr::select(ecoregion_id, ecoregion_name, binomial)

# and updating name to work with function
realm_biome_species <- 
  hotspots_species

# running function
specsens_limits <-
  table_function(sdm_threshold = 'specsens',
                 climate_migration_limits = 'Yes')

specsens_nolimits <-
  table_function(sdm_threshold = 'specsens',
                 climate_migration_limits = 'No')

prevalence_limits <-
  table_function(sdm_threshold = 'prevalence',
                 climate_migration_limits = 'Yes')

prevalence_nolimits <-
  table_function(sdm_threshold = 'prevalence',
                 climate_migration_limits = 'No')

out_df_hotspots <-
  rbind(prevalence_limits,
        prevalence_nolimits,
        specsens_limits,
        specsens_nolimits) %>%
  dplyr::rename(biodiversity_hotspot = ecoregion_name) %>%
  dplyr::select(-ecoregion_id)


###
# And now for countries
coo_dat <-
  rbind(read_csv('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/COO_Data/NewCountryOfOccurrenceData.csv'),
        read_csv('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/COO_Data/coo_dat_reptiles.csv')) %>%
  mutate(binomial = paste0(taxon,'/',species)) %>%
  dplyr::select(ecoregion_id = un_code, ecoregion_name = ISO3, binomial)

# and updating name to work with function
realm_biome_species <- 
  coo_dat


specsens_limits <-
  table_function(sdm_threshold = 'specsens',
                 climate_migration_limits = 'Yes')

specsens_nolimits <-
  table_function(sdm_threshold = 'specsens',
                 climate_migration_limits = 'No')

prevalence_limits <-
  table_function(sdm_threshold = 'prevalence',
                 climate_migration_limits = 'Yes')

prevalence_nolimits <-
  table_function(sdm_threshold = 'prevalence',
                 climate_migration_limits = 'No')

out_df_countries <-
  rbind(prevalence_limits,
        prevalence_nolimits,
        specsens_limits,
        specsens_nolimits) %>%
  dplyr::rename(country_iso3c = ecoregion_name) %>%
  dplyr::select(-ecoregion_id)


###
# And saving files
write.csv(out_df_ecoregion,
          paste0('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Figures and Tables/Supplemental Data/Extended_Data_Table1_Ecoregions.csv',
                 row.names = FALSE))

write.csv(out_df_hotspots,
          paste0('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Figures and Tables/Supplemental Data/Extended_Data_Table2_biodiversity_hotspot.csv',
                 row.names = FALSE))

write.csv(out_df_countries,
          paste0('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Figures and Tables/Supplemental Data/Extended_Data_Table3_country.csv',
                 row.names = FALSE))
