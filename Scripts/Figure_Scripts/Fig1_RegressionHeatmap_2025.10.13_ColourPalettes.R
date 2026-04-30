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
# Realm / biome combinations with >25% species habitat
realm_biome_species <- 
  read_csv('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/species_by_ecoregion_esh_maps.csv') %>%
  mutate(species = gsub('_[A-Z].*|_[0-9].*','',species)) %>%
  mutate(binomial = paste0(taxa,'/',species)) %>%
  dplyr::select(-c(species,taxa)) %>%
  left_join(.,realms_biomes) %>%
  dplyr::group_by(binomial, realm, biome) %>%
  dplyr::summarise(prop_cells = sum(prop_cells)) %>%
  filter(!is.na(realm),
         !is.na(biome)) %>%
  filter(prop_cells > .1)

###
# outcomes in 2050
out_df_2050 <-
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
  distinct()

### Abs outcomes in 2020
file_list <-
  list.files(paste0(getwd(),'/Aggregated_CSV_Files/'),
             pattern = 'Abs',
             full.names = TRUE) %>%
  .[grepl('prevalence',.)] %>%
  .[grepl('nopatch',.)] %>%
  .[!grepl('Climate_Suitable',.)]

# importing files
out_df_abs <-
  do.call(rbind, lapply(file_list, read_csv)) %>%
  mutate(taxon = gsub('/.*','',binomial)) %>%
  mutate(species = gsub('.*/','',binomial))

out_df_2020_2050_abs <-
  out_df_abs %>%
  filter(paste0(taxon,species) %in% paste0(species_list$taxon,species_list$species)) %>%
  filter(year %in% c(2020,2050)) %>%
  filter(grepl('.int.urb',scenario)) %>% # Getting only total outcomes
  distinct()

###
# Adding in habitat area in 2020
out_df_2050_check <-
  left_join(out_df_2050,
            out_df_2020_2050_abs %>%
              filter(year %in% 2020) %>%
              dplyr::select(binomial,
                            tot_pop_2020 = tot_pop,
                            tot_area_2020 = tot_area,
                            area_largest_patch_2020 = area_largest_patch,
                            pop_largest_patch_2020 = pop_largest_patch) %>%
              distinct())

###
# Adding in species characteristics (rarity, etc)
rarity <- read_csv('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/hab_pref_rarity_19Sep2023.csv')
out_df_2050_check <-
  left_join(out_df_2050_check,
            rarity %>%
              dplyr::select(species, Hab_Rarity, Area_Rarity) %>%
              distinct()) %>%
  mutate(Area_Rarity = case_when(is.na(Area_Rarity) ~ 'common', .default = Area_Rarity),
         Hab_Rarity = case_when(is.na(Hab_Rarity) ~ 'common',.default = Hab_Rarity))

###
# Adding in species body mass
body_mass <- 
  read_csv('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/file_transfers_aug212024/Body Mass Estimates 20February2020 Updated Taxonomy.csv') %>%
  dplyr::select(taxon = Class,
                binomial,
                est_mass_kg) %>%
  distinct() %>%
  mutate(taxon = 
           case_when(grepl('AMPHIBIA',taxon) ~ 'Amphibians',
                     grepl('AVES',taxon) ~ 'Birds',
                     grepl('MAMMALIA',taxon) ~ 'Mammals'))

rep_mass <- 
  read_csv('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/file_transfers_aug212024/reptile_body_masses.csv') %>%
  mutate(taxon = 'Reptiles') %>%
  mutate(binomial = gsub(' ','_',binomial)) %>%
  mutate(est_mass_kg = `mass (g)` / 1000) %>%
  dplyr::select(taxon, binomial, est_mass_kg) %>%
  distinct()

body_mass <-
  rbind(body_mass, 
        rep_mass) %>% 
  distinct() %>%
  mutate(binomial = paste0(taxon,'/',binomial))

out_df_2050_check <-
  left_join(out_df_2050_check,
            body_mass) %>%
  distinct()

###
# Developing glm model not controlling for realms
realms <- 
  sort(unique(out_df_2050_check$realm)) %>%
  .[!grepl('Antar',.)]

biomes <-
  c(unique(out_df_2050_check$biome) %>% .[grepl('Temperate',.)],
    unique(out_df_2050_check$biome) %>% .[grepl('Montane',.)],
    unique(out_df_2050_check$biome) %>% .[grepl('Medi',.)],
    unique(out_df_2050_check$biome) %>% .[grepl('Tropical',.)],
    unique(out_df_2050_check$biome) %>% .[grepl('Mangrove',.)],
    unique(out_df_2050_check$biome) %>% .[grepl('Boreal|Ice|Desert|Flooded',.)]
  )

#####
###
# Getting data for extent of habitat fragmentation
files <-
  list.files(path = '/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Aggregated_CSV_Files_Migration/',
             pattern = 'Abs',
             full.names = TRUE) %>%
  .[grepl('prevalence',.)]

out_df <- 
  do.call(rbind,lapply(files,read_csv)) %>%
  dplyr::select(pop_dispersal, binomial, year, num_pops = num_patches,
                area_largest_pop = area_largest_patch, pop_largest_pop = pop_largest_patch)

no_mig_files <-
  list.files(path = '/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Aggregated_CSV_Files',
             pattern = 'Abs',
             full.names = TRUE) %>%
  .[grepl('prevalence',.)] %>%
  .[!grepl('Climate',.)]

out_df_nomig <- do.call(rbind,lapply(no_mig_files,read_csv)) %>%
  filter(grepl('.int.urb',scenario)) %>%
  dplyr::select(binomial, year, num_patches,
                area_largest_patch, pop_largest_patch,
                tot_area, tot_pop) %>%
  distinct()


out_merge <-
  left_join(out_df %>% distinct(),
            out_df_nomig %>% distinct()) %>%
  mutate(num_pops = case_when(num_pops < 0 ~ 0,
                              .default = num_pops)) %>%
  mutate(pop_ratio = num_patches / num_pops) %>%
  mutate(taxa = gsub('/.*','',binomial))


# ggplot(out_merge, aes(x = taxa, y = log10(pop_ratio))) +
#   geom_point() +
#   geom_boxplot() +
#   # geom_violin() +
#   facet_wrap(.~year) +
#   coord_cartesian(ylim = c(0,2))


sum_df <-
  out_merge %>%
  dplyr::group_by(taxa,year) %>%
  dplyr::summarise(median_pop_ratio = median(pop_ratio,na.rm=TRUE),
                   mean_pop_ratio = mean(pop_ratio,na.rm=TRUE))


dat_2020 <-
  out_merge %>%
  filter(year %in% 2020)

dat_2050 <-
  out_merge %>%
  filter(year %in% 2050)

dat_ratio <-
  left_join(dat_2050 %>% dplyr::select(pop_dispersal, taxa, binomial, num_pops_2050 = num_pops, pop_largest_pop_2050 = pop_largest_pop, pop_ratio_2050 = pop_ratio, tot_area_2050 = tot_area, tot_pop_2050 = tot_pop, area_largest_pop_2050 = area_largest_pop, num_patches_2050 = num_patches),
            dat_2020 %>% dplyr::select(pop_dispersal, taxa, binomial, num_pops_2020 = num_pops, pop_largest_pop_2020 = pop_largest_pop, pop_ratio_2020 = pop_ratio, tot_area_2020 = tot_area, tot_pop_2020 = tot_pop, area_largest_pop_2020 = area_largest_pop, num_patches_2020 = num_patches)) %>%
  filter(binomial %in% paste0(species_list$taxon,'/',species_list$species)) %>%
  mutate(pop_ratio = pop_ratio_2050 / pop_ratio_2020,
         largest_pop_ratio = pop_largest_pop_2050 / pop_largest_pop_2020,
         area_largest_pop_ratio = area_largest_pop_2050 / area_largest_pop_2020, 
         num_pops_ratio = num_pops_2050 / num_pops_2020,
         num_patches_ratio = num_patches_2050 / num_patches_2020,
         tot_area_ratio = tot_area_2050 / tot_area_2020,
         tot_pop_ratio = tot_pop_2050 / tot_pop_2020) %>%
  mutate(avg_area_pop_2020 = tot_area_2020 / num_pops_2020,
         avg_area_pop_2050 = tot_area_2050 / num_pops_2050) %>%
  mutate(avg_area_pop_ratio = avg_area_pop_2050 / avg_area_pop_2020)


###
# Adding dispersal distances to the data frame
# Dispersal distance function
amp_dispersal_distance_frame <- read.csv('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Amphibian Migration Distance 29Nov2019.csv')
mam_dispersal_distance_frame <- read.csv('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/DispersalDistances_Estimated.csv')

body_mass_frame <- read.csv('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/file_transfers_aug212024/Body Mass Estimates 20February2020 Updated Taxonomy.csv')

dispersal_distance_function <-
  function(ss) {
    # Identifying mammals/birds/amphibians
    taxon <- species.frame$taxa[ss]
    
    if(taxon %in% 'Amphibians') {
      # Matching family/species/genus/etc as much as possible
      if(species.frame$binomial[ss] %in% paste0('Amphibians/',amp_dispersal_distance_frame$Binomial)) {
        dispersal_distance = mean(amp_dispersal_distance_frame$Distance[paste0('Amphibians/', amp_dispersal_distance_frame$Binomial) %in% species.frame$binomial[ss]])
      } else if(species.frame$Genus[ss] %in% amp_dispersal_distance_frame$Genus) {
        dispersal_distance = mean(amp_dispersal_distance_frame$Distance[amp_dispersal_distance_frame$Genus %in% species.frame$Genus[ss]])
      } else if(species.frame$Family[ss] %in% amp_dispersal_distance_frame$Family) {
        dispersal_distance = mean(amp_dispersal_distance_frame$Distance[amp_dispersal_distance_frame$Family %in% species.frame$Family[ss]])
      } else {
        dispersal_distance = mean(amp_dispersal_distance_frame$Distance)
      }
    } else if(taxon %in% 'Mammals') {
      # Getting dispersal distance for the species
      # This is currently here to make sure the code is working
      # Will ultimately be taken from a database
      # If body mass exists, reclaculating here using coefficients and body mass used for consistency with density estimates
      dispersal_distance = 10^(-1.810227 + 0.6628845 * log10(species.frame$est_mass_kg[ss] * 1000)) * 1000
      
      # If body mass does not exist, trying to match with binomial, genus, family, and then order
      if(is.na(species.frame$est_mass_kg[ss]) |
         max(dispersal_distance) < 0) {
        dispersal_distance = mam_dispersal_distance_frame$estimated_median_dispersal_m[paste0('Mammals/',mam_dispersal_distance_frame$binomial) %in% species.frame$binomial[ss]]
      }
      # If statement to catch pspecies not in the data frame
      if(is.na(species.frame$est_mass_kg[ss]) |
         max(dispersal_distance) < 0) {
        dispersal_distance = mean(mam_dispersal_distance_frame$estimated_median_dispersal_m[mam_dispersal_distance_frame$genus %in% species.frame$Genus[ss]], na.rm=TRUE)
      }
      if(is.na(species.frame$est_mass_kg[ss]) |
         max(dispersal_distance) < 0) {
        dispersal_distance = mean(mam_dispersal_distance_frame$estimated_median_dispersal_m[mam_dispersal_distance_frame$genus %in% species.frame$Genus[ss]], na.rm=TRUE)
      }
      if(is.na(species.frame$est_mass_kg[ss]) |
         max(dispersal_distance) < 0) {
        dispersal_distance = mean(mam_dispersal_distance_frame$estimated_median_dispersal_m[mam_dispersal_distance_frame$family %in% species.frame$Family[ss]], na.rm=TRUE)
      }
      if(is.na(species.frame$est_mass_kg[ss]) |
         max(dispersal_distance) < 0) {
        dispersal_distance = mean(mam_dispersal_distance_frame$estimated_median_dispersal_m[mam_dispersal_distance_frame$order %in% species.frame$Order[ss]], na.rm=TRUE)
      }
      if(is.na(species.frame$est_mass_kg[ss]) |
         max(dispersal_distance) < 0) {
        # if nothing else work,s then matching with 90th percentile dispersal distance across all mammals
        dispersal_distance = 4803
      }
    } else if (taxon %in% 'Birds') {
      
      # Based on taxonomies
      # And filling gaps, where needed
      dispersal_distance <-
        read.csv('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/bird_dispersal_distances.csv', stringsAsFactors = FALSE) %>%
        mutate(Scientific_Name = gsub(' ','_', Scientific_Name)) %>%
        filter(paste0('Birds/',Scientific_Name) %in% species.frame$binomial[ss]) %>%
        dplyr::select(dispersal_distance) %>%
        as.numeric()
      
      if(is.na(dispersal_distance)) {
        dispersal_distance <- 
          quantile(read.csv('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/bird_dispersal_distances.csv', stringsAsFactors = FALSE)$dispersal_distance,.9,na.rm=TRUE) 
      }
      
      if(max(dispersal_distance) < 0) {
        # using 90th percentile to fill gap if species not matched
        dispersal_distance <- 
          quantile(read.csv('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/bird_dispersal_distances.csv', stringsAsFactors = FALSE)$dispersal_distance,.9,na.rm=TRUE) 
      }
    } else if (taxon %in% 'Reptiles') {
      
      # data from: Southwood and Avens, 2010, Journal of Comparative Physiology, 'Physiological, behavioural, and ecological aspects of migration in reptiles'
      if(grepl('Thamnophis',species.frame$binomial[ss])) {
        dispersal_distance <- 17000
      } else if(grepl('Liasis',species.frame$binomial[ss])) {
        dispersal_distance <- 12000
      } else if(grepl('Iguana',species.frame$binomial[ss])) {
        dispersal_distance <- 3000
      } else if(grepl('Conolophus',species.frame$binomial[ss])) {
        dispersal_distance <- 10000
      } else if(grepl('Crocodylus',species.frame$binomial[ss])) {
        dispersal_distance <- 10000
      } else {dispersal_distance <- 0}
    }# end reptiles
    
    # returning dispersal distance
    return(dispersal_distance)
  } 


disperse_distance_frame <-
  left_join(dat_ratio,
            body_mass_frame %>%
              mutate(class = 
                       case_when(grepl('AMPHI',Class) ~ 'Amphibians',
                                 grepl('MAMM',Class) ~ 'Mammals',
                                 grepl('AVES',Class) ~ 'Birds')) %>%
              mutate(binomial = paste0(class,'/',binomial))) %>%
  mutate(dispersal_distance = NA) %>%
  dplyr::select(taxa,binomial,
                Class, Order, Family, Genus, Species,
                est_mass_kg, class, dispersal_distance) %>%
  distinct()

species.frame <- disperse_distance_frame


if(grepl('Dispersal_',list.files('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Dispersal_Distances/'))) {
  disperse_distance_frame <-
    read.csv('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Dispersal_Distances/Dispersal_Distances_2025.04.28.csv',
             stringsAsFactors = FALSE)
} else {
  
  for(ss in 1:nrow(disperse_distance_frame)) {
    disperse_distance_frame$dispersal_distance[ss] <-
      dispersal_distance_function(ss)
  }
  
  ###
  # Saving csv
  write.csv(disperse_distance_frame,
            '/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Dispersal_Distances/Dispersal_Distances_2025.04.28.csv',
            row.names = FALSE)
}


###
# Merging in dispersal distances back in
dat_ratio <-
  left_join(dat_ratio,
            disperse_distance_frame %>%
              dplyr::select(binomial, dispersal_distance))

#####
###
# Setting factors for the model - setting reference groups, etc
biomes_list <- sort(unique(out_df_2050_check$biome))
biomes_list <-
  c(biomes_list[grepl('Temperate.*Con',biomes_list)],
    biomes_list[!grepl('Temperate.*Con',biomes_list)])

realms_list <- sort(unique(out_df_2050_check$realm))
realms_list <-
  c(realms_list[grepl('Neotropic',realms_list)],
    realms_list[!grepl('Neotropic',realms_list)])


mod_dat <-
  out_df_2050_check %>% 
  filter(grepl('.int.urb',scenario)) %>%
  left_join(., 
            dat_ratio %>%
              filter(grepl('full.*migration',pop_dispersal)) %>%
              dplyr::select(-c(tot_pop_2020,tot_area_2020)) %>%
              distinct()) %>%
  filter(!is.na(tot_area)) %>% 
  filter(is.finite(tot_area)) %>%
  filter(!is.na(tot_pop)) %>%
  filter(is.finite(tot_pop)) %>%
  distinct() %>%
  mutate(extinction_debt = 
           case_when(tot_area %in% 0 ~ 1,
                     tot_area > 0 ~ 0),
         iucn_listing = 
           case_when(area_largest_patch < 25 | pop_largest_patch < 1000 ~ 0,
                     .default = 1)) %>%
  mutate(winner = 
           case_when(tot_area > 1 & tot_pop > 1 ~ 1,
                     .default = 0),
         big_loser = 
           case_when(tot_area <= .5 & tot_pop <= .5 ~ 1,
                     .default = 0)) %>%
  mutate(Hab_Rarity = case_when(Hab_Rarity %in% 'common' ~ 'Common', .default = 'Rare'),
         Area_Rarity = case_when(Area_Rarity %in% 'common' ~ 'Common', .default = 'Rare')) %>%
  mutate(rarity_combined =
           paste0('Habitat Rarity: ',Hab_Rarity,
                  ' / Area Rarity: ',Area_Rarity)) %>%
  transform(realm = factor(realm, levels = realms_list)) %>%
  transform(biome = factor(biome, levels = biomes_list)) %>%
  transform(taxon = factor(taxon, levels = c('Mammals','Amphibians','Birds','Reptiles')))

outcome_list <-
  c('tot_pop','tot_area','num_patches','largest_pop_ratio','area_largest_pop_ratio', # linear glm gaussian models
    'num_patches_ratio','num_pops_ratio','tot_pop_ratio',# more linear glm gaussian models
    'extinction_debt','winner','big_loser') # binomial models

# outcome_list <-
#   c('extinction_debt','iucn_listing', 'winner','big_loser')


#####
###
# Building the models
out_df_models <- data.frame()

for(ii in outcome_list) {
  out_lm_mod <- as.data.frame(mod_dat)
  out_lm_mod$mod_outcome <- (out_lm_mod[,names(out_lm_mod) %in% ii])
  if(ii %in% c('extinction_debt','winner','big_loser')) {
    
    if(ii %in% 'extinction_debt') {
      out_lm_realm <-
        glm(
          formula = 
            mod_outcome ~ 
            # Hab_Rarity * Area_Rarity +
            log10(est_mass_kg)*taxon + 
            rarity_combined + 
            # log2(tot_area_2020+1e-9) + 
            # log2(tot_pop_2020+1e-9) +
            realm + biome,#  - 1,
          family = 'binomial',
          data = out_lm_mod %>% 
            filter(!is.na(tot_area)) %>% 
            filter(is.finite(tot_area)) %>%
            filter(!is.na(tot_pop)) %>%
            filter(is.finite(tot_pop)) %>%
            filter(grepl('.int.urb',scenario)) %>%
            distinct()) 
    } else {
      out_lm_realm <-
        glm(
          formula = 
            mod_outcome ~ 
            # Hab_Rarity * Area_Rarity +
            log10(est_mass_kg)*taxon + 
            rarity_combined + 
            # log2(tot_area_2020+1e-9) + 
            # log2(tot_pop_2020+1e-9) +
            realm + biome,#  - 1,
          family = 'binomial',
          data = out_lm_mod %>% 
            filter(!is.na(tot_area)) %>% 
            filter(is.finite(tot_area)) %>%
            filter(tot_area > 0) %>%
            filter(!is.na(tot_pop)) %>%
            filter(is.finite(tot_pop)) %>%
            filter(grepl('.int.urb',scenario)) %>%
            distinct()) 
    }
    
    
    summary_df <-
      data.frame(summary(out_lm_realm)[['coefficients']]) %>%
      mutate(model_outcome = ii,
             estimate = Estimate, 
             std_error = `Std..Error`,
             p_value = `Pr...z..`) %>%
      dplyr::select(model_outcome,
                    estimate,
                    std_error,
                    p_value) #%>%
    # mutate(estimate = (exp(estimate) - 1))
    
    summary_df$variable_name <- row.names(summary_df)
  } else {
    out_lm_realm <-
      glm(
        formula = 
          log2(mod_outcome + 1e-9) ~ 
          # Hab_Rarity * Area_Rarity +
          log10(est_mass_kg)*taxon + 
          rarity_combined + 
          # log2(tot_area_2020+1e-9) + 
          # log2(tot_pop_2020+1e-9) +
          realm + biome,# - 1,
        data = out_lm_mod %>% 
          filter(!is.na(tot_area)) %>% 
          filter(is.finite(tot_area)) %>%
          filter(!is.na(tot_pop)) %>%
          filter(is.finite(tot_pop)) %>%
          filter(tot_area>0) %>%
          filter(!is.na(tot_area)) %>%
          filter(is.finite(log2(mod_outcome + 1e-9))) %>%
          filter(grepl('.int.urb',scenario)) %>%
          distinct()) 
    
    summary_df <-
      data.frame(summary(out_lm_realm)[['coefficients']]) %>%
      mutate(model_outcome = ii,
             estimate = Estimate, 
             std_error = `Std..Error`,
             p_value = `Pr...t..`) %>%
      dplyr::select(model_outcome,
                    estimate,
                    std_error,
                    p_value)
    
    summary_df$variable_name <- row.names(summary_df)
  }
  
  
  
  out_df_models <-
    rbind(out_df_models,
          summary_df)
}

#####
###
# And visualising the models

###
# Adding variables for reference group, etc
variables <- unique(out_df_models$variable_name)
variable_order <-
  c('Reference Taxon: Mammals',
    variables[grepl('taxon',variables)] %>% .[!grepl('est_mass_kg',.)],
    '',
    variables[grepl('est_mass',variables)],
    '','Reference Rarity: Habitat Common / Area Common',
    variables[grepl('Rarity',variables)],
    '','Reference Realm: Neotropic',
    variables[grepl('realm',variables)],
    '','Reference Biome: Temperate Coniferous Forests',
    variables[grepl('biome',variables)])

###
# Adding colour palette


# col_palette <-
#   c(brewer.pal('Reds',n=7)[5:1],
#     '#ffffff',
#     brewer.pal('RdYlBu',n=11)[7:11])


#     )

col_palette <-
  rev(c(colorRampPalette(c('#3581b9','#99e9fc','#cbf1f7'))(5),
    '#dddddd',
    colorRampPalette(c('#f2bf91','#f28379','#d95252'))(5)))

# col_palette <- brewer.pal('RdBu',n=11)

out_df_models_colour <-
  out_df_models %>%
  mutate(fill = 
           case_when(estimate <= -0.8 & grepl('tot_pop|tot_area|num_patches|prop_pop_largest_patch|prop_area_largest_patch|ratio',model_outcome) ~ col_palette[1],
                     estimate < -0.6 & grepl('tot_pop|tot_area|num_patches|prop_pop_largest_patch|prop_area_largest_patch|ratio',model_outcome) ~ col_palette[2],
                     estimate < -0.4 & grepl('tot_pop|tot_area|num_patches|prop_pop_largest_patch|prop_area_largest_patch|ratio',model_outcome) ~ col_palette[3],
                     estimate < -0.2 & grepl('tot_pop|tot_area|num_patches|prop_pop_largest_patch|prop_area_largest_patch|ratio',model_outcome) ~ col_palette[4],
                     estimate < 0 & grepl('tot_pop|tot_area|num_patches|prop_pop_largest_patch|prop_area_largest_patch|ratio',model_outcome) ~ col_palette[5],
                     # estimate < .1 & grepl('tot_pop|tot_area|num_patches|prop_pop_largest_patch|prop_area_largest_patch',model_outcome) ~ col_palette[6],
                     estimate >= .8 & grepl('tot_pop|tot_area|num_patches|prop_pop_largest_patch|prop_area_largest_patch|ratio',model_outcome) ~ col_palette[11],
                     estimate > .6 & grepl('tot_pop|tot_area|num_patches|prop_pop_largest_patch|prop_area_largest_patch|ratio',model_outcome) ~ col_palette[10],
                     estimate > .4 & grepl('tot_pop|tot_area|num_patches|prop_pop_largest_patch|prop_area_largest_patch|ratio',model_outcome) ~ col_palette[9],
                     estimate > .2 & grepl('tot_pop|tot_area|num_patches|prop_pop_largest_patch|prop_area_largest_patch|ratio',model_outcome) ~ col_palette[8],
                     estimate > 0 & grepl('tot_pop|tot_area|num_patches|prop_pop_largest_patch|prop_area_largest_patch|ratio',model_outcome) ~ col_palette[7],
                     estimate <= -2 & grepl('winner',model_outcome) ~ col_palette[1],
                     estimate < -1.5 & grepl('winner',model_outcome) ~ col_palette[2],
                     estimate < -1.0 & grepl('winner',model_outcome) ~ col_palette[3],
                     estimate < -0.5 & grepl('winner',model_outcome) ~ col_palette[4],
                     estimate < 0 & grepl('winner',model_outcome) ~ col_palette[5],
                     estimate > 2 & grepl('winner',model_outcome) ~ col_palette[11],
                     estimate > 1.5 & grepl('winner',model_outcome) ~ col_palette[10],
                     estimate > 1.0 & grepl('winner',model_outcome) ~ col_palette[9],
                     estimate > .5 & grepl('winner',model_outcome) ~ col_palette[8],
                     estimate > 0 & grepl('winner',model_outcome) ~ col_palette[7],
                     estimate <= -2 & grepl('big_loser|extinction',model_outcome) ~ col_palette[11],
                     estimate < -1.5 & grepl('big_loser|extinction',model_outcome) ~ col_palette[10],
                     estimate < -1.0 & grepl('big_loser|extinction',model_outcome) ~ col_palette[9],
                     estimate < -0.5 & grepl('big_loser|extinction',model_outcome) ~ col_palette[8],
                     estimate < 0 & grepl('big_loser|extinction',model_outcome) ~ col_palette[7],
                     estimate > 2 & grepl('big_loser|extinction',model_outcome) ~ col_palette[1],
                     estimate > 1.5 & grepl('big_loser|extinction',model_outcome) ~ col_palette[2],
                     estimate > 1.0 & grepl('big_loser|extinction',model_outcome) ~ col_palette[3],
                     estimate > .5 & grepl('big_loser|extinction',model_outcome) ~ col_palette[4],
                     estimate > 0 & grepl('big_loser|extinction',model_outcome) ~ col_palette[5])) %>%
  # Updating for big loser and extinction debt
  # mutate(fill = 
  #          case_when(estimate <= -2 & grepl('big_loser|extinction',model_outcome) ~ col_palette[11],
  #                    estimate < -1.5 & grepl('big_loser|extinction',model_outcome) ~ col_palette[10],
  #                    estimate < -1.0 & grepl('big_loser|extinction',model_outcome) ~ col_palette[9],
  #                    estimate < -0.5 & grepl('big_loser|extinction',model_outcome) ~ col_palette[8],
  #                    estimate < 0 & grepl('big_loser|extinction',model_outcome) ~ col_palette[7],
  #                    estimate > 2 & grepl('big_loser|extinction',model_outcome) ~ col_palette[1],
  #                    estimate > 1.5 & grepl('big_loser|extinction',model_outcome) ~ col_palette[2],
  #                    estimate > 1.0 & grepl('big_loser|extinction',model_outcome) ~ col_palette[3],
  #                    estimate > .5 & grepl('big_loser|extinction',model_outcome) ~ col_palette[4],
  #                    estimate > 0 & grepl('big_loser|extinction',model_outcome) ~ col_palette[5])) %>%
  # estimate >= 2.25 & !grepl('tot_pop|tot_area|num_patches|prop_pop_largest_patch|prop_area_largest_patch',model_outcome) ~ col_palette[11])) %>%
  mutate(fill = 
           case_when(p_value > 0.05 ~ '#ffffff',
                     .default = fill)) %>%
  # mutate(estimate = 
  #          case_when(p_value > 0.05 ~ NA,
  #                    .default = estimate)) %>%
  # mutate(std_error = 
  #          case_when(p_value > 0.05 ~ NA,
  #                    .default = std_error)) %>%
  mutate(colour = 'black') %>%
  mutate(colour = case_when(p_value > 0.05 ~ '#bbbbbb',.default = colour)) %>%
  mutate(significance = case_when(p_value < .001 ~ '***',
                                  p_value < .01 ~ '**',
                                  p_value < 0.05 ~ '*',
                                  .default = ''))


###
# X location
x_vals <-
  data.frame(model_outcome = 
               c('tot_area',
                 # 'tot_pop',
                 'area_largest_pop_ratio',
                 # 'largest_pop_ratio',
                 # 'num_patches_ratio',
                 'num_pops_ratio',
                 'tot_pop_ratio',
                 '',
                 'winner','big_loser','extinction_debt'),
             x_labels = 
               c('Total Area',
                 # 'Total Abundance',
                 'Area Largest\nPopulation',
                 # 'Abundance Largest\nPopulation', 
                 # 'Number of\nPatches',
                 'Number of\nPopulations',
                 'Connectivity of\nPopulations',
                 '',
                 'Winning Species',
                 'Losing Species',
                 'Potential\nExtinction Debt'),
             x_vals = NA) %>%
  mutate(x_vals = 1:nrow(.) * 2 - 1) %>%
  mutate(x_vals =
           case_when(x_vals > 16 ~ x_vals - 1,
                     .default = x_vals)) %>%
  mutate(xmin = x_vals - 1,
         xmax = x_vals + 1)

y_vals <-
  data.frame(variable_name = variable_order,
             y_labels = variable_order %>% str_remove(.,'taxon|realm|biome|rarity_combined'),
             y_vals = NA) %>%
  mutate(y_vals = 1:nrow(.))


out_df_models_colour <-
  left_join(out_df_models_colour,
            x_vals) %>%
  left_join(.,y_vals)
  
  
###
#
# cont_plot <-
#   ggplot(dat = out_df_models_colour,
#          aes(x = model_outcome, y = variable_name, fill = fill)) +
#   geom_tile(fill = out_df_models_colour$fill) +
#   # scale_fill_gradient2(limits = c(-15,15), low = 'red',high = 'blue', mid = 'white',midpoint = 0) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1)) +
#   scale_x_discrete(limits = c('tot_area','tot_pop','area_largest_pop_ratio','largest_pop_ratio',
#                               'num_patches_ratio','num_pops_ratio','tot_pop_ratio',
#                               '',
#                               'winner','big_loser','extinction_debt'),
#                    labels = c('Total Area','Total Abundance','Area of Largest Population', 'Abundance of Largest Population', 
#                               'Number of Patches','Number of Populations','Connectivity of Populations',
#                               '',
#                               'Winning Species','Losing Species','Potential Extinction Debt')) +
#   scale_y_discrete(limits = rev(variable_order),
#                    labels = rev(variable_order) %>% str_remove(.,'taxon|realm|biome|rarity_combined')) +
#   labs(y = NULL, x = NULL) +
#   geom_text(aes(label = paste0(round(estimate,digits = 2),' (',round(std_error, digits = 2),')', significance)), size = 6 / .pt, colour = out_df_models_colour$colour) +
#   theme(text = element_text(size = 6)) +
#   theme(axis.text.y = element_text(angle = 45))

mod_type_text <-
  data.frame(x = c(mean(x_vals$x_vals[!grepl('Winning|Losing|Extinction',x_vals$x_labels)]) - 1,
                   mean(x_vals$x_vals[grepl('Winning|Losing|Extinction',x_vals$x_labels)])),
             y = rep(max(y_vals$y_vals),2) + 2,
             label = c('Continuous GLMs','Binomial GLMs'))

# 
# cont_plot <-
#   ggplot(dat = out_df_models_colour,
#          aes(fill = fill)) +
#   geom_rect(aes(xmin = xmin, xmax = xmax, ymin = y_vals - .5, ymax = y_vals + .5),fill = out_df_models_colour$fill) +
#   # scale_fill_gradient2(limits = c(-15,15), low = 'red',high = 'blue', mid = 'white',midpoint = 0) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 0, vjust = 1,hjust = .5)) +
#   scale_x_continuous(breaks = unique(x_vals$x_vals),
#                      labels = unique(x_vals$x_labels),
#                      expand = c(0,0),
#                      position = 'top') +
#   scale_y_continuous(breaks = unique(y_vals$y_vals),
#                    labels = (y_vals$y_labels),
#                    transform = 'reverse',
#                    # limits = c(-5,max(y_vals$y_vals)),
#                    # limits = c(45,2),
#                    expand = expansion(mult = c(.1,0))) +
#   labs(y = NULL, x = NULL) +
#   geom_text(aes(x = x_vals, y = y_vals, label = paste0(round(estimate,digits = 2),' (',round(std_error, digits = 2),')', significance)), size = 6 / .pt, colour = out_df_models_colour$colour) +
#   theme(text = element_text(size = 6)) +
#   theme(axis.text.y = element_text(angle = 45)) +
#   theme(axis.ticks.x = element_blank(),
#         axis.ticks.y = element_blank()) +
#   theme(axis.line = element_blank()) +
#   geom_vline(xintercept = x_vals$x_vals[x_vals$model_outcome %in% ''] - .5,
#              colour = 'grey',
#              linetype = 2) + 
#   geom_text(dat = mod_type_text, aes(x = x, y = y, label = label, fill = NULL, colour = NULL), size = 6/.pt) +
#   geom_segment(aes(x = 0,
#                    xend = 0,
#                    y = 40.5,
#                    yend = 41.5)) +
#   geom_segment(aes(x = 14,
#                    xend = 14,
#                    y = 40.5,
#                    yend = 41.5)) +
#   geom_segment(aes(x = 0,
#                    xend = 14,
#                    y = 41.5,
#                    yend = 41.5)) +
#   geom_segment(aes(x = 15,
#                    xend = 15,
#                    y = 40.5,
#                    yend = 41.5)) +
#   geom_segment(aes(x = 21,
#                    xend = 21,
#                    y = 40.5,
#                    yend = 41.5)) +
#   geom_segment(aes(x = 15,
#                    xend = 21,
#                    y = 41.5,
#                    yend = 41.5))
  # coord_cartesian(ylim = c(-5,40))
  # geom_hline(yintercept = y_vals$y_vals[y_vals$variable_name %in% ''] + .5,
  #            colour = 'grey',
  #            linetype = 2)


# cont_plot
# 
# ggsave('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Figures and Tables/Figures/Heatmaps/Regression_Results_2025.05.07.pdf',
#        width = 18.2,
#        height = 25,
#        units = 'cm')




sum(dat_ratio$area_largest_pop_ratio[grepl('full.*mig',dat_ratio$pop_dispersal)] <=0.5)

###
# setting up data structure for plotting in base R
out_df_models_colour <-
  out_df_models_colour %>%
  mutate(x_position =
           case_when(grepl('tot_area', model_outcome) ~ 1,
                     grepl('area_largest_pop',model_outcome) ~ 3,
                     grepl('num_pops_ratio',model_outcome) ~ 5,
                     grepl('tot_pop_ratio',model_outcome) ~ 7,
                     grepl('winner',model_outcome) ~ 10,
                     grepl('loser',model_outcome) ~ 12,
                     grepl('extinction',model_outcome) ~ 14)) %>%
  # mutate(y_vals = case_when(!grepl('\\btaxon[A-Z]',variable_name) ~ y_vals + 1,
  #                           .default = y_vals)) %>%
  mutate(x_left = x_position - 1,
         x_right = x_position + 1,
         y_bottom = y_vals - 1,
         y_top = y_vals) %>%
  filter(!is.na(x_left),
         !is.na(y_top))


pdf(paste0(getwd(),'/Figures and Tables/Figures/Table1_TableHeatMap_',Sys.Date(),'.pdf'),
    width = 25 / 2.54,
    height = 8)

# Setting margins
par(mar = c(1,11,3,0))

# Creating empty plot
plot(NA,
     xlim = c(min(out_df_models_colour$x_left, na.rm = TRUE),
              max(out_df_models_colour$x_right, na.rm = TRUE)),
     ylim = rev(c(-min(out_df_models_colour$y_bottom, na.rm = TRUE),
              -max(out_df_models_colour$y_top, na.rm = TRUE))),
     xlab = '',
     ylab = '',
     axes = FALSE)

# Adding rectangles
for(ii in 1:nrow(out_df_models_colour)) {
  rect(xleft = out_df_models_colour$x_left[ii],
       xright = out_df_models_colour$x_right[ii],
       ybottom = -out_df_models_colour$y_bottom[ii],
       ytop = -out_df_models_colour$y_top[ii],
       col = out_df_models_colour$fill[ii],
       border = out_df_models_colour$fill[ii])
}

# Adding text
text(x = out_df_models_colour$x_position,
     y = -out_df_models_colour$y_vals + .5,
     label = paste0(round(out_df_models_colour$estimate,digits = 2),' (',round(out_df_models_colour$std_error, digits = 2),')', out_df_models_colour$significance),
     col = out_df_models_colour$colour,
     cex = .6,
     adj = c(.5,.5))

axis(side = 3,
     at = c(1,3,5,7),
     labels = x_vals$x_labels[1:4],
     cex.axis = .6)

axis(side = 3,
     at = c(10,12,14),
     labels = x_vals$x_labels[6:8],
     cex.axis = .6,
     xpd = TRUE)

# axis(side = 2,
#      at = -seq(from = min(out_df_models_colour$y_vals) - 1.5,
#                to = max(out_df_models_colour$y_vals) - .5,
#                by = 1),
#      # labels = unique(out_df_models_colour$y_labels),
#      labels = (variable_order) %>% str_remove(.,'taxon|realm|biome|rarity_combined'),
#      cex.axis = .6,
#      xpd = TRUE,
#      las = 1)

axis(side = 2,
     at = 
       -seq(from = min(out_df_models_colour$y_vals[grepl('^taxon[A-Z]',out_df_models_colour$variable_name)]) - 1,
           to = max(out_df_models_colour$y_vals[grepl('^taxon[A-Z]',out_df_models_colour$variable_name)]),
           by = 1) + .5,
     # labels = unique(out_df_models_colour$y_labels),
     labels = (variable_order) %>% str_remove(.,'taxon|realm|biome|rarity_combined') %>% .[1:4] %>% gsub(' Taxon','',.),
     cex.axis = .6,
     xpd = TRUE,
     las = 1)

axis(side = 2,
     at = 
       -seq(from = min(out_df_models_colour$y_vals[grepl('est_mass_kg',out_df_models_colour$variable_name)]),
            to = max(out_df_models_colour$y_vals[grepl('est_mass_kg',out_df_models_colour$variable_name)]),
            by = 1) + .5,
     # labels = unique(out_df_models_colour$y_labels),
     labels = (variable_order) %>% str_remove(.,'taxon|realm|biome|rarity_combined') %>% .[grepl('est_mass',.)],
     cex.axis = .6,
     xpd = TRUE,
     las = 1)

axis(side = 2,
     at = 
       -seq(from = min(out_df_models_colour$y_vals[grepl('Rarity',out_df_models_colour$variable_name)]) - 1,
            to = max(out_df_models_colour$y_vals[grepl('Rarity',out_df_models_colour$variable_name)]),
            by = 1) + .5,
     # labels = unique(out_df_models_colour$y_labels),
     labels = (variable_order) %>% str_remove(.,'taxon|realm|biome|rarity_combined') %>% .[grepl('Rarity',.)] %>% gsub(' Rarity','',.) %>% gsub(':','',.),
     cex.axis = .6,
     xpd = TRUE,
     las = 1)

axis(side = 2,
     at = 
       -seq(from = min(out_df_models_colour$y_vals[grepl('Realm',out_df_models_colour$variable_name, ignore.case = TRUE)]) - 1,
            to = max(out_df_models_colour$y_vals[grepl('Realm',out_df_models_colour$variable_name, ignore.case = TRUE)]),
            by = 1) + .5,
     # labels = unique(out_df_models_colour$y_labels),
     labels = (variable_order) %>% .[grepl('Realm',.,ignore.case = TRUE)]  %>% str_remove(.,'taxon|realm|biome|rarity_combined') %>% gsub(' Realm','',.),
     cex.axis = .6,
     xpd = TRUE,
     las = 1)

axis(side = 2,
     at = 
       -seq(from = min(out_df_models_colour$y_vals[grepl('biome',out_df_models_colour$variable_name, ignore.case = TRUE)]) - 1,
            to = max(out_df_models_colour$y_vals[grepl('biome',out_df_models_colour$variable_name, ignore.case = TRUE)]),
            by = 1) + .5,
     # labels = unique(out_df_models_colour$y_labels),
     labels = 
       (variable_order) %>% 
       .[grepl('biome',.,ignore.case = TRUE)] %>% 
       str_remove(.,'taxon|realm|biome|rarity_combined') %>% 
       gsub(' Biome','',.) %>%
       str_replace('Temperate','Temp.') %>%
       str_replace('Tropical and Subtropical','Trop. & Subtrop.') %>%
       str_replace('Grasslands, Savannas and Shrublands','Grass, Savanna, and Shrub') %>%
       str_replace('Mediterranean Forests, Woodlands and Scrub','Mediterranean Forests and Scrub'),
     cex.axis = .6,
     xpd = TRUE,
     las = 1)

segments(y0 = -max(out_df_models_colour$y_vals) - 1,
         y1 = -max(out_df_models_colour$y_vals) - 1,
         x0 = 0,
         x1 = 8,
         xpd = TRUE)

segments(y0 = -max(out_df_models_colour$y_vals) - 1,
         y1 = -max(out_df_models_colour$y_vals) - .5,
         x0 = 0,
         x1 = 0,
         xpd = TRUE)

segments(y0 = -max(out_df_models_colour$y_vals) - 1,
         y1 = -max(out_df_models_colour$y_vals) - .5,
         x0 = 8,
         x1 = 8,
         xpd = TRUE)


segments(y0 = -max(out_df_models_colour$y_vals) - 1,
         y1 = -max(out_df_models_colour$y_vals) - 1,
         x0 = 9,
         x1 = 15,
         xpd = TRUE)

segments(y0 = -max(out_df_models_colour$y_vals) -1,
         y1 = -max(out_df_models_colour$y_vals) - .5,
         x0 = 9,
         x1 = 9,
         xpd = TRUE)

segments(y0 = -max(out_df_models_colour$y_vals) -1,
         y1 = -max(out_df_models_colour$y_vals) - .5,
         x0 = 15,
         x1 = 15,
         xpd = TRUE)

text(x = 4,
     y = -max(out_df_models_colour$y_vals) - 1.5,
     label = 'Continuous Models',
     cex = .6,
     xpd = TRUE,
     adj = c(.5,1))

text(x = 12,
     y = -max(out_df_models_colour$y_vals) - 1.5,
     label = 'Binomial Models',
     cex = .6,
     xpd = TRUE,
     adj = c(.5,1))

dev.off()
