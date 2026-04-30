#####
###
# Making data for supplemental table 4

###
# Libraries
library(plyr)
library(dplyr)
library(readr)
library(readxl)

###
# Setting working directory
setwd("/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity")

#####
###
# Background data management

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
# getting files containing migration estimates
file_list <-
  list.files(paste0(getwd(),'/Aggregated_CSV_Files_Migration/'),
             pattern = 'Prop',
             full.names = TRUE) %>%
  # .[grepl('prevalence',.)] %>%
  .[grepl('nopatch',.)] %>%
  .[!grepl('Climate_Suitable',.)]

# function for importing files
import_fun <-
  function(ii) {
    read_csv(ii) %>%
      mutate(threshold = str_extract(ii,'prevalence|specsens'))
  }

# importing files
out_df <-
  do.call(rbind,lapply(file_list,import_fun))

# summarising outcomes
out_df_sum <-
  out_df %>%
  filter(year %in% 2050) %>%
  mutate(num_patches = 
           case_when(is.finite(num_patches) ~ num_patches,
                     !is.finite(num_patches) ~ NA)) %>% 
  mutate(taxon = str_remove(binomial,'/.*')) %>%
  filter(binomial %in% paste0(species_list$taxon,'/',species_list$species)) %>%
  dplyr::group_by(taxon, threshold, pop_dispersal) %>%
  dplyr::summarise(#fifth_patches = quantile(num_patches, .05, na.rm = TRUE),
                   twentyfifth_patches = quantile(num_patches, .25, na.rm = TRUE),
                   median_patches = median(num_patches,na.rm=TRUE),
                   seventyfifth_patches = quantile(num_patches, .75, na.rm = TRUE)) %>%
  mutate(pop_dispersal = str_remove(pop_dispersal,'_patch.*'))
                   # ninetyfifth_patches = quantile(num_patches, .95, na.rm = TRUE),
                   # mean_pathces = mean(num_patches,na.rm=TRUE))

write.xlsx(out_df_sum,
          paste0(getwd(),'/Figures and Tables/Tables/TableS4_',Sys.Date(),'.xlsx'),
          rowNames = FALSE)
