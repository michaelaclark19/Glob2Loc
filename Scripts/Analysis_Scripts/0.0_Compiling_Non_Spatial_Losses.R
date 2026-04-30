#!/usr/bin/env Rscript

### Script for calculating biodiversity loss in 2020 relative to 2050
# Does this in three parts
# 1) Import files
# 2) Merge files from 2020 to 2050
# 3) Calculate outcomes
# 4) Visualising outcomes


# Libraries
library(plyr)
library(dplyr)
library(parallel)
library(stringr)

# setting working directory
setwd('/data/pubh-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity')

# number of cores
num_cores = 10

# creating function to import files
# this identifies taxa, year, assumptions on migration limit, and sdm threshold used
patch_csv_import <-
  function(i) {
    # importing csv
    out_csv <- 
      read.csv(i) %>%
      mutate(taxa = gsub(".*BAU/","",i)) %>% # Identifying taxa
      mutate(taxa = gsub("/Spatial.*","",taxa)) %>%
      mutate(year = str_extract(i,"_[0-9]{4,4}")) %>%
      mutate(year = gsub("_","", year)) %>%
      mutate(climate_migration_limits = ifelse(grepl("nomigrationlimits",i),"No_Migration_Limits","Migration_Limits")) %>%
      mutate(sdm_threshold = ifelse(grepl("prevalence", i), "prevalence", "specsens")) %>%
      mutate(binomial = gsub(".*_csv_","",i)) %>%
      mutate(binomial = gsub("[0-9]{4,4}.*","",binomial))
  }


# function for lapply across species in a year
species_lapply_function <-
  function(s) {

  cat(s,'\n')
  # subsetting full list of files to only this species
  species_files <- 
    out_files[grepl(s, out_files)] %>%
    .[grepl(y,.)]
  
  # now importing files for the species
  species_df <- 
    do.call(rbind,
            lapply(species_files,
                   patch_csv_import)) %>%
    unique() %>%
    mutate(estimated_population = as.numeric(estimated_population))
  
  # species_df_sum <-
  #   species_df %>%
  #   dplyr::group_by(species, year, land_scen, sdm_threshold, climate_migration_limits, scenario) %>%
  #   dplyr::summarise(habitat_availability_sqkm = sum(habitat_availability_sqkm),
  #                    estimated_population = sum(estimated_population),
  #                    patch_migrate = max(patch_migrate),
  #                    patch_no_migrate = max(patch_no_migrate))

  # rbinding to taxa df
  return(species_df)
} # End function to import files for a species

# List of taxa - will be looping across species within each taxa
taxa <- list.files(paste0(getwd(),'/Outputs/CSV_File_Outputs/BAU'),full.names=TRUE)

# creating out files
dir.create(paste0(getwd(),'/Outputs/Aggregated_CSV_File_Outputs'))


# now looping across taxa, amphibians, reptiles, and mammals
# need to make an exception for birds because of memory issues
# for(i in taxa) {
for(i in taxa[2]) {
  # getting files for the taxa
  out_files <- list.files(i, full.names = TRUE, pattern = 'Spatial_Pop')


for(y in c(2020,2050)) {

  # getting list of species in the year
  species_year <-
    out_files[grepl(y,out_files)] %>%
    .[grepl('nomigration',.)]

  if(grepl('Birds',i)) { # exception for birbs - they take too much memory...
    # chunking, because birds are a bit of a pain in terms of memory limits
    chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
    species_year_chunks <- chunk2(species_year,num_cores)
    # species_year_chunks <- chunk2(species_year,10)

    # now looping across species via lapply and rbind
    for(ii in 1:length(species_year_chunks)) {

      cat('Birds: Chunk ',ii,'\n')

      # and looping across species via lapply and do.call(rbind)
      out_df <- do.call(rbind,mclapply(species_year_chunks[[ii]],species_lapply_function, mc.cores=num_cores))

      # and writing file
      write.csv(out_df,
                paste0(getwd(),'/Outputs/Aggregated_CSV_File_Outputs/Stacked_File_',
                       gsub('.*BAU/','',i),
                       '_',y,'_chunk_',ii,'_nomigrationlimits.csv'),
                row.names = FALSE)
    } # end loop for bird chunks
  } else { # for all other taxa
    out_df <- do.call(rbind,mclapply(species_year,species_lapply_function, mc.cores=num_cores))

    # and writing file
    write.csv(out_df,
              paste0(getwd(),'/Outputs/Aggregated_CSV_File_Outputs/Stacked_File_',
                     gsub('.*BAU/','',i),
                     '_',y,'_nomigrationlimits.csv.csv'),
              row.names = FALSE)
  } # end else statement
} # end loop for year
} # End loop for taxa

