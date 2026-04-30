#!/usr/bin/env Rscript

#####
# script to rbind species-level csv files
#####

# libraries
library(plyr)
library(dplyr)
library(parallel)
library(readr)

# working directory
setwd('/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/Outputs/Managed_CSV_Files/')
dir.create(gsub('/Managed_CSV_Files','/Aggregated_CSV_Files',getwd()))

# list of taxa
taxa_list <- c('Amphibians','Mammals','Reptiles','Birds')

# thresholds
thresh_list <- c('prevalence','specsens')

# impact_type
impact_type <- c('Abs_Impacts','Prop_Impacts')

read_csv_function <- 
	function(iii) {

		tmp <- read_csv(iii)
		tmp <- tmp[,names(names_df)]
		return(tmp)
	}


# looping through threshold types
for(thresh in thresh_list) {
  # looping through taxa types
  for(tt in taxa_list) {
    # looping through impact types (absolute impacts or prop impacts)
    for(ii in impact_type) {
	    cat('Starting:',thresh, tt, ii, '\n')
     # file_list <- 
     #   list.files(paste0(getwd(),'/',ii,'/',tt), full.names = TRUE) %>% 
     #   .[grepl(thresh,.)] %>% 
     #   .[!grepl('nopatch',.)]
      
     # out_df_migrate <- do.call(rbind, mclapply(file_list, read_csv, mc.cores = 3))
      
      file_list <- 
        list.files(paste0(getwd(),'/',ii,'/',tt), full.names = TRUE) %>% 
        .[grepl(thresh,.)] %>% 
        .[grepl('nopatch',.)] %>%
	.[grepl('Climate_Suitable',.)]

names_df <- read_csv(file_list[1])
      
     out_df_nomigrate <- do.call(rbind, mclapply(file_list, read_csv_function, mc.cores = 10))
      

      write.csv(out_df_nomigrate,
                paste0(gsub('/Managed_CSV_Files','/Aggregated_CSV_Files',getwd()),'/',ii,'_',tt,'_',thresh,'_Climate_Suitable_nopatchmigration.csv'),
                row.names = FALSE)
     # write.csv(out_df_migrate,
     #           paste0(gsub('/Managed_CSV_Files','/Aggregated_CSV_Files',getwd()),'/',ii,'_',tt,'_',thresh,'_patchmigration.csv'),
     #           row.names = FALSE)
    } # end impact type
  } # end taxa type
} # end threshold type

###
# end script!







