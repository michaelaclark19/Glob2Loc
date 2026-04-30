###
# Script to identify number of species in analysis
# Based on quality of SDM
# Accuracy of SDMS
# And reasons for species exclusion

# libraries
library(plyr)
library(dplyr)


# list of files for overall model accuracy
files_list <- 
  list.files(path = '/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/ESH_RCPs/Weighted_Threshold_Accuracy',
             full.names = TRUE)
files_out <-
  do.call(c,lapply(files_list, list.files, full.names = TRUE))

# function to import
import_fun <-
  function(ii) {
    return(read.csv(ii) %>%
             dplyr::select(taxa, species, threshold, region, tot_accuracy,
                           in_both, esh_not_sdm, sdm_not_esh, out_both))
  }

# importing functions
stacked_files <-
  do.call(rbind, lapply(files_out %>% .[!grepl('error',.)], import_fun))

###
# getting summary stats

# by taxa
stacked_out <-
  stacked_files %>%
  dplyr::group_by(taxa, threshold, region) %>%
  dplyr::summarise(mean_acc = mean(tot_accuracy,na.rm=TRUE),
                   median_acc = median(tot_accuracy,na.rm=TRUE),
                   sd_acc = sd(tot_accuracy,na.rm=TRUE),
                   count = n()) %>%
  mutate(se_acc = sd_acc / sqrt(count))

# overall
stacked_files %>%
  dplyr::group_by(threshold, region) %>%
  dplyr::summarise(mean_acc = mean(tot_accuracy,na.rm=TRUE),
                   median_acc = median(tot_accuracy,na.rm=TRUE),
                   sd_acc = sd(tot_accuracy,na.rm=TRUE),
                   count = n()) %>%
  mutate(se_acc = sd_acc / sqrt(count))

# writing file
write.csv(stacked_files,
          '/data/ouce-glob2loc/pubh0329/stacked_sdm_accuracy.csv',
          row.names = FALSE)

tmp <- read.csv('/data/ouce-glob2loc/pubh0329/stacked_sdm_accuracy.csv')

stacked_out <-
  tmp %>%
  filter(!grepl('eco',region)) %>%
  dplyr::group_by(taxa, threshold) %>%
  dplyr::summarise(mean_acc = mean(tot_accuracy,na.rm=TRUE) * 100,
                   median_acc = median(tot_accuracy,na.rm=TRUE) * 100,
                   sd_acc = sd(tot_accuracy,na.rm=TRUE) * 100,
                   se_acc = (sd(tot_accuracy,na.rm=TRUE)/ sqrt(n()))* 100,
                   number_species = n())



write.xlsx(stacked_out,
           '/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Checking_SDM_Accuracy/Table_S1.xlsx',
           row.names = FALSE)


