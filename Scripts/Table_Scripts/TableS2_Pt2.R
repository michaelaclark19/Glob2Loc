#!/usr/bin/env Rscript

# Libraries
library(plyr)
library(dplyr)
library(parallel)
library(data.table)

# Setting wd
setwd('...')

file_list <- list.files(paste0(getwd(),'/Analyses/SDM_Analyses/Species_Results'),full.names=TRUE)
cat(length(file_list),'\n')

tmp <- rbindlist(mclapply(file_list, fread,mc.cores=9), use.names = TRUE, fill = TRUE)
# tmp <- do.call(rbind,mclapply(file_list,read.csv,mc.cores=9))


tmp_2020 <- 
  tmp %>% 
  filter(year %in% 2020) %>%
  dplyr::rename(num_cells_2020 = num_cells)

tmp_merge <- 
  left_join(tmp,
            tmp_2020 %>%
              dplyr::select(-year)) %>%
  mutate(ratio = num_cells / num_cells_2020,
         log_ratio = log2(num_cells / num_cells_2020))


tmp_summary <-
  tmp_merge %>%
  dplyr::group_by(year, taxon, threshold, climate_migration_limits) %>%
  #filter(is.finite(ratio)) %>%
  dplyr::summarise(median_ratio = median(ratio, na.rm = TRUE),
                   ratio_25th_percentile = quantile(ratio,.25,na.rm=TRUE),
                   ratio_75th_percentile = quantile(ratio,.75,na.rm=TRUE),
                   ratio_5th_percentile = quantile(ratio,.05,na.rm=TRUE),
                   ratio_95th_percentile = quantile(ratio,.95,na.rm=TRUE),
                   se_median_ratio = sd(ratio, na.rm = TRUE),
                   #mean_ratio = mean(ratio, na.rm = TRUE),
                   median_log_ratio = median(log_ratio, na.rm = TRUE),
                   log_ratio_25th = quantile(log_ratio,.25,na.rm=TRUE),
                   log_ratio_75th = quantile(log_ratio,.75,na.rm=TRUE),
                   num_species = n()) %>%
                   #mean_log_ratio = mean(log_ratio, na.rm = TRUE)) %>%
  filter(year %in% 2050) %>%
  mutate(median_log_ratio = 2^median_log_ratio) %>%
  mutate(se_median_ratio = se_median_ratio / sqrt(num_species))

tmp_summary <-
  as.data.frame(tmp_summary %>% 
                  ungroup() %>%
                  dplyr::select(taxon, threshold, climate_migration_limits,
                                median_ratio, ratio_25th_percentile, ratio_75th_percentile, ratio_5th_percentile, ratio_95th_percentile,se_median_ratio,num_species) %>%
                  mutate(median_ratio = round(median_ratio,digits=2),
                         ratio_25th_percentile = round(ratio_25th_percentile,digits=2),
                         ratio_75th_percentile = round(ratio_75th_percentile,digits=2),
                         ratio_5th_percentile = round(ratio_5th_percentile,digits=2),
                         ratio_95th_percentile = round(ratio_95th_percentile,digits=2),
                         se_median_ratio = round(se_median_ratio, digits = 3))) %>%
  #mutate(median_ratio = paste0(median_ratio, '(',ratio_25th_percentile, ' - ', ratio_75th_percentile, ')')) %>%
  #mutate(percentiles_5th_to_95th = paste0(ratio_5th_percentile, ' - ',ratio_95th_percentile)) %>%
  dplyr::select(taxon, threshold, climate_migration_limits,
                median_ratio, ratio_25th_percentile, ratio_75th_percentile, num_species)
  

tmp_summary



write.csv(tmp_summary,
          paste0(getwd(),'/Analyses/SDM_Analyses/Species_Results_Aggregated_2026.04.16.csv'),
          row.names = FALSE)
