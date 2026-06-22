#!/usr/bin/env Rscript

### reptile species richness

# libraries
library(raster)
library(plyr)
library(dplyr)
library(parallel)

# setting directory
setwd('/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity')

# getting maps
reptile_maps <- list.files(paste0(getwd(),'/ESH_Tifs_12Oct/Mammals'), full.names = TRUE)

# getting template map
template <- raster('Ecoregions_Feb2023/biomes_realms_raster.tif')
template[!is.na(template)] <- 0
# template_df <- data.frame(richness = getValues(template))
# 
# # making function
# richness_fun <-
#   function(i) {
#     tmp <- raster(i)
#     tmp[is.na(tmp)] <- 0
#     tmp <- raster::extend(tmp, template, value = 0)
#     template <- template+tmp
#     return(template)
#   }
# 
# # testing with chameleons
# reptile_maps_check <- reptile_maps[grepl('Calumma_',reptile_maps)]
# 
# # checking time
# # about 3 species/min
# t1 = Sys.time()
# for(i in reptile_maps_check[1:3]) {
#   tmp <- raster(i)
#   tmp[is.na(tmp)] <- 0
#   tmp <- raster::extend(tmp, template, value = 0)
#   template_df <-
#     template_df %>%
#     mutate(tmp_richness = getValues(tmp)) %>%
#     mutate(richness = richness + tmp_richness) %>%
#     dplyr::select(-tmp_richness)
# }
# Sys.time() - t1
# 
# # converting into function that can be run in parallel
# 
# t1 = Sys.time()
# 
# # getting template map
# template <- raster('Ecoregions_Feb2023/biomes_realms_raster.tif')
# template[!is.na(template)] <- 0
# template_df <- data.frame(richness = getValues(template))

# function
richness_fun <- 
  function(z) {
    
    # creating template df
    template_df <- data.frame(richness = getValues(template))
    
    # updating species_list
    species_list <- chunk_list[[z]]
    
    # # now looping through species
    for(i in species_list) {
      cat(i,'\n')
      tmp <- raster(i)
      tmp[is.na(tmp)] <- 0
      tmp <- raster::extend(tmp, template, value = 0)
      template_df <-
        template_df %>%
        mutate(tmp_richness = getValues(tmp)) %>%
        mutate(richness = richness + tmp_richness) %>%
        dplyr::select(-tmp_richness)
    }
    # returning
    # out_raster <- raster(matrix(template_df$richness, byrow = TRUE, nrow = template@nrows))
    # extent(out_raster) <- extent(template)
    # crs(out_raster) <- crs(template)
    
    # saving data frame
    write.csv(template_df,
              paste0(getwd(),'/Other Data Inputs/Species Richness Rasters/mammal_richness_df_',z,'.csv'),
              row.names = FALSE)
    
    # return(template_df)
  }

# tmp <- lapply(reptile_maps_check[1:3],richness_fun)

# checking with chameleons
# it works...
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
# chunk_list <- chunk2(reptile_maps_check,5)
# tmp <- mclapply(chunk_list,richness_fun,mc.cores=5)

# chunking into 5 parts


# chunk_list <- chunk2(reptile_maps,20)
chunk_list <- chunk2(reptile_maps,20)

# and checking time of parallel
t1 = Sys.time()
out_list <- mclapply(chunk_list,richness_fun,mc.cores = length(chunk_list))
Sys.time() - t1

