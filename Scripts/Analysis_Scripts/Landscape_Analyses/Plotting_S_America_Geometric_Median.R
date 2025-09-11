#####
###
# Looking at geometric mean losses across taxa in the amazon

###
# libraries
library(raster)
library(terra)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(scales)

###
# getting shapefile
template <- raster('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Landscape_Analyses/S_America_prevalence_Amphibians_2050_future.extent.esh.exp.int.urb.tif_pop_density_geometric_median_2024-09-15.tif')
my_sf <- read_sf('/Users/michael/Downloads/TM_WORLD_BORDERS-0.3/TM_WORLD_BORDERS-0.3.shp')
my_sf_reprojected <- st_transform(my_sf$geometry, crs(template))

# country
country_map <- 
  raster('/Users/michael/Desktop/MollweideCountryID_1.5km.tif') %>%
  raster::crop(.,template)

country_map[!(country_map %in% 76)] <- NA

###
# function for plotting
plot_fun <-
  function(yy) {
    # raster
    tmp_raster <- raster(yy)
    tmp_raster[is.na(country_map)] <- NA
    
    # getting quantiles
    quants <- max(abs(quantile(getValues(tmp_raster),c(.075,.925),na.rm=TRUE)))
    quants_cutoff <- (quantile(getValues(tmp_raster),c(.1,.9),na.rm=TRUE))
    
    # converting
    # tmp_raster <- tmp_raster < quants_cutoff[1]
    tmp_raster <- squish(tmp_raster, range = c(-quants, quants))
    
    # plotting
    plot(tmp_raster,
         col = brewer.pal('PiYG', n = 11),
         # zlim = c(-quants,quants),
         box = FALSE, 
         axes = FALSE,
         main = str_extract(yy,'Birds|Mammals|Amphibians|Reptiles'))
    
    # adding shape file
    plot(my_sf_reprojected, add = TRUE, lwd = .5)
  } # End function



###
# importing and plotting files
dev.off()
par(mfrow = c(2,2))
par(mar = c(1,1,1,1))

# raster list
raster_list <-
  list.files(path = '/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Landscape_Analyses/',
             pattern = 'geometri',
             full.names = TRUE) %>%
  .[grepl('median',.)] %>%
  .[grepl('.exp.int.urb',.)]

# plotting
lapply(raster_list, plot_fun)


###
# getting total loss
library(raster)
library(plyr)
library(dplyr)
library(terra)
library(stringr)


###
# files
file_list <-
  list.files(path = '/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/Analyses/Landscape_Analyses/',
             full.names = TRUE) %>%
  .[!grepl('geometri',.)] %>%
  .[!grepl('richness',.)] %>%
  .[grepl('hab',.)]

### scenarios
scen_list <-
  c('current.extent.esh.tif',
    'future.extent.esh.exp.tif',
    'future.extent.esh.urb.tif',
    'future.extent.esh.int.tif')

###
# 2020
rast_2020 <-
  rast(file_list %>%
         .[grepl(2020,.)] %>%
         .[grepl('prevalence',.)] %>%
         .[grepl('future',.)])

rast_2020 <- terra::app(rast_2020,sum)

### looping
out_list <- list()
for(ss in scen_list) {
  cat(ss,'\n')
  # files for scenario
  scen_files <- 
    file_list %>%
    .[grepl(ss,.)] %>%
    .[grepl('prevalence',.)] %>%
    .[grepl(2050,.)]
  
  if(length(scen_files) > 4) {
    scen_files <-
      scen_files %>%
      .[!grepl('Bird.*9-10',.)]
  }
  
  cat(length(scen_files),'\\n')
  
  # stacking
  raster_stack <- rast(scen_files)
  raster_out <- app(raster_stack, sum)
  raster_out <- raster_out - rast_2020
  
  # adding to list
  out_list[[ss]] <- raster_out
}

# country map
country_map <- rast('/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/Global Mollweide Maps/MollweideCountryID_1.5km.tif')
country_map <- terra::crop(country_map, raster_out)
country_map[!(country_map %in% 76)] <- NA

out_df <- data.frame()
country_raster <- raster(country_map)
# getting data frame
for(rr in names(out_list) %>% .[!grepl('.int.tif|.exp.int.urb',.)]) {
    tmp_rast <- raster(out_list[[rr]])
    tmp_rast[is.na(country_raster)] <- NA
    
    tmp_values <- getValues(tmp_rast)
    quants <- quantile(tmp_values, .2, na.rm = TRUE)
    tmp_values <- tmp_values < quants
    
    if(rr %in% names(out_list)[1]) {
      out_df <- data.frame(rr = tmp_values)
      names(out_df) <- rr
    } else {
      out_df[,rr] <- tmp_values
    }
}


out_df_keep <- out_df
out_df <- out_df_keep

out_df <-
  out_df %>%
  mutate(current.extent.esh.tif = case_when(current.extent.esh.tif %in% 1 ~ 'Climate'),
         future.extent.esh.exp.tif = case_when(future.extent.esh.exp.tif %in% 1 ~ 'Ag Exp'),
         future.extent.esh.urb.tif = case_when(future.extent.esh.urb.tif %in% 1 ~ 'Urb Exp')) %>%
  mutate(tot_stress = paste(current.extent.esh.tif,future.extent.esh.exp.tif,future.extent.esh.urb.tif))
  

#
tmp_table <- 
  table(out_df$tot_stress) %>% 
  data.frame(.) %>%
  mutate(num_stressors = 3 - str_count(Var1,'NA')) %>%
  filter(num_stressors > 0) %>%
  mutate(Freq = Freq / sum(Freq)) %>%
  mutate(Freq = Freq * 100) %>%
  arrange(Freq)

tmp_table_sum <-
  tmp_table %>%
  dplyr::group_by(num_stressors) %>%
  dplyr::summarise(Freq = sum(Freq))

tmp_table_sum

tmp_table %>%
  filter(num_stressors > 1) %>%
  mutate(Freq = Freq / sum(Freq)) %>%
  mutate(Freq = Freq * 100) %>%
  arrange(Freq)
  
