# Script for forecasting spatial pattern of urbanization ----
# 2050 projection from Huang et al 2019 ERL (SSP2 data from Seto Lab)
# 2010 data taken from GlobCov
# Same change in each 5-year interval for each country
# Cells flip based on proximity to current urban areas
# Urban 2050 available from here: https://figshare.com/articles/Global_Urban_Land_Expansion_by_2050/7897010

# Libraries ----
library(raster)
library(plyr)
library(dplyr)

# (0) setting working directory
setwd("/Users/maclark/Desktop/Multiple Stresses of Biodiversity")

# (1) Importing  rasters
# Urban rasters
urban.2010 <- raster(paste0(getwd(),'/ESA LandCov Maps/Urban2010_Corrected_GlobCov.tif'))
urban.2050 <- raster(paste0(getwd(),'/Raw Urban 2050 Projections/urban-ssp2.tif'))
# Travel to urban areas
travel.time <- raster(paste0(getwd(),'/Global Mollweide Maps/Accessibility_to_UrbanAreas.tif'))
# Country id
country.id <- raster(paste0(getwd(),'/Global Mollweide Maps/MollweideCountryID_1.5km.tif'))

# Reprojecting 2050 map
urban.2050 <- projectRaster(urban.2050, urban.2010, method = 'bilinear')

# Creating data frame
df <- 
  data.frame(urban.2010 = getValues(urban.2010),
             urban.2050 = getValues(urban.2050),
             travel.time = getValues(travel.time),
             country.id = getValues(country.id),
             cell.id = 1:(country.id@ncols*country.id@nrows)) %>%
  mutate(count = 1)

# converting nas to 0s in urban.2050
df$urban.2050[!is.na(df$urban.2010) & is.na(df$urban.2050)] <- 0
# Dropping nas
df <- df[!is.na(df$urban.2010),]
df <- df[!is.na(df$urban.2050),]
df <- df[!is.na(df$country.id),]
df <- df[!is.na(df$travel.time),]

# Getting difference in land
df$urban.dif <- df$urban.2050 - df$urban.2010

# Urban 2050 cannot be less than urban 2010
df$urban.2050[df$urban.dif < 0] <- 
  df$urban.2010[df$urban.dif < 0]

# And recallculating to check
df$urban.dif <- df$urban.2050 - df$urban.2010

# (2) Summing difference by country -----
df.sum <-
  rowsum(df[,c('urban.2010','urban.2050','country.id','count')], group = df[,'country.id'], na.rm = TRUE) %>% # Summing by country
  mutate(five.year.change = (urban.2050 - urban.2010) / 8) %>% # Change per five year period
  mutate(country.id = country.id/count) %>% # updating country id values
  filter(five.year.change > 0) %>% # Only keeping countries with positive change
  filter(country.id > 0) # Negative values indicate the ocean...

# And merigng in targets
df.targets <-
  left_join(df,
            df.sum %>% dplyr::select(country.id, five.year.change))

# (3) Updating urban cells by year -----

# Ordering data frame first
df.targets <-
  df.targets[order(df.targets$country.id, df.targets$travel.time, -df.targets$urban.2010),]

# Next getting rid of rows where difference in urban area == 0
# Doing this to speed up processing
df.targets <-
  df.targets %>%
  filter(urban.dif > 0)

# Getting cumulative sum by country
df.targets <-
  df.targets %>%
  group_by(country.id) %>%
  mutate(urban.exp_cumsum = cumsum(urban.dif))
  

# And year for urban expansion
df.targets <-
  df.targets %>%
  mutate(year.exp = ifelse(urban.exp_cumsum <= five.year.change, 2015,
                           ifelse(urban.exp_cumsum <= five.year.change*2, 2020,
                                  ifelse(urban.exp_cumsum <= five.year.change*3, 2025,
                                         ifelse(urban.exp_cumsum <= five.year.change*4, 2030,
                                                ifelse(urban.exp_cumsum <= five.year.change*5, 2035,
                                                       ifelse(urban.exp_cumsum <= five.year.change*6, 2040,
                                                              ifelse(urban.exp_cumsum <= five.year.change*7, 2045, 2050))))))))

# And checking to make sure this has worked
# Looks alright!
df.targets.check <-
  df.targets %>%
  group_by(country.id, year.exp) %>%
  summarise(urban.dif = sum(urban.dif))

###
# Alrighty, now making rasters out of these again
# limiting columns for processing time
df.rasters <-
  df.targets %>% 
  as.data.frame(.) %>%
  dplyr::select(cell.id, urban.2010, urban.2050, year.exp)

# Rows to add
add.df <-
  data.frame(cell.id = 1:(country.id@nrows * country.id@ncols),
             urban.2010 = getValues(urban.2010),
             urban.2050 = getValues(urban.2010),
             year.exp = NA) %>%
  filter(!(cell.id %in% df.rasters$cell.id))

# Adding back in missing cell ids
df.rasters1 <-
  rbind(df.rasters,
        data.frame(cell.id = 1:(country.id@nrows * country.id@ncols),
                   urban.2010 = getValues(urban.2010),
                   urban.2050 = getValues(urban.2010),
                   year.exp = 2010) %>%
          filter(!(cell.id %in% df.rasters$cell.id)))

# (4) Making rasters ----

# Sorting df rasters to save time
df.rasters1 <-
  df.rasters1[order(df.rasters1$cell.id),]
# Creating directory for the urbanisation projections
if('Urbanization_Forecasts' %in% list.files(getwd())) {
  # Do nothing
} else {
  dir.create(paste0(getwd(),'/Urbanization_Forecasts/'))
}

# Looping through rasters
for(y in seq(2010, 2050, length.out = 9)) {
 if(y %in% 2010) {
   # Creating raster in new location
   writeRaster(urban.2010,
               paste0(getwd(),'/Urbanization_Forecasts/UrbanProjection_',y,'.tif'))
 } else {
   # Updating values
   tmp.df <-
     df.rasters1 %>%
     mutate(urban.extent = ifelse(year.exp <= y, urban.2050, urban.2010)) %>%
     mutate()
   # Making raster
   tmp.raster <-
     raster(matrix(tmp.df$urban.extent, nrow = country.id@nrows, ncol = country.id@ncols, byrow = TRUE))
   crs(tmp.raster) = crs(country.id)
   extent(tmp.raster) = extent(country.id)
   
   
 }
}

for(i in countries.loop) {
  # Getting country id
  tmp.df <- 
    df %>% 
    filter(country.id %in% i) %>%
    filter(urban.2010 > 0 | urban.2050 > 0) %>% # only need cells with some expansion
    filter(urban.dif > 0)

  # Ordering by travel time to city and current urban extent
  tmp.df <-
    tmp.df[order(tmp.df$travel.time, -tmp.df$urban.2010),]
  
  # Getting cumulative sum
  tmp.df <-
    tmp.df %>%
    mutate(cumsum_urban = cumsum(urban.dif))
  
  # Assigning
  }





