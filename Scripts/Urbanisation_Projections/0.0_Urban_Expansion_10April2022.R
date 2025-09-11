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
library(rgdal)

# (0) setting working directory
setwd("/data/pubh-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity")

# Raw urban projections
urban.projections <- list.files(paste0(getwd(),'/Raw Urban 2050 Projections'), full.names = TRUE)

# Function to convert the raw urban projections into urban forecasts

urban.project.function <-
  function(raster.file) {
    # (1) Importing rasters and raster management ----
    # 2010 raster from ESA
    # 2050 Raster from huang et al 2019
    urban.2010 <- raster(paste0(getwd(),'/ESA LandCov Maps/Urban2010_Corrected_GlobCov.tif'))
    urban.2050 <- raster(raster.file)
    # 2015 raster from globmod GHS_SMOD
    urban.2015.globmod <- readOGR(paste0(getwd(),"/Global Mollweide Maps/CurrentUrbanGlobMod/GHS_SMOD_POP2015_GLOBE_R2019A_54009_1K_labelHDC_V2_0.gpkg"))
    urban.2015.raster <- rasterize(urban.2015.globmod, urban.2050, 'BU_2015')
    # Projecting all to 5x5km mollweide for consistency
    urban.2010.tmp <-
      projectRaster(urban.2010, urban.2050, method = 'bilinear')
    
    # (2) Comparing 2010 esa and 2015 from globmod ----
    # Converting to numeric
    urban.2015.raster.keep <- urban.2015.raster
    urban.2015.raster[is.na(urban.2015.raster)] <- 0
    urban.2015.raster[urban.2015.raster>0]<-1
    
    # Total urban extent
    # These are close
    sum(getValues(urban.2010.tmp), na.rm = TRUE)
    sum(getValues(urban.2015.raster))
    
    # Difference between rasters
    # THese are decently close
    # BUt definite differences between the ssp map is binary, whereas the esa map is continuous
    urban.delta.ssp.esa <- urban.2015.raster - urban.2010.tmp
    plot(urban.delta.ssp.esa)
    
    # (3) Getting cells where urban extent expands ----
    urban.2050[is.na(urban.2050)] <- 0
    urban.delta.esa <- urban.2050-urban.2010.tmp
    urban.delta.ssp <- urban.2050 - urban.2015.raster
    
    # Maps are different enough that we are:
    # (a) getting delta between the ssp maps for consisntency, and summing by country
    # (b) applying this delta to 2010 esa map
    # This means urban increase is the same, but exact location of future urban land is not
    
    # Don't want negative urban expansion
    # This results from mismatched maps
    urban.delta.ssp[urban.delta.ssp < 0] <- 0
    
    # (4) Urban expansion by country
    # country id raster
    country.id <- raster(paste0(getwd(),"/Global Mollweide Maps/MollweideCountryID_1.5km.tif"))
    
    # projecting to 5x5km
    country.id.rough <- projectRaster(country.id, urban.2050, method = 'ngb')
    
    # Getting 5 year expansion by country
    exp.country.df <-
      data.frame(country = getValues(country.id.rough),
                 urban.2015 = getValues(urban.2015.raster),
                 urban.2050 = getValues(urban.2050)) %>%
      mutate(count = 1) %>%
      filter(!is.na(country)) %>%
      rowsum(.,group = .$country) %>%
      mutate(country = country / count) %>%
      mutate(exp_5_years = (urban.2050 - urban.2015) / 8)
    
    # (5) Reprojecting urban 2050 to 1.5 x 1.5km
    urban.2050.fine <- projectRaster(urban.2050, urban.2010, method = 'bilinear')
    urban.delta.fine <- projectRaster(urban.delta.ssp, urban.2010, method = 'bilinear')
    # plot(urban.2050.fine)
    
    # (6) Creating data frame ----
    # This will be used to select cells that experience an expansion in urban extent
    
    # Importing rasters on travel to cities
    travel.time <- raster(paste0(getwd(),'/Global Mollweide Maps/Accessibility_to_UrbanAreas.tif'))
    
    # Creating data frame
    df.fine <-
      data.frame(country = getValues(country.id),
                 urban.2010 = getValues(urban.2010),
                 urban.2050 = getValues(urban.2050.fine),
                 urban.delta = getValues(urban.delta.fine),
                 travel.time = getValues(travel.time),
                 cell.id = 1:(urban.2010@nrows * urban.2010@ncols)) %>%
      mutate(urban.delta = ifelse(!is.na(urban.delta) & urban.delta < 0,0,urban.delta)) # No urban abandonment
    
    # Urban delta + urban 2010 cannot be bigger than 1
    df.fine$check <- df.fine$urban.2010 + df.fine$urban.delta
    df.fine$urban.delta[df.fine$check > 1 & !is.na(df.fine$check)] <-
      1 - df.fine$urban.2010[df.fine$check > 1 & !is.na(df.fine$check)]
    
    # # Filter out oceans and cells with no urban area
    # df.fine <-
    #   df.fine %>%
    #   filter(!is.na(country)) %>%
    #   filter(country > 0) %>%
    #   filter(!is.na(urban.2010)) %>%
    #   filter(!is.na(urban.2050)) %>%
    #   # filter(urban.2010 > 0) %>%
    #   filter(urban.2050 > 0 | urban.2010 > 0) %>%
    #   filter(urban.delta > 0)
    
    # (7) And expanding ----
    # Assuming that each cell expands 1/8 of their total expansion in each 5 year interval
    df.fine <-
      df.fine %>%
      mutate(urban.2015 = urban.2010 + urban.delta * (1/8)) %>%
      mutate(urban.2020 = urban.2010 + urban.delta * (2/8)) %>%
      mutate(urban.2025 = urban.2010 + urban.delta * (3/8)) %>%
      mutate(urban.2030 = urban.2010 + urban.delta * (4/8)) %>%
      mutate(urban.2035 = urban.2010 + urban.delta * (5/8)) %>%
      mutate(urban.2040 = urban.2010 + urban.delta * (6/8)) %>%
      mutate(urban.2045 = urban.2010 + urban.delta * (7/8)) %>%
      mutate(urban.2050 = urban.2010 + urban.delta * (8/8))
    
    # (8) Making rasters ----
    
    # Creating directory to save
    if('Urbanization_Forecasts' %in% list.files(getwd())) {
      # Do nothing
    } else {
      # Create folder
      dir.create(paste0(getwd(),'/Urbanization_Forecasts'))
    }
    
    # Creating directory for the SSPs
    ssp.scenario <- str_extract(raster.file,'ssp[0-9].tif') %>% gsub(".tif","",.)
    dir.create(paste0(getwd(),'/Urbanization_Forecasts/',ssp.scenario))
    # And looping through to make rasters
    for(i in seq(2010, 2050, by = 5)) {
      # Extracting values
      tmp.vector = df.fine[,grepl(i,names(df.fine))]
      # Making raster
      tmp.raster = raster(matrix(tmp.vector, byrow = TRUE, nrow = urban.2010@nrows))
      crs(tmp.raster) = crs(urban.2010)
      extent(tmp.raster) = extent(urban.2010)
      # Saving rasters
      writeRaster(tmp.raster,
                  paste0(getwd(),'/Urbanization_Forecasts/',ssp.scenario,'Urban_Extent_', i,'.tif'),
                  overwrite = TRUE)
    } # End loop through years
  } # End function


# And now running the forecasts
lapply(urban.projections, urban.project.function)

# END