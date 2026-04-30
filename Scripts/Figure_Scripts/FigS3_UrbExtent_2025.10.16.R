###
# Plotting change in urbanisation from 2020 to 2050

# libraries
library(raster)
library(RColorBrewer)
library(viridis)
library(terra)

# loading rasters
urb_2020 <- rast('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Urbanization_Forecasts/Urban Extent 2020.tif')
urb_2050 <- rast('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Urbanization_Forecasts/Urban Extent 2050.tif')

# Removing non land cells
crop_2020 <- rast('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Land_Forecasting_Outputs/crop_mean2015_2020.tif')

urb_2020[is.na(crop_2020)] <- NA
urb_2050[is.na(crop_2020)] <- NA
urb_2020[is.na(urb_2020) & !is.na(crop_2020)] <- 0
urb_2050[is.na(urb_2050) & !is.na(crop_2020)] <- 0

# getting change
urb_delta <- urb_2050 - urb_2020

# getting map of countries
borders <- vect('/Users/michael/Downloads/TM_WORLD_BORDERS-0.3/TM_WORLD_BORDERS-0.3.shp')
borders <- project(borders,crop_2020)

# Colour palette
red_palette <- colorRampPalette(c('#ffffff','#f2bf91','#f28379','#d95252'))(100)
blue_palette <- colorRampPalette(rev(c('#3581b9','#99e9fc','#cbf1f7','#ffffff')))(100)

my_palette <- 
  c(rev(blue_palette),
    '#ffffff',
    red_palette)

###
# Making function for plot
plot_function <-
  function(raster_input,
           palette,
           limits,
           title) {
    # cropland
    # plot(r_oval, 
    #      # col = 'light blue',
    #      col = brewer.pal(n = 9, name = 'Blues')[2],
    #      # col = '#eeeeee',
    #      box = FALSE,
    #      axes = FALSE,
    #      legend = FALSE,
    #      mar = c(0,0,0,0))
    
    plot(raster_input,
         plg = list(loc = 'left', title = 'Proportion\nof Cell', cex = 1, title.cex = 1, legend.size = .1),
         colNA = brewer.pal(n = 9, name = 'Blues')[2],
         col = palette,
         # col = my_palette,
         main = title,
         cex.main = 1,
         axes = FALSE,
         # add = TRUE,
         range = limits)
    
    plot(borders,
         add = TRUE,
         border = '#dddddd',
         lwd = .5)
  }


###
# plotting
# dev.off()
pdf(paste0('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Figures and Tables/Figures/FigS3_UrbExpansion_',Sys.Date(),'.pdf'),
    width = 12, height = 4)
par(mfrow = c(1,3))

# cropland
plot_function(raster_input = urb_2020,
              palette = viridis(99),
              limits = c(0,1),
              title = 'Urban Extent in 2020')

plot_function(raster_input = urb_2050,
              palette = viridis(99),
              limits = c(0,1),
              title = 'Urban Extent in 2050')

plot_function(raster_input = urb_delta,
              palette = my_palette[101:201],
              limits = c(0,1),
              title = 'Change in Urban Extent 2020 - 2050')

dev.off()