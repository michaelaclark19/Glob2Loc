###
# Plotting change in urbanisation from 2020 to 2050

# libraries
library(raster)
library(RColorBrewer)
library(viridis)
library(terra)

# loading rasters
gdd_binary <- rast('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/GDD_Binary_2025.tif')

# 0 corresponds to not suitable in 2020 or 2050 (by growing length)
# 1 corresponds to suitable under both 5 and 10c
# 2 corresponds to suitable under 5 not 10
# 3 corresponds to suitable under 10 not 5

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
dev.off()
pdf(paste0('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Figures and Tables/Figures/FigS3_AgExpansion_',Sys.Date(),'.pdf'),
    width = 10, height = 8)
par(mfrow = c(3,2))

# cropland
plot_function(raster_input = crop_2020,
              palette = viridis(99),
              limits = c(0,1),
              title = 'Cropland Extent in 2020')

plot_function(raster_input = crop_delta,
              palette = my_palette,
              limits = c(-1,1),
              title = 'Change in Cropland Extent 2020 - 2050')

plot_function(raster_input = past_2020,
              palette = viridis(99),
              limits = c(0,1),
              title = 'Pastureland Extent in 2020')

plot_function(raster_input = past_delta,
              palette = my_palette,
              limits = c(-1,1),
              title = 'Change in Pastureland Extent 2020 - 2050')

plot_function(raster_input = ag_2020,
              palette = viridis(99),
              limits = c(0,1),
              title = 'Agricultural Extent in 2020')

plot_function(raster_input = ag_delta,
              palette = my_palette,
              limits = c(-1,1),
              title = 'Change in Agricultural Extent 2020 - 2050')

dev.off()