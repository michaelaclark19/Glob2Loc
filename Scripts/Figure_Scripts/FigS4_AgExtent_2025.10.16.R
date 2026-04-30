###
# Plotting change in urbanisation from 2020 to 2050

# libraries
library(raster)
library(RColorBrewer)
library(viridis)
library(terra)

# loading rasters
crop_2020 <- rast('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Land_Forecasting_Outputs/crop_mean2015_2020.tif')
crop_2050 <- rast('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Land_Forecasting_Outputs/crop_mean2045_2050.tif')

past_2020 <- rast('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Land_Forecasting_Outputs/pasture_mean2015_2020.tif')
past_2050 <- rast('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Land_Forecasting_Outputs/pasture_mean2045_2050.tif')

ag_2020 <- crop_2020 + past_2020
ag_2050 <- crop_2050 + past_2050

# getting change
crop_delta <- crop_2050 - crop_2020
past_delta <- past_2050 - past_2020
ag_delta <- ag_2050 - ag_2020

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
# Creating a mask for the ocean
r <- rast(nrows=180, ncols=360, xmin=-180, xmax=180, ymin=-90, ymax=90, vals=runif(180*360))
crs(r) <- "EPSG:4326"

# Project raster to Mollweide (you can skip this if yours is already projected)
r_moll <- project(r, "+proj=moll +datum=WGS84")

# Create graticule points to approximate the globe
lons <- seq(-180, 180, by=1)
lats <- seq(-90, 90, by=1)
grd <- expand.grid(lon = lons, lat = lats)
pts <- vect(grd, crs = "EPSG:4326")  # Vector points in lat/lon

# Project graticule points to Mollweide
pts_moll <- project(pts, crs(r_moll))

# Convert to a convex hull polygon that covers the oval shape
hull <- convHull(pts_moll)

# Rasterize the convex hull to match the raster
mask_raster <- rasterize(hull, r_moll, field=1, background=NA)

# Mask the raster to keep only the oval
r_oval <- mask(r_moll, mask_raster)

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