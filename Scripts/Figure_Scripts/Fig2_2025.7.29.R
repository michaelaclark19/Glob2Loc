#####
###
# Making figures to show changes in habitat area and population abundance by ecoregion

###
# Libraries
library(raster)
library(terra)
library(plyr)
library(dplyr)
library(data.table)
library(readr)
library(RColorBrewer)

###
# Setting working directory
setwd('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity')

###
# Getting habitat loss
file_list <-
  list.files(paste0(getwd(),'/Outputs/Aggregated_CSV_Files/'),
             pattern = 'Prop',
             full.names = TRUE) %>%
  .[grepl('prevalence',.)] %>%
  .[grepl('nopatch',.)] %>%
  .[!grepl('Climate_Suitable',.)]

# importing files
out_df <-
  do.call(rbind, lapply(file_list, read_csv)) %>%
  mutate(taxon = gsub('/.*','',binomial)) %>%
  mutate(species = gsub('.*/','',binomial))

# checking size
out_df_unique <-
  out_df %>%
  dplyr::select(taxon,species) %>%
  distinct()

# Species by ecoregion list
prop_biodiv_hotspot <-
  read_csv(paste0(getwd(),'/Ecoregions_Feb2023/species_by_ecoregion_esh_maps.csv'))

###
# Data management for function
hotspots_raster <- 
  raster(paste0(getwd(),'/Ecoregions_Feb2023/TNC_Ecoregions_Map.tif'))

biome_id_map <-
  raster(paste0(getwd(),'/Ecoregions_Feb2023/biome_id__mollweide.tif'))

# Making data table
hotspots_df <- 
  data.table(ecoregion_id = getValues(hotspots_raster),
             biome_id = getValues(biome_id_map)) %>%
  mutate(ecoregion_id = ifelse(biome_id %in% 99,NA,ecoregion_id)) %>%
  dplyr::select(-biome_id) %>%
  mutate(row_index = 1:nrow(.))

setkey(hotspots_df,ecoregion_id)

# Function to convert outputs to map
map_function <- 
  function(ii) {
    cat(ii)
    # Getting scenario data
    tmp_df <-
      out_df_2050 %>%
      filter(scenario %in% ii) %>%
      filter(ecoregion_id > 0)
    
    tmp_df <- data.table(tmp_df)
    
    # Merging outcomes by ecoregion to the data frame containing ids for each ecoregion
    
    # Checking whether this merges quickly
    tmp_df <- 
      tmp_df[hotspots_df, on = 'ecoregion_id', nomatch = NA] 
    setorder(tmp_df, row_index)
    
    # Converting to raster
    tmp_raster <-
      raster(matrix(tmp_df$tot_pop,
                    nrow = hotspots_raster@nrows,
                    byrow = TRUE),
             crs = crs(hotspots_raster))
    
    # Updating extent
    extent(tmp_raster) <- extent(hotspots_raster)
    names(tmp_raster) <- ii
    
    rm(tmp_df)
    gc()
    
    # And returning
    writeRaster(tmp_raster,
                paste0(getwd(),'/Files_For_Fig3/',ii,'.tif'),
                overwrite = TRUE)
    rm(tmp_raster)
  }


###
# Getting summary data for different combos of stressors
out_df_2050 <-
  rbind(out_df %>%
          filter(year %in% 2050) %>%
          left_join(.,
                    prop_biodiv_hotspot %>%
                      mutate(species = gsub('_[A-Z].*|_[0-9].*','',species)) %>%
                      mutate(binomial = paste0(taxa,'/',species)) %>%
                      dplyr::select(-c(species,taxa))) %>%
          dplyr::group_by(scenario,ecoregion_id) %>%
          dplyr::summarise(median_pop = median(log2(tot_pop), na.rm = TRUE),
                           median_area = median(log2(tot_area), na.rm = TRUE),
                           lower_pop = quantile(log2(tot_pop),.05,na.rm=TRUE),
                           median_patches = median(log2(num_patches), na.rm = TRUE),
                           num_species = n()) %>%
          filter(!is.na(ecoregion_id)) %>%
          mutate(median_area = 2^(median_area),
                 median_pop = 2^(median_pop),
                 median_patches = 2^(median_patches),
                 lower_pop = 2^(lower_pop)) %>%
          mutate(tot_pop = median_pop) %>%
          dplyr::select(scenario, ecoregion_id, tot_pop),
        out_df %>%
          filter(year %in% 2050,
                 grepl('.int.urb',scenario)) %>%
          left_join(.,
                    prop_biodiv_hotspot %>%
                      mutate(species = gsub('_[A-Z].*|_[0-9].*','',species)) %>%
                      mutate(binomial = paste0(taxa,'/',species)) %>%
                      dplyr::select(-c(species,taxa))) %>%
          dplyr::group_by(scenario,taxon,ecoregion_id) %>%
          dplyr::summarise(median_pop = median(log2(tot_pop), na.rm = TRUE),
                           median_area = median(log2(tot_area), na.rm = TRUE),
                           lower_pop = quantile(log2(tot_pop),.05,na.rm=TRUE),
                           median_patches = median(log2(num_patches), na.rm = TRUE),
                           num_species = n()) %>%
          filter(!is.na(ecoregion_id)) %>%
          mutate(median_area = 2^(median_area),
                 median_pop = 2^(median_pop),
                 median_patches = 2^(median_patches),
                 lower_pop = 2^(lower_pop)) %>%
          mutate(tot_pop = median_pop) %>%
          dplyr::select(scenario = taxon, ecoregion_id, tot_pop))

for(ii in unique(out_df_2050$scenario)) { map_function(ii)}

out_raster_list <- lapply(list.files(path = paste0(getwd(),'/Files_For_Fig3/'), full.names = TRUE) %>% .[!grepl('lower',.)],rast)

plot_order <-
  c('Amphibians','Birds','Mammals','Reptiles',
    '.int.urb',
    'current','esh.exp.tif','esh.urb','esh.int.tif')

title_order <-
  c('Amphibians','Birds','Mammals','Reptiles',
    'All\nStressors',
    'Climate\nChange','Agricultural\nExpansion','Urban\nExpansion','Agricultural\nIntensification')


# matrix_layout <-
#   matrix(c(0,2,3,0,
#            1,5,5,4,
#            6,5,5,9,
#            0,7,8,0),
#          byrow = TRUE, nrow = 4)

###
# Getting country borders
borders <- vect(paste0(getwd(),'/TM_WORLD_BORDERS-0.3/TM_WORLD_BORDERS-0.3.shp'))
borders <- project(borders,rast(out_raster_list[1]))

# matrix_layout <-
#   matrix(c(2,1,1,1,1,6,
#            3,1,1,1,1,7,
#            4,1,1,1,1,8,
#            5,1,1,1,1,9),
#          byrow = TRUE, nrow = 4)


# plot_order <-
#   c('.int.urb',
#     'Amphibians','Birds','Mammals','Reptiles',
#     'current','esh.exp.tif','esh.urb','esh.int.tif')
# 
# title_order <-
#   c('All\nStressors',
#     'Amphibians','Birds','Mammals','Reptiles',
#     'Climate\nChange','Agricultural\nExpansion','Agricultural\nIntensification', 'Urban\nExpansion')

#####
###
# Creating a mask for the ocean
# Example raster (replace with your actual raster)
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

# Plot
# plot(r_oval, col = 'light blue', colNA = "white")

#####
###
# Getting this to work with base R
for(xx in c(1:9)) {
# for(xx in 1) {
  
  if(xx %in% 5) {
    pdf(paste0(getwd(),'/Figures and Tables/Figures/Fig1_Results_',title_order[xx],'_',Sys.Date(),'.pdf'),
        width = 15 / 2.54,
        height = 3.5)
    par(mar = c(0,0,0,0))
  } else {
    pdf(paste0(getwd(),'/Figures and Tables/Figures/Fig1_Results_',title_order[xx],'_',Sys.Date(),'.pdf'),
        width = 10 / 2.54,
        height = 2.3333)
    par(mar = c(0,0,0,0))
  }
  tmp <- list.files(path = paste0(getwd(),'/Files_For_Fig2/'),
                    pattern = plot_order[xx],
                    full.names = TRUE) %>%
    .[!grepl('_lower',.)] %>%
    rast()
  
  # tmp <- mask(tmp, mask_raster)
  
  # Getting range for plot
  # And updating range for plot
  tmp_range_df <- 
    out_df_2050 %>%
    filter(scenario %in% names(tmp))
  
  # Getting range for plot
  tmp_range <- c(0,2) 
  if(min(tmp_range_df$tot_pop,na.rm=TRUE) < 0.7) {tmp_range[1] <- .7}
  if(max(tmp_range_df$tot_pop,na.rm=TRUE) > 1.1) {tmp_range[2] <- 1.1}
  if(tmp_range[1] %in% 0) {tmp_range[1] <- quantile(tmp_range_df$tot_pop,c(0.01), na.rm = TRUE)}
  if(tmp_range[2] %in% 2) {tmp_range[2] <- quantile(tmp_range_df$tot_pop,c(0.99), na.rm = TRUE)}
  if(tmp_range[2] < 1){tmp_range[2]<-1}
  
  # Adjusting for plotting
  tmp_range[1] <- floor(min(tmp_range)*100)/100
  tmp_range[2] <- ceiling(max(tmp_range)*100)/100
  # tmp_range[2] <- 1
  if(tmp_range[2] %in% 1.11) {tmp_range[2] <- 1.1}
  if(tmp_range[2] > 1) {tmp_range[2] <- 1}
  
  # Upper and lower limits to the range
  tmp[tmp<min(tmp_range)] <- min(tmp_range)+1e-5
  # tmp[tmp>max(tmp_range,1)] <- max(tmp_range,1)
  
  ###
  # Making palette - centered on 1, but reds < 1 and blues > 1
  # Have two options coded below - just need to uncomment the one you would like to use
  
  # Getting reds and blues
  
  red_palette <- colorRampPalette(c('#dddddd','#f2bf91','#f28379','#d95252'))
  blue_palette <- colorRampPalette(rev(c('#3581b9','#99e9fc','#cbf1f7','#dddddd')))
  
  # red_palette <- colorRampPalette(brewer.pal(n = 11, 'RdBu')[6:1])
  # blue_palette <- colorRampPalette(brewer.pal(n = 11, 'RdBu')[6:11])
  
  # Getting relative lengths of the reds and blues
  red_length <- round((1 - min(tmp_range)) / (max(tmp_range) - min(tmp_range)) * 100, digits = 0)
  blue_length <- 100 - red_length
  
  if(red_length >= 100) {
    red_length <- 100
    blue_length <- 0
  }
  
  # Making palette based on this skew
  my_palette_plot <-
    c(rev(red_palette(red_length)),
      blue_palette(blue_length))
  
  my_palette_plot <- rev(red_palette(99))
  
  ###
  # Creating plot
  plot(r_oval, 
       # col = 'light blue',
       col = brewer.pal(n = 9, name = 'Blues')[2],
       # col = '#eeeeee',
       box = FALSE,
       axes = FALSE,
       legend = FALSE,
       mar = c(0,0,0,0))
  plot(tmp,
       add = TRUE,
       box = FALSE,
       axes = FALSE,
       legend = FALSE,
       range = tmp_range,
       col = my_palette_plot,
       mar = c(0,0,0,0))
  
  if(max(minmax(tmp)) > 1) {
    tmp_gains <- tmp
    tmp_gains[tmp_gains <= 1] <- NA
    
    plot(tmp_gains,
         add = TRUE,
         box = FALSE,
         axes = FALSE,
         legend = FALSE,
         mar = c(0,0,0,0),
         col = brewer.pal(n = 9, name = 'Blues')[7])
         # col = 'pinxk')
  }
 
  
  plot(borders,
       add = TRUE,
       lwd = .5)
  
  map_boundaries <- boundaries(x = r_oval)
  map_boundaries[map_boundaries %in% 0] <- NA
  
  plot(map_boundaries, add = TRUE, col = 'black')
  
  
  ### Adding legend - location changes by plot
  if(xx %in% 5) {
    ###
    # Adding legend
    # Define rectangle bounds
    xleft <- .1
    xright <- .12
    ybottom <- .35
    ytop <- .55
    
    # Convert to user coordinates
    xleft <- grconvertX(xleft, from = "ndc", to = "user")
    xrightlabel <- grconvertX(xright+.025, from = "ndc", to = "user")
    xright <- grconvertX(xright, from = "ndc", to = "user")
    ybottom <- grconvertY(ybottom, from = "ndc", to = "user")
    ytitlelabel <- grconvertY(ytop+.1, from = "ndc", to = "user")
    ytop_pos <- grconvertY(ytop + .05, from = "ndc", to = "user")
    ybottom_pos <- grconvertY(ytop + .00, from = "ndc", to = "user")
    ytop <- grconvertY(ytop, from = "ndc", to = "user")
    
    n_colors <- length(my_palette_plot)
    
    # Compute width of each sub-rectangle
    y_steps <- seq(ybottom, ytop, length.out = n_colors + 1)
    
    # Draw gradient using thin vertical rectangles
    for (i in 1:n_colors) {
      rect(xleft = xleft, ybottom = y_steps[i],
           xright = xright, ytop = y_steps[i+1],
           col = my_palette_plot[i], border = NA,
           xpd = TRUE)
    }
    
    rect(xleft = xleft,
         xright = xright,
         ybottom = ybottom_pos,
         ytop = ytop_pos,
         col = brewer.pal(n = 9, name = 'Blues')[7],
         xpd = TRUE,
         border = NA)
    
    # Optional: draw border around the whole area
    rect(xleft, ybottom, xright, ytop, border = "black", xpd = TRUE)
    rect(xleft, ybottom_pos, xright, ytop_pos, border = "black", xpd = TRUE)
    
    
    ###
    # Adding legend text
    text(x = xrightlabel, 
         y = ytop, 
         label = round(max(tmp_range),digits = 2),
         xpd = TRUE,
         adj = c(0,.5),
         cex = .6)
    
    text(x = xrightlabel, 
         y = mean(c(ytop,ybottom)), 
         label = round(mean(tmp_range),digits = 2),
         xpd = TRUE,
         adj = c(0,.5),
         cex = .6)
    
    text(x = xrightlabel, 
         y = ybottom, 
         label = round(min(tmp_range),digits = 2),
         xpd = TRUE,
         adj = c(0,.5),
         cex = .6)
    
    text(x = xrightlabel, 
         y = mean(c(ytop_pos)), 
         label = '> 1',
         xpd = TRUE,
         adj = c(0,.5),
         cex = .6)
    
    
    text(x = mean(xleft,xright),
         y = ytitlelabel,
         label = 'Change in\nPop. Abundance',
         adj = c(.5,.5),
         cex = .6,
         xpd = TRUE)
    
    text(x = 0,
         y = ext(tmp)[4] - (ext(tmp)[4] * 1.78),
         label = title_order[xx],
         adj = c(.5,0),
         xpd = TRUE)
  } else {
    ###
    # Adding legend
    # Define rectangle bounds
    if(max(minmax(tmp)) >= 1) {
      xleft <- .1
      xright <- .13
      ybottom <- .25
      ytop <- .45
      
      # Convert to user coordinates
      xleft <- grconvertX(xleft, from = "ndc", to = "user")
      xrightlabel <- grconvertX(xright+.025, from = "ndc", to = "user")
      xright <- grconvertX(xright, from = "ndc", to = "user")
      ybottom <- grconvertY(ybottom, from = "ndc", to = "user")
      ytitlelabel <- grconvertY(ytop+.1, from = "ndc", to = "user")
      ytop_pos <- grconvertY(ytop + .05, from = "ndc", to = "user")
      ybottom_pos <- grconvertY(ytop + .00, from = "ndc", to = "user")
      ytop <- grconvertY(ytop, from = "ndc", to = "user")
      
      
      n_colors <- length(my_palette_plot)
      
      # Compute width of each sub-rectangle
      # Compute width of each sub-rectangle
      y_steps <- seq(ybottom, ytop, length.out = n_colors + 1)
      
      # Draw gradient using thin vertical rectangles
      for (i in 1:n_colors) {
        rect(xleft = xleft, ybottom = y_steps[i],
             xright = xright, ytop = y_steps[i+1],
             col = my_palette_plot[i], border = NA,
             xpd = TRUE)
      }
      
      rect(xleft = xleft,
           xright = xright,
           ybottom = ybottom_pos,
           ytop = ytop_pos,
           col = brewer.pal(n = 9, name = 'Blues')[7],
           xpd = TRUE,
           border = NA)
      
      # Optional: draw border around the whole area
      rect(xleft, ybottom, xright, ytop, border = "black", xpd = TRUE)
      rect(xleft, ybottom_pos, xright, ytop_pos, border = "black", xpd = TRUE)
      
      
      ###
      # Adding legend text
      text(x = xrightlabel, 
           y = ytop, 
           label = round(max(tmp_range),digits = 2),
           xpd = TRUE,
           adj = c(0,.5),
           cex = .6)
      
      text(x =xrightlabel, 
           y = mean(c(ybottom,ytop)), 
           label = round(mean(tmp_range),digits = 2),
           xpd = TRUE,
           adj = c(0,.5),
           cex = .6)
      
      text(x = xrightlabel, 
           y = mean(c(ybottom)), 
           label = round(min(tmp_range),digits = 3),
           xpd = TRUE,
           adj = c(0,.5),
           cex = .6)
      
      text(x = xrightlabel,
           y = ytop_pos,
           label = '>1',
           xpd = TRUE,
           adj = c(.5,1),
           cex = .6)
      
      
      text(x = mean(c(xleft,xright)),
           y = ytitlelabel,
           label = 'Change in\nPop. Abundance',
           adj = c(.5,0),
           cex = .6,
           xpd = TRUE)
      
      text(x = 0,
           y = ext(tmp)[4] - (ext(tmp)[4] * 1.78),
           label = title_order[xx],
           adj = c(.5,0),
           xpd = TRUE)
      
    } else {
      xleft <- .08
      xright <- .11
      ybottom <- .28
      ytop <- .53
      
      # Convert to user coordinates
      xleft <- grconvertX(xleft, from = "ndc", to = "user")
      xrightlabel <- grconvertX(xright+.025, from = "ndc", to = "user")
      xright <- grconvertX(xright, from = "ndc", to = "user")
      ybottom_label <- grconvertY(ybottom-.025, from = "ndc", to = "user")
      ybottom <- grconvertY(ybottom, from = "ndc", to = "user")
      ytitlelabel <- grconvertY(ytop+.025, from = "ndc", to = "user")
      ytop <- grconvertY(ytop, from = "ndc", to = "user")
      
      n_colors <- length(my_palette_plot)
      
      # Compute width of each sub-rectangle
      n_colors <- length(my_palette_plot)
      
      # Compute width of each sub-rectangle
      # Compute width of each sub-rectangle
      y_steps <- seq(ybottom, ytop, length.out = n_colors + 1)
      
      # Draw gradient using thin vertical rectangles
      for (i in 1:n_colors) {
        rect(xleft = xleft, ybottom = y_steps[i],
             xright = xright, ytop = y_steps[i+1],
             col = my_palette_plot[i], border = NA,
             xpd = TRUE)
      }
      
      # Optional: draw border around the whole area
      rect(xleft, ybottom, xright, ytop, border = "black", xpd = TRUE)
      
      
      ###
      # Adding legend text
      text(x = xrightlabel, 
           y = ytop, 
           label = round(max(tmp_range),digits = 2),
           xpd = TRUE,
           adj = c(0,.5),
           cex = .6)
      
      text(x =xrightlabel, 
           y = mean(c(ybottom,ytop)), 
           label = round(mean(tmp_range),digits = 2),
           xpd = TRUE,
           adj = c(0,.5),
           cex = .6)
      
      text(x = xrightlabel, 
           y = mean(c(ybottom)), 
           label = round(min(tmp_range),digits = 3),
           xpd = TRUE,
           adj = c(0,.5),
           cex = .6)
      
      
      text(x = mean(c(xleft,xright)),
           y = ytitlelabel,
           label = 'Change in\nPop. Abundance',
           adj = c(.5,0),
           cex = .6,
           xpd = TRUE)
      
      text(x = 0,
           y = ext(tmp)[4] - (ext(tmp)[4] * 1.78),
           label = title_order[xx],
           adj = c(.5,0),
           xpd = TRUE)
    }
    
  }
  dev.off()
}

