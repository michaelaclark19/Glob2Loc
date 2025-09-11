#####
# Plotting location of populations map
#####


###
# Function to identify overlapping pops
pops_2020 <- raster('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Species_Analyses_New/hab_availability_Phasmahyla_guttata_prevalence_2020_esh.future.extent.esh.exp.int.urb.tif')
pops_2050 <- raster('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Species_Analyses_New/hab_availability_Phasmahyla_guttata_prevalence_2050_esh.future.extent.esh.exp.int.urb.tif')


pops_2020 <- 
  raster(outcome_rasters_prevalence %>%
           .[grepl(ss,.)] %>%
           .[grepl('2020',.)] %>%
           .[grepl('hab_avail.*.int.urb.tif',.)])
# Getting 2050 map to extend raster
pops_2050 <- 
  raster(outcome_rasters_prevalence %>%
           .[grepl(ss,.)] %>%
           .[grepl('2050',.)] %>%
           .[grepl('hab_avail.*.int.urb.tif',.)])

pops_2020 <- raster::extend(pops_2020,pops_2050)
pops_2050 <- raster::extend(pops_2050,pops_2020)

# converting to matrix
vals_pops_2020 <- getValues(pops_2020)
trial_ncols = pops_2020@ncols

# now getting migration distance between patches
py_run_string("d = r.dispersal_distance") # Dispersal distance
py_run_string("dp = d / 1000 / 1.5") # Numbers of cells
py_run_string("im = r.vals_pops_2020") # Creating rasters with ones and 0s
py_run_string("im = np.reshape(im,(-1,r.trial_ncols))") # Creating raster
py_run_string("im_migrate = copy.copy(im)") # Copying
py_run_string("im_migrate[np.isnan(im_migrate)] = 0") # NAs to 0s so convolve works
py_run_string("im_migrate[im_migrate>=.2] = 1") # Adjusting for migration ability
py_run_string("im_migrate[im_migrate<.2] = 0") # Assumes cell needs >=20% suitable habitat for species for a species to migrate through it

# Assuming there is no migration between habitat patches
py_run_string("nim = im_migrate != im_migrate.min()")
py_run_string("bim_no_migrate = morphology.binary_dilation(nim, morphology.disk(0, dtype=bool))")
py_run_string("labim_no_migrate = measure.label(bim_no_migrate)")
py_run_string("labim_no_migrate[~nim] = 0")

# Assuming migration between patches
# With an if statement to save processing time
# E.g. if dispersal less than a cell, then can't disperse
py_run_string("if(r.dispersal_distance < 1500): 
                      labim = copy.copy(labim_no_migrate)")
py_run_string("if(r.dispersal_distance >= 1500):
                      bim = morphology.binary_dilation(nim, morphology.disk(dp, dtype=bool))
                      labim = measure.label(bim)
                      labim[~nim] = 0")

# returning to R for pops in 2020
pops_2020_plot = raster(py$labim, crs = crs(pops_2020))
extent(pops_2020_plot) = extent(pops_2020)
# removing non habitat areas
pops_2020_plot[is.na(pops_2020)] <- NA
# getting background
pops_2020_background <- !is.na(pops_2020)
pops_2020_background[pops_2020_background %in% 0] <- NA

# returning to R for pops in 2050
vals_pops_2050 <- getValues(pops_2050)
trial_ncols = pops_2050@ncols

# now getting migration distance between patches
py_run_string("d = r.dispersal_distance") # Dispersal distance
py_run_string("dp = d / 1000 / 1.5") # Numbers of cells
py_run_string("im = r.vals_pops_2050") # Creating rasters with ones and 0s
py_run_string("im = np.reshape(im,(-1,r.trial_ncols))") # Creating raster
py_run_string("im_migrate = copy.copy(im)") # Copying
py_run_string("im_migrate[np.isnan(im_migrate)] = 0") # NAs to 0s so convolve works
py_run_string("im_migrate[im_migrate>=.2] = 1") # Adjusting for migration ability
py_run_string("im_migrate[im_migrate<.2] = 0") # Assumes cell needs >=20% suitable habitat for species for a species to migrate through it

# Assuming there is no migration between habitat patches
py_run_string("nim = im_migrate != im_migrate.min()")
py_run_string("bim_no_migrate = morphology.binary_dilation(nim, morphology.disk(0, dtype=bool))")
py_run_string("labim_no_migrate = measure.label(bim_no_migrate)")
py_run_string("labim_no_migrate[~nim] = 0")

# Assuming migration between patches
# With an if statement to save processing time
# E.g. if dispersal less than a cell, then can't disperse
py_run_string("if(r.dispersal_distance < 1500): 
                      labim = copy.copy(labim_no_migrate)")
py_run_string("if(r.dispersal_distance >= 1500):
                      bim = morphology.binary_dilation(nim, morphology.disk(dp, dtype=bool))
                      labim = measure.label(bim)
                      labim[~nim] = 0")

# and returning to R
pops_2050_plot = raster(py$labim, crs = crs(pops_2050))
extent(pops_2050_plot) = extent(pops_2050)
# removing non habitat areas
pops_2050_plot[is.na(pops_2050)] <- NA
pops_2050_background <- !is.na(pops_2050)
pops_2050_background[pops_2050_background %in% 0] <- NA
#
pop_df <-
  data.frame(pops_2020 = getValues(pops_2020_plot),
             pops_2050 = getValues(pops_2050_plot)) %>%
  distinct() %>%
  arrange(pops_2020) %>%
  mutate(pops_2020_check = !is.na(pops_2020),
         pops_2050_check = !is.na(pops_2050)) 

pops_only_2020_df <-
  pop_df %>%
  dplyr::group_by(pops_2020) %>%
  dplyr::summarise(check_2020 = sum(pops_2050_check)) %>%
  filter(check_2020 %in% 0) 

pops_only_2050_df <-
  pop_df %>%
  dplyr::group_by(pops_2050) %>%
  dplyr::summarise(check_2050 = sum(pops_2020_check)) %>%
  filter(check_2050 %in% 0) 


# Updating ID numbers
# But only if this needs to be done

# Updating ids of populations only in 2020
if(nrow(pops_only_2020_df) >= 1) {
  pops_only_2020_df <- 
    pops_only_2020_df %>%
    mutate(id_num_2020 = -1:-nrow(.))
  for(nn in 1:nrow(pops_only_2020_df)) {
    pops_2020_plot[pops_2020_plot %in% pops_only_2020_df$pops_2020[nn]] <- pops_only_2020$id_num_2020[nn]
  }
}

# Updating ids of populations only in 2050
if(nrow(pops_only_2050_df) >= 1) {
  pops_only_2050_df <- 
    pops_only_2050_df %>%
    mutate(id_num_2050 = -1:-nrow(.))
  for(nn in 1:nrow(pops_only_2050_df)) {
    pops_2050_plot[pops_2050_plot %in% pops_only_2050_df$pops_2050[nn]] <- pops_only_2050$id_num_2050[nn]
  }
}

# Updating ids of populations in both 2020 and 2050
if(nrow(pop_df) >= 1) {
  pops_2050_both_out <-
    pops_2050_plot
  pops_2050_both_out[!is.na(pops_2050_both_out)] <- NA
  
  
  pop_df_update <-
    pop_df %>%
    filter(!is.na(pops_2020)) %>%
    filter(!is.na(pops_2050)) %>%
    dplyr::group_by(pops_2050) %>%
    mutate(pops_2050_update = max(pops_2020))
  
  for(nn in unique(pop_df_update$pops_2050_update)) {
    # getting pop ids to update
    ids_update <- unique(pop_df_update$pops_2050[pop_df_update$pops_2050_update %in% nn])
    pops_2050_both_out[pops_2050_plot %in% ids_update] <- nn
  }
}



# dev.off()
par(mar = c(0,0,0,0))
layout.matrix <- matrix(c(1,1,1,1,
                          2,2,3,3,
                          2,2,3,3,
                          2,2,3,3,
                          4,4,4,4),
                        byrow = TRUE,nrow = 5)
layout(mat = layout.matrix)
layout.show(4)

plot.new()
text(x=.5,y=.5,'Estimated Location of Populations',font = 2,cex = 2)

# 
if(nrow(pops_only_2020_df) >= 1) {
  pops_2020_plot_both <- pops_2020_plot
  pops_2020_plot_both[pops_2020_plot_both%in%pops_only_2020$id_num_2020] <- NA
  
  pops_2020_plot_only <- pops_2020_plot
  pops_2020_plot_only[!(pops_2020_plot_only%in%pops_only_2020$id_num_2020)] <- NA
} else {
  pops_2020_plot_both <- pops_2020_plot
  pops_2020_plot_only <- pops_2020_plot
}

if(nrow(pops_only_2050_df) >= 1) {
  pops_2050_plot_both <- pops_2050_plot
  pops_2050_plot_both[pops_2020_plot_both%in%pops_only_2050_df$id_num_2050] <- NA
  
  pops_2050_plot_only <- pops_2050_plot
  pops_2050_plot_only[!(pops_2050_plot_only%in%pops_only_2050_df$id_num_2050)] <- NA
} else {
  pops_2050_plot_both <- pops_2050_plot
  pops_2050_plot_only <- pops_2050_plot
}


plot(rast(pops_2020_plot_both),
     axes = FALSE, 
     type = 'continuous',
     col = viridis::rocket(99),
     # legend = 'bottom',
     legend = FALSE,
     main = '2020',
     range = c(1,max(maxValue(pops_2020_plot_both), maxValue(pops_2050_both_out))))

if(nrow(pops_only_2020_df) >= 1) {
  ### Getting data to create legend
  # Number of populations
  num_pops_2020 <- nrow(pops_only_2020_df)
  
  # Colour palettes
  col_palette_2020_only <- 
    colorRampPalette(brewer.pal(n=9, 'Blues'))(99)[seq(from = 25, to = 99, l = num_pops_2020)]
  
  plot(rast(pops_2020_plot_only),
       axes = FALSE,
       legend = FALSE,
       add = TRUE,
       col = col_palette_2020_only)
}


plot(rast(pops_2050_both_out),
     axes = FALSE, 
     # legend = 'bottom',
     col = viridis::rocket(99),
     type = 'continuous',
     legend = FALSE,
     # add = TRUE,
     main = '2050',
     # range = c(1,100))
     range = c(1,max(maxValue(pops_2020_plot_both), maxValue(pops_2050_both_out))))

if(nrow(pops_only_2050_df) >= 1) {
  ### Getting data to create legend
  # Number of populations
  num_pops_2050 <- nrow(pops_only_2050_df)
  
  # Colour palettes
  col_palette_2050_only <- 
    colorRampPalette(brewer.pal(n=9, 'Greens'))(99)[seq(from = 25, to = 99, l = num_pops_2050)]
  
  plot(rast(pops_2050_plot_only),
       axes = FALSE,
       legend = FALSE,
       add = TRUE,
       col = col_palette_2050_only)
}


plot(c(-3,3),c(0-legend_height,1+legend_height),type = 'n', axes = F,xlab = '', ylab = '')

# Creating legends - for populations in both 2020 and 2050
# text(x = 0,y=4,'Population ID Number',font = 2)

###
# Getting starting point for legend text
# Labels for legend

# for pops only in 2020
labels_2020_pop_legend <-
  seq(from = 1,
      to = num_pops_2020,
      l = min(num_pops_2020, 3))

# for pops only in both 2020 and 2050
labels_start_both <- max(labels_2020_pop_legend) + 1

# for pops only in 2050
labels_start_2050 <- 
  labels_start_both + 
  max(minmax(rast(pops_2020_plot_both)),
      minmax(rast(pops_2050_plot_both)))

# Labels for legend
labels_2050_pop_legend <-
  seq(from = labels_start_2050,
      to = labels_start_2050 + nrow(pops_only_2050_df),
      l = 2)

legend_image <- as.raster(matrix(viridis::rocket(99)[seq(from = 1, to = 99, l = max(maxValue(pops_2020_plot_both), maxValue(pops_2050_both_out)))], nrow=1))
text(y=-.5, 
     x = seq(-.5,.5,l=2), 
     labels = seq(label_start_both,
                  (max(minmax(rast(pops_2020_plot_both)),
                       minmax(rast(pops_2050_plot_both))) + 
                     (labels_start_both - 1)),
                  l=2))
text(y=3, x = 0, labels = 'Population ID Number\n(Populations in both 2020 and 2050)\n',adj = .5)
rasterImage(legend_image, -.5,0, .5,1)

# Creating legend - for population in 2020 not 2050
if(nrow(pops_only_2020_df) >= 1) {
  legend_image <- as.raster(matrix(col_palette_2020_only, nrow=1))
  text(y=-.5, x = seq(-2.5,-1.5,l=2), labels = labels_2020_pop_legend)
  text(y=3, x = -2, labels = 'Population ID Number\n(Populations in 2020\nnot present in 2050)')
  rasterImage(legend_image, -2.5,0,-1.5,1)
} else {
  text(y=3, x = -2, labels = 'No Populations in 2020\nthat are not present in 2050\n')
}


# Creating legend - for population in 2050 not 2020
if(nrow(pops_only_2050_df) >= 1) {
  legend_image <- as.raster(matrix(col_palette_2050_only, nrow=1))
  text(y=-.5, x = seq(1.5,2.5,l = 2), labels = labels_2050_pop_legend)
  text(y=3, x = 2, labels = 'Population ID Number\n(Populations in 2050\nnot present in 2020)')
  rasterImage(legend_image, 1.5,0,2.5,1)
} else {
  text(y=3, x = 2, labels = 'No Populations in 2050\nthat are not present in 2020\n')
}



