  #####
  ###
  # Fig 5 - species level outcomes
  # Four panels
  # (a) Total outcome in 2050 including climate suitable habitat, habitat remaining, patches, and populations
  # (b) - (d) Change in outcomes from 2020 - 2050 for ag exp, urb, and climate change
  
  
  
  ###
  # Libraries
  library(terra)
  library(sp)
  library(RColorBrewer)
  library(plyr)
  library(dplyr)
  library(viridis)
  library(viridisLite)
  
  ###
  # Background data management
  
  # Setting colour palette
  my_pal <-
    rev(c(colorRampPalette(c(brewer.pal('GnBu',n=6)[6:3],'#eeeeee'))(round(99 * 3.5/7)),
      colorRampPalette(c('#eeeeee',
                         '#F8EBA3',
                         '#F0B082',
                         '#C94A6B',
                         '#aa00aa'))(round(99*3.5/7))))
  
  # Importing files
  clim_suit <- rast('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Phylloscartes_oustaleti_SSP2_prevalence_threshold__2050.tif')
  baseline_2020 <- rast('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/TempFigs_14Nov2024/hab_availability_Phylloscartes_oustaleti_prevalence_2020_esh.future.extent.esh.exp.int.urb.tif')
  all_2050 <- rast('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/TempFigs_14Nov2024/hab_availability_Phylloscartes_oustaleti_prevalence_2050_esh.future.extent.esh.exp.int.urb.tif')
  abund_2050 <- rast('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/TempFigs_14Nov2024/pop_density_Phylloscartes_oustaleti_prevalence_2050_esh.future.extent.esh.exp.int.urb.tif')
  ag_2050 <- rast('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/TempFigs_14Nov2024/hab_availability_Phylloscartes_oustaleti_prevalence_2050_esh.future.extent.esh.exp.tif')
  urb_2050 <- rast('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/TempFigs_14Nov2024/hab_availability_Phylloscartes_oustaleti_prevalence_2050_esh.future.extent.esh.urb.tif')
  clim_2050 <- rast('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/TempFigs_14Nov2024/hab_availability_Phylloscartes_oustaleti_prevalence_2050_esh.current.extent.esh.tif')
  patches_2050 <- rast('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/TempFigs_14Nov2024/patches_2050_Phylloscartes_oustaleti.tif')
  pops_2050 <- rast('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/TempFigs_14Nov2024/pops_2050_Phylloscartes_oustaleti.tif')
  
  # Getting states of brazil
  brazil_states <- ne_states(country = 'Brazil', returnclass = 'sv')
  brazil_states <- project(brazil_states, clim_2050)
  
  # Updating extents
  max_extent <-
    ext(c(min(ext(baseline_2020)[1],ext(all_2050)[1]),
      max(ext(baseline_2020)[2],ext(all_2050)[2]),
      min(ext(baseline_2020)[3],ext(all_2050)[3]),
      max(ext(baseline_2020)[4],ext(all_2050)[4])))
  
  baseline_2020 <- extend(baseline_2020, max_extent)
  all_2050 <- extend(all_2050,max_extent)
  ag_2050 <- extend(ag_2050,max_extent)
  urb_2050 <- extend(urb_2050,max_extent)
  clim_2050 <- extend(clim_2050,max_extent)
  patches_2050 <- extend(patches_2050,max_extent)
  pops_2050 <- extend(pops_2050, max_extent)
  mollweide_id_full <- rast('/Users/michael/Desktop/MollweideCountryID_1.5km.tif')
  
  # Updating clim suit map to remove oceans, etc
  mollweide_id <- crop(mollweide_id_full, clim_suit)
  clim_suit[is.na(mollweide_id)] <- NA
  clim_vect <- as.polygons(clim_suit)
  clim_sp <- as(clim_vect,'Spatial')
  
  # Getting change in available habitat - 2020 to 2050 for individual stressors
  delta_ag <- ag_2050 - baseline_2020
  delta_urb <- urb_2050 - baseline_2020
  
  clim_2050_delta <- clim_2050
  clim_2050_delta[is.na(clim_2050_delta)] <- 0
  baseline_2020_clim <- baseline_2020
  baseline_2020_clim[is.na(baseline_2020_clim)] <- 0
  delta_clim <- clim_2050_delta - baseline_2020_clim
  delta_clim[is.na(baseline_2020) & is.na(clim_2050)] <- NA
  patches_2050[is.na(pops_2050)] <- NA
  
  ###
  # Panel a
  pops_2050_poly <- as.polygons(pops_2050)
  pops_sf_1 <- as(pops_2050_poly[pops_2050_poly$pops_2050_Phylloscartes_oustaleti %in% 1], 'Spatial')
  pops_sf_2 <- as(pops_2050_poly[pops_2050_poly$pops_2050_Phylloscartes_oustaleti %in% 2], 'Spatial')
  clim_sp <- crop(clim_vect, 
                  round(all_2050,digits = 2) %>% extend(.,450))
  clim_sp <- as(clim_sp,'Spatial')
  
  # brazil_states_rast <- rasterize(brazil_states, mollweide_id_full)
  # dev.off()
  pdf('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Figures and Tables/Figures/Fig5a_Species_HabAvail_2025-05-17.pdf',
      width = 9 / 2.54,
      height = 9 / 2.54)
  par(mar = c(3,1,1,1))
  
  
  terra::plot(clim_sp,
       box = FALSE,
       legend = FALSE,
       axes = FALSE)
       # add = TRUE,
       # background = '#eeeeee',
       # density = 15,
       # angle = 45)
  
  plot(brazil_states,
       add = TRUE,
       box = FALSE,
       axes = FALSE,
       legend = FALSE,
       col = '#eeeeee',
       border = 'grey',
       xpd = TRUE)
  
  terra::plot(clim_sp,
       box = FALSE,
       legend = FALSE,
       axes = FALSE,
       add = TRUE)
       # background = '#eeeeee',
       # density = 15,
       # angle = 45)
  
  plot(round(all_2050,digits = 2) %>% extend(.,150),
       box = FALSE,
       legend = FALSE,
       axes = FALSE,
       range = c(0,1),
       add = TRUE,
       col = my_pal[51:100])
  
  plot(brazil_states,
       add = TRUE,
       box = FALSE,
       axes = FALSE,
       legend = FALSE,
       border = 'grey')
  
  
  plot(clim_sp,
       # density = 10,
       # angle = 45,
       # col = '#dddddd',
       box = FALSE,
       axes = FALSE,
       legend = FALSE,
       add = TRUE)
  
  
  # Define rectangle bounds
  xleft <- .53
  xright <- .93
  ybottom <- .05
  ytop <- .1
  
  # Convert to user coordinates
  xleft <- grconvertX(xleft, from = "ndc", to = "user")
  xrightlabel <- grconvertX(xright+.025, from = "ndc", to = "user")
  xright <- grconvertX(xright, from = "ndc", to = "user")
  ybottom_label <- grconvertY(ybottom - .025, from = "ndc", to = "user")
  ybottom <- grconvertY(ybottom, from = "ndc", to = "user")
  ytitlelabel <- grconvertY(ytop+.05, from = "ndc", to = "user")
  ytop <- grconvertY(ytop, from = "ndc", to = "user")
  
  # Compute width of each sub-rectangle
  my_pal_tmp <- my_pal[51:100]
  x_steps <- seq(xleft, xright, length.out = length(my_pal_tmp) + 1)
  
  # Draw gradient using thin vertical rectangles
  for (i in 1:length(my_pal_tmp)) {
    rect(xleft = x_steps[i], ybottom = ybottom,
         xright = x_steps[i+1], ytop = ytop,
         col = my_pal_tmp[i], border = NA,
         xpd = TRUE)
  }
  
  rect(xleft, ybottom, xright, ytop, border = "black", xpd = TRUE)
  
  text(x = xright, 
       y = ybottom_label, 
       label = '1',
       xpd = TRUE,
       adj = c(.5,.5),
       cex = .6)
  
  text(x = xleft, 
       y = ybottom_label, 
       label = '0',
       xpd = TRUE,
       adj = c(.5,.5),
       cex = .6)
  
  text(x = mean(c(xright,xleft)), 
       y = ybottom_label,
       label = '.5',
       xpd = TRUE,
       adj = c(.5,.5),
       cex = .6)
  
  text(x = mean(c(xleft,xright)),
       y = ytitlelabel,
       label = 'Available Habitat in 2050\n(Proportion of Cell)',
       adj = c(0.5,.5),
       cex = .6,
       xpd = TRUE)
  
  text(x = grconvertX(.00001, from = 'ndc', to = 'user'),
       y = grconvertY(.99999, from = 'ndc', to = 'user'),
       label = '(a)',
       adj = c(0,1),
       xpd = TRUE,
       cex = .6)
  
  rect(xleft = grconvertX(.53, from = 'ndc', to = 'user'),
           xright = grconvertX(.58, from = 'ndc', to = 'user'),
           ybottom = grconvertY(.225, from = 'ndc', to = 'user'),
           ytop = grconvertY(.275, from = 'ndc', to = 'user'),
           lwd = 1,
           # col = 'white',
           xpd = TRUE)
       # density = 15,
       # angle = 45)
  
  text(x = grconvertX(.7, from = 'ndc', to = 'user'),
       y = grconvertY(.25, from = 'ndc', to = 'user'),
       label = 'Climate Suitable\nHabitat Area',
       adj = c(.5,.5),
       xpd = TRUE,
       cex = .6)
  
  # Adding distance scale
  tmp_df <- 
    round(all_2050,digits = 2) %>% 
    extend(.,150) %>%
    as.data.frame(xy = TRUE) %>%
    dplyr::select(x,y) %>% 
    distinct()
  
  tmp_df_x_min <- quantile(tmp_df$x,.73) - 150000
  tmp_df_x_max <- quantile(tmp_df$x,.73) + 150000
  
  y_pos <- quantile(unique(tmp_df$y),.05)
  
  y_pos_new <- grconvertY(y = 0.01, from = 'ndc', to = 'user')
  tmp_df_x_min_user <- grconvertX(tmp_df_x_min, from = 'user', to = 'ndc')
  tmp_df_x_max_user <- grconvertX(tmp_df_x_max, from = 'user', to = 'ndc')
  tmp_df_x_min_dif <- tmp_df_x_max_user - tmp_df_x_min_user
  tmp_df_x_min_new <- grconvertX(x = 0.025, from = 'ndc', to = 'user')
  tmp_df_x_max_new <- grconvertX(x = 0.025 + tmp_df_x_min_dif, from = 'ndc', to = 'user')
  
  segments(x0 = tmp_df_x_min_new,
           x = tmp_df_x_max_new,
           y0 = y_pos_new,
           y = y_pos_new,
           lwd = 1,
           xpd = TRUE)
  
  segments(x0 = tmp_df_x_min_new,
           x = tmp_df_x_min_new,
           y0 = y_pos_new,
           y = y_pos_new + 15000,
           lwd = 1,
           xpd = TRUE)
  
  segments(x0 = tmp_df_x_max_new,
           x = tmp_df_x_max_new,
           y0 = y_pos_new,
           y = y_pos_new + 15000,
           lwd = 1,
           xpd = TRUE)
  
  segments(x0 = mean(c(tmp_df_x_max_new,tmp_df_x_min_new)),
           x = mean(c(tmp_df_x_max_new,tmp_df_x_min_new)),
           y0 = y_pos_new,
           y = y_pos_new + 15000,
           lwd = 1,
           xpd = TRUE)
  
  text(x = tmp_df_x_min_new,
       y = y_pos_new + 25000,
       label = '0km',
       adj = c(0,0),
       cex = .6,
       xpd = TRUE,
       srt = 45)
  
  text(x = mean(c(tmp_df_x_max_new,tmp_df_x_min_new)),
       y = y_pos_new + 25000,
       label = '150km',
       adj = c(0,0),
       cex = .6,
       xpd = TRUE, 
       srt = 45)
  
  text(x = tmp_df_x_max_new,
       y = y_pos_new + 25000,
       label = '300km',
       adj = c(0,0),
       cex = .6,
       xpd = TRUE, 
       srt = 45)
  
  
  
  xleft_inset <- grconvertX(0,from = 'ndc', to = 'user')
  xright_inset <- grconvertX(1,from = 'ndc', to = 'user')
  ybottom_inset <- grconvertY(0,from = 'ndc', to = 'user')
  ytop_inset <- grconvertY(1,from = 'ndc', to = 'user')
  
  # Creating inset map, showing where this is in Brazil
  par(fig = c(0.05, 0.3, 0.7, 0.95), new = TRUE)
  par(mar = c(0,0,0,0))
  
  plot(brazil_states,
       box = FALSE,
       background = 'white',
       col = '#eeeeee',
       axes = FALSE,
       legend = FALSE,
       # border = 'white',
       mar = c(0,0,0,0))
  
  rect(xleft = xleft_inset,
       xright = xright_inset,
       ybottom = ybottom_inset,
       ytop = ytop_inset,
       lwd = 2)
  
  
  dev.off()

###
# Panel b - pop abundance by quintiles
num_breaks <- 5
breaks <- quantile(values(abund_2050), probs = seq(0, 1, 1 / num_breaks), na.rm = TRUE, type = 7)

# Reclassify raster into quintiles
# The 'classify' function needs a matrix of from-to-newvalue
rc_matrix <- cbind(
  breaks[-length(breaks)],
  breaks[-1],
  1:num_breaks
)

# Handle edge cases to include the maximum value
rc_matrix[num_breaks, 2] <- rc_matrix[num_breaks, 2] + 1e-10  # tiny increment to include max

# Apply the classification
quintile_raster <- classify(abund_2050, rc_matrix, include.lowest = TRUE)
# plot(quintile_raster) # Checking

pdf('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Figures and Tables/Figures/Fig5b_Species_PopAbund_2025-05-17.pdf',
    width = 9 / 2.54,
    height = 9 / 2.54)
par(mar = c(3,1,1,1))


plot(clim_sp,
     box = FALSE,
     legend = FALSE,
     axes = FALSE)
     # add = TRUE,
     # density = 15,
     # angle = 45)

plot(brazil_states,
     add = TRUE,
     box = FALSE,
     axes = FALSE,
     legend = FALSE,
     col = '#eeeeee',
     border = 'grey',
     xpd = TRUE)

plot(clim_sp,
     box = FALSE,
     legend = FALSE,
     axes = FALSE,
     add = TRUE)
     # density = 15,
     # angle = 45)

plot(quintile_raster,
     box = FALSE,
     legend = FALSE,
     axes = FALSE,
     range = c(0,1),
     add = TRUE)


plot(brazil_states,
     add = TRUE,
     box = FALSE,
     axes = FALSE,
     legend = FALSE,
     border = 'grey')

# plot(brazil_states,
#      add = TRUE,
#      box = FALSE,
#      axes = FALSE,
#      legend = FALSE,
#      border = 'grey')


plot(clim_sp,
     # density = 10,
     # angle = 45,
     # col = '#dddddd',
     box = FALSE,
     axes = FALSE,
     legend = FALSE,
     add = TRUE)


# Define rectangle bounds
legend(x = grconvertX(.51, from = 'ndc', to = 'user'),
       y = grconvertY(.01, from = 'ndc', to = 'user'), 
       legend = 1:5, fill = viridis(n = 5),
       bty = "n", cex = 0.6, inset = c(0, 0),
       title = 'Population Density Quintile\n(1 = Lowest; 5 = Highest Density)',
       hor = TRUE,
       # y.intersp = .8,
       xpd = TRUE,
       xjust = 0,
       yjust = 0)

rect(xleft = grconvertX(.53, from = 'ndc', to = 'user'),
     xright = grconvertX(.58, from = 'ndc', to = 'user'),
     ybottom = grconvertY(.225, from = 'ndc', to = 'user'),
     ytop = grconvertY(.275, from = 'ndc', to = 'user'),
     lwd = 1,
     # col = 'white',
     xpd = TRUE)#,
     # density = 15,
     # angle = 45)

text(x = grconvertX(.7, from = 'ndc', to = 'user'),
     y = grconvertY(.25, from = 'ndc', to = 'user'),
     label = 'Climate Suitable\nHabitat Area',
     adj = c(.5,.5),
     xpd = TRUE,
     cex = .6)


# Adding distance scale
tmp_df <- 
  round(all_2050,digits = 2) %>% 
  extend(.,150) %>%
  as.data.frame(xy = TRUE) %>%
  dplyr::select(x,y) %>% 
  distinct()

tmp_df_x_min <- quantile(tmp_df$x,.73) - 150000
tmp_df_x_max <- quantile(tmp_df$x,.73) + 150000

y_pos <- quantile(unique(tmp_df$y),.05)

y_pos_new <- grconvertY(y = 0.01, from = 'ndc', to = 'user')
tmp_df_x_min_user <- grconvertX(tmp_df_x_min, from = 'user', to = 'ndc')
tmp_df_x_max_user <- grconvertX(tmp_df_x_max, from = 'user', to = 'ndc')
tmp_df_x_min_dif <- tmp_df_x_max_user - tmp_df_x_min_user
tmp_df_x_min_new <- grconvertX(x = 0.025, from = 'ndc', to = 'user')
tmp_df_x_max_new <- grconvertX(x = 0.025 + tmp_df_x_min_dif, from = 'ndc', to = 'user')

segments(x0 = tmp_df_x_min_new,
         x = tmp_df_x_max_new,
         y0 = y_pos_new,
         y = y_pos_new,
         lwd = 1,
         xpd = TRUE)

segments(x0 = tmp_df_x_min_new,
         x = tmp_df_x_min_new,
         y0 = y_pos_new,
         y = y_pos_new + 15000,
         lwd = 1,
         xpd = TRUE)

segments(x0 = tmp_df_x_max_new,
         x = tmp_df_x_max_new,
         y0 = y_pos_new,
         y = y_pos_new + 15000,
         lwd = 1,
         xpd = TRUE)

segments(x0 = mean(c(tmp_df_x_max_new,tmp_df_x_min_new)),
         x = mean(c(tmp_df_x_max_new,tmp_df_x_min_new)),
         y0 = y_pos_new,
         y = y_pos_new + 15000,
         lwd = 1,
         xpd = TRUE)

text(x = tmp_df_x_min_new,
     y = y_pos_new + 25000,
     label = '0km',
     adj = c(0,0),
     cex = .6,
     xpd = TRUE,
     srt = 45)

text(x = mean(c(tmp_df_x_max_new,tmp_df_x_min_new)),
     y = y_pos_new + 25000,
     label = '150km',
     adj = c(0,0),
     cex = .6,
     xpd = TRUE, 
     srt = 45)

text(x = tmp_df_x_max_new,
     y = y_pos_new + 25000,
     label = '300km',
     adj = c(0,0),
     cex = .6,
     xpd = TRUE, 
     srt = 45)

text(x = grconvertX(.00001, from = 'ndc', to = 'user'),
     y = grconvertY(.99999, from = 'ndc', to = 'user'),
     label = '(b)',
     adj = c(0,1),
     xpd = TRUE,
     cex = .6)


dev.off()


###
# Panel c
# Updating y position for new format (scale)



pdf('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Figures and Tables/Figures/Fig5c_Species_Clim_2025-05-17.pdf',
    width = 9 / 2.54 * 2/3,
    height = 9 / 2.54 * 2/3)
par(mar = c(2,1,1,1))
mollweide_id <- crop(mollweide_id, delta_clim)
delta_clim[is.na(mollweide_id)] <- NA

plot(clim_sp,
     box = FALSE,
     legend = FALSE,
     axes = FALSE,
     col = 'white',
     border = 'white',
     # add = TRUE,
     density = 15,
     angle = 45)

y_pos_new <- grconvertY(y = 0.01, from = 'ndc', to = 'user')
tmp_df_x_min_user <- grconvertX(tmp_df_x_min, from = 'user', to = 'ndc')
tmp_df_x_max_user <- grconvertX(tmp_df_x_max, from = 'user', to = 'ndc')
tmp_df_x_min_dif <- tmp_df_x_max_user - tmp_df_x_min_user
tmp_df_x_min_new <- grconvertX(x = 0.025, from = 'ndc', to = 'user')
tmp_df_x_max_new <- grconvertX(x = 0.025 + tmp_df_x_min_dif, from = 'ndc', to = 'user')

plot(brazil_states,
     add = TRUE,
     box = FALSE,
     axes = FALSE,
     legend = FALSE,
     col = '#eeeeee',
     border = 'grey',
     xpd = TRUE)

plot(round(delta_clim,digits = 2),
     box = FALSE,
     legend = FALSE,
     axes = FALSE,
     add = TRUE,
     range = c(-1,1),
     col = (my_pal))


plot(brazil_states,
     add = TRUE,
     box = FALSE,
     axes = FALSE,
     legend = FALSE,
     border = 'grey')


# Define rectangle bounds
xleft <- .53
xright <- .93
ybottom <- .05
ytop <- .1

# Convert to user coordinates
xleft <- grconvertX(xleft, from = "ndc", to = "user")
xrightlabel <- grconvertX(xright+.025, from = "ndc", to = "user")
xright <- grconvertX(xright, from = "ndc", to = "user")
ybottom_label <- grconvertY(ybottom - .025, from = "ndc", to = "user")
ybottom <- grconvertY(ybottom, from = "ndc", to = "user")
ytitlelabel <- grconvertY(ytop+.1, from = "ndc", to = "user")
ytop <- grconvertY(ytop, from = "ndc", to = "user")

# Compute width of each sub-rectangle
my_pal_tmp <- my_pal
x_steps <- seq(xleft, xright, length.out = length(my_pal_tmp) + 1)

# Draw gradient using thin vertical rectangles
for (i in 1:length(my_pal_tmp)) {
  rect(xleft = x_steps[i], ybottom = ybottom,
       xright = x_steps[i+1], ytop = ytop,
       col = my_pal_tmp[i], border = NA,
       xpd = TRUE)
}

rect(xleft, ybottom, xright, ytop, border = "black", xpd = TRUE)



text(x = xright, 
     y = ybottom_label, 
     label = '1',
     xpd = TRUE,
     adj = c(.5,.5),
     cex = .6)

text(x = xleft, 
     y = ybottom_label, 
     label = '-1',
     xpd = TRUE,
     adj = c(.5,.5),
     cex = .6)

text(x = mean(c(xright,xleft)), 
     y = ybottom_label,
     label = '0',
     xpd = TRUE,
     adj = c(.5,.5),
     cex = .6)

text(x = mean(c(xleft,xright)),
     y = ytitlelabel,
     label = 'Change in Available\nHabitat from 2020 - 2050\n(Proportion of Cell)',
     adj = c(0.5,.5),
     cex = .6,
     xpd = TRUE)

text(x = grconvertX(.00001, from = 'ndc', to = 'user'),
     y = grconvertY(.99999, from = 'ndc', to = 'user'),
     label = '(e)',
     adj = c(0,1),
     xpd = TRUE,
     cex = .6)


segments(x0 = tmp_df_x_min_new,
         x = tmp_df_x_max_new,
         y0 = y_pos_new,
         y = y_pos_new,
         lwd = 1,
         xpd = TRUE)

segments(x0 = tmp_df_x_min_new,
         x = tmp_df_x_min_new,
         y0 = y_pos_new,
         y = y_pos_new+15000,
         lwd = 1,
         xpd = TRUE)

segments(x0 = tmp_df_x_max_new,
         x = tmp_df_x_max_new,
         y0 = y_pos_new,
         y = y_pos_new+15000,
         lwd = 1,
         xpd = TRUE)

segments(x0 = mean(c(tmp_df_x_max_new,tmp_df_x_min_new)),
         x = mean(c(tmp_df_x_max_new,tmp_df_x_min_new)),
         y0 = y_pos_new,
         y = y_pos_new+15000,
         lwd = 1,
         xpd = TRUE)

text(x = tmp_df_x_min_new,
     y = y_pos_new+25000,
     label = '0km',
     adj = c(0,0),
     cex = .6,
     srt = 45,
     xpd = TRUE)

text(x = mean(c(tmp_df_x_max_new,tmp_df_x_min_new)),
     y = y_pos_new + 25000,
     label = '150km',
     adj = c(0,0),
     cex = .6,
     srt = 45,
     xpd = TRUE)

text(x = tmp_df_x_max_new,
     y = y_pos_new + 25000,
     label = '300km',
     adj = c(0,0),
     cex = .6,
     srt = 45,
     xpd = TRUE)

dev.off()



###
# Panel d
pdf('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Figures and Tables/Figures/Fig5d_Species_AgExp_2025-05-17.pdf',
    width = 9 / 2.54 * 2/3,
    height = 9 / 2.54 * 2/3)
par(mar = c(2,1,1,1))
mollweide_id <- crop(mollweide_id, delta_clim)
delta_clim[is.na(mollweide_id)] <- NA

plot(clim_sp,
     box = FALSE,
     legend = FALSE,
     axes = FALSE,
     col = 'white',
     border = 'white',
     # add = TRUE,
     density = 15,
     angle = 45)

y_pos_new <- grconvertY(y = 0.01, from = 'ndc', to = 'user')
tmp_df_x_min_user <- grconvertX(tmp_df_x_min, from = 'user', to = 'ndc')
tmp_df_x_max_user <- grconvertX(tmp_df_x_max, from = 'user', to = 'ndc')
tmp_df_x_min_dif <- tmp_df_x_max_user - tmp_df_x_min_user
tmp_df_x_min_new <- grconvertX(x = 0.025, from = 'ndc', to = 'user')
tmp_df_x_max_new <- grconvertX(x = 0.025 + tmp_df_x_min_dif, from = 'ndc', to = 'user')

plot(brazil_states,
     add = TRUE,
     box = FALSE,
     axes = FALSE,
     legend = FALSE,
     col = '#eeeeee',
     border = 'grey',
     xpd = TRUE)

plot(round(delta_ag,digits = 2),
     box = FALSE,
     legend = FALSE,
     axes = FALSE,
     add = TRUE,
     range = c(-1,1),
     col = (my_pal))


plot(brazil_states,
     add = TRUE,
     box = FALSE,
     axes = FALSE,
     legend = FALSE,
     border = 'grey')


# Define rectangle bounds
xleft <- .53
xright <- .93
ybottom <- .05
ytop <- .1

# Convert to user coordinates
xleft <- grconvertX(xleft, from = "ndc", to = "user")
xrightlabel <- grconvertX(xright+.025, from = "ndc", to = "user")
xright <- grconvertX(xright, from = "ndc", to = "user")
ybottom_label <- grconvertY(ybottom - .025, from = "ndc", to = "user")
ybottom <- grconvertY(ybottom, from = "ndc", to = "user")
ytitlelabel <- grconvertY(ytop+.1, from = "ndc", to = "user")
ytop <- grconvertY(ytop, from = "ndc", to = "user")

# Compute width of each sub-rectangle
my_pal_tmp <- my_pal
x_steps <- seq(xleft, xright, length.out = length(my_pal_tmp) + 1)

# Draw gradient using thin vertical rectangles
for (i in 1:length(my_pal_tmp)) {
  rect(xleft = x_steps[i], ybottom = ybottom,
       xright = x_steps[i+1], ytop = ytop,
       col = my_pal_tmp[i], border = NA,
       xpd = TRUE)
}

rect(xleft, ybottom, xright, ytop, border = "black", xpd = TRUE)



text(x = xright, 
     y = ybottom_label, 
     label = '1',
     xpd = TRUE,
     adj = c(.5,.5),
     cex = .6)

text(x = xleft, 
     y = ybottom_label, 
     label = '-1',
     xpd = TRUE,
     adj = c(.5,.5),
     cex = .6)

text(x = mean(c(xright,xleft)), 
     y = ybottom_label,
     label = '0',
     xpd = TRUE,
     adj = c(.5,.5),
     cex = .6)

text(x = mean(c(xleft,xright)),
     y = ytitlelabel,
     label = 'Change in Available\nHabitat from 2020 - 2050\n(Proportion of Cell)',
     adj = c(0.5,.5),
     cex = .6,
     xpd = TRUE)

text(x = grconvertX(.00001, from = 'ndc', to = 'user'),
     y = grconvertY(.99999, from = 'ndc', to = 'user'),
     label = '(c)',
     adj = c(0,1),
     xpd = TRUE,
     cex = .6)


segments(x0 = tmp_df_x_min_new,
         x = tmp_df_x_max_new,
         y0 = y_pos_new,
         y = y_pos_new,
         lwd = 1,
         xpd = TRUE)

segments(x0 = tmp_df_x_min_new,
         x = tmp_df_x_min_new,
         y0 = y_pos_new,
         y = y_pos_new+15000,
         lwd = 1,
         xpd = TRUE)

segments(x0 = tmp_df_x_max_new,
         x = tmp_df_x_max_new,
         y0 = y_pos_new,
         y = y_pos_new+15000,
         lwd = 1,
         xpd = TRUE)

segments(x0 = mean(c(tmp_df_x_max_new,tmp_df_x_min_new)),
         x = mean(c(tmp_df_x_max_new,tmp_df_x_min_new)),
         y0 = y_pos_new,
         y = y_pos_new+15000,
         lwd = 1,
         xpd = TRUE)

text(x = tmp_df_x_min_new,
     y = y_pos_new+25000,
     label = '0km',
     adj = c(0,0),
     cex = .6,
     srt = 45,
     xpd = TRUE)

text(x = mean(c(tmp_df_x_max_new,tmp_df_x_min_new)),
     y = y_pos_new + 25000,
     label = '150km',
     adj = c(0,0),
     cex = .6,
     srt = 45,
     xpd = TRUE)

text(x = tmp_df_x_max_new,
     y = y_pos_new + 25000,
     label = '300km',
     adj = c(0,0),
     cex = .6,
     srt = 45,
     xpd = TRUE)


dev.off()

###
# Panel e

pdf('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Figures and Tables/Figures/Fig5e_Species_Urb_2025-05-17.pdf',
    width = 9 / 2.54 * 2/3,
    height = 9 / 2.54 * 2/3)
par(mar = c(2,1,1,1))
mollweide_id <- crop(mollweide_id, delta_clim)
delta_clim[is.na(mollweide_id)] <- NA

plot(clim_sp,
     box = FALSE,
     legend = FALSE,
     axes = FALSE,
     col = 'white',
     border = 'white',
     # add = TRUE,
     density = 15,
     angle = 45)

y_pos_new <- grconvertY(y = 0.01, from = 'ndc', to = 'user')
tmp_df_x_min_user <- grconvertX(tmp_df_x_min, from = 'user', to = 'ndc')
tmp_df_x_max_user <- grconvertX(tmp_df_x_max, from = 'user', to = 'ndc')
tmp_df_x_min_dif <- tmp_df_x_max_user - tmp_df_x_min_user
tmp_df_x_min_new <- grconvertX(x = 0.025, from = 'ndc', to = 'user')
tmp_df_x_max_new <- grconvertX(x = 0.025 + tmp_df_x_min_dif, from = 'ndc', to = 'user')

plot(brazil_states,
     add = TRUE,
     box = FALSE,
     axes = FALSE,
     legend = FALSE,
     col = '#eeeeee',
     border = 'grey',
     xpd = TRUE)

plot(round(delta_urb,digits = 2),
     box = FALSE,
     legend = FALSE,
     axes = FALSE,
     add = TRUE,
     range = c(-1,1),
     col = (my_pal))


plot(brazil_states,
     add = TRUE,
     box = FALSE,
     axes = FALSE,
     legend = FALSE,
     border = 'grey')


# Define rectangle bounds
xleft <- .53
xright <- .93
ybottom <- .05
ytop <- .1

# Convert to user coordinates
xleft <- grconvertX(xleft, from = "ndc", to = "user")
xrightlabel <- grconvertX(xright+.025, from = "ndc", to = "user")
xright <- grconvertX(xright, from = "ndc", to = "user")
ybottom_label <- grconvertY(ybottom - .025, from = "ndc", to = "user")
ybottom <- grconvertY(ybottom, from = "ndc", to = "user")
ytitlelabel <- grconvertY(ytop+.1, from = "ndc", to = "user")
ytop <- grconvertY(ytop, from = "ndc", to = "user")

# Compute width of each sub-rectangle
my_pal_tmp <- my_pal
x_steps <- seq(xleft, xright, length.out = length(my_pal_tmp) + 1)

# Draw gradient using thin vertical rectangles
for (i in 1:length(my_pal_tmp)) {
  rect(xleft = x_steps[i], ybottom = ybottom,
       xright = x_steps[i+1], ytop = ytop,
       col = my_pal_tmp[i], border = NA,
       xpd = TRUE)
}

rect(xleft, ybottom, xright, ytop, border = "black", xpd = TRUE)



text(x = xright, 
     y = ybottom_label, 
     label = '1',
     xpd = TRUE,
     adj = c(.5,.5),
     cex = .6)

text(x = xleft, 
     y = ybottom_label, 
     label = '-1',
     xpd = TRUE,
     adj = c(.5,.5),
     cex = .6)

text(x = mean(c(xright,xleft)), 
     y = ybottom_label,
     label = '0',
     xpd = TRUE,
     adj = c(.5,.5),
     cex = .6)

text(x = mean(c(xleft,xright)),
     y = ytitlelabel,
     label = 'Change in Available\nHabitat from 2020 - 2050\n(Proportion of Cell)',
     adj = c(0.5,.5),
     cex = .6,
     xpd = TRUE)

text(x = grconvertX(.00001, from = 'ndc', to = 'user'),
     y = grconvertY(.99999, from = 'ndc', to = 'user'),
     label = '(d)',
     adj = c(0,1),
     xpd = TRUE,
     cex = .6)


segments(x0 = tmp_df_x_min_new,
         x = tmp_df_x_max_new,
         y0 = y_pos_new,
         y = y_pos_new,
         lwd = 1,
         xpd = TRUE)

segments(x0 = tmp_df_x_min_new,
         x = tmp_df_x_min_new,
         y0 = y_pos_new,
         y = y_pos_new+15000,
         lwd = 1,
         xpd = TRUE)

segments(x0 = tmp_df_x_max_new,
         x = tmp_df_x_max_new,
         y0 = y_pos_new,
         y = y_pos_new+15000,
         lwd = 1,
         xpd = TRUE)

segments(x0 = mean(c(tmp_df_x_max_new,tmp_df_x_min_new)),
         x = mean(c(tmp_df_x_max_new,tmp_df_x_min_new)),
         y0 = y_pos_new,
         y = y_pos_new+15000,
         lwd = 1,
         xpd = TRUE)

text(x = tmp_df_x_min_new,
     y = y_pos_new+25000,
     label = '0km',
     adj = c(0,0),
     cex = .6,
     srt = 45,
     xpd = TRUE)

text(x = mean(c(tmp_df_x_max_new,tmp_df_x_min_new)),
     y = y_pos_new + 25000,
     label = '150km',
     adj = c(0,0),
     cex = .6,
     srt = 45,
     xpd = TRUE)

text(x = tmp_df_x_max_new,
     y = y_pos_new + 25000,
     label = '300km',
     adj = c(0,0),
     cex = .6,
     srt = 45,
     xpd = TRUE)


dev.off()
