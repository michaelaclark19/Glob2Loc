#####
# HABITAT AVAIL AND POP ABUNDANCE TEMPLATE IN 2050
#
#####

###
# Rasters
# For trials
pop_abund <- rast('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Species_Analyses_New/pop_density_Phasmahyla_guttata_prevalence_2050_esh.future.extent.esh.exp.int.urb.tif')
hab_avail <- rast('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Species_Analyses_New/hab_availability_Phasmahyla_guttata_prevalence_2050_esh.future.extent.esh.exp.int.urb.tif')

# (1) hab avail, 2050
hab_avail <- 
  rast(outcome_rasters_prevalence %>%
         .[grepl(ss,.)] %>%
         .[grepl('2050',.)] %>%
         .[grepl('hab_avail.*.int.urb.tif',.)])
hab_avail[hab_avail < 0] <- 0

# (2) pop abund, 2050
pop_abund <- 
  rast(outcome_rasters_prevalence %>%
         .[grepl(ss,.)] %>%
         .[grepl('2050',.)] %>%
         .[grepl('pop_dens.*.int.urb.tif',.)])
pop_abund[pop_abund < 0] <- 0


### 
# Making plot template
layout.matrix <- matrix(c(rep(1,8),
                          rep(c(2,2,2,2,3,3,3,3),4)),
                        byrow = TRUE,
                        nrow = 5)
layout(mat = layout.matrix)
layout.show(3)
par(mar = c(0,0,0,0))

# (0) Title for section
plot.new()
text(x=.5,y=.5,'Projected Habitat and Abundance in 2050',font = 2,cex = 2)

# (1) Hab avail in 2050
plot(hab_avail,
     axes = FALSE, 
     col = viridis::viridis(100),
     range = c(0,1),
     main = 'Habitat Availability',
     plg = list(title = NULL, x = 'bottom', title.adj = .5))
# ,
# col = (viridis::mako(99)))
# col = rev(colorRampPalette(brewer.pal(n = 9, 'YlGnBu'))(99)))
plot(my_sf_reprojected,
     add = TRUE,
     lwd = .5)

text(x = mean(ext(hab_avail)[1:2]),
     y = min(ext(hab_avail)[3:4]) - (ext(hab_avail)[4] - ext(hab_avail)[3])*.05,
     'Proportion of Cell',
     adj = .5, 
     xpd = TRUE)



plot(pop_abund,
     axes = FALSE, 
     col = viridis::viridis(100),
     main = 'Population Abundance',
     plg = list(title = NULL, x = 'bottom'))

plot(my_sf_reprojected,
     add = TRUE,
     lwd = .5)

text(x = mean(ext(pop_abund)[1:2]),
     y = min(ext(pop_abund)[3:4]) - (ext(pop_abund)[4] - ext(pop_abund)[3])*.05,
     'Number of Individuals per sq km',
     adj = .5, 
     xpd = TRUE)
