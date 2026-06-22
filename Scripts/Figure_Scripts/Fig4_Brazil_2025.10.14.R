#####
# Fig 3
# Landscape analysis for Brazil
# Panels:
# (a) change of available habitat area in Brazil
# (b) bar charts showing overlap / no overlap between stressors
# (c) losses avoided by biodiversity hotspot
#####

###
# Working directory
setwd('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity')

###
# Libraries
library(terra)
library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(readr)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(sf)

###
# data management

###
# Panel (a) - total change habitat available 

# Habitat 2020
tot_hab_2020 <-
  list.files(path = paste0(getwd(),'/Files_For_Fig4'),
         pattern = 'hab_avail',
         full.names = TRUE) %>%
  .[grepl('_2020',.)] %>%
  .[!grepl('richness',.)] %>%
  rast(.) %>%
  app(.,sum,na.rm=TRUE)

# Habitat 2050
tot_hab_2050 <-
  list.files(path = paste0(getwd(),'/Files_For_Fig4'),
             pattern = 'hab_avail',
             full.names = TRUE) %>%
  .[grepl('_2050',.)] %>%
  .[!grepl('richness',.)] %>%
  .[grepl('.int.urb',.)] %>% # Limiting to all stressor scenario
  rast(.) %>%
  app(.,sum,na.rm=TRUE)

# Richness 2020
tot_richness_2020 <-
  list.files(path = paste0(getwd(),'/Files_For_Fig4/'),
             pattern = 'hab_avail',
             full.names = TRUE) %>%
  .[grepl('_2020',.)] %>%
  .[grepl('richness',.)] %>% 
  rast(.) %>%
  app(.,sum,na.rm=TRUE)

# Richness 2050
tot_richness_2050 <-
  list.files(path = paste0(getwd(),'/Files_For_Fig4/'),
             pattern = 'hab_avail',
             full.names = TRUE) %>%
  .[grepl('_2020',.)] %>%
  .[grepl('richness',.)] %>%
  .[grepl('.int.urb',.)] %>% # Limiting to all stressor scenario
  rast(.) %>%
  app(.,sum,na.rm=TRUE)

# Map --- difference in total habitat
dif_tot_hab <- tot_hab_2050 - tot_hab_2020

country_id <- rast(paste0(getwd(),'/Global Mollweide Maps/MollweideCountryID_1.5km.tif'))
country_id <- crop(country_id, dif_tot_hab)
country_id[country_id != 76] <- NA

dif_tot_hab[is.na(country_id)] <- NA

# Palette for fig
# my_pal <-
  # c(colorRampPalette(rev(brewer.pal('Reds',n=9)))(round((99 * (5/7)))),
  #   colorRampPalette(brewer.pal('Blues',n=9))(round((99 * (2/7)))))

my_pal <-
  # Palette for consistency with other figs
  # c(colorRampPalette(c('#3581b9','#99e9fc','#cbf1f7','#fefefe'))(round(99 * 2/7)),
  #   colorRampPalette(c('#fefefe','#f2bf91','#f28379','#d95252'))(round(99*3/7)))
  
  # Palette to make the cols pop out
  c(colorRampPalette(c(brewer.pal('GnBu',n=6)[6:3],'#fefefe'))(round((99 * (2/7)))),
    colorRampPalette(c('#fefefe','#f8eba3','#f0b082','#c94a6b','#aa00aa'))(round((99 * 3/7))))
                       
#   c(brewer.pal('GnBu',n=6)[6:3],
#     '#d6d6d6',
#     '#F8EBA3',
#     # '#D7E8EF',
#     '#F0B082',
#     '#C94A6B',
#     '#aa00aa')




###
# Setting plot margins
# par(mar = c(0,0,0,0))

###
# Layout for fig
# matrix_layout <-
#   matrix(c(1,1,2,
#            1,1,3,
#            4,4,4),
#          nrow = 3, byrow = TRUE)
# 
# 
# layout(matrix_layout,
#        heights = c(1,1,1.7),
#        widths = c(1,1,1))
# layout.show(4)

###
# Getting shapefile for brazilian states
brazil_states <- ne_states(country = 'Brazil', returnclass = 'sv')
brazil_states <- project(brazil_states, dif_tot_hab)

# Getting shapefile for cerrado
ecoregions_vect <- vect(paste0(getwd(),'/Ecoregions_Feb2023/Terrestrial_Ecoregions/Terrestrial_Ecoregions.shp'))
cerrado_vect <- ecoregions_vect[ecoregions_vect$ECO_ID_U %in% 10438]
cerrado_vect <- project(cerrado_vect, dif_tot_hab)

### Plotting
pdf(paste0(getwd(),'/Figures and Tables/Figures/Fig4_a_Brazil_',Sys.Date(),'.pdf'),
    width = 18.2 / 2.54 * 2/3,
    height = 5)

par(mar = c(0,0,0,0))

dif_tot_hab[dif_tot_hab <=-300] <- -300
plot(dif_tot_hab,
     box = FALSE,
     # legend = 'bottom',
     axes = FALSE,
     range = c(-300,200),
     col = rev(my_pal),
     mar = c(0,0,0,0),
     legend = FALSE)

### Adding shape files

plot(cerrado_vect, add = TRUE, border = '#333333', lty = 5)
plot(brazil_states, add = TRUE, border = '#aaaaaa',background = 'white')



###
# Adding mask raster
brazil_mask <- dif_tot_hab
brazil_mask[is.na(brazil_mask)] <- -999
brazil_mask[brazil_mask > -999] <- NA

plot(brazil_mask,
     box = FALSE,
     # legend = 'bottom',
     axes = FALSE,
     # range = c(-300,200),
     col = 'white',
     mar = c(0,0,0,0),
     legend = FALSE,
     add = TRUE)

# text(x = ext(dif_tot_hab)[1] + 
#        (ext(dif_tot_hab)[2] - ext(dif_tot_hab)[1]) * .025, 
#      y = ext(dif_tot_hab)[4] - 
#        (ext(dif_tot_hab)[4] - ext(dif_tot_hab)[3]) * .025, 
#      label = '(a)',
#      adj = c(0,1),
#      cex = .6)
text(x = par('usr')[1],
     y = par('usr')[4],
     label = '(a)',
     adj = c(0,1),
     cex = .6,
     xpd = TRUE)

# Adding manual legend
# Define the color gradient
n_colors <- length(my_pal)  # Number of gradient steps
gradient_colors <- rev(my_pal)

# Define rectangle bounds
xleft <- .2
xright <- .23
ybottom <- .1
ytop <- .4

# Convert to user coordinates
xleft <- grconvertX(xleft, from = "ndc", to = "user")
xrightlabel <- grconvertX(xright+.025, from = "ndc", to = "user")
xright <- grconvertX(xright, from = "ndc", to = "user")
ybottom <- grconvertY(ybottom, from = "ndc", to = "user")
ytitlelabel <- grconvertY(ytop+.05, from = "ndc", to = "user")
ytop <- grconvertY(ytop, from = "ndc", to = "user")

# Compute width of each sub-rectangle
y_steps <- seq(ybottom, ytop, length.out = n_colors + 1)

# Draw gradient using thin vertical rectangles
for (i in 1:n_colors) {
  rect(xleft = xleft, ybottom = y_steps[i],
       xright = xright, ytop = y_steps[i+1],
       col = gradient_colors[i], border = NA,
       xpd = TRUE)
}

# Optional: draw border around the whole area
rect(xleft, ybottom, xright, ytop, border = "black", xpd = TRUE)

text(x = xrightlabel, 
     y = ytop, 
     label = '200',
     xpd = TRUE,
     adj = c(0,.5),
     cex = .6)

text(x = xrightlabel, 
     y = ybottom, 
     label = '-300',
     xpd = TRUE,
     adj = c(0,.5),
     cex = .6)

text(x = xrightlabel, 
     y = ybottom*.4+ytop*.6, 
     label = '0',
     xpd = TRUE,
     adj = c(0,.5),
     cex = .6)

text(x = mean(xleft,xright),
     y = ytitlelabel,
     label = 'Change in Total\nHabitat Availability',
     adj = c(0.5,.5),
     cex = .6,
     xpd = TRUE)




dev.off()

### with ggplot
# plot_df <- as.data.frame(dif_tot_hab, xy = TRUE)
# 
# full_brazil_plot <-
#   ggplot(plot_df, aes(x = x, y = y, fill = sum)) +
#   geom_raster() +
#   theme_classic() +
#   theme(axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.line = element_blank()) +
#   theme(legend.position = c(.2,.2)) +
#   labs(x = '', y = '', fill = 'Change in Total\nHabitat Area') +
#   scale_fill_gradientn(colors = my_pal,
#                        limits = c(-500,200))

###
# Panel b --- overlap / no overlap between stressors
# Abundance in 2020
abundance_2020 <-
  list.files(path = paste0(getwd(),'/Files_For_Fig4/'), full.names = TRUE) %>%
  .[grepl(pattern = '2020',.)] %>%
  .[grepl(pattern = 'pop_dens',.)] %>%
  .[!grepl(pattern = 'richness',.)] %>%
  rast(.) %>%
  app(sum,na.rm=TRUE)

# Plotting differences by scenario
scen_list <- 
  c('current',
    'future.extent.esh.exp.tif',
    'future.extent.esh.urb.tif',
    'future.extent.esh.int.tif')


###
country_id <- rast(paste0(getwd(),'/Global Mollweide Maps/MollweideCountryID_1.5km.tif'))
country_id <- crop(country_id, abundance_2020)
country_id[country_id != 76] <- NA

col_list <- 
  c('#4169E1',
    '#DC143C',
    '#FFD700',
    '#008080')

# dev.off()
# par(mfrow = c(1,1))
for(ss in scen_list) {
  tmp_hab <-
    list.files(path = paste0(getwd(),'/Files_For_Fig4/'), full.names = TRUE) %>%
    .[grepl(pattern = '2050',.)] %>%
    .[grepl(pattern = 'pop_',.)] %>%
    .[!grepl(pattern = 'richness',.)] %>%
    .[grepl(ss,.)] %>%
    rast(.) %>%
    app(sum,na.rm=TRUE)
  
  tmp_hab[is.na(country_id)] <- NA
  
  tmp_plot <- tmp_hab / abundance_2020
  tmp_quant <- quantile(values(tmp_plot),.1,na.rm=TRUE)
  if(tmp_quant >= 1) {tmp_quant <- .999}
  tmp_plot_quant <- tmp_plot
  tmp_plot_quant[tmp_plot_quant > tmp_quant] <- NA
  tmp_plot_quant[!is.na(tmp_plot_quant)] <- 1
  
  # plot(tmp_plot)
  # plot(tmp_plot_quant, 
  #      main = paste0(ss,': ',round(tmp_quant,digits=2)),
  #      col = col_list[which(scen_list %in% ss)],
  #      alpha = .5)
  
  if(ss %in% scen_list[1]) {
    out_raster <- tmp_plot
    out_raster[out_raster <= tmp_quant] <- -1
    out_raster[out_raster >= 0] <- 0
    out_raster_out_stressor <- out_raster
    
    out_raster_pos <- tmp_plot
    out_raster_pos[out_raster_pos <= 1.001] <- 0
    out_raster_pos[out_raster_pos > 1.001] <- -1
    out_raster_pos_stressor <- out_raster_pos
    
    # plot(tmp_plot_quant, 
    #      # main = paste0(ss,': ',round(tmp_quant,digits=2)),
    #      col = col_list[which(scen_list %in% ss)],
    #      alpha = .5,
    #      main = 'High Risk Areas By Stressor',
    #      box = FALSE,
    #      axes = FALSE,
    #      legend = FALSE)
  } else {
    out_raster <- tmp_plot
    out_raster[out_raster <= tmp_quant] <- -1
    out_raster[out_raster >= 0] <- 0
    out_raster_out_stressor <- out_raster + out_raster_out_stressor
    
    out_raster_pos <- tmp_plot
    out_raster_pos[out_raster_pos <= 1.001] <- 0
    out_raster_pos[out_raster_pos > 1.001] <- -1
    out_raster_pos_stressor <- out_raster_pos + out_raster_pos_stressor
    
    # plot(tmp_plot_quant, 
    #      # main = paste0(ss,': ',round(tmp_quant,digits=2)),
    #      col = col_list[which(scen_list %in% ss)],
    #      alpha = .5,
    #      add = TRUE)
  }
}

# dev.off()

my_pal <-
  c(brewer.pal('GnBu',n=6)[6:3],
    '#fefefe',
    '#F8EBA3',
    # '#D7E8EF',
    '#F0B082',
    '#C94A6B',
    '#aa00aa')

###
# 
# plot(out_raster_pos_stressor - out_raster_pos_stressor,
#      col = rev(my_pal),
#      type = 'continuous',
#      range = c(-4,4),
#      alpha = 1,
#      legend = 'bottom',
#      box = FALSE,
#      axes = FALSE,
#      plg = list(cex = .5, width = .1))

pdf(paste0(getwd(),'/Figures and Tables/Figures/Fig4_c_stressors_Brazil_',Sys.Date(),'.pdf'),
    width = 18.2 / 2.54 * 1/3,
    height = 2.5)

par(mar = c(0,0,0,0))

plot(-out_raster_out_stressor,
     col = c('#fefefe',
             '#F8EBA3',
             # '#D7E8EF',
             '#F0B082',
             '#C94A6B'),
     box = FALSE,
     axes = FALSE,
     legend = FALSE,
     mar = c(0,0,0,0))#,
     # legend = 'bottomleft')

out_raster_pos_stressor[out_raster_pos_stressor %in% 0] <- NA
plot(-out_raster_pos_stressor,
     col = brewer.pal('GnBu',n=6)[3:6],
     add = TRUE,
     box = FALSE,
     axes = FALSE,
     alpha = .75,
     legend = FALSE,
     mar = c(0,0,0,0))#,
     # legend = 'bottomright')


plot(brazil_states, add = TRUE, border = '#aaaaaa')

# legend("bottomleft", legend = 4:0, fill = (my_pal), 
#        bty = "n", cex = 0.6, inset = c(0.1, 0.01),
#        title = 'Number of Stressors\nWith Positive Impacts',
#        y.intersp = .8)
# 
# legend("bottomleft", legend = 4:0, fill = (my_pal)[9:5], 
#        bty = "n", cex = 0.6, inset = c(0.3, 0.01),
#        title = 'Number of Stressors\nWith Negative Impacts',
#        y.intersp = .8)

lgd <-
  legend("bottomleft", legend = c(4:0,4:0), fill = c(my_pal[1:5],my_pal[9:5]), ncol = 2,
         bty = "n", cex = 0.6, inset = c(0.01, 0.01),
         # title = 'Number of Taxa With:\nAny Positive Impact (Blues)\nLarge Negative Impact (Reds)',
         # adj = 0,
         y.intersp = .8,
         xpd = TRUE)

text(x = lgd$rect$left,   # left edge of legend box
     y = lgd$rect$top + 0.5,  # slightly above the legend
     labels = 'Number Stressors With:\nAny Benefit (Blue)\nLarge Negative Impact (Red)',
     adj = c(0,0),  # left-justify text
     cex = .6,
     xpd = TRUE)

# text(x = ext(out_raster_out_stressor)[1] + 
#        (ext(out_raster_out_stressor)[2] - ext(out_raster_out_stressor)[1]) * .025, 
#      y = ext(out_raster_out_stressor)[4] - 
#        (ext(out_raster_out_stressor)[4] - ext(out_raster_out_stressor)[3]) * .025, 
#      label = '(b)',
#      adj = c(0,1),
#      cex = .6)

text(x = par('usr')[1],
     y = par('usr')[4],
     label = '(c)',
     adj = c(0,1),
     cex = .6,
     xpd = TRUE)

dev.off()


###
# Plotting
# df_pos <- as.data.frame(-out_raster_pos_stressor, xy = TRUE)
# df_neg <- as.data.frame(out_raster_out_stressor, xy = TRUE)
# 
# stressor_plot <-
#   ggplot() +
#   geom_raster(data = df_neg, aes(x = x, y = y, fill = sum)) +
#   geom_raster(data = df_pos, aes(x = x, y = y, fill = sum)) +
#   theme_classic() +
#   theme(axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.line = element_blank()) +
#   scale_fill_gradientn(colors = rev(my_pal),
#                        limits = c(-4,4),
#                        guide = guide_colourbar(direction = 'horizontal',
#                                                title.position = 'top',
#                                                title.hjust = .5),
#                        name = 'Number of Stressors') +
#   theme(legend.position = c(.2,.2)) +
#   labs(x = '', y = '')

###
# Repeating by taxa
# Plotting differences by scenario
scen_list <- 
  c('Amphibians',
    'Birds',
    'Mammals',
    'Reptiles')


###
country_id <- rast(paste0(getwd(),'/Global Mollweide Maps/MollweideCountryID_1.5km.tif'))
country_id <- crop(country_id, abundance_2020)
country_id[country_id != 76] <- NA

# dev.off()
# par(mfrow = c(1,1))
for(ss in scen_list) {
  tmp_hab <-
    list.files(path = paste0(getwd(),'/Files_For_Fig4/'), full.names = TRUE) %>%
    .[grepl(pattern = '2050',.)] %>%
    .[grepl(pattern = 'pop_',.)] %>%
    .[!grepl(pattern = 'richness',.)] %>%
    .[grepl(ss,.)] %>%
    .[grepl('.int.urb',.)] %>%
    rast(.) #%>%
    # app(sum,na.rm=TRUE)
  
  hab_2020 <-
    list.files(path = paste0(getwd(),'/Files_For_Fig4/'), full.names = TRUE) %>%
    .[grepl(pattern = '2020',.)] %>%
    .[grepl(pattern = 'pop_',.)] %>%
    .[!grepl(pattern = 'richness',.)] %>%
    .[grepl(ss,.)] %>%
    .[grepl('.int.urb',.)] %>%
    rast(.) #%>%
  
  tmp_hab[is.na(country_id)] <- NA
  hab_2020[is.na(country_id)] <- NA
  
  tmp_plot <- tmp_hab / hab_2020
  tmp_quant <- quantile(values(tmp_plot),.1,na.rm=TRUE)
  if(tmp_quant >= 1) {tmp_quant <- .999}
  
  tmp_plot_quant <- tmp_plot
  tmp_plot_quant[tmp_plot_quant > tmp_quant] <- NA
  tmp_plot_quant[!is.na(tmp_plot_quant)] <- 1
  
  if(ss %in% scen_list[1]) {
    out_raster <- tmp_plot
    out_raster[out_raster <= tmp_quant] <- -1
    out_raster[out_raster >= 0] <- 0
    out_raster_out_taxa <- out_raster
    
    out_raster_pos <- tmp_plot
    out_raster_pos[out_raster_pos <= 1.001] <- 0
    out_raster_pos[out_raster_pos > 1.001] <- -1
    out_raster_pos_taxa <- out_raster_pos
    
    # plot(tmp_plot_quant, 
    #      # main = paste0(ss,': ',round(tmp_quant,digits=2)),
    #      col = col_list[which(scen_list %in% ss)],
    #      alpha = .5,
    #      main = 'High Risk Areas By Taxa',
    #      box = FALSE,
    #      axes = FALSE)
  } else {
    out_raster <- tmp_plot
    out_raster[out_raster <= tmp_quant] <- -1
    out_raster[out_raster >= 0] <- 0
    out_raster_out_taxa <- out_raster + out_raster_out_taxa
    
    out_raster_pos <- tmp_plot
    out_raster_pos[out_raster_pos <= 1.001] <- 0
    out_raster_pos[out_raster_pos > 1.001] <- -1
    out_raster_pos_taxa <- out_raster_pos + out_raster_pos_taxa
    # plot(tmp_plot_quant, 
    #      # main = paste0(ss,': ',round(tmp_quant,digits=2)),
    #      col = col_list[which(scen_list %in% ss)],
    #      alpha = .5,
    #      add = TRUE)
  }
}

# dev.off()

###
# 
pdf(paste0(getwd(),'/Figures and Tables/Figures/Fig4_b_taxa_Brazil_',Sys.Date(),'.pdf'),
    width = 18.2 / 2.54 * 1/3,
    height = 2.5)

par(mar = c(0,0,0,0))

plot(-out_raster_out_taxa,
     col = c('#fefefe',
             '#F8EBA3',
             # '#D7E8EF',
             '#F0B082',
             '#C94A6B',
             '#9B2171'),
     box = FALSE,
     axes = FALSE,
     legend = FALSE,
     mar = c(0,0,0,0))

out_raster_pos_taxa[out_raster_pos_taxa %in% 0] <- NA
plot(-out_raster_pos_taxa,
     col = brewer.pal('GnBu',n=6)[3:6],
     add = TRUE,
     alpha = .75,
     legend = FALSE,
     mar = c(0,0,0,0))
# Adding lines for brazil states
plot(brazil_states, add = TRUE, border = '#aaaaaa')
# Adding legend
# legend("bottomleft", legend = 4:0, fill = (my_pal), 
#        bty = "n", cex = 0.6, inset = c(0.05, 0.01),
#        title = 'Number of Taxa Experiencing\nPositive Impacts (Blues)\nNegative Impacts (Reds)',
#        y.intersp = .8)
# 
# legend("bottomright", legend = 4:0, fill = (my_pal)[9:5], 
#        bty = "n", cex = 0.6, inset = c(0.05, 0.01),
#        title = 'Number of Taxa\nWith Negative Impacts',
#        y.intersp = .8)


lgd <-
  legend("bottomleft", legend = c(4:0,4:0), fill = c(my_pal[1:5],my_pal[9:5]), ncol = 2,
       bty = "n", cex = 0.6, inset = c(0.01, 0.01),
       # title = 'Number of Taxa With:\nAny Positive Impact (Blues)\nLarge Negative Impact (Reds)',
       # adj = 0,
       y.intersp = .8,
       xpd = TRUE)

text(x = lgd$rect$left,   # left edge of legend box
     y = lgd$rect$top + 0.5,  # slightly above the legend
     labels = 'Number Taxa With:\nAny Benefit (Blue)\nLarge Negative Impact (Red)',
     adj = c(0,0),  # left-justify text
     cex = .6,
     xpd = TRUE)

# text(x = ext(out_raster_pos_taxa)[1] + 
#        (ext(out_raster_pos_taxa)[2] - ext(out_raster_pos_taxa)[1]) * .025, 
#      y = ext(out_raster_pos_taxa)[4] - 
#        (ext(out_raster_pos_taxa)[4] - ext(out_raster_pos_taxa)[3]) * .025, 
#      label = '(c)',
#      adj = c(0,1),
#      cex = .6)

text(x = par('usr')[1],
     y = par('usr')[4],
     label = '(b)',
     adj = c(0,1),
     cex = .6,
     xpd = TRUE)

dev.off()



###
# Panel d --- avoided impacts by biodiversity hotspot


###
# Getting list of unique species
# list of files
file_list <-
  list.files(paste0(getwd(),'/Aggregated_CSV_Files/'),
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

table(out_df_unique$taxon)

# checking it worked -- good
head(out_df)

###
# which species have estimated area in 2020
species_list <-
  out_df %>%
  filter(year %in% 2020) %>%
  dplyr::select(taxon, species, tot_pop, tot_area) %>%
  unique() %>%
  mutate(tot_pop = as.numeric(tot_pop),
         tot_area = as.numeric(tot_area)) %>%
  # filter(!is.na(tot_pop) &
  #        !is.na(tot_area)) %>%
  filter(!is.na(tot_area)) %>%
  filter(tot_area > 0)

###
# Getting outcomes for bau
file_list <-
  list.files(paste0(getwd(),'/Aggregated_CSV_Files/'),
             pattern = 'Prop',
             full.names = TRUE) %>%
  .[grepl('prevalence',.)] %>%
  .[grepl('nopatch',.)] %>%
  .[!grepl('Climate_Suitable',.)]

# importing files
out_df_bau <-
  do.call(rbind, lapply(file_list, read_csv)) %>%
  mutate(taxon = gsub('/.*','',binomial)) %>%
  mutate(species = gsub('.*/','',binomial)) %>%
  filter(paste0(taxon,species) %in% paste0(species_list$taxon,species_list$species)) %>%
  filter(year %in% 2050) %>%
  filter(grepl('.int.urb',scenario)) %>%
  dplyr::select(taxon,species,binomial,
                tot_pop_bau = tot_pop,
                tot_area_bau = tot_area,
                num_patches_bau = num_patches)

###
# Getting outcomes for eat lancet
file_list <-
  list.files(paste0(getwd(),'/Aggregated_CSV_Files_EAT_Lancet/'),
             pattern = 'Prop',
             full.names = TRUE) %>%
  .[grepl('prevalence',.)] %>%
  .[grepl('nopatch',.)] %>%
  .[!grepl('Climate_Suitable',.)]

# importing files
out_df_eat <-
  do.call(rbind, lapply(file_list, read_csv)) %>%
  mutate(taxon = gsub('/.*','',binomial)) %>%
  mutate(species = gsub('.*/','',binomial)) %>%
  filter(paste0(taxon,species) %in% paste0(species_list$taxon,species_list$species)) %>%
  filter(year %in% 2050) %>%
  filter(grepl('.int.urb',scenario)) %>%
  dplyr::select(taxon,species,binomial,
                tot_pop_eat = tot_pop,
                tot_area_eat = tot_area,
                num_patches_eat = num_patches)

###
# Getting relative outcomes by biodiv hotspot
###
# For specices endemic to the atlantic forest
biodiv_hotspots <-
  read.csv(paste0(getwd(),'/Analyses/species_by_biodiv_hotspot_esh_maps.csv')) %>%
  filter(ecoregion_id > 0) %>%
  dplyr::group_by(taxa,species) %>%
  dplyr::summarise(prop_cells = sum(prop_cells,na.rm = TRUE)) %>%
  filter(prop_cells %in% 1) %>%
  mutate(binomial = paste0(taxa,'/',species)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(taxa,species))

hotspots_vect <- vect(paste0(getwd(),'/Ecoregions_Feb2023/hotspots_2016_1/hotspots_2016_1.shp'))
ecoregions_vect <- vect(paste0(getwd(),'/Ecoregions_Feb2023/Terrestrial_Ecoregions/Terrestrial_Ecoregions.shp'))
hotspots_map <- raster(paste0(getwd(),'/Ecoregions_Feb2023/Biodiversity_Hotspots_Raster.tif'))
realms_map <- raster(paste0(getwd(),'/Ecoregions_Feb2023/realms_raster.tif'))

hotspots_realms_df <-
  data.frame(hotspot = getValues(hotspots_map),
             realms = getValues(realms_map)) %>%
  filter(!is.na(hotspot)) %>%
  filter(!is.na(realms)) %>%
  distinct() %>%
  left_join(.,
            data.frame(realm_name = ecoregions_vect$WWF_REALM2) %>%
              distinct() %>%
              mutate(realms = 
                       case_when(grepl('Neotropic',realm_name) ~ 1,
                                 grepl('Palearctic',realm_name) ~ 2,
                                 grepl('Nearctic',realm_name) ~ 3,
                                 grepl('Indo-Malay',realm_name) ~ 5,
                                 grepl('Afrotropic',realm_name) ~ 6,
                                 grepl('Oceania',realm_name) ~ 7,
                                 grepl('Austral',realm_name) ~ 8,
                                 grepl('Antarctic',realm_name) ~ 9))) %>%
  mutate(hotspots_name = 
           case_when(hotspot %in% 1 ~ 'Atlantic Forest',
                     hotspot %in% 2 ~ 'California Floristic Province',
                     hotspot %in% 3 ~ 'California Floristic Province (Outer Limits)',
                     hotspot %in% 4 ~ 'Cape Floristic Region',
                     hotspot %in% c(5) ~ 'Caribbean Islands',
                     hotspot %in% c(6) ~ 'Caribbean Islands (Outer Limits',
                     hotspot %in% 7 ~ 'Caucasus',
                     hotspot %in% 8 ~ 'Cerrado',
                     hotspot %in% 9 ~ 'Chilean Winter Rainfall and Valdivian Forests',
                     hotspot %in% 10 ~ 'Chilean Winter Rainfall and Valdivian Forests (Outer Limits)',
                     hotspot %in% 11 ~ 'Coastal Forests of Eastern Africa',
                     hotspot %in% 12 ~ 'East Melanesian Islands',
                     hotspot %in% 13 ~ 'East Melanesian Islands (Outer Limits)',
                     hotspot %in% 14 ~ 'Eastern Afromontane',
                     hotspot %in% 15 ~ 'Guinean Forests of West Africa',
                     hotspot %in% 16 ~ 'Himalaya',
                     hotspot %in% 17 ~ 'Horn of Africa',
                     hotspot %in% 18 ~ 'Horn of Africa (Outer Limits)',
                     hotspot %in% 19 ~ 'Indo-Burma',
                     hotspot %in% 20 ~ 'Indo-Burma (Outer Limits)',
                     hotspot %in% 21 ~ 'Irano-Anatolian',
                     hotspot %in% 22 ~ 'Japan',
                     hotspot %in% 23 ~ 'Japan (Outer Limits)',
                     hotspot %in% 24 ~ 'Madagascar +',
                     hotspot %in% 25 ~ 'Madagascar +',
                     hotspot %in% 26 ~ 'Madrean Pine-Oak Woodlands',
                     hotspot %in% 27 ~ 'Maputaland-Pondoland-Albany',
                     hotspot %in% 28 ~ 'Mediterranean Basin',
                     hotspot %in% 29 ~ 'Mediterranean Basin (Outer Limits)',
                     hotspot %in% 30 ~ 'Mesoamerica',
                     hotspot %in% 31 ~ 'Mesoamerica (Outer Limits)',
                     hotspot %in% 32 ~ 'Mountains of Central Asia',
                     hotspot %in% 33 ~ 'Mountains of Southwest China',
                     hotspot %in% 34 ~ 'New Caledonia',
                     hotspot %in% 35 ~ 'New Caledonia (Outer Limits)',
                     hotspot %in% 36 ~ 'New Zealand',
                     hotspot %in% 37 ~ 'New Zealand (Outer Limits)',
                     hotspot %in% 38 ~ 'Philippines',
                     hotspot %in% 39 ~ 'Philippines (Outer Limit)',
                     hotspot %in% 40 ~ 'Polynesia-Micronesia',
                     hotspot %in% 41 ~ 'Polynesia-Micronesia (Outer Limit)',
                     hotspot %in% 42 ~ 'Southwest Australia',
                     hotspot %in% 43 ~ 'Succulent Karoo',
                     hotspot %in% 44 ~ 'Sundaland',
                     hotspot %in% 45 ~ 'Sundaland (Outer Limits)',
                     hotspot %in% 46 ~ 'Tropical Andes',
                     hotspot %in% 47 ~ 'Tumbes-Choco-Magdalena',
                     hotspot %in% 48 ~ 'Tumbes-Choco-Magdalena (Outer Limits)',
                     hotspot %in% 49 ~ 'Wallacea',
                     hotspot %in% 50 ~ 'Wallacea (Outer Limits)',
                     hotspot %in% 51 ~ 'Western Ghats and Sri Lanka',
                     hotspot %in% 52 ~ 'Forests of East Australia',
                     hotspot %in% 53 ~ 'North American Coastal Plain')) %>%
  filter(grepl('Atlantic Forest|Cerrado',hotspots_name))

# Getting relative outcomes
out_df_merged <- 
  left_join(out_df_bau,
            out_df_eat) %>%
  mutate(loss_pop_bau = log2(tot_pop_bau + 1e-9), # Calculating loss in each outcome by scenario
         loss_area_bau = log2(tot_area_bau + 1e-9),
         frag_bau = log2(num_patches_bau + 1e-9),
         loss_pop_eat = log2(tot_pop_eat + 1e-9),
         loss_area_eat = log2(tot_area_eat + 1e-9),
         frag_eat = log2(num_patches_eat + 1e-9)) %>%
  mutate(relative_avoided_loss_pop = (loss_pop_eat - loss_pop_bau) / loss_pop_bau,
         relative_avoided_loss_area = (loss_area_eat - loss_area_bau) / loss_area_bau,
         relative_avoided_loss_patches = (frag_eat - frag_bau) / frag_bau) %>%
  # mutate(relative_avoided_loss_pop = log2(tot_pop_eat / tot_pop_bau),
  #        relative_avoided_loss_area = log2(tot_area_eat / tot_area_bau),
  #        relative_avoided_loss_patches = log2(num_patches_eat / num_patches_bau)) %>%
  filter(!is.na(tot_area_eat)) %>% # Removing species not endemic to biodiversity hotspots
  filter(binomial %in% biodiv_hotspots$binomial) %>%
  left_join(.,
            read.csv(paste0(getwd(),'/Ecoregions_Feb2023/species_by_ecoregion_esh_maps.csv')) %>%
              mutate(species = gsub('_[A-Z].*|_[0-9].*','',species)) %>%
              mutate(binomial = paste0(taxa,'/',species)) %>%
              dplyr::select(binomial, ecoregion_id, prop_cells)) %>%
  left_join(.,
            data.frame(ecoregion_id = ecoregions_vect$ECO_ID_U,
                       realm_name = ecoregions_vect$WWF_REALM2)) %>%
  dplyr::select(-c(ecoregion_id,prop_cells)) %>%
  distinct() %>%
  left_join(.,
            read.csv(paste0(getwd(),'/Analyses/species_by_biodiv_hotspot_esh_maps.csv')) %>%
              mutate(binomial = paste0(taxa,'/',species)) %>%
              dplyr::select(binomial, hotspot_id = ecoregion_id, prop_cells_hotspot = prop_cells) %>%
              filter(hotspot_id > -99) %>%
              distinct()) %>%
  left_join(.,
            hotspots_realms_df %>%
              dplyr::rename(hotspot_id = hotspot) %>%
              distinct()) %>%
  filter(grepl('Atlantic Forest|Cerrado', hotspots_name))

###
# How many species have net gains...
cat('Species with gains to population: ',nrow(out_df_eat %>% filter(tot_pop_eat > 1)))
cat('Species with gains to habitat area: ',nrow(out_df_eat %>% filter(tot_area_eat > 1)))
cat('Species with gains to patches: ',nrow(out_df_eat %>% filter(num_patches_eat >= 1)))


winners_ssp1 <-
  out_df_eat %>%
  filter(tot_pop_eat > 1 &
           tot_area_eat > 1)


# Plot data - need to order by hotspots
dat_plot = 
  rbind(out_df_merged %>%
          dplyr::select(binomial,hotspot_id,realm_name, biodiv_hotspot = hotspots_name, outcome = relative_avoided_loss_pop, outcome_bau = loss_pop_bau, outcome_eat = loss_pop_eat) %>%
          mutate(impact = 'Avoided Population Loss'),
        out_df_merged %>%
          dplyr::select(binomial,hotspot_id,realm_name, biodiv_hotspot = hotspots_name,outcome = relative_avoided_loss_area, outcome_bau = loss_area_bau, outcome_eat = loss_area_eat) %>%
          mutate(impact = 'Avoided Habitat Loss'),
        out_df_merged %>%
          dplyr::select(binomial,hotspot_id,realm_name,biodiv_hotspot = hotspots_name,outcome = relative_avoided_loss_patches, outcome_bau = frag_bau, outcome_eat = frag_eat) %>%
          mutate(impact = 'Avoided Habitat Fragmentation')) 

# Rbinding empty rows for placeholder between realms
dat_plot <-
  rbind(dat_plot,
        dat_plot[1,] %>% 
          mutate_all(~NA) %>%
          mutate(realm_name = 'Afrotropic',
                 biodiv_hotspot = 'ZZZ'),
        dat_plot[1,] %>% 
          mutate_all(~NA) %>%
          mutate(realm_name = 'Australasia',
                 biodiv_hotspot = 'ZZZ'),
        dat_plot[1,] %>% 
          mutate_all(~NA) %>%
          mutate(realm_name = 'Indo-Malay',
                 biodiv_hotspot = 'ZZZ'),
        dat_plot[1,] %>% 
          mutate_all(~NA) %>%
          mutate(realm_name = 'Nearctic',
                 biodiv_hotspot = 'ZZZ'),
        dat_plot[1,] %>% 
          mutate_all(~NA) %>%
          mutate(realm_name = 'Neotropic',
                 biodiv_hotspot = 'ZZZ'),
        # dat_plot[1,] %>% 
        #   mutate_all(~NA) %>%
        #   mutate(realm_name = 'Oceania',
        #          biodiv_hotspot = 'ZZZ'),
        dat_plot[1,] %>% 
          mutate_all(~NA) %>%
          mutate(realm_name = 'Palearctic',
                 biodiv_hotspot = 'ZZZ')) %>%
  arrange(realm_name, biodiv_hotspot) %>%
  mutate(realm_hotspot = paste0(realm_name, ': ', biodiv_hotspot)) %>%
  dplyr::group_by(realm_hotspot) %>%
  mutate(num_species = round(n() / 3),digits=0) %>%
  filter(num_species > 10 | grepl('ZZZ',realm_hotspot)) %>%
  mutate(realm_hotspot = paste0(realm_hotspot, ' (n = ',num_species,')')) %>%
  transform(realm_hotspot = 
              factor(realm_hotspot, 
                     levels = unique(.$realm_hotspot)))


###
# Getting order for axis
dat_order <-
  dat_plot %>%
  dplyr::group_by(realm_name, realm_hotspot,impact) %>%
  dplyr::summarise(median_outcome = median(outcome,na.rm=TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(realm_name, realm_hotspot) %>%
  dplyr::summarise(mean_outcome = mean(median_outcome,na.rm=TRUE)) %>%
  arrange(realm_name, mean_outcome)

###
# Adding labels for axis plots
plot_labels <-
  dat_order %>%
  dplyr::ungroup() %>%
  dplyr::select(realm_hotspot) %>%
  unique() %>%
  mutate(axis_labels = 
           case_when(grepl('ZZZ',realm_hotspot) ~ '',
                     .default = realm_hotspot)) %>%
  filter(!grepl('NA',realm_hotspot))


# Prepare data
filtered_data <- 
  subset(dat_plot %>%
           filter(!is.na(impact)), 
         !grepl("NA", realm_hotspot)) %>%
  filter(!is.na(outcome)) %>%
  filter(is.finite(outcome))

# Adjust outcome as in ggplot (multiply by -100)
filtered_data$outcome_adj <- 
  -filtered_data$outcome * 100

# Set colors manually for impacts
impact_colors <- c('#00cd6c', '#ffc61e', '#009ade')
names(impact_colors) <- unique(filtered_data$impact) 

# Create a combined group identifier for boxplot positions
filtered_data$group <- interaction(filtered_data$impact, filtered_data$realm_hotspot, drop = TRUE)

###
# Ordering for boxplots
dat_order_baseplot <-
  rbind(dat_order %>% mutate(impact = unique(filtered_data$impact)[1]),
        dat_order %>% mutate(impact = unique(filtered_data$impact)[2]),
        dat_order %>% mutate(impact = unique(filtered_data$impact)[3])) %>%
  mutate(group_levels = paste0(impact,'.',realm_hotspot)) %>%
  transform(impact = factor(impact,
                            levels = c('Avoided Population Loss','Avoided Habitat Loss','Avoided Habitat Fragmentation'))) %>%
  arrange(realm_name, mean_outcome, realm_hotspot, (impact)) %>%
  mutate(x_axis_loc = NA) %>%
  filter(!grepl('Pop.*Loss.*ZZZ|Habitat Loss.*ZZZ',group_levels)) %>%
  filter(!is.na(realm_name)) %>%
  filter(!grepl('NA',realm_hotspot)) %>%
  filter(!grepl('Palear.*ZZZ',realm_hotspot))

###
# Adding location
tmp_x <- 0
for(ii in 1:nrow(dat_order_baseplot)) {
  
  if(ii > 1) {
    if(dat_order_baseplot$realm_hotspot[ii] !=
       dat_order_baseplot$realm_hotspot[ii-1]) {# &
      # !grepl('ZZZ',dat_order_baseplot$realm_hotspot[ii])) {
      tmp_x <- tmp_x + 1
    } else {
      
    }
  }
  tmp_x <- tmp_x + 1
  if(grepl('ZZZ',dat_order_baseplot$realm_hotspot[ii])) {
    tmp_x <- tmp_x - 1
  }
  
  dat_order_baseplot$x_axis_loc[ii] <- tmp_x
}

# Create a numeric position for each unique group
group_levels <- dat_order_baseplot$group_levels
group_positions <- dat_order_baseplot$x_axis_loc

# Map group names to positions
group_map <- setNames(group_positions, group_levels)

# Calculate x-axis positions and labels for realm_hotspot
realms <- unique(filtered_data$realm_hotspot)
x_labels <- plot_labels$axis_labels
names(x_labels) <- plot_labels$realm_hotspot
realm_positions <- dat_order_baseplot$x_axis_loc

# Start plotting


# names(filtered_data)

filtered_data_plot <-
  filtered_data %>%
  mutate(taxon = gsub('/.*','',binomial)) %>%
  dplyr::select(-binomial) %>%
  ungroup() %>%
  filter(grepl('Atlantic',biodiv_hotspot)) %>%
  # dplyr::group_by(biodiv_hotspot, taxon) %>%
  arrange(biodiv_hotspot,impact, taxon, -outcome_eat) %>%
  # transform(binomial = factor(binomial, levels = unique(.$binomial))) %>%
  mutate(x = 1:nrow(.),
         xend = 1:nrow(.)) %>%
  mutate(alpha_segment = 
           case_when(outcome_bau < -4  ~ .5,
                     .default = 1)) %>%
  mutate(alpha_point = 
           case_when(outcome_eat < -4 ~ .5,
                     .default = 1)) %>%
  mutate(outcome_bau =
           case_when(outcome_bau < -4 ~ -4,
                     .default = outcome_bau)) %>%
  mutate(outcome_eat =
           case_when(outcome_eat < -4 ~ -4,
                     .default = outcome_eat))


# Filter data
data_plot <- subset(filtered_data_plot,
                    grepl("Atlantic", biodiv_hotspot) &
                      is.finite(outcome) &
                      !is.na(outcome)) %>%
  filter(grepl('Habitat Loss',impact)) %>%
  mutate(point = 
           case_when(alpha_point < 1 ~ 1,
                     .default = 19))

x_count <- 1
for(ii in 1:nrow(data_plot)) {
  
  
  if(ii > 1) {
    x_count <- x_count + 1
    if(data_plot$taxon[ii] != data_plot$taxon[ii - 1]) {
      x_count <- x_count + 2
    }
  }
  
  data_plot$x[ii] <- x_count
}

# Transform y values
y_bau <- data_plot$outcome_bau
y_eat <- data_plot$outcome_eat

# Add jitter to x positions (like ggplot's position = 'jitter')
set.seed(1)
# x_jittered <- jitter(data_plot$x, amount = 0.2)
# xend_jittered <- jitter(data_plot$xend, amount = 0.2)

# Define y-axis limits
ylim <- c(-4, 4)

# Assign colors by taxon
taxa <- unique(data_plot$taxon)
taxon_colors <- setNames(c('#009e73',
                           '#d55e00',
                           '#0072b2',
                           '#cc79a7'), 
                         taxa)

pdf(paste0(getwd(),'/Figures and Tables/Figures/Fig4_d_boxplots_Brazil_',Sys.Date(),'.pdf'),
    width = 18.2 / 2.54 * 3/3,
    height = 2.5)

par(mar = c(2,5,1,.5))
# Create empty plot
plot(NA, xlim = c(min(data_plot$x),max(data_plot$x)),
     ylim = c(-4,2.5), 
     xlab = "", 
     ylab = "Habitat Loss (2020 - 2050)", 
     axes = FALSE,
     cex.lab = .6)

# Horizontal reference line at y = 1
abline(h = 0, col = "black")

# Plot segments (from outcome_bau to outcome_eat)
for (i in seq_len(nrow(data_plot))) {
  lines(c(data_plot$x[i], data_plot$x[i]),
        c(data_plot$outcome_bau[i], data_plot$outcome_eat[i]),
        # col = taxon_colors[data_plot$taxon[i]],
        col = adjustcolor(taxon_colors[data_plot$taxon[i]], 
                          alpha.f = data_plot$alpha_segment[i]),
        lend = "square")
}

# Overlay open circles at outcome_eat
for(ii in 1:nrow(data_plot)) {
  points(data_plot$x[ii], 
         data_plot$outcome_eat[ii], 
         pch = 19,
         # pch = data_plot$point[ii],
         cex = 0.5, 
         # col = adjustcolor('#000000', 
         #   alpha.f = data_plot$alpha_point[ii]))
         col = adjustcolor(taxon_colors[data_plot$taxon[ii]], 
                           alpha.f = data_plot$alpha_point[ii]))
}


# Y-axis with log2 labels
axis(2, 
     at = c(-3:3), 
     labels = c("-87.5%", '-75%',"-50%", '0%',"+50%",'+75%', "+87.5%"), 
     las = 1,
     cex.axis = .6)

text_df <-
  data_plot %>%
  dplyr::group_by(taxon) %>%
  dplyr::summarise(x_loc = mean(x))

text(x = text_df$x_loc,
     y = 2,
     labels = text_df$taxon,
     adj = c(.5,0),
     cex = .6) 

###
# Adding separators by realm
for(ii in unique(data_plot$taxon)[1:3]) {
  xloc = max(data_plot$x[data_plot$taxon %in% ii]) + 1.5
  abline(v = xloc, col = 'grey')
}

# Box around plot
# rect(xleft = min(data_plot$x)-2,
#      xright = max(data_plot$x)+2,
#      ybottom = -4.5,
#      ytop = 4.5,
#      xpd = TRUE)
box()

x_fig <- grconvertX(0.00001, from = "ndc", to = "user")  # near left edge
y_fig <- grconvertY(0.99999, from = "ndc", to = "user")  # near top edge

text(#x = par('usr')[1],
  x = x_fig,
  #y = par('usr')[4],
  y = y_fig,
  label = '(d)',
  adj = c(0,1),
  cex = .6,
  xpd = TRUE)

grconvertX(0.00001, from = "ndc", to = "user")  # near left edge
y_fig <- grconvertY(0.99999, from = "ndc", to = "user")  # near top edge
# 
text(x = grconvertX(.25, from = "ndc", to = "user"),
     y = grconvertY(0.1, from = "ndc", to = "user"),
     label = 'Segment Start: BAU Scenario',
     adj = c(0,0),
     cex = .6,
     xpd = TRUE)
# 
text(x = grconvertX(.5, from = "ndc", to = "user"),
     y = grconvertY(0.1, from = "ndc", to = "user"),
     label = 'Circle: Mitigation Scenario',
     adj = c(0,0),
     cex = .6,
     xpd = TRUE)
#
text(x = grconvertX(.75, from = "ndc", to = "user"),
     y = grconvertY(0.1, from = "ndc", to = "user"),
     label = 'Shading: >95% Habitat Loss',
     adj = c(0,0),
     cex = .6,
     xpd = TRUE)
# text(x = data_plot$x[1]-1,
#      y = data_plot$outcome_bau[1],
#      label = 'Segment\nStart:\nBAU Outcome',
#      adj = c(1,0),
#      cex = .6,
#      xpd = TRUE)

dev.off()