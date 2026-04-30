###
# Libraries
library(fmsb)
library(plyr)
library(dplyr)
library(terra)
library(data.table)
library(readr)
library(RColorBrewer)
library(ggplot2)
library(scales)

###
# Setting working directory
setwd("/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity")

###
# Plotting raster map
# This is the big plot
winners_losers_rast <- rast('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Temp_Raster_Folders/ecoregions_winners_losers.tif')


levels(winners_losers_rast) <- 
  data.frame(value = 1:9, 
             name = c('<5% Winners; <5% Losers',
    '5-15% Winners; <5% Losers',
    '>15% Winners; <5% Losers',
    '<5% Winners; 5-15% Losers',
    '5-15% Winners; 5-15% Losers',
    '>15% Winners; 5-15% Losers',
    '<5% Winners; >15% Losers',
    '5-15% Winners; >15% Losers',
    '>15% Winners; >15% Losers'))

# Colours for raster
# cols <-
#   c('#DDDDDD', 
#     '#9BC8EB', 
#     '#64B9FA', 
#     '#EBDC82',
#     '#96BE87',
#     '#509B8C',
#     '#FAD74B', 
#     '#96AF3C', 
#     '#418232')

cols <-
  # c('#cdcdcd','#7f9cba','#396da6',
  #   '#b27c85','#6d6078','#33456d',
  #   '#932e3d','#5d2639','#2b1b33')


  # c('#cdcdcd','#85af72','#3c9120',
  #   '#c479a6','#7e695e','#3b551c',
  #   '#b41f72','#751441','#371215')
  
  # c('#cdcdcd','#8eb0ba','#5291a4',
  #   '#bd8b8b','#83777c','#4b6470',
  #   '#ae4848','#793e42','#45353a')
  # 
  # c('#3581b9','#99e9fc','#cbf1f7','#dddddd')
  # c('#dddddd','#f2bf91','#f28379','#d95252')
  
  c('#dddddd','#99e9fc','#3581b9',
    '#f2bf91','#BAAEB9','#66729A',
    '#d95252','#A76070','#876985')
  




###
# Getting country borders
country_borders <- vect('/Users/michael/Downloads/TM_WORLD_BORDERS-0.3/TM_WORLD_BORDERS-0.3.shp')
country_borders <- project(country_borders, winners_losers_rast)

pdf(paste0(getwd(),'/Figures and Tables/Figures/Fig3_a_winners_losers',Sys.Date(),'.pdf'),
    width = 12 / 2.54,
    height = 3)

# layout(matrix(c(1,2),nrow = 2))
# layout.show(2)

par(mar = c(2,0,0,0))
plot(winners_losers_rast,
     box = FALSE,
     axes = FALSE,
     col = cols,
     legend = FALSE)#,
     # colNA = 'lightgrey')#,
     # mar = c(0,0,0,0))

legend(x = grconvertX(.01, from = "ndc", to = "user"),
       y = grconvertY(0.15, from = "ndc", to = "user"),
       legend = levels(winners_losers_rast)[[1]]$name[1:3], fill = (cols[1:3]), 
       bty = "n", cex = 0.6, 
       # inset = c(0.1, 0.01),
       adj = c(0,.5),
       y.intersp = .8,
       # horiz = 'TRUE',
       xpd = TRUE)

legend(x = grconvertX(.32, from = "ndc", to = "user"),
       y = grconvertY(0.15, from = "ndc", to = "user"),
       legend = levels(winners_losers_rast)[[1]]$name[4:6], fill = (cols[4:6]), 
       bty = "n", cex = 0.6, 
       # inset = c(0.1, 0.01),
       adj = c(0,.5),
       y.intersp = .8,
       # horiz = 'TRUE',
       xpd = TRUE)

legend(x = grconvertX(.65, from = "ndc", to = "user"),
       y = grconvertY(0.15, from = "ndc", to = "user"),
       legend = levels(winners_losers_rast)[[1]]$name[7:9], fill = (cols[7:9]), 
       bty = "n", 
       cex = 0.6, 
       # inset = c(0.1, 0.01),
       adj = c(0,.5),
       y.intersp = .8,
       # horiz = 'TRUE',
       xpd = TRUE)


text(x = grconvertX(0.00001, from = "ndc", to = "user"),
       y = grconvertY(0.99999, from = "ndc", to = "user"),
       label = '(a)',
       # bty = "n", 
       cex = 0.6, 
       # inset = c(0.1, 0.01),
       adj = c(0,1),
       # horiz = 'TRUE',
       xpd = TRUE)

plot(country_borders,
     add = TRUE,
     colNA = 'lightgrey',
     lwd = .5)


dev.off()

###
# Getting spearmans correlations for heatmaps
###
# Getting habitat loss
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
# Filtering to species in analysis
# Rearranging to wide format
out_df_heatmap <-
  out_df %>%
  filter(year %in% 2050) %>%
  filter(paste0(taxon,species) %in%
           paste0(species_list$taxon,species_list$species))

###
# Begging plot
# Setting plotting parameters



###
# Now getting plotting for ecoregions by taxon
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
# Filtering to species in analysis
# Rearranging to wide format
out_df_heatmap <-
  out_df %>%
  filter(year %in% 2050) %>%
  filter(paste0(taxon,species) %in%
           paste0(species_list$taxon,species_list$species))

###
# Adding ecoregions
out_df_ecoregion <- 
  left_join(out_df_heatmap,
            read_csv('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/species_by_ecoregion_esh_maps.csv') %>%
              mutate(species = gsub('_[A-Z].*|_[0-9].*','',species)) %>%
              mutate(binomial = paste0(taxa,'/',species)) %>%
              distinct()) %>%
  distinct()

###
# Now looping to get outcomes by stressor
outcome_list <- unique(out_df_heatmap$scenario)
out_df_eco_taxon <- data.frame()


for(tt in unique(out_df_heatmap$taxon)) {
  for(ii in outcome_list) {
    # Getting outcome
    tmp_outcome_ii <-
      out_df_ecoregion %>%
      filter(scenario %in% ii,
             taxon %in% tt) %>%
      dplyr::group_by(taxon, ecoregion_id) %>%
      dplyr::summarise(tot_pop_ii = median(tot_pop, na.rm = TRUE)) %>%
      mutate(tot_pop_ii = rank(tot_pop_ii))
    
    # Looping through other outcome
    for(yy in outcome_list[!(outcome_list %in% ii)]) {
      # Getting other outcome
      tmp_outcome_yy <-
        out_df_ecoregion %>%
        filter(scenario %in% yy,
               taxon %in% tt) %>%
        dplyr::group_by(taxon, ecoregion_id) %>%
        dplyr::summarise(tot_pop_yy = median(tot_pop, na.rm = TRUE)) %>%
        mutate(tot_pop_yy = rank(tot_pop_yy))
      
      # Merging data frames
      tmp_outcome_out <-
        left_join(tmp_outcome_ii %>% dplyr::select(taxon, ecoregion_id, tot_pop_ii),
                  tmp_outcome_yy %>% dplyr::select(taxon, ecoregion_id, tot_pop_yy))
      
      # Getting spearmans correlation
      cor_out <-
        cor.test(x = tmp_outcome_out$tot_pop_ii, 
                 y = tmp_outcome_out$tot_pop_yy,
                 method = 'spearman')
      
      out_df_eco_taxon <-
        rbind(out_df_eco_taxon,
              data.frame(taxon = tt,
                         outcome_ii = ii,
                         outcome_yy = yy,
                         rho = cor_out$estimate,
                         p_value = cor_out$p.value))
      
      out_lm <- lm(tmp_outcome_out$tot_pop_ii ~ tmp_outcome_out$tot_pop_yy)
    }
  }
}


# Changing names
out_df_eco_taxon_plot <-
  out_df_eco_taxon %>%
  mutate(outcome_ii = 
           case_when(grepl('.int.urb',outcome_ii) ~ 'All Stressors',
                     grepl('urb',outcome_ii) ~ 'Urban Expansion',
                     grepl('.int',outcome_ii) ~ 'Agricultural Intensification',
                     grepl('.exp',outcome_ii) ~ 'Agricultural Expansion',
                     grepl('current',outcome_ii) ~ 'Climate Change')) %>%
  mutate(outcome_yy = 
           case_when(grepl('.int.urb',outcome_yy) ~ 'All Stressors',
                     grepl('urb',outcome_yy) ~ 'Urban Expansion',
                     grepl('.int',outcome_yy) ~ 'Agricultural Intensification',
                     grepl('.exp',outcome_yy) ~ 'Agricultural Expansion',
                     grepl('current',outcome_yy) ~ 'Climate Change')) %>%
  transform(outcome_ii = factor(outcome_ii, levels = c('All Stressors','Climate Change','Agricultural Expansion','Agricultural Intensification', 'Urban Expansion'))) %>%
  transform(outcome_yy = factor(outcome_yy, levels = rev(c('All Stressors','Climate Change','Agricultural Expansion','Agricultural Intensification', 'Urban Expansion')))) %>%
  mutate(outcome_ii_numeric = as.numeric(outcome_ii),
         outcome_yy_numeric = 
           case_when(grepl('Urban',outcome_yy) ~ 1,
                     grepl('Ag.*Int',outcome_yy) ~ 2,
                     grepl('Ag.*Exp',outcome_yy) ~ 3,
                     grepl('Climate',outcome_yy) ~ 4)) %>%
  mutate(outcome_ii_numeric =
           case_when(taxon %in% 'Birds' ~ outcome_ii_numeric + 5,
                     taxon %in% 'Mammals' ~ outcome_ii_numeric + 5,
                     taxon %in% 'Reptiles' ~ outcome_ii_numeric + 0,
                     .default = outcome_ii_numeric)) %>%
  mutate(outcome_yy_numeric =
           case_when(taxon %in% 'Amphibians' ~ outcome_yy_numeric + 5,
                     taxon %in% 'Birds' ~ outcome_yy_numeric + 5,
                     taxon %in% 'Reptiles' ~ outcome_yy_numeric + 0,
                     .default = outcome_yy_numeric)) %>%
  mutate(x_start = outcome_ii_numeric - 1,
         x_end = outcome_ii_numeric,
         y_start = outcome_yy_numeric - 1,
         y_end = outcome_yy_numeric) %>%
  # filter(outcome_ii != outcome_yy) %>%
  # filter(!grepl('All',outcome_yy)) %>%
  filter(!(grepl('Climate',outcome_yy) & grepl('Agricul',outcome_ii))) %>%
  filter(!(grepl('Agr.*Exp',outcome_yy) & grepl('Agri.*Inte',outcome_ii))) %>%
  filter(!(grepl('Urban',outcome_ii))) %>%
  mutate(rho_merge = round(rho,digits = 2)) %>%
  mutate(text_label = 
           case_when(p_value < 0.001 ~ paste0(round(rho,digits = 2),'***'),
                     p_value < 0.01 ~ paste0(round(rho,digits = 2),'**'),
                     p_value <= 0.05 ~ paste0(round(rho, digits = 2),'*'),
                     p_value > 0.05 ~ as.character(round(rho,digits = 2))))

col_df <-
  data.frame(col = 
              c(colorRampPalette(c('#3581b9','#99e9fc','#cbf1f7','#dddddd'))(100),
                '#dddddd',
                colorRampPalette(c('#dddddd','#f2bf91','#f28379','#d95252'))(100)),
  # data.frame(col = rev(colorRampPalette(brewer.pal('RdBu',n=11))(201)),
             rho_merge = seq(from = -1, to = 1, length.out = 201)) %>%
  mutate(rho_merge = round(rho_merge, digits = 2))
  

out_df_eco_taxon_plot <-
  left_join(out_df_eco_taxon_plot,
            col_df)

###
# Plotting heatmap
pdf(paste0(getwd(),'/Figures and Tables/Figures/Fig3_b_heatmaps_',Sys.Date(),'.pdf'),
    width = 12 / 2.54,
    height = 3)

par(mar = c(4.5,5.5,0.5,.5))
plot(NA, 
     xlim = c(0,max(out_df_eco_taxon_plot$x_end)),
     ylim = c(0,9), 
     xlab = "", 
     ylab = "", 
     axes = FALSE,
     cex.lab = .6)

axis(side = 2, 
     at = c(.5:3.5), 
     labels = rep(levels(out_df_eco_taxon_plot$outcome_yy) %>% .[!grepl('All',.)],1) %>%
       gsub('Agricultural','Ag.',.), 
     las = 1,
     cex.axis = .6)

axis(side = 2, 
     at = c(5.5:8.5), 
     labels = rep(levels(out_df_eco_taxon_plot$outcome_yy) %>% .[!grepl('All',.)],1) %>%
       gsub('Agricultural','Ag.',.), 
     las = 1,
     cex.axis = .6)

###
# x axis labels
x_axis_labels <-
  out_df_eco_taxon_plot %>%
  dplyr::select(x_start, 
                outcome_ii) %>%
  distinct() %>%
  arrange(x_start)

axis(side = 1, 
     at = x_axis_labels$x_start[x_axis_labels$x_start<5] + .5, 
     labels = FALSE,
     las = 1,
     cex.axis = .6,
     xpd = TRUE)

axis(side = 1, 
     at = x_axis_labels$x_start[x_axis_labels$x_start>=5] + .5, 
     labels = FALSE,
     las = 1,
     cex.axis = .6,
     xpd = TRUE)

# x axis label text
text(x = x_axis_labels$x_start + .5,
     y = -1,
     label = x_axis_labels$outcome_ii %>% gsub('Agricultural','Ag.',.),
     xpd = TRUE,
     cex = .6,
     srt = 45,
     adj = c(1,.5))

for(ii in 1:nrow(out_df_eco_taxon_plot)) {
  rect(xleft = out_df_eco_taxon_plot$x_start[ii],
       xright = out_df_eco_taxon_plot$x_end[ii],
       ybottom = out_df_eco_taxon_plot$y_start[ii],
       ytop = out_df_eco_taxon_plot$y_end[ii],
       col = out_df_eco_taxon_plot$col[ii])
}

text(x = out_df_eco_taxon_plot$x_start + .5,
     y = out_df_eco_taxon_plot$y_start + .5,
     label = out_df_eco_taxon_plot$text_label,
     adj = c(.5,.5),
     cex = .6)

text(x = grconvertX(0.00001, from = "ndc", to = "user"),
     y = grconvertY(0.99999, from = "ndc", to = "user"),
     label = '(b)',
     adj = c(0,1),
     cex = .6,
     xpd = TRUE)

###
# Adding labels for taxa
taxa_labels <-
  out_df_eco_taxon_plot %>%
  dplyr::select(taxon,x_start,y_end) %>%
  distinct() %>%
  dplyr::group_by(taxon) %>%
  dplyr::summarise(x_mean = mean(x_start, na.rm = TRUE),
                   y_mean = max(y_end, na.rm = TRUE))

text(x = taxa_labels$x_mean + .75,
     y = taxa_labels$y_mean,
     label = taxa_labels$taxon,
     adj = c(.5,.5),
     cex = .6,
     xpd = TRUE)


# box()

dev.off()
# 
# 
# ###
# # Plotting heatmaps
# ###
# # Plotting correlation between stressors across all species
# heatmap_plot <- 
#   ggplot(out_df_eco_taxon %>%
#          filter(outcome_ii != 'Urban Expansion',
#                 outcome_yy != 'All Stressors') %>%
#          filter(!(grepl('Climate',outcome_yy) & grepl('Agricul',outcome_ii))) %>%
#          filter(!(grepl('Agr.*Exp',outcome_yy) & grepl('Agri.*Inte',outcome_ii))), 
#        aes(x = outcome_ii, y = outcome_yy, fill = rho)) +
#   geom_tile() +
#   geom_text(aes(label = round(rho,digits = 3)), size = 6 / .pt) + 
#   scale_fill_distiller(palette = 'RdBu') +
#   # scale_fill_continuous(palette = brewer.pal('RdBu',n=11)) +
#   labs(x = '', y = '', fill = "Spearman's Rho") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   facet_wrap(.~taxon, nrow = 1) +
#   theme(text = element_text(size = 6)) +
#   theme(legend.position = 'none')
# 
# 
# 
# ###
# # Saving outputs
# out_plot_1 <- 
#   plot_grid(raster_plot,legend_plot,
#             nrow = 1,
#             rel_widths = c(3,1))
# out_plot <-
#   plot_grid(out_plot_1,
#             heatmap_plot,
#             ncol = 1,
#             labels = c('(a)','(b)'),
#             label_size = 6)# ,
#             # rel_heights = c(1,1),
#             # rel_widths = c(3,2),
#             # align = 'hv')
# 
# out_plot
# 
# ggsave('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Figures and Tables/Figures/Fig3_2025.05.07.pdf',
#        width = 180,
#        height = 120,
#        units = 'mm')
