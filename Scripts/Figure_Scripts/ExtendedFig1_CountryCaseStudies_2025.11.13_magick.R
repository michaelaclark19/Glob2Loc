###
# libraries
library(plyr)
library(dplyr)
library(parallel)
library(terra)
library(countrycode)

###
# files
file_list <-
  list.files('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Country_Analyses_Transfer/',
             full.names = TRUE)

###
# countries
countries <- c('AUS','UZB','TZA',
               'ECU',
               #'EGY',
               'SWE','JPN','GHA')

###
# stressors
stressor_list <-
  data.frame(title = 
               c('Agricultural Expansion',
                 'Agricultural Intensification',
                 'Climate Change',
                 'Urban Expansion'),
             stressor = c('esh.exp.tif',
                          'esh.int.tif',
                          'current',
                          'esh.urb.tif')) 

###
# taxa
taxa_list <-
  c('Amphibians','Birds','Mammals','Reptiles')

###
# Palette
# red_palette <- c('#ffffff','#f2bf91','#f28379','#d95252')
# blue_palette <- c('#3581b9','#99e9fc','#cbf1f7','#ffffff')
# 
# my_palette <- 
#   c(blue_palette,
#     red_palette) %>%
#   unique()

red_palette <- colorRampPalette(c('#ffffff','#f2bf91','#f28379','#d95252'))(100)
blue_palette <- colorRampPalette(rev(c('#3581b9','#99e9fc','#cbf1f7','#ffffff')))(100)

my_palette <- 
  rev(c(rev(blue_palette),
        '#ffffff',
        red_palette))

###
# Making plot function
# plot_function <-
#   function(raster_input,
#            # palette,
#            #limits,
#            legend_title,
#            title) {
#     
#     plot_vals <- values(raster_input) 
#     plot_vals <- plot_vals[plot_vals != 0]
#     plot_limits <-
#       c(quantile(plot_vals,.05, na.rm=TRUE),
#         max(plot_vals, na.rm = TRUE))
#     
#     plot_limits <-
#       as.vector(c(floor(plot_limits[1]/10) * 10,
#                 ceiling(plot_limits[2]/10) * 10))
#     
#     plot_limits <- c(quantile(plot_vals,c(.05,.95),na.rm=TRUE))
#     plot_limits <- c(-max(abs(plot_limits)),
#                      max(abs(plot_limits)))
#     
#     if(max(plot_limits) < 0) {
#       my_palette <- colorRampPalette(my_palette[7:4])(100)
#       plot_limits <- c(plot_limits[1],0)
#       raster_input[raster_input < plot_limits[1]] <- plot_limits[1]
#     } else {
#       palette_neg <- colorRampPalette(my_palette[7:4])(50)
#       palette_pos <- colorRampPalette(my_palette[1:4])(50)
#       my_palette_plot <- 
#         c(palette_neg,
#           '#ffffff',
#           rev(palette_pos))
#       
#       raster_input[raster_input < plot_limits[1]] <- plot_limits[1]
#       raster_input[raster_input > plot_limits[2]] <- plot_limits[2]
#     }
# 
#     plot(raster_input,
#          plg = list(loc = 'left', title = legend_title, cex = 1, title.cex = 1, legend.size = .1),
#          # colNA = brewer.pal(n = 9, name = 'Blues')[2],
#          col = my_palette_plot,
#          # col = my_palette,
#          main = title,
#          cex.main = 1,
#          axes = FALSE,
#          # mar = c(0,0,0,0),
#          # add = TRUE,
#          range = plot_limits)
#     
#     plot(country_shapes,
#          add = TRUE,
#          border = 'black',
#          lwd = .5)
#   }


plot_function <-
  function(raster_input,
           palette,
           # limits,
           legend.title,
           title,
           multiplier,
           file_name) {
    # cropland
    # plot(r_oval, 
    #      # col = 'light blue',
    #      col = brewer.pal(n = 9, name = 'Blues')[2],
    #      # col = '#eeeeee',
    #      box = FALSE,
    #      axes = FALSE,
    #      legend = FALSE,
    #      mar = c(0,0,0,0))
    
    pdf(paste0(folder_save_tmp_images,file_name),
        height = tot_height * multiplier,
        width = tot_width * multiplier)
    
    # par(mar = c(3,1,1,1))
    
    plot_vals <- values(raster_input) 
    plot_vals <- plot_vals[plot_vals != 0]
    plot_limits <-
      c(quantile(plot_vals,.05, na.rm=TRUE),
        max(plot_vals, na.rm = TRUE))
    
    plot_limits <-
      as.vector(c(floor(plot_limits[1]/10) * 10,
                  ceiling(plot_limits[2]/10) * 10))
    
    plot_limits <- c(quantile(plot_vals,c(.05,.95),na.rm=TRUE))
    plot_limits <- c(-max(abs(plot_limits)),
                     max(abs(plot_limits)))
    
    raster_input[raster_input < plot_limits[1]] <- plot_limits[1]
    raster_input[raster_input > plot_limits[2]] <- plot_limits[2]
    
    plot(raster_input,
         plg = list(loc = 'left', title = 'Proportion\nof Cell', cex = 1, title.cex = 1, legend.size = .1),
         # colNA = brewer.pal(n = 9, name = 'Blues')[2],
         col = palette,
         # col = my_palette,
         main = title,
         cex.main = 1,
         axes = FALSE,
         legend = FALSE,
         mar = c(2,1,2,1),
         range = plot_limits)
    
    plot(country_shapes,
         add = TRUE,
         border = 'black',
         lwd = .5)
    
    ### 
    # Adding legend 
    
    # Define rectangle bounds for the legend colour bar
    xleft <- .3
    xright <- .7
    ybottom <- .05
    ytop <- .1
    
    if(multiplier == .5) {
      ybottom <- 0.035
      ytop <- 0.06
    }
    
    # Convert to user coordinates
    xleft <- grconvertX(xleft, from = "ndc", to = "user")
    xrightlabel <- grconvertX(xright+.025, from = "ndc", to = "user")
    ybottomlabel <- grconvertY(ybottom-.0125, from = "ndc", to = "user")
    xright <- grconvertX(xright, from = "ndc", to = "user")
    ybottom <- grconvertY(ybottom, from = "ndc", to = "user")
    ytitlelabel <- grconvertY(ytop+.02, from = "ndc", to = "user")
    ytop_pos <- grconvertY(ytop + .05, from = "ndc", to = "user")
    ybottom_pos <- grconvertY(ytop + .00, from = "ndc", to = "user")
    ytop <- grconvertY(ytop, from = "ndc", to = "user")
    
    
    
    # Number of colours
    n_colors <- length(palette)
    
    # Compute width of each sub-rectangle
    x_steps <- seq(xleft, xright, length.out = n_colors + 1)
    
    # Adding colour bar
    # Do this in a loop
    for (i in 1:n_colors) {
      rect(xleft = x_steps[i], ybottom = ybottom,
           xright = x_steps[i+1], ytop = ytop,
           col = palette[i], border = NA,
           xpd = TRUE)
    }
    
    # Adding black rectngle around colour bar
    rect(xleft, ybottom, xright, ytop, border = "black", xpd = TRUE)
    
    # Adding text for the colour bar
    text(x = xleft, 
         y = ybottomlabel, 
         label = round(plot_limits[2],digits = 2),
         xpd = TRUE,
         adj = c(.5,1),
         cex = .6)
    
    
    text(x = mean(c(xleft,xright)), 
         y = ybottomlabel, 
         label = 0,
         xpd = TRUE,
         adj = c(.5,1),
         cex = .6)
    
    text(x = xrightlabel, 
         y = ybottomlabel, 
         label = round(plot_limits[1], digits = 2),
         xpd = TRUE,
         adj = c(.5,1),
         cex = .6)
    
    # Adding legend title
    text(x = mean(c(xleft,xright)),
         y = ytitlelabel,
         label = legend.title,
         adj = c(.5,.5),
         cex = .6,
         xpd = TRUE)
    
    
    if(multiplier == .5) {
      
      text(x = xrightlabel - .25 * (xright - xleft), 
           y = ybottomlabel, 
           label = round(plot_limits[1]/2, digits = 2),
           xpd = TRUE,
           adj = c(.5,1),
           cex = .6)
      
      
      text(x = xleft + .25 * (xright - xleft), 
           y = ybottomlabel, 
           label = round(plot_limits[2]/2, digits = 2),
           xpd = TRUE,
           adj = c(.5,1),
           cex = .6)
    }
    
    
    dev.off()
    
  }


###
# Figure layout:
# Top row - all species all stressors: (left) = abs change in hab area; (right) = % change in hab area
# Second row - by taxa, all stressors: amphibians, birds, mammals, reptiles
# Third row - by stressor, all taxa: ag exp, ag int, climate change, urbanisation


###
# Folder for tmp pdf files
folder_save_tmp_images <- '/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Figures and Tables/Figures/Supplemental_Country_Figures/Temp_Country_Panels/'

###
# Looping through countries
for(cc in countries) {
  
  ###
  # getting rast files for the country
  file_countries <-
    file_list %>%
    .[grepl(cc,.)]
  
  tmp_rast <-
    rast(file_countries[1])
  
  # fig size ratio...
  ratio <- nrow(tmp_rast) / ncol(tmp_rast)
  
  tot_width = 10
  tot_height = tot_width * ratio
  
  if(tot_height < 10) {
    tot_height = 10
    tot_width = tot_height / ratio
  }
  
  
  # par(mar = c(1,1,1,1))
  
  ###
  # getting shapes 
  country_shapes <-
    ne_states(iso_a2 = 
                countrycode(cc, origin = 'iso3c', destination = 'iso2c'), 
              returnclass = 'sv')
  
  

  ### 
  # Change for all taxa all stressors
  
  # 2020 files
  file_list_2020 <-
    file_countries %>%
    .[grepl('2020',.)] %>%
    .[grepl('hab_availability',.)] %>%
    .[!grepl('richness',.)]
  
  # baseline 2020
  rast_2020 <-
    app(rast(file_list_2020), fun = sum)
  
  # projecting countryshapes
  country_shapes <-
    project(country_shapes, rast_2020)
  
  # all stressors 2050
  file_list_2050 <- 
    file_countries %>%
    .[grepl('2050',.)] %>%
    .[grepl('hab_availability',.)] %>%
    .[!grepl('richness',.)] %>%
    .[grepl('.int.urb',.)]
  
  rast_2050 <-
    app(rast(file_list_2050),fun = sum)
  
  ### getting rasters
  # abs change
  abs_change_2050 <- rast_2050 - rast_2020
  prop_change_2050 <- log2(rast_2050 / rast_2020)
  
  
  ###
  # plotting first two
  plot_function(
    raster_input = abs_change_2050,
    palette = my_palette,
    # limits = c(-200,200),
    legend.title = 'Number of Cells',
    title = 'Absolute Change in Habitat Area',
    multiplier=.5,
    file_name = paste0(cc,'_Abs_Change_Hab_',Sys.Date(),'.pdf'))
  
  plot_function(
    raster_input = prop_change_2050,
    palette = (my_palette),
    # limits = c(-200,200),
    legend.title = 'Proportional Change (log(2))',
    title = 'Proportional Change in Habitat Area',
    multiplier=.5,
    file_name = paste0(cc,'_Prop_Change_Hab_',Sys.Date(),'.pdf'))
  ###
  # Change by taxa all stressors
  for(tt in taxa_list) {
    files_taxa <-
      file_countries %>%
      .[grepl(tt,.)] %>%
      .[grepl('.int.urb',.)] %>%
      .[!grepl('richness',.)] %>%
      .[grepl('hab_avail',.)]
    
    taxa_rast <-
      rast(files_taxa %>% .[grepl('2050',.)]) -
      rast(files_taxa %>% .[grepl('2020',.)])
    
    taxa_rast <- 
      log2((rast(files_taxa %>% .[grepl('2050',.)]) + .0001) /
      (rast(files_taxa %>% .[grepl('2020',.)]) + .0001))
    
    plot_function(raster_input = taxa_rast,
                  palette = my_palette,
                  legend.title = 'Prop. Change (log(2))',
                  title = tt,
                  multiplier = .25,
                  file_name = paste0(cc,'_Taxa_',tt,'_',Sys.Date(),'.pdf'))
  }
  
  ###
  # Change by stressors across all taxa
  for(ss in 1:nrow(stressor_list)) {
    files_stressor <-
      file_countries %>%
      .[grepl(stressor_list$stressor[ss],.)] %>%
      # .[grepl('.int.urb',.)] %>%
      .[!grepl('richness',.)] %>%
      .[grepl('hab_avail',.)] %>%
      .[grepl(2050,.)]
    
    stressor_rast <-
      app(rast(files_stressor), fun = sum) -
      rast_2020
    
    stressor_rast <-
      log2(app(rast(files_stressor), fun = sum) /
      rast_2020)
    
    if(stressor_list$title[ss] %in% 'Agricultural Intensification') {
      files_stressor <-
        file_countries %>%
        .[grepl(stressor_list$stressor[ss],.)] %>%
        # .[grepl('.int.urb',.)] %>%
        .[!grepl('richness',.)] %>%
        .[grepl('pop',.)] %>%
        .[grepl(2050,.)]
      
      rast_2020_pop_abund <-
        file_countries %>%
        # .[grepl(stressor_list$stressor[ss],.)] %>%
        # .[grepl('.int.urb',.)] %>%
        .[!grepl('richness',.)] %>%
        .[grepl('pop',.)] %>%
        .[grepl(2020,.)]
      
      stressor_rast <-
        log2((app(rast(files_stressor), fun = sum) + .0001)/
        (app(rast(rast_2020_pop_abund), fun = sum) + .0001))
      
      plot_function(raster_input = stressor_rast,
                    legend.title = 'Prop. Change (log(2))',
                    title = stressor_list$title[ss] %>% gsub(' ','\n',.),
                    palette = my_palette,
                    multiplier = .25,
                    file_name = paste0(cc,'_Stressor_',stressor_list$title[ss],'_',Sys.Date(),'.pdf'))
      
      
    } else {
      plot_function(raster_input = stressor_rast,
                  legend.title = 'Prop. Change (log(2))',
                  palette = my_palette,
                  title = stressor_list$title[ss],
                  multiplier = .25,
                  file_name = paste0(cc,'_Stressor_',stressor_list$title[ss],'_',Sys.Date(),'.pdf'))
    }
  }
  
  ###
  # And appending together
  
  # List of pdfs for the country
  image_list <-
    list.files(folder_save_tmp_images, 
               full.names = TRUE,
               pattern = cc) %>%
    .[grepl(Sys.Date(),.)]
  
  top_row <-
    image_append(
      c(image_read_pdf(image_list %>% .[grepl('Abs',.)]),
        image_read_pdf(image_list %>% .[grepl('Prop',.)]))
    )
  
  middle_row <-
    image_append(
      c(image_read_pdf(image_list %>% .[grepl('Amp',.)]),
        image_read_pdf(image_list %>% .[grepl('Bird',.)]),
        image_read_pdf(image_list %>% .[grepl('Mamm',.)]),
        image_read_pdf(image_list %>% .[grepl('Rep',.)]))
    )
  
  bottom_row <-
    image_append(
      c(image_read_pdf(image_list %>% .[grepl('Ag.*Exp',.)]),
        image_read_pdf(image_list %>% .[grepl('Ag.*Int',.)]),
        image_read_pdf(image_list %>% .[grepl('Clim',.)]),
        image_read_pdf(image_list %>% .[grepl('Urb',.)]))
    )
  
  image_out <-
    image_append(c(top_row,
                 middle_row,
                 bottom_row),
                 stack = TRUE)
  
  # And saving
  image_write(image_out,
              path = paste0('/Users/michael/Desktop/Research/Multiple_Stresses_of_Biodiversity/Figures and Tables/Figures/Supplemental_Country_Figures/',
                            cc,'_',Sys.Date(),'.pdf'),
              format = 'pdf')
    
  
  
  # dev.off()
}
