###
# Plotting change in urbanisation from 2020 to 2050

# libraries
library(raster)
library(RColorBrewer)
library(viridis)
library(terra)

# loading rasters
int_2020 <- rast(paste0(getwd(),'/Ag Intensity Outputs/Cropland Intensity Forecasts/BAU/Cropland_Intensity_2015_2020.tif'))
int_2050 <- rast(paste0(getwd(),'/Ag Intensity Outputs/Cropland Intensity Forecasts/BAU/Cropland_Intensity_2045_2050.tif'))

# Getting other land areas
crop_2020 <- rast(paste0(getwd(),'/Land Forecast Outputs/BAU/Global tifs/crop_mean2015_2020.tif'))

int_2020[!is.na(crop_2020) & is.na(int_2020)] <- 0
int_2050[!is.na(crop_2020) & is.na(int_2050)] <- 0

delta_int <- int_2050 - int_2020

# 0 corresponds to not suitable in 2020 or 2050 (by growing length)
# 1 corresponds to suitable under both 5 and 10c
# 2 corresponds to suitable under 5 not 10
# 3 corresponds to suitable under 10 not 5

# Colour palette
red_palette <- colorRampPalette(c('#ffffff','#f2bf91','#f28379','#d95252'))(4)
blue_palette <- colorRampPalette(rev(c('#3581b9','#99e9fc','#cbf1f7','#ffffff')))(4)

my_palette <- 
  c(rev(blue_palette),
    red_palette[2:4])


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
         legend = 'bottomleft',
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

# Changing levels on the raster
lvl <- data.frame(
  value = 0:3,  # or 1:4, depending on your data
  label = c(
    "No Agriculture",
    "Low Intensity Agriculture",
    "Moderate Intensity Agriculture",
    "High Intensity Agriculture"
  )
)

levels(int_2020) <- lvl
levels(int_2050) <- lvl


###
# plotting
dev.off()
pdf(paste0(getwd(),'/Figures and Tables/Figures/FigS6_AgIntensity_',Sys.Date(),'.pdf'),
    width = 12, height = 4)
par(mfrow = c(1,3))

# cropland
plot_function(raster_input = int_2020,
              palette = viridis(4),
              limits = c(0,3),
              title = 'Agricultural Intensity in 2020')

plot_function(raster_input = int_2050,
              palette = viridis(4),
              # limits = c(-1,1),
              title = 'Agricultural Intensity in 2050')

plot_function(raster_input = delta_int,
              palette = my_palette[2:7],
              limits = c(-3,3),
              title = 'Change in Agricultural Intensity 2020 - 2050')


dev.off()



###
# Getting values for table S3 - which contains number of cells in no ag / low / mod / high intensity ag in 2020 and 2050
tmp_df <-
  data.frame(int_2020 = values(int_2020),
             int_2050 = values(int_2050)) 

tmp_df_2 <-
  tmp_df %>%
  mutate(tot_int = label + label.1) %>%
  dplyr::rename(int_2020 = label,
                int_2050 = label.1)

tmp_df_3 <-
  tmp_df_2 %>%
  filter(tot_int > 0) %>%
  # filter(!is.na(int_2020) & !is.na(int_2050)) %>%
  mutate(int_2020_class = 
           case_when(int_2020 %in% 0 | is.na(int_2020) ~ 'No Cropland in 2020',
                     int_2020 %in% 1 ~ 'Low Intensity in 2020',
                     int_2020 %in% 2 ~ 'Moderate Intensity in 2020',
                     int_2020 %in% 3 ~ 'High Intensity in 2020')) %>%
  mutate(int_2050_class = 
           case_when(int_2050 %in% 0 | is.na(int_2050) ~ 'No Cropland in 2050',
                     int_2050 %in% 1 ~ 'Low Intensity in 2050',
                     int_2050 %in% 2 ~ 'Moderate Intensity in 2050',
                     int_2050 %in% 3 ~ 'High Intensity in 2050')) %>%
  transform(int_2020_class = 
              factor(int_2020_class,
                     levels = c('No Cropland in 2020',
                                'Low Intensity in 2020',
                                'Moderate Intensity in 2020',
                                'High Intensity in 2020'))) %>%
  transform(int_2050_class = 
              factor(int_2050_class,
                     levels = c('No Cropland in 2050',
                                'Low Intensity in 2050',
                                'Moderate Intensity in 2050',
                                'High Intensity in 2050')))
out_df <- 
  tmp_df_3 %>%
  dplyr::group_by(int_2020_class,
                  int_2050_class) %>%
  dplyr::summarise(count = n()) %>%
  pivot_wider(names_from = int_2050_class,
              values_from = count)

write.csv(out_df,
          paste0(getwd(),'/Figures and Tables/Tables/TableS3_CroplandIntensity_',Sys.Date(),'.csv'),
          row.names = FALSE)
