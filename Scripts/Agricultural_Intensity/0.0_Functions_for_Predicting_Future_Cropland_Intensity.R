# Functions for predicting future cropland intensity


# Creating function to weight models based on probability of each outcome (e.g. minimal, moderate, and high cropland intensity)
weight.fun <-
  function(mod.list, # list of models
           weight.list, # list of model weights
           df, # data frame used to predict model
           mod.type) { # type of model
    for(i in 1:length(mod.list)) { # Looping through model
      probs.tmp = predict(mod.list[[i]], df, 'prob') %>% # predicted probabilities
        as.data.frame(.) # Converting to data frame
      names(probs.tmp) = c('Minimal','Moderate','Intense') # Adding names
      
      probs.keep <- # Data frame with weighted probability for model i
        probs.tmp %>%
        mutate(Minimal = Minimal * weight.list[[i]], # Weighting outcome by the model weight
               Moderate = Moderate * weight.list[[i]],
               Intense = Intense * weight.list[[i]])
      
      if(i == 1) {
        probs.out = probs.keep # And data frame with weighted probability across all models
      } else {
        probs.out = probs.out + probs.keep
      }
    }
    
    probs.out <-
      probs.out %>%
      mutate(Predicted_Classification = # Predicted classification is an ordered factor
               ifelse(Minimal > Moderate & Minimal > Intense, 'Minimal',
                      ifelse(Moderate > Intense, 'Moderate', 'Intense')))
    
    # Changing names to be reflective of the model type
    names(probs.out) <- 
      paste(names(probs.out),
            "_",
            mod.type,
            sep = "")
    
    # Adding columns to data frame
    df <- cbind(df, probs.out) 
    
    return(df)
  }








# Now making this look good
# Making confusion matrix in a way that meshes with ggplot
conf.plot.function <-
  function(df, model_type) {
    df[,'tmp_outcome'] <- df[,grep(paste0('predicted.*', model_type), names(df), ignore.case = TRUE)]
    
    tmp.sum <- 
      data.frame(land_cover_classification = c(rep('Minimal',3),rep('Moderate',3),rep('Intense',3)),
                 tmp_outcome = rep(c('Minimal','Moderate','Intense'),3)) %>%
      left_join(., df %>% group_by(land_cover_classification, tmp_outcome) %>% dplyr::summarise(Count = n())) %>%
      left_join(., df %>% group_by(land_cover_classification) %>% dplyr::summarise(Count_Actual = n())) %>%
      mutate(Count = ifelse(is.na(Count), 0 , Count)) %>%
      mutate(plot.x = ifelse(land_cover_classification %in% 'Minimal',1,
                             ifelse(land_cover_classification %in% 'Moderate',2,3))) %>%
      mutate(plot.y = ifelse(tmp_outcome %in% 'Minimal',1,
                             ifelse(tmp_outcome %in% 'Moderate',2,3))) %>%
      transform(tmp_outcome = factor(tmp_outcome, levels = c('Minimal','Moderate','Intense'), ordered = TRUE)) %>%
      mutate(percent_predicted = Count/Count_Actual) %>%
      as.data.frame(.) %>%
      mutate(label = paste('Predicted: ',Count,'\n',
                           'Actual: ', Count_Actual,'\n',
                           'Percent Predicted: ',round(percent_predicted, digits = 2)*100,"%",
                           sep = ""))
    
    # Setting colour palette
    # Color scale
    my.palette = c('#E6E1B9','#5BC2AE','#0A6164') %>%
      colorRampPalette(.) 
    palette.table = 
      data.frame(value = 0:100,
                 color = my.palette(101))
    
    # Merging in colour palette
    tmp.sum1 <-
      tmp.sum %>%
      mutate(value = round(percent_predicted, digits = 2) * 100) %>%
      left_join(., palette.table)
    
    # Manging data to make the plot
    tmp.plot <- 
      tmp.sum1 %>%
      dplyr::select(land_cover_classification, tmp_outcome,plot.x,plot.y,percent_predicted,color, label) %>%
      mutate(plot.y = ifelse(land_cover_classification %in% 'Minimal',1,
                             ifelse(land_cover_classification %in% 'Moderate',2,3))) %>%
      mutate(plot.x = ifelse(tmp_outcome %in% 'Minimal',1,
                             ifelse(tmp_outcome %in% 'Moderate',2,3))) %>%
      mutate(color = ifelse(land_cover_classification == tmp_outcome, 'white','black'))
    
    # And plotting
    out.plot = ggplot(tmp.plot) +
      geom_tile(aes(x = plot.x, y = plot.y, fill = percent_predicted * 100)) +
      geom_text(aes(x = plot.x, y = plot.y, label = label, colour = color), size = 2.5, show.legend = FALSE) +
      scale_fill_gradientn(colours=my.palette(100), limits = c(0,100), expand = c(0,0)) +
      scale_colour_manual(values = c('black','white')) +
      scale_x_continuous(breaks = c(1,2,3),labels = (c('Minimal','Moderate','Intense'))) +
      scale_y_reverse(breaks = c(1,2,3),labels = rev(c('Intense','Moderate','Minimal'))) +
      labs(x = 'Predicted Crop Intensity', y = 'Observed Crop Intensity', fill = 'Relative\nPercent') +
      theme(axis.text = element_text(size = 7.5)) +
      theme(axis.title = element_text(size = 10)) +
      ggtitle(paste0("Weighted ", model_type, " Model"))
    
    return(out.plot)
    
  }


# This projects future cropland intensity
# And is designed to be used with mclapply
# y is a list of years over which to loop
project.intensity.function <-
  function(y) {
    # Getting raster for target year
    crop.forecasts.years <- 
      crop.forecasts[grepl(y, crop.forecasts)]
    
    # Importing raster
    tmp.raster.keep <- raster(crop.forecasts.years)
    
    # Cropping and reprojecting
    tmp.raster <- projectRaster(tmp.raster.keep, intens, method = 'bilinear')
    
    # Getting adjacent cropland extent
    # Converting to a matrix for convolution
    input_array <- matrix(getValues(tmp.raster), nrow = tmp.raster@nrows, byrow = TRUE)
    # Crop convolved
    py_run_string("input_array = copy.copy(r.input_array)")
    py_run_string("input_array[np.isnan(input_array)] = 0")
    # py_run_string("input_array[input_array > 0] = 1")
    py_run_string("crop_convolved = convolve(input_array, circle, mode='constant')")
    # Recreating raster
    crop_convolved <- raster(matrix(py$crop_convolved,byrow = FALSE, nrow = tmp.raster@nrows),
                             crs = crs(tmp.raster))
    extent(crop_convolved) <- extent(tmp.raster)
    
    # Creating data frame and predicting
    cover.df <-
      data.frame(intens = raster::getValues(intens), 
                 crop = raster::getValues(tmp.raster),
                 adj_crop = raster::getValues(crop_convolved),
                 iso = as.numeric(getValues(country.iso))) %>% 
      mutate(adj_crop_by_intens = adj_crop * intens) %>%
      left_join(., yields %>% mutate(iso = countrycode(ISO3, origin = 'iso3c', destination = 'iso3n')) %>% mutate(iso = ifelse(ISO3 %in% 'ZZZ',728,iso)) %>% unique(.)) %>%
      mutate(yields_increase = .[,which(names(.) %in% paste0('X',gsub(".*_","",y)))]) %>%
      mutate(yields_increase = ifelse(iso %in% 729, unique(.$yields_increase[.$iso %in% 728]),yields_increase)) %>% # updating yields increase for sudan and south sudan
      mutate(yields_increase = ifelse(is.na(yields_increase),1,yields_increase)) %>% # Converting countires with no projections to a yields increase of 0%. This is conservative, but not a terrible issue beause almost all of these countries are island nations
      mutate(intens = intens * yields_increase)
    
    # Predicting intensity Random forest
    cover.df <-
      weight.fun(mod.list = mod.list.rf,
                 weight.list = model.acc.rf1,
                 df = cover.df,
                 mod.type = 'RandomForest')
    
    # Classifying into low medium and high
    # ANd making sure intensity predictions only exist in cells where there is cropland
    cover.df <-
      cover.df %>%
      mutate(predicted_intens = ifelse(Predicted_Classification_RandomForest %in% 'Minimal',1,
                                       ifelse(Predicted_Classification_RandomForest %in% 'Moderate',2,
                                              ifelse(Predicted_Classification_RandomForest %in% 'Intense',3,0)))) %>%
      mutate(predicted_intens = ifelse(crop %in% 0 | is.na(crop), 0, predicted_intens))
    
    
    # Making intensity raster
    intens.cov.raster <-
      raster(matrix(cover.df$predicted_intens, byrow = TRUE, nrow = intens@nrows), crs = crs(intens))
    extent(intens.cov.raster) <- extent(intens)
    
    # and reprojecting
    intens.cov.raster <- projectRaster(intens.cov.raster, tmp.raster.keep, method = 'ngb')
    
    return(intens.cov.raster)
  }

