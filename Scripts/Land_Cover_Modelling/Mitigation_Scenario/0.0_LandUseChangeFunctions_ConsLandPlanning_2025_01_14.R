# Read me ----
# Script to hold functions
#
# Function for sampling -----
# Much faster than sample()
fastSampleReject <- function(all, n, w){
  out <- numeric(0)
  while(length(out) < n)
    out <- unique(c(out, sample(all, size = n, replace = TRUE, prob = w)))
  out[1:n]
}

# Proportion of surrounding ag or pasture land ----
convolution_fun <- function(m){
  m_rows <- nrow(m)
  m_cols <- ncol(m)
  m_tmp <- rbind(matrix(0, nrow = 1, ncol = m_cols), m)[-(m_rows + 1),]
  m_tmp <- m_tmp + 
    rbind(matrix(0, nrow = 1, ncol = m_cols),
          cbind(m, matrix(0, nrow = m_rows, ncol = 1))[,-1])[-(m_rows + 1),]
  m_tmp <- m_tmp + 
    cbind(m, matrix(0, nrow = m_rows, ncol = 1))[,-1]
  m_tmp <- m_tmp + 
    rbind(cbind(m, matrix(0, nrow = m_rows, ncol = 1))[,-1],
          matrix(0, nrow = 1, ncol = m_cols))[-1,]
  m_tmp <- m_tmp + 
    rbind(m, matrix(0, nrow = 1, ncol = m_cols))[-1,]
  m_tmp <- m_tmp + 
    rbind(cbind(matrix(0, nrow = m_rows, ncol = 1), m)[,-(m_cols + 1)],
          matrix(0, nrow = 1, ncol = m_cols))[-1,]
  m_tmp <- m_tmp + cbind(matrix(0, nrow = m_rows, ncol = 1), m)[,-(m_cols + 1)]
  m_tmp <- m_tmp + 
    rbind(matrix(0, nrow = 1, ncol = m_cols),
          cbind(matrix(0, nrow = m_rows, ncol = 1), m)[,-(m_cols + 1)])[-(m_rows + 1),]
  return(m_tmp)
}

adj.ag.f <- function(m, r, c, crop){
  if(crop == "crop"){
    # Convert the column to a matrix
    tmp_mat <- matrix(m$prop_to_predict_from_crop, 
                      nrow = r, 
                      ncol = c, 
                      byrow = FALSE)
    # Replace NAs with 0s to allow convolution
    tmp_mat[is.na(tmp_mat)] <- 0
    # Run the convolution
    tmp_mat <- convolution_fun(tmp_mat)
    # Put data in and put NAs back where they should be
    m$prop_adj_to_predict_from_crop <- as.vector(tmp_mat)
    m[is.na(m$prop_to_predict_from_crop), "prop_adj_to_predict_from_crop"] <- NA
    return(m)
  } else {
    # Convert the column to a matrix
    tmp_mat <- matrix(m$prop_to_predict_from_pasture, 
                      nrow = r, 
                      ncol = c, 
                      byrow = FALSE)
    # Replace NAs with 0s to allow convolution
    tmp_mat[is.na(tmp_mat)] <- 0
    # Run the convolution
    tmp_mat <- convolution_fun(tmp_mat)
    # Put data in and put NAs back where they should be
    m$prop_adj_to_predict_from_pasture <- as.vector(tmp_mat)
    m[is.na(m$prop_to_predict_from_pasture), "prop_adj_to_predict_from_pasture"] <- NA
    return(m)
  }
}

# Predicting probabilites and amount of conversion ----
# for troubleshooting, p = model_data_current
pred_fun <- function(p, crop){
  if(crop == "crop"){
    
    # And adding in new predictor variables
    # log(dist), log(dist)^2, etc
    p$dist_50k_log <- log(p$dist_50k + .001) # Adding .001 to avoid NAs in calculations
    p$dist_50k_log_sq <- p$dist_50k_log^2
    
    p$prop_crop_sq <- p$prop_to_predict_from_crop^2
    p$prop_pasture_sq <- p$prop_to_predict_from_pasture^2
    
    p$prop_adj_crop_sq <- p$prop_adj_to_predict_from_crop^2
    p$prop_adj_pasture_sq <- p$prop_adj_to_predict_from_pasture^2
    
    p$prop_cell_delta_crop_sq <- p$prop_to_predict_from_cell_delta_crop^2
    p$prop_cell_delta_pasture_sq <- p$prop_to_predict_from_cell_delta_pasture^2
    
    p$ISO_numeric = factor(p$ISO_numeric, ordered = FALSE)
    
    # Cropland -----
    # First: probability of change for crop
    # pred_tmp are temporary predictor values. Gives data frame with probability of either no change, increase, or decrease in crop extent
    # Manual model predictions
    iso3s = unique(p$ISO_numeric) %>% as.character(.) %>% .[!is.na(.)]
    
    if(length(iso3s) > 1) {
      dat.out = data.frame()
      for(i in iso3s) {
        pred_tmp <- multinom.predict.2b.trial(dat = p %>% filter(ISO_numeric %in% i),
                                              coefs = coefs.crop.prob,
                                              vars = predictor_vars_crop[1:(length(predictor_vars_crop) - 1)],
                                              factor_vars = 'ISO_numeric')
        dat.out <- 
          rbind(dat.out,
                p %>% filter(ISO_numeric %in% i) %>% cbind(pred_tmp))
      }
      
      p <- 
        dat.out %>%
        dplyr::rename(a_pred_to_no_change_crop = aNoChange,
                      pred_to_increase_crop = Increase,
                      pred_to_decrease_crop = Decrease)
    } else {
      tmp <-
        multinom.predict.2b.trial(dat = p,
                                  coefs = coefs.crop.prob,
                                  vars = predictor_vars_crop[1:(length(predictor_vars_crop) - 1)],
                                  factor_vars = 'ISO_numeric')
      # And updating names
      p[,'a_pred_to_no_change_crop'] <-
        tmp$aNoChange
      p[,'pred_to_increase_crop'] <-
        tmp$Increase
      p[,'pred_to_decrease_crop'] <-
        tmp$Decrease
      
      # 
      # p <- cbind(p,
      #            multinom.predict.2b.trial(dat = p,
      #                                      coefs = coefs.crop.prob,
      #                                      vars = predictor_vars_crop[1:(length(predictor_vars_crop) - 1)],
      #                                      factor_vars = 'ISO_numeric'))
    }
    
    
    
    # Repeat for decreasing models
    # Proportion to decrease
    # And manually predicting
    # Have to take negative values because the undelrying models can only predict values between 0 and 1
    
    # Also looping through isos
    if(length(iso3s) > 1) {
      dat.out = data.frame()
      for(i in iso3s) {
        dat.out <- rbind(dat.out,
                         p %>% filter(ISO_numeric %in% i) %>%
                           mutate(pred_change_cell_prop_crop_decrease = 
                                    as.vector(-glm.predict.2b(dat = p %>% filter(ISO_numeric %in% i),
                                                              coefs = coefs.crop.decrease,
                                                              vars = predictor_vars_crop[1:(length(predictor_vars_crop) - 1)],
                                                              factor_vars = 'ISO_numeric'))))
        
        
      }

      p <- dat.out
    } else {
      p[,'pred_change_cell_prop_crop_decrease'] =
        as.vector(-glm.predict.2b(dat = p,
                                  coefs = coefs.crop.decrease,
                                  vars = predictor_vars_crop[1:(length(predictor_vars_crop) - 1)],
                                  factor_vars = 'ISO_numeric'))
        
    }
    # Avoiding NAs
    p$pred_change_cell_prop_crop_decrease[is.na(p$pred_change_cell_prop_crop_decrease)] <- 0
    
    # Repeat for increasing models
    # Proportion to increase
    
    # Also looping through isos
    if(length(iso3s) > 1) {
      dat.out = data.frame()
      for(i in iso3s) {
        dat.out <- rbind(dat.out,
                         p %>% filter(ISO_numeric %in% i) %>%
                           mutate(pred_change_cell_prop_crop_increase = 
                                    as.vector(glm.predict.2b(dat = p %>% filter(ISO_numeric %in% i),
                                                             coefs = coefs.crop.increase,
                                                             vars = predictor_vars_crop[1:(length(predictor_vars_crop) - 1)],
                                                             factor_vars = 'ISO_numeric'))))
        
        
      }
      p <- dat.out
    } else {
      p[,'pred_change_cell_prop_crop_increase'] <- 
        as.vector(glm.predict.2b(dat = p,
                                 coefs = coefs.crop.increase,
                                 vars = predictor_vars_crop[1:(length(predictor_vars_crop) - 1)],
                                 factor_vars = 'ISO_numeric'))
    }
   
    # Avoiding NAs
    p$pred_change_cell_prop_crop_increase[is.na(p$pred_change_cell_prop_crop_increase)] <- 0
    
    # Check that the amount of cropland plus change plus urban land plus nonarable land is not greater than the area of the cell.
    # Creates temporary cropland to make sure potential future area of cell after crop expansion does not exceed 1
    # Creating temp data frame and converting NAs to 0s to avoid issues summing across columns
    p.tmp <- dplyr::select(p, 
                           urban_prop, 
                           Prop.NonArable.2010, 
                           prop_to_predict_from_crop, 
                           pred_change_cell_prop_crop_increase, 
                           pred_change_cell_prop_crop_decrease)
    # Converting NAs to 0s
    p.tmp[is.na(p.tmp)] <- 0
    # Creating new column to adjust values of predictions so that they cannot exceed 1
    p.tmp$tmp <- rowSums(p.tmp[,c('urban_prop',
                                  'Prop.NonArable.2010',
                                  'prop_to_predict_from_crop',
                                  'pred_change_cell_prop_crop_increase')],
                         na.rm = TRUE)



    # Adjust predictions if so
    # If area does exceed 1, then adjust so that predicted cropland change does not cause cell value to exceed 1
    p.tmp[p.tmp$tmp > 1, "pred_change_cell_prop_crop_increase"] <-
      1 - (p.tmp[p.tmp$tmp > 1, "urban_prop"] + 
             p.tmp[p.tmp$tmp > 1, "Prop.NonArable.2010"] + 
             p.tmp[p.tmp$tmp > 1, "prop_to_predict_from_crop"])
    
    # Repeat so that the area can't be negative
    p.tmp$tmp <-     p.tmp$tmp <- rowSums(p.tmp[,c('pred_change_cell_prop_crop_increase',
                                                   'prop_to_predict_from_crop')],
                                          na.rm = TRUE)
    
    p.tmp[p.tmp$tmp < 0, "pred_change_cell_prop_crop_increase"] <-
      0 - p.tmp[p.tmp$tmp < 0, "prop_to_predict_from_crop"]
    
    # Repeat for decreasing models
    # Check that the amount of cropland plus change plus urban land plus non arable land is not greater than the area of the cell.
    p.tmp$tmp <- rowSums(p.tmp[,c('urban_prop',
                                  'Prop.NonArable.2010',
                                  'prop_to_predict_from_crop',
                                  'pred_change_cell_prop_crop_decrease')],
                         na.rm = TRUE)
    
    # Adjust predictions if so
    p.tmp[p.tmp$tmp > 1, "pred_change_cell_prop_crop_decrease"] <-
      1 - rowSums(p.tmp[p.tmp$tmp>1,c('urban_prop',
                           'Prop.NonArable.2010',
                           'prop_to_predict_from_crop')],
                  na.rm = TRUE)
    
    # Repeat so that the area can't be negative
    p.tmp$tmp <- rowSums(p.tmp[,c('pred_change_cell_prop_crop_decrease',
                                             'prop_to_predict_from_crop')],
                         na.rm = TRUE)

    
    p.tmp[p.tmp$tmp < 0, "pred_change_cell_prop_crop_decrease"] <-
      0 - p.tmp[p.tmp$tmp < 0, "prop_to_predict_from_crop"]
    
    # Finally, check that fully urbanised areas can't be converted at all
    # If urban proportion is one, then cropland extent cannot change in these cells
    p[p$urban_prop %in% 1, 
      c("pred_to_increase_crop",
        "pred_to_decrease_crop",
        "pred_change_cell_prop_crop_decrease", 
        "pred_change_cell_prop_crop_increase")] <- 0
    
    # Inserting updated estimates for prop increase and prop decrease into original data frame
    p$pred_change_cell_prop_crop_decrease[!is.na(p$pred_change_cell_prop_crop_decrease)] <- 
      p.tmp$pred_change_cell_prop_crop_decrease[!is.na(p$pred_change_cell_prop_crop_decrease)]
    p$pred_change_cell_prop_crop_increase[!is.na(p$pred_change_cell_prop_crop_increase)] <- 
      p.tmp$pred_change_cell_prop_crop_increase[!is.na(p$pred_change_cell_prop_crop_increase)]
    # p$tmp <- p.tmp$tmp
    # Removing p.tmp
    rm(p.tmp)
    
    # p <- subset(p, select = -c(tmp))
    
    # Make sure areas with a suitability of 0 can't have agricultural expansion, but CAN have contraction
    # Don't need these anymore. 
    # JK. We do

    p[which(p$ag_suitability == 0 & p$prop_to_predict_from_crop < .025), 
      "pred_to_increase_crop"] <- 0
    
    # Preventing expansion in cells where GDD == 0 and there is no cropland
    # Expansion can occur if GDD == 0 and there is cropland in a cell
    # Contraction can occur in these cells
    p[which(p$GDD_binary == 0 & p$prop_to_predict_from_crop < .025), 
      "pred_to_increase_crop"] <- 0
    
    # Get the area converted
    p$area_converted_crop_increase <- p$pred_change_cell_prop_crop_increase * cell_area
    p$area_converted_crop_decrease <- p$pred_change_cell_prop_crop_decrease * cell_area
    
    # Round all the numbers to avoid rounding errors
    # Will not meet country-level targets if we don't do this
    p[, c("prop_adj_to_predict_from_crop", "prop_adj_to_predict_from_pasture",
          "prop_to_predict_from_crop", "prop_to_predict_from_pasture",
          "a_pred_to_no_change_crop", "pred_to_increase_crop", "pred_to_decrease_crop", 
          "pred_change_cell_prop_crop_increase", "pred_change_cell_prop_crop_decrease")] <-
      round(p[, c("prop_adj_to_predict_from_crop", "prop_adj_to_predict_from_pasture",
                  "prop_to_predict_from_crop", "prop_to_predict_from_pasture",
                  "a_pred_to_no_change_crop", "pred_to_increase_crop", "pred_to_decrease_crop", 
                  # Re-add this line later
                  "pred_change_cell_prop_crop_increase", "pred_change_cell_prop_crop_decrease")],
            digits = 8)
    
    # Changing so data frames can stack again
    p$ISO_numeric <- as.numeric(p$ISO_numeric)
    return(p)
  } else {
    # Pastureland -----
    # First: probability of change
    # pred_tmp are temporary predictor values. 
    # Gives data frame with probability of either no change, increase, or decrease in pasture extent
    
    iso3s = unique(p$ISO_numeric) %>% as.character(.) %>% .[!is.na(.)]
    
    if(length(iso3s) > 1) {
      dat.out = data.frame()
      for(i in iso3s) {
        pred_tmp <- multinom.predict.2b.trial(dat = p %>% filter(ISO_numeric %in% i),
                                              coefs = coefs.pasture.prob,
                                              vars = predictor_vars_pasture[1:(length(predictor_vars_pasture) - 1)],
                                              factor_vars = 'ISO_numeric')
        dat.out <- 
          rbind(dat.out,
                p %>% filter(ISO_numeric %in% i) %>% cbind(pred_tmp))
      }
      
      p <- 
        dat.out %>%
        dplyr::rename(a_pred_to_no_change_pasture = aNoChange,
                      pred_to_decrease_pasture = Decrease,
                      pred_to_increase_pasture = Increase)
    } else {
      pred_tmp = 
        multinom.predict.2b.trial(dat = p,
                                  coefs = coefs.pasture.prob,
                                  vars = predictor_vars_pasture[1:(length(predictor_vars_pasture) - 1)],
                                  factor_vars = 'ISO_numeric')
      
      # Put predictions from pred_tmp into the dataframe
      p$a_pred_to_no_change_pasture <- pred_tmp$aNoChange
      p$pred_to_decrease_pasture <- pred_tmp$Decrease
      p$pred_to_increase_pasture <- pred_tmp$Increase 
    }

    
    
    # Repeat for increasing models
    # Proportion to increase
    
    # Also looping through isos
    if(length(iso3s) > 1) {
      dat.out = data.frame()
      for(i in iso3s) {
        dat.out <- rbind(dat.out,
                         p %>% filter(ISO_numeric %in% i) %>%
                           mutate(pred_change_cell_prop_pasture_increase = 
                                    as.vector(glm.predict.2b(dat = p %>% filter(ISO_numeric %in% i),
                                                             coefs = coefs.pasture.increase,
                                                             vars = predictor_vars_pasture[1:(length(predictor_vars_pasture) - 1)],
                                                             factor_vars = 'ISO_numeric'))))
        
        
      }
      p <- dat.out
    } else {
      p[,'pred_change_cell_prop_pasture_increase'] <- 
        as.vector(glm.predict.2b(dat = p,
                                 coefs = coefs.pasture.increase,
                                 vars = predictor_vars_pasture[1:(length(predictor_vars_pasture) - 1)],
                                 factor_vars = 'ISO_numeric'))
    }
    
    # This might cause NAs, converting these to 0s
    p$pred_change_cell_prop_pasture_increase[is.na(p$pred_change_cell_prop_pasture_increase)] <- 0
    
    # Check that the amount of pasture, crop, urban land and nonarable land is not greater 
    #     than the area of the cell.
    # Creates temporary cropland to make sure potential future area of cell after crop expansion 
    #     does not exceed 1
    # Creating temp data frame and converting NAs to 0s to avoid issues summing across columns
    p.tmp <- dplyr::select(p, 
                           urban_prop, 
                           Prop.NonArable.2010, 
                           prop_to_predict_from_crop, 
                           prop_to_predict_from_pasture,
                           pred_change_cell_prop_pasture_increase)
    # Converting NAs to 0s
    # Creating new column to adjust values of predictions so that they cannot exceed 1
    p.tmp$tmp <- rowSums(p.tmp[,c('urban_prop',
                                  'Prop.NonArable.2010',
                                  'prop_to_predict_from_crop',
                                  'prop_to_predict_from_pasture',
                                  'pred_change_cell_prop_pasture_increase')],
                         na.rm = TRUE)

    # Adjust predictions if so
    # If area does exceed 1, then adjust so that predicted pasture change does not cause cell value to exceed 1
    p.tmp[p.tmp$tmp > 1, "pred_change_cell_prop_pasture_increase"] <-
      1 - rowSums(p.tmp[p.tmp$tmp > 1,c('urban_prop',
                           'Prop.NonArable.2010',
                           'prop_to_predict_from_crop',
                           'prop_to_predict_from_pasture')],
                  na.rm = TRUE)
    
    # Repeat so that the area can't be negative
    p.tmp$tmp <- rowSums(p.tmp[,
                               c('pred_change_cell_prop_pasture_increase',
                                 'prop_to_predict_from_pasture')],
                         na.rm = TRUE)
    
    p.tmp[p.tmp$tmp < 0, "pred_change_cell_prop_pasture_increase"] <-
      0 - p.tmp[p.tmp$tmp < 0, "prop_to_predict_from_pasture"]
    
    # Finally, check that fully urbanised areas can't be converted at all
    # If urban proportion is one, then cropland extent cannot change in these cells
    p[p$urban_prop == 1 & !is.na(p$urban_prop), c("pred_to_increase_pasture",
                                                  "pred_to_decrease_pasture",
                                                  "pred_change_cell_prop_pasture_increase")] <- 0
    
    # Inserting updated estimates for prop increase into original data frame
    p$pred_change_cell_prop_pasture_increase[!is.na(p$pred_change_cell_prop_pasture_increase)] <- 
      p.tmp$pred_change_cell_prop_pasture_increase[!is.na(p$pred_change_cell_prop_pasture_increase)]
    # Removing p.tmp
    rm(p.tmp)
    
    # Make sure areas with a suitability of 0 can't have agricultural expansion, but CAN have contraction
    # Unless more than .025 portions of a cell currently have crop
    # Don't need these anymore. 
    # JK we do
    # p[which(p$ag_suitability == 0 & p$prop_to_predict_from_crop < .025), 
    #   "pred_to_increase_pasture"] <- 0
  
    # Preventing expansion in cells where GDD == 0 and there is no cropland
    # Expansion can occur if GDD == 0 and there is more than .025 cropland in the cell
    # Contraction can occur in these cells
    # p[which(p$GDD_binary == 0 & p$prop_to_predict_from_crop < .025), 
    #   "pred_to_increase_pasture"] <- 0
    
    # Get the area converted
    p$area_converted_pasture_increase <- p$pred_change_cell_prop_pasture_increase * cell_area
    
    # Round all the numbers to avoid rounding errors
    # Will not meet country-level targets if we don't do this
    p[, c("prop_adj_to_predict_from_crop", "prop_adj_to_predict_from_pasture",
          "prop_to_predict_from_crop", "prop_to_predict_from_pasture",
          # "a_pred_to_no_change_crop", "pred_to_increase_crop", "pred_to_decrease_crop", 
          "pred_to_increase_pasture", "pred_change_cell_prop_pasture_increase")] <- 
    # "pred_change_cell_prop_crop_increase", "pred_change_cell_prop_crop_decrease")] <-
    round(p[, c("prop_adj_to_predict_from_crop", "prop_adj_to_predict_from_pasture",
                "prop_to_predict_from_crop", "prop_to_predict_from_pasture",
                #         "a_pred_to_no_change_crop", "pred_to_increase_crop", "pred_to_decrease_crop", 
                # Re-add this line later
                "pred_to_increase_pasture", "pred_change_cell_prop_pasture_increase")],
          #         "pred_change_cell_prop_crop_increase", "pred_change_cell_prop_crop_decrease")],
          digits = 8)
  }
  # Converting back to numeric so data frames can stack
  p$ISO_numeric <- as.numeric(p$ISO_numeric)
  return(p)
}

# Function to pick cells to convert -----
crop_conv_fun <- function(target, data, cons_planning_prop_to_reduce){
  # Get just country's cells, and those with non-NA values
  df <- filter(data, 
               ISO3 == target$ISO3[1]) %>%
    as.data.frame()
  # Exclude NAs
  df <- df %>% 
    filter(!is.na(area_converted_crop_increase),
           !is.na(area_converted_crop_decrease))
  
  # If conservation land use planning is in place
  # Converting probability for cells to increase in crop or pasture extent to 0
  # In locations where conesrvation actions are implemented
  # These areas can still decrease in extent
  if(conservation_land_planning %in% 'yes' &
     cons_planning_prop_to_reduce %in% 1) {
    # Which rows in the data frame include conservation land use planning
    which_rows_land_planning <- which(df$cons_planning_areas >= 1)
    
    # Now forcing columns to 0
    df[which_rows_land_planning,
       c('pred_to_increase_pasture',
          'pred_to_increase_crop',
          'pred_change_cell_prop_pasture_increase',
          'pred_change_cell_prop_crop_increase',
          'ag_suitability')] <- 0
  }
  
  # Remove cells where ag_suitability == 0 and prop_crop == 0
  # Prevents expansion in cells where ag_suitability == 0 and there is no current crop production
  df.tmp <- df[which(df$ag_suitability != 0 & df$prop_to_predict_from_crop != 0),]
  
  # Set target. Creating vector from data_frame targets
  t <- target$target_temp
  # If troubleshooting, t = target_decade$target[1]
  
  target_orig <- t # set to target_orig to have something to compare to
  
  # Setting logical value. Used in if loop below
  # If true, multiplies values by -1 
  # Does this so that below script works with both positive and negative targets
  switch = FALSE
  
  # Logical statement to only create prob and area vector if target != 0
  if(t > 0) {
    probs <- as.vector(df$pred_to_increase_crop)
    area <- as.vector(df$pred_change_cell_prop_crop_increase * cell_area)
    Row.Index <- as.vector(df$row_index)
    # GDD <- as.vector(df$GDD_binary)
  }
  
  # Logical to make sure target and direction of cell_area are in same direction
  if(t < 0) {
    probs <- as.vector(df$pred_to_decrease_crop)
    area <- -1 * as.vector(df$pred_change_cell_prop_crop_decrease * cell_area)
    Row.Index <- as.vector(df$row_index)
    # GDD <- as.vector(df$GDD_binary)
    t = t * -1
    switch = TRUE
  }
  
  # Drop cells with ag_suitability == 0 and prop_crop == 0
  
  # Drop zero probability cells and NA cells
  # Don't want them to be picked
  area <- area[probs > 0 & !is.na(probs)]
  Row.Index <- Row.Index[probs > 0 & !is.na(probs)]
  # GDD <- GDD[probs > 0]
  probs <- probs[probs > 0 & !is.na(probs)]
  
  # Also drop zero area cells
  # Creates shorter vectors and (should) make code faster
  Row.Index <- Row.Index[area!= 0 & !is.na(area)]
  probs <- probs[area!= 0 & !is.na(area)]
  # GDD <- GDD[area!= 0]
  area <- area[area!= 0 & !is.na(area)]
  
  
  # Counter to see if we're double counting cells
  n1 = 0
  # Selecting all cells if sum of cells with correct direction of change is less than target
  if(sum(area[which(area > 0)]) <= t) {
    # And adding row.nums to cells_tmp
    row.num.list <- 1:(length(area))
    # Select all cells with positive values
    cells_tmp <- row.num.list[which(area > 0)]
    # Creating list used to 
    # And removing first value of cell
    # Adjust target
    t = t - sum(area[which(area>0)])
    n1 = n1 + length(which(area > 0))
  } else {
    
    # Getting estimate of number of cells needed
    n.cells = round(t/mean(area))
    if(n.cells > length(Row.Index)) {
      Row.Index <- Row.Index[which(area > 0)]
      probs <- probs[which(area > 0)]
      # GDD <- GDD[which(area > 0)]
      area <- area[which(area > 0)]
    }
    
    # n.cells can't be negative
    # Removing cells where area < 0
    # Avoids getting stuck in a while loop
    if(n.cells < 0) {
      Row.Index <- Row.Index[which(area > 0)]
      probs <- probs[which(area > 0)]
      # GDD <- GDD[which(area > 0)]
      area <- area[which(area > 0)]
      # Updating estimate of cells needed for conversion
      n.cells = round(t/mean(area) * .95)
    }
    # Creating list of place in the Row.Index list of cells that have been selected
    # Is used so that these cells are not selected again
    row.num.list <- 1:(length(area) + 1)
    # Need to append values onto probs so that vectors is same length as row.num.list
    # Not doing this causes the fastSampleReject function to crash
    probs <- c(probs, 0)
    
    # Creating list of row indices of cells that were selected for a change in cropland extent
    # First value of cells_tmp is equivalent to last value of row.num.list
    # This value will be dropped
    cells_tmp <- row.num.list[length(row.num.list)]
    # Selecting many cells to meet target rapidly
    while(t > 3) {
      # Getting estimate of number of cells needed
      n.cells = round(t/mean(area) * .95)
      # Adding in exception if n.cells is greater than the number of remaining cells that can be picked
      # Just assuming all of the remaining cells are selected
      # This could result in an overestimate of crop conversion...
      # But this will almost always be a small overestimate
      if(n.cells > length(area[-cells_tmp])) {
        n.cells <- length(area[-cells_tmp])
      }
      
      # Selecting rows
      row.nums <- fastSampleReject(all = row.num.list[-cells_tmp], n = n.cells, w = probs[-cells_tmp])
      # Cutting off cells if sum of selected cells is greater than target
      # Put in a while loop so that it loops below code until condition is satisfied
      while(sum(area[row.nums]) > t) {
        # If sum area[row.nums] > t, then reduce number of cells selected and resample until this is met
        n.cells = round(n.cells * .9)
        # Reselect cells
        row.nums <- fastSampleReject(all = row.num.list[-cells_tmp], n = n.cells, w = probs[-cells_tmp])
      }
      
      # Appending selected rows 
      cells_tmp <- c(cells_tmp, row.nums)
      # Adjusting target
      t = t - sum(area[row.nums])
      # Adding counter to verify we're not double-selecting cells
      n1 = n1 + length(row.nums)
    }
    
    # Selecting individual cells to approach target slowly
    while(t > 0) {
      # Selecting rows
      # Need to select first object out because everything else is an NA
      # Putting this in an if statement to prevent the code from freaking out
      if(length(probs[-cells_tmp] < 1) | length(row.num.list[-cells_tmp] < 1)) {
        t = 0
      } else {
        row.nums <- fastSampleReject(all = row.num.list[-cells_tmp], n = 1, w = probs[-cells_tmp])[1]
        # Appending selected rows 
        cells_tmp <- c(cells_tmp, row.nums)
        # Adjusting target
        t = t - area[row.nums]
        # Adding counter to verify we're not double-selecting cells
        n1 = n1 + 1
      }
    }
  }
  
  # Resetting switch
  # Also resetting values of area and t
  # Because these were multiplied by -1 above
  if(switch == TRUE) {
    t = t * -1
    area = area * -1
    switch = FALSE
  }
  # Dropping first value from cells_tmp because this is a place holder
  # But only need to do this if counter is not equal to length of cells_tmp
  if(length(cells_tmp) != n1) {
    cells_tmp <- cells_tmp[-1]
  }
  
  # Getting vector of row.indices to tell us what cells to convert
  cells_tmp <- Row.Index[cells_tmp]
  # Data frame containing all of the cells converted, how much area has been converted
  # Data frame can be quite large. A single row for every cell picked, but single value for area converted
  # Data frame can also contain no cells. Need an if statement if this is the case
  if(length(cells_tmp != 0)) {
    out <- data.frame(index = cells_tmp,
                      area_converted = target_orig - t,
                      target_remaining = t,
                      row.names = 1:length(cells_tmp))
    
  }
  # And creating a data frame if no cells were picked for conversion
  # Need to do this so that script continues running
  # Picking an arbitrarily large index for cells_tmp
  # R will return the entirevectors if the index is outside of their range
  if(length(cells_tmp) == 0) {
    out <- data.frame(index = 1e12,
                      area_converted = 0,
                      target_remaining = t,
                      row.names = 1e12)
  }
  # I have no idea why I need to rename things, but this works...
  names(out) <- c("index", "area_converted","target_remaining")
  return(out)
}

# Function to pick cells to convert for pasture (could probably combine this with the one above) ----
pasture_conv_fun <- function(target, data, cons_planning_prop_to_reduce){
  # Filter to country's cells ----
  # Get just country's cells, and those with non-NA values
  df <- filter(data, 
               ISO3 == target$ISO3[1]) %>%
    as.data.frame()
  # Exclude NAs
  df <- df %>% 
    filter(!is.na(area_converted_pasture_increase))
  #!is.na(area_converted_pasture_decrease))
  
  # If conservation land use planning is in place
  # Converting probability for cells to increase in crop or pasture extent to 0
  # In locations where conesrvation actions are implemented
  # These areas can still decrease in extent
  if(conservation_land_planning %in% 'yes' &
     cons_planning_prop_to_reduce %in% 1) {
    # Which rows in the data frame include conservation land use planning
    which_rows_land_planning <- which(df$cons_planning_areas >= 1)
    
    # Now forcing columns to 0
    df[which_rows_land_planning,
       c('pred_to_increase_pasture',
          'pred_to_increase_crop',
          'pred_change_cell_prop_pasture_increase',
          'pred_change_cell_prop_crop_increase',
          'ag_suitability')] <- 0
  }
  
  # We DON'T remove cells with a suitability of 0 for pasture
  # Set target. ----
  # Creating vector from data_frame targets
  t <- target$target
  # If troubleshooting, t = target_decade$target[1]
  target_orig <- t # set to target_orig to have something to compare to
  
  # Setting logical value. Used in if loop below
  # If true, multiplies values by -1 
  # Does this so that below script works with both positive and negative target
  switch = FALSE
  
  # Logical statement to only create prob and area vector if target != 0
  if(t > 0) {
    probs <- as.vector(df$pred_to_increase_pasture)
    area <- as.vector(df$pred_change_cell_prop_pasture_increase * cell_area)
    Row.Index <- as.vector(df$row_index)
    # GDD <- as.vector(df$GDD_binary)
  }
  
  # Drop zero probability cells
  # Don't want them to be picked
  area <- area[probs > 0 & !is.na(probs)]
  Row.Index <- Row.Index[probs > 0 & !is.na(probs)]
  # GDD <- GDD[probs > 0]
  probs <- probs[probs > 0 & !is.na(probs)]
  
  # Also drop zero area cells
  # Creates shorter vectors and (should) make code faster
  Row.Index <- Row.Index[area!= 0 & !is.na(area)]
  probs <- probs[area!= 0 & !is.na(area)]
  # GDD <- GDD[area!= 0]
  area <- area[area!= 0 & !is.na(area)]
  
  # Counter to see if we're double counting cells
  n1 = 0
  # Selecting all cells if sum of cells with correct direction of change is less than target
  if(sum(area[which(area > 0)]) <= t) {
    # And adding row.nums to cells_tmp
    row.num.list <- 1:(length(area))
    # Select all cells with positive values
    cells_tmp <- row.num.list[which(area > 0)]
    # Creating list used to 
    # And removing first value of cell
    # Adjust target
    t = t - sum(area[which(area>0)])
    n1 = n1 + length(which(area > 0))
  } else {
    
    # Getting estimate of number of cells needed
    n.cells = round(t/mean(area))
    if(n.cells > length(Row.Index)) {
      Row.Index <- Row.Index[which(area > 0)]
      probs <- probs[which(area > 0)]
      # GDD <- GDD[which(area > 0)]
      area <- area[which(area > 0)]
    }
    
    
    # n.cells can't be negative
    # Removing cells where area < 0
    # Avoids getting stuck in a while loop
    if(n.cells < 0 | median(area) < 0) {
      Row.Index <- Row.Index[which(area > 0)]
      probs <- probs[which(area > 0)]
      # GDD <- GDD[which(area > 0)]
      area <- area[which(area > 0)]
      # Updating estimate of cells needed for conversion
      n.cells = round(t/mean(area) * .95)
    }
    
    # Creating list of place in the Row.Index list of cells that have been selected
    # Is used so that these cells are not selected again
    row.num.list <- 1:(length(area) + 1)
    # Need to append values onto probs so that vectors is same length as row.num.list
    # Not doing this causes the fastSampleReject function to crash
    probs <- c(probs, 0)
    
    # Creating list of row indices of cells that were selected for a change in cropland extent
    # First value of cells_tmp is equivalent to last value of row.num.list
    # This value will be dropped
    cells_tmp <- row.num.list[length(row.num.list)]
    # Selecting many cells to meet target rapidly
    while(t > 3) {
      # Getting estimate of number of cells needed
      n.cells = round(t/mean(area[-cells_tmp]) * .95)
      
      # Adding in exception if n.cells is greater than the number of remaining cells that can be picked
      # Just assuming all of the remaining cells are selected
      if(n.cells > length(area[-cells_tmp])) {
        n.cells <- length(area[-cells_tmp])
      }
      # Selecting rows
      row.nums <- fastSampleReject(all = row.num.list[-cells_tmp], n = n.cells, w = probs[-cells_tmp])
      # Cutting off cells if sum of selected cells is greater than target
      # Put in a while loop so that it loops below code until condition is satisfied
      while(sum(area[row.nums]) > t) {
        # If sum area[row.nums] > t, then reduce number of cells selected and resample until this is met
        n.cells = round(n.cells * .95)
        # Reselect cells
        row.nums <- fastSampleReject(all = row.num.list[-cells_tmp], n = n.cells, w = probs[-cells_tmp])
      }
      
      # Appending selected rows 
      cells_tmp <- c(cells_tmp, row.nums)
      # Adjusting target
      t = t - sum(area[row.nums])
      # Adding counter to verify we're not double-selecting cells
      n1 = n1 + length(row.nums)
    }
    
  }
  
  # Resetting switch
  # Also resetting values of area and t
  # Because these were multiplied by -1 above
  if(switch == TRUE) {
    t = t * -1
    area = area * -1
    switch = FALSE
  }
  # Dropping first value from cells_tmp because this is a place holder
  # But only need to do this if counter is not equal to length of cells_tmp
  if(length(cells_tmp) != n1) {
    cells_tmp <- cells_tmp[-1]
  }
  
  # Getting vector of row.indices to tell us what cells to convert
  cells_tmp <- Row.Index[cells_tmp]
  # Data frame containing all of the cells converted, how much area has been converted
  # Data frame can be quite large. A single row for every cell picked, but single value for area converted
  # Data frame can also contain no cells. Need an if statement if this is the case
  if(length(cells_tmp != 0)) {
    out <- data.frame(index = cells_tmp,
                      area_converted = target_orig - t,
                      target_remaining = t,
                      row.names = 1:length(cells_tmp))
    
  }
  # And creating a data frame if no cells were picked for conversion
  # Need to do this so that script continues running
  # Picking an arbitrarily large index for cells_tmp
  # R will return the entirevectors if the index is outside of their range
  if(length(cells_tmp) == 0) {
    out <- data.frame(index = 1e12,
                      area_converted = 0,
                      target_remaining = t,
                      row.names = 1e12)
  }
  # I have no idea why I need to rename things, but this works...
  names(out) <- c("index", "area_converted","target_remaining")
  return(out)
}

# Function to run through each year, convert cells and update predictions ----
year_fun <- function(x){
  mod_df <- model_data_current

  # Really 5-yr time intervals. Decs because originally forecasting at 10-yr intervals
  decs <- gsub("target_", 
               "",
               names(mod_df)[grep("target_", names(mod_df))])
  
  # Getting years for urban expansion
  urban.years <- 
    decs %>%
    gsub(".*_","",.)
  
  # Set up list for messages ----
  # Messages because warnings
  messages.list <- list()
  # Loop through decades -----
  # for(i in 1:1){
  for(i in 1:length(decs)){
    # Set up the decade to use
    dec <- decs[i]
    
    #cat('Decades:',decs[i])
    # Residual target for crop and pasture based on urban expansion
    # Assuming it expands into non arable or other lands before crop and pasture
    # e.g. if urban_tmp + crop + pasture > 1, then need to reduce crop and pasture
    # otherwise, reduce non.arable and other equally
    
    # Doing this in 3 steps
    
    # (1) Which cells have urban expansion
    # Doing this to check if we need to worry about urban expansion
    cells.expansion <- which(mod_df[,paste0('UrbanExtent',urban.years[i])] > mod_df[,'urban_prop'])
    
    # If yes
    if(length(cells.expansion) >= 1) {
      # (2) For cells where urban + crop + pasture < 1, adjust non arable and other lands
      tmp <- rowSums(mod_df[,c(paste0('UrbanExtent',urban.years[i]),'prop_to_predict_from_crop','prop_to_predict_from_pasture')], na.rm = TRUE)
      tmp.index.urb.crop.past <- which(tmp <= 1 & !is.na(mod_df$ISO_numeric))
      tmp.index.urb.crop.past <- cells.expansion[cells.expansion %in% tmp.index.urb.crop.past]
      
      # Another if statement to hopefully save time
      if(length(tmp.index.urb.crop.past) >= 1) {
        # Land needed to lose
        tmp.loss.urb.crop.past <- 
          mod_df[tmp.index.urb.crop.past,c('Prop.NonArable.2010','other_2010')] / # This gets ratio of nonarable and other lands
          rowSums(mod_df[tmp.index.urb.crop.past,c('Prop.NonArable.2010','other_2010')], na.rm = TRUE) *
          (mod_df[tmp.index.urb.crop.past,paste0('UrbanExtent',urban.years[i])] - mod_df[tmp.index.urb.crop.past,'urban_prop']) # How much land is lost to urban expansion

        
      }
   cat('L917\n')   
      
      # (3) For cells where urban + crop + pasture > 1, adjust non arable and other lands first
      # And then adjust crop and pasture
      tmp <- rowSums(mod_df[,c(paste0('UrbanExtent',urban.years[i]),'prop_to_predict_from_crop','prop_to_predict_from_pasture')], na.rm = TRUE)
      tmp.index.all <- which(tmp > 1 & !is.na(mod_df$ISO_numeric))
      tmp.index.all <- cells.expansion[cells.expansion %in% tmp.index.all]
      
      # If statement to save time
      if(length(tmp.index.all) >= 1) {
        # This is calculated as how much crop and pasture is left
        # Noting that non-arable and other lands go to 0 automatically
        tmp.loss.all <- 
          mod_df[tmp.index.all,c('prop_to_predict_from_crop','prop_to_predict_from_pasture')] / # This gets ratio of existing crop and pasture
          rowSums(mod_df[tmp.index.all,c('prop_to_predict_from_crop','prop_to_predict_from_pasture')], na.rm = TRUE) *
          (1 - mod_df[tmp.index.all,paste0('UrbanExtent',urban.years[i])]) # How much land is left for crop and pasture
      }
    }
    
    
    # And updating urban area
    # And getting residual target for crop and pasture
    mod_df$urban_prop <-
      mod_df[,paste0('UrbanExtent',urban.years[i])] 
    
    # And adjusting other land covers
    # First if urban + crop + pasture < 1
    #Making sure it exists - need an exception for certain countries without project urban expansion
    if(exists('tmp.index.urb.crop.past')) {
      mod_df[tmp.index.urb.crop.past,c('Prop.NonArable.2010','other_2010')] <-
        mod_df[tmp.index.urb.crop.past,c('Prop.NonArable.2010','other_2010')] -
        tmp.loss.urb.crop.past
    }
    
    # Second if urban + crop + pasture > 1
    # Likewise, ned exceptions for countries without urban expansion
    if(exists('tmp.index.all')) {
      mod_df[tmp.index.all,c('Prop.NonArable.2010','other_2010')] <- 0
      
      # Calculating land targets from urban expansion
      # These are in sq km
      crop.target.tmp <- 
        sum(mod_df[tmp.index.all,c('prop_to_predict_from_crop')] -
              tmp.loss.all[,'prop_to_predict_from_crop'], na.rm = TRUE) *
        cell_area
      
      pasture.target.tmp <- 
        sum(mod_df[tmp.index.all,c('prop_to_predict_from_pasture')] -
              tmp.loss.all[,'prop_to_predict_from_pasture'], na.rm = TRUE) *
        cell_area
      
      # And updating values in the data from
      mod_df[tmp.index.all,c('prop_to_predict_from_crop','prop_to_predict_from_pasture')] <-
        tmp.loss.all
    } else {
      # No ag land displaced from urban expansion
      crop.target.tmp <- 0
      pasture.target.tmp <- 0
    }
    cat('L976\n')

    # And a quick check to make sure these sum correctly
    # min(rowSums(mod_df[,c('prop_to_predict_from_crop','prop_to_predict_from_pasture','urban_prop','Prop.NonArable.2010','other_2010')]), na.rm = TRUE)
    # max(rowSums(mod_df[,c('prop_to_predict_from_crop','prop_to_predict_from_pasture','urban_prop','Prop.NonArable.2010','other_2010')]), na.rm = TRUE)
    
    # And updating GDD_binary
    mod_df$GDD_binary <- mod_df[,paste0('GDD_Binary_',urban.years[i])]
    
    # Set demand to correct decade
    mod_df$change_demand_crop <- mod_df[, paste("prop_change", dec, sep = "_")] + 1
    
    # Get targets -----
    target_decade <- as.data.frame(mod_df[!duplicated(mod_df$ISO3),
                                          c("ISO3", paste("target", dec, sep = "_"))])
    names(target_decade)[2] <- "target"
    target_decade <- target_decade[!is.na(target_decade$ISO3),]
    
   cat(decs[i],'L995\n') 
    # And manually updating target_decade$target if needed
    # For some reasons this can result in an NA for some, but not all, countries
    if(is.na(target_decade$target)) {
      # First getting value for appropriate decade
      tmp.target.dec <- targets[targets$ISO3 == target_decade$ISO3,paste('target_',
                                                                         dec,
                                                                         sep = '')]
      
      # And second dropping NA value
      tmp.target.dec <- tmp.target.dec[!is.na(tmp.target.dec)]
      cat('L1005\n')
      
      if(is.finite(max(tmp.target.dec))) {
        target_decade$target <- tmp.target.dec
      } else {
        target_decade$target[is.na(target_decade$target)] <- rep(targets[targets$ISO3 %in% iso3c[c], paste0('target_',dec)], length(is.na(target_decade$target)))
      }
    }
    cat(decs[i],'L1013\n')
    
    # And adding crop target
    target_decade$target <-
      target_decade$target +
      crop.target.tmp
    
    # Adding checker for whether conservation land needs to be reduced
    cons_reduce_true <- 0 # 0 for false
    
    # Adding if statement, to add crop target for land that is displaced from conservation land use planning
    # Only implemented if...
    # (1) Conservation land use planning is used as a scenario
    # (2) Year is after 2020 (e.g. this starts being implemented in 2020 to 2025)
    if(conservation_land_planning %in% 'yes' &
       grepl('_2025|203|204',decs[i])) {

cat('Cons Land Planning Check:', decs[i])
	    #cat('Cons Land Planning Check 1\n')
      decs_remaining <- decs[grepl('_2025|203|204',decs)]
      decs_remaining <- str_extract(decs_remaining,'_[0-9]{4,4}\\b')
      
      our_dec <- decs[i]
      our_dec <- str_extract(our_dec,'_[0-9]{4,4}\\b')
      
      decs_remaining <- decs_remaining[decs_remaining >= our_dec]
      decs_remaining <- length(decs_remaining)
      #cat('Cons Land Planning Check 2\n',our_dec,decs_remaining,'\n')
      
      # Getting which rows need to be replaced
      which_rows_cont_areas <- which(mod_df$contraction_areas >= 1)
      
      # Displaced ag areas
      displaced_ag_areas <-
        mod_df %>%
        filter(contraction_areas >= 1) %>%
        filter(ISO3 %in% iso3c[c]) %>%
        dplyr::group_by('displaced_ag_area') %>%
        dplyr::summarise(displaced_crop = sum(prop_to_predict_from_crop, na.rm = TRUE), # This gets total area
                         displaced_past = sum(prop_to_predict_from_pasture, na.rm = TRUE)) %>%
        mutate(displaced_crop = displaced_crop / decs_remaining, # And amount removed in this time period
               displaced_past = displaced_past / decs_remaining)
      
      # Total displaced areas in square kilometers
      displaced_ag_areas_crop <- displaced_ag_areas$displaced_crop
      displaced_ag_areas_pasture <- displaced_ag_areas$displaced_past

      cat(exists('displaced_ag_areas_crop'),'\n')
      cat(is.numeric(displaced_ag_areas_crop),'\n')
      cat(displaced_ag_areas_crop,'\n')
      displaced_ag_areas_crop
      print(displaced_ag_areas_crop)
      if(max(displaced_ag_areas_crop)<0) {displaced_ag_areas_crop <- 0}
      if(max(displaced_ag_areas_pasture)<0) {displaced_ag_areas_pasture <- 0}
      cat('Cons Check 3:',decs[i],displaced_ag_areas_crop,displaced_ag_areas_pasture,'\n')
      # Reducing amount of ag land in the cell
      # This is done as current_amount * (1 - 1/remaining decades)
      # So e.g., if there are 6 decades left, then this will reduce amount of ag in the cell by 1/6th
      # If 4 decades left, it will reduce amoung of ag in the cell by 1/4th, etc
      # This is updated here, because prop_to_predict_from_crop, and prop_to_predict_from_pasture are later used to update the amount of ag land in the 'current' decade
      mod_df <-
        mod_df %>%
        mutate(prop_to_predict_from_crop = ifelse(contraction_areas >= 1,
						  prop_to_predict_from_crop * (1 - 1/decs_remaining),
						  prop_to_predict_from_crop),
               prop_to_predict_from_pasture = ifelse(contraction_areas >= 1,
						     prop_to_predict_from_pasture * (1 - 1/decs_remaining),
						     prop_to_predict_from_pasture))
     cat('Target:',target_decade$target,'\n')
    cat('Target Remaining:',target_remaining,'\n')
    cat('Data Frame:',is.data.frame(displaced_ag_areas_crop))
    #cat('Null::',is.null(displaced_ag_areas_crop))
    #cat('Numeric:',is.numeric(displaced_ag_areas_crop))
    #cat('Character:',is.character(displaced_ag_areas_crop))
   cat('Displaced Ag Areas:',displaced_ag_areas_crop,'\n') 
      # And updating cropland target, for elsewhere in the country
      target_decade$target <-
        target_decade$target +
        displaced_ag_areas_crop

cat('New Target:',target_decade$target,'\n')
      
      # And updating checker
      cons_reduce_true <- 1 # Updating to 1 for true
#cat('Cons Check Completed\n')
    }
    
    
    # Put warnings into a list
    messages.list[[dec]] <- list()
    
    # Adding remaining target from previous decade to this decade
    if(i != 1) {
      target_decade$target <- target_decade$target + target_remaining
    }
    
    # Get the absolute value of the targets, so I can run a while() loop
    # I don't think is needed anymore
    target_decade$abs.target <- abs(target_decade[,"target"])
    
    # Split target into parts -----
    # Looping through the below while loop n amount of times
    # Each loop hits 1/n part of the total target fo rthe time interval
    # This makes it more likely to for cropland expansion to be concentrated
    # Setting n
    n = 5
    # Logical to fix n if it is less than one
    if(n < 1) {
      n <- 1
    }
    # Setting split target (previously "target_temp" - this is so that 1/5th of the 5-year target is met each 'year'
    # Now doing this outside of the for() loop
    target_decade$split_target <- target_decade$target / n
    # Reset remaining target
    target_remaining <- 0
    
    # Re-run predictions for the decade -----
    # NB. We are doing this only once at the start of each decade, rather than in every
    #     loop of the while() loop, or y-loop.
    # This is because it reflects the logic of our model: we are predicting how things
    #     change in a five-year period, based on the situation at the start of that period,
    #     not how it changes throughout the period.
    # First get the adjacent agricultural proportions
    mod_df <- adj.ag.f(m = mod_df,
                       r = rows,
                       c = cols,
                       crop = "crop")  %>%
      as.data.frame()
    mod_df <- adj.ag.f(m = mod_df,
                       r = rows,
                       c = cols,
                       crop = "pasture")  %>%
      as.data.frame()
    
    # Then get new predictions
    na.data <- mod_df[is.na(mod_df$ISO3),]
    mod_df <- mod_df[!is.na(mod_df$ISO3),]
    # Dropping columns from previous predictions
    # 
    
    mod_df <- pred_fun(mod_df, crop = "crop")
    mod_df <- pred_fun(mod_df, crop = "pasture")
    mod_df <- bind_rows(mod_df, na.data)
    mod_df <- mod_df[order(mod_df$row_index),]
    
    # Run through split target ----
    for(y in 1:n) {
      #cat('Decades:',i,'\n',
      #    'Ag Loop:',y,'\n')
      # Adding counter to loop through multiple rounds of conversion
      counter <- 1
      # If statement to avoid looping over countries that do not have any cropland expansion
      if(target_decade$abs.target < .1) {
        # Set counter to 19 to avoid while loop below
        counter <- 19
        # And update target_remaining so that it does not catch on above logical statement
        # Not doing this results in an error
        target_remaining <- target_decade$target
        
        
        
        # And updating values for extent conversion and absolute crop for countries that have no expansion
        mod_df[[paste("prop_crop", dec, sep = "_")]] <- mod_df$prop_to_predict_from_crop
      }  
      # Get a temporary target: the portion of the overall target that we're trying to meet 
      #     (at the moment 1/5) plus any remaining target from the previous loop
      target_decade$target_temp <- target_decade$split_target + target_remaining

      # Pick cells for crop expansion ----
      while(abs(target_decade$target_temp) > 1 & counter <= 15){
        # Make sure that crop predictions can't be too high or low -----
        # This is lifted from pred_fun()
        mod_df$tmp <- rowSums(mod_df[,c('urban_prop',
                                      'Prop.NonArable.2010',
                                      'prop_to_predict_from_crop',
                                      'pred_change_cell_prop_crop_increase')],
                             na.rm = TRUE)

        # Adjust predictions if so
        # If area does exceed 1, then adjust so that predicted cropland change does not cause cell value to exceed 1
        mod_df[mod_df$tmp > 1, "pred_change_cell_prop_crop_increase"] <-
          1 - rowSums(mod_df[mod_df$tmp>1,c('urban_prop',
                                'Prop.NonArable.2010',
                                'prop_to_predict_from_crop')],
                      na.rm = TRUE)

        
        # Repeat so that the area can't be negative
        mod_df$tmp <- rowSums(mod_df[,c('pred_change_cell_prop_crop_increase',
                                        'prop_to_predict_from_crop')],
                              na.rm = TRUE)
        
        mod_df[mod_df$tmp < 0 & !is.na(mod_df$tmp), "pred_change_cell_prop_crop_increase"] <-
          0 - mod_df[mod_df$tmp < 0 & !is.na(mod_df$tmp), "prop_to_predict_from_crop"]
        
        # Repeat for decreasing models
        # Check that the amount of cropland plus change plus urban land plus non arable land is not greater than the area of the cell.
        mod_df$tmp <- rowSums(mod_df[,c('urban_prop',
                                        'Prop.NonArable.2010',
                                        'prop_to_predict_from_crop',
                                        'pred_change_cell_prop_crop_decrease')],
                              na.rm = TRUE)
        
        # Adjust predictions if so
        mod_df[mod_df$tmp > 1, "pred_change_cell_prop_crop_decrease"] <-
          1 - rowSums(mod_df[mod_df$tmp>1, c('urban_prop',
                                'Prop.NonArable.2010',
                                'prop_to_predict_from_crop')],
                      na.rm = TRUE)
        
        # Repeat so that the area can't be negative
        mod_df$tmp <- rowSums(mod_df[,c('pred_change_cell_prop_crop_decrease',
                                        'prop_to_predict_from_crop')],
                              na.rm = TRUE)
        
        mod_df[mod_df$tmp < 0 & !is.na(mod_df$tmp), "pred_change_cell_prop_crop_decrease"] <-
          0 - mod_df[mod_df$tmp < 0 & !is.na(mod_df$tmp), "prop_to_predict_from_crop"]
        
        # Pick cells ----
        if(nrow(target_decade[target_decade$target_temp > 0,]) > 0){
		#cat('Crop Cells Expansion:',i,y,target_decade$target,'\n')
          cells_expansion <- crop_conv_fun(target = target_decade, 
                                           data = mod_df,
                                           cons_planning_prop_to_reduce = cons_reduce_true)
        } else {
          cells_expansion <- data.frame(index = NULL,
                                        area_converted = NULL)
        }
        if(nrow(target_decade[target_decade$target_temp < 0,]) > 0){
		#cat('Crop Cells Contraction:',i,y,target_decade$target,'\n')
          cells_contraction <- crop_conv_fun(target = target_decade, 
                                             data = mod_df,
                                             cons_planning_prop_to_reduce = cons_reduce_true) 
        } else {
          cells_contraction <- data.frame(index = NULL,
                                          area_converted = NULL)
        }
        
        # Convert chosen cells ----
        # Expansion
        mod_df[mod_df$row_index %in% cells_expansion$index, "prop_to_predict_from_crop"] <- 
          mod_df[mod_df$row_index %in% cells_expansion$index, "prop_to_predict_from_crop"] + 
          mod_df[mod_df$row_index %in% cells_expansion$index, "pred_change_cell_prop_crop_increase"]
        # Contraction
        mod_df[mod_df$row_index %in% cells_contraction$index, "prop_to_predict_from_crop"] <- 
          mod_df[mod_df$row_index %in% cells_contraction$index, "prop_to_predict_from_crop"] + 
          mod_df[mod_df$row_index %in% cells_contraction$index, "pred_change_cell_prop_crop_decrease"]
        
        # Update the change in the previous time-step
        mod_df[mod_df$row_index %in% cells_expansion$index, "prop_to_predict_from_cell_delta_crop"] <- 
          mod_df[mod_df$row_index %in% cells_expansion$index, "pred_change_cell_prop_crop_increase"]
        # Contraction
        mod_df[mod_df$row_index %in% cells_contraction$index, "prop_to_predict_from_cell_delta_crop"] <- 
          mod_df[mod_df$row_index %in% cells_contraction$index, "pred_change_cell_prop_crop_decrease"]
        
        # Getting remaining target. e.g. land that could not be converted ----
        # This will be added onto target for future decades
        if(target_decade$target > 0) {
          target_remaining = unique(cells_expansion$target_remaining)
        }
        if(target_decade$target < 0) {
          target_remaining = unique(cells_contraction$target_remaining)
        }
        if(target_decade$target == 0) {
          target_remaining = 0
        }
        
        # Adjust targets -----
        cells <- bind_rows(cells_contraction, 
                           cells_expansion)
        target_decade$converted <- cells$area_converted[1]
        # Take the area converted off the target
        # In an if loop because target_decade$converted is undefined if target_decade$target == 0
        if(target_decade$target != 0) {
          target_decade$target_temp <- target_decade$target_temp - target_decade$converted
        }
        
        # }
        # Adding in logical to prevent script from looping through while loop if there are no cells left to expand
        # Also prevents from looping if target has been hit
        if(nrow(cells) == 0) {
          counter <- 16
        } else {if(abs(target_decade$split_target - target_decade$converted) < 1) {
          counter <- 17
        } else {
          counter <- counter + 1
        }
        }
        # Below bracket marks end of while loop that meets 1/n proportion of time interval's crop target
      }
      
      # Pasture target setting -----
      # # Work out how much pasture has been lost:
      # # Area available for pasture is the area without urban or cropland
      mod_df$prop_for_pasture <- 0
      # Getting row index of cells that changed in crop extent

      # Calculating amount for pasture
      mod_df$prop_for_pasture <- 1 - rowSums(mod_df[,c('urban_prop',
                                                                  'prop_to_predict_from_crop',
                                                                  'Prop.NonArable.2010')],
                                                        na.rm = TRUE)
      
      # Calculating amount of pasture lost
      mod_df$pasture_loss_tmp <- 0
      mod_df$pasture_loss_tmp <- mod_df$prop_to_predict_from_pasture -
        mod_df$prop_for_pasture
      
      # Dropping value
      mod_df <- dplyr::select(mod_df,
                              -prop_for_pasture)
      
      # Pasture is not necessarily gained in a cell if cropland extent decreases
      # Negative values correspond with pasture area being gained by a cell
      # So Limiting range of pasture loss values from 0 to 1
      # Negative values in pasture_loss_tmp would indicate a net increase in pasture extent in the cell
      mod_df$pasture_loss_tmp <- squish(mod_df$pasture_loss_tmp, 
                                                   range = c(0,1))
      
      # Adjust the amount of pastureland so it is equal to that available
      mod_df$prop_to_predict_from_pasture <-
        mod_df$prop_to_predict_from_pasture - mod_df$pasture_loss_tmp
      # 
      
      # Round this value. NOTE: I am rounding this to one fewer digit (7 vs 8). This is 
      #     because taking two 8-digit decimals from each other and rounding it to 8 d.p. 
      #     can still give a floating point error
      mod_df$prop_to_predict_from_pasture <- round(mod_df$prop_to_predict_from_pasture, 
                                                   digits = 7)
      
      # # Sum up the total area lost for each country
      # Need to make exception for first loop to avoid repeatedly adding
      # Pasture displaced by urban expansion
      if(y %in% 1) {
        pasture_target <- data.frame(ISO3 = target_decade$ISO3[1],
                                     target = sum(mod_df$pasture_loss_tmp * cell_area,
                                                  na.rm = TRUE)) %>%
          mutate(target = target + pasture.target.tmp) %>%
          mutate(abs.target = abs(target))
      } else {
        pasture_target <- data.frame(ISO3 = target_decade$ISO3[1],
                                     target = sum(mod_df$pasture_loss_tmp * cell_area,
                                                  na.rm = TRUE)) %>%
          mutate(abs.target = abs(target))
      }
      
      
      # Manually updating values of pasture_target in the case that this data frame is not created
      # Need to make exception for first loop to avoid repeatedly adding
      # Pasture displaced by urban expansion
      if(nrow(pasture_target) == 0) {
        if(y %in% 1) {
          pasture_target <- data.frame(ISO3 = target_decade$ISO3[1],
                                       target = 0,
                                       abs.target = 0) +
            pasture.target.tmp
        } else {
          pasture_target <- data.frame(ISO3 = target_decade$ISO3[1],
                                       target = 0,
                                       abs.target = 0)
        }
        
      }
      
      # And updating pasture target, based on amount of ag land displaced by conservation land use planning
      if(conservation_land_planning %in% 'yes' &
         grepl('_2025|203|204',decs[i])) {
        pasture_target <-
          pasture_target %>%
          mutate(target = target + displaced_ag_areas_pasture / n) %>%
          mutate(abs.target = abs(target))
      }
      
      # If statement to avoid looping over countries that do not have any cropland expansion
      if(pasture_target$abs.target < .1) {
        # Setting counter to 19
        counter <- 19
        # and updating pasture extent in next time period if no cells are changed
        # And updating values for extent conversion and absolute crop for countries that have no expansion
        mod_df[[paste("prop_pasture", dec, sep = "_")]] <- mod_df$prop_to_predict_from_pasture
      } else {counter <- 1}
      
      # Pick cells for pasture expansion ----
      while(abs(pasture_target$target) > 1 & counter <= 15){
        # Make sure that pasture predictions can't be too high or low ----
        mod_df$tmp <- rowSums(mod_df[,c('urban_prop',
                                        'Prop.NonArable.2010',
                                        'prop_to_predict_from_crop',
                                        'prop_to_predict_from_pasture',
                                        'pred_change_cell_prop_pasture_increase')],
                              na.rm = TRUE)
        
        # Adjust predictions if so
        # If area does exceed 1, then adjust so that predicted pasture change does not cause cell value to exceed 1
        mod_df[mod_df$tmp > 1, "pred_change_cell_prop_pasture_increase"] <-
          1 - rowSums(mod_df[mod_df$tmp > 1,c('urban_prop',
                                'Prop.NonArable.2010',
                                'prop_to_predict_from_crop',
                                'prop_to_predict_from_pasture')],
                      na.rm = TRUE)
        
        # Repeat so that the area can't be negative
        mod_df$tmp <- rowSums(mod_df[,c('pred_change_cell_prop_pasture_increase',
                                        'prop_to_predict_from_pasture')],
                              na.rm = TRUE)
        
        mod_df[mod_df$tmp < 0 & !is.na(mod_df$tmp), "pred_change_cell_prop_pasture_increase"] <-
          0 - mod_df[mod_df$tmp < 0 & !is.na(mod_df$tmp), "prop_to_predict_from_pasture"]
        
        mod_df <- dplyr::select(mod_df, -tmp)
        
        # Pick cells ----
        cells_expansion <- pasture_conv_fun(target = pasture_target, 
                                            data = mod_df,
                                            cons_planning_prop_to_reduce = cons_reduce_true)
        # I'm (naively?) assuming that pasture targets will be met. If they aren't then they aren't, and that's
        #     tough
        # Convert chosen cells ----
        # Expansion
        mod_df[mod_df$row_index %in% cells_expansion$index, "prop_to_predict_from_pasture"] <- 
          mod_df[mod_df$row_index %in% cells_expansion$index, "prop_to_predict_from_pasture"] + 
          mod_df[mod_df$row_index %in% cells_expansion$index, "pred_change_cell_prop_pasture_increase"]
        
        # Update the amount changed in the previous time step    
        mod_df[mod_df$row_index %in% cells_expansion$index, "prop_to_predict_from_cell_delta_pasture"] <- 
          mod_df[mod_df$row_index %in% cells_expansion$index, "pred_change_cell_prop_pasture_increase"]
        
        # Adjust targets -----
        pasture_target$converted <- cells_expansion$area_converted[1]
        
        # cells <- rbind(cells_expansion)
        # Take the area converted off the target
        # In an if loop because pasture_target$converted is undefined if pasture_target$target == 0
        if(pasture_target$target != 0) {
          pasture_target$target <- pasture_target$target - pasture_target$converted
        }
        
        # Subset the targets to only those that remain unmet
        met <- pasture_target[pasture_target$abs.target <= 0,]

        # Adding in logical to prevent script from looping through while loop if there are no cells left to expand
        if(!exists('cells')) {
          counter <- 16
        } else  if(nrow(cells) == 0) {
          counter <- 16
        } else {if(abs(pasture_target$target) < 1) {
          counter <- 17
        } else {
          counter <- counter + 1
        }
       
        }
        # Below bracket marks end of while loop for pasture conversion
      }
      
      # Need to reset the pasture target for each of mini-loops within each time interval
      mod_df$pasture_loss_tmp <- 0
      
      # Small logical that creates target_decade$converted if this column does not already exist
      if(!("converted" %in% colnames(target_decade))) {
        target_decade$converted <- 0
      }
      
      # Also need to adjust the absolute target into make sure that the messages reflect
      #     whether the target is met
      target_decade$abs.target <- target_decade$abs.target - abs(target_decade$converted)
      
      # Below bracket marks end for loop that crop and pasture conversion loops are embedded in
      # The for loop divides up the time interval target into n parts of equal size (assuming each part's target is met in full)
    }
    
    # Print the countries where targets are met
    if(target_decade$abs.target <= 0) {
      messages.list[[dec]][["Messages.Crop"]] <- c("Target Met", NA)
    }
    # Print countries where targets are not met 
    if(target_decade$abs.target > 0) {
      messages.list[[dec]][["Messages.Crop"]] <- c("Target Not Met - Land Remaining", target_decade$abs.target)
    }
    
    # Print the countries where targets are met for pasture
    if(abs(pasture_target$target) < 1) {
      messages.list[[dec]][["Messages.Pasture"]] <- "Pasture target Met"
    }
    # Print countries where targets are not met for pasture
    if(abs(pasture_target$target) >= 1) {
      messages.list[[dec]][["Messages.Pasture"]] <- "Pasture Target Not Met - Not Enough Available Land"
    }
    
    # Get the final area of cropland and pastureland in each cell at the end of the 
    #     time period
    # If statement. These columns are created above if target == 0
    if(!(paste("prop_crop", dec, sep = "_") %in% colnames(mod_df))) {
      mod_df[[paste("prop_crop", dec, sep = "_")]] <- mod_df$prop_to_predict_from_crop
      mod_df[[paste("prop_pasture", dec, sep = "_")]] <- mod_df$prop_to_predict_from_pasture
    }
    # The bracket below marks the end of the for() loop going through each time period
    
    
    # Dropping columns in mod_df to clear space
    # Getting index of columns
    col.identifiers <- c(grep(paste("prop_change_",
                                  dec,
                                  sep = ""),
                            names(mod_df)),
                         grep(paste("target_",
                                    dec,
                                    sep = ""),
                              names(mod_df)))

    mod_df <- mod_df[,-col.identifiers]

  }
  
  # If we only want to get the meaningful cells:
  mod_df <- mod_df[!is.na(mod_df$ISO3),]
  out.list <- list()
  out.list[["LandUse"]] <- mod_df
  out.list[["messages"]] <- messages.list
  return(out.list)
}
# out.list is a dataframe with all the messages created
# Above function returns a list of lists. One item for each model run, and each model run has warning messages and data frames


# Function for making post-conversion rasters ----
mean.sd.fun <- function(column.names) {
  
  # Puting row index into first column in matrix
  plot.table <- data.table(Row.Index = test.list[[1]]$LandUse$row_index)
  # Dropping any random columns that may have ended up in the data frame
  column.names <- column.names[grep('_20', column.names)]
  
  for (i in 1:length(column.names))  {
    # Creating empty matrix
    # Used to hold crop or pasture extent values
    matrix.holder <- matrix(NA, nrow = nrow(plot.table), ncol = length(test.list))
    # Updating values in matrix.holder
    for(rep in 1:length(test.list)){
      matrix.holder[,rep] <- test.list[[rep]]$LandUse[[column.names[i]]]
    }
    
    # And taking row means and row Sds of matrix.holder
    plot.table[[paste(column.names[i],
                      "mean",
                      sep = "_")]] <-  rowMeans(matrix.holder,
                                                na.rm = TRUE)
    
    # And putting sds into a column
     plot.table[[paste(column.names[i],
                      "sd",
                      sep = "_")]] <- rowSds(matrix.holder,
                                             na.rm = TRUE)
  }
  # Alright, now have a matrix that contains row sds and row means
  # Next need to join this matrix with a matrix that contains row index values for the entire country
  # This includes water area and ocean area
  row.index.table <- data.table(Row.Index = getValues(row.index.raster))

  # Changing row names. Don't know why I need to do this
  names(row.index.table)[1] <- "Row.Index"
  names(plot.table)[grep('row.*index', names(plot.table),ignore.case = TRUE)] <- "Row.Index"
  # And now joining two matrices together
  plot.table <- left_join(row.index.table,
                          plot.table,
			  by = setNames(names(plot.table)[1], names(row.index.table)[1]))  
  # And ordering by row.index
  plot.table <- setorder(plot.table, 
			  Row.Index)
  raster.list <- list()
  # And now making a list of rasters
  names.loop <- 
    names(plot.table)[grep('_20', names(plot.table))]
  
  for(i in names.loop) {
    raster.list[[i]] <- raster(matrix(plot.table[[i]],
                                      nrow = row.index.raster@nrows,
                                      byrow = TRUE))
    crs(raster.list[[i]]) <- crs(row.index.raster)
    extent(raster.list[[i]]) <- extent(row.index.raster)
  }
  return(raster.list)
}

# Bootstrapping function ----
bootstrap.sd.function <- function(column.names) {
  # First looking at sd across all iterations
  # Doing this to limit sample to only cells that experience a change in crop extent in any iteration
  # Don't want to artifically deflate the SD that we're sampling
  # Puting row index into first column in matrix
  plot.table <- data.table(Row.Index = test.list[[1]]$LandUse$row_index)
  
  # Looping through iterations to make a list of tables
  iteration.list <- list()
  for(n.its in 2:length(test.list)) {
    # Puting row index into first column in matrix
    plot.table <- data.table(Row.Index = test.list[[1]]$LandUse$row_index)
    
    # Dropping any random columns that may have ended up in the data frame
    column.names <- column.names[grep('_20', column.names)]
    for (i in 1:length(column.names))  {
      # Creating empty matrix
      # Used to hold crop or pasture extent values
      matrix.holder <- matrix(NA, nrow = nrow(plot.table), ncol = n.its)
      # Updating values in matrix.holder
      for(rep in 1:n.its){
        matrix.holder[,rep] <- test.list[[rep]]$LandUse[,column.names[i]]
      }
      
      # And putting sds into a column
      plot.table[[paste(column.names[i],
                        "sd",
                        sep = "_")]] <- rowSds(matrix.holder,
                                               na.rm = TRUE)
    }
    # Adding plot table to list of tables
    iteration.list[[paste("Iteration_",
                          n.its,
                          sep = "")]] <- plot.table
  }
  
  
  
  sample.vector <- 1:nrow(plot.table)
  sample.size = round(length(sample.vector)/100)
  
  # Looping through bootstrap
  for(bootstrap in 1:100) {
    # Getting subset of cells to sample
    tmp.cells <- fastSampleReject(all = sample.vector,
                                  n = sample.size,
                                  w = rep(1, length(sample.vector)))
    for(z in 1:length(iteration.list)) {
      # Taking subset of cells
      bootstrap.frame <- iteration.list[[z]][tmp.cells,]
      # Adding indicators for the bootstrap and number of iterations
      bootstrap.frame$N.Iterations <- z
      bootstrap.frame$Bootstrap <- bootstrap
      if(z == 1) {
        master.frame <- bootstrap.frame
      } else {
        master.frame <- rbind(master.frame,
                              bootstrap.frame)
      }
    }
    
    if(bootstrap == 1) {
      complete.frame <- master.frame
    }
    if(bootstrap > 1) {
      complete.frame <- rbind(complete.frame,
                              master.frame)
    }
    # Below bracket is end of bootstrap loop
  }
  return(complete.frame)
}



# Fast multinomial prdict model
multinom.predict.2b <- function(dat, coefs, vars, factor_vars) {
  # Ordering coefs data frame to make sure it matches order of the predictor vars
  # There's got to be a simple way of doing this, but I don't know this off the top of my head!
  dat <- as.data.frame(dat)
  coefs.tmp <- coefs
  coefs.tmp[,2:ncol(coefs.tmp)] <- NA
  # Dropping coefficients for factor variables
  # Will be added in later
  coefs.tmp <- coefs.tmp[,1:(length(vars) + 1)]
  # And looping through to reorder coefficients
  for(i in 2:ncol(coefs.tmp)) {
    coefs.tmp[,i] <- coefs[,colnames(coefs) == vars[i-1]]
  }
  
  
  # Adding in factor coefficients
  # First getting any coefficient for a factor variable
  coefs_factor <- coefs[,grep(factor_vars, colnames(coefs))]
  
  # Changing ot matrix if needed
  if(!is.data.frame(coefs_factor)) {
    coefs_factor = data.frame(tmp = coefs_factor)
    names(coefs_factor) = gsub('.*)','',colnames(coefs)[ncol(coefs)])
  }
  
  # Changing names of coefs_factor
  colnames(coefs_factor) <- gsub(".*)","",colnames(coefs_factor))
  
  # Getting iso numeric value
  tmp.iso = (unique(dat[,factor_vars])) %>% as.data.frame(.)
  tmp.iso = tmp.iso[!is.na(tmp.iso)]
  
  # And getting coefficients for the specific factor being fed to the data frame
  coefs_factor <- coefs_factor[,which(colnames(coefs_factor) %in% tmp.iso)]
  
  # If cateogrical variable is the baseline categorical variable, then need to create a matrix with 0s
  if(length(coefs_factor) == 0) {
    coefs_factor <- matrix(c(0,0), ncol = 1)
    row.names(coefs_factor) <- row.names(coefs)
  }
  
  # creating empty list 
  # Used to store model estimates
  vectors.list <- list()
  
  # Creating names to store in the list
  cat.names <- c('aNoChange',
                 row.names(coefs))
  
  # First category has coefficients of 0
  # Creating empty vector with these values
  vectors.list[[cat.names[1]]] <- rep(exp(0), nrow(dat))
  
  # And looping through other categories to get variable estimates
  for(n.cats in 1:nrow(coefs)) {
    # Creating matrix of coefs
    # I have no idea why as.matrix(coefs) isn't working
    # Giving an error indicating that coefs.matrix is not a matrix or vector
    # So doing this instead
    coefs.matrix = as.matrix(coefs[n.cats,2:ncol(coefs.tmp)], nrow = 12)
    coefs.matrix <- matrix(coefs.matrix[1,],nrow = 12)
    
    # Doing matrix multiplication
    # And then exponentiating
    # To get vector of relative model probabilities
    vectors.list[[cat.names[n.cats + 1]]] <- exp(data.matrix(dat[,vars]) %*% coefs.matrix + 
                                                   coefs.tmp[n.cats, 1] + coefs_factor[n.cats])
  }
  
  # Creating matrix with relative probabilities
  prob.matrix <- matrix(c(vectors.list[[1]],
                          vectors.list[[2]],
                          vectors.list[[3]]),
                        ncol = (nrow(coefs) + 1),
                        byrow = FALSE)
  
  # And then getting relative probabilities
  prob.matrix <- prob.matrix / rowSums(prob.matrix)
  # And adding column names
  colnames(prob.matrix) <- cat.names
  # And returing prob.matrix
  return(as.data.frame(prob.matrix))
}




multinom.predict.2b.trial <- function(dat, coefs, vars, factor_vars) {
  # Ordering coefs data frame to make sure it matches order of the predictor vars
  # There's got to be a simple way of doing this, but I don't know this off the top of my head!
  dat <- as.data.frame(dat)
  coefs.tmp <- coefs
  coefs.tmp[,2:ncol(coefs.tmp)] <- NA
  # Dropping coefficients for factor variables
  # Will be added in later
  coefs.tmp <- coefs.tmp[,1:(length(vars) + 1)]
  # And looping through to reorder coefficients
  for(i in 2:ncol(coefs.tmp)) {
    coefs.tmp[,i] <- coefs[,colnames(coefs) == vars[i-1]]
  }
  
  
  # Adding in factor coefficients
  # First getting any coefficient for a factor variable
  coefs_factor <- coefs[,grep(factor_vars, colnames(coefs))]
  
  # Changing ot matrix if needed
  if(!is.data.frame(coefs_factor)) {
    coefs_factor = data.frame(tmp = coefs_factor)
    names(coefs_factor) = gsub('.*)','',colnames(coefs)[ncol(coefs)])
  }
  
  # Changing names of coefs_factor
  colnames(coefs_factor) <- gsub(".*)","",colnames(coefs_factor))
  
  # Getting iso numeric value
  tmp.iso = (unique(dat[,factor_vars])) %>% as.data.frame(.)
  tmp.iso = tmp.iso[!is.na(tmp.iso)]
  
  # And getting coefficients for the specific factor being fed to the data frame
  coefs_factor <- coefs_factor[,which(colnames(coefs_factor) %in% tmp.iso)]
  
  # If cateogrical variable is the baseline categorical variable, then need to create a matrix with 0s
  if(length(coefs_factor) == 0) {
    coefs_factor <- matrix(c(0,0), ncol = 1)
    row.names(coefs_factor) <- row.names(coefs)
  }
  
  # creating empty list 
  # Used to store model estimates
  vectors.list <- list()
  
  # Creating names to store in the list
  cat.names <- c('aNoChange',
                 row.names(coefs))
  
  # First category has coefficients of 0
  # Creating empty vector with these values
  vectors.list[[cat.names[1]]] <- rep(exp(0), nrow(dat))
  
  # And looping through other categories to get variable estimates
  for(n.cats in 1:nrow(coefs)) {
    # Creating matrix of coefs
    # I have no idea why as.matrix(coefs) isn't working
    # Giving an error indicating that coefs.matrix is not a matrix or vector
    # So doing this instead
    coefs.matrix = as.matrix(coefs[n.cats,2:ncol(coefs.tmp)], nrow = length(vars))
    coefs.matrix <- matrix(coefs.matrix[1,],nrow = length(vars))
    
    # Doing matrix multiplication
    # And then exponentiating
    # To get vector of relative model probabilities
    vectors.list[[cat.names[n.cats + 1]]] <- exp(data.matrix(dat[,vars]) %*% coefs.matrix + 
                                                   coefs.tmp[n.cats, 1] + coefs_factor[n.cats])
  }
  
  # Creating matrix with relative probabilities
  prob.matrix <- matrix(c(vectors.list[[1]],
                          vectors.list[[2]],
                          vectors.list[[3]]),
                        ncol = (nrow(coefs) + 1),
                        byrow = FALSE)
  
  # And then getting relative probabilities
  prob.matrix <- prob.matrix / rowSums(prob.matrix)
  # And adding column names
  colnames(prob.matrix) <- cat.names
  # And returing prob.matrix
  return(as.data.frame(prob.matrix))
}





# Second function part two
# Does not require that variable inputs are in the same order as the model inputs
# Can deal with any number of inputs
# Assumes coefs[1] is the coefficient for the intercept
# Works with factors, but only if data frames that contain a single factor are fed into the function
glm.predict.2b <- function(dat, coefs, vars, factor_vars) {
  # Ordering coefs data frame to make sure it matches order of the predictor vars
  # There's got to be a simple way of doing this, but I don't know this off the top of my head!
  coefs.tmp <- rep(NA, length(vars) + 1)
  coefs.tmp[1] <- coefs$Estimate[coefs$variable == '(Intercept)']
  # And looping through to reorder coefficients
  for(i in 1:length(vars)) {
    coefs.tmp[i+1] <- coefs$Estimate[coefs$variable == vars[i]]
  }
  
  # Adding in factor coefficients
  # First getting any coefficient for a factor variable
  coefs_factor <- coefs$Estimate[grep(factor_vars, coefs$variable)]
  # Second getting specific coefficient for a factor variable
  # Doing this in two steps
  # First dropping name of the factor variable from the column name
  names(coefs_factor) <- coefs$variable[grep('factor',coefs$variable)] %>%gsub('.*ordered = FALSE)', "", .)
  # And second only getting the estimate for the specific factor variable
  coefs_factor <- coefs_factor[names(coefs_factor) %in% unique(dat[,factor_vars])]
  
  # And adding in factor variable to the coefficient matrix
  if(length(coefs_factor) == 0) {
    coefs.tmp <- c(coefs.tmp,
                   0)
  }
  if(length(coefs_factor) >= 1) {
    coefs.tmp <- c(coefs.tmp,
                   coefs_factor)
  }
  
  # Conerting levels to 1 for everything
  # Needed to make this work
  dat[,factor_vars] <- 1
  
  
  # And doing maths!
  return(exp(as.matrix(dat[,c(vars, factor_vars)]) %*% matrix(coefs.tmp[2:length(coefs.tmp)]) +
               coefs.tmp[1]))
}



###
# Function to identify areas of either:
# (1) Areas of agricultural contraction (e.g. no ag land in these areas by 2050, can be urban areas)
# (2) Areas of no agricultural expansion (e.g. no additional ag land in these areas, can be urban expansion)

# These input layers need to be: 0s for outside of areas, >= 1 for inside no expansion / contraction areas
# NAs are assumed to be non-land
cons_area_function <-
  function(contraction_areas,
           no_expansion_areas) {
    
    # Getting single layer for areas of ag contraction
    if(length(contraction_areas) >= 1) {
      contraction_areas <- rast(contraction_areas)
      contraction_areas <- app(contraction_areas, fun = 'sum', na.rm = TRUE)
    } else {
      contraction_areas <- 'none'
    }
    
    # Getting single layer for areas of no ag expansion
    if(length(contraction_areas) >= 1) {
      no_expansion_areas <- rast(no_expansion_areas)
      no_expansion_areas <- app(no_expansion_areas, fun = 'sum', na.rm = TRUE)
    } else {
      no_expansion_areas <- 'none'
    }
    
    # # Updating no expansion areas, if they overlap with contraction areas
    # # E.g. want to remove these, to avoid odd things happening
    # # This converts areas in the no_expansion_areas raster to 0, if areas in the contraction_areas raster are >= 1
    # if(!(is.character(contraction_areas)) &
    #    !(is.character(no_expansion_areas))) {
    #   no_expansion_areas[contraction_areas >= 1] <- 0
    # }
    # 
    # # For contraction areas, getting the amount of displaced cropland and pastureland
    # if(!is.character(contraction_areas)) {
    #   # Extracting values
    #   displaced_ag_areas <-
    #     data.frame(country_id = values(tmp_country),
    #                contraction_areas = values(contraction_areas),
    #                crop_area = values(tmp_crop),
    #                past_area = values(tmp_past))
    #   # Changing names of the data frame
    #   names(displaced_ag_areas) <-
    #     c('country_id','contraction_areas','crop_area','past_area')
    #   # Getting summary by country
    #   # Where total area for contraction is number of cells
    #   displaced_ag_areas_sum <-
    #     displaced_ag_areas %>%
    #     filter(!is.na(country_id)) %>%
    #     filter(contraction_areas >= 1) %>%
    #     dplyr::group_by(country_id) %>%
    #     dplyr::summarise(crop_area = sum(crop_area, na.rm = TRUE),
    #                      past_area = sum(past_area, na.rm = TRUE)) %>%
    #     mutate(crop_area = ifelse(is.na(crop_area),0,crop_area),
    #            past_area = ifelse(is.na(past_area),0,crop_area))
    # } # End if statement calculating how much ag land per country needs to be removed
    # 
    # And returning info
    out_list <-
      list(raster(contraction_areas),
           raster(no_expansion_areas))
    
    names(out_list) <- 
      c('contraction_areas_raster',
        'no_expansion_areas_raster')
    
    rm(contraction_areas,
       no_expansion_areas)
    
    return(out_list)
  }
