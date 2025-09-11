#####
# Function to estimate loss of esh from different stresses
#####

esh.loss.function.old <-
  function(current,
           stresses,
           species.df,
           tmp.crop.pref,
           tmp.past.pref,
           tmp.urban.pref) {
    
    ###
    # Creating temporary vectors
    # Converting 0s to NAs to make script work
    store.crop.2010 <- species.df$crop.2010
    store.crop.values <- species.df$crop.values
    
    store.past.2010 <- species.df$past.2010
    store.past.values <- species.df$pasture.values
    
    species.df$crop.2010[is.na(species.df$crop.2010)] <- 0
    species.df$crop.values[is.na(species.df$crop.values)] <- 0
    
    species.df$past.2010[is.na(species.df$past.2010)] <- 0
    species.df$pasture.values[is.na(species.df$pasture.values)] <- 0
    
    
    if('climate' %in% stresses) {
      ### Creating temporary vectors for esh 45 and 85
      tmp.esh.sdm <- species.df$esh.values
      tmp.esh.sdm[!is.na(species.df$esh.values) & is.na(species.df$esh.values.sdm)] <- NA
      
      tmp.esh.45 <- species.df$esh.values
      tmp.esh.45[!is.na(species.df$esh.values) & is.na(species.df$esh.values.45)] <- NA
      
      # Second for RCP 85
      tmp.esh.85 <- species.df$esh.values
      tmp.esh.85[!is.na(species.df$esh.values) & is.na(species.df$esh.values.85)] <- NA
    }
    
    
    
    ###
    # Getting current esh loss
    if(current == 'yes') {
      ###
      # Adding esh current based on current ag extent and intensity
      
      ### 
      # Calculating ESH loss based on actual ESH raster from Team Europe
      # Based on ESH raster
      # This is for extent of ag
      species.df$esh.current.extent.esh <-
        species.df$esh.values - 
        species.df$crop.2010 * species.df$crop.intens.2010 -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.2010.values * (tmp.urban.pref)
      
      
      ###
      # And repeating for SDMs
      # If statement here to save time
      if('climate' %in% stresses) {
        species.df$esh.current.extent.sdm <-
          species.df$esh.values.sdm - 
          species.df$crop.2010 * species.df$crop.intens.2010 -
          species.df$past.2010 * species.df$past.intens.2010 -
          species.df$urban.2010.values * (tmp.urban.pref)
        
        
        
        species.df$esh.current.extent.esh.sdm <-
          tmp.esh.sdm -
          species.df$crop.2010 * species.df$crop.intens.2010 -
          species.df$past.2010 * species.df$past.intens.2010 -
          species.df$urban.2010.values * (tmp.urban.pref)
      }
    } # End of estimate of current esh extent
    
    ###
    # Loop for ag expansion
    # This looks at ESH from ag expansion relative to current ESH
    # e.g., extent and intensity of current ag production is accounted for
    # But only change in extent of future ag production is incorporate into this
    # Unless extent decreases, in which intensity penalty in those cells is also ecreased
    if('expansion' %in% stresses) {
      # Calculating loss from hab expansion
      species.df$esh.future.extent.esh.exp <-
        species.df$esh.values - 
        species.df$crop.values * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.2010 -
        species.df$urban.2010.values * (tmp.urban.pref)
      
      ###
      # Repeating for SDM model under current climate
      if('climate' %in% stresses) {
        # Account for habitat loss from changes in crop extent
        species.df$esh.future.extent.sdm.exp <-
          species.df$esh.values.sdm - 
          species.df$crop.values * species.df$crop.intens.2010 -
          species.df$pasture.values * species.df$past.intens.2010 -
          species.df$urban.2010.values * (tmp.urban.pref)
        
        species.df$esh.future.extent.esh.sdm.exp <-
          tmp.esh.sdm - 
          species.df$crop.values * species.df$crop.intens.2010 -
          species.df$pasture.values * species.df$past.intens.2010 -
          species.df$urban.2010.values * (tmp.urban.pref)
        
      } 
    } # End of function for expansion
    
    ###
    # Loop for ag intensification
    # This looks at ESH loss from ag intensification relative to current ESH
    # e.g. extent and intensity of current ag production is accounted for
    # But only change in intensification on existing cells is accounted for
    # Unless ag extent decreases in cells currently occupied by ag
    if('intensification' %in% stresses) {
      # Isolating the impact of increased intensification on currently existing ag landscapes
      # This is for intensity of ag
      
      ###
      # First calculating habitat loss from on current ag extent
      species.df$esh.future.extent.esh.int <-
        species.df$esh.values - 
        species.df$crop.2010* species.df$crop.intens.tmp -
        species.df$past.2010 * species.df$past.intens.tmp -
        species.df$urban.2010.values * (tmp.urban.pref)
      
      ###
      # Repeating for SDM model under current climate
      if('climate' %in% stresses) {
        # Isolating the impact of increased intensification on currently existing ag landscapes
        # This is for intensity of ag
        
        ###
        # First calculating esh loss from current cropland extent
        species.df$esh.future.extent.sdm.int <-
          species.df$esh.values.sdm - 
          species.df$crop.2010* species.df$crop.intens.tmp -
          species.df$past.2010 * species.df$past.intens.tmp -
          species.df$urban.2010.values * (tmp.urban.pref)
        
        species.df$esh.future.extent.esh.sdm.int <-
          tmp.esh.sdm - 
          species.df$crop.2010* species.df$crop.intens.tmp -
          species.df$past.2010 * species.df$past.intens.tmp -
          species.df$urban.2010.values * (tmp.urban.pref)
      }
    }
    
    ###
    # Loop for climate change
    # Shirnks/expands habitat based on outputs from SDM model
    if('climate' %in% stresses) {
      ###
      # For current esh models, limit ag extent with current extent
      # First for RCP 45
      
      ###
      # Accounting for hab loss within historical crop extent
      species.df$esh.future.extent.esh.sdm.rcp45 <-
        tmp.esh.45 - 
        species.df$crop.2010* species.df$crop.intens.2010 -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.2010.values * (tmp.urban.pref)
      
      species.df$esh.future.extent.esh.sdm.rcp85 <-
        tmp.esh.85 -
        species.df$crop.2010* species.df$crop.intens.2010 -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.2010.values * (tmp.urban.pref)
      
      
      
      ###
      # And repeating for SDM models
      # Don't need to clip the rasters beause we've already done this
      # And have already calculated esh loss from ag extent in 2010 under the current SDM
      species.df$esh.future.extent.sdm.rcp45 <-
        species.df$esh.values.45 - 
        species.df$crop.2010* species.df$crop.intens.2010 -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.2010.values * (tmp.urban.pref)
      
      species.df$esh.future.extent.sdm.rcp85 <-
        species.df$esh.values.85 - 
        species.df$crop.2010* species.df$crop.intens.2010 -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.2010.values * (tmp.urban.pref)
    }
    
    ###
    # Loop for climate change and expansion
    if('climate' %in% stresses &
       'expansion' %in% stresses) {
      ###
      # We've already accounted for climate in the above if statement
      # Now accounting for expansion on top of climate
      
      ###
      # Now accounting for expansion of ag
      species.df$esh.future.extent.esh.sdm.rcp45.exp <-
        tmp.esh.45 -
        species.df$crop.values * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.2010 -
        species.df$urban.2010.values * (tmp.urban.pref)
      
      species.df$esh.future.extent.esh.sdm.rcp85.exp <-
        tmp.esh.85 - 
        species.df$crop.values * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.2010 -
        species.df$urban.2010.values * (tmp.urban.pref)
      
      ###
      # Repeating for the SDM rasters
      # Second for the SDM rasters
      species.df$esh.future.extent.sdm.rcp45.exp <-
        species.df$esh.values.45 - 
        species.df$crop.values * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.2010 -
        species.df$urban.2010.values * (tmp.urban.pref)
      
      species.df$esh.future.extent.sdm.rcp85.exp <-
        species.df$esh.values.85 - 
        species.df$crop.values * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.2010 -
        species.df$urban.2010.values * (tmp.urban.pref)
    }
    
    ###
    # Loop for climate change and intensification
    if('climate' %in% stresses &
       'intensification' %in% stresses) {
      ###
      # Updating ESH rasters for climate change
      
      ###
      # First for esh rasters
      # First step is to account for hab loss from ag extent in 2010
      species.df$esh.future.extent.esh.sdm.rcp45.int <-
        tmp.esh.45 - 
        species.df$crop.2010 * species.df$crop.intens.tmp -
        species.df$past.2010 * species.df$past.intens.tmp -
        species.df$urban.2010.values * (tmp.urban.pref)
      
      species.df$esh.future.extent.esh.sdm.rcp85.int <-
        tmp.esh.85 - 
        species.df$crop.2010 * species.df$crop.intens.tmp -
        species.df$past.2010 * species.df$past.intens.tmp -
        species.df$urban.2010.values * (tmp.urban.pref)
      
      ###
      # Second for sdm rasters
      # Don't need to adjust for climate because this is already done in the SDM models
      
      ###
      # First step is to account for hab loss from ag extent in 2010
      species.df$esh.future.extent.sdm.rcp45.int <-
        species.df$esh.values.45 - 
        species.df$crop.2010 * species.df$crop.intens.tmp -
        species.df$past.2010 * species.df$past.intens.tmp -
        species.df$urban.2010.values * (tmp.urban.pref)
      
      species.df$esh.future.extent.sdm.rcp85.int <-
        species.df$esh.values.85 - 
        species.df$crop.2010 * species.df$crop.intens.tmp -
        species.df$past.2010 * species.df$past.intens.tmp -
        species.df$urban.2010.values * (tmp.urban.pref)
    }
    
    ###
    # Loop for expansion and intensification
    if('expansion' %in% stresses &
       'intensification' %in% stresses) {
      
      ###
      # ESH models first
      
      ###
      # First account for hab loss
      species.df$esh.future.extent.esh.exp.int <-
        species.df$esh.values - 
        species.df$crop.values * species.df$crop.intens.tmp -
        species.df$pasture.values * species.df$past.intens.tmp -
        species.df$urban.2010.values * (tmp.urban.pref)
      
      ###
      # Repeating for SDM model under current climate
      if(c('climate' %in% stresses)) {
        species.df$esh.future.extent.esh.sdm.exp.int <-
          tmp.esh.sdm - 
          species.df$crop.values * species.df$crop.intens.tmp -
          species.df$pasture.values * species.df$past.intens.tmp -
          species.df$urban.2010.values * (tmp.urban.pref)
        
        species.df$esh.future.extent.sdm.exp.int <-
          species.df$esh.values.sdm - 
          species.df$crop.values * species.df$crop.intens.tmp -
          species.df$pasture.values * species.df$past.intens.tmp -
          species.df$urban.2010.values * (tmp.urban.pref)
      }
    }
    
    ###
    # Loop for all three
    if('expansion' %in% stresses &
       'intensification' %in% stresses &
       'climate' %in% stresses) {
      
      ###
      # First ESH models
      
      ###
      # Accounting for climate
      # Updating ESH rasters for climate change
      
      ###
      # Accounting for extensification
      species.df$esh.future.extent.esh.sdm.rcp45.exp.int <-
        tmp.esh.45 - 
        species.df$crop.values * species.df$crop.intens.tmp -
        species.df$pasture.values * species.df$past.intens.tmp -
        species.df$urban.2010.values * (tmp.urban.pref)
      
      species.df$esh.future.extent.esh.sdm.rcp85.exp.int <-
        tmp.esh.85 - 
        species.df$crop.values * species.df$crop.intens.tmp -
        species.df$pasture.values * species.df$past.intens.tmp -
        species.df$urban.2010.values * (tmp.urban.pref)
      
      
      ###
      # And repeating for the SDMs
      
      ###
      # Accounting for climate
      # This is already done in the SDM models
      
      ###
      # Accounting for extensification
      species.df$esh.future.extent.sdm.rcp45.exp.int <-
        species.df$esh.values.45 - 
        species.df$crop.values * species.df$crop.intens.tmp -
        species.df$pasture.values * species.df$past.intens.tmp -
        species.df$urban.2010.values * (tmp.urban.pref)
      
      species.df$esh.future.extent.sdm.rcp85.exp.int <-
        species.df$esh.values.85 - 
        species.df$crop.values * species.df$crop.intens.tmp -
        species.df$pasture.values * species.df$past.intens.tmp -
        species.df$urban.2010.values * (tmp.urban.pref)
    }
    
    ###
    # Adding urban extent in
    if('urbanization' %in% stresses) {
      # Calculating loss from hab expansion
      species.df$esh.future.extent.esh.urb <-
        species.df$esh.values - 
        species.df$crop.2010 * species.df$crop.intens.2010 -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.values * (tmp.urban.pref)
      
      ###
      # Repeating for SDM model under current climate
      if('climate' %in% stresses) {
        # Account for habitat loss from changes in crop extent
        species.df$esh.future.extent.sdm.urb <-
          species.df$esh.values.sdm - 
          species.df$crop.2010 * species.df$crop.intens.2010 -
          species.df$past.2010 * species.df$past.intens.2010 -
          species.df$urban.values * (tmp.urban.pref)
        
        species.df$esh.future.extent.esh.sdm.urb <-
          tmp.esh.sdm - 
          species.df$crop.2010 * species.df$crop.intens.2010 -
          species.df$past.2010 * species.df$past.intens.2010 -
          species.df$urban.values * (tmp.urban.pref)
      } 
    } # End of if statement for urbanization
    
    # Unless extent decreases, in which intensity penalty in those cells is also ecreased
    if('expansion' %in% stresses &
       'urbanization' %in% stresses) {
      # Calculating loss from hab expansion
      species.df$esh.future.extent.esh.exp.urb <-
        species.df$esh.values - 
        species.df$crop.values * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.2010 -
        species.df$urban.values * (tmp.urban.pref)
      
      ###
      # Repeating for SDM model under current climate
      if('climate' %in% stresses) {
        # Account for habitat loss from changes in crop extent
        species.df$esh.future.extent.sdm.exp.urb <-
          species.df$esh.values.sdm - 
          species.df$crop.values * species.df$crop.intens.2010 -
          species.df$pasture.values * species.df$past.intens.2010 -
          species.df$urban.values * (tmp.urban.pref)
        
        species.df$esh.future.extent.esh.sdm.exp.urb <-
          tmp.esh.sdm - 
          species.df$crop.values * species.df$crop.intens.2010 -
          species.df$pasture.values * species.df$past.intens.2010 -
          species.df$urban.values * (tmp.urban.pref)
        
      } 
    } # End of function for expansion
    
    ###
    # Loop for ag intensification
    # This looks at ESH loss from ag intensification relative to current ESH
    # e.g. extent and intensity of current ag production is accounted for
    # But only change in intensification on existing cells is accounted for
    # Unless ag extent decreases in cells currently occupied by ag
    if('intensification' %in% stresses &
       'urbanization' %in% stresses) {
      # Isolating the impact of increased intensification on currently existing ag landscapes
      # This is for intensity of ag
      
      ###
      # First calculating habitat loss from on current ag extent
      species.df$esh.future.extent.esh.int.urb <-
        species.df$esh.values - 
        species.df$crop.2010* species.df$crop.intens.tmp -
        species.df$past.2010 * species.df$past.intens.tmp -
        species.df$urban.values * (tmp.urban.pref)
      
      ###
      # Repeating for SDM model under current climate
      if('climate' %in% stresses) {
        # Isolating the impact of increased intensification on currently existing ag landscapes
        # This is for intensity of ag
        
        ###
        # First calculating esh loss from current cropland extent
        species.df$esh.future.extent.sdm.int.urb <-
          species.df$esh.values.sdm - 
          species.df$crop.2010* species.df$crop.intens.tmp -
          species.df$past.2010 * species.df$past.intens.tmp -
          species.df$urban.values * (tmp.urban.pref)
        
        species.df$esh.future.extent.esh.sdm.int.urb <-
          tmp.esh.sdm - 
          species.df$crop.2010* species.df$crop.intens.tmp -
          species.df$past.2010 * species.df$past.intens.tmp -
          species.df$urban.values * (tmp.urban.pref)
      }
    }
    
    ###
    # Loop for climate change
    # Shirnks/expands habitat based on outputs from SDM model
    if('climate' %in% stresses &
       'urbanization' %in% stresses) {
      ###
      # For current esh models, limit ag extent with current extent
      # First for RCP 45
      
      ###
      # Accounting for hab loss within historical crop extent
      species.df$esh.future.extent.esh.sdm.rcp45.urb <-
        tmp.esh.45 - 
        species.df$crop.2010* species.df$crop.intens.2010 -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.values * (tmp.urban.pref)
      
      species.df$esh.future.extent.esh.sdm.rcp85.urb <-
        tmp.esh.85 -
        species.df$crop.2010* species.df$crop.intens.2010 -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.values * (tmp.urban.pref)
      
      
      
      ###
      # And repeating for SDM models
      # Don't need to clip the rasters beause we've already done this
      # And have already calculated esh loss from ag extent in 2010 under the current SDM
      species.df$esh.future.extent.sdm.rcp45.urb <-
        species.df$esh.values.45 - 
        species.df$crop.2010* species.df$crop.intens.2010 -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.values * (tmp.urban.pref)
      
      species.df$esh.future.extent.sdm.rcp85.urb <-
        species.df$esh.values.85 - 
        species.df$crop.2010* species.df$crop.intens.2010 -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.values * (tmp.urban.pref)
    }
    
    ###
    # Loop for climate change and expansion
    if('climate' %in% stresses &
       'expansion' %in% stresses &
       'urbanization' %in% stresses) {
      ###
      # We've already accounted for climate in the above if statement
      # Now accounting for expansion on top of climate
      
      ###
      # Now accounting for expansion of ag
      species.df$esh.future.extent.esh.sdm.rcp45.exp.urb <-
        tmp.esh.45 -
        species.df$crop.values * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.2010 -
        species.df$urban.values * (tmp.urban.pref)
      
      species.df$esh.future.extent.esh.sdm.rcp85.exp.urb <-
        tmp.esh.85 - 
        species.df$crop.values * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.2010 -
        species.df$urban.values * (tmp.urban.pref)
      
      ###
      # Repeating for the SDM rasters
      # Second for the SDM rasters
      species.df$esh.future.extent.sdm.rcp45.exp.urb <-
        species.df$esh.values.45 - 
        species.df$crop.values * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.2010 -
        species.df$urban.values * (tmp.urban.pref)
      
      species.df$esh.future.extent.sdm.rcp85.exp.urb <-
        species.df$esh.values.85 - 
        species.df$crop.values * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.2010 -
        species.df$urban.values * (tmp.urban.pref)
    }
    
    ###
    # Loop for climate change and intensification
    if('climate' %in% stresses &
       'intensification' %in% stresses &
       'urbanization' %in% stresses) {
      ###
      # Updating ESH rasters for climate change
      
      ###
      # First for esh rasters
      # First step is to account for hab loss from ag extent in 2010
      species.df$esh.future.extent.esh.sdm.rcp45.int.urb <-
        tmp.esh.45 - 
        species.df$crop.2010 * species.df$crop.intens.tmp -
        species.df$past.2010 * species.df$past.intens.tmp -
        species.df$urban.values * (tmp.urban.pref)
      
      species.df$esh.future.extent.esh.sdm.rcp85.int.urb <-
        tmp.esh.85 - 
        species.df$crop.2010 * species.df$crop.intens.tmp -
        species.df$past.2010 * species.df$past.intens.tmp -
        species.df$urban.values * (tmp.urban.pref)
      
      ###
      # Second for sdm rasters
      # Don't need to adjust for climate because this is already done in the SDM models
      
      ###
      # First step is to account for hab loss from ag extent in 2010
      species.df$esh.future.extent.sdm.rcp45.int.urb <-
        species.df$esh.values.45 - 
        species.df$crop.2010 * species.df$crop.intens.tmp -
        species.df$past.2010 * species.df$past.intens.tmp -
        species.df$urban.values * (tmp.urban.pref)
      
      species.df$esh.future.extent.sdm.rcp85.int.urb <-
        species.df$esh.values.85 - 
        species.df$crop.2010 * species.df$crop.intens.tmp -
        species.df$past.2010 * species.df$past.intens.tmp -
        species.df$urban.values * (tmp.urban.pref)
    }
    
    ###
    # Loop for expansion and intensification
    if('expansion' %in% stresses &
       'intensification' %in% stresses &
       'urbanization' %in% stresses) {
      
      ###
      # ESH models first
      
      ###
      # First account for hab loss
      species.df$esh.future.extent.esh.exp.int.urb <-
        species.df$esh.values - 
        species.df$crop.values * species.df$crop.intens.tmp -
        species.df$pasture.values * species.df$past.intens.tmp -
        species.df$urban.values * (tmp.urban.pref)
      
      ###
      # Repeating for SDM model under current climate
      if(c('climate' %in% stresses)) {
        species.df$esh.future.extent.esh.sdm.exp.int.urb <-
          tmp.esh.sdm - 
          species.df$crop.values * species.df$crop.intens.tmp -
          species.df$pasture.values * species.df$past.intens.tmp -
          species.df$urban.values * (tmp.urban.pref)
        
        species.df$esh.future.extent.sdm.exp.int.urb <-
          species.df$esh.values.sdm - 
          species.df$crop.values * species.df$crop.intens.tmp -
          species.df$pasture.values * species.df$past.intens.tmp -
          species.df$urban.values * (tmp.urban.pref)
      }
    }
    
    ###
    # Loop for all three
    if('expansion' %in% stresses &
       'intensification' %in% stresses &
       'climate' %in% stresses &
       'urbanization' %in% stresses) {
      
      ###
      # First ESH models
      
      ###
      # Accounting for climate
      # Updating ESH rasters for climate change
      
      ###
      # Accounting for extensification
      species.df$esh.future.extent.esh.sdm.rcp45.exp.int.urb <-
        tmp.esh.45 - 
        species.df$crop.values * species.df$crop.intens.tmp -
        species.df$pasture.values * species.df$past.intens.tmp -
        species.df$urban.values * (tmp.urban.pref)
      
      species.df$esh.future.extent.esh.sdm.rcp85.exp.int.urb <-
        tmp.esh.85 - 
        species.df$crop.values * species.df$crop.intens.tmp -
        species.df$pasture.values * species.df$past.intens.tmp -
        species.df$urban.values * (tmp.urban.pref)
      
      
      ###
      # And repeating for the SDMs
      
      ###
      # Accounting for climate
      # This is already done in the SDM models
      
      ###
      # Accounting for extensification
      species.df$esh.future.extent.sdm.rcp45.exp.int.urb <-
        species.df$esh.values.45 - 
        species.df$crop.values * species.df$crop.intens.tmp -
        species.df$pasture.values * species.df$past.intens.tmp -
        species.df$urban.values * (tmp.urban.pref)
      
      species.df$esh.future.extent.sdm.rcp85.exp.int.urb <-
        species.df$esh.values.85 - 
        species.df$crop.values * species.df$crop.intens.tmp -
        species.df$pasture.values * species.df$past.intens.tmp -
        species.df$urban.values * (tmp.urban.pref)
    }
    
    ###
    # Restoring crop and pasture values
    species.df$crop.2010 <- store.crop.2010
    species.df$crop.values <- store.crop.values
    
    species.df$past.2010 <- store.past.2010
    species.df$pasture.values <- store.past.values
    
    
    # Removing objects to save space...
    # Exception for climate - don't want to throw a warning
    if('climate' %in% stresses) {
      rm(tmp.esh.45)
      rm(tmp.esh.85)  
    }
    
    rm(store.crop.2010)
    rm(store.crop.values)
    rm(store.past.2010)
    rm(store.past.values)
    
    ###
    # And returning the data frame
    return(species.df)
  }



# New AOH function
#####
# Function to estimate loss of esh from different stresses
#####

esh.loss.function <-
  function(current,
           stresses,
           species.df,
           tmp.crop.pref,
           tmp.past.pref,
           tmp.urban.pref) {

    ###
    # Creating temporary vectors
    # Converting 0s to NAs to make script work
    store.crop.2010 <- species.df$crop.2010
    store.crop.values <- species.df$crop.values

    store.past.2010 <- species.df$past.2010
    store.past.values <- species.df$pasture.values

    species.df$crop.2010[is.na(species.df$crop.2010)] <- 0
    species.df$crop.values[is.na(species.df$crop.values)] <- 0

    species.df$past.2010[is.na(species.df$past.2010)] <- 0
    species.df$pasture.values[is.na(species.df$pasture.values)] <- 0

    # Calculating 'absolute' impacts
    # Assuming only one stressor exists
    # Cannot do this for intensification
    # As this cannot exist in the absence of ag expansion

    # Ag expansion
    if('expansion' %in% stresses) {
      species.df$esh.agexp.abs <-
        species.df$esh.values -
        species.df$crop.values * species.df$crop.intens.tmp -
        species.df$pasture.values * species.df$past.intens.tmp

      species.df$esh.cropexp.abs <-
        species.df$esh.values -
        species.df$crop.values * species.df$crop.intens.tmp

      species.df$esh.pastexp.abs <-
        species.df$esh.values -
        species.df$pasture.values * species.df$past.intens.tmp
    }

    # Urb expansion
    if('urbanization' %in% stresses) {
      species.df$esh.urb.abs <-
        species.df$esh.values -
        species.df$urban.values * (tmp.urban.pref)
    }



    # Calculating 'marginal' impacts
    # Assuming only one stressor changes from 2010 baseline
    # Coutnerfactual = 2010, all stressors

    if('climate' %in% stresses) {
      ### Creating temporary vectors for esh 45 and 85
      tmp.esh.sdm <- species.df$esh.values
      tmp.esh.sdm[!is.na(species.df$esh.values) & is.na(species.df$esh.values.sdm)] <- NA

      tmp.esh.45 <- species.df$esh.values
      tmp.esh.45[!is.na(species.df$esh.values) & is.na(species.df$esh.values.45)] <- NA

      # Second for RCP 85
      tmp.esh.85 <- species.df$esh.values
      tmp.esh.85[!is.na(species.df$esh.values) & is.na(species.df$esh.values.85)] <- NA
    }



    ###
    # Getting current esh loss
    if(current == 'yes') {
      ###
      # Adding esh current based on current ag extent and intensity

      ###
      # Calculating ESH loss based on actual ESH raster from Team Europe
      # Based on ESH raster
      # This is for extent of ag
      species.df$esh.current.extent.esh <-
        species.df$esh.values -
        species.df$crop.2010 * species.df$crop.intens.2010 -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.2010.values * (tmp.urban.pref)


      ###
      # And repeating for SDMs
      # If statement here to save time
      if('climate' %in% stresses) {
        species.df$esh.current.extent.sdm <-
          species.df$esh.values.sdm -
          species.df$crop.2010 * species.df$crop.intens.2010 -
          species.df$past.2010 * species.df$past.intens.2010 -
          species.df$urban.2010.values * (tmp.urban.pref)



        species.df$esh.current.extent.esh.sdm <-
          tmp.esh.sdm -
          species.df$crop.2010 * species.df$crop.intens.2010 -
          species.df$past.2010 * species.df$past.intens.2010 -
          species.df$urban.2010.values * (tmp.urban.pref)
      }
    } # End of estimate of current esh extent

    ###
    # Loop for ag expansion
    # This looks at ESH from ag expansion relative to current ESH
    # e.g., extent and intensity of current ag production is accounted for
    # But only change in extent of future ag production is incorporate into this
    # Unless extent decreases, in which intensity penalty in those cells is also ecreased
    if('expansion' %in% stresses) {
      # Calculating loss from hab expansion

      # All ag changes
      species.df$esh.future.extent.esh.exp <-
        species.df$esh.values -
        species.df$crop.values * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.2010 -
        species.df$urban.2010.values * (tmp.urban.pref)

      # Crop only changes
      species.df$esh.future.extent.esh.cropexp <-
        species.df$esh.values -
        species.df$crop.values * species.df$crop.intens.2010 -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.2010.values * (tmp.urban.pref)

      # Pasture only changes
      species.df$esh.future.extent.esh.pastexp <-
        species.df$esh.values -
        species.df$crop.2010 * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.2010 -
        species.df$urban.2010.values * (tmp.urban.pref)

      ###
      # Repeating for SDM model under current climate
      if('climate' %in% stresses) {
        # Account for habitat loss from changes in crop extent
        species.df$esh.future.extent.sdm.exp <-
          species.df$esh.values.sdm -
          species.df$crop.values * species.df$crop.intens.2010 -
          species.df$pasture.values * species.df$past.intens.2010 -
          species.df$urban.2010.values * (tmp.urban.pref)

        species.df$esh.future.extent.esh.sdm.exp <-
          tmp.esh.sdm -
          species.df$crop.values * species.df$crop.intens.2010 -
          species.df$pasture.values * species.df$past.intens.2010 -
          species.df$urban.2010.values * (tmp.urban.pref)

      }
    } # End of function for expansion

    ###
    # Loop for ag intensification
    # This looks at ESH loss from ag intensification relative to current ESH
    # e.g. extent and intensity of current ag production is accounted for
    # But only change in intensification on existing cells is accounted for
    # Unless ag extent decreases in cells currently occupied by ag
    if('intensification' %in% stresses) {
      # Isolating the impact of increased intensification on currently existing ag landscapes
      # This is for intensity of ag

      ###
      # First calculating habitat loss from on current ag extent
      species.df$esh.future.extent.esh.int <-
        species.df$esh.values -
        species.df$crop.2010* species.df$crop.intens.tmp -
        species.df$past.2010 * species.df$past.intens.tmp -
        species.df$urban.2010.values * (tmp.urban.pref)

      ###
      # Repeating for SDM model under current climate
      if('climate' %in% stresses) {
        # Isolating the impact of increased intensification on currently existing ag landscapes
        # This is for intensity of ag

        ###
        # First calculating esh loss from current cropland extent
        species.df$esh.future.extent.sdm.int <-
          species.df$esh.values.sdm -
          species.df$crop.2010* species.df$crop.intens.tmp -
          species.df$past.2010 * species.df$past.intens.tmp -
          species.df$urban.2010.values * (tmp.urban.pref)

        species.df$esh.future.extent.esh.sdm.int <-
          tmp.esh.sdm -
          species.df$crop.2010* species.df$crop.intens.tmp -
          species.df$past.2010 * species.df$past.intens.tmp -
          species.df$urban.2010.values * (tmp.urban.pref)
      }
    }

    ###
    # Loop for climate change
    # Shirnks/expands habitat based on outputs from SDM model
    if('climate' %in% stresses) {
      ###
      # For current esh models, limit ag extent with current extent
      # First for RCP 45

      ###
      # Accounting for hab loss within historical crop extent
      species.df$esh.future.extent.esh.sdm.rcp45 <-
        tmp.esh.45 -
        species.df$crop.2010* species.df$crop.intens.2010 -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.2010.values * (tmp.urban.pref)

      species.df$esh.future.extent.esh.sdm.rcp85 <-
        tmp.esh.85 -
        species.df$crop.2010* species.df$crop.intens.2010 -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.2010.values * (tmp.urban.pref)



      ###
      # And repeating for SDM models
      # Don't need to clip the rasters beause we've already done this
      # And have already calculated esh loss from ag extent in 2010 under the current SDM
      species.df$esh.future.extent.sdm.rcp45 <-
        species.df$esh.values.45 -
        species.df$crop.2010* species.df$crop.intens.2010 -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.2010.values * (tmp.urban.pref)

      species.df$esh.future.extent.sdm.rcp85 <-
        species.df$esh.values.85 -
        species.df$crop.2010* species.df$crop.intens.2010 -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.2010.values * (tmp.urban.pref)
    }

    ###
    # Loop for climate change and expansion
    if('climate' %in% stresses &
       'expansion' %in% stresses) {
      ###
      # We've already accounted for climate in the above if statement
      # Now accounting for expansion on top of climate

      ###
      # Now accounting for expansion of ag
      species.df$esh.future.extent.esh.sdm.rcp45.exp <-
        tmp.esh.45 -
        species.df$crop.values * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.2010 -
        species.df$urban.2010.values * (tmp.urban.pref)

      species.df$esh.future.extent.esh.sdm.rcp85.exp <-
        tmp.esh.85 -
        species.df$crop.values * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.2010 -
        species.df$urban.2010.values * (tmp.urban.pref)

      ###
      # Repeating for the SDM rasters
      # Second for the SDM rasters
      species.df$esh.future.extent.sdm.rcp45.exp <-
        species.df$esh.values.45 -
        species.df$crop.values * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.2010 -
        species.df$urban.2010.values * (tmp.urban.pref)

      species.df$esh.future.extent.sdm.rcp85.exp <-
        species.df$esh.values.85 -
        species.df$crop.values * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.2010 -
        species.df$urban.2010.values * (tmp.urban.pref)
    }

    ###
    # Loop for climate change and intensification
    if('climate' %in% stresses &
       'intensification' %in% stresses) {
      ###
      # Updating ESH rasters for climate change

      ###
      # First for esh rasters
      # First step is to account for hab loss from ag extent in 2010
      species.df$esh.future.extent.esh.sdm.rcp45.int <-
        tmp.esh.45 -
        species.df$crop.2010 * species.df$crop.intens.tmp -
        species.df$past.2010 * species.df$past.intens.tmp -
        species.df$urban.2010.values * (tmp.urban.pref)

      species.df$esh.future.extent.esh.sdm.rcp85.int <-
        tmp.esh.85 -
        species.df$crop.2010 * species.df$crop.intens.tmp -
        species.df$past.2010 * species.df$past.intens.tmp -
        species.df$urban.2010.values * (tmp.urban.pref)

      ###
      # Second for sdm rasters
      # Don't need to adjust for climate because this is already done in the SDM models

      ###
      # First step is to account for hab loss from ag extent in 2010
      species.df$esh.future.extent.sdm.rcp45.int <-
        species.df$esh.values.45 -
        species.df$crop.2010 * species.df$crop.intens.tmp -
        species.df$past.2010 * species.df$past.intens.tmp -
        species.df$urban.2010.values * (tmp.urban.pref)

      species.df$esh.future.extent.sdm.rcp85.int <-
        species.df$esh.values.85 -
        species.df$crop.2010 * species.df$crop.intens.tmp -
        species.df$past.2010 * species.df$past.intens.tmp -
        species.df$urban.2010.values * (tmp.urban.pref)
    }

    ###
    # Loop for expansion and intensification
    if('expansion' %in% stresses &
       'intensification' %in% stresses) {

      ###
      # ESH models first

      ###
      # First account for hab loss

      # All ag changes
      species.df$esh.future.extent.esh.exp.int <-
        species.df$esh.values -
        species.df$crop.values * species.df$crop.intens.tmp -
        species.df$pasture.values * species.df$past.intens.tmp -
        species.df$urban.2010.values * (tmp.urban.pref)

      # Crop only changes
      species.df$esh.future.extent.esh.cropexp.int <-
        species.df$esh.values -
        species.df$crop.values * species.df$crop.intens.tmp -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.2010.values * (tmp.urban.pref)

      # Pasture only changes
      species.df$esh.future.extent.esh.pastexp.int <-
        species.df$esh.values -
        species.df$crop.2010 * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.tmp -
        species.df$urban.2010.values * (tmp.urban.pref)

      ###
      # Repeating for SDM model under current climate
      if(c('climate' %in% stresses)) {
        species.df$esh.future.extent.esh.sdm.exp.int <-
          tmp.esh.sdm -
          species.df$crop.values * species.df$crop.intens.tmp -
          species.df$pasture.values * species.df$past.intens.tmp -
          species.df$urban.2010.values * (tmp.urban.pref)

        species.df$esh.future.extent.sdm.exp.int <-
          species.df$esh.values.sdm -
          species.df$crop.values * species.df$crop.intens.tmp -
          species.df$pasture.values * species.df$past.intens.tmp -
          species.df$urban.2010.values * (tmp.urban.pref)
      }
    }

    ###
    # Loop for all three
    if('expansion' %in% stresses &
       'intensification' %in% stresses &
       'climate' %in% stresses) {

      ###
      # First ESH models

      ###
      # Accounting for climate
      # Updating ESH rasters for climate change

      ###
      # Accounting for extensification
      species.df$esh.future.extent.esh.sdm.rcp45.exp.int <-
        tmp.esh.45 -
        species.df$crop.values * species.df$crop.intens.tmp -
        species.df$pasture.values * species.df$past.intens.tmp -
        species.df$urban.2010.values * (tmp.urban.pref)

      species.df$esh.future.extent.esh.sdm.rcp85.exp.int <-
        tmp.esh.85 -
        species.df$crop.values * species.df$crop.intens.tmp -
        species.df$pasture.values * species.df$past.intens.tmp -
        species.df$urban.2010.values * (tmp.urban.pref)


      ###
      # And repeating for the SDMs

      ###
      # Accounting for climate
      # This is already done in the SDM models

      ###
      # Accounting for extensification
      species.df$esh.future.extent.sdm.rcp45.exp.int <-
        species.df$esh.values.45 -
        species.df$crop.values * species.df$crop.intens.tmp -
        species.df$pasture.values * species.df$past.intens.tmp -
        species.df$urban.2010.values * (tmp.urban.pref)

      species.df$esh.future.extent.sdm.rcp85.exp.int <-
        species.df$esh.values.85 -
        species.df$crop.values * species.df$crop.intens.tmp -
        species.df$pasture.values * species.df$past.intens.tmp -
        species.df$urban.2010.values * (tmp.urban.pref)
    }

    ###
    # Adding urban extent in
    if('urbanization' %in% stresses) {
      # Calculating loss from hab expansion
      species.df$esh.future.extent.esh.urb <-
        species.df$esh.values -
        species.df$crop.2010 * species.df$crop.intens.2010 -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.values * (tmp.urban.pref)

      ###
      # Repeating for SDM model under current climate
      if('climate' %in% stresses) {
        # Account for habitat loss from changes in crop extent
        species.df$esh.future.extent.sdm.urb <-
          species.df$esh.values.sdm -
          species.df$crop.2010 * species.df$crop.intens.2010 -
          species.df$past.2010 * species.df$past.intens.2010 -
          species.df$urban.values * (tmp.urban.pref)

        species.df$esh.future.extent.esh.sdm.urb <-
          tmp.esh.sdm -
          species.df$crop.2010 * species.df$crop.intens.2010 -
          species.df$past.2010 * species.df$past.intens.2010 -
          species.df$urban.values * (tmp.urban.pref)
      }
    } # End of if statement for urbanization

    # Unless extent decreases, in which intensity penalty in those cells is also ecreased
    if('expansion' %in% stresses &
       'urbanization' %in% stresses) {

      # Calculating loss from hab expansion
      # All ag changes
      species.df$esh.future.extent.esh.exp.urb <-
        species.df$esh.values -
        species.df$crop.values * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.2010 -
        species.df$urban.values * (tmp.urban.pref)

      # Crop only changes
      species.df$esh.future.extent.esh.cropexp.urb <-
        species.df$esh.values -
        species.df$crop.values * species.df$crop.intens.2010 -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.values * (tmp.urban.pref)

      # Pasture only changes
      species.df$esh.future.extent.esh.pastexp.urb <-
        species.df$esh.values -
        species.df$crop.2010 * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.2010 -
        species.df$urban.values * (tmp.urban.pref)

      ###
      # Repeating for SDM model under current climate
      if('climate' %in% stresses) {
        # Account for habitat loss from changes in crop extent
        species.df$esh.future.extent.sdm.exp.urb <-
          species.df$esh.values.sdm -
          species.df$crop.values * species.df$crop.intens.2010 -
          species.df$pasture.values * species.df$past.intens.2010 -
          species.df$urban.values * (tmp.urban.pref)

        species.df$esh.future.extent.esh.sdm.exp.urb <-
          tmp.esh.sdm -
          species.df$crop.values * species.df$crop.intens.2010 -
          species.df$pasture.values * species.df$past.intens.2010 -
          species.df$urban.values * (tmp.urban.pref)

      }
    } # End of function for expansion

    ###
    # Loop for ag intensification
    # This looks at ESH loss from ag intensification relative to current ESH
    # e.g. extent and intensity of current ag production is accounted for
    # But only change in intensification on existing cells is accounted for
    # Unless ag extent decreases in cells currently occupied by ag
    if('intensification' %in% stresses &
       'urbanization' %in% stresses) {
      # Isolating the impact of increased intensification on currently existing ag landscapes
      # This is for intensity of ag

      ###
      # First calculating habitat loss from on current ag extent
      species.df$esh.future.extent.esh.int.urb <-
        species.df$esh.values -
        species.df$crop.2010* species.df$crop.intens.tmp -
        species.df$past.2010 * species.df$past.intens.tmp -
        species.df$urban.values * (tmp.urban.pref)

      ###
      # Repeating for SDM model under current climate
      if('climate' %in% stresses) {
        # Isolating the impact of increased intensification on currently existing ag landscapes
        # This is for intensity of ag

        ###
        # First calculating esh loss from current cropland extent
        species.df$esh.future.extent.sdm.int.urb <-
          species.df$esh.values.sdm -
          species.df$crop.2010* species.df$crop.intens.tmp -
          species.df$past.2010 * species.df$past.intens.tmp -
          species.df$urban.values * (tmp.urban.pref)

        species.df$esh.future.extent.esh.sdm.int.urb <-
          tmp.esh.sdm -
          species.df$crop.2010* species.df$crop.intens.tmp -
          species.df$past.2010 * species.df$past.intens.tmp -
          species.df$urban.values * (tmp.urban.pref)
      }
    }

    ###
    # Loop for climate change
    # Shirnks/expands habitat based on outputs from SDM model
    if('climate' %in% stresses &
       'urbanization' %in% stresses) {
      ###
      # For current esh models, limit ag extent with current extent
      # First for RCP 45

      ###
      # Accounting for hab loss within historical crop extent
      species.df$esh.future.extent.esh.sdm.rcp45.urb <-
        tmp.esh.45 -
        species.df$crop.2010* species.df$crop.intens.2010 -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.values * (tmp.urban.pref)

      species.df$esh.future.extent.esh.sdm.rcp85.urb <-
        tmp.esh.85 -
        species.df$crop.2010* species.df$crop.intens.2010 -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.values * (tmp.urban.pref)



      ###
      # And repeating for SDM models
      # Don't need to clip the rasters beause we've already done this
      # And have already calculated esh loss from ag extent in 2010 under the current SDM
      species.df$esh.future.extent.sdm.rcp45.urb <-
        species.df$esh.values.45 -
        species.df$crop.2010* species.df$crop.intens.2010 -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.values * (tmp.urban.pref)

      species.df$esh.future.extent.sdm.rcp85.urb <-
        species.df$esh.values.85 -
        species.df$crop.2010* species.df$crop.intens.2010 -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.values * (tmp.urban.pref)
    }

    ###
    # Loop for climate change and expansion
    if('climate' %in% stresses &
       'expansion' %in% stresses &
       'urbanization' %in% stresses) {
      ###
      # We've already accounted for climate in the above if statement
      # Now accounting for expansion on top of climate

      ###
      # Now accounting for expansion of ag
      species.df$esh.future.extent.esh.sdm.rcp45.exp.urb <-
        tmp.esh.45 -
        species.df$crop.values * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.2010 -
        species.df$urban.values * (tmp.urban.pref)

      species.df$esh.future.extent.esh.sdm.rcp85.exp.urb <-
        tmp.esh.85 -
        species.df$crop.values * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.2010 -
        species.df$urban.values * (tmp.urban.pref)

      ###
      # Repeating for the SDM rasters
      # Second for the SDM rasters
      species.df$esh.future.extent.sdm.rcp45.exp.urb <-
        species.df$esh.values.45 -
        species.df$crop.values * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.2010 -
        species.df$urban.values * (tmp.urban.pref)

      species.df$esh.future.extent.sdm.rcp85.exp.urb <-
        species.df$esh.values.85 -
        species.df$crop.values * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.2010 -
        species.df$urban.values * (tmp.urban.pref)
    }

    ###
    # Loop for climate change and intensification
    if('climate' %in% stresses &
       'intensification' %in% stresses &
       'urbanization' %in% stresses) {
      ###
      # Updating ESH rasters for climate change

      ###
      # First for esh rasters
      # First step is to account for hab loss from ag extent in 2010
      species.df$esh.future.extent.esh.sdm.rcp45.int.urb <-
        tmp.esh.45 -
        species.df$crop.2010 * species.df$crop.intens.tmp -
        species.df$past.2010 * species.df$past.intens.tmp -
        species.df$urban.values * (tmp.urban.pref)

      species.df$esh.future.extent.esh.sdm.rcp85.int.urb <-
        tmp.esh.85 -
        species.df$crop.2010 * species.df$crop.intens.tmp -
        species.df$past.2010 * species.df$past.intens.tmp -
        species.df$urban.values * (tmp.urban.pref)

      ###
      # Second for sdm rasters
      # Don't need to adjust for climate because this is already done in the SDM models

      ###
      # First step is to account for hab loss from ag extent in 2010
      species.df$esh.future.extent.sdm.rcp45.int.urb <-
        species.df$esh.values.45 -
        species.df$crop.2010 * species.df$crop.intens.tmp -
        species.df$past.2010 * species.df$past.intens.tmp -
        species.df$urban.values * (tmp.urban.pref)

      species.df$esh.future.extent.sdm.rcp85.int.urb <-
        species.df$esh.values.85 -
        species.df$crop.2010 * species.df$crop.intens.tmp -
        species.df$past.2010 * species.df$past.intens.tmp -
        species.df$urban.values * (tmp.urban.pref)
    }

    ###
    # Loop for expansion and intensification
    if('expansion' %in% stresses &
       'intensification' %in% stresses &
       'urbanization' %in% stresses) {

      ###
      # ESH models first

      ###
      # First account for hab loss
      # All ag changes
      species.df$esh.future.extent.esh.exp.int.urb <-
        species.df$esh.values -
        species.df$crop.values * species.df$crop.intens.tmp -
        species.df$pasture.values * species.df$past.intens.tmp -
        species.df$urban.values * (tmp.urban.pref)

      # Crop only changes
      species.df$esh.future.extent.esh.cropexp.int.urb <-
        species.df$esh.values -
        species.df$crop.values * species.df$crop.intens.tmp -
        species.df$past.2010 * species.df$past.intens.2010 -
        species.df$urban.values * (tmp.urban.pref)

      # Pasture only changes
      species.df$esh.future.extent.esh.pastexp.int.urb <-
        species.df$esh.values -
        species.df$crop.2010 * species.df$crop.intens.2010 -
        species.df$pasture.values * species.df$past.intens.tmp -
        species.df$urban.values * (tmp.urban.pref)

      ###
      # Repeating for SDM model under current climate
      if(c('climate' %in% stresses)) {
        species.df$esh.future.extent.esh.sdm.exp.int.urb <-
          tmp.esh.sdm -
          species.df$crop.values * species.df$crop.intens.tmp -
          species.df$pasture.values * species.df$past.intens.tmp -
          species.df$urban.values * (tmp.urban.pref)

        species.df$esh.future.extent.sdm.exp.int.urb <-
          species.df$esh.values.sdm -
          species.df$crop.values * species.df$crop.intens.tmp -
          species.df$pasture.values * species.df$past.intens.tmp -
          species.df$urban.values * (tmp.urban.pref)
      }
    }

    ###
    # Loop for all three
    if('expansion' %in% stresses &
       'intensification' %in% stresses &
       'climate' %in% stresses &
       'urbanization' %in% stresses) {

      ###
      # First ESH models

      ###
      # Accounting for climate
      # Updating ESH rasters for climate change

      ###
      # Accounting for extensification
      species.df$esh.future.extent.esh.sdm.rcp45.exp.int.urb <-
        tmp.esh.45 -
        species.df$crop.values * species.df$crop.intens.tmp -
        species.df$pasture.values * species.df$past.intens.tmp -
        species.df$urban.values * (tmp.urban.pref)

      species.df$esh.future.extent.esh.sdm.rcp85.exp.int.urb <-
        tmp.esh.85 -
        species.df$crop.values * species.df$crop.intens.tmp -
        species.df$pasture.values * species.df$past.intens.tmp -
        species.df$urban.values * (tmp.urban.pref)


      ###
      # And repeating for the SDMs

      ###
      # Accounting for climate
      # This is already done in the SDM models

      ###
      # Accounting for extensification
      species.df$esh.future.extent.sdm.rcp45.exp.int.urb <-
        species.df$esh.values.45 -
        species.df$crop.values * species.df$crop.intens.tmp -
        species.df$pasture.values * species.df$past.intens.tmp -
        species.df$urban.values * (tmp.urban.pref)

      species.df$esh.future.extent.sdm.rcp85.exp.int.urb <-
        species.df$esh.values.85 -
        species.df$crop.values * species.df$crop.intens.tmp -
        species.df$pasture.values * species.df$past.intens.tmp -
        species.df$urban.values * (tmp.urban.pref)
    }

    ###
    # Restoring crop and pasture values
    species.df$crop.2010 <- store.crop.2010
    species.df$crop.values <- store.crop.values

    species.df$past.2010 <- store.past.2010
    species.df$pasture.values <- store.past.values


    # Removing objects to save space...
    # Exception for climate - don't want to throw a warning
    if('climate' %in% stresses) {
      rm(tmp.esh.45)
      rm(tmp.esh.85)
    }

    rm(store.crop.2010)
    rm(store.crop.values)
    rm(store.past.2010)
    rm(store.past.values)

    ###
    # And returning the data frame
    return(species.df)
  }


#####
###
# Writing a function to getpopulation density estimates based on climate variables and species life trait characteristics
density_estimate_function <- function(esh.raster_stack,
                                      coefs,
                                      coef_order,
                                      coef_family,
                                      coef_binomial,
                                      body_mass_grams,
                                      npp_raster,
                                      pcv_raster_stack,
                                      pwarmest_raster_stack,
                                      richness_raster,
                                      scenarios) {
  
  # Limiting stacks to only max values
  max_values = maxValue(esh.raster_stack)
  if(sum(max_values, na.rm = TRUE) >= 1) {
    
    # Coef family cannot be NA
    if(is.na(coef_family)) {
      coef_family <- 0
    }
    
    # Extracting values from richness raster
    tmp.richness = getValues(crop(richness_raster, esh.raster_stack))
    
    # and npp raster
    tmp.npp = getValues(npp_raster)
    
    
    esh.raster_stack = esh.raster_stack[[which(max_values %in% 1)]]
    pcv_raster_stack = pcv_raster_stack[[which(max_values %in% 1)]]
    pwarmest_raster_stack = pwarmest_raster_stack[[which(max_values %in% 1)]]
    scenarios = scenarios[which(max_values %in% 1)]
    
    # Cropping other rasters
    pcv_raster_stack <- crop(pcv_raster_stack, esh.raster_stack)
    pwarmest_raster_stack <- crop(pwarmest_raster_stack, esh.raster_stack)
    
    return.frame = data.frame()
    
    # Loop to get mean values across the scenarios
    for(z in 1:length(scenarios)) {
      # Extracting values from the rasters
      tmp.pcv = getValues(pcv_raster_stack[[z]])
      tmp.pwarmest = getValues(pwarmest_raster_stack[[z]])
      tmp.esh = getValues(esh.raster_stack[[z]])
      
      # Getting mean estimates
      npp.mean.tmp = mean(tmp.npp[!is.na(tmp.esh)], na.rm = TRUE)
      pcv.mean.tmp = mean(tmp.pcv[!is.na(tmp.esh)], na.rm = TRUE)
      pwarmest.mean.tmp = mean(pwarmest_raster_stack[!is.na(tmp.esh)], na.rm = TRUE)
      richness.mean.tmp = mean(tmp.richness[!is.na(tmp.esh)], na.rm = TRUE)
      
      # Getting median estimates
      npp.median.tmp = median(tmp.npp[!is.na(tmp.esh)], na.rm = TRUE)
      pcv.median.tmp = median(tmp.pcv[!is.na(tmp.esh)], na.rm = TRUE)
      pwarmest.median.tmp = median(pwarmest_raster_stack[!is.na(tmp.esh)], na.rm = TRUE)
      richness.median.tmp = median(tmp.richness[!is.na(tmp.esh)], na.rm = TRUE)
      
      # Getting standard errors
      npp.mean.se = sd(tmp.npp[!is.na(tmp.esh)], na.rm = TRUE) / length(tmp.npp[!is.na(tmp.esh)])
      pcv.mean.se = sd(tmp.pcv[!is.na(tmp.esh)], na.rm = TRUE) / length(tmp.pcv[!is.na(tmp.esh)])
      pwarmest.mean.se = sd(tmp.pwarmest[!is.na(tmp.esh)], na.rm = TRUE) / length(tmp.pwarmest[!is.na(tmp.esh)])
      richness.mean.se = sd(tmp.richness[!is.na(tmp.esh)], na.rm = TRUE) / length(tmp.richness[!is.na(tmp.esh)])
      
      # Putting estimates into a data frame
      if(z == 1 & !is.na(pcv.mean.tmp)) {
        return.frame = 
          data.frame(Scenario = scenarios[z],
                     Estimate = c('Mean','Median','se'),
                     npp = c(npp.mean.tmp, npp.median.tmp, npp.mean.se),
                     pcv = c(pcv.mean.tmp, pcv.median.tmp, pcv.mean.se),
                     pwarmest = c(pwarmest.mean.tmp, pwarmest.median.tmp, pwarmest.mean.se),
                     richness = c(richness.mean.tmp, richness.median.tmp, richness.mean.se))
      } else if(z >1 & !is.na(pcv.mean.tmp)) {
        return.frame = 
          rbind(return.frame,
                data.frame(Scenario = scenarios[z],
                           Estimate = c('Mean','Median','se'),
                           npp = c(npp.mean.tmp, npp.median.tmp, npp.mean.se),
                           pcv = c(pcv.mean.tmp, pcv.median.tmp, pcv.mean.se),
                           pwarmest = c(pwarmest.mean.tmp, pwarmest.median.tmp, pwarmest.mean.se),
                           richness = c(richness.mean.tmp, richness.median.tmp, richness.mean.se)))
      } else {
        return.frame = rbind(return.frame,
                             data.frame(Scenario = scenarios[z],
                                        Estimate = c('Mean','Median','se'),
                                        npp = NA,
                                        pcv = NA,
                                        pwarmest = NA,
                                        richness = NA))
      }
      
      
      
    } # End of for loop
    
    # And now estimating densities based on every possible combination of the variables
    # Within the different climate scenarios
    
    
    ###
    # Mean and median estimates
    return.frame <-
      return.frame %>%
      merge(.,coefs) %>%
      mutate(order_intercept_adj_weighted = coef_order,
             family_intercept_adj_weighted = coef_family,
             binomial_intercept_adj_weighted = coef_binomial,
             log10_body_mass_g = log10(body_mass_grams)) %>%
      mutate(predicted_coefficients = intercept_weighted_coef + 
               binomial_intercept_adj_weighted +
               family_intercept_adj_weighted + 
               order_intercept_adj_weighted +
               log10_body_mass_g * log10_body_mass_g_weighted_coef +
               log10_body_mass_g^2 * log10_body_mass_g2_weighted_coef + 
               log10_body_mass_g^3 * log10_body_mass_g3_weighted_coef +
               npp * npp_weighted_coef +
               npp^2 * npp2_weighted_coef +
               pcv * pcv_weighted_coef + 
               pcv^2 * pcv2_weighted_coef +
               pwarmest * pwarmest_weighted_coef +
               pwarmest^2 * pwarmest2_weighted_coef +
               richness * richness_weighted_coef) %>%
      mutate(Density = 10^(predicted_coefficients))
    
    
    mean.median.estimates <-
      return.frame %>%
      filter(Estimate %in% c('Mean','Median'))
    
    coef.cols <-
      which(names(return.frame) == 'pcv') :
      which(names(return.frame) == 'richness')
    
    # Looping to get every possible combination
    for(j in 1:length(scenarios)) {
      tmp.coefs <- 
        rbind(data.frame(return.frame %>% filter(Scenario %in% scenarios[j]) %>% filter(Estimate %in% c('Mean','se')) %>%
                           mutate(npp = npp[1] + 2*npp[2], pcv = pcv[1] + 2*pcv[2], pwarmest = pwarmest[1]+2*pwarmest[2], richness = richness[1] + 2*richness[2]) %>% .[1,]),
              data.frame(return.frame %>% filter(Scenario %in% scenarios[j]) %>% filter(Estimate %in% c('Mean','se')) %>%
                           mutate(npp = npp[1] - 2*npp[2], pcv = pcv[1] - 2*pcv[2], pwarmest = pwarmest[1]-2*pwarmest[2], richness = richness[1] - 2*richness[2]) %>% .[1,]),
              data.frame(return.frame %>% filter(Scenario %in% scenarios[j]) %>% filter(Estimate %in% c('Median','se')) %>%
                           mutate(npp = npp[1] + 2*npp[2], pcv = pcv[1] + 2*pcv[2], pwarmest = pwarmest[1]+2*pwarmest[2], richness = richness[1] + 2*richness[2]) %>% .[1,]),
              data.frame(return.frame %>% filter(Scenario %in% scenarios[j]) %>% filter(Estimate %in% c('Median','se')) %>%
                           mutate(npp = npp[1] - 2*npp[2], pcv = pcv[1] - 2*pcv[2], pwarmest = pwarmest[1]-2*pwarmest[2], richness = richness[1] - 2*richness[2]) %>% .[1,])) %>%
        mutate(Density_pcv = pcv_weighted_coef * pcv + pcv2_weighted_coef * pcv^2,
               Density_warmest = pwarmest_weighted_coef * pwarmest + pwarmest2_weighted_coef * pwarmest^2,
               Density_Richness = richness_weighted_coef * richness,
               Density_npp = npp_weighted_coef * npp + npp2_weighted_coef * npp^2) %>%
        mutate(Scenario_Est = c('Mean','Mean','Median','Median'))
      
      # And getting min and max density estimates
      estimate.frame <-
        data.frame(Scenario = scenarios[j],
                   Estimate = c('Min_Mean','Max_Mean','Min_Median','Max_Median'),
                   Density = c(10^(coefs$intercept_weighted_coef + coef_family +
                                     coefs$log10_body_mass_g_weighted_coef * log10(body_mass_grams) +
                                     coefs$log10_body_mass_g2_weighted_coef * log10((body_mass_grams))^2 +
                                     coefs$log10_body_mass_g3_weighted_coef * log10((body_mass_grams))^3 +
                                     min(tmp.coefs$Density_pcv[tmp.coefs$Scenario_Est %in% 'Mean']) + min(tmp.coefs$Density_Richness[tmp.coefs$Scenario_Est %in% 'Mean']) + min(tmp.coefs$Density_warmest[tmp.coefs$Scenario_Est %in% 'Mean']) + min(tmp.coefs$Density_npp[tmp.coefs$Scenario_Est %in% 'Mean'])),
                               10^(coefs$intercept_weighted_coef + coef_family +
                                     coefs$log10_body_mass_g_weighted_coef * log10(body_mass_grams) +
                                     coefs$log10_body_mass_g2_weighted_coef * log10((body_mass_grams))^2 +
                                     coefs$log10_body_mass_g3_weighted_coef * log10((body_mass_grams))^3 +
                                     max(tmp.coefs$Density_pcv[tmp.coefs$Scenario_Est %in% 'Mean']) + max(tmp.coefs$Density_Richness[tmp.coefs$Scenario_Est %in% 'Mean']) + max(tmp.coefs$Density_warmest[tmp.coefs$Scenario_Est %in% 'Mean'] + max(tmp.coefs$Density_npp[tmp.coefs$Scenario_Est %in% 'Mean']))),
                               10^(coefs$intercept_weighted_coef + coef_family +
                                     coefs$log10_body_mass_g_weighted_coef * log10(body_mass_grams) +
                                     coefs$log10_body_mass_g2_weighted_coef * log10((body_mass_grams))^2 +
                                     coefs$log10_body_mass_g3_weighted_coef * log10((body_mass_grams))^3 +
                                     min(tmp.coefs$Density_pcv[tmp.coefs$Scenario_Est %in% 'Median']) + min(tmp.coefs$Density_Richness[tmp.coefs$Scenario_Est %in% 'Median']) + min(tmp.coefs$Density_warmest[tmp.coefs$Scenario_Est %in% 'Median']) + min(tmp.coefs$Density_npp[tmp.coefs$Scenario_Est %in% 'Median'])),
                               10^(coefs$intercept_weighted_coef + coef_family +
                                     coefs$log10_body_mass_g_weighted_coef * log10(body_mass_grams) +
                                     coefs$log10_body_mass_g2_weighted_coef * log10((body_mass_grams))^2 +
                                     coefs$log10_body_mass_g3_weighted_coef * log10((body_mass_grams))^3 +
                                     max(tmp.coefs$Density_pcv[tmp.coefs$Scenario_Est %in% 'Median']) + max(tmp.coefs$Density_Richness[tmp.coefs$Scenario_Est %in% 'Median']) + max(tmp.coefs$Density_warmest[tmp.coefs$Scenario_Est %in% 'Median']) + max(tmp.coefs$Density_npp[tmp.coefs$Scenario_Est %in% 'Median']))))
      
      if(j == 1) {
        estimate.out <- estimate.frame
      } else {
        estimate.out <- 
          rbind(estimate.out,estimate.frame)
      }
    }
    
    
    density.frame <-
      rbind(dplyr::select(mean.median.estimates, Scenario, Estimate, Density),
            estimate.out)
  } else {
    density.frame = data.frame(Scenario = NA, Estimate = NA, Density = NA)
  }
  
  
  
  # 
  return(density.frame)
}


# For mammals
density_estimate_function_mammals <- function(esh.raster_stack,
                                              coefs,
                                              coef_order,
                                              coef_family,
                                              coef_binomial,
                                              body_mass_grams,
                                              npp_raster,
                                              pcv_raster_stack,
                                              pwarmest_raster_stack,
                                              richness_raster,
                                              scenarios) {
  
  # Limiting stacks to only max values
  max_values = maxValue(esh.raster_stack)
  if(sum(max_values, na.rm = TRUE) >= 1) {
    
    # Coefficients cannot be NA
    if(is.na(coef_order)) {
      coef_order <- 0
    }
    if(is.na(coef_family)) {
      coef_family <- 0
    }
    if(is.na(coef_binomial)) {
      coef_binomial <- 0
    }
    
    # Extracting values from richness raster
    tmp.richness = getValues(crop(richness_raster, esh.raster_stack))
    
    # and npp raster
    tmp.npp = getValues(npp_raster)
    
    
    esh.raster_stack = esh.raster_stack[[which(max_values %in% 1)]]
    pcv_raster_stack = pcv_raster_stack[[which(max_values %in% 1)]]
    pwarmest_raster_stack = pwarmest_raster_stack[[which(max_values %in% 1)]]
    scenarios = scenarios[which(max_values %in% 1)]
    
    # Cropping other rasters
    pcv_raster_stack <- crop(pcv_raster_stack, esh.raster_stack)
    pwarmest_raster_stack <- crop(pwarmest_raster_stack, esh.raster_stack)
    
    return.frame = data.frame()
    
    # Loop to get mean values across the scenarios
    for(z in 1:length(scenarios)) {
      # Extracting values from the rasters
      tmp.pcv = getValues(pcv_raster_stack[[z]])
      tmp.pwarmest = getValues(pwarmest_raster_stack[[z]])
      tmp.esh = getValues(esh.raster_stack[[z]])
      
      # Getting mean estimates
      npp.mean.tmp = mean(tmp.npp[!is.na(tmp.esh)], na.rm = TRUE)
      pcv.mean.tmp = mean(tmp.pcv[!is.na(tmp.esh)], na.rm = TRUE)
      pwarmest.mean.tmp = mean(pwarmest_raster_stack[!is.na(tmp.esh)], na.rm = TRUE)
      richness.mean.tmp = mean(tmp.richness[!is.na(tmp.esh)], na.rm = TRUE)
      
      # Getting median estimates
      npp.median.tmp = median(tmp.npp[!is.na(tmp.esh)], na.rm = TRUE)
      pcv.median.tmp = median(tmp.pcv[!is.na(tmp.esh)], na.rm = TRUE)
      pwarmest.median.tmp = median(pwarmest_raster_stack[!is.na(tmp.esh)], na.rm = TRUE)
      richness.median.tmp = median(tmp.richness[!is.na(tmp.esh)], na.rm = TRUE)
      
      # Getting standard errors
      npp.mean.se = sd(tmp.npp[!is.na(tmp.esh)], na.rm = TRUE) / length(tmp.npp[!is.na(tmp.esh)])
      pcv.mean.se = sd(tmp.pcv[!is.na(tmp.esh)], na.rm = TRUE) / length(tmp.pcv[!is.na(tmp.esh)])
      pwarmest.mean.se = sd(tmp.pwarmest[!is.na(tmp.esh)], na.rm = TRUE) / length(tmp.pwarmest[!is.na(tmp.esh)])
      richness.mean.se = sd(tmp.richness[!is.na(tmp.esh)], na.rm = TRUE) / length(tmp.richness[!is.na(tmp.esh)])
      
      # Putting estimates into a data frame
      if(z == 1 & !is.na(pcv.mean.tmp)) {
        return.frame = 
          data.frame(Scenario = scenarios[z],
                     Estimate = c('Mean','Median','se'),
                     npp = c(npp.mean.tmp, npp.median.tmp, npp.mean.se),
                     pcv = c(pcv.mean.tmp, pcv.median.tmp, pcv.mean.se),
                     pwarmest = c(pwarmest.mean.tmp, pwarmest.median.tmp, pwarmest.mean.se),
                     richness = c(richness.mean.tmp, richness.median.tmp, richness.mean.se))
      } else if(z >1 & !is.na(pcv.mean.tmp)) {
        return.frame = 
          rbind(return.frame,
                data.frame(Scenario = scenarios[z],
                           Estimate = c('Mean','Median','se'),
                           npp = c(npp.mean.tmp, npp.median.tmp, npp.mean.se),
                           pcv = c(pcv.mean.tmp, pcv.median.tmp, pcv.mean.se),
                           pwarmest = c(pwarmest.mean.tmp, pwarmest.median.tmp, pwarmest.mean.se),
                           richness = c(richness.mean.tmp, richness.median.tmp, richness.mean.se)))
      } else {
        return.frame = rbind(return.frame,
                             data.frame(Scenario = scenarios[z],
                                        Estimate = c('Mean','Median','se'),
                                        npp = NA,
                                        pcv = NA,
                                        pwarmest = NA,
                                        richness = NA))
      }
      
      
      
    } # End of for loop
    
    # And now estimating densities based on every possible combination of the variables
    # Within the different climate scenarios
    
    
    ###
    # Mean and median estimates
    return.frame <-
      return.frame %>%
      merge(.,coefs) %>%
      mutate(family_intercept_adj_weighted = coef_family,
             order_intercept_adj_weighted = coef_order,
             binomial_intercept_adj_weighted = coef_binomial,
             log10_body_mass_g = log10(body_mass_grams)) %>%
      mutate(predicted_coefficients = intercept_weighted_coef + 
               binomial_intercept_adj_weighted +
               family_intercept_adj_weighted + 
               order_intercept_adj_weighted +
               log10_body_mass_g * log10_body_mass_g_weighted_coef +
               log10_body_mass_g^2 * log10_body_mass_g2_weighted_coef + 
               log10_body_mass_g^3 * log10_body_mass_g3_weighted_coef +
               npp * npp_weighted_coef +
               npp^2 * npp2_weighted_coef +
               pcv * pcv_weighted_coef + 
               pcv^2 * pcv2_weighted_coef +
               pwarmest * pwarmest_weighted_coef +
               pwarmest^2 * pwarmest2_weighted_coef +
               richness * richness_mammals_weighted_coef) %>%
      mutate(Density = 10^(predicted_coefficients))
    
    
    mean.median.estimates <-
      return.frame %>%
      filter(Estimate %in% c('Mean','Median'))
    
    coef.cols <-
      which(names(return.frame) == 'pcv') :
      which(names(return.frame) == 'richness')
    
    # Looping to get every possible combination
    for(j in 1:length(scenarios)) {
      tmp.coefs <- 
        rbind(data.frame(return.frame %>% filter(Scenario %in% scenarios[j]) %>% filter(Estimate %in% c('Mean','se')) %>%
                           mutate(npp = npp[1] + 2*npp[2], pcv = pcv[1] + 2*pcv[2], pwarmest = pwarmest[1]+2*pwarmest[2], richness = richness[1] + 2*richness[2]) %>% .[1,]),
              data.frame(return.frame %>% filter(Scenario %in% scenarios[j]) %>% filter(Estimate %in% c('Mean','se')) %>%
                           mutate(npp = npp[1] - 2*npp[2], pcv = pcv[1] - 2*pcv[2], pwarmest = pwarmest[1]-2*pwarmest[2], richness = richness[1] - 2*richness[2]) %>% .[1,]),
              data.frame(return.frame %>% filter(Scenario %in% scenarios[j]) %>% filter(Estimate %in% c('Median','se')) %>%
                           mutate(npp = npp[1] + 2*npp[2], pcv = pcv[1] + 2*pcv[2], pwarmest = pwarmest[1]+2*pwarmest[2], richness = richness[1] + 2*richness[2]) %>% .[1,]),
              data.frame(return.frame %>% filter(Scenario %in% scenarios[j]) %>% filter(Estimate %in% c('Median','se')) %>%
                           mutate(npp = npp[1] - 2*npp[2], pcv = pcv[1] - 2*pcv[2], pwarmest = pwarmest[1]-2*pwarmest[2], richness = richness[1] - 2*richness[2]) %>% .[1,])) %>%
        mutate(Density_pcv = pcv_weighted_coef * pcv + pcv2_weighted_coef * pcv^2,
               Density_warmest = pwarmest_weighted_coef * pwarmest + pwarmest2_weighted_coef * pwarmest^2,
               Density_Richness = richness_mammals_weighted_coef * richness,
               Density_npp = npp_weighted_coef * npp + npp2_weighted_coef * npp^2) %>%
        mutate(Scenario_Est = c('Mean','Mean','Median','Median'))
      
      # And getting min and max density estimates
      estimate.frame <-
        data.frame(Scenario = scenarios[j],
                   Estimate = c('Min_Mean','Max_Mean','Min_Median','Max_Median'),
                   Density = c(10^(coefs$intercept_weighted_coef + coef_family + coef_order + coef_binomial +
                                     coefs$log10_body_mass_g_weighted_coef * log10(body_mass_grams) +
                                     coefs$log10_body_mass_g2_weighted_coef * log10((body_mass_grams))^2 +
                                     coefs$log10_body_mass_g3_weighted_coef * log10((body_mass_grams))^3 +
                                     min(tmp.coefs$Density_pcv[tmp.coefs$Scenario_Est %in% 'Mean']) + min(tmp.coefs$Density_Richness[tmp.coefs$Scenario_Est %in% 'Mean']) + min(tmp.coefs$Density_warmest[tmp.coefs$Scenario_Est %in% 'Mean']) + min(tmp.coefs$Density_npp[tmp.coefs$Scenario_Est %in% 'Mean'])),
                               10^(coefs$intercept_weighted_coef + coef_family + coef_order + coef_binomial +
                                     coefs$log10_body_mass_g_weighted_coef * log10(body_mass_grams) +
                                     coefs$log10_body_mass_g2_weighted_coef * log10((body_mass_grams))^2 +
                                     coefs$log10_body_mass_g3_weighted_coef * log10((body_mass_grams))^3 +
                                     max(tmp.coefs$Density_pcv[tmp.coefs$Scenario_Est %in% 'Mean']) + max(tmp.coefs$Density_Richness[tmp.coefs$Scenario_Est %in% 'Mean']) + max(tmp.coefs$Density_warmest[tmp.coefs$Scenario_Est %in% 'Mean'] + max(tmp.coefs$Density_npp[tmp.coefs$Scenario_Est %in% 'Mean']))),
                               10^(coefs$intercept_weighted_coef + coef_family + coef_order + coef_binomial +
                                     coefs$log10_body_mass_g_weighted_coef * log10(body_mass_grams) +
                                     coefs$log10_body_mass_g2_weighted_coef * log10((body_mass_grams))^2 +
                                     coefs$log10_body_mass_g3_weighted_coef * log10((body_mass_grams))^3 +
                                     min(tmp.coefs$Density_pcv[tmp.coefs$Scenario_Est %in% 'Median']) + min(tmp.coefs$Density_Richness[tmp.coefs$Scenario_Est %in% 'Median']) + min(tmp.coefs$Density_warmest[tmp.coefs$Scenario_Est %in% 'Median']) + min(tmp.coefs$Density_npp[tmp.coefs$Scenario_Est %in% 'Median'])),
                               10^(coefs$intercept_weighted_coef + coef_family + coef_order + coef_binomial +
                                     coefs$log10_body_mass_g_weighted_coef * log10(body_mass_grams) +
                                     coefs$log10_body_mass_g2_weighted_coef * log10((body_mass_grams))^2 +
                                     coefs$log10_body_mass_g3_weighted_coef * log10((body_mass_grams))^3 +
                                     max(tmp.coefs$Density_pcv[tmp.coefs$Scenario_Est %in% 'Median']) + max(tmp.coefs$Density_Richness[tmp.coefs$Scenario_Est %in% 'Median']) + max(tmp.coefs$Density_warmest[tmp.coefs$Scenario_Est %in% 'Median']) + max(tmp.coefs$Density_npp[tmp.coefs$Scenario_Est %in% 'Median']))))
      
      if(j == 1) {
        estimate.out <- estimate.frame
      } else {
        estimate.out <- 
          rbind(estimate.out,estimate.frame)
      }
    }
    
    
    density.frame <-
      rbind(dplyr::select(mean.median.estimates, Scenario, Estimate, Density),
            estimate.out)
  } else {
    density.frame = data.frame(Scenario = NA, Estimate = NA, Density = NA)
  }
  
  
  
  # 
  return(density.frame)
}

# For amphibians
density_estimate_function_amphibians <- function(esh.raster_stack,
                                                 coefs,
                                                 coef_order,
                                                 coef_family,
                                                 coef_binomial,
                                                 body_mass_grams,
                                                 npp_raster,
                                                 pcv_raster_stack,
                                                 pwarmest_raster_stack,
                                                 richness_raster,
                                                 scenarios) {
  
  # Limiting stacks to only max values
  max_values = maxValue(esh.raster_stack)
  if(sum(max_values, na.rm = TRUE) >= 1) {
    
    # Coefficients cannot be NA
    if(is.na(coef_order)) {
      coef_order <- 0
    }
    if(is.na(coef_family)) {
      coef_family <- 0
    }
    if(is.na(coef_binomial)) {
      coef_binomial <- 0
    }
    
    # Extracting values from richness raster
    tmp.richness = getValues(crop(richness_raster, esh.raster_stack))
    
    # and npp raster
    tmp.npp = getValues(npp_raster)
    
    
    esh.raster_stack = esh.raster_stack[[which(max_values %in% 1)]]
    pcv_raster_stack = pcv_raster_stack[[which(max_values %in% 1)]]
    pwarmest_raster_stack = pwarmest_raster_stack[[which(max_values %in% 1)]]
    scenarios = scenarios[which(max_values %in% 1)]
    
    # Cropping other rasters
    pcv_raster_stack <- crop(pcv_raster_stack, esh.raster_stack)
    pwarmest_raster_stack <- crop(pwarmest_raster_stack, esh.raster_stack)
    
    return.frame = data.frame()
    
    # Loop to get mean values across the scenarios
    for(z in 1:length(scenarios)) {
      # Extracting values from the rasters
      tmp.pcv = getValues(pcv_raster_stack[[z]])
      tmp.pwarmest = getValues(pwarmest_raster_stack[[z]])
      tmp.esh = getValues(esh.raster_stack[[z]])
      
      # Getting mean estimates
      npp.mean.tmp = mean(tmp.npp[!is.na(tmp.esh)], na.rm = TRUE)
      pcv.mean.tmp = mean(tmp.pcv[!is.na(tmp.esh)], na.rm = TRUE)
      pwarmest.mean.tmp = mean(pwarmest_raster_stack[!is.na(tmp.esh)], na.rm = TRUE)
      richness.mean.tmp = mean(tmp.richness[!is.na(tmp.esh)], na.rm = TRUE)
      
      # Getting median estimates
      npp.median.tmp = median(tmp.npp[!is.na(tmp.esh)], na.rm = TRUE)
      pcv.median.tmp = median(tmp.pcv[!is.na(tmp.esh)], na.rm = TRUE)
      pwarmest.median.tmp = median(pwarmest_raster_stack[!is.na(tmp.esh)], na.rm = TRUE)
      richness.median.tmp = median(tmp.richness[!is.na(tmp.esh)], na.rm = TRUE)
      
      # Getting standard deviations
      # Ignore the variable naming...
      npp.mean.se = sd(tmp.npp[!is.na(tmp.esh)], na.rm = TRUE) #/ length(tmp.npp[!is.na(tmp.esh)])
      pcv.mean.se = sd(tmp.pcv[!is.na(tmp.esh)], na.rm = TRUE) #/ length(tmp.pcv[!is.na(tmp.esh)])
      pwarmest.mean.se = sd(tmp.pwarmest[!is.na(tmp.esh)], na.rm = TRUE) #/ length(tmp.pwarmest[!is.na(tmp.esh)])
      richness.mean.se = sd(tmp.richness[!is.na(tmp.esh)], na.rm = TRUE) #/ length(tmp.richness[!is.na(tmp.esh)])
      
      # Putting estimates into a data frame
      if(z == 1 & !is.na(pcv.mean.tmp)) {
        return.frame = 
          data.frame(Scenario = scenarios[z],
                     Estimate = c('Mean','Median','se'),
                     npp = c(npp.mean.tmp, npp.median.tmp, npp.mean.se),
                     pcv = c(pcv.mean.tmp, pcv.median.tmp, pcv.mean.se),
                     pwarmest = c(pwarmest.mean.tmp, pwarmest.median.tmp, pwarmest.mean.se),
                     richness = c(richness.mean.tmp, richness.median.tmp, richness.mean.se))
      } else if(z >1 & !is.na(pcv.mean.tmp)) {
        return.frame = 
          rbind(return.frame,
                data.frame(Scenario = scenarios[z],
                           Estimate = c('Mean','Median','se'),
                           npp = c(npp.mean.tmp, npp.median.tmp, npp.mean.se),
                           pcv = c(pcv.mean.tmp, pcv.median.tmp, pcv.mean.se),
                           pwarmest = c(pwarmest.mean.tmp, pwarmest.median.tmp, pwarmest.mean.se),
                           richness = c(richness.mean.tmp, richness.median.tmp, richness.mean.se)))
      } else {
        return.frame = rbind(return.frame,
                             data.frame(Scenario = scenarios[z],
                                        Estimate = c('Mean','Median','se'),
                                        npp = NA,
                                        pcv = NA,
                                        pwarmest = NA,
                                        richness = NA))
      }
      
      
      
    } # End of for loop
    
    # And now estimating densities based on every possible combination of the variables
    # Within the different climate scenarios
    
    
    ###
    # Mean and median estimates
    return.frame <-
      return.frame %>%
      merge(.,coefs) %>%
      mutate(family_intercept_adj_weighted = coef_family,
             order_intercept_adj_weighted = coef_order,
             binomial_intercept_adj_weighted = coef_binomial,
             log10_body_mass_g = log10(body_mass_grams)) %>%
      mutate(predicted_coefficients = intercept_weighted_coef + 
               binomial_intercept_adj_weighted +
               family_intercept_adj_weighted + 
               order_intercept_adj_weighted +
               log10_body_mass_g * log10_body_mass_g_weighted_coef +
               log10_body_mass_g^2 * log10_body_mass_g2_weighted_coef + 
               log10_body_mass_g^3 * log10_body_mass_g3_weighted_coef +
               npp * npp_weighted_coef +
               npp^2 * npp2_weighted_coef +
               pcv * pcv_weighted_coef + 
               pcv^2 * pcv2_weighted_coef +
               pwarmest * pwarmest_weighted_coef +
               pwarmest^2 * pwarmest2_weighted_coef +
               richness * richness_amps_weighted_coef) %>%
      mutate(Density = 10^(predicted_coefficients))
    
    
    mean.median.estimates <-
      return.frame %>%
      filter(Estimate %in% c('Mean','Median'))
    
    coef.cols <-
      which(names(return.frame) == 'pcv') :
      which(names(return.frame) == 'richness')
    
    # Looping to get every possible combination
    for(j in 1:length(scenarios)) {
      tmp.coefs <- 
        rbind(data.frame(return.frame %>% filter(Scenario %in% scenarios[j]) %>% filter(Estimate %in% c('Mean','se')) %>%
                           mutate(npp = npp[1] + 2*npp[2], pcv = pcv[1] + 2*pcv[2], pwarmest = pwarmest[1]+2*pwarmest[2], richness = richness[1] + 2*richness[2]) %>% .[1,]),
              data.frame(return.frame %>% filter(Scenario %in% scenarios[j]) %>% filter(Estimate %in% c('Mean','se')) %>%
                           mutate(npp = npp[1] - 2*npp[2], pcv = pcv[1] - 2*pcv[2], pwarmest = pwarmest[1]-2*pwarmest[2], richness = richness[1] - 2*richness[2]) %>% .[1,]),
              data.frame(return.frame %>% filter(Scenario %in% scenarios[j]) %>% filter(Estimate %in% c('Median','se')) %>%
                           mutate(npp = npp[1] + 2*npp[2], pcv = pcv[1] + 2*pcv[2], pwarmest = pwarmest[1]+2*pwarmest[2], richness = richness[1] + 2*richness[2]) %>% .[1,]),
              data.frame(return.frame %>% filter(Scenario %in% scenarios[j]) %>% filter(Estimate %in% c('Median','se')) %>%
                           mutate(npp = npp[1] - 2*npp[2], pcv = pcv[1] - 2*pcv[2], pwarmest = pwarmest[1]-2*pwarmest[2], richness = richness[1] - 2*richness[2]) %>% .[1,])) %>%
        mutate(Density_pcv = pcv_weighted_coef * pcv + pcv2_weighted_coef * pcv^2,
               Density_warmest = pwarmest_weighted_coef * pwarmest + pwarmest2_weighted_coef * pwarmest^2,
               Density_Richness = richness_amps_weighted_coef * richness,
               Density_npp = npp_weighted_coef * npp + npp2_weighted_coef * npp^2) %>%
        mutate(Scenario_Est = c('Mean','Mean','Median','Median'))
      
      # And getting min and max density estimates
      estimate.frame <-
        data.frame(Scenario = scenarios[j],
                   Estimate = c('Min_Mean','Max_Mean','Min_Median','Max_Median'),
                   Density = c(10^(coefs$intercept_weighted_coef + coef_family + coef_order + coef_binomial +
                                     coefs$log10_body_mass_g_weighted_coef * log10(body_mass_grams) +
                                     coefs$log10_body_mass_g2_weighted_coef * log10((body_mass_grams))^2 +
                                     coefs$log10_body_mass_g3_weighted_coef * log10((body_mass_grams))^3 +
                                     min(tmp.coefs$Density_pcv[tmp.coefs$Scenario_Est %in% 'Mean']) + min(tmp.coefs$Density_Richness[tmp.coefs$Scenario_Est %in% 'Mean']) + min(tmp.coefs$Density_warmest[tmp.coefs$Scenario_Est %in% 'Mean']) + min(tmp.coefs$Density_npp[tmp.coefs$Scenario_Est %in% 'Mean'])),
                               10^(coefs$intercept_weighted_coef + coef_family + coef_order + coef_binomial +
                                     coefs$log10_body_mass_g_weighted_coef * log10(body_mass_grams) +
                                     coefs$log10_body_mass_g2_weighted_coef * log10((body_mass_grams))^2 +
                                     coefs$log10_body_mass_g3_weighted_coef * log10((body_mass_grams))^3 +
                                     max(tmp.coefs$Density_pcv[tmp.coefs$Scenario_Est %in% 'Mean']) + max(tmp.coefs$Density_Richness[tmp.coefs$Scenario_Est %in% 'Mean']) + max(tmp.coefs$Density_warmest[tmp.coefs$Scenario_Est %in% 'Mean'] + max(tmp.coefs$Density_npp[tmp.coefs$Scenario_Est %in% 'Mean']))),
                               10^(coefs$intercept_weighted_coef + coef_family + coef_order + coef_binomial +
                                     coefs$log10_body_mass_g_weighted_coef * log10(body_mass_grams) +
                                     coefs$log10_body_mass_g2_weighted_coef * log10((body_mass_grams))^2 +
                                     coefs$log10_body_mass_g3_weighted_coef * log10((body_mass_grams))^3 +
                                     min(tmp.coefs$Density_pcv[tmp.coefs$Scenario_Est %in% 'Median']) + min(tmp.coefs$Density_Richness[tmp.coefs$Scenario_Est %in% 'Median']) + min(tmp.coefs$Density_warmest[tmp.coefs$Scenario_Est %in% 'Median']) + min(tmp.coefs$Density_npp[tmp.coefs$Scenario_Est %in% 'Median'])),
                               10^(coefs$intercept_weighted_coef + coef_family + coef_order + coef_binomial +
                                     coefs$log10_body_mass_g_weighted_coef * log10(body_mass_grams) +
                                     coefs$log10_body_mass_g2_weighted_coef * log10((body_mass_grams))^2 +
                                     coefs$log10_body_mass_g3_weighted_coef * log10((body_mass_grams))^3 +
                                     max(tmp.coefs$Density_pcv[tmp.coefs$Scenario_Est %in% 'Median']) + max(tmp.coefs$Density_Richness[tmp.coefs$Scenario_Est %in% 'Median']) + max(tmp.coefs$Density_warmest[tmp.coefs$Scenario_Est %in% 'Median']) + max(tmp.coefs$Density_npp[tmp.coefs$Scenario_Est %in% 'Median']))))
      
      if(j == 1) {
        estimate.out <- estimate.frame
      } else {
        estimate.out <- 
          rbind(estimate.out,estimate.frame)
      }
    }
    
    
    density.frame <-
      rbind(dplyr::select(mean.median.estimates, Scenario, Estimate, Density),
            estimate.out)
  } else {
    density.frame = data.frame(Scenario = NA, Estimate = NA, Density = NA)
  }
  
  
  
  # 
  return(density.frame)
}

new_density_function = function(coefs,
                                coef_order,
                                coef_family,
                                coef_binomial,
                                body_mass_grams,
                                npp_input,
                                pcv_input,
                                pwarmest_input,
                                richness_input,
                                habitat_input,
                                hab_min,
                                hab_max,
                                patches_migrate,
                                patches_no_migrate) {
  
  tmp.df = 
    data.frame(npp = getValues(npp_input),
               pcv = getValues(pcv_input),
               pwarmest = getValues(pwarmest_input),
               richness = getValues(richness_input),
               habitat = getValues(habitat_input)) %>%
    mutate(Density = 
             10^(coefs$intercept_weighted_coef + 
                   coef_order + coef_family + coef_binomial + 
                   log10(body_mass_grams) * coefs$log10_body_mass_g_weighted_coef +
                   log10(body_mass_grams)^2 * coefs$log10_body_mass_g2_weighted_coef + 
                   log10(body_mass_grams)^3 * coefs$log10_body_mass_g3_weighted_coef +
                   npp * coefs$npp_weighted_coef +
                   npp^2 * coefs$npp2_weighted_coef +
                   pcv * coefs$pcv_weighted_coef + 
                   pcv^2 * coefs$pcv2_weighted_coef +
                   pwarmest * coefs$pwarmest_weighted_coef +
                   pwarmest^2 * coefs$pwarmest2_weighted_coef +
                   richness * coefs$richness_amps_weighted_coef)) %>%
    mutate(Density = squish(Density, range = c(hab_min, hab_max))) %>%
    mutate(Pop_estimate = Density * habitat * 2.25) %>%
    dplyr::select(.,Pop_estimate, habitat) %>%
    mutate(patch_migrate = getValues(patches_migrate),
           patch_no_migrate = getValues(patches_no_migrate),
           count = 1) %>%
    rowsum(.,group = .$patch_no_migrate,na.rm=TRUE) %>%
    mutate(patch_no_migrate = patch_no_migrate / count,
           patch_migrate = patch_migrate/count) %>%
    dplyr::select(.,
                  patch_no_migrate,patch_migrate,
                  habitat, estimated_population = Pop_estimate,
                  n_cells = count) %>%
    filter(patch_no_migrate > 0)
  
  
  return(tmp.df)
  # Inputs no longer need to be a raster stack, can just pass individual rasters
  # Need to pass patch raster (migrate + no migrate) and species loss raster
  # Create data frame, multiply env variables by patch size
  # Row sums group by patch to get total value of varaibles, then divide by total patch to get weighted average by patch size
  # Run these estimates through density coefs by patch
  # And return this data frame
  
  #
  
  #
  
  # Output is data frame with three coumns with n rows, where n = number of patches
  # Where columns indicate mean estimate, low estimate, and high density estimate
  # Each row indicates estimate for different patch
  
  
}





# Function to make raster with pop density estimates----
pop_density_function = 
  function(coefs,
           coef_order,
           coef_family,
           coef_binomial,
           body_mass_grams,
           npp_input,
           pcv_input,
           pwarmest_input,
           richness_input,
           hab_min,
           hab_max,
           nrows) {
    tmp.df = 
      data.frame(npp = getValues(npp_input),
                 pcv = getValues(pcv_input),
                 pwarmest = getValues(pwarmest_input),
                 richness = getValues(richness_input)) %>%
      mutate(Density = 
               10^(coefs$intercept_weighted_coef + 
                     coef_order + coef_family + coef_binomial + 
                     log10(body_mass_grams) * coefs$log10_body_mass_g_weighted_coef +
                     log10(body_mass_grams)^2 * coefs$log10_body_mass_g2_weighted_coef + 
                     log10(body_mass_grams)^3 * coefs$log10_body_mass_g3_weighted_coef +
                     npp * coefs$npp_weighted_coef +
                     npp^2 * coefs$npp2_weighted_coef +
                     pcv * coefs$pcv_weighted_coef + 
                     pcv^2 * coefs$pcv2_weighted_coef +
                     pwarmest * coefs$pwarmest_weighted_coef +
                     pwarmest^2 * coefs$pwarmest2_weighted_coef +
                     richness * coefs$richness_weighted_coef)) %>%
      mutate(Density = squish(Density, range = c(hab_min, hab_max)))
    
    density_raster = raster(matrix(tmp.df$Density, byrow = TRUE, nrow = nrows))
    
    return(density_raster)
  }

# For amphibians
pop_density_function_amps = 
  function(coefs,
           coef_order,
           coef_family,
           coef_binomial,
           body_mass_grams,
           npp_input,
           pcv_input,
           pwarmest_input,
           richness_input,
           hab_min,
           hab_max,
           nrows) {
    tmp.df = 
      data.frame(npp = getValues(npp_input),
                 pcv = getValues(pcv_input),
                 pwarmest = getValues(pwarmest_input),
                 richness = getValues(richness_input)) %>%
      mutate(Density = 
               10^(coefs$intercept_weighted_coef + 
                     coef_order + coef_family + coef_binomial + 
                     log10(body_mass_grams) * coefs$log10_body_mass_g_weighted_coef +
                     log10(body_mass_grams)^2 * coefs$log10_body_mass_g2_weighted_coef + 
                     log10(body_mass_grams)^3 * coefs$log10_body_mass_g3_weighted_coef +
                     npp * coefs$npp_weighted_coef +
                     npp^2 * coefs$npp2_weighted_coef +
                     pcv * coefs$pcv_weighted_coef + 
                     pcv^2 * coefs$pcv2_weighted_coef +
                     pwarmest * coefs$pwarmest_weighted_coef +
                     pwarmest^2 * coefs$pwarmest2_weighted_coef +
                     richness * coefs$richness_amps_weighted_coef)) %>%
      mutate(Density = squish(Density, range = c(hab_min, hab_max)))
    
    density_raster = raster(matrix(tmp.df$Density, byrow = TRUE, nrow = nrows))
    
    return(density_raster)
  }

# For mammals
pop_density_function_mammals = 
  function(coefs,
           coef_order,
           coef_family,
           coef_binomial,
           body_mass_grams,
           npp_input,
           pcv_input,
           pwarmest_input,
           richness_input,
           hab_min,
           hab_max,
           nrows) {
    tmp.df = 
      data.frame(npp = getValues(npp_input),
                 pcv = getValues(pcv_input),
                 pwarmest = getValues(pwarmest_input),
                 richness = getValues(richness_input)) %>%
      mutate(Density = 
               10^(coefs$intercept_weighted_coef + 
                     coef_order + coef_family + coef_binomial + 
                     log10(body_mass_grams) * coefs$log10_body_mass_g_weighted_coef +
                     log10(body_mass_grams)^2 * coefs$log10_body_mass_g2_weighted_coef + 
                     log10(body_mass_grams)^3 * coefs$log10_body_mass_g3_weighted_coef +
                     npp * coefs$npp_weighted_coef +
                     npp^2 * coefs$npp2_weighted_coef +
                     pcv * coefs$pcv_weighted_coef + 
                     pcv^2 * coefs$pcv2_weighted_coef +
                     pwarmest * coefs$pwarmest_weighted_coef +
                     pwarmest^2 * coefs$pwarmest2_weighted_coef +
                     richness * coefs$richness_mammals_weighted_coef)) %>%
      mutate(Density = squish(Density, range = c(hab_min, hab_max)))
    
    density_raster = raster(matrix(tmp.df$Density, byrow = TRUE, nrow = nrows))
    
    return(density_raster)
  }


pop_estimate_spatial = 
  function(density_raster,
           patch_no_migrate,
           patch_migrate,
           habitat_availability) {
    tmp.df = data.frame(density = getValues(density_raster),
                        patch_no_migrate = getValues(patch_no_migrate),
                        patch_migrate = getValues(patch_migrate),
                        habitat_availability_sqkm = getValues(habitat_availability)) %>%
      mutate(estimated_population = habitat_availability_sqkm * density * 2.25) %>%
      mutate(count = 1) %>%
      rowsum(., group = .$patch_no_migrate, na.rm=TRUE) %>%
      dplyr::select(.,patch_no_migrate, patch_migrate, estimated_population, count, habitat_availability_sqkm) %>%
      mutate(patch_no_migrate = patch_no_migrate / count,
             patch_migrate = patch_migrate / count) %>%
      filter(patch_no_migrate > 0) %>%
      mutate(habitat_availability_sqkm = habitat_availability_sqkm * 2.25) %>%
      dplyr::select(-count)
    
    return(tmp.df)
  }


# Dispersal distance functino
dispersal_distance_function <-
  function(species, species.frame) {
    # Identifying mammals/birds/amphibians
    taxon <- species.frame$taxon[species.frame$Species %in% species] 
    
    if(taxon %in% 'Amphibians') {
      # Matching family/species/genus/etc as much as possible
      if(species %in% amp_dispersal_distance_frame) {
        dispersal_distance = mean(amp_dispersal_distance_frame$Distance[paste0(amp_dispersal_distance_frame$Genus, amp_dispersal_distance_frame$Species, sep = "_") %in% species])
      } else if(species.frame$genus[species.frame$Species %in% species] %in% amp_dispersal_distance_frame$Genus) {
        dispersal_distance = mean(amp_dispersal_distance_frame$Distance[amp_dispersal_distance_frame$Genus %in% species.frame$genus[species.frame$Species %in% species]])
      } else if(species.frame$family[species.frame$Species %in% species] %in% amp_dispersal_distance_frame$Family) {
        dispersal_distance = mean(amp_dispersal_distance_frame$Distance[amp_dispersal_distance_frame$Family %in% species.frame$family[species.frame$Species %in% species]])
      } else {
        dispersal_distance = mean(amp_dispersal_distance_frame$Distance)
      }
    } else if(taxon %in% 'Mammals') {
      # Getting dispersal distance for the species
      # This is currently here to make sure the code is working
      # Will ultimately be taken from a database
      # If body mass exists, reclaculating here using coefficients and body mass used for consistency with density estimates
      dispersal_distance = 10^(-1.810227 + 0.6628845 * log10(species.frame$est_mass_kg[k] * 1000)) * 1000
      
      # If body mass does not exist, trying to match with binomial, genus, family, and then order
      if(is.na(species.frame$est_mass_kg[k])) {
        dispersal_distance = mam_dispersal_distance_frame$estimated_median_dispersal_m[mam_dispersal_distance_frame$binomial %in% species.frame$Species[species.frame$Species %in% species]]
      }
      # If statement to catch pspecies not in the data frame
      if(max(dispersal_distance) < 0) {
        dispersal_distance = mean(mam_dispersal_distance_frame$estimated_median_dispersal_m[mam_dispersal_distance_frame$genus %in% species.frame$genus[species.frame$Species %in% species]], na.rm=TRUE)
      }
      if(max(dispersal_distance)<0) {
        dispersal_distance = mean(mam_dispersal_distance_frame$estimated_median_dispersal_m[mam_dispersal_distance_frame$genus %in% species.frame$genus[species.frame$Species %in% species]], na.rm=TRUE)
      }
      if(max(dispersal_distance)<0) {
        dispersal_distance = mean(mam_dispersal_distance_frame$estimated_median_dispersal_m[mam_dispersal_distance_frame$family %in% species.frame$family[species.frame$Species %in% species]], na.rm=TRUE)
      }
      if(max(dispersal_distance < 0)) {
        dispersal_distance = mean(mam_dispersal_distance_frame$estimated_median_dispersal_m[mam_dispersal_distance_frame$order %in% species.frame$order[species.frame$Species %in% species]], na.rm=TRUE)
      }
      if(max(dispersal_distance < 0)) {
        # if nothing else work,s then matching with 90th percentile dispersal distance across all mammals
        dispersal_distance = 4803
      }
    } else if (taxon %in% 'Birds') {

	    # Based on taxonomies
	    # And filling gaps, where needed
	    dispersal_distance <-
		    read.csv(paste0(getwd(),'/Other Data Inputs/Patch_Dispersal_Distances/bird_dispersal_distances.csv'), stringsAsFactors = FALSE) %>%
		    mutate(Scientific_Name = gsub(' ','_', Scientific_Name)) %>%
		    filter(Scientific_Name %in% species.frame$Species[k]) %>%
		    dplyr::select(dispersal_distance) %>%
		    as.numeric()

	    if(max(dispersal_distance < 0)) {
		    # using 90th percentile to fill gap if species not matched
		    dispersal_distance <- 
			    quantile(read.csv(paste0(getwd(),'/Other Data Inputs/Patch_Dispersal_Distances/bird_dispersal_distances.csv'), stringsAsFactors = FALSE)$dispersal_distance,.9,na.rm=TRUE) 
	    }
    } else if (taxon %in% 'Reptiles') {

	    # data from: Southwood and Avens, 2010, Journal of Comparative Physiology, 'Physiological, behavioural, and ecological aspects of migration in reptiles'
	    if(grepl('Thamnophis',species.frame$Species[k])) {
	    dispersal_distance <- 17000
	    } else if(grepl('Liasis',species.frame$Species[k])) {
	    dispersal_distance <- 12000
	    } else if(grepl('Iguana',species.frame$Species[k])) {
            dispersal_distance <- 3000
            } else if(grepl('Conolophus',species.frame$Species[k])) {
            dispersal_distance <- 10000
            } else if(grepl('Crocodylus',species.frame$Species[k])) {
            dispersal_distance <- 10000
            } else {dispersal_distance <- 0}
    }# end reptiles

    # returning dispersal distance
    return(dispersal_distance)
  } 









density.function.new <-
  function(esh.list,density.list) {
    
    climate_scens = c('esh','rcp_current','rcp_45','rcp_85') # list of climate scenarios
    out.df <- data.frame() # Data frame used to save density outputs
    for(d in 1:length(names(esh.list))) { # Looping through habitat maps
      tmp.df <- # Extracting habitat range and density estimates
        data.frame(esh = getValues(esh.list[[d]]),
                   density = getValues(density.list[[d]]))
      tmp.df <- tmp.df[tmp.df$esh %in% 1,] #Limiting only to species' habitat
      tmp.df <- tmp.df[!is.na(tmp.df$density),] # And limiting only where there are density estimates
      
      if(nrow(tmp.df) > 0) {
        out.df <- # Can only calculate these if there's remaining habitat...
          rbind(out.df,
                data.frame(climate_scen = climate_scens[d],
                           mean_density = mean(tmp.df$density),
                           median_density = median(tmp.df$density),
                           sd_density = sd(tmp.df$density),
                           density_5th = quantile(tmp.df$density,.05),
                           density_95th = quantile(tmp.df$density,.95)))
      } else {
        out.df <- # And an exception if there is no remaining habitat...
          rbind(out.df,
                data.frame(climate_scen = climate_scens[d],
                           mean_density = 'No habitat remaining',
                           median_density = 'No habitat remaining',
                           sd_density = 'No habitat remaining',
                           density_5th = 'No habitat remaining',
                           density_95th = 'No habitat remaining'))
      }
      
      
    }
    
    return(out.df)
  }



# Manages data indicating each species' rarity
species.traits.function <-
  function(i) {
    
    # Ag intensity coefficients
    intensity.coefs <- read.csv(paste0(getwd(),"/Ag Intensity Outputs/Response to Habitat Intensification/Coefficients_Rarity_Models_13Nov2023.csv"))
    
    # Importing species traits
    species.traits <-
      read.csv(paste0(getwd(),'/Ag Intensity Outputs/Response to Habitat Intensification/hab_pref_rarity_19Sep2023.csv')) %>%
      mutate(Class = ifelse(taxon %in% 'Amphibians','Amphibia',
                            ifelse(taxon %in% 'Birds','Aves',
                                   ifelse(taxon %in% 'Mammals','Mammalia',
                                          ifelse(taxon %in% 'Reptiles','Reptilia', NA))))) %>%
      mutate(Hab_Rarity = ifelse(is.na(Hab_Rarity),'common',Hab_Rarity)) %>%
      mutate(Area_Rarity = ifelse(is.na(Area_Rarity),'common',Area_Rarity)) %>%
      mutate(Density_Rarity = ifelse(is.na(Density_Rarity),'common',Density_Rarity)) %>%
      mutate(Habitat_Breadth = str_to_title(Hab_Rarity),
             Habitat_Range = str_to_title(Area_Rarity)) %>%
      dplyr::select(Class, binomial = species, Habitat_Breadth, Habitat_Range) 
    
    # Mering traits with response estimates
    species.traits <-
      left_join(species.traits,
                intensity.coefs %>%
                  dplyr::select(Class,
                                Habitat_Range, Habitat_Breadth, estimate_type,
                                Urban = LandUseUrban,
                                Low_Cropland = LandUseLow.Intensity.Cropland,
                                Mod_High_Cropland = LandUseModerate_or_High.Intensity.Cropland,
                                Pasture = LandUsePasture))
    
    # Updating estimates so everything is in response to 0
    # E.g. 0 = no effect of land use on species
    # > 0 = positive effect of land use on species
    # < 0 = negative effect on species
    species.traits <-
      species.traits %>%
      mutate(Urban = Urban - 1,
             Low_Cropland = Low_Cropland - 1,
             Mod_High_Cropland = Mod_High_Cropland - 1,
             Pasture = Pasture - 1)
      
    
    # Returning data frame
    return(species.traits)
  } # End species traits function


species.traits.function_old <-
  function(i) {

	  # Ag intensity coefficients
    intensity.coefs <- read.csv(paste0(getwd(),"/Ag Intensity Outputs/Response to Habitat Intensification/Coefficients_Rarity_Models_13Nov2023.csv"))
    
    
    # Importing species traits
    species.traits <-
      read.csv(paste0(getwd(),"/Ag Intensity Outputs/Response to Habitat Intensification/hab_pref_rarity_19Sep2023.csv"), stringsAsFactors = FALSE) %>%
      dplyr::select(taxon, binomial,Density_Common, Niche_Breadth_Common, Range_Size_Common) %>% # Limiting to only needed columns
      mutate(Rarity_Abundance = ifelse(is.na(Density_Common),'Common',Density_Common)) %>% # Filling in gaps
      # Assumes species are common if we don't have data for them. 
      # This avoids overestimating future losses
      # Common species typically more tolerant of human activities than are less common species
      mutate(Rarity_Habitat = ifelse(is.na(Niche_Breadth_Common),'Common',Niche_Breadth_Common)) %>% # Filling in gaps
      mutate(Rarity_Range = ifelse(is.na(Range_Size_Common),'Common',Range_Size_Common)) %>% # Filling in gaps
      mutate(abund.quants = ifelse(Density_Common %in% 'Common',1,2), # Filling in gaps - assuming 
             habitat.quants = ifelse(Niche_Breadth_Common %in% 'Common',1,2), # Filling in gaps
             range.quants = ifelse(Range_Size_Common %in% 'Common',1,2)) %>% # Filling in gaps
      mutate(combined.quants = paste(paste(abund.quants, habitat.quants, range.quants))) %>%
      left_join(.,intensity.coefs) %>%
      dplyr::select(taxon, binomial,Rarity_Abundance,Rarity_Habitat,Rarity_Range,combined.quants,
                    LandUse,Estimate,Lower,Upper) %>%
      unique(.)
    
    # Mean estimate
    species.estimate <- 
      species.traits %>%
      dplyr::select(taxon, binomial, Rarity_Habitat, Rarity_Range, Rarity_Abundance, combined.quants, LandUse, Estimate) %>%
      tidyr::spread(.,LandUse,Estimate)
    
    # Lower estimate, if needed for sensitivity
    species.lower <-
      species.traits %>%
      dplyr::select(taxon, binomial, Rarity_Habitat, Rarity_Range, Rarity_Abundance, combined.quants, LandUse, Lower) %>%
      tidyr::spread(.,LandUse, Lower) %>%
      dplyr::select(taxon, binomial, 
                    lower_Amphibia = Amphibia, lower_Aves = Aves, lower_Mammalia = Mammalia,
                    `lower_Minimal Cropland` = `Minimal Cropland`, `lower_Mod. or Int. Cropland` = `Mod. or Int. Cropland`, lower_Pasture = Pasture)
    
    # Upper estimate, if needed for sensitivity
    species.upper <-
      species.traits %>%
      dplyr::select(taxon, binomial, Rarity_Habitat, Rarity_Range, Rarity_Abundance, combined.quants, LandUse, Upper) %>%
      tidyr::spread(.,LandUse, Upper) %>%
      dplyr::select(taxon, binomial, 
                    upper_Amphibia = Amphibia, upper_Aves = Aves, upper_Mammalia = Mammalia,
                    `upper_Minimal Cropland` = `Minimal Cropland`, `upper_Mod. or Int. Cropland` = `Mod. or Int. Cropland`, upper_Pasture = Pasture)
    
    # Manging output data frame
    species.traits <-
      species.estimate %>%
      mutate(Class = ifelse(taxon %in% 'Aves',Aves + Amphibia, 
                            ifelse(taxon %in% 'Mammalia',Mammalia + Amphibia,Amphibia))) %>%
      mutate(`Minimal Cropland` = `Minimal Cropland` + Class) %>%
      mutate(`Mod. or Int. Cropland` = `Mod. or Int. Cropland` + Class) %>%
      mutate(Pasture = Pasture + Class) %>%
      mutate(Class = 10^(Class + .001)) %>%
      mutate(`Minimal Cropland` = 10^(`Minimal Cropland` + .001)) %>%
      mutate(`Mod. or Int. Cropland` = 10^(`Mod. or Int. Cropland` + .001)) %>%
      mutate(Pasture = 10^(Pasture + .001)) %>%
      mutate(`Minimal Cropland` = `Minimal Cropland` / Class,
             `Mod. or Int. Cropland` = `Mod. or Int. Cropland` / Class,
             Pasture = Pasture / Class,
             Class = Class / Class)
    
    # Returning data frame
    return(species.traits)
  } # End species traits function



# Wrapper to get species' IUCN habitat preferences
hab.prefs.function <-
  function(i) {
    # Number habitat preferences
    amp.pref <- read.csv(paste0(getwd(),"/Other Data Inputs/Habitat_Preferences/AmphibianHabitats2018.csv"), stringsAsFactors = FALSE) %>% # Data
      filter(Suitability %in% c('Suitable','Marginal')) %>%
      filter(Habitat %in% c('14.1. Artificial/Terrestrial - Arable Land','14.3. Artificial/Terrestrial - Plantations','15.8. Artificial/Aquatic - Seasonally Flooded Agricultural Land',
                            '14.2. Artificial/Terrestrial - Pastureland', "14.5. Artificial/Terrestrial - Urban Areas")) %>%
      mutate(Habitat = ifelse(Habitat == '14.2. Artificial/Terrestrial - Pastureland','Pastureland',
                              ifelse(Habitat %in% "14.5. Artificial/Terrestrial - Urban Areas",'Urban','Cropland'))) %>%
      mutate(Pasture_Tolerance = ifelse(Habitat == 'Pastureland',1,0)) %>%
      mutate(Crop_Tolerance = ifelse(Habitat %in% 'Cropland',1,0)) %>%
      mutate(Urban_Tolerance = ifelse(Habitat %in% 'Urban' & Suitability %in% 'Marginal',.5,
                                      ifelse(Habitat %in% 'Urban' & Suitability %in% 'Suitable',1,0))) %>%
      mutate(Taxon = 'Amphibia') %>%
      dplyr::select(Taxon, Species_ID, species, Habitat, Suitability, Pasture_Tolerance, Crop_Tolerance, Urban_Tolerance)
    
    
    
    # Repeating for birds
    bird.pref <- read.csv(paste0(getwd(),"/Other Data Inputs/Habitat_Preferences/BirdHabitats2017.csv"), stringsAsFactors = FALSE) %>% # Data
      mutate(Habitat = Habitats.Classification.Scheme.Level.2) %>%
      filter(Suitability %in% c('Suitable','Marginal')) %>%
      filter(Habitat %in% c('Artificial/Terrestrial - Arable Land','Artificial/Terrestrial - Plantations','Artificial/Aquatic - Seasonally Flooded Agricultural Land',
                            'Artificial/Terrestrial - Pastureland', "Artificial/Terrestrial - Urban Areas")) %>%
      mutate(Habitat = ifelse(Habitat == 'Artificial/Terrestrial - Pastureland','Pastureland',
                              ifelse(Habitat %in% "Artificial/Terrestrial - Urban Areas",'Urban','Cropland'))) %>%
      mutate(Pasture_Tolerance = ifelse(Habitat == 'Pastureland',1,0)) %>%
      mutate(Crop_Tolerance = ifelse(Habitat %in% 'Cropland',1,0)) %>%
      mutate(Urban_Tolerance = ifelse(Habitat %in% 'Urban' & Suitability %in% 'Marginal',.5,
                                      ifelse(Habitat %in% 'Urban' & Suitability %in% 'Suitable',1,0))) %>%
      mutate(Taxon = 'Aves') %>%
      dplyr::select(Taxon, Species_ID = SIS.ID, species = Scientific.name, Habitat, Suitability, Pasture_Tolerance, Crop_Tolerance, Urban_Tolerance)
    
    
    
    # Repeating for mammals
    mam.pref <- read.csv(paste0(getwd(),"/Other Data Inputs/Habitat_Preferences/MammalHabitats2017.csv"), stringsAsFactors = FALSE) %>% # Data
      mutate(Habitat = Habitat.description) %>%
      filter(Suitability %in% c('Suitable','Marginal')) %>%
      filter(Habitat %in% c('Artificial/Terrestrial - Arable Land','Artificial/Terrestrial - Plantations','Artificial/Aquatic - Seasonally Flooded Agricultural Land',
                            'Artificial/Terrestrial - Pastureland', "Artificial/Terrestrial - Urban Areas")) %>%
      mutate(Habitat = ifelse(Habitat == 'Artificial/Terrestrial - Pastureland','Pastureland',
                              ifelse(Habitat %in% "Artificial/Terrestrial - Urban Areas",'Urban','Cropland'))) %>%
      mutate(Pasture_Tolerance = ifelse(Habitat == 'Pastureland',1,0)) %>%
      mutate(Crop_Tolerance = ifelse(Habitat %in% 'Cropland',1,0)) %>%
      mutate(Urban_Tolerance = ifelse(Habitat %in% 'Urban' & Suitability %in% 'Marginal',.5,
                                      ifelse(Habitat %in% 'Urban' & Suitability %in% 'Suitable',1,0))) %>%
      mutate(Taxon = 'Mammalia') %>%
      dplyr::select(Taxon,Species_ID = assessment.id, species = Species, Habitat, Suitability, Pasture_Tolerance, Crop_Tolerance, Urban_Tolerance)
   
    # Repeating for amphibians
    mam.pref <- read.csv(paste0(getwd(),"/Other Data Inputs/Habitat_Preferences/MammalHabitats2017.csv"), stringsAsFactors = FALSE) %>% # Data
      mutate(Habitat = Habitat.description) %>%
      filter(Suitability %in% c('Suitable','Marginal')) %>%
      filter(Habitat %in% c('Artificial/Terrestrial - Arable Land','Artificial/Terrestrial - Plantations','Artificial/Aquatic - Seasonally Flooded Agricultural Land',
                            'Artificial/Terrestrial - Pastureland', "Artificial/Terrestrial - Urban Areas")) %>%
      mutate(Habitat = ifelse(Habitat == 'Artificial/Terrestrial - Pastureland','Pastureland',
                              ifelse(Habitat %in% "Artificial/Terrestrial - Urban Areas",'Urban','Cropland'))) %>%
      mutate(Pasture_Tolerance = ifelse(Habitat == 'Pastureland',1,0)) %>%
      mutate(Crop_Tolerance = ifelse(Habitat %in% 'Cropland',1,0)) %>%
      mutate(Urban_Tolerance = ifelse(Habitat %in% 'Urban' & Suitability %in% 'Marginal',.5,
                                      ifelse(Habitat %in% 'Urban' & Suitability %in% 'Suitable',1,0))) %>%
      mutate(Taxon = 'Mammalia') %>%
      dplyr::select(Taxon,Species_ID = assessment.id, species = Species, Habitat, Suitability, Pasture_Tolerance, Crop_Tolerance, Urban_Tolerance)

# reptiles
rep.pref <- 
  read.csv(paste0(getwd(),"/Other Data Inputs/Habitat_Preferences/ReptileHabitats2023.csv"), stringsAsFactors = FALSE) %>% # data
  dplyr::rename(Habitat = name,
                Suitability = suitability,
                assessment.id = assessmentId,
                Species = scientificName) %>%
  filter(Suitability %in% c('Suitable','Marginal')) %>%
  filter(Habitat %in% c('Artificial/Terrestrial - Arable Land','Artificial/Terrestrial - Plantations','Artificial/Aquatic - Seasonally Flooded Agricultural Land',
                        'Artificial/Terrestrial - Pastureland', "Artificial/Terrestrial - Urban Areas")) %>%
  mutate(Habitat = ifelse(Habitat == 'Artificial/Terrestrial - Pastureland','Pastureland',
                          ifelse(Habitat %in% "Artificial/Terrestrial - Urban Areas",'Urban','Cropland'))) %>%
  mutate(Pasture_Tolerance = ifelse(Habitat == 'Pastureland',1,0)) %>%
  mutate(Crop_Tolerance = ifelse(Habitat %in% 'Cropland',1,0)) %>%
  mutate(Urban_Tolerance = ifelse(Habitat %in% 'Urban' & Suitability %in% 'Marginal',.5,
                                  ifelse(Habitat %in% 'Urban' & Suitability %in% 'Suitable',1,0))) %>%
  mutate(Taxon = 'Reptilia') %>%
  dplyr::select(Taxon,Species_ID = assessment.id, species = Species, Habitat, Suitability, Pasture_Tolerance, Crop_Tolerance, Urban_Tolerance)
    

    # And rbinding
    aoh.habs <-
      rbind(amp.pref,
            bird.pref,
            mam.pref,
	    rep.pref) %>%
      mutate(species = gsub(" ","_",species))
    
    # And returning
    return(aoh.habs)
    
  } # End hab prefs function





# Function to upload and manage population density coefficients
pop.dens.coef.function <-
  function(species.frame) {
    
    # Uploading coefficients for pop density estimates
    # Order
    coef.order <- 
      rbind(read.csv(paste0(getwd() %>% gsub('pubh-glob2loc','ouce-glob2loc',.),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_order_Mammals.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Mammals'),
            # read.csv(paste0(getwd(),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_order_Amphibians.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Amphibians'),
            read.csv(paste0(getwd() %>% gsub('pubh-glob2loc','ouce-glob2loc',.),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_order_Birds.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Birds'))
    # Family
    coef.family <- 
      rbind(read.csv(paste0(getwd() %>% gsub('pubh-glob2loc','ouce-glob2loc',.),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_family_Mammals.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Mammals'),
            read.csv(paste0(getwd() %>% gsub('pubh-glob2loc','ouce-glob2loc',.),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_family_Amphibians.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Amphibians'),
            read.csv(paste0(getwd() %>% gsub('pubh-glob2loc','ouce-glob2loc',.),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_family_Birds.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Birds'),
	    read.csv(paste0(getwd() %>% gsub('pubh-glob2loc','ouce-glob2loc',.),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_family_Reptiles.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Reptiles'))
    # Binomial
    coef.binomial <- 
      rbind(read.csv(paste0(getwd() %>% gsub('pubh-glob2loc','ouce-glob2loc',.),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_binomial_Mammals.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Mammals'),
            read.csv(paste0(getwd() %>% gsub('pubh-glob2loc','ouce-glob2loc',.),"/Population Density Estimates/RefittedSantiniModelCoefficients_random_binomial_Birds.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Birds'))
    # Climate
    coef.climate <- 
      rbind(read.csv(paste0(getwd() %>% gsub('pubh-glob2loc','ouce-glob2loc',.),"/Population Density Estimates/RefittedSantiniModelCoefficients_fixed_Mammals.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Mammals'),
            read.csv(paste0(getwd() %>% gsub('pubh-glob2loc','ouce-glob2loc',.),"/Population Density Estimates/RefittedSantiniModelCoefficients_fixed_Amphibians.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Amphibians'),
            read.csv(paste0(getwd() %>% gsub('pubh-glob2loc','ouce-glob2loc',.),"/Population Density Estimates/RefittedSantiniModelCoefficients_fixed_Birds.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Birds'),
	    read.csv(paste0(getwd() %>% gsub('pubh-glob2loc','ouce-glob2loc',.),"/Population Density Estimates/RefittedSantiniModelCoefficients_fixed_Reptiles.csv"), stringsAsFactors = FALSE) %>% mutate(taxon = 'Reptiles'))

    # Managing body mass data
    # And then merging into the data set

    # bird/amp/mammal body mass
body_mass <-
  read.csv(paste0(getwd() %>% gsub('pubh-glob2loc','ouce-glob2loc',.),"/Other Data Inputs/Pop Density Inputs/Body Mass Estimates 20February2020 Updated Taxonomy.csv"), stringsAsFactors = FALSE) %>% # Importing body mass file
  #read.csv(paste0(getwd(),"/Other Data Inputs/Pop Density Inputs/Body Mass Estimates 14November2019.csv"), stringsAsFactors = FALSE) %>% # Importing body mass file
  dplyr::select(order = Order, family = Family, genus = Genus, species_merge = binomial, est_mass_kg, Family_Mass_kg, Genus_Mass_kg) %>%  # Only keeping necessary columns
  mutate(binomial = species_merge)

# reptile body mass
rep_body_mass <-
  read.csv('/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/Other Data Inputs/Pop Density Inputs/reptile_body_masses.csv',
           stringsAsFactors = FALSE) %>%
  dplyr::rename(order = Order,
                family = Family) %>%
  mutate(est_mass_kg = `mass..g.` / 1000) %>%
  mutate(Family_Mass_kg = est_mass_kg,
         Genus_Mass_kg = est_mass_kg) %>%
  mutate(binomial = gsub(' ','_', `X...binomial`)) %>%
  mutate(genus = gsub('_.*','',binomial),
         species_merge = binomial) 
  
  # rbinding the body mass estimates for reps and non reps
  body_mass_merge <-
  rbind(body_mass,
        rep_body_mass[,names(body_mass)])
    
    # Merging body mass and family coefficients into the species data set ----
    species.frame <-
      left_join(species.frame %>%
		mutate(species_merge = gsub('.*/','',Species)),
	body_mass_merge) %>%
      left_join(.,coef.order) %>%
      mutate(family = str_to_title(family)) %>%
      left_join(., coef.family) %>%
      left_join(.,coef.binomial)


    # Converting NAs to 0s in the order, family, and binomial coefficients
    species.frame$order_intercept_adj_weighted[is.na(species.frame$order_intercept_adj_weighted)] <- 0
    species.frame$family_intercept_adj_weighted[is.na(species.frame$family_intercept_adj_weighted)] <- 0
    species.frame$binomial_intercept_adj_weighted[is.na(species.frame$binomial_intercept_adj_weighted)] <- 0
    
    # Returning data frame
    return(species.frame)
    
  } # End function






