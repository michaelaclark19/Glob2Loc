#!/usr/bin/env Rscript

#######
# Testing accuracy based on 2010 SDM outputs
# Doing this based on model thresholds of the different threshold types
# Is there one that is consistently best?
#######

###
# Libraries
library(plyr)
library(dplyr)
library(raster)
library(stringr)
library(parallel)

# functions for the SDMs
weight.raw.rasters <-
  function(file.list) {
    # Accuracy csv
    auc.csv <-
      read.csv(auc.thresh.list[grepl(tmp.name,auc.thresh.list)] %>%
                 .[grepl('_Accuracy_',.)])
    
    # Threshold csv
    thresh.csv <-
      read.csv(auc.thresh.list[grepl(tmp.name,auc.thresh.list)] %>%
                 .[grepl('_Thresholds_',.)])
    
    # Stacking models
    models <- stack()
    for(k in file.list) {
      models <- raster::stack(models, raster(k))
    }
    
    # Changing names of models
    names(models) <- gsub('.*[0-9]{4,4}_','',names(models))
    
    # Accuracy for weighting
    names(auc.csv)[names(auc.csv) %in% 'rf'] <- 'random.forest'
    auc.weight <- as.numeric(auc.csv[,names(models)])
    
    ###
    # Next getting model average based on auc values
    names(auc.weight) <- names(models)
    # Getting weighted model means
    # List of model aucs
    auc <- auc.weight
    # Weighting based on accuracy
    # AUC = 0.5 means model is as good as random
    # Weighting based on square of the distance fomr auc = 0.5
    # So more weighting on more accurate models
    w <- (auc-0.5)^2
    
    
    # Weighting the models
    # mods.weighted <- weighted.mean(models, w)
    # This is slow, the below returns the same estimates, with the exception of rounding errors on the order of 1e-16
    w <- w/sum(w)
    
    # Speeding up if only one SDM fit well
    if(length(names(models)) %in% 1) {
      mods.weighted <- models
    } else {
      mods.weighted <- weight.fun.raster(w, models)
    }
    
    # Getting weighted threshold
    tr.weight <- as.numeric(thresh.csv[thresh.csv$model_name %in% names(models),
                                       gsub('.*2010_threshold_','',i) %>% gsub('.tif','',.)])
    
    # Getting weighted thresholds
    tr.weight <- weighted.mean(tr.weight, w)
    
    # Returning weighted models and weighted threshold
    return(list(mods.weighted,tr.weight))
  }

weight.fun.raster <-
  function(weights,
           models) {
    tmp.mat <-
      matrix(ncol = length(weights),
             nrow = models@nrows * models@ncols, NA)
    
    for(m in 1:ncol(tmp.mat)) {
      tmp.mat[,m] <- getValues(models[[which(names(models) %in% names(weights)[m])]]) * weights[[m]]
    }
    
    out.raster <- raster(matrix(base::rowSums(tmp.mat),byrow=TRUE,nrow = models@nrows), crs = crs(models))
    extent(out.raster) <- extent(models)
    
    return(out.raster)
  }


# Installing functions
setwd("/data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity")
source(paste0(getwd(),'/Scripts/SDM Scripts/0.0_SDM_Functions_13Sept2022.R'))

# Ecoregion map
tnc.eco <- raster(paste0(getwd(),'/Ecoregions_Feb2023/biomes_realms_raster.tif'))

# Have multiple working directories because of memory issues
# Needed to get creative!
wd.list <- getwd()
  
# List of species we have
species.list <-
  list.dirs(paste0(getwd(),'/ESH_RCPs/Weighted_Threshold_Accuracy')) %>% 
  lapply(.,list.files,full.names=TRUE) %>%
  do.call(c,.) %>%
  gsub('\\.csv','',.) %>%
  gsub('.*Weighted_Threshold_Accuracy/','',.)

# Getting list of species with outputs in 2090
species.list.2050 <-
  do.call(c,
          lapply(paste0(wd.list,'/ESH_RCPs/Raw/') %>% paste0(.,sort(rep(c('Amphibians','Mammals','Birds','Reptiles'),2))),
                 list.files,
                 full.names = TRUE)) %>%
  .[grepl('SSP2.*204',.)] %>%
  gsub('__SSP2.*','',.) %>%
  gsub('_SSP2.*','',.) %>%
  unique()

# Shortening to check overlaps...
species.list.check <- species.list.2050 %>% gsub('.*Raw/','',.) %>% unique()

# And filtering to species with data in both time periods
species.list <- species.list[species.list %in% species.list.check]

# List of sdm outputs
sdms.list <-
  do.call(c,
          lapply(paste0(wd.list,'/ESH_RCPs/Raw/') %>% paste0(.,sort(rep(c('Amphibians','Mammals','Birds','Reptiles'),2))),
                 list.files,
                 full.names = TRUE))

# list of accuracy and threshold files
auc.list <-
  list.files(paste0(getwd(),'/ESH_RCPs/Mod_Accuracy'), full.names = TRUE)

# list of AOH rasters
aoh.list <-
  list.dirs(paste0(getwd(),'/ESH_Tifs_12Oct')) %>%
  list.files(.,full.names=TRUE)

###
# Removing species we already have...
species.have <-
  c(list.files(paste0(getwd(),'/ESH_RCPs/SSP5-8.5/Amphibians'), full.names = TRUE),
    list.files(paste0(getwd(),'/ESH_RCPs/SSP5-8.5/Birds'), full.names = TRUE),
    list.files(paste0(getwd(),'/ESH_RCPs/SSP5-8.5/Mammals'), full.names = TRUE),
    list.files(paste0(getwd(),'/ESH_RCPs/SSP5-8.5/Reptiles'), full.names = TRUE)) %>%
  .[grepl('.csv',.)] %>% gsub('_completed.*','',.) %>% gsub('.*SSP5-8.5/','',.)

# Overlap between species list and species already have

if(length(species.have) > 0) {
  # Filter species
  species.need <- which(!(gsub('.*Raw/','',species.list) %in% species.have))
  species.list <- species.list[species.need]
} else {
  # Don't filter species
  species.list <- species.list
}

# And removing species we already have

# List of years
years <- 
  sdms.list %>% 
  str_extract(.,'[0-9]{4,4}.[0-9]{4,4}') %>% unique() %>% .[!is.na(.)] %>%
  .[!grepl('1850|196|198|208',.)]

# List of SSPs
ssps <- list.files(paste0(getwd(),'/ESH_RCPs'), pattern = 'SSP') %>% str_extract(.,'SSP[0-9]')
ssps.write <- list.files(paste0(getwd(),'/ESH_RCPs'), pattern = 'SSP')


cat('\n', ssps[c(2)],'\n')

# Checking it made it this far...
head(species.list)

# Function that loops through species
species_wrap_fun <-
  function(i) {
    
    # Species taxa
    tmp.taxa <-
      gsub('.*/ESH_RCPs/Raw/','',i) %>%
      gsub('\\/.*','',.)
    
    # tmp.taxa <- 'Amphibians'
    
    # Species name
    tmp.name <-
      gsub(paste0('.*',tmp.taxa),'',i) %>%
      gsub('\\/','',.) %>%
      gsub('__His.*','',.)
    
    # 
    # tmp.names <- 'Brachycephalus_ephippium'
    
    # Getting esh file
    esh.raster <- raster(aoh.list[grep(paste0(tmp.taxa,'.*',tmp.name), aoh.list)])
    
    # Which realms_biomes are in the species esh?
    # This gets those unique values, and then clips the realm biome raster to only include those
    # This is used as a masking layer later on
    tmp.realm.biome <- raster::crop(tnc.eco, esh.raster)
    tmp.realm.values <- unique(tmp.realm.biome[esh.raster %in% 1]) %>% .[!is.na(.)] %>% .[!(.%in%1)]
    tmp.realm.biome[!(tmp.realm.biome %in% tmp.realm.values)] <- NA
    
    # Auc csv
    auc.csv <-
      tryCatch(read.csv(auc.list[grepl(tmp.name,auc.list)] %>% .[grepl('_Accuracy_',.)]),
               error = function(e) {return('Whoops')})
    
    # Threshold csv
    thresh.csv <-
      tryCatch(read.csv(auc.list[grepl(tmp.name,auc.list)] %>% .[grepl('_Thresholds_',.)]), 
               error = function(e) {return('Whoops')})
    
    # If statement making sure these exist
    if(is.data.frame(auc.csv) & is.data.frame(thresh.csv)) {
      
      # Looping through SSPs
      for(ss in ssps[c(2)]) {
        cat(ss,'\n')
        
        # Empty raster stack
        out.year.stack <- raster::stack()
        
        # Loop through years
        for(yy in years) {
          cat(ss,yy,'\n')
          
          # Getting models
          if(yy %in% years[1]) {
            model.files <-
              sdms.list[grepl(paste0(i,'.*',yy),sdms.list)]
          } else {
            model.files <-
              sdms.list[grepl(paste0(i,'.*',ss,'.*',yy),sdms.list)]
          }
          
          if(sum(grepl('__Historic|__SSP',model.files)) > 1 & sum(grepl('[a-z]_Historic|[a-z]_SSP',model.files)) > 1) {
            model.files <- model.files[grepl('[a-z]_Historic|[a-z]_SSP',model.files)]
          }
          cat(model.files)
          
          
          
          
          if(sum(grepl(getwd(),model.files)) > 0 & sum(grepl(getwd(),model.files)) > 0) { # Checking if there are files in both the ouce and pubh directories. If so, take the ouce files
            model.files <- model.files[grepl(getwd(),model.files)]
          } # End if statement
          
          # Stacking if needed
          models <- stack()
          for(k in model.files) {
            cat('Importing',k,'\n')
            raster(k)
            cat('End Import',k,'\n')
            models <- raster::stack(models, raster(k))
          }
          
          cat('End Raw SDMs Import')
          # Names of models
          names(models) <-
            gsub('.*[0-9]{4,4}_','',names(models))
          
          # Checking
          if(length(grep(tmp.name,auc.list)) %in% 0) {
            # Do nothing
          } else {
            cat('Importing Thresholds and Accuracy \n')
            # Getting thresholds and accuracy
            # Accuracy
            auc.csv <-
              read.csv(auc.list[grepl(tmp.name,auc.list)] %>%
                         .[grepl('_Accuracy_',.)])
            # And another check
            if(length(auc.list[grepl(tmp.name,auc.list)] %>%
                      .[grepl('_Thresholds_',.)]) > 0) {
              
              # Thresholds
              thresh.csv <-
                read.csv(auc.list[grepl(tmp.name,auc.list)] %>%
                           .[grepl('_Thresholds_',.)])
              
              # Accuracy for weighting
              names(auc.csv)[names(auc.csv) %in% 'rf'] <- 'random.forest'
              auc.weight <- as.numeric(auc.csv[,names(models)])
              
              ###
              # Next getting model average based on auc values
              names(auc.weight) <- names(models)
              # Getting weighted model means
              # List of model aucs
              auc <- auc.weight
              # Weighting based on accuracy
              # AUC = 0.5 means model is as good as random
              # Weighting based on square of the distance fomr auc = 0.5
              # So more weighting on more accurate models
              w <- (auc-0.5)^2
              
              
              # Weighting the models
              # mods.weighted <- weighted.mean(models, w)
              # This is slow, the below returns the same estimates, with the exception of rounding errors on the order of 1e-16
              w <- w/sum(w)
              
              # Getting thresholds - both prevalence and equal spec sens
              # These are the two approaches normally used for weighting SDMs
              # Prevalence
              thresh.prev <- 
                sum(as.numeric(thresh.csv$prevalence[match(names(w),thresh.csv$model_name)]) * # Ordering thresholds
                      w) # Multiplying by weight
              # Spec sens
              thresh.specsens <- 
                sum(as.numeric(thresh.csv$equal_sens_spec[match(names(w),thresh.csv$model_name)]) * # Ordering thresholds
                      w) # Multiplying by weight
              
              cat('Weighting Rasters \n') 
              # Speeding up if only one SDM fit well
              if(length(names(models)) %in% 1) {
                mods.weighted <- models
              } else {
                mods.weighted <- weight.fun.raster(w, models)
              }
              
              # Adding to raster stack
              out.year.stack <- stack(out.year.stack, mods.weighted)
              
            } # End check for thresholds
          } # End else statement checking whether accuracy data for species SDM exists
        } # End of year loop
        
        
        cat('End Stacking Rasters','\n')
        
        # Updating names
        names(out.year.stack) <- years
        
        # Empty raster stack
        out.save.stack.prev <- raster::stack()
        out.save.stack.specsens <- raster::stack()
        
        # Looping through years to interpolate
        for(yy in 1:(length(years) - 1)) {
          cat(ss,yy,'\n')
          years.interp <- 
            seq(from = round_any(mean(as.numeric(unlist(str_extract_all(years[yy],'[0-9]{4,4}')))),5),
                to = round_any(mean(as.numeric(unlist(str_extract_all(years[yy+1],'[0-9]{4,4}')))),5),
                by = 5)
          
          for(y in seq_along(years.interp)) {
            
            # Interp factors
            interp.factors <- seq(from = 1, to = 0, length = length(years.interp))
            if(y %in% 1) {
              out.raster <- out.year.stack[[yy]]
            } else if (y %in% length(years.interp)) {
              out.raster <- out.year.stack[[yy + 1]]
            } else {
              out.raster <-
                out.year.stack[[yy]] * interp.factors[y] +
                out.year.stack[[yy+1]] * (1 - interp.factors[y])
            }
            # Species can only be in the realm/biome combination in which it is currently found
            out.raster[is.na(tmp.realm.biome)] <- NA
            
            # Templates for the two types of thresholds
            out.raster.prev <- out.raster
            out.raster.specsens <- out.raster
            
            # Updating based on thresholds
            out.raster.prev[out.raster.prev < thresh.prev] <- NA
            out.raster.prev[!is.na(out.raster.prev)] <- 1
            
            out.raster.specsens[out.raster.specsens < thresh.specsens] <- NA
            out.raster.specsens[!is.na(out.raster.specsens)] <- 1
            
            # removing water
            
            # Updating names
            names(out.raster.prev) <- years.interp[y]
            names(out.raster.specsens) <- years.interp[y]
            # And stacking, with a logic check
            if(names(out.raster.prev) %in% names(out.save.stack.prev)) {
              # Do nothing
              # Write raster to save space
              writeRaster(out.raster.prev, 
                          paste0(getwd(),'/ESH_RCPs/',ssps.write[grep(ss,ssps.write)],'/',tmp.taxa,'/',tmp.name,'_',ss,'_','prevalence_threshold__',years.interp[y],'.tif'),
                          # paste0('/data/pubh-glob2loc/pubh0329/trial_transfer_21Sept2022/',tmp.name,'_',ss,'_','prevalence_threshold'),
                          # suffix = gsub('X','',names(out.save.stack.prev)),
                          format = 'GTiff',
                          # bylayer = TRUE,
                          overwrite = TRUE)
              
              writeRaster(out.raster.specsens, 
                          paste0(getwd(),'/ESH_RCPs/',ssps.write[grep(ss,ssps.write)],'/',tmp.taxa,'/',tmp.name,'_',ss,'_','specsens_threshold__',years.interp[y],'.tif'),
                          # paste0('/data/pubh-glob2loc/pubh0329/trial_transfer_21Sept2022/',tmp.name,'_',ss,'_','prevalence_threshold'),
                          # suffix = gsub('X','',names(out.save.stack.prev)),
                          format = 'GTiff',
                          # bylayer = TRUE,
                          overwrite = TRUE)
              
            } else { # Stack
              # out.save.stack.prev <- raster::stack(out.save.stack.prev, out.raster.prev)
              # out.save.stack.specsens <- raster::stack(out.save.stack.specsens, out.raster.specsens)
              writeRaster(out.raster.prev, 
                          paste0(getwd(),'/ESH_RCPs/',ssps.write[grep(ss,ssps.write)],'/',tmp.taxa,'/',tmp.name,'_',ss,'_','prevalence_threshold__',years.interp[y],'.tif'),
                          # paste0('/data/pubh-glob2loc/pubh0329/trial_transfer_21Sept2022/',tmp.name,'_',ss,'_','prevalence_threshold'),
                          # suffix = gsub('X','',names(out.save.stack.prev)),
                          format = 'GTiff',
                          # bylayer = TRUE,
                          overwrite = TRUE)
              
              writeRaster(out.raster.specsens, 
                          paste0(getwd(),'/ESH_RCPs/',ssps.write[grep(ss,ssps.write)],'/',tmp.taxa,'/',tmp.name,'_',ss,'_','specsens_threshold__',years.interp[y],'.tif'),
                          # paste0('/data/pubh-glob2loc/pubh0329/trial_transfer_21Sept2022/',tmp.name,'_',ss,'_','prevalence_threshold'),
                          # suffix = gsub('X','',names(out.save.stack.prev)),
                          format = 'GTiff',
                          # bylayer = TRUE,
                          overwrite = TRUE)
            }
          } # End loop through years to interpolate (these are the 5 year periods)
          # names(out.save.stack) <- years.interp
        } # End loop through year periods to interpolate (1995, 2030, etc)
        
        # # And saving raster stack
        # writeRaster(out.save.stack.prev, 
        #             paste0(getwd(),'/ESH_RCPs/',ssps.write[grep(ss,ssps.write)],'/',tmp.taxa,'/',tmp.name,'_',ss,'_','prevalence_threshold_'),
        #             # paste0('/data/pubh-glob2loc/pubh0329/trial_transfer_21Sept2022/',tmp.name,'_',ss,'_','prevalence_threshold'),
        #             suffix = gsub('X','',names(out.save.stack.prev)),
        #             format = 'GTiff',
        #             bylayer = TRUE,
        #             overwrite = TRUE)
        # 
        # writeRaster(out.save.stack.specsens, 
        #             paste0(getwd(),'/ESH_RCPs/',ssps.write[grep(ss,ssps.write)],'/',tmp.taxa,'/',tmp.name,'_',ss,'_','specsens_threshold_'),
        #             # paste0('/data/pubh-glob2loc/pubh0329/trial_transfer_21Sept2022/',tmp.name,'_',ss,'_','equalspecsens_threshold'),
        #             suffix = gsub('X','',names(out.save.stack.specsens)),
        #             format = 'GTiff',
        #             bylayer = TRUE,
        #             overwrite = TRUE)
        
        gc()
        rm(out.save.stack.specsens)
        rm(out.save.stack.prev)
        rm(out.raster.prev)
        rm(out.raster.specsens)
        rm(out.year.stack)
        rm(thresh.csv)
        rm(auc.csv)
      } # End loop through SSPS
    } # End if statement checking for accuracy and thresholds
    
    
    # And pinging a csv file to save progress - using these to check progress ensures rasters have been properly saved...
    write.csv(data.frame(complete = 'complete'),
              paste0(getwd(),'/ESH_RCPs/SSP5-8.5/',tmp.taxa,'/',tmp.name,'_completed.csv'))  
    
    
    # Updating based on threshold
    # out.save.stack[out.save.stack < thresh] <- NA
    # out.save.stack[!is.na(out.save.stack)] <- 1
  } # End function

species.have <- species.have[!grepl('.csv',species.have)] %>% unique(.)

cat('Amphibians: ',sum(grepl('Amphibian',species.have)),'Species Completed\n')
cat('Birds: ',sum(grepl('Bird',species.have)),'Species Completed\n')
cat('Mammals: ',sum(grepl('Mammal',species.have)),'Species Completed\n')
cat('Reptiles: ',sum(grepl('Reptile',species.have)),'Species Completed\n')

cat('Amphibians: ',sum(grepl('Amphibian',species.list)),'Species Left\n')
cat('Birds: ',sum(grepl('Bird',species.list)),'Species Left\n')
cat('Mammals: ',sum(grepl('Mammal',species.list)),'Species Left\n')
cat('Reptiles: ',sum(grepl('Reptile',species.list)),'Species Left\n')


# lapply(species.list,species_wrap_fun)
# Getting location to split species list
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
# species.list.amps <- species.list[grepl('Mammals',species.list)]

species.list <- chunk2(species.list,5)[[1]]
#species_wrap_fun(species.list[[1]][1])
# and running

# Adding a try catch statement around the full function
# This catches species with errors
# But ensures that script will keep running
# And documents species where the error doesn't work

# And building trycatch around the accuracy checking function
try_catch_function <-
  function(k) {
    tryCatch(species_wrap_fun(k),
             #k*5,
             error = function(e) { # Try catch error function
               tmp_species <- gsub('.*ESH_RCPs/Raw','',k)
               
               # which species did it not work for?
               cat('Whoopsies:',tmp_species, '\n')
               # tracking species it did not work for
               
               write.csv(data.frame(complete = 'complete'),
                         paste0(getwd(),'/ESH_RCPs/SSP5-8.5/',tmp_species, '_error_completed.csv'))
             }) # End try catch 
  }

length(species.list)
#lapply(species.list,try_catch_function)
mclapply(species.list, try_catch_function, mc.cores = 5)

# cat(species.list)
# lapply(rev(species.list),species_wrap_fun)
# 
# 
# plan(multisession, workers = 3)
# # future_lapply(rep(5,5),Sys.sleep)
# future_lapply((species.list),species_wrap_fun)
# plan(sequential)



