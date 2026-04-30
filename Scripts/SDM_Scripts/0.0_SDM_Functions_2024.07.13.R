###
# Writing functions for SDM modelling
###

###
# First function runs models
sdm.mod.fun <-
  function(raster.stack,
           raster.stack.bc.xm,
           df.dat,
           df.dat.bc.xm,
           esh.raster,
           realm,
           biome) {
    #####
    # Getting data ready to run bioclim, domain, and maxent
    # Extracting coordinates of cells that contain the species
    tmp.coords <-
      xyFromCell(esh.raster,
                 df.dat$cell_num)
    
    df.dat <-
      cbind(df.dat,
            as.data.frame(tmp.coords))
    
    df.dat.bc.xm <-
      cbind(df.dat.bc.xm,
            as.data.frame(tmp.coords))
    
    # Limiting data frames to only the ecoregions and land cover types
    # The species currently exists in
    eco.dat <- data.frame(biome = getValues(biome), realm = getValues(realm), esh.values = getValues(esh.raster))
    # Getting summary by ecoregion
    biome.sum.dat <- eco.dat %>% filter(esh.values %in% 1) %>% dplyr::select(-realm) %>% as.matrix() %>% rowsum(.,group = .[,1], na.rm = TRUE) %>%
      as.data.frame() %>% mutate(biome = biome/esh.values) %>% mutate(esh.values = esh.values/sum(esh.values,na.rm=TRUE))
    # Getting summary by biome
    realm.sum.dat <- eco.dat %>% filter(esh.values %in% 1) %>% dplyr::select(-biome) %>% as.matrix() %>% rowsum(.,group = .[,1], na.rm = TRUE) %>%
      as.data.frame() %>% mutate(realm = realm/esh.values) %>% mutate(esh.values = esh.values/sum(esh.values,na.rm=TRUE))
    
    # Sense checks to make sure
    # (a) there is some overlap between esh raster and biome and realm rasters
    # (b) at least 25% of species habitat range is on land
    
    
    # This first makes sure there is some overlap between the rasters
    if(!(max(unique(biome.sum.dat$biome)) %in% 0 | max(unique(realm.sum.dat$realm)) %in% 0)) {
      # This second makes sure that at least 25% of the habitat range is on land
      if(sum(biome.sum.dat$esh.values[biome.sum.dat$biome != 0]) >= .75 &
         sum(realm.sum.dat$esh.values[realm.sum.dat$realm != 0]) >= .75) {
        ### Filtering data frames by 
        # Land cover type
        esh.df <- df.dat[eco.dat$biome %in% biome.sum.dat$biome &
                           eco.dat$realm %in% realm.sum.dat$realm,]
        esh.df.bc.xm <- df.dat.bc.xm[eco.dat$biome %in% biome.sum.dat$biome &
                                       eco.dat$realm %in% realm.sum.dat$realm,]
        
        ###
        # Limiting only to cells that contain the species
        esh.df <-
          esh.df[!is.na(esh.df$esh),]
        
        esh.df.bc.xm <- 
          esh.df.bc.xm[!is.na(esh.df.bc.xm$esh),]
        
        ###
        # And now running models
        # This chunk is for setting up data
        
        ###
        # Bioclim model
        # Bioclim uses only presence values
        # pred_nf <- dropLayer(predictors,'biome')
        # group <- kfold(bradypus, 5)
        # pres_train <- bradypus[group != 1, ]
        # pres_test <- bradypus[group == 1, ]
        # bc <- bioclim(pred_nf, pres_train
        
        # Extracting values in these cells
        # tmp.test <-
        #   extract(esh.raster,
        #           tmp.coords,
        #           method = 'simple')
        
        # Dividng into groups for cross validation
        # Making another exception!
        # Need to make sure there are more than 50 data points in the presence data
        # Doing this here
        # Because this is after clipping land cover and ecoregion
        
        # Creating empty object
        # Doing this for the boolean later on
        backg <- c()
        
        if(nrow(esh.df)>50 & !is.na(minValue(raster.stack$GDD_2010))) {
          esh.df$group <-
            kfold(esh.df,
                  k = 5)
          if(nrow(esh.df) == nrow(esh.df.bc.xm)) {
            esh.df.bc.xm$group <-
              esh.df$group
          } else {
            esh.df.bc.xm$group <-
              kfold(esh.df.bc.xm,
                    k = 5)
          }
          
          
          ####
          # Splitting into test and training groups
          # Test group is 80% of the data set
          # Training group is 20% of the data set
          pres_train <-
            esh.df[esh.df$group != 1,]
          pres_test <-
            esh.df[esh.df$group == 1,]
          
          pres_train.rf.glm <-
            esh.df.bc.xm[esh.df.bc.xm$group != 1,]
          pres_test.rf.glm <-
            esh.df.bc.xm[esh.df.bc.xm$group == 1,]
          
          ###
          # Repeating for background data
          # For maxent, this is 10,000 background points
          # First just randomly selecting all the points
          # Note that this randomizes the order
          # This might also throw a message saying it did not extract as many points as wanted
          # This is ok
          # backg <- randomPoints(raster.stack, (esh.raster@nrows*esh.raster@ncols - sum(!is.na(getValues(esh.raster)))),tryf=10)
          backg <- randomPoints(raster.stack, min(round((esh.raster@nrows*esh.raster@ncols - sum(!is.na(getValues(esh.raster))))/4),(nrow(pres_train) + nrow(pres_test))*4),tryf=10)
          
          backg <- cbind(backg,1:nrow(backg))
          colnames(backg) = c('lon','lat','cell_num')
          
          # Limiting to data points where esh = 0 
          # Getting ESH values at the backg points
          # Ideally these will all be NAs
          esh.vals.backg <- raster::extract(x = esh.raster,
                                            y = backg[,c('lon','lat')],
                                            method = 'simple')
          # Dropping coordinates that don't contain NAs in the esh raster
          backg <- backg[which(is.na(esh.vals.backg)),]
          
          # Getting ecoregion and biome of the selected points
          back.land <- data.frame(biome.tmp = raster::extract(biome,y = backg[,c('lon','lat')], method = 'simple'),
                                  realm.tmp = raster::extract(realm,y = backg[,c('lon','lat')], method = 'simple'),
                                  gdd = raster::extract(raster.stack$GDD_2010,y = backg[,c('lon','lat')], method = 'simple'),
                                  temp = raster::extract(raster.stack$Min_Temp_2010,y = backg[,c('lon','lat')], method = 'simple'),
                                  water_bal = raster::extract(raster.stack$Water_2010,y = backg[,c('lon','lat')], method = 'simple'),
                                  precip = raster::extract(raster.stack$Precip_2010,y = backg[,c('lon','lat')], method = 'simple'))
          back.land$cell_num <- backg[,colnames(backg) == 'cell_num']
          
          # And filtering
          # Note that this also gets rid of all non-land areas
          # Because the eco regions and biome data stes we are using are only for terrestrial regions
          backg <- backg[back.land$biome.tmp %in% biome.sum.dat$biome & back.land$realm.tmp %in% realm.sum.dat$realm,]
          
          
        } # End of if statement
        
        
        
        
        if(!is.matrix(backg)) {
          backg <- matrix(0,0,0,nrow=1)
        }
        
        ### For island-dwelling species, this can result in no out-of-sample predictions
        # Preventing the script from running on these species (will return an error)
        # This means that these species will be excluded from our analysis
        # Or at least from the climate change part of our analysis
        # Using an if statement to do this
        if(nrow(backg) > 50) {
          
          # Uncomment if you want to check to be sure
          # na.check <- extract(tmp.stack,
          #                     backg[,c('lon','lat')],
          #                     method = 'simple')
          # na.check.list <- which(!complete.cases(na.check))
          
          # Randomly selecting 10,000 cells for maxent 
          backg_xm <- 
            sample_n(as.data.frame(backg),
                     size = min(10000,nrow(backg)),
                     replace = FALSE)
          
          # Counters to keep the loop running
          row.check = nrow(backg)
          row.check.1 = nrow(backg)
          if(nrow(backg) < (nrow(pres_train) + nrow(pres_test))) {
            row.check.1 = -1
          }
          # Repeating if sample size is not adequately large
          counter = 1
          # Looping to add cells to the data frame
          # Note that this is used for random forest and glm
          while(nrow(backg) < (nrow(pres_train) + nrow(pres_test)) &
                counter < 5 &
                row.check != row.check.1)  {
            # Random points
            backg.tmp <- randomPoints(raster.stack, min(round((esh.raster@nrows*esh.raster@ncols - sum(!is.na(getValues(esh.raster))))/4),(nrow(pres_train) + nrow(pres_test))*4),tryf=10)
            colnames(backg.tmp) = c('lon','lat')
            
            # Limiting to data points where esh = 0 
            # Getting ESH values at the backg points
            # Ideally these will all be NAs
            esh.vals.backg <- raster::extract(x = esh.raster,
                                              y = backg.tmp,
                                              method = 'simple')
            
            # Dropping coordinates that don't contain NAs in the esh raster
            backg.tmp <- backg.tmp[which(is.na(esh.vals.backg)),]
            
            # Getting ecoregion and biome of the selected points
            back.land.tmp <- data.frame(biome.tmp = raster::extract(biome,y = backg.tmp[,c('lon','lat')],method = 'simple'),
                                        realm.tmp = raster::extract(realm,y = backg.tmp[,c('lon','lat')], method = 'simple'))
            
            # And filtering
            # Have to do this a silly way because this is throwing errors
            backg.tmp <- backg.tmp[back.land.tmp$biome.tmp %in% biome.sum.dat$biome & back.land.tmp$realm.tmp %in% realm.sum.dat$realm,]
            
            backg.tmp <- cbind(backg.tmp,-1)
            
            # Adding to backg
            backg <- rbind(backg, backg.tmp)
            
            backg <- unique(backg[,1:2])
            backg <- cbind(backg,-1)
            
            # Making sure loop doesn't get stuck
            # You really want to comment this out when running the scripts!
            # Like really really want to!
            # backg <- backg[1:5,]
            
            # Adding to counter
            counter = counter + 1
            row.check.1
          }
          
          # Now limiting to same number of rows as in pres test + pres train
          backg <- 
            sample_n(as.data.frame(backg),
                     size = min((nrow(pres_test) + nrow(pres_train)),nrow(backg)),
                     replace = FALSE)
          
          # And getting the background data for random forest and the glm model
          backg.rf.glm <-
            raster::extract(raster.stack.bc.xm,
                            backg[,c('lon','lat')])
          
          
          # Dividing into groups for testing and training
          group <- kfold(backg.rf.glm, 5)
          backg_train.rf.glm <- backg.rf.glm[group != 1, ]
          backg_test.rf.glm <- backg.rf.glm[group == 1, ]
          
          # Adding point locations
          backg_train.rf.glm <- cbind(backg_train.rf.glm, backg[group != 1,c('lon','lat')])
          backg_test.rf.glm <- cbind(backg_test.rf.glm, backg[group == 1,c('lon','lat')])
          
          backg_train.rf.glm <- as.data.frame(backg_train.rf.glm)
          backg_test.rf.glm <- as.data.frame(backg_test.rf.glm)
          
          ###
          # Testing accuracy + processing time of bioclim
          # n.points = c()
          # t.vector <- c()
          # auc.vector <- c()
          # for(i in seq(from = 5000, to = 150000, by = 5000)) {
          #   # for(y in seq(from = 50, to = 500, by = 50)) {
          #   for(z in 1:10) {
          #     
          #     
          #     pres.train <- sample_n(pres_train.rf.glm,
          #                            i,
          #                            replace = FALSE)
          #     pres.test <- sample_n(pres_test,
          #                             i/5,
          #                             replace = FALSE)
          #     backg.test <- sample_n(backg_test.rf.glm,
          #                            i/5,
          #                            replace = FALSE)
          #     t1 = Sys.time()
          #     bc <- bioclim(raster.stack, 
          #                   pres.train[,c('x','y')])
          #     
          #     t1 = Sys.time() - t1
          #     # Getting auc
          #     auc_bc <-
          #       evaluate(pres_test[,c('x','y')], backg_test.rf.glm[,c('lon','lat')], bc, raster.stack)
          #     
          #     pres.test <- sample_n(pres_test.rf.glm,
          #                           round(i/8,digits=0),
          #                           replace=FALSE)
          #     abs.test <- sample_n(backg_test.rf.glm,
          #                          round(i/8,digits=0),
          #                          replace=FALSE)
          #     
          #     # Getting AUC
          #     auc_bc <-
          #       evaluate(pres.test[,c('x','y')], abs.test[,c('lon','lat')], xm, raster.stack.bc.xm)
          #     
          #     # and data frame
          #     n.points = c(n.points,i)
          #     t.vector <- c(t.vector,t1)
          #     auc.vector <- c(auc.vector,auc_bc@auc)
          #     # n.trees = c(n.trees, y)
          #   }
          #   # }
          #   
          #   
          #   # And checking progress
          #   print(i)
          #   
          # }
          # tmp.frame <- data.frame(n.points = n.points,
          #                         time = t.vector,
          #                         auc = auc.vector)
          
          # And plotting
          # par(mfrow = c(1,1))
          # 
          # tmp1 = ggplot(dat = tmp.frame, aes(x = n.points, y = time, group = n.points)) +
          #   geom_boxplot() +
          #   labs(title = 'points vs time')
          # 
          # tmp2 = ggplot(dat = tmp.frame, aes(x = n.points, y = auc, group = n.points)) +
          #   geom_boxplot() +
          #   labs(title = 'points vs accuracy')
          # 
          # tmp3 = ggplot(dat = tmp.frame, aes(x = time, y = auc, colour = factor(n.points))) +
          #   geom_point() +
          #   labs(title = 'time vs accuracy')
          # 
          # library(cowplot)
          # tmp3 = plot_grid(tmp1,tmp2,tmp3,nrow=3)
          # tmp3
          # 
          # ggplot(dat = tmp.frame, aes(x = n.points, y = auc)) +
          #   geom_point() +
          #   geom_smooth()
          # 
          # ggsave("/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Performance Testing/Aciononyx jubatus Bioclim Model.pdf",
          #        width=10,height=10)
          # write.csv(tmp.frame,
          #           "/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Performance Testing/Aciononyx jubatus Bioclim Model.csv")
          # 
          
          ###
          # Bioclim
          # For performance vs accuracy tradeoff, limiting to 15,000 presence points (e.g. nrow = 10,000)
          # This means that the testing data sets each contain 3,000 data points
          
          if(nrow(pres_train) <= 15000) {
            pres.train.bc <- pres_train
            pres.test.bc <- pres_test
            backg.test.bc <- backg_test.rf.glm
          }
          if(nrow(pres_train) > 15000) {
            pres.train.bc <- sample_n(pres_train,
                                      15000,
                                      replace = FALSE)
            pres.test.bc <- sample_n(pres_test,
                                     min(3000,nrow(pres_test)),
                                     replace = FALSE)
            backg.test.bc <- sample_n(backg_test.rf.glm,
                                      min(3000,nrow(backg_test.rf.glm)),
                                      replace = FALSE)
          }
          
          bc <- bioclim(raster.stack, 
                        pres.train.bc[,c('x','y')])
          
          # Getting auc
          auc_bc <-
            evaluate(pres.test.bc[,c('x','y')], backg.test.bc[,c('lon','lat')], bc, raster.stack)
          # Getting threshold of values within the existing ESH
          bc.tr <- threshold(auc_bc)
          
          
          # bc.tmp <-
          #   bc.pred
          # bc.tmp[bc.tmp < bc.tr] <- NA
          # 
          # plot(bc.pred)
          # plot(esh.raster, add = TRUE)
          # plot(bc.tmp, col = 'black', add = TRUE)
          
          ###
          # Domain next
          # dm <- domain(raster.stack, 
          #              pres_train[,c('x','y')])
          # 
          # # Getting AUC
          # auc_dm <-
          #   evaluate(pres_test[,c('x','y')], backg_test, dm, raster.stack)
          # # And getting threshold
          # dm.tr <-
          #   threshold(auc_dm,'spec_sens')
          
          ###
          # Maxent
          # THe commented out code immediately below tests processing time vs accuracy vs number of points
          # Do not run unless specifically testing for this
          
          # 
          # n.points = c()
          # t.vector <- c()
          # auc.vector <- c()
          # for(i in seq(from = 500, to = 15000, by = 500)) {
          #   # for(y in seq(from = 50, to = 500, by = 50)) {
          #     for(z in 1:10) {
          #       
          # 
          #       pres.train <- sample_n(pres_train.rf.glm,
          #                           i/2,
          #                           replace = FALSE)
          #       backg.train <- sample_n(backg_train.rf.glm,
          #                               i/2,
          #                               replace = FALSE)
          #       t1 = Sys.time()
          #       xm <- maxent(x = raster.stack.bc.xm,
          #                    p = pres.train[,c('x','y')],
          #                    a = backg.train[,c('lon','lat')])
          #       
          #       
          #       t1 = Sys.time() - t1
          # 
          #       pres.test <- sample_n(pres_test.rf.glm,
          #                                 round(i/8,digits=0),
          #                                 replace=FALSE)
          #       abs.test <- sample_n(backg_test.rf.glm,
          #                                round(i/8,digits=0),
          #                                replace=FALSE)
          # 
          #       # Getting AUC
          #       auc_xm <-
          #         evaluate(pres.test[,c('x','y')], abs.test[,c('lon','lat')], xm, raster.stack.bc.xm)
          # 
          #       # and data frame
          #       n.points = c(n.points,i)
          #       t.vector <- c(t.vector,t1)
          #       auc.vector <- c(auc.vector,auc_xm@auc)
          #       # n.trees = c(n.trees, y)
          #     }
          #   # }
          # 
          # 
          #   # And checking progress
          #   print(i)
          # 
          # }
          # tmp.frame <- data.frame(n.points = n.points,
          #                         time = t.vector,
          #                         auc = auc.vector)
          # 
          # 
          # # And plotting
          # par(mfrow = c(1,1))
          # 
          # tmp1 = ggplot(dat = tmp.frame, aes(x = n.points, y = time, group = n.points)) +
          #   geom_boxplot() +
          #   labs(title = 'points vs time')
          # 
          # tmp2 = ggplot(dat = tmp.frame, aes(x = n.points, y = auc, group = n.points)) +
          #   geom_boxplot() +
          #   labs(title = 'points vs accuracy')
          # 
          # tmp3 = ggplot(dat = tmp.frame, aes(x = time, y = auc, colour = factor(n.points))) +
          #   geom_point() +
          #   labs(title = 'time vs accuracy')
          # 
          # library(cowplot)
          # tmp3 = plot_grid(tmp1,tmp2,tmp3,nrow=3)
          # tmp3
          # 
          # ggsave("/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Performance Testing/Maxent Model.pdf",
          #        width=10,height=10)
          # 
          # write.csv(tmp.frame,
          #           "/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Performance Testing/Aciononyx jubatus Maxent Model.csv")
          # 
          
          
          
          # Limiting for processing time and performance
          if(nrow(pres_train.rf.glm) <= 7500) {
            pres.train.xm <- pres_train.rf.glm
            backg.train.xm <- backg_train.rf.glm
            pres.test.xm <- pres_test.rf.glm
            backg.test.xm <- backg_test.rf.glm
          }
          if(nrow(pres_train.rf.glm) > 7500) {
            pres.train.xm <- sample_n(pres_train.rf.glm,min(7500,nrow(pres_train.rf.glm)),replace = FALSE)
            backg.train.xm <- sample_n(backg_train.rf.glm,min(7500,nrow(backg_train.rf.glm)),replace=FALSE)
            pres.test.xm <- sample_n(pres_test.rf.glm,min(1875,nrow(pres_test.rf.glm)),replace=FALSE)
            backg.test.xm <- sample_n(backg_test.rf.glm,min(1875,nrow(backg_test.rf.glm)),replace=FALSE)
          }
          # xm <- maxent(x = raster.stack.bc.xm,
          #              p = pres.train.xm[,c('x','y')],
          #              a = backg.train.xm[,c('lon','lat')])
          
          xm <- maxnet(p=c(rep(1, nrow(pres.train.xm)),rep(0, nrow(backg.train.xm))), 
                       data=rbind(pres.train.xm[,names(raster.stack.bc.xm)],backg.train.xm[,names(raster.stack.bc.xm)]), 
                       f=maxnet.formula(p=p, data=data, classes='lpq'), 
                       addsamplestobackground = FALSE)
          
          # Getting AUC
          auc_xm <- dismo::evaluate(pres.test.xm[,c('x','y')], backg.test.xm[,c('lon','lat')], xm, raster.stack.bc.xm, type='cloglog')
          # auc_xm <-
          #   evaluate(pres.test.xm[,c('x','y')], backg.test.xm[,c('lon','lat')], xm, raster.stack.bc.xm)
          # And getting threshold
          xm.tr <- threshold(auc_xm)
          
          ###
          # Random Forest
          # Setting up training data
          # Converting NAs to 0s in the esh data
          backg_train.rf.glm$esh <- 0
          backg_test.rf.glm$esh <- 0
          
          # Binding data inside and outside ESH
          # dplyr::select to reorganize column orders
          train <- rbind(dplyr::select(pres_train.rf.glm,esh,GDD,Precip,Water_Bal,Min_temp,GDD_sq,Precip_sq,Water_Bal_sq,Min_temp_sq,x,y),
                         dplyr::select(backg_train.rf.glm,
                                       esh,
                                       GDD = GDD_2010,
                                       Precip = Precip_2010,
                                       Water_Bal = Water_2010,
                                       Min_temp = Min_Temp_2010,
                                       GDD_sq = GDD_2010_sq,
                                       Precip_sq = Precip_2010_sq,
                                       Water_Bal_sq = Water_2010_sq,
                                       Min_temp_sq = Min_Temp_2010_sq,
                                       x = lon,
                                       y = lat))
          test <- rbind(dplyr::select(pres_test.rf.glm,esh,GDD,Precip,Water_Bal,Min_temp,GDD_sq,Precip_sq,Water_Bal_sq,Min_temp_sq,x,y),
                        dplyr::select(backg_test.rf.glm,
                                      esh,
                                      GDD = GDD_2010,
                                      Precip = Precip_2010,
                                      Water_Bal = Water_2010,
                                      Min_temp = Min_Temp_2010,
                                      GDD_sq = GDD_2010_sq,
                                      Precip_sq = Precip_2010_sq,
                                      Water_Bal_sq = Water_2010_sq,
                                      Min_temp_sq = Min_Temp_2010_sq,
                                      x = lon,
                                      y = lat))
          
          # Getting complete cases
          # Models throw an error if we don't do this
          train <-
            train[complete.cases(train),]
          test <- test[complete.cases(test),]
          
          # Running the model
          model = esh ~ GDD + Precip + Water_Bal + Min_temp
          
          ###
          # Chunk of code below tests performance
          # 
          # n.points = c()
          # t.vector <- c()
          # auc.vector <- c()
          # n.trees = c()
          # for(i in seq(from = 500, to = 15000, by = 500)) {
          #   for(y in seq(from = 50, to = 500, by = 50)) {
          #     for(z in 1:10) {
          #     
          #     
          #     tmp.train <- sample_n(train,
          #                            i,
          #                            replace = FALSE)
          #     tmp.test <- sample_n(test,
          #                             i/4,
          #                             replace = FALSE)
          #     t1 = Sys.time()
          #     # Getting AUC
          #     rf <- randomForest(model,
          #                        data = tmp.train[,c('esh','GDD','Precip','Water_Bal','Min_temp')],
          #                        ntree = y)
          #     t1 = Sys.time() - t1
          #     # Getting accuracy
          #     tmp.test.pres <- sample_n(pres_test.rf.glm,
          #                              round(i/8),
          #                              replace = FALSE)
          #     tmp.test.backg <- sample_n(backg_test.rf.glm,
          #                               round(i/8),
          #                               replace = FALSE)
          #     
          #     tmp.test.backg <- dplyr::rename(tmp.test.backg,
          #                                    GDD = GDD_2010,
          #                                    Min_temp = Min_Temp_2010,
          #                                    Precip = Precip_2010,
          #                                    Water_Bal = Water_2010,
          #                                    GDD_sq = GDD_2010_sq,
          #                                    Min_temp_sq = Min_Temp_2010_sq,
          #                                    Precip_sq = Precip_2010_sq,
          #                                    Water_Bal_sq = Water_2010_sq,
          #                                    x = lon,
          #                                    y = lat)
          #     
          #     auc_rf <-
          #       evaluate(p = tmp.test.pres[complete.cases(tmp.test.pres),],
          #                a = tmp.test.backg[complete.cases(tmp.test.backg),], 
          #                model = rf)
          #     
          #     
          #     
          #     # and data frame
          #     n.points = c(n.points,i)
          #     t.vector <- c(t.vector,t1)
          #     auc.vector <- c(auc.vector,auc_rf@auc)
          #     n.trees = c(n.trees, y)
          #     }
          #   }
          #   
          #   
          #   # And checking progress
          #   print(i)
          #   
          # }
          # tmp.frame <- data.frame(n.points = n.points,
          #                         time = t.vector,
          #                         auc = auc.vector,
          #                         n.trees = n.trees)
          # 
          # 
          # # And plotting
          # par(mfrow = c(1,1))
          # 
          # tmp1 = ggplot(dat = tmp.frame, aes(x = n.points, y = time, group = n.points)) +
          #   geom_boxplot() +
          #   facet_wrap(.~n.trees) +
          #   labs(title = 'points vs time facet n.trees')
          # 
          # tmp2 = ggplot(dat = tmp.frame, aes(x = n.points, y = auc, group = n.points)) +
          #   geom_boxplot() +
          #   facet_wrap(.~n.trees) +
          #   labs(title = 'points vs accuracy facet n.trees')
          #   
          #   tmp3 = ggplot(dat = tmp.frame, aes(x = n.trees, y = time, group = n.trees)) +
          #     geom_boxplot() +
          #     facet_wrap(.~n.points) +
          #     labs(title = 'n.trees vs time facet points')
          #   
          #   tmp4 = ggplot(dat = tmp.frame, aes(x = n.trees, y = auc, group = n.trees)) +
          #     geom_boxplot() +
          #     facet_wrap(.~n.points) +
          #     labs(title = 'n.trees vs accuracy facet points')
          # 
          # library(cowplot)
          # tmp5 = plot_grid(tmp1,tmp2,tmp3,tmp4,nrow=2,ncol=2)
          # tmp5
          # 
          # ggsave("/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Performance Testing/Aciononyx jubatus Random Forest Model.pdf",
          #        width=20,height=10)
          # 
          # tmp1 = ggplot(dat = tmp.frame, aes(x = time, y = auc, group = n.trees)) +
          #   geom_boxplot() +
          #   facet_wrap(.~n.points)+
          #   labs(title = 'time vs accuracy facet points')
          # 
          # tmp2 = ggplot(dat = tmp.frame, aes(x = time, y = auc, group = n.points)) +
          #   geom_boxplot() +
          #   facet_wrap(.~n.trees)+
          #   labs(title = 'time vs accuracy facet n.trees')
          # 
          # tmp3 = plot_grid(tmp1,tmp2,nrow=2)
          # tmp3
          # 
          # 
          # 
          # write.csv(tmp.frame,
          #           "/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Performance Testing/Aciononyx jubatus Random Forest Model.csv")
          # 
          # 
          
          
          
          ###
          # And for performance/processing time tradeoff of random forest
          # Randomly selecting 5000 absence and presence points
          # Creating dummy variables
          if(nrow(train)<=10000) {
            train.rf <- train
            test.rf <- test
            pres_test.rf <- pres_test.rf.glm
            backg_test.rf <- backg_test.rf.glm
            
            backg_test.rf <- dplyr::rename(backg_test.rf,
                                           GDD = GDD_2010,
                                           Min_temp = Min_Temp_2010,
                                           Precip = Precip_2010,
                                           Water_Bal = Water_2010,
                                           GDD_sq = GDD_2010_sq,
                                           Min_temp_sq = Min_Temp_2010_sq,
                                           Precip_sq = Precip_2010_sq,
                                           Water_Bal_sq = Water_2010_sq,
                                           x = lon,
                                           y = lat)
          }
          if(nrow(train)>10000){
            train.rf <- rbind(sample_n(train[train$esh == 0,],min(5000,sum(train$esh==0)),replace = FALSE),
                              sample_n(train[train$esh == 1,],min(5000,sum(train$esh==1)),replace = FALSE))
            
            test.rf <- rbind(sample_n(test[test$esh==0,],min(1250,test$esh==0),replace = FALSE),
                             sample_n(test[test$esh==1,],min(1250,test$esh==1),replace = FALSE))
            
            pres_test.rf <- sample_n(pres_test.rf.glm,
                                     min(1250,nrow(pres_test.rf.glm)),
                                     replace = FALSE)
            backg_test.rf <- sample_n(backg_test.rf.glm,
                                      min(1250,nrow(backg_test.rf.glm)),
                                      replace = FALSE)
            
            backg_test.rf <- dplyr::rename(backg_test.rf,
                                           GDD = GDD_2010,
                                           Min_temp = Min_Temp_2010,
                                           Precip = Precip_2010,
                                           Water_Bal = Water_2010,
                                           GDD_sq = GDD_2010_sq,
                                           Min_temp_sq = Min_Temp_2010_sq,
                                           Precip_sq = Precip_2010_sq,
                                           Water_Bal_sq = Water_2010_sq,
                                           x = lon,
                                           y = lat)
          }
          
          
          
          # This will return a warning about the output having < 5 unique responses
          # This is ok, seeing as our response variable is binary.
          rf <- randomForest(model,
                             data = train.rf[,c('esh','GDD','Precip','Water_Bal','Min_temp')],
                             ntree = 500)
          
          # Changing names to test the model
          backg_test.rf.glm <- dplyr::rename(backg_test.rf.glm,
                                             GDD = GDD_2010,
                                             Min_temp = Min_Temp_2010,
                                             Precip = Precip_2010,
                                             Water_Bal = Water_2010,
                                             GDD_sq = GDD_2010_sq,
                                             Min_temp_sq = Min_Temp_2010_sq,
                                             Precip_sq = Precip_2010_sq,
                                             Water_Bal_sq = Water_2010_sq,
                                             x = lon,
                                             y = lat)
          
          
          # Now gettin auc for the model
          auc_rf <-
            evaluate(p = pres_test.rf[complete.cases(pres_test.rf),],
                     a = backg_test.rf[complete.cases(backg_test.rf),], 
                     model = rf)
          # Getting threshold
          rf.tr <- threshold(auc_rf)
          # plot(rf.pred>rf.tr)
          
          
          
          
          ###
          # GLM Model
          # THe commented out code immediately below tests processing time vs accuracy vs number of points
          # Do not run unless specifically testing for this
          
          
          # n.points = c()
          # t.vector <- c()
          # auc.vector <- c()
          # # n.trees <- c()
          # for(i in seq(from = 5000, to = 100000, by = 5000)) {
          #   # for(y in seq(from = 50, to = 500, by = 50)) {
          #     for(z in 1:10) {
          #       tmp.dat <- sample_n(train,
          #                           i,
          #                           replace = FALSE)
          #       t1 = Sys.time()
          #     
          #       gm.gaus <-
          #         glm(esh ~ GDD + Precip + Water_Bal + Min_temp + GDD_sq + Precip_sq + Water_Bal_sq + Min_temp_sq,
          #             family = gaussian(link = "identity"),
          #             data=tmp.dat)
          #       t1 = Sys.time() - t1
          # 
          #       tmp.pres.test <- sample_n(pres_test.rf.glm,
          #                                 round(i/8,digits=0),
          #                                 replace=FALSE)
          #       tmp.abs.test <- sample_n(backg_test.rf.glm,
          #                                round(i/8,digits=0),
          #                                replace=FALSE)
          #       tmp.abs.test <- dplyr::select(tmp.abs.test,
          #                                     esh,
          #                                     GDD=GDD_2010,
          #                                     Precip = Precip_2010,
          #                                     Water_Bal = Water_2010,
          #                                     Min_temp = Min_Temp_2010,
          #                                     GDD_sq=GDD_2010_sq,
          #                                     Precip_sq = Precip_2010_sq,
          #                                     Water_Bal_sq = Water_2010_sq,
          #                                     Min_temp_sq = Min_Temp_2010_sq,
          #                                     x = lon,
          #                                     y=lat)
          # 
          #       auc_gm <-
          #         evaluate(p = tmp.pres.test[complete.cases(tmp.pres.test),],
          #                  a = tmp.abs.test[complete.cases(tmp.abs.test),],
          #                  model = gm.gaus)
          # 
          #       # and data frame
          #       n.points = c(n.points,i)
          #       t.vector <- c(t.vector,t1)
          #       auc.vector <- c(auc.vector,auc_gm@auc)
          #       # n.trees = c(n.trees, y)
          #     }
          #   # }
          # 
          # 
          #   # And checking progress
          #   print(i)
          # 
          # }
          # tmp.frame <- data.frame(n.points = n.points,
          #                         time = t.vector,
          #                         auc = auc.vector)
          # 
          # write.csv(tmp.frame,
          #           "/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Performance Testing/Aciononyx jubatus GLM Model.csv")
          # 
          # 
          # # And plotting
          # par(mfrow = c(1,1))
          # 
          # tmp1 = ggplot(dat = tmp.frame, aes(x = n.points, y = time, group = n.points)) +
          #   geom_boxplot() +
          #   labs(title = 'points vs time')
          # 
          # tmp2 = ggplot(dat = tmp.frame, aes(x = n.points, y = auc, group = n.points)) +
          #   geom_boxplot() +
          #   labs(title = 'points vs accuracy')
          # 
          # tmp3 = ggplot(dat = tmp.frame, aes(x = time, y = auc, colour = factor(n.points))) +
          #   geom_point() +
          #   labs(title = 'time vs accuracy')
          # 
          # library(cowplot)
          # tmp3 = plot_grid(tmp1,tmp2,tmp3,nrow=3)
          # tmp3
          # 
          # ggsave("/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Performance Testing/Aciononyx jubatus GLM Model Big Data.pdf",
          #        width=10,height=15)
          # 
          
          
          # Limiting data set for processing and performance
          if(nrow(train)<=150000) {
            train.glm <- train
            pres.test.glm <- pres_test.rf.glm
            backg.test.glm <- backg_test.rf.glm
          }
          if(nrow(train) > 150000) {
            train.glm <- rbind(sample_n(train[train$esh==0,],min(75000,nrow(train[train$esh==0,])),replace=FALSE),
                               sample_n(train[train$esh==1,],min(75000,nrow(train[train$esh==1,])),replace=FALSE))
            pres.test.glm <- sample_n(pres_test.rf.glm,min(18750,nrow(pres_test.rf.glm)),replace=FALSE)
            backg.test.glm <- sample_n(backg_test.rf.glm,min(18750,nrow(backg_test.rf.glm)),replace=FALSE)
          }
          
          gm.gaus <-
            glm(esh ~ GDD + Precip + Water_Bal + Min_temp + GDD_sq + Precip_sq + Water_Bal_sq + Min_temp_sq,             
                family = gaussian(link = "identity"), 
                data=train.glm)
          # Testing model evaluation
          auc_gm <-
            evaluate(pres.test.glm, 
                     backg.test.glm,  
                     gm.gaus)
          # Getting threshold
          gm.tr <- threshold(auc_gm)
          
          ###
          # And list of things to return
          mod.list <-
            list(bc, xm, rf, gm.gaus)
          names(mod.list) <- c('bioclim','maxent','random.forest','glm.gaussian')
          
          auc.list <-
            list(auc_bc, auc_xm, auc_rf, auc_gm)
          names(auc.list) <- c('bioclim','maxent','random.forest','glm.gaussian')
          
          tr.list <-
            list(bc.tr, xm.tr, rf.tr, gm.tr)
          names(tr.list) <- c('bioclim','maxent','random.forest','glm.gaussian')
          
          return.list <-
            list(mod.list,
                 auc.list,
                 tr.list)
          
          names(return.list) <-
            c('Models','AUC','Thresholds')
          
          # And returning this list
          return(return.list) 
        } # End of if statement preventing code from running on species with small out-of-sample data sets
      } # End of if statement preventing code from running on species with no range overlap between biome and realm  
    }
  }

###
# Second function predicts models based on RCP 4.5 and RCP 8.5
###
# Second function predicts models based on RCP 4.5 and RCP 8.5
sdm.pred.fun <-
  function(#raster.stack.pred,
    raster.stack.bc.xm.pred,
    #df.dat,
    df.dat.bc.xm,
    mod.list,
    path.sdms.write.tmp,
    path.sdms.write.mosaiced.tmp,
    split.current.extents,
    species_name_save) {
    
    ### 
    # Limiting to models with an auc > 0.8 or < .99
    # Models with auc < .8 are poorly fitting models
    # Models with auc >= .99 are over fitting models
    tmp.list <- c()
    for(z in 1:4) {
      if(mod.list$AUC[[z]]@auc >= 0.75 &
         mod.list$AUC[[z]]@auc < 1) {
        tmp.list<-c(tmp.list,z)
      }
    }
    
    # if(length(tmp.list >= 1)) {
    #   mod.list$AUC <- mod.list$AUC[tmp.list]
    #   mod.list$Thresholds <- mod.list$Thresholds[tmp.list]
    #   mod.list$Models <- mod.list$Models[tmp.list]
    # }
    
    # Below lines are herefor troubleshooting
    # mod.list$AUC <- mod.list$AUC[1]
    #   mod.list$Thresholds <- mod.list$Thresholds[1]
    #   mod.list$Models <- mod.list$Models[1]
    
    if(length(tmp.list) >= 1) {
      # Creating empty list
      # Used for a logical test
      models<-c()
      # Adding auc objects to this list to be used later
      auc.list <- list()
      # Adding threshold objects to this list to be used later
      threshold.list <- list()
      
      ### Splitting main raster into multiple smaller ones
      # Doing this so that the computer doesn't run out of RAM when predicting future SDMs
      # Each raster has a max of ~3,000,000 cells
      
      ###
      # First making rasters for the individual predictions
      if('bioclim' %in% names(mod.list$Models)) {
        # Bioclim
        # t1 = Sys.time()
        
        for(i in 1:length(split.current.extents)) {
          # Logic check - 
          # Predict the SDM, or create a dummy SDM
          # Doing this based on AUC of the models
          # And to save time
          if('bioclim' %in% names(mod.list$Models)[tmp.list]) { # Predict
            pred.raster.stack <- crop(raster.stack.bc.xm.pred[[c('GDD_2010','Min_Temp_2010','Precip_2010','Water_2010')]], split.current.extents[[i]])
            
            bc.pred.tmp <- predict(pred.raster.stack, 
                                   mod.list$Models[['bioclim']], 
                                   ext = extent(pred.raster.stack))
          } else { # Else create a raster full of NAs based on split current extent
            bc.pred.tmp <- split.current.extents[[i]]
            bc.pred.tmp[!is.na(bc.pred.tmp)] <- NA
          } # End if else statement check whether to forecast bioclim
          
          # Writing if the species has more than one tile in their habitat range
          if(length(split.current.extents) %in% 1) {
            # And saving
            writeRaster(bc.pred.tmp,
                        paste0(path.sdms.write.tmp,
                               '/',
                               species_name_save,
                               '_bioclim.tif'),
                        overwrite = TRUE)
            
            # rm(out.raster)
            rm(bc.pred.tmp)
          } else {
            # Else write the raster tile
            writeRaster(bc.pred.tmp,
                        paste0(path.sdms.write.tmp,
                               '/',
                               species_name_save,
                               '_tmp_sdm_tile_',
                               i,
                               ".tif"),
                        overwrite = TRUE)
            rm(bc.pred.tmp)
            
            if(i %in% length(split.current.extents)) {
              # If on the last tile of the habitat range, then mosaic all other tiles
              
              # Getting list of rasters to mosaic
              gdalfile <- 
                list.files(path = path.sdms.write.tmp,
                           pattern = species_name_save,
                           full.names = TRUE) %>% .[grepl('_tmp',.)]
              
              # Place to write the sdm
              dst_dataset = paste0(path.sdms.write.tmp,'/',species_name_save,'_bioclim.tif')
              
              # And mosaicing
              mosaic_rasters(gdalfile, dst_dataset)
            } # End statement to merge
          } # End statement looping through tiles
        } # End bioclim
        
        # Removing files in the folders
        tmp.files <- list.files(path.sdms.write.tmp,full.names=TRUE, pattern = paste0(species_name_save,'_tmp'))
        file.remove(tmp.files)
        
        tmp.files <- list.files(path.sdms.write.mosaiced.tmp,full.names=TRUE, pattern = species_name_save)
        file.remove(tmp.files)
        
        # models <- stack(bc.merge)
        auc.list['bioclim'] <- mod.list$AUC$bioclim
        threshold.list['bioclim'] <- mod.list$Thresholds$bioclim
        
        # rm(bc.merge)
      }
      
      
      
      
      # Domain
      # dm.pred <- predict(raster.stack, 
      #                    mod.list$Models$domain, 
      #                    ext = extent(raster.stack), 
      #                    progress = 'text')
      if('maxent' %in% names(mod.list$Models)) {
        
        # Predicting maxent
        # Unfortunately cannot do this in parallel through R... =|
        # Alright, so can do this in parallel, but need to have it loop through mclapply
        for(i in 1:length(split.current.extents)) {
          # Logic check - 
          # Predict the SDM, or create a dummy SDM
          # Doing this based on AUC of the models
          # And to save time
          if('maxent' %in% names(mod.list$Models)[tmp.list]) { # Predict
            pred.raster.stack <- crop(raster.stack.bc.xm.pred, split.current.extents[[i]])
            
            xm.pred.tmp <- predict(pred.raster.stack, 
                                   mod.list$Models[['maxent']], 
                                   ext = extent(pred.raster.stack),
                                   type = 'cloglog')
          } else { # Else create a raster full of NAs based on split current extent
            xm.pred.tmp <- split.current.extents[[i]]
            xm.pred.tmp[!is.na(xm.pred.tmp)] <- NA
          } # End if else statement check whether to forecast maxent
          
          # Writing if the species has more than one tile in their habitat range
          if(length(split.current.extents) %in% 1) {
            # And saving
            writeRaster(xm.pred.tmp,
                        paste0(path.sdms.write.tmp,
                               '/',
                               species_name_save,
                               '_maxent.tif'),
                        overwrite = TRUE)
            
            # rm(out.raster)
            rm(xm.pred.tmp)
          } else {
            # Else write the raster tile
            writeRaster(xm.pred.tmp,
                        paste0(path.sdms.write.tmp,
                               '/',
                               species_name_save,
                               '_tmp_sdm_tile_',
                               i,
                               ".tif"),
                        overwrite = TRUE)
            rm(xm.pred.tmp)
            
            if(i %in% length(split.current.extents)) {
              # If on the last tile of the habitat range, then mosaic all other tiles
              
              # Getting list of rasters to mosaic
              gdalfile <- 
                list.files(path = path.sdms.write.tmp,
                           pattern = species_name_save,
                           full.names = TRUE) %>% .[grepl('_tmp',.)]
              
              # Place to write the sdm
              dst_dataset = paste0(path.sdms.write.tmp,'/',species_name_save,'_maxent.tif')
              
              # And mosaicing
              mosaic_rasters(gdalfile, dst_dataset)
            } # End statement to merge
          } # End of if loop across tiles
        }
        
        # xm.list <-
        #   chunk2(1:length(split.current.extents),5)
        # 
        # tryCatch(
        #   expr = {
        #     withTimeout({lapply(xm.list,xm.lapply)}, 
        #                 timeout = (30*length(xm.list[[1]])+5))
        #   }, 
        #   TimeoutException = function(ex) cat("Timeout. Skipping.\n")
        # )
        # 
        # xm.merge <-
        #   raster(paste0(path.sdms.write.tmp,
        #                 '/',
        #                 species_name_save,
        #                 '_maxent.tif'))
        # crs(xm.merge) <- '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'
        
        
        # Removing files in the folders
        tmp.files <- list.files(path.sdms.write.tmp,full.names=TRUE, pattern = paste0(species_name_save, '_tmp'))
        file.remove(tmp.files)
        
        tmp.files <- list.files(path.sdms.write.mosaiced.tmp,full.names=TRUE, pattern = paste0(species_name_save,'_tmp'))
        file.remove(tmp.files)
        
        # And saving
        auc.list['maxent'] <- mod.list$AUC$maxent
        threshold.list['maxent'] <- mod.list$Thresholds$maxent
        
        # rm(xm.merge)
        # rm(out.raster)
        # rm(xm.merge)
      }
      
      if('random.forest' %in% names(mod.list$Models)) {
        # Need to create a new data frame with the new extents
        df.dat = list()
        
        # split.stack.bc.xm <- list()
        # split.stack.bc.xm[[1]] <- crop(raster.stack.bc.xm.pred, split.current.extents[[1]])
        # # Getting column indices for these
        # tmp.names = names(split.stack.bc.xm[[1]])
        # gdd <- grep("GDD",tmp.names)[!(grep("GDD",tmp.names) %in% grep("_sq",tmp.names))]
        # precip <- grep("Precip",tmp.names)[!(grep("Precip",tmp.names) %in% grep("_sq",tmp.names))]
        # water <- grep("Water",tmp.names)[!(grep("Water",tmp.names) %in% grep("_sq",tmp.names))]
        # temp <- grep("emp",tmp.names)[!(grep("emp",tmp.names) %in% grep("_sq",tmp.names))]
        # gdd_sq <- grep("GDD",tmp.names)[(grep("GDD",tmp.names) %in% grep("_sq",tmp.names))]
        # precip_sq <- grep("Precip",tmp.names)[(grep("Precip",tmp.names) %in% grep("_sq",tmp.names))]
        # water_sq <- grep("Water",tmp.names)[(grep("Water",tmp.names) %in% grep("_sq",tmp.names))]
        # temp_sq <- grep("emp",tmp.names)[(grep("emp",tmp.names) %in% grep("_sq",tmp.names))]
        
        for(i in 1:length(split.current.extents)) {
          # Logic check - 
          # Predict the SDM, or create a dummy SDM
          # Doing this based on AUC of the models
          # And to save time
          if('random.forest' %in% names(mod.list$Models)[tmp.list]) { # Predict
            pred.raster.stack <- crop(raster.stack.bc.xm.pred, split.current.extents[[i]])
            
            df.dat.tmp <- 
              data.frame(GDD = getValues(pred.raster.stack[[grep("GDD",names(pred.raster.stack))[!(grep("GDD",names(pred.raster.stack)) %in% grep("_sq",names(pred.raster.stack)))]]]),
                         Precip = getValues(pred.raster.stack[[grep("Precip",names(pred.raster.stack), ignore.case = TRUE)[!(grep("Precip",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                         Water_Bal = getValues(pred.raster.stack[[grep("Water",names(pred.raster.stack), ignore.case = TRUE)[!(grep("Water",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                         Min_temp = getValues(pred.raster.stack[[grep("Temp",names(pred.raster.stack), ignore.case = TRUE)[!(grep("Temp",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                         GDD_sq = getValues(pred.raster.stack[[grep("GDD",names(pred.raster.stack))[(grep("GDD",names(pred.raster.stack)) %in% grep("_sq",names(pred.raster.stack)))]]]),
                         Precip_sq = getValues(pred.raster.stack[[grep("Precip",names(pred.raster.stack), ignore.case = TRUE)[(grep("Precip",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                         Water_Bal_sq = getValues(pred.raster.stack[[grep("Water",names(pred.raster.stack), ignore.case = TRUE)[(grep("Water",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                         Min_temp_sq = getValues(pred.raster.stack[[grep("Temp",names(pred.raster.stack), ignore.case = TRUE)[(grep("Temp",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]))
            # Creating tmp raster
            rf.pred.tmp <-
              raster(matrix(predict(mod.list$Models[['random.forest']],
                                    df.dat.tmp[,c('GDD','Precip','Water_Bal','Min_temp')]),
                            nrow = pred.raster.stack@nrows,
                            ncol = pred.raster.stack@ncols,
                            byrow = TRUE),
                     crs = crs(pred.raster.stack))
            extent(rf.pred.tmp) <- extent(pred.raster.stack)
          } else { # Else create a raster full of NAs based on split current extent
            rf.pred.tmp <- split.current.extents[[i]]
            rf.pred.tmp[!is.na(rf.pred.tmp)] <- NA
          } # End if else statement check whether to forecast random forest
          
          
          
          # Writing if the species has more than one tile in their habitat range
          if(length(split.current.extents) %in% 1) {
            # And saving
            writeRaster(rf.pred.tmp,
                        paste0(path.sdms.write.tmp,
                               '/',
                               species_name_save,
                               '_random.forest.tif'),
                        overwrite = TRUE)
            
            # rm(out.raster)
            rm(rf.pred.tmp)
          } else {
            # Else write the raster tile
            writeRaster(rf.pred.tmp,
                        paste0(path.sdms.write.tmp,
                               '/',
                               species_name_save,
                               '_tmp_sdm_tile_',
                               i,
                               ".tif"),
                        overwrite = TRUE)
            rm(rf.pred.tmp)
            
            if(i %in% length(split.current.extents)) {
              # If on the last tile of the habitat range, then mosaic all other tiles
              
              # Getting list of rasters to mosaic
              gdalfile <- 
                list.files(path = path.sdms.write.tmp,
                           pattern = species_name_save,
                           full.names = TRUE) %>% .[grepl('_tmp',.)]
              
              # Place to write the sdm
              dst_dataset = paste0(path.sdms.write.tmp,'/',species_name_save,'_random.forest.tif')
              
              # And mosaicing
              mosaic_rasters(gdalfile, dst_dataset)
            } # End statement to merge
          } # End of if statement for each tile
        } # End of loop across tiles
        
        
        
        # rf.merge <- 
        #   raster(paste0(path.sdms.write.tmp,
        #                 '/',
        #                 species_name_save,
        #                 '_random.forest.tif'))
        # crs(rf.merge) <- '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'
        # Removing files in the folders
        tmp.files <- list.files(path.sdms.write.tmp,full.names=TRUE, pattern = paste0(species_name_save,'_tmp'))
        file.remove(tmp.files)
        
        tmp.files <- list.files(path.sdms.write.mosaiced.tmp,full.names=TRUE, pattern = paste0(species_name_save,'_tmp'))
        file.remove(tmp.files)
        
        # And saving outputs
        # if(!is.Raster(models)) {
        #   models <- stack(rf.merge)
        # } else {
        #   models <- stack(models,rf.merge)
        # }
        auc.list['random.forest'] <- mod.list$AUC$random.forest
        threshold.list['random.forest'] <- mod.list$Thresholds$random.forest
        
        # rm(rf.merge)
      }
      
      if('glm.gaussian' %in% names(mod.list$Models)) {
        # Need to create a new data frame with the new extents
        df.dat = list()
        
        split.stack.bc.xm <- list()
        split.stack.bc.xm[[1]] <- crop(raster.stack.bc.xm.pred, split.current.extents[[1]])
        # Getting column indices for these
        tmp.names = names(split.stack.bc.xm[[1]])
        gdd <- grep("GDD",tmp.names)[!(grep("GDD",tmp.names) %in% grep("_sq",tmp.names))]
        precip <- grep("Precip",tmp.names)[!(grep("Precip",tmp.names) %in% grep("_sq",tmp.names))]
        water <- grep("Water",tmp.names)[!(grep("Water",tmp.names) %in% grep("_sq",tmp.names))]
        temp <- grep("emp",tmp.names)[!(grep("emp",tmp.names) %in% grep("_sq",tmp.names))]
        gdd_sq <- grep("GDD",tmp.names)[(grep("GDD",tmp.names) %in% grep("_sq",tmp.names))]
        precip_sq <- grep("Precip",tmp.names)[(grep("Precip",tmp.names) %in% grep("_sq",tmp.names))]
        water_sq <- grep("Water",tmp.names)[(grep("Water",tmp.names) %in% grep("_sq",tmp.names))]
        temp_sq <- grep("emp",tmp.names)[(grep("emp",tmp.names) %in% grep("_sq",tmp.names))]
        
        ### 
        # Looping through files
        # lapply(1:length(split.current.extents), gl.pred.fun2)#, mc.cores = min(c(1,length(split.current.extents))))
        
        for(i in 1:length(split.current.extents)) {
          # Logic check - 
          # Predict the SDM, or create a dummy SDM
          # Doing this based on AUC of the models
          # And to save time
          if('glm.gaussian' %in% names(mod.list$Models)[tmp.list]) { # Predict
            pred.raster.stack <- crop(raster.stack.bc.xm.pred, split.current.extents[[i]])
            
            df.dat.tmp <- 
              data.frame(GDD = getValues(pred.raster.stack[[grep("GDD",names(pred.raster.stack))[!(grep("GDD",names(pred.raster.stack)) %in% grep("_sq",names(pred.raster.stack)))]]]),
                         Precip = getValues(pred.raster.stack[[grep("Precip",names(pred.raster.stack), ignore.case = TRUE)[!(grep("Precip",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                         Water_Bal = getValues(pred.raster.stack[[grep("Water",names(pred.raster.stack), ignore.case = TRUE)[!(grep("Water",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                         Min_temp = getValues(pred.raster.stack[[grep("Temp",names(pred.raster.stack), ignore.case = TRUE)[!(grep("Temp",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                         GDD_sq = getValues(pred.raster.stack[[grep("GDD",names(pred.raster.stack))[(grep("GDD",names(pred.raster.stack)) %in% grep("_sq",names(pred.raster.stack)))]]]),
                         Precip_sq = getValues(pred.raster.stack[[grep("Precip",names(pred.raster.stack), ignore.case = TRUE)[(grep("Precip",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                         Water_Bal_sq = getValues(pred.raster.stack[[grep("Water",names(pred.raster.stack), ignore.case = TRUE)[(grep("Water",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                         Min_temp_sq = getValues(pred.raster.stack[[grep("Temp",names(pred.raster.stack), ignore.case = TRUE)[(grep("Temp",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]))
            # Creating tmp raster
            gl.pred.tmp <-
              raster(matrix(predict(mod.list$Models[['glm.gaussian']],
                                    df.dat.tmp, type = 'response'),
                            nrow = pred.raster.stack@nrows,
                            ncol = pred.raster.stack@ncols,
                            byrow = TRUE),
                     crs = crs(pred.raster.stack))
            extent(gl.pred.tmp) <- extent(pred.raster.stack)
          } else { # Else create a raster full of NAs based on split current extent
            gl.pred.tmp <- split.current.extents[[i]]
            gl.pred.tmp[!is.na(gl.pred.tmp)] <- NA
          } # End if else statement check whether to forecast random forest
          
          # Writing if the species has more than one tile in their habitat range
          if(length(split.current.extents) %in% 1) {
            # And saving
            writeRaster(gl.pred.tmp,
                        paste0(path.sdms.write.tmp,
                               '/',
                               species_name_save,
                               '_glm.gaussian.tif'),
                        overwrite = TRUE)
            
            # rm(out.raster)
            rm(gl.pred.tmp)
          } else {
            # Else write the raster tile
            writeRaster(gl.pred.tmp,
                        paste0(path.sdms.write.tmp,
                               '/',
                               species_name_save,
                               '_tmp_sdm_tile_',
                               i,
                               ".tif"),
                        overwrite = TRUE)
            rm(gl.pred.tmp)
            
            if(i %in% length(split.current.extents)) {
              # If on the last tile of the habitat range, then mosaic all other tiles
              
              # Getting list of rasters to mosaic
              gdalfile <- 
                list.files(path = path.sdms.write.tmp,
                           pattern = species_name_save,
                           full.names = TRUE) %>% .[grepl('_tmp',.)]
              
              # Place to write the sdm
              dst_dataset = paste0(path.sdms.write.tmp,'/',species_name_save,'_glm.gaussian.tif')
              
              # And mosaicing
              mosaic_rasters(gdalfile, dst_dataset)
            } # End statement to merge
          }
        }
        
        # gl.merge <-
        #   raster(paste0(path.sdms.write.tmp,
        #                 '/',
        #                 species_name_save,
        #                 '_glm.gaussian.tif'))
        # crs(gl.merge) = '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'
        
        # Removing files in the folders
        tmp.files <- list.files(path.sdms.write.tmp,full.names=TRUE, pattern = paste0(species_name_save,'_tmp'))
        file.remove(tmp.files)
        
        tmp.files <- list.files(path.sdms.write.mosaiced.tmp,full.names=TRUE, pattern = species_name_save)
        file.remove(tmp.files)
        
        # And saving outputs
        # if(!is.Raster(models)) {
        #   models <- stack(gl.merge)
        # } else {
        #   models <- stack(models,gl.merge)
        # }
        auc.list['glm.gaussian'] <- mod.list$AUC$glm.gaussian
        threshold.list['glm.gaussian'] <- mod.list$Thresholds$glm.gaussian
        # rm(gl.merge)
      }
      
      
      # Importing models
      models <- stack()
      for(i in tmp.list) {
        models <- stack(models,
                        raster(list.files(path = path.sdms.write.tmp,
                                          full.names = TRUE,
                                          pattern = paste0(species_name_save,'_')) %>%
                                 .[grep(names(threshold.list)[i], .)] %>%
				 .[!grepl('.aux.xml',.)]))
        # auc.weight <- auc.list[[i]]
        # tr.weight <- mod.list$Thresholds[i]
      }
      names(models) <- names(mod.list$Models)[tmp.list]
      # Returning models
      return(models)
      
      
      ## Accuracy for weighting
      # auc.weight <- lapply(tmp.list, function(i) return(auc.list[[i]]))
      # auc.weight <- auc.list[[tmp.list]]
      # Thresholds for weighting
      # tr.weight <- as.numeric(mod.list$Thresholds[tmp.list])
      
      ###
      # Next getting model average based on auc values
      # names(models) <- names(mod.list$Models)[tmp.list]
      # names(auc.weight) <- names(mod.list$Models)[tmp.list]
      # Getting weighted model means
      # List of model aucs
      # auc <- sapply(auc.weight, function(x) x@auc) 
      # Weighting based on accuracy
      # AUC = 0.5 means model is as good as random
      # Weighting based on square of the distance fomr auc = 0.5
      # So more weighting on more accurate models
      # w <- (auc-0.5)^2
      
      
      
      # And weighting the thresholds
      # Getting weighted threshold
      # tr.weight <-
      #  weighted.mean(tr.weight,
      #                w)
      
      # Weighting the models
      # mods.weighted <- weighted.mean(models, w)
      # This is slow, the below returns the same estimates, with the exception of rounding errors on the order of 1e-16
      # w <- w/sum(w)
      
      # Saving time - don't need to weight models if only one is used!
      # if(length(names(models)) %in% 1) {
      #   mods.weighted <- models
      # } else {
      #   mods.weighted <- weight.fun.raster(w, models)
      # }
      ###
      # stacking weighted mod with the individual mdoels
      # models <-
      #   stack(mods.weighted,
      #        models)
      
      # names(models)[1] <- 'weighted'
      
      
      # Creating a list of thresholds
      # threshold.list[['weighted']] <- tr.weight
      # And reorganizing to be the same order as the model list and auc list
      # threshold.list<-threshold.list[c(length(threshold.list),1:(length(threshold.list)-1))]
      
      # Creating list to return
      # return.list <-
      #  list(models,
      #       threshold.list)
      # names(return.list) <-
      #  c("Model.Predictions",
      #    "Thresholds")
      
      ###
      # And returning
      # return(return.list) 
    }
    
    if(length(tmp.list) == 0) {
      return(NA)
    }
  }

###
# And function for updating SDM predictions based on habitat preference and TNC ecoregions
# sdm.hab.fun <-
#   function(pred.rasters,
#            glob_cover.raster,
#            tnc_ecoregions.raster,
#            esh.raster) {
#     ###
#     # Cropping glob cover and tnc rasters
#     glob.cover.tmp <-
#       crop(glob_cover.raster,
#            esh.raster)
#     tnc.ecoregions.tmp <-
#       crop(ecoregions.map,
#            esh.raster)
#     
#     ###
#     # Limiting extent by glob cover raster values
#     glob.cover.values <-
#       getValues(glob.cover.tmp)
#     
#     esh.values <-
#       getValues(esh.raster)
#     
#     # Dropping NAs
#     glob.cover.values <-
#       glob.cover.values[esh.values == 1 &
#                           !is.na(esh.values)]
#     
#     esh.values <-
#       esh.values[esh.values == 1 &
#                    !is.na(esh.values)]
#     
#     # Count by habitat cover
#     tmp.matrix <-
#       matrix(c(glob.cover.values,
#                esh.values),
#              ncol = 2,
#              nrow = length(glob.cover.values),
#              byrow = FALSE)
#     
#     glob.cover.sums <-
#       rowsum(tmp.matrix,
#              group = tmp.matrix[,1])
#     
#     glob.cover.sums[,1] <-
#       glob.cover.sums[,1] / glob.cover.sums[,2]
#     
#     glob.cover.sums <-
#       as.data.frame(glob.cover.sums)
#     # And proportion of cells within each land cover type
#     glob.cover.sums$prop.cover <-
#       glob.cover.sums[,2] /
#       sum(glob.cover.sums[,2])
#     # And dropping land cover types with <= 5% of habitat cover
#     # Unless they are also in the IUCN habitat list
#     glob.cover.sums <-
#       glob.cover.sums[glob.cover.sums$prop.cover >= .025,]
#     
#     # Adding to list
#     # glob.cover.list <-
#     #   unique(c(glob.cover.list,
#     #            glob.cover.sums$V1))
#     
#     # Cropping out land cover types that are not suitable for the species
#     m3 <- pred.rasters$Model.Predictions$weighted
#     m3[m3 > pred.rasters$Thresholds$weighted &
#          !is.na(m3)] <- 1
#     m3[m3 <= pred.rasters$Thresholds$weighted |
#          is.na(m3)] <- NA
#     
#     m3[!(glob.cover.tmp %in% glob.cover.sums$V1) &
#          !is.na(m3)] <- NA
#     m3[m3 == 0] <- NA
#     
#     ###
#     # Repeating for TNC Ecoregions
#     # Getting unique values of habitat covers where the species exists
#     ecoregion.values <-
#       getValues(tnc.ecoregions.tmp)
#     
#     esh.values <-
#       getValues(esh.raster)
#     
#     # Dropping NAs
#     ecoregion.values <-
#       ecoregion.values[esh.values == 1 &
#                          !is.na(esh.values)]
#     
#     esh.values <-
#       esh.values[esh.values == 1 &
#                    !is.na(esh.values)]
#     
#     # Count by habitat cover
#     tmp.matrix <-
#       matrix(c(ecoregion.values,
#                esh.values),
#              ncol = 2,
#              nrow = length(glob.cover.values),
#              byrow = FALSE)
#     
#     ecoregion.sums <-
#       rowsum(tmp.matrix,
#              group = tmp.matrix[,1])
#     
#     ecoregion.sums[,1] <-
#       ecoregion.sums[,1] / ecoregion.sums[,2]
#     
#     ecoregion.sums <-
#       as.data.frame(ecoregion.sums)
#     # And proportion of cells within each land cover type
#     ecoregion.sums$prop.cover <-
#       ecoregion.sums[,2] /
#       sum(ecoregion.sums[,2])
#     # And dropping land cover types with <= 5% of habitat cover
#     # If there are also < 100 cells in that habitat type
#     ecoregion.sums <-
#       ecoregion.sums[ecoregion.sums$prop.cover >= .025,]
#     
#     # Cropping out land cover types that are not suitable for the species
#     m3[!(tnc.ecoregions.tmp %in% ecoregion.sums$V1) &
#          !is.na(m3)] <- NA
#     
#     ###
#     # And returning the map
#     return(m3)
#   }


###
# And function for updating SDM predictions based on habitat preference and TNC ecoregions
sdm.hab.fun <-
  function(pred.rasters,
           esh.raster) {
    
    # Converting from raster, to vector, then back to raster
    # This might not make a huge difference on small rasters
    # But should make a large difference on larger species
    
    # Getting values
    m3 <- getValues(pred.rasters$Model.Predictions$weighted)
    
    # Converting to NAs and 1s
    m3[!is.na(m3) & m3 < pred.rasters$Thresholds$weighted] <- NA
    m3[!is.na(m3)] <- 1
    
    # Converting back to rasters
    m3 <- raster(matrix(m3, byrow = TRUE, nrow = pred.rasters$Model.Predictions$weighted@nrows))
    crs(m3) <- crs(pred.rasters$Model.Predictions$weighted)
    extent(m3) <- extent(pred.rasters$Model.Predictions$weighted)
    
    # And returning the map
    return(m3)
  }

###
# And one more function for clipping habitat based on distance from current esh
sdm.migrate.fun <-
  function(esh.raster,
           hab.clip.raster,
           migration.distance) {
    
    ###
    # Getting raster containing distance from current ESH
    distance.raster <-
      distance(esh.raster)
    
    ###
    # Limiting habitat raster to within a certain migration distance
    hab.clip.raster[distance.raster > migration.distance] <- NA
    
    ###
    # And returning
    return(hab.clip.raster)
  }
# 
# 
# 
# test.mods <-
#   sdm.mod.fun(raster.stack = tmp.stack,
#               df.dat = tmp.df,
#               esh.raster = esh.raster)
# 
# test.preds <-
#   sdm.pred.fun(raster.stack = tmp.stack,
#                df.dat = tmp.df,
#                mod.list = test.mods)
# 
# test.clip.habs <-
#   sdm.hab.fun(pred.rasters = test.preds,
#               glob_cover.raster = glob.cover,
#               tnc_ecoregions.raster = ecoregions.map,
#               esh.raster = esh.raster)
# 
# test.migration <-
#   sdm.migrate.fun(esh.raster = esh.raster,
#                   hab.clip.raster = test.clip.habs,
#                   migration.distance = 15000)






# First function runs models
sdm.mod.fun.eco <-
  function(#raster.stack,
    raster.stack.bc.xm,
    #df.dat,
    df.dat.bc.xm,
    esh.raster,
    realm,
    biome,
    eco) {
    #####
    # Getting data ready to run bioclim, domain, and maxent
    # Extracting coordinates of cells that contain the species
    tmp.coords <-
      xyFromCell(esh.raster,
                 df.dat.bc.xm$cell_num)
    
    df.dat.bc.xm <-
      cbind(df.dat.bc.xm,
            as.data.frame(tmp.coords))
    # 
    # df.dat.bc.xm <-
    #   cbind(df.dat.bc.xm,
    #         as.data.frame(tmp.coords))
    
    # Limiting data frames to only the ecoregions and land cover types
    # The species currently exists in
    eco.dat <- data.frame(biome = getValues(biome), realm = getValues(realm), eco = getValues(eco), esh.values = getValues(esh.raster))
    # Getting summary by biome
    # biome.sum.dat <- eco.dat %>% filter(esh.values %in% 1) %>% filter(!is.na(biome)) %>% dplyr::select(-realm) %>% dplyr::select(-eco) %>% as.matrix() %>% rowsum(.,group = .[,1], na.rm = TRUE) %>%
    #   as.data.frame() %>% mutate(biome = biome/esh.values) %>% mutate(esh.values = esh.values/sum(esh.values,na.rm=TRUE))
    
    biome.sum.dat <-
      tryCatch(eco.dat %>% filter(esh.values %in% 1) %>% filter(!is.na(biome)) %>% dplyr::select(-realm) %>% dplyr::select(-eco) %>% as.matrix() %>% rowsum(.,group = .[,1], na.rm = TRUE) %>%
                 as.data.frame() %>% mutate(biome = biome/esh.values) %>% mutate(esh.values = esh.values/sum(esh.values,na.rm=TRUE)),
               error = function(err) {return(data.frame(biome = 0))}
      ) %>%
      filter(esh.values >= .01)
    # Getting summary by realm
    # realm.sum.dat <- eco.dat %>% filter(esh.values %in% 1) %>% filter(!is.na(realm)) %>% dplyr::select(-biome) %>% dplyr::select(-eco) %>% as.matrix() %>% rowsum(.,group = .[,1], na.rm = TRUE) %>%
    #   as.data.frame() %>% mutate(realm = realm/esh.values) %>% mutate(esh.values = esh.values/sum(esh.values,na.rm=TRUE))
    
    realm.sum.dat <-
      tryCatch(eco.dat %>% filter(esh.values %in% 1) %>% filter(!is.na(realm)) %>% dplyr::select(-biome) %>% dplyr::select(-eco) %>% as.matrix() %>% rowsum(.,group = .[,1], na.rm = TRUE) %>%
                 as.data.frame() %>% mutate(realm = realm/esh.values) %>% mutate(esh.values = esh.values/sum(esh.values,na.rm=TRUE)),
               error = function(err) {return(data.frame(realm = 0))}
      ) %>%
      filter(esh.values >= .01)
    # Getting summary by ecoregion
    # eco.sum.dat <- eco.dat %>% filter(esh.values %in% 1) %>% filter(!is.na(eco)) %>% dplyr::select(-biome) %>% dplyr::select(-realm) %>% as.matrix() %>% rowsum(.,group = .[,1], na.rm = TRUE) %>%
    #   as.data.frame() %>% mutate(eco = eco/esh.values) %>% mutate(esh.values = esh.values/sum(esh.values,na.rm=TRUE))
    
    eco.sum.dat <-
      tryCatch(eco.dat %>% filter(esh.values %in% 1) %>% filter(!is.na(eco)) %>% dplyr::select(-biome) %>% dplyr::select(-realm) %>% as.matrix() %>% rowsum(.,group = .[,1], na.rm = TRUE) %>%
                 as.data.frame() %>% mutate(eco = eco/esh.values) %>% mutate(esh.values = esh.values/sum(esh.values,na.rm=TRUE)),
               error = function(err) {return(data.frame(eco = 0))}
      ) %>%
      filter(esh.values >= .01) %>%
      filter(eco > 1)
    
    # Sense checks to make sure
    # (a) there is some overlap between esh raster and biome and realm rasters
    # (b) at least 25% of species habitat range is on land
    
    
    # This first makes sure there is some overlap between the rasters
    if(!(max(unique(biome.sum.dat$biome)) %in% 0 | max(unique(realm.sum.dat$realm)) %in% 0 | max(unique(eco.sum.dat$eco) %in% 0))) {
      # This second makes sure that at least 25% of the habitat range is on land
      if(sum(biome.sum.dat$esh.values[biome.sum.dat$biome != 0]) >= .75 &
         sum(realm.sum.dat$esh.values[realm.sum.dat$realm != 0]) >= .75) {
        ### Filtering data frames by 
        # Land cover type
        # esh.df <- df.dat[eco.dat$biome %in% biome.sum.dat$biome &
        #                    eco.dat$realm %in% realm.sum.dat$realm &
        #                    eco.dat$eco %in% eco.sum.dat$eco,]
        esh.df.bc.xm <- df.dat.bc.xm[eco.dat$biome %in% biome.sum.dat$biome &
                                       eco.dat$realm %in% realm.sum.dat$realm &
                                       eco.dat$eco %in% eco.sum.dat$eco,]
        
        ###
        # Limiting only to cells that contain the species
        # esh.df <-
        #   esh.df[!is.na(esh.df$esh),]
        
        esh.df.bc.xm <- 
          esh.df.bc.xm[!is.na(esh.df.bc.xm$esh),]
        
        ###
        # And now running models
        # This chunk is for setting up data
        
        ###
        # Bioclim model
        # Bioclim uses only presence values
        # pred_nf <- dropLayer(predictors,'biome')
        # group <- kfold(bradypus, 5)
        # pres_train <- bradypus[group != 1, ]
        # pres_test <- bradypus[group == 1, ]
        # bc <- bioclim(pred_nf, pres_train
        
        # Extracting values in these cells
        # tmp.test <-
        #   extract(esh.raster,
        #           tmp.coords,
        #           method = 'simple')
        
        # Dividng into groups for cross validation
        # Making another exception!
        # Need to make sure there are more than 50 data points in the presence data
        # Doing this here
        # Because this is after clipping land cover and ecoregion
        
        # Creating empty object
        # Doing this for the boolean later on
        backg <- c()
        
        if(nrow(esh.df.bc.xm)>25 & !is.na(minValue(raster.stack.bc.xm$GDD_2010))) { # if statement to only run on species w/ enough cells and that overlap w/ land
          
          # SPlitting into groups for testing and training
          esh.df.bc.xm$group <-
            kfold(esh.df.bc.xm,
                  k = 5)
          
          
          
          ####
          # Splitting into test and training groups
          # Test group is 80% of the data set
          # Training group is 20% of the data set
          # pres_train <-
          #   esh.df[esh.df$group != 1,]
          # pres_test <-
          #   esh.df[esh.df$group == 1,]
          
          pres_train.rf.glm <-
            esh.df.bc.xm[esh.df.bc.xm$group != 1,]
          pres_test.rf.glm <-
            esh.df.bc.xm[esh.df.bc.xm$group == 1,]
          
          ###
          # Repeating for background data
          # For maxent, this is 10,000 background points
          # First just randomly selecting all the points
          # Note that this randomizes the order
          # This might also throw a message saying it did not extract as many points as wanted
          # This is ok
          # backg <- randomPoints(raster.stack, (esh.raster@nrows*esh.raster@ncols - sum(!is.na(getValues(esh.raster)))),tryf=10)
          backg <- randomPoints(raster.stack.bc.xm, min(round((esh.raster@nrows*esh.raster@ncols - sum(!is.na(getValues(esh.raster))))/4),(nrow(pres_train.rf.glm) + nrow(pres_test.rf.glm))*4),tryf=10)
          
          backg <- cbind(backg,1:nrow(backg))
          colnames(backg) = c('lon','lat','cell_num')
          
          # Limiting to data points where esh = 0 
          # Getting ESH values at the backg points
          # Ideally these will all be NAs
          esh.vals.backg <- raster::extract(x = esh.raster,
                                            y = backg[,c('lon','lat')],
                                            method = 'simple')
          # Dropping coordinates that don't contain NAs in the esh raster
          backg <- backg[which(is.na(esh.vals.backg)),]
          
          # Getting ecoregion and biome of the selected points
          back.land <- data.frame(biome.tmp = raster::extract(biome,y = backg[,c('lon','lat')], method = 'simple'),
                                  realm.tmp = raster::extract(realm,y = backg[,c('lon','lat')], method = 'simple'),
                                  eco.tmp = raster::extract(eco, y = backg[,c('lon','lat')], method = 'simple'),
                                  gdd = raster::extract(raster.stack.bc.xm$GDD_2010,y = backg[,c('lon','lat')], method = 'simple'),
                                  temp = raster::extract(raster.stack.bc.xm$Min_Temp_2010,y = backg[,c('lon','lat')], method = 'simple'),
                                  water_bal = raster::extract(raster.stack.bc.xm$Water_2010,y = backg[,c('lon','lat')], method = 'simple'),
                                  precip = raster::extract(raster.stack.bc.xm$Precip_2010,y = backg[,c('lon','lat')], method = 'simple'))
          back.land$cell_num <- backg[,colnames(backg) == 'cell_num']
          
          # And filtering
          # Note that this also gets rid of all non-land areas
          # Because the eco regions and biome data stes we are using are only for terrestrial regions
          backg <- backg[back.land$biome.tmp %in% biome.sum.dat$biome & back.land$realm.tmp %in% realm.sum.dat$realm & back.land$eco.tmp %in% eco.sum.dat$eco,]
          
          
        } # End of if statement
        
        
        
        
        if(!is.matrix(backg)) {
          backg <- matrix(0,0,0,nrow=1)
        }
        
        ### For island-dwelling species, this can result in no out-of-sample predictions
        # Preventing the script from running on these species (will return an error)
        # This means that these species will be excluded from our analysis
        # Or at least from the climate change part of our analysis
        # Using an if statement to do this
        if(nrow(backg) > 25) {
          
          # Uncomment if you want to check to be sure
          # na.check <- extract(tmp.stack,
          #                     backg[,c('lon','lat')],
          #                     method = 'simple')
          # na.check.list <- which(!complete.cases(na.check))
          
          # Randomly selecting 10,000 cells for maxent 
          backg_xm <- 
            sample_n(as.data.frame(backg),
                     size = min(10000,nrow(backg)),
                     replace = FALSE)
          
          # Counters to keep the loop running
          row.check = nrow(backg)
          row.check.1 = nrow(backg)
          if(nrow(backg) < (nrow(pres_train.rf.glm) + nrow(pres_test.rf.glm))) {
            row.check.1 = -1
          }
          # Repeating if sample size is not adequately large
          counter = 1
          # Looping to add cells to the data frame
          # Note that this is used for random forest and glm
          while(nrow(backg) < (nrow(pres_train.rf.glm) + nrow(pres_test.rf.glm)) &
                counter < 5 &
                row.check != row.check.1)  {
            # Random points
            backg.tmp <- randomPoints(raster.stack.bc.xm, min(round((esh.raster@nrows*esh.raster@ncols - sum(!is.na(getValues(esh.raster))))/4),(nrow(pres_train.rf.glm) + nrow(pres_test.rf.glm))*4),tryf=10)
            colnames(backg.tmp) = c('lon','lat')
            
            # Limiting to data points where esh = 0 
            # Getting ESH values at the backg points
            # Ideally these will all be NAs
            esh.vals.backg <- raster::extract(x = esh.raster,
                                              y = backg.tmp,
                                              method = 'simple')
            
            # Dropping coordinates that don't contain NAs in the esh raster
            backg.tmp <- backg.tmp[which(is.na(esh.vals.backg)),]
            
            # Getting ecoregion and biome of the selected points
            back.land.tmp <- data.frame(biome.tmp = raster::extract(biome,y = backg.tmp[,c('lon','lat')],method = 'simple'),
                                        realm.tmp = raster::extract(realm,y = backg.tmp[,c('lon','lat')], method = 'simple'),
                                        eco.tmp = raster::extract(eco, y = backg.tmp[,c('lon','lat')], method = 'simple'))
            
            # And filtering
            # Have to do this a silly way because this is throwing errors
            backg.tmp <- backg.tmp[back.land.tmp$biome.tmp %in% biome.sum.dat$biome & back.land.tmp$realm.tmp %in% realm.sum.dat$realm & back.land.tmp$eco.tmp %in% eco.sum.dat$eco,]
            
            backg.tmp <- cbind(backg.tmp,-1)
            
            # Adding to backg
            backg <- rbind(backg, backg.tmp)
            
            backg <- unique(backg[,1:2])
            backg <- cbind(backg,-1)
            
            # Making sure loop doesn't get stuck
            # You really want to comment this out when running the scripts!
            # Like really really want to!
            # backg <- backg[1:5,]
            
            # Adding to counter
            counter = counter + 1
            row.check.1
          }
          
          # Now limiting to same number of rows as in pres test + pres train
          backg <- 
            sample_n(as.data.frame(backg),
                     size = min((nrow(pres_test.rf.glm) + nrow(pres_train.rf.glm)),nrow(backg)),
                     replace = FALSE)
          
          # And getting the background data for random forest and the glm model
          backg.rf.glm <-
            raster::extract(raster.stack.bc.xm,
                            backg[,c('lon','lat')])
          
          
          # Dividing into groups for testing and training
          group <- kfold(backg.rf.glm, 5)
          backg_train.rf.glm <- backg.rf.glm[group != 1, ]
          backg_test.rf.glm <- backg.rf.glm[group == 1, ]
          
          # Adding point locations
          backg_train.rf.glm <- cbind(backg_train.rf.glm, backg[group != 1,c('lon','lat')])
          backg_test.rf.glm <- cbind(backg_test.rf.glm, backg[group == 1,c('lon','lat')])
          
          backg_train.rf.glm <- as.data.frame(backg_train.rf.glm)
          backg_test.rf.glm <- as.data.frame(backg_test.rf.glm)
          
          ###
          # Testing accuracy + processing time of bioclim
          # n.points = c()
          # t.vector <- c()
          # auc.vector <- c()
          # for(i in seq(from = 5000, to = 150000, by = 5000)) {
          #   # for(y in seq(from = 50, to = 500, by = 50)) {
          #   for(z in 1:10) {
          #     
          #     
          #     pres.train <- sample_n(pres_train.rf.glm,
          #                            i,
          #                            replace = FALSE)
          #     pres.test <- sample_n(pres_test,
          #                             i/5,
          #                             replace = FALSE)
          #     backg.test <- sample_n(backg_test.rf.glm,
          #                            i/5,
          #                            replace = FALSE)
          #     t1 = Sys.time()
          #     bc <- bioclim(raster.stack, 
          #                   pres.train[,c('x','y')])
          #     
          #     t1 = Sys.time() - t1
          #     # Getting auc
          #     auc_bc <-
          #       evaluate(pres_test[,c('x','y')], backg_test.rf.glm[,c('lon','lat')], bc, raster.stack)
          #     
          #     pres.test <- sample_n(pres_test.rf.glm,
          #                           round(i/8,digits=0),
          #                           replace=FALSE)
          #     abs.test <- sample_n(backg_test.rf.glm,
          #                          round(i/8,digits=0),
          #                          replace=FALSE)
          #     
          #     # Getting AUC
          #     auc_bc <-
          #       evaluate(pres.test[,c('x','y')], abs.test[,c('lon','lat')], xm, raster.stack.bc.xm)
          #     
          #     # and data frame
          #     n.points = c(n.points,i)
          #     t.vector <- c(t.vector,t1)
          #     auc.vector <- c(auc.vector,auc_bc@auc)
          #     # n.trees = c(n.trees, y)
          #   }
          #   # }
          #   
          #   
          #   # And checking progress
          #   print(i)
          #   
          # }
          # tmp.frame <- data.frame(n.points = n.points,
          #                         time = t.vector,
          #                         auc = auc.vector)
          
          # And plotting
          # par(mfrow = c(1,1))
          # 
          # tmp1 = ggplot(dat = tmp.frame, aes(x = n.points, y = time, group = n.points)) +
          #   geom_boxplot() +
          #   labs(title = 'points vs time')
          # 
          # tmp2 = ggplot(dat = tmp.frame, aes(x = n.points, y = auc, group = n.points)) +
          #   geom_boxplot() +
          #   labs(title = 'points vs accuracy')
          # 
          # tmp3 = ggplot(dat = tmp.frame, aes(x = time, y = auc, colour = factor(n.points))) +
          #   geom_point() +
          #   labs(title = 'time vs accuracy')
          # 
          # library(cowplot)
          # tmp3 = plot_grid(tmp1,tmp2,tmp3,nrow=3)
          # tmp3
          # 
          # ggplot(dat = tmp.frame, aes(x = n.points, y = auc)) +
          #   geom_point() +
          #   geom_smooth()
          # 
          # ggsave("/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Performance Testing/Aciononyx jubatus Bioclim Model.pdf",
          #        width=10,height=10)
          # write.csv(tmp.frame,
          #           "/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Performance Testing/Aciononyx jubatus Bioclim Model.csv")
          # 
          
          ###
          # Bioclim
          # For performance vs accuracy tradeoff, limiting to 15,000 presence points (e.g. nrow = 10,000)
          # This means that the testing data sets each contain 3,000 data points
          
          if(nrow(pres_train.rf.glm) <= 15000) {
            pres.train.bc <- pres_train.rf.glm#[,c('cell_num','esh','GDD','Precip','Water_Bal','Min_temp','x','y','group')]
            pres.test.bc <- pres_test.rf.glm#[,c('cell_num','esh','GDD','Precip','Water_Bal','Min_temp','x','y','group')]
            backg.test.bc <- backg_test.rf.glm
          }
          if(nrow(pres_train.rf.glm) > 15000) {
            pres.train.bc <- sample_n(pres_train.rf.glm,
                                      15000,
                                      replace = FALSE)
            pres.test.bc <- sample_n(pres_test.rf.glm,
                                     min(3000,nrow(pres_test.rf.glm)),
                                     replace = FALSE)
            backg.test.bc <- sample_n(backg_test.rf.glm,
                                      min(3000,nrow(backg_test.rf.glm)),
                                      replace = FALSE)
          }
          
          bc <- dismo::bioclim(raster.stack.bc.xm[[c('GDD_2010','Min_Temp_2010','Precip_2010','Water_2010')]], 
                               pres.train.bc[,c('x','y')])
          
          # Getting auc
          auc_bc <-
            dismo::evaluate(p = pres.test.bc[,c('x','y')], 
			    a = backg.test.bc[,c('lon','lat')], 
			    model = bc, 
			    x = raster.stack.bc.xm[[c('GDD_2010','Min_Temp_2010','Precip_2010','Water_2010')]])
          # Getting threshold of values within the existing ESH
          # bc.tr <-
          #  threshold(auc_bc,'spec_sens')
    	bc.tr <- threshold(auc_bc)      

          # Testing accuracy
          # If accuracy is close to 1, then model is almost certainly overfitting
          # if(auc_bc@auc >= .99) {
          #   
          #   # Separating into 4 chunks
          #   pres.train.bc$partition <- dismo::kfold(pres.train.bc,4)
          #   
          #     # Additional cross validations
          #     tmp.bc <- dismo::bioclim(raster.stack.bc.xm[[c('GDD_2010','Min_Temp_2010','Precip_2010','Water_2010')]], 
          #                              pres.train.bc[pres.train.bc$partition %in% 1:2,c('x','y')])
          #     
          #     # Getting auc
          #     tmp.auc <-
          #       dismo::evaluate(pres.train.bc[pres.train.bc$partition %in% 3:4,c('x','y')], backg.test.bc[,c('lon','lat')], tmp.bc, raster.stack.bc.xm[[c('GDD_2010','Min_Temp_2010','Precip_2010','Water_2010')]])
          #     # Getting threshold of values within the existing ESH
          #     tmp.tr <-
          #       threshold(auc_bc,'spec_sens')
          #     
          #     # And updating if needed
          #     if(tmp.auc@auc < .99) { # If accuracy is not super high any more, then update values
          #       bc <- tmp.bc
          #       auc_bc <- tmp.auc
          #       bc.tr <- tmp.tr
          #       
          #       rm(tmp.bc)
          #       rm(tmp.auc)
          #       rm(tmp.tr)
          #     } # Else do nothing
          #   } # End check for model overfitting
          
          
          # bc.tmp <-
          #   bc.pred
          # bc.tmp[bc.tmp < bc.tr] <- NA
          # 
          # plot(bc.pred)
          # plot(esh.raster, add = TRUE)
          # plot(bc.tmp, col = 'black', add = TRUE)
          
          ###
          # Domain next
          # dm <- domain(raster.stack, 
          #              pres_train[,c('x','y')])
          # 
          # # Getting AUC
          # auc_dm <-
          #   evaluate(pres_test[,c('x','y')], backg_test, dm, raster.stack)
          # # And getting threshold
          # dm.tr <-
          #   threshold(auc_dm,'spec_sens')
          
          ###
          # Maxent
          # THe commented out code immediately below tests processing time vs accuracy vs number of points
          # Do not run unless specifically testing for this
          
          # 
          # n.points = c()
          # t.vector <- c()
          # auc.vector <- c()
          # for(i in seq(from = 500, to = 15000, by = 500)) {
          #   # for(y in seq(from = 50, to = 500, by = 50)) {
          #     for(z in 1:10) {
          #       
          # 
          #       pres.train <- sample_n(pres_train.rf.glm,
          #                           i/2,
          #                           replace = FALSE)
          #       backg.train <- sample_n(backg_train.rf.glm,
          #                               i/2,
          #                               replace = FALSE)
          #       t1 = Sys.time()
          #       xm <- maxent(x = raster.stack.bc.xm,
          #                    p = pres.train[,c('x','y')],
          #                    a = backg.train[,c('lon','lat')])
          #       
          #       
          #       t1 = Sys.time() - t1
          # 
          #       pres.test <- sample_n(pres_test.rf.glm,
          #                                 round(i/8,digits=0),
          #                                 replace=FALSE)
          #       abs.test <- sample_n(backg_test.rf.glm,
          #                                round(i/8,digits=0),
          #                                replace=FALSE)
          # 
          #       # Getting AUC
          #       auc_xm <-
          #         evaluate(pres.test[,c('x','y')], abs.test[,c('lon','lat')], xm, raster.stack.bc.xm)
          # 
          #       # and data frame
          #       n.points = c(n.points,i)
          #       t.vector <- c(t.vector,t1)
          #       auc.vector <- c(auc.vector,auc_xm@auc)
          #       # n.trees = c(n.trees, y)
          #     }
          #   # }
          # 
          # 
          #   # And checking progress
          #   print(i)
          # 
          # }
          # tmp.frame <- data.frame(n.points = n.points,
          #                         time = t.vector,
          #                         auc = auc.vector)
          # 
          # 
          # # And plotting
          # par(mfrow = c(1,1))
          # 
          # tmp1 = ggplot(dat = tmp.frame, aes(x = n.points, y = time, group = n.points)) +
          #   geom_boxplot() +
          #   labs(title = 'points vs time')
          # 
          # tmp2 = ggplot(dat = tmp.frame, aes(x = n.points, y = auc, group = n.points)) +
          #   geom_boxplot() +
          #   labs(title = 'points vs accuracy')
          # 
          # tmp3 = ggplot(dat = tmp.frame, aes(x = time, y = auc, colour = factor(n.points))) +
          #   geom_point() +
          #   labs(title = 'time vs accuracy')
          # 
          # library(cowplot)
          # tmp3 = plot_grid(tmp1,tmp2,tmp3,nrow=3)
          # tmp3
          # 
          # ggsave("/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Performance Testing/Maxent Model.pdf",
          #        width=10,height=10)
          # 
          # write.csv(tmp.frame,
          #           "/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Performance Testing/Aciononyx jubatus Maxent Model.csv")
          # 
          
          
          
          # Limiting for processing time and performance
          if(nrow(pres_train.rf.glm) <= 7500) {
            pres.train.xm <- pres_train.rf.glm
            backg.train.xm <- backg_train.rf.glm
            pres.test.xm <- pres_test.rf.glm
            backg.test.xm <- backg_test.rf.glm
          }
          if(nrow(pres_train.rf.glm) > 7500) {
            pres.train.xm <- sample_n(pres_train.rf.glm,min(7500,nrow(pres_train.rf.glm)),replace = FALSE)
            backg.train.xm <- sample_n(backg_train.rf.glm,min(7500,nrow(backg_train.rf.glm)),replace=FALSE)
            pres.test.xm <- sample_n(pres_test.rf.glm,min(1875,nrow(pres_test.rf.glm)),replace=FALSE)
            backg.test.xm <- sample_n(backg_test.rf.glm,min(1875,nrow(backg_test.rf.glm)),replace=FALSE)
          }
          # xm <- maxent(x = raster.stack.bc.xm,
          #              p = pres.train.xm[,c('x','y')],
          #              a = backg.train.xm[,c('lon','lat')])
          
          pres.train.xm <-
            pres.train.xm %>%
            dplyr::rename(GDD_2010 = GDD,
                          GDD_2010_sq = GDD_sq,
                          Min_Temp_2010 = Min_temp,
                          Min_Temp_2010_sq = Min_temp_sq,
                          Water_2010 = Water_Bal,
                          Water_2010_sq = Water_Bal_sq,
                          Precip_2010 = Precip,
                          Precip_2010_sq = Precip_sq) %>%
            filter(complete.cases(.))
          
          # names(backg.train.xm) <- gsub("_2010",'',names(backg.train.xm))
          
          backg.train.xm <- 
            backg.train.xm %>%
            filter(complete.cases(.))

          xm <- tryCatch(maxnet(p=c(rep(1, nrow(pres.train.xm)),rep(0, nrow(backg.train.xm))),
                       data=rbind(pres.train.xm[,names(raster.stack.bc.xm)],backg.train.xm[,names(raster.stack.bc.xm)]),
                       f=maxnet.formula(p=c(rep(1, nrow(pres.train.xm)),rep(0, nrow(backg.train.xm))),
                                        data=rbind(pres.train.xm[,names(raster.stack.bc.xm)],backg.train.xm[,names(raster.stack.bc.xm)]),
                                        classes='lpq'),
                       addsamplestobackground = FALSE),
         error = function(e) 'error')
          # Getting AUC
	  if(xm[[1]][1] %in% 'error') {
	  auc_xm <- auc_bc
	  auc_xm@auc <- 1.5
	  auc_xm@cor <- 1.5
	  xm.tr <- bc.tr
	  xm.tr[,1:ncol(xm.tr)] <- 1.5
	  } else {
	   auc_xm <- dismo::evaluate(p = pres.test.xm[,c('x','y')], 
                                    a = backg.test.xm[,c('lon','lat')], 
                                    model = xm, 
                                    x = raster.stack.bc.xm, type='cloglog')
          # auc_xm <-
          #   evaluate(pres.test.xm[,c('x','y')], backg.test.xm[,c('lon','lat')], xm, raster.stack.bc.xm)
          # And getting threshold
          xm.tr <-
            threshold(auc_xm)
	  }
          
          # Checking on accuracy
          # If accuracy is close to 1
          # Then need to rerun model using less data
          # if(auc_xm@auc >=.99) {
          #   pres.train.xm$partition <- dismo::kfold(pres.train.xm,4)
          #   backg.train.xm$partition <- dismo::kfold(backg.train.xm,4)
          #   tmp.xm <- maxnet(p=c(rep(1, nrow(pres.train.xm %>% filter(partition %in% 1:2))),rep(0, nrow(backg.train.xm %>% filter(partition %in% 1:2)))), 
          #                data=rbind(pres.train.xm[pres.train.xm$partition %in% 1:2,names(raster.stack.bc.xm)],backg.train.xm[backg.train.xm$partition %in% 1:2,names(raster.stack.bc.xm)]), 
          #                f=maxnet.formula(p=c(rep(1, nrow(pres.train.xm %>% filter(partition %in% 1:2))),rep(0, nrow(backg.train.xm %>% filter(partition %in% 1:2)))), 
          #                                 data=rbind(pres.train.xm[pres.train.xm$partition %in% 1:2,names(raster.stack.bc.xm)],backg.train.xm[backg.train.xm$partition %in% 1:2,names(raster.stack.bc.xm)]),
          #                                 classes='lpq'), 
          #                addsamplestobackground = FALSE)
          #   
          #   # Getting AUC
          #   tmp.auc_xm <- dismo::evaluate(pres.test.xm[,c('x','y')], backg.test.xm[,c('lon','lat')], tmp.xm, raster.stack.bc.xm, type='cloglog')
          #   # auc_xm <-
          #   #   evaluate(pres.test.xm[,c('x','y')], backg.test.xm[,c('lon','lat')], xm, raster.stack.bc.xm)
          #   # And getting threshold
          #   tmp.xm.tr <-
          #     threshold(auc_xm,'spec_sens')
          #   
          #   if(tmp.auc_xm@auc < .99) {
          #     xm <- tmp.xm
          #     auc_xm <- tmp.auc_xm
          #     xm.tr <- tmp.xm.tr
          #     
          #     rm(tmp.xm)
          #     rm(tmp.auc_xm)
          #     rm(tmp.xm.tr)
          #   }
          # } # End accuracy check for max ent
          
          # rm(raster.stack.tmp)
          
          ###
          # Random Forest
          # Setting up training data
          # Converting NAs to 0s in the esh data
          backg_train.rf.glm$esh <- 0
          backg_test.rf.glm$esh <- 0
          
          # Binding data inside and outside ESH
          # dplyr::select to reorganize column orders
          train <- rbind(dplyr::select(pres_train.rf.glm,esh,GDD,Precip,Water_Bal,Min_temp,GDD_sq,Precip_sq,Water_Bal_sq,Min_temp_sq,x,y),
                         dplyr::select(backg_train.rf.glm,
                                       esh,
                                       GDD = GDD_2010,
                                       Precip = Precip_2010,
                                       Water_Bal = Water_2010,
                                       Min_temp = Min_Temp_2010,
                                       GDD_sq = GDD_2010_sq,
                                       Precip_sq = Precip_2010_sq,
                                       Water_Bal_sq = Water_2010_sq,
                                       Min_temp_sq = Min_Temp_2010_sq,
                                       x = lon,
                                       y = lat))
          test <- rbind(dplyr::select(pres_test.rf.glm,esh,GDD,Precip,Water_Bal,Min_temp,GDD_sq,Precip_sq,Water_Bal_sq,Min_temp_sq,x,y),
                        dplyr::select(backg_test.rf.glm,
                                      esh,
                                      GDD = GDD_2010,
                                      Precip = Precip_2010,
                                      Water_Bal = Water_2010,
                                      Min_temp = Min_Temp_2010,
                                      GDD_sq = GDD_2010_sq,
                                      Precip_sq = Precip_2010_sq,
                                      Water_Bal_sq = Water_2010_sq,
                                      Min_temp_sq = Min_Temp_2010_sq,
                                      x = lon,
                                      y = lat))
          
          # Getting complete cases
          # Models throw an error if we don't do this
          train <-
            train[complete.cases(train),]
          test <- test[complete.cases(test),]
          
          # Running the model
          model = esh ~ GDD + Precip + Water_Bal + Min_temp
          
          ###
          # Chunk of code below tests performance
          # 
          # n.points = c()
          # t.vector <- c()
          # auc.vector <- c()
          # n.trees = c()
          # for(i in seq(from = 500, to = 15000, by = 500)) {
          #   for(y in seq(from = 50, to = 500, by = 50)) {
          #     for(z in 1:10) {
          #     
          #     
          #     tmp.train <- sample_n(train,
          #                            i,
          #                            replace = FALSE)
          #     tmp.test <- sample_n(test,
          #                             i/4,
          #                             replace = FALSE)
          #     t1 = Sys.time()
          #     # Getting AUC
          #     rf <- randomForest(model,
          #                        data = tmp.train[,c('esh','GDD','Precip','Water_Bal','Min_temp')],
          #                        ntree = y)
          #     t1 = Sys.time() - t1
          #     # Getting accuracy
          #     tmp.test.pres <- sample_n(pres_test.rf.glm,
          #                              round(i/8),
          #                              replace = FALSE)
          #     tmp.test.backg <- sample_n(backg_test.rf.glm,
          #                               round(i/8),
          #                               replace = FALSE)
          #     
          #     tmp.test.backg <- dplyr::rename(tmp.test.backg,
          #                                    GDD = GDD_2010,
          #                                    Min_temp = Min_Temp_2010,
          #                                    Precip = Precip_2010,
          #                                    Water_Bal = Water_2010,
          #                                    GDD_sq = GDD_2010_sq,
          #                                    Min_temp_sq = Min_Temp_2010_sq,
          #                                    Precip_sq = Precip_2010_sq,
          #                                    Water_Bal_sq = Water_2010_sq,
          #                                    x = lon,
          #                                    y = lat)
          #     
          #     auc_rf <-
          #       evaluate(p = tmp.test.pres[complete.cases(tmp.test.pres),],
          #                a = tmp.test.backg[complete.cases(tmp.test.backg),], 
          #                model = rf)
          #     
          #     
          #     
          #     # and data frame
          #     n.points = c(n.points,i)
          #     t.vector <- c(t.vector,t1)
          #     auc.vector <- c(auc.vector,auc_rf@auc)
          #     n.trees = c(n.trees, y)
          #     }
          #   }
          #   
          #   
          #   # And checking progress
          #   print(i)
          #   
          # }
          # tmp.frame <- data.frame(n.points = n.points,
          #                         time = t.vector,
          #                         auc = auc.vector,
          #                         n.trees = n.trees)
          # 
          # 
          # # And plotting
          # par(mfrow = c(1,1))
          # 
          # tmp1 = ggplot(dat = tmp.frame, aes(x = n.points, y = time, group = n.points)) +
          #   geom_boxplot() +
          #   facet_wrap(.~n.trees) +
          #   labs(title = 'points vs time facet n.trees')
          # 
          # tmp2 = ggplot(dat = tmp.frame, aes(x = n.points, y = auc, group = n.points)) +
          #   geom_boxplot() +
          #   facet_wrap(.~n.trees) +
          #   labs(title = 'points vs accuracy facet n.trees')
          #   
          #   tmp3 = ggplot(dat = tmp.frame, aes(x = n.trees, y = time, group = n.trees)) +
          #     geom_boxplot() +
          #     facet_wrap(.~n.points) +
          #     labs(title = 'n.trees vs time facet points')
          #   
          #   tmp4 = ggplot(dat = tmp.frame, aes(x = n.trees, y = auc, group = n.trees)) +
          #     geom_boxplot() +
          #     facet_wrap(.~n.points) +
          #     labs(title = 'n.trees vs accuracy facet points')
          # 
          # library(cowplot)
          # tmp5 = plot_grid(tmp1,tmp2,tmp3,tmp4,nrow=2,ncol=2)
          # tmp5
          # 
          # ggsave("/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Performance Testing/Aciononyx jubatus Random Forest Model.pdf",
          #        width=20,height=10)
          # 
          # tmp1 = ggplot(dat = tmp.frame, aes(x = time, y = auc, group = n.trees)) +
          #   geom_boxplot() +
          #   facet_wrap(.~n.points)+
          #   labs(title = 'time vs accuracy facet points')
          # 
          # tmp2 = ggplot(dat = tmp.frame, aes(x = time, y = auc, group = n.points)) +
          #   geom_boxplot() +
          #   facet_wrap(.~n.trees)+
          #   labs(title = 'time vs accuracy facet n.trees')
          # 
          # tmp3 = plot_grid(tmp1,tmp2,nrow=2)
          # tmp3
          # 
          # 
          # 
          # write.csv(tmp.frame,
          #           "/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Performance Testing/Aciononyx jubatus Random Forest Model.csv")
          # 
          # 
          
          
          
          ###
          # And for performance/processing time tradeoff of random forest
          # Randomly selecting 5000 absence and presence points
          # Creating dummy variables
          if(nrow(train)<=10000) {
            train.rf <- train
            test.rf <- test
            pres_test.rf <- pres_test.rf.glm
            backg_test.rf <- backg_test.rf.glm
            
            backg_test.rf <- dplyr::rename(backg_test.rf,
                                           GDD = GDD_2010,
                                           Min_temp = Min_Temp_2010,
                                           Precip = Precip_2010,
                                           Water_Bal = Water_2010,
                                           GDD_sq = GDD_2010_sq,
                                           Min_temp_sq = Min_Temp_2010_sq,
                                           Precip_sq = Precip_2010_sq,
                                           Water_Bal_sq = Water_2010_sq,
                                           x = lon,
                                           y = lat)
          }
          if(nrow(train)>10000){
            train.rf <- rbind(sample_n(train[train$esh == 0,],min(5000,sum(train$esh==0)),replace = FALSE),
                              sample_n(train[train$esh == 1,],min(5000,sum(train$esh==1)),replace = FALSE))
            
            test.rf <- rbind(sample_n(test[test$esh==0,],min(1250,test$esh==0),replace = FALSE),
                             sample_n(test[test$esh==1,],min(1250,test$esh==1),replace = FALSE))
            
            pres_test.rf <- sample_n(pres_test.rf.glm,
                                     min(1250,nrow(pres_test.rf.glm)),
                                     replace = FALSE)
            backg_test.rf <- sample_n(backg_test.rf.glm,
                                      min(1250,nrow(backg_test.rf.glm)),
                                      replace = FALSE)
            
            backg_test.rf <- dplyr::rename(backg_test.rf,
                                           GDD = GDD_2010,
                                           Min_temp = Min_Temp_2010,
                                           Precip = Precip_2010,
                                           Water_Bal = Water_2010,
                                           GDD_sq = GDD_2010_sq,
                                           Min_temp_sq = Min_Temp_2010_sq,
                                           Precip_sq = Precip_2010_sq,
                                           Water_Bal_sq = Water_2010_sq,
                                           x = lon,
                                           y = lat)
          }
          
          
          
          # This will return a warning about the output having < 5 unique responses
          # This is ok, seeing as our response variable is binary.
          rf <- randomForest(model,
                             data = train.rf[,c('esh','GDD','Precip','Water_Bal','Min_temp')],
                             ntree = 500,
                             type = 'classification')
          
          # Changing names to test the model
          backg_test.rf.glm <- dplyr::rename(backg_test.rf.glm,
                                             GDD = GDD_2010,
                                             Min_temp = Min_Temp_2010,
                                             Precip = Precip_2010,
                                             Water_Bal = Water_2010,
                                             GDD_sq = GDD_2010_sq,
                                             Min_temp_sq = Min_Temp_2010_sq,
                                             Precip_sq = Precip_2010_sq,
                                             Water_Bal_sq = Water_2010_sq,
                                             x = lon,
                                             y = lat)
          
          
          # Now gettin auc for the model
          auc_rf <-
            dismo::evaluate(p = pres_test.rf[complete.cases(pres_test.rf),],
                            a = backg_test.rf[complete.cases(backg_test.rf),], 
                            model = rf)
          # Getting threshold
          rf.tr <-
            threshold(auc_rf)
          
          # Checking for overfitting
          # if(auc_rf@auc >= .99) {
          #   # Adding partition
          #   train.rf$partition <- dismo::kfold(train.rf, 4)
          #   
          #   # This will return a warning about the output having < 5 unique responses
          #   # This is ok, seeing as our response variable is binary.
          #   tmp.rf <- randomForest(model,
          #                      data = train.rf[train.rf$partition %in% 1:2,c('esh','GDD','Precip','Water_Bal','Min_temp')],
          #                      ntree = 500,
          #                      type = 'classification')
          #   
          #   # Changing names to test the model
          #   # backg_test.rf.glm <- dplyr::rename(backg_test.rf.glm,
          #   #                                    GDD = GDD_2010,
          #   #                                    Min_temp = Min_Temp_2010,
          #   #                                    Precip = Precip_2010,
          #   #                                    Water_Bal = Water_2010,
          #   #                                    GDD_sq = GDD_2010_sq,
          #   #                                    Min_temp_sq = Min_Temp_2010_sq,
          #   #                                    Precip_sq = Precip_2010_sq,
          #   #                                    Water_Bal_sq = Water_2010_sq,
          #   #                                    x = lon,
          #   #                                    y = lat)
          #   
          #   
          #   # Now gettin auc for the model
          #   tmp.auc_rf <-
          #     dismo::evaluate(p = pres_test.rf[complete.cases(pres_test.rf),],
          #                     a = backg_test.rf[complete.cases(backg_test.rf),], 
          #                     model = tmp.rf)
          #   # Getting threshold
          #   tmp.rf.tr <-
          #     threshold(auc_rf,'spec_sens')
          #   
          #   # Updating outputs if models are overfit
          #   if(tmp.auc_rf@auc < .99) {
          #     rf <- tmp.rf
          #     auc_rf <- tmp.auc_rf
          #     rf.tr <- tmp.rf.tr
          #     
          #     rm(rf)
          #   }
          # } # End check for overfitting in random forest model

          
          
          
          
          ###
          # GLM Model
          # THe commented out code immediately below tests processing time vs accuracy vs number of points
          # Do not run unless specifically testing for this
          
          
          # n.points = c()
          # t.vector <- c()
          # auc.vector <- c()
          # # n.trees <- c()
          # for(i in seq(from = 5000, to = 100000, by = 5000)) {
          #   # for(y in seq(from = 50, to = 500, by = 50)) {
          #     for(z in 1:10) {
          #       tmp.dat <- sample_n(train,
          #                           i,
          #                           replace = FALSE)
          #       t1 = Sys.time()
          #     
          #       gm.gaus <-
          #         glm(esh ~ GDD + Precip + Water_Bal + Min_temp + GDD_sq + Precip_sq + Water_Bal_sq + Min_temp_sq,
          #             family = gaussian(link = "identity"),
          #             data=tmp.dat)
          #       t1 = Sys.time() - t1
          # 
          #       tmp.pres.test <- sample_n(pres_test.rf.glm,
          #                                 round(i/8,digits=0),
          #                                 replace=FALSE)
          #       tmp.abs.test <- sample_n(backg_test.rf.glm,
          #                                round(i/8,digits=0),
          #                                replace=FALSE)
          #       tmp.abs.test <- dplyr::select(tmp.abs.test,
          #                                     esh,
          #                                     GDD=GDD_2010,
          #                                     Precip = Precip_2010,
          #                                     Water_Bal = Water_2010,
          #                                     Min_temp = Min_Temp_2010,
          #                                     GDD_sq=GDD_2010_sq,
          #                                     Precip_sq = Precip_2010_sq,
          #                                     Water_Bal_sq = Water_2010_sq,
          #                                     Min_temp_sq = Min_Temp_2010_sq,
          #                                     x = lon,
          #                                     y=lat)
          # 
          #       auc_gm <-
          #         evaluate(p = tmp.pres.test[complete.cases(tmp.pres.test),],
          #                  a = tmp.abs.test[complete.cases(tmp.abs.test),],
          #                  model = gm.gaus)
          # 
          #       # and data frame
          #       n.points = c(n.points,i)
          #       t.vector <- c(t.vector,t1)
          #       auc.vector <- c(auc.vector,auc_gm@auc)
          #       # n.trees = c(n.trees, y)
          #     }
          #   # }
          # 
          # 
          #   # And checking progress
          #   print(i)
          # 
          # }
          # tmp.frame <- data.frame(n.points = n.points,
          #                         time = t.vector,
          #                         auc = auc.vector)
          # 
          # write.csv(tmp.frame,
          #           "/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Performance Testing/Aciononyx jubatus GLM Model.csv")
          # 
          # 
          # # And plotting
          # par(mfrow = c(1,1))
          # 
          # tmp1 = ggplot(dat = tmp.frame, aes(x = n.points, y = time, group = n.points)) +
          #   geom_boxplot() +
          #   labs(title = 'points vs time')
          # 
          # tmp2 = ggplot(dat = tmp.frame, aes(x = n.points, y = auc, group = n.points)) +
          #   geom_boxplot() +
          #   labs(title = 'points vs accuracy')
          # 
          # tmp3 = ggplot(dat = tmp.frame, aes(x = time, y = auc, colour = factor(n.points))) +
          #   geom_point() +
          #   labs(title = 'time vs accuracy')
          # 
          # library(cowplot)
          # tmp3 = plot_grid(tmp1,tmp2,tmp3,nrow=3)
          # tmp3
          # 
          # ggsave("/Users/maclark/Desktop/Multiple Drivers of Biodiversity/Performance Testing/Aciononyx jubatus GLM Model Big Data.pdf",
          #        width=10,height=15)
          # 
          
          
          # Limiting data set for processing and performance
          if(nrow(train)<=150000) {
            train.glm <- train
            pres.test.glm <- pres_test.rf.glm
            backg.test.glm <- backg_test.rf.glm
          }
          if(nrow(train) > 150000) {
            train.glm <- rbind(sample_n(train[train$esh==0,],min(75000,nrow(train[train$esh==0,])),replace=FALSE),
                               sample_n(train[train$esh==1,],min(75000,nrow(train[train$esh==1,])),replace=FALSE))
            pres.test.glm <- sample_n(pres_test.rf.glm,min(18750,nrow(pres_test.rf.glm)),replace=FALSE)
            backg.test.glm <- sample_n(backg_test.rf.glm,min(18750,nrow(backg_test.rf.glm)),replace=FALSE)
          }
          
          # gm.gaus <-
          #   glm(esh ~ GDD + Precip + Water_Bal + Min_temp + GDD_sq + Precip_sq + Water_Bal_sq + Min_temp_sq,
          #       family = gaussian(link = "identity"),
          #       data=train.glm)
          
          gm.gaus <-
            glm(esh ~ GDD + Precip + Water_Bal + Min_temp + GDD_sq + Precip_sq + Water_Bal_sq + Min_temp_sq,             
                family = binomial(link = "logit"), 
                data=train.glm)
          # Testing model evaluation
          auc_gm <-
            dismo::evaluate(p = pres.test.glm, 
                            a = backg.test.glm,  
                            model = gm.gaus, 
                            type = 'response')
          # Getting threshold
          gm.tr <-
            threshold(auc_gm)
          
          # Testing accuracy
          # If accuracy high - fit with smaller subset of data to reduce chance of overfitting
          # if(auc_gm@auc >= .99) {
          #   train.glm$partition <- dismo::kfold(train.glm, 4)
          #   tmp.gm.gaus <-
          #     glm(esh ~ GDD + Precip + Water_Bal + Min_temp + GDD_sq + Precip_sq + Water_Bal_sq + Min_temp_sq,             
          #         family = binomial(link = "logit"), 
          #         data=train.glm %>% filter(partition %in% 1:2))
          #   # Testing model evaluation
          #   tmp.auc_gm <-
          #     dismo::evaluate(pres.test.glm, 
          #                     backg.test.glm,  
          #                     tmp.gm.gaus, 
          #                     type = 'response')
          #   # Getting threshold
          #   tmp.gm.tr <-
          #     threshold(auc_gm,'spec_sens')
          #   
          #   # Updating if accuracy is below .99
          #   if(tmp.auc_gm@auc < .99) {
          #     gm.gaus <- tmp.gm.gaus
          #     auc_gm <- tmp.auc_gm
          #     gm.tr <- tmp.gm.tr
          #     
          #     rm(gm.gaus)
          #   }
          # }
          
          ###
          # And list of things to return
          mod.list <-
            list(bc, xm, rf, gm.gaus)
          names(mod.list) <- c('bioclim','maxent','random.forest','glm.gaussian')
          
          auc.list <-
            list(auc_bc, auc_xm, auc_rf, auc_gm)
          names(auc.list) <- c('bioclim','maxent','random.forest','glm.gaussian')
          
          tr.list <-
            list(bc.tr, xm.tr, rf.tr, gm.tr)
          names(tr.list) <- c('bioclim','maxent','random.forest','glm.gaussian')
          
          return.list <-
            list(mod.list,
                 auc.list,
                 tr.list)
          
          names(return.list) <-
            c('Models','AUC','Thresholds')
          
          # And returning this list
          return(return.list) 
        } # End of if statement preventing code from running on species with small out-of-sample data sets
      } # End of if statement preventing code from running on species with no range overlap between biome and realm  
    }
  }








bioclim.pred.fun2 <- 
  function(index.list) {
    
    # for(i in index.list) {
    pred.raster.stack <- crop(raster.stack.bc.xm.pred[[c('GDD_2010','Min_Temp_2010','Precip_2010','Water_2010')]], split.current.extents[[index.list]])
    
    bc.pred.tmp <- predict(pred.raster.stack, 
                           test.mods$Models[['bioclim']], 
                           ext = extent(pred.raster.stack))
    writeRaster(bc.pred.tmp,
                paste0(path.sdms.write.tmp,
                       '/',
                       df$species_name[k],
                       '_tmp_sdm_tile_',
                       index.list,
                       ".tif"),
                overwrite = TRUE)
    # }
  }


bioclim.pred.fun.loop <- 
  function(index.list,
           raster.stack.bc.xm.pred,
           split.current.extents,
           test.mods,
           species_name_save) {
    
    # for(i in index.list) {
    pred.raster.stack <- crop(raster.stack.bc.xm.pred[[c('GDD_2010','Min_Temp_2010','Precip_2010','Water_2010')]], split.current.extents[[index.list]])
    
    bc.pred.tmp <- predict(pred.raster.stack, 
                           test.mods$Models[['bioclim']], 
                           ext = extent(pred.raster.stack))
    writeRaster(bc.pred.tmp,
                paste0(path.sdms.write.tmp,
                       '/',
                       species_name_save,
                       '_tmp_sdm_tile_',
                       index.list,
                       ".tif"),
                overwrite = TRUE)
    # }
  }

bioclim.pred.fun <- 
  function(index.list,
           extent.raster.list,
           predict.raster.stack,
           model.name) {
    
    for(i in index.list) {
      pred.raster.stack <- crop(raster.stack.bc.xm.pred[[c('GDD_2010','Min_Temp_2010','Precip_2010','Water_2010')]], extent.raster.list[[i]])
      
      bc.pred.tmp <- predict(pred.raster.stack, 
                             test.mods$Models[[model.name]], 
                             ext = extent(pred.raster.stack))
      writeRaster(bc.pred.tmp,
                  paste0(path.sdms.write.tmp,
                         '/',
                         df$species_name[k],
                         '_tmp_sdm_tile_',
                         i,
                         ".tif"),
                  overwrite = TRUE)
    }
  }


# Function 2 to work with mclapply
maxent.pred.fun2 <- 
  function(index.num, raster.stack.bc.xm.pred, split.current.extents, test.mods, species_name_save) {
    
    # for(i in index.list) {
    pred.raster.stack <- crop(raster.stack.bc.xm.pred, split.current.extents[[index.num]])
    
    # setTimeLimit(elapse = 500,
    #              transient = TRUE)
    
    xm.pred.tmp <- predict(pred.raster.stack, 
                           test.mods$Models[['maxent']], 
                           ext = extent(pred.raster.stack),
                           type = 'cloglog')
    
    writeRaster(xm.pred.tmp,
                paste0(path.sdms.write.tmp,
                       '/',
                       species_name_save,
                       '_tmp_sdm_tile_',
                       index.num,
                       ".tif"),
                overwrite = TRUE)
    
    Sys.sleep(10)
    
    removeTmpFiles()
    gc()
    xlcFreeMemory()
    # }
  }

xm.lapply <-
  function(split.current.extents, test.mods, raster.stack.bc.xm.pred, species_name_save) {
    
    loop.list <- 1:length(split.current.extents)
    for(raster.chunk in loop.list) {
      pred.raster.stack <- crop(raster.stack.bc.xm.pred, split.current.extents[[raster.chunk]])  
      
      
      
      xm.pred.tmp <- predict(pred.raster.stack, 
                             test.mods$Models[['maxent']], 
                             ext = extent(pred.raster.stack),
                             type = 'cloglog')
      
      writeRaster(xm.pred.tmp,
                  paste0(path.sdms.write.tmp,
                         '/',
                         species_name_save,
                         '_tmp_sdm_tile_',
                         raster.chunk,
                         ".tif"),
                  overwrite = TRUE)
      
      Sys.sleep(3)
    }
    
    
    Sys.sleep(10)
    
    removeTmpFiles()
    gc()
    xlcFreeMemory()
    
  }


rf.pred.fun.loop <-
  function(split.current.extents,
           raster.stack.bc.xm.pred,
           test.mods,
           index.list,
           species_name_save) {
    
    pred.raster.stack <- crop(raster.stack.bc.xm.pred, split.current.extents[[index.list]])
    
    
    
    df.dat.tmp <- 
      data.frame(GDD = getValues(pred.raster.stack[[grep("GDD",names(pred.raster.stack))[!(grep("GDD",names(pred.raster.stack)) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Precip = getValues(pred.raster.stack[[grep("Precip",names(pred.raster.stack), ignore.case = TRUE)[!(grep("Precip",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Water_Bal = getValues(pred.raster.stack[[grep("Water",names(pred.raster.stack), ignore.case = TRUE)[!(grep("Water",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Min_temp = getValues(pred.raster.stack[[grep("Temp",names(pred.raster.stack), ignore.case = TRUE)[!(grep("Temp",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 GDD_sq = getValues(pred.raster.stack[[grep("GDD",names(pred.raster.stack))[(grep("GDD",names(pred.raster.stack)) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Precip_sq = getValues(pred.raster.stack[[grep("Precip",names(pred.raster.stack), ignore.case = TRUE)[(grep("Precip",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Water_Bal_sq = getValues(pred.raster.stack[[grep("Water",names(pred.raster.stack), ignore.case = TRUE)[(grep("Water",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Min_temp_sq = getValues(pred.raster.stack[[grep("Temp",names(pred.raster.stack), ignore.case = TRUE)[(grep("Temp",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]))
    # Creating tmp raster
    rf.pred.tmp <-
      raster(matrix(predict(test.mods$Models[['random.forest']],
                            df.dat.tmp[,c('GDD','Precip','Water_Bal','Min_temp')]),
                    nrow = pred.raster.stack@nrows,
                    ncol = pred.raster.stack@ncols,
                    byrow = TRUE),
             crs = crs(pred.raster.stack))
    extent(rf.pred.tmp) <- extent(pred.raster.stack)
    # Saving tmp raster
    writeRaster(rf.pred.tmp,
                paste0(path.sdms.write.tmp,
                       '/',
                       species_name_save,
                       '_tmp_sdm_tile_',
                       index.list,
                       ".tif"),
                overwrite = TRUE)
    
  }


gl.pred.fun.loop <-
  function(split.current.extents,
           raster.stack.bc.xm.pred,
           test.mods,
           index.list,
           species_name_save) {
    
    pred.raster.stack <- crop(raster.stack.bc.xm.pred, split.current.extents[[index.list]])
    
    df.dat.tmp <- 
      data.frame(GDD = getValues(pred.raster.stack[[grep("GDD",names(pred.raster.stack))[!(grep("GDD",names(pred.raster.stack)) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Precip = getValues(pred.raster.stack[[grep("Precip",names(pred.raster.stack), ignore.case = TRUE)[!(grep("Precip",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Water_Bal = getValues(pred.raster.stack[[grep("Water",names(pred.raster.stack), ignore.case = TRUE)[!(grep("Water",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Min_temp = getValues(pred.raster.stack[[grep("Temp",names(pred.raster.stack), ignore.case = TRUE)[!(grep("Temp",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 GDD_sq = getValues(pred.raster.stack[[grep("GDD",names(pred.raster.stack))[(grep("GDD",names(pred.raster.stack)) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Precip_sq = getValues(pred.raster.stack[[grep("Precip",names(pred.raster.stack), ignore.case = TRUE)[(grep("Precip",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Water_Bal_sq = getValues(pred.raster.stack[[grep("Water",names(pred.raster.stack), ignore.case = TRUE)[(grep("Water",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Min_temp_sq = getValues(pred.raster.stack[[grep("Temp",names(pred.raster.stack), ignore.case = TRUE)[(grep("Temp",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]))
    # Creating tmp raster
    gl.pred.tmp <-
      raster(matrix(predict(test.mods$Models[['glm.gaussian']],
                            df.dat.tmp, type = 'response'),
                    nrow = pred.raster.stack@nrows,
                    ncol = pred.raster.stack@ncols,
                    byrow = TRUE),
             crs = crs(pred.raster.stack))
    extent(gl.pred.tmp) <- extent(pred.raster.stack)
    # Saving tmp raster
    writeRaster(gl.pred.tmp,
                paste0(path.sdms.write.tmp,
                       '/',
                       species_name_save,
                       '_tmp_sdm_tile_',
                       index.list,
                       ".tif"),
                overwrite = TRUE)
    
  }


rf.pred.fun2 <- 
  function(index.list) {
    pred.raster.stack <- crop(raster.stack.bc.xm.pred, split.current.extents[[index.list]])
    
    
    
    df.dat.tmp <- 
      data.frame(GDD = getValues(pred.raster.stack[[grep("GDD",names(pred.raster.stack))[!(grep("GDD",names(pred.raster.stack)) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Precip = getValues(pred.raster.stack[[grep("Precip",names(pred.raster.stack), ignore.case = TRUE)[!(grep("Precip",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Water_Bal = getValues(pred.raster.stack[[grep("Water",names(pred.raster.stack), ignore.case = TRUE)[!(grep("Water",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Min_temp = getValues(pred.raster.stack[[grep("Temp",names(pred.raster.stack), ignore.case = TRUE)[!(grep("Temp",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 GDD_sq = getValues(pred.raster.stack[[grep("GDD",names(pred.raster.stack))[(grep("GDD",names(pred.raster.stack)) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Precip_sq = getValues(pred.raster.stack[[grep("Precip",names(pred.raster.stack), ignore.case = TRUE)[(grep("Precip",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Water_Bal_sq = getValues(pred.raster.stack[[grep("Water",names(pred.raster.stack), ignore.case = TRUE)[(grep("Water",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Min_temp_sq = getValues(pred.raster.stack[[grep("Temp",names(pred.raster.stack), ignore.case = TRUE)[(grep("Temp",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]))
    # Creating tmp raster
    rf.pred.tmp <-
      raster(matrix(predict(test.mods$Models[['random.forest']],
                            df.dat.tmp[,c('GDD','Precip','Water_Bal','Min_temp')]),
                    nrow = pred.raster.stack@nrows,
                    ncol = pred.raster.stack@ncols,
                    byrow = TRUE),
             crs = crs(pred.raster.stack))
    extent(rf.pred.tmp) <- extent(pred.raster.stack)
    # Saving tmp raster
    writeRaster(rf.pred.tmp,
                paste0(path.sdms.write.tmp,
                       '/',
                       df$species_name[k],
                       '_tmp_sdm_tile_',
                       index.list,
                       ".tif"),
                overwrite = TRUE)
  }




gl.pred.fun2 <- 
  function(index.list) {
    pred.raster.stack <- crop(raster.stack.bc.xm.pred, split.current.extents[[index.list]])
    
    df.dat.tmp <- 
      data.frame(GDD = getValues(pred.raster.stack[[grep("GDD",names(pred.raster.stack))[!(grep("GDD",names(pred.raster.stack)) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Precip = getValues(pred.raster.stack[[grep("Precip",names(pred.raster.stack), ignore.case = TRUE)[!(grep("Precip",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Water_Bal = getValues(pred.raster.stack[[grep("Water",names(pred.raster.stack), ignore.case = TRUE)[!(grep("Water",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Min_temp = getValues(pred.raster.stack[[grep("Temp",names(pred.raster.stack), ignore.case = TRUE)[!(grep("Temp",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 GDD_sq = getValues(pred.raster.stack[[grep("GDD",names(pred.raster.stack))[(grep("GDD",names(pred.raster.stack)) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Precip_sq = getValues(pred.raster.stack[[grep("Precip",names(pred.raster.stack), ignore.case = TRUE)[(grep("Precip",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Water_Bal_sq = getValues(pred.raster.stack[[grep("Water",names(pred.raster.stack), ignore.case = TRUE)[(grep("Water",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]),
                 Min_temp_sq = getValues(pred.raster.stack[[grep("Temp",names(pred.raster.stack), ignore.case = TRUE)[(grep("Temp",names(pred.raster.stack), ignore.case = TRUE) %in% grep("_sq",names(pred.raster.stack)))]]]))
    # Creating tmp raster
    gl.pred.tmp <-
      raster(matrix(predict(test.mods$Models[['glm.gaussian']],
                            df.dat.tmp, type = 'response'),
                    nrow = pred.raster.stack@nrows,
                    ncol = pred.raster.stack@ncols,
                    byrow = TRUE),
             crs = crs(pred.raster.stack))
    extent(gl.pred.tmp) <- extent(pred.raster.stack)
    # Saving tmp raster
    writeRaster(gl.pred.tmp,
                paste0(path.sdms.write.tmp,
                       '/',
                       df$species_name[k],
                       '_tmp_sdm_tile_',
                       index.list,
                       ".tif"),
                overwrite = TRUE)
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



split.raster.function <-
  function(raster.stack, esh.raster) {
    # Getting number of cells in the raster
    n.cells <- raster.stack@nrows * raster.stack@ncols
    n.rasters <- round(n.cells/3000000,digits=0)
    
    if(n.rasters %in% 0) {
      n.rasters <- 1
    }
    
    # Getting number of rasters in the raster stack
    div <- seq_len(n.rasters) 
    div <- div[n.rasters %% div == 0]
    # Making sure it isn't a prime number
    while(length(div) <= 2 & n.rasters >= 5) {
      n.rasters <- n.rasters + 1
      div <- seq_len(n.rasters)
      div <- div[n.rasters %% div == 0]
    }
    # And getting number of sides to split into
    if(length(div) %% 2 == 0 & n.rasters != 3) {
      sides.x = div[length(div)/2]
      sides.y = div[length(div)/2 + 1]
    }
    if(length(div) %% 2 != 0 & n.rasters != 3) {
      sides.x = div[ceiling(length(div)/2)]
      sides.y = div[ceiling(length(div)/2)]
    }
    if(n.rasters == 3) {
      if(raster.stack@nrows > raster.stack@ncols) {
        sides.x = 3
        sides.y = 1
      } else{
        sides.y = 3
        sides.x = 1
      }
      
    }
    
    # Splitting climate rasters
    # This gets extents
    # split.current.extents = splitRaster(raster.stack$GDD_2010, nx = sides.x, ny = sides.y, buffer = c(1,1))
    split.current.extents = splitRaster(esh.raster, nx = sides.x, ny = sides.y, buffer = c(1,1))
    
    # Which values have NAs
    max.values <- c()
    split.current.extents.tmp <- list()
    counter = 1
    for(i in 1:length(split.current.extents)) {
      # Getting max value of raster
      tmp.max.values <- maxValue(split.current.extents[[i]])
      # Adding to list if true
      if(!is.na(tmp.max.values)) {
        split.current.extents.tmp[[counter]] <- split.current.extents[[i]]
        counter = counter + 1
      }
    }
    
    # Updating name of object
    split.current.extents <- split.current.extents.tmp
    # 
    # # Now repeating to limit to rasters within appropriate ecoregion
    # split.current.extents.tmp <- list()
    # counter = 1
    # for(i in 1:length(split.current.extents)) {
    #   tmp.ecoregion = crop(ecoregions.map,split.current.extents[[i]])
    #   if(sum(unique(getValues(tmp.ecoregion)) %in% mod.list$Ecoregions$ecoregion[mod.list$Ecoregions$ecoregion != 0]) >= 1) {
    #     split.current.extents.tmp[[counter]] = split.current.extents[[i]]
    #     counter = counter + 1
    #   }
    # }
    # 
    # # Updating name of object
    # split.current.extents <- split.current.extents.tmp
    # 
    # # Now repeating to limit to rasters within appropriate biome
    # split.current.extents.tmp <- list()
    # counter = 1
    # for(i in 1:length(split.current.extents)) {
    #   tmp.ecoregion = crop(glob.cover,split.current.extents[[i]])
    #   if(sum(unique(getValues(tmp.ecoregion)) %in% mod.list$GlobBiomes$glob_cov[mod.list$GlobBiomes$glob_cov != 0]) >= 1) {
    #     split.current.extents.tmp[[counter]] = split.current.extents[[i]]
    #     counter = counter + 1
    #   }
    # }
    # 
    # split.current.extents <- split.current.extents.tmp
    rm(split.current.extents.tmp)
    # 
    removeTmpFiles()
    gc()
    
    return(split.current.extents)
  }



chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 


