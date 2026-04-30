#####
# Trial code
#####

###
# Libraries
library(raster)

###
# Creating a temporary matrix
# tmp.matrix <- 
#   matrix(c(rep(1, 30),
#          rep(0,20),
#          rep(2,20),
#          rep(0,10),
#          rep(3,20),
#          rep(0,10)),
#              nrow = 10,
#              ncol = 11)
# 
# tmp.matrix <- 
#   matrix(rep(c(rep(1, 30),
#            rep(0,20),
#            rep(2,20),
#            rep(0,10),
#            rep(3,20),
#            rep(0,10)),
#            1600),
#          nrow = 400,
#          ncol = 440)
# 
# vector.matrix <-
#   matrix(1:110,
#          byrow = TRUE,
#          nrow = 10,
#          ncol = 11)
# 
# tmp.raster <-
#   raster(tmp.matrix)

#####
# Creating function
# Goal of this function is to look at which patches are 
# Within a certain buffer of each other
# Where a buffer is a rectangle around the target cell
# E.g. a buffer of 2 creates a rectangle with that has
# Dimensions that are 5 x 5, or 25 total cells
#####

connectivity.function <-
  function(input.raster,
           buffer.cells) {
    
    # ###
    # # Converting vector into a raster
    # input.raster <-
    #   raster(matrix(data.vector,
    #                 nrow = template.raster@nrows,
    #                 ncol = template.raster@ncols,
    #                 byrow = TRUE),
    #          crs = crs(template.raster))
    # extent(input.raster) <- extent(template.raster)
    
    n.row = input.raster@nrows
    n.col = input.raster@ncols
    
    ### 
    # Converting into a vector
    tmp.vector <-
      getValues(input.raster)
    
    #####
    # Add exception here so that script only loops over cells that 
    # Don't correspond with NA or 0
    #####
    
    ###
    # Getting location of cells that contain values > 0
    tmp.vector2 <-
      which(tmp.vector >= 1 &
                   !is.na(tmp.vector))
    
    ###
    # Creating empty df
    # This will be appended to later
    stacked.df <-
      c()
    
    ###
    # Looping through cells in the raster
    for(i in 1:length(tmp.vector2)) {
      ###
      # Getting column number of the point
      rownum = 
        ceiling(tmp.vector2[i] / n.col)
      
      ###
      # Getting row number of the point
      colnum = 
        tmp.vector2[i] -
        (rownum - 1) * n.col
      
      # Creating vector that contains colnums of the buffer
      cols.buffer <-
        (colnum - buffer.cells) : (colnum + buffer.cells)
      
      cols.buffer <-
        cols.buffer[cols.buffer >= 1 &
                      cols.buffer <= n.col]
      
      # Creating vector that contains rownums of the buffuer
      rows.buffer <-
        (rownum - buffer.cells) : (rownum + buffer.cells)
      
      rows.buffer <-
        rows.buffer[rows.buffer >= 1 &
                      rows.buffer <= n.row]
      
      ### 
      # Getting value of corners of buffer cells
      
      # Upper left corner
      up.left <-
        n.col *
        (rows.buffer[1] - 1) +
        cols.buffer[1]
      
      # Lower left corner
      up.right <-
        n.col *
        (rows.buffer[1] - 1) +
        cols.buffer[length(cols.buffer)]
      
      # Upper right corner
      low.left <-
        n.col *
        (rows.buffer[length(rows.buffer)] - 1) +
        cols.buffer[1]
      
      # Lower right corner
      low.right <-
        n.col *
        (rows.buffer[length(rows.buffer)] - 1) +
        cols.buffer[length(cols.buffer)]
      
      ###
      # Creating vector that contains value of cell indices in the buffer
      
      # First getting values of the top of the box
      left.vals <-
        seq(up.left,
            low.left,
            by = n.col)
      
      # Next getting values of the bottom of the box
      right.vals <-
        seq(up.right,
            low.right,
            by = n.col)
      
      # Saving empty vector
      index.values = c()
      
      # And appending these two lists together
      for(j in 1:length(left.vals)) {
        index.values <-
          c(index.values,
            left.vals[j]:right.vals[j])
      }
      
      
      ###
      # Getting values in the matrix in these locations
      tmp.values <-
        unique(tmp.vector[index.values])
      
      
      ###
      # Creating data frame
      stacked.df <-
        rbind(stacked.df,
              data.frame(Patch_1 = rep(tmp.vector[tmp.vector2[i]],
                                       length(tmp.values)),
                         Patch_2 = tmp.values))
    }
    
    stacked.df <-
      unique(stacked.df)
    stacked.df <-
      stacked.df[stacked.df$Patch_2 != 0,]
    stacked.df <-
      stacked.df[stacked.df$Patch_1 < stacked.df$Patch_2,]
    
    return(stacked.df)
}

#####
# Second function
# This moves a matrix rather than moving a box around a vector
#####


patch.connectivity.function <-
  function(input.raster,
           buffer.cells) {
    
    # ###
    # # Converting vector into a raster
    # input.raster <-
    #   raster(matrix(data.vector,
    #                 nrow = template.raster@nrows,
    #                 ncol = template.raster@ncols,
    #                 byrow = TRUE),
    #          crs = crs(template.raster))
    # extent(input.raster) <- extent(template.raster)
    
    
    n.row = input.raster@nrows
    n.col = input.raster@ncols
    
    ### 
    # Converting into a vector
    tmp.matrix <-
      as.matrix(input.raster)
    
    ###
    # Converting 0s to NAs
    tmp.matrix[tmp.matrix == 0] <-
      NA
    
    #####
    # Add exception here so that script only loops over cells that 
    # Don't correspond with NA or 0
    #####
    
    connect.df <-
      c()
    
    ###
    # Looping to move second matrix in x direction
    for(i in 0:buffer.cells) {
      
      # Moving the second matrix
      # First moving it left
      # Which means the columns on the left get cut off
      tmp.matrix.left <-
        tmp.matrix[, (i + 1) : n.col]
      
      # And repositioning the first matrix
      # Need to cut off the cells on the right
      tmp.matrix1 <-
        tmp.matrix[, 1 : (n.col - i)]
      
      # And getting list of cells where these values do not match
      connect.df <-
        rbind(connect.df,
        unique(data.frame(Patch_1 = tmp.matrix1[tmp.matrix1 != tmp.matrix.left],
                          Patch_2 = tmp.matrix.left[tmp.matrix.left != tmp.matrix1])))
      
      
      
      ###
      # Looping through to move matrix y direction
      for(j in 0:buffer.cells) {
        ###
        # Only need to look at overlap if i or j > 1
        # If i and j are both < 1
        # This says to look in cells immediately adjacent to target cell
        # Which is already looked at by frag stats
        
        if(i > 1 | j > 1) {
          ###
          # Moving second matrix up
          tmp.matrix.left.up <-
            tmp.matrix.left[(j + 1) : n.row,]
          # Repositioning first matrix
          tmp.matrix2 <-
            tmp.matrix1[1 : (n.row - j),]
          
          ###
          # Moving second matrix down
          tmp.matrix.left.down <-
            tmp.matrix.left[1 : (n.row - j), ]
          # Repositioning second matrix
          tmp.matrix3 <-
            tmp.matrix1[(j + 1) : n.row,]
          
          ###
          # Getting data frame of values where these don't match
          connect.df <-
            rbind(connect.df,
                  unique(data.frame(Patch_1 = tmp.matrix2[tmp.matrix2 != tmp.matrix.left.up],
                                    Patch_2 = tmp.matrix.left.up[tmp.matrix.left.up != tmp.matrix2])))
          
          connect.df <-
            rbind(connect.df,
                  unique(data.frame(Patch_1 = tmp.matrix3[tmp.matrix3 != tmp.matrix.left.down],
                                    Patch_2 = tmp.matrix.left.down[tmp.matrix.left.down != tmp.matrix3])))
        }
      }

      # Now repeating if the second matrix moves to the right
      # Which means the cells on the end get cut off
      tmp.matrix.right <-
        tmp.matrix[, 1 : (n.col - i)]
      
      # And repositioning the first matrix
      tmp.matrix1 <-
        tmp.matrix[, (i+1) : n.col]
      
      # And getting connectivity between these
      connect.df <-
        rbind(connect.df,
              unique(data.frame(Patch_1 = tmp.matrix1[tmp.matrix1 != tmp.matrix.right],
                                Patch_2 = tmp.matrix.right[tmp.matrix.right != tmp.matrix1])))
      
      ###
      # Looping to move second matrix in y direction
      for(j in 0:buffer.cells) {
        
        ###
        # Only need to look at overlap if i or j > 1
        # If i and j are both < 1
        # This says to look in cells immediately adjacent to target cell
        # Which is already looked at by frag stats
        if(i > 1 |
           j > 1) {
          ###
          # Moving second matrix up
          tmp.matrix.right.up <-
            tmp.matrix.right[(j + 1) : n.row,]
          # Repositioning first matrix
          tmp.matrix2 <-
            tmp.matrix1[1 : (n.row - j),]
          
          ###
          # Moving second matrix down
          tmp.matrix.right.down <-
            tmp.matrix.right[1 : (n.row - j), ]
          # Repositioning second matrix
          tmp.matrix3 <-
            tmp.matrix1[(j + 1) : n.row,]
          
          ###
          # Getting data frame of values where these don't match
          connect.df <-
            rbind(connect.df,
                  unique(data.frame(Patch_1 = tmp.matrix2[tmp.matrix2 != tmp.matrix.right.up],
                                    Patch_2 = tmp.matrix.right.up[tmp.matrix.right.up != tmp.matrix2])))
          
          connect.df <-
            rbind(connect.df,
                  unique(data.frame(Patch_1 = tmp.matrix3[tmp.matrix3 != tmp.matrix.right.down],
                                    Patch_2 = tmp.matrix.right.down[tmp.matrix.right.down != tmp.matrix3])))
        }
      }
    }
        
  # Only need unique combinations
  connect.df <-
    unique(connect.df)
  # And don't need transitivie combinations
  # E.g. patch 2 connecting to patch 3 is the same as patch 3 connecting to patch 2
  connect.df <-
    connect.df[connect.df$Patch_1 < connect.df$Patch_2,]
  
  # And removing NAs
  connect.df <-
    connect.df[complete.cases(connect.df),]

    
  return(connect.df)
}

###
# Note that the output of this does not identify patches connected to themselves
# Eg will not return a row where patch_1 == 1 and patch_2 == 1

# t1 = Sys.time()
# trial.matrices <-
#   patch.connectivity.function(input.raster = tmp.raster,
#                               buffer.cells = 3)
# Sys.time() - t1


#####
# Now create function that gets
# (a) unique combinations from the patches, and 
# (b) updates values in the matrix so that they are considered the same patch
#####
# dat <-
#   trial.matrices


###
# Updating patch numbers
# The output of this function is a raster that has updated patch numbers
# Which will then be overlaid with a raster that contains
# Habitat area by cell to get total size of each patch
update.patch.numbers.function <-
  function(dat,
           input.raster) {
    
    # Creating new column
    dat$New_Patch <-
      dat$Patch_2
    # Getting list of unique patches
    unique.patches <-
      sort(unique(dat$Patch_2))
    
    # And second list of unique pathces
    unique.patches.tmp <-
      sort(unique(dat$Patch_1))
    
    ###
    # Need to add a loop to allow for connectivity between intermediate patches
    for(i in length(unique.patches.tmp):1) {
      tmp.vector <- sort(dat$Patch_2[dat$Patch_1 == unique.patches.tmp[i]])
      
      tmp.vector2 <- sort(dat$Patch_1[dat$Patch_2 %in% tmp.vector])
      tmp.vector2 <- sort(c(tmp.vector,tmp.vector2))
      while(length(tmp.vector) != length(tmp.vector2)) {
        tmp.vector <- tmp.vector2
        tmp.vector2 <-
          sort(unique(c(tmp.vector2, dat$Patch_2[dat$Patch_1 %in% tmp.vector2])))
        tmp.vector2 <-
          sort(unique(c(tmp.vector2, dat$Patch_1[dat$Patch_2 %in% tmp.vector2])))
      }
      
      dat$New_Patch[dat$Patch_1 %in% tmp.vector] <- i
      dat$New_Patch[dat$Patch_2 %in% tmp.vector] <- i
    }
      
    
    # ###
    # # First condensing into single patches
    # for(i in length(unique.patches) : 1) {
    #   # Updating values
    #   dat$New_Patch[dat$New_Patch == unique.patches[i]] <-
    #     min(dat$Patch_1[dat$Patch_2 == unique.patches[i]])
    # }
    
    ### 
    # And now updating values in rasters
    # Doing this by converting to raster to a vector
    # Updating values
    # And then reconverting back into a raster
    # This is much faster than updating values in a raster format
    tmp.vector <-
      getValues(input.raster)
    
    # Getting unique values of the new patches
    unique.patches <-
      unique(dat$New_Patch)
    
    # Updating values in the vector
    for(i in 1:length(unique.patches)) {
      tmp.vector[tmp.vector %in% unique(dat$Patch_2[dat$New_Patch == unique.patches[i]])] <-
        unique.patches[i]
      tmp.vector[tmp.vector %in% unique(dat$Patch_1[dat$New_Patch == unique.patches[i]])] <-
        unique.patches[i]
    }
    
    # Reconverting back into a raster
    remade.raster <-
      raster(matrix(tmp.vector,
                    nrow = input.raster@nrows,
                    ncol = input.raster@ncols,
                    byrow = TRUE),
             crs = crs(input.raster))
    # Setting resolution
    extent(remade.raster) <-
      extent(input.raster)
    
    # And returning the raster
    return(remade.raster)
  }


#####
# And now creating a function that gets sizes of patches
# That are larger than the size reestimated to contain a minimum viable population
#####

mvp.patch.function <-
  function(patch.raster,
           loss.raster,
           mvp_size) {
    
    ###
    # Converting to vectors
    # Because R is really slow with rasters
    
    # Vector with patch numbers
    patch.vector <-
      getValues(patch.raster)
    
    # Vector that contains proportion of each cell that is habitable for the species
    loss.vector <- 
      getValues(loss.raster)
    
    # Converting to a matrix
    df <-
      matrix(c(patch.vector,
               loss.vector,
               rep(1, length(patch.vector))),
             nrow = length(patch.vector),
             ncol = 3,
             byrow = FALSE)
    
    # Removing NAs
    # Because this freaks out the function
    df <-
      df[!is.na(df[,1]),]
    
    # And now taking row sums
    sum.df <-
      as.data.frame(rowsum(df,
            group = df[,1],
            na.rm = TRUE))
    
    # Getting patch number
    sum.df[,1] <-
      sum.df[,1] /
      sum.df[,3]
    
    # Limiting to patches that are larger than the esimated minimum viable population
    sum.df <-
      sum.df[sum.df[,2] >= mvp_size,]
    
    # Removing patches that aren't adequately large
    patches.keep <-
      as.vector(unique(sum.df[,1]))
    
    patch.vector[!(patch.vector %in% unique(sum.df[,1])) &
                   !is.na(patch.vector)] <-
      NA
    
    # Recreating the patch raster
    patch.raster.new <-
      raster(matrix(patch.vector,
                    nrow = patch.raster@nrows,
                    ncol = patch.raster@ncols,
                    byrow = TRUE),
             crs = crs(patch.raster))
    
    # Updating extent
    extent(patch.raster.new) <-
      extent(patch.raster)
    
    # And returning the raster
    return(patch.raster.new)
    
    ###
    # And checking to make sure this works
    # Do not need to run this chunk of code
    # This takes a while on larger rasters
    # patch.raster.new[patch.raster.new == 0] <- NA
    # plot(patch.raster.new)

  }


###
# And testing the functions
# trial.matrices <-
#   patch.connectivity.function(input.raster = tmp.raster,
#                               buffer.cells = 2)
# 
# updated.patch.raster <-
#   update.patch.numbers.function(dat = trial.matrices,
#                                 input.raster = tmp.raster)

