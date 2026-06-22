# 0.0
# Functions for forecasting crop and pasture models

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


# 0.0
# Functions for forecasting crop and pasture models

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
