library(plyr)
library(dplyr)
library(stringr)

# Setting working directory
setwd('/Users/macuser/Desktop')

# List of raw data files
file.list <- 
  list.files(path = '/Users/macuser/Downloads/',
             pattern = 'CMIP6',
             full.names = TRUE) %>%
  unique(.)

# List of years and ssp scenarios
years <- str_extract(file.list,'[0-9]{4,4}-[0-9]{4,4}') %>% unique(.) %>% sort(.)
ssps <- str_extract(file.list,'Historic|SSP.{5,5}') %>% unique(.) %>% sort(.) %>% c('Historic',.)
indicators <- c('Total_precip','Mean_temp','Minimum_temp')

# Main directory
dir.create(paste0(getwd(),'/CMIP6_Climate_Data'))

# Checking that all directories have 12 files
directory.append <- c()

# Appending other directories
for(s in ssps) {
  # Creating ssp directory
  dir.create(paste0(getwd(),'/CMIP6_Climate_Data/',s))
  for(y in years) {
    # Creating year directory
    dir.create(paste0(getwd(),'/CMIP6_Climate_Data/',s,'/',y))
    # Creating indicator directory
    for(i in indicators) {
      dir.create(paste0(getwd(),'/CMIP6_Climate_Data/',s,'/',y,'/',i))
      
      # List of files to move
      files.move <- 
        file.list %>%
        .[grepl(y,.,ignore.case = TRUE)] %>%
        .[grepl(gsub('_',' ',i),.,ignore.case = TRUE)]
      
      if(s != 'Historic') {
        files.move <- 
          files.move %>%
          .[grepl(s,.,ignore.case = TRUE)]
      }
      
      # Directory to move files to...
      files.move.to <-
        gsub('/Users/macuser/Downloads/',
             paste0(getwd(),'/CMIP6_Climate_Data/',s,'/',y,'/',i),
             files.move)
      
      # And moving files
      file.copy(from = files.move,
                to = files.move.to)
      
      # Checking number of files - should be twelve here
      # If not twelve files, then saving the directory
      if(length(list.files(paste0(getwd(),'/CMIP6_Climate_Data/',s,'/',y,'/',i))) != 12) {
        directory.append <- c(directory.append,paste0(getwd(),'/CMIP6_Climate_Data/',s,'/',y,'/',i))
      } # Else do nothing
      
    } # End loop through indicators
  } # End loop through years
} # End loop through SSPs

# Checking directories that do not have 12 files
directory.append

# Filtering out the ones that should not have twelve files
directory.append <-
  directory.append %>%
  .[!(grepl('2[0-9]{3,3}-',.) & grepl('Historic',.))] %>%
  .[!(grepl('1[0-9]{3,3}-',.) & grepl('SSP',.))]

lapply(directory.append,function(i) {list.files(i) %>% sort()})

# END
