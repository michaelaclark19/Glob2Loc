##### General
This GitHub folder contains the scripts used in the paper titled: "Cumulative Impact of Multiple Anthropogenic Stressors to Terrestrial Vertebrate Biodiversity".
Each folder contains scripts for different purposes e.g. land cover forecasting, estimating agricultural intensity, calculating biodiversity impacts, constructing distribution models, etc.
The names of scripts in each folder start with a number that indicate the order in which the scripts should be run, starting with the script named "1.0...". The "0.0..." scripts contain functions that are called in other scripts and do not need to be run. 


##### Running Order of Different Parts of the Analyses
This list below indicates the order in which the different types of scripts (e.g. land cover forecasting, estimating agricultural intensity, etc) should be run to complete the analysis.
1. Climate Change Data
2. Creating GDD Rasters
3. Species Richness
4. Population Density Estimates
5. Urbanisation Projections
6. Land Cover Modelling
7. Agricultural Intensity
8. SDM Scripts
9. Calculating Biodiversity Scripts
10. Analysis Scripts
11. Figure Scripts
12. Table Scripts



##### High-level description of the scripts in each folder.

###
1. Climate Change Data: This contains scripts that manage raw climate change data (from CMIP6), and reproject/rescale/etc it to the projection used in this analysis. The output from these scripts are global raster maps of current and projected future climate data that will be used throughout this analysis.
"0.0...": Moves climate maps from one directory to another. Does not need to be run if the climate raster maps are already in the correct folder.
"1.0...": Calculates precipitation in the warmest 3 months of the year
"1.1...": Calculates mean temperature across the year
"2.0...": Calculates growing degree days and variance in precipitation (Across months)
"3.0...": Reprojects rasters
"4.0...": Reprojects the minimum temperature rasters
"5.0...": Calculates water balance by cell
"6.0...": Calculates annual precipitation by cell

###
2. Creating GDD Rasters: This creates GDD rasters - this is an old script and does not need to be run.

###
3. Species Richness: This contains scripts that calculates species richness using area of habitat maps (amphibians, birds, mammals) and IUCN range maps (reptiles). The output from this script is a series of maps that contain current estimated species richness for amphibians, birds, mammals, and reptiles.
"1.0...": There are four versions of this script. Each version is used to estimate species richness for a different taxa (amphibians, birds, mammals, reptiles).

###
4. Population Density: This gets coefficients used to estimate population density of terrestrial vertebrates. It is based off Santini et al (2017). https://onlinelibrary.wiley.com/doi/pdfdirect/10.1111/geb.12758. The output from these scripts are correlation coefficients that are used to estimate population density, measured in the number of mature breeding individuals per sq km, for each species.
"1.0...": There are two versions of this script. One version is for amphibians, birds, and mammals. The other is for reptiles. Both scripts recreate the method from Santini et al (2017), but using the CMIP6 climate data for internal consistency in this analysis.

###
5. Urbanisation Projections: This projects urban land cover from 2010 to 2050 using data from Seto et al (2019). The output from these scripts are projected urban land cover from 2010 to 2050 at 5-year intervals.
"0.0...": This contains the functions used to project urban land cover to 2050, and then runs these functions. 

###
6. Land Cover Modelling: This contains scripts that allocate projected changes in agricultural land use from national levels (i.e. non-spatial) to subnational (i.e. spatial). E.g., where in Brazil an additional 100 sq km of land would be allocated.
"0.0...": This contains the functions used to allocate land use.
"1.0...": There are two versions. These two scripts estimate the correlation between historic land cover change and variables that are known to be associated with land cover change (e.g. proximity to markets, location in national parks and other protected areas, density of agricultural production, etc). These two scripts can be run in any order.
"2.0...": This script runs the projections of agricultural land use for each country. It creates rasters of projected cropland and pastureland use from 2010 to 2050 at 5-year intervals for each country.
"3.0...": This script stitches the country maps back into a global map.
Mitigation Scenario Folder: This contains scripts used to run the mitigation scenario. The structure of the scripts is as indicated above, with the exception that there is no "1.0..." script.

###
7. Agricultural Intensity: This contains a series of scripts that are used to project cropland intensity at 5-year intervals from 2010 to 2050. It does so by first creating a multinomial classification model that estimates cropland intensity, and then applying this model to the cropland projections created in the 'Land Cover Modelling' folder.
"0.0...": Contains functions used to weight different regression models, and then plot/analyse the accuracy of these regression models.
"1.0...": This develops and then saves the multinomial classification model that is used to project pastureland intensity. It tests two model forms - a multinomial classification model and a random forest classification model - and ultimately uses the multinomial classification model due to higher overall accuracy. 
"2.0...": This projects cropland intensity using the outputs from the "Land Cover Modelling" agricultural land use projections.
"3.0...": This classifies species into different rarity traits. This is then used to understand how species respond to different types of human modified land cover.
Mitigation Scenario Folder: Contains the script used to project cropland intensity for the mitigation scenario.

###
8. SDM Scripts: This contains a series of scripts that are used to project each species' potential habitat under different scenarios of climate change. The outputs are, for each species, an estimate of each species' climate-suitable and climate-accessible habitat 5-year intervals from 2010 to 2050.
"0.0...": This contains the functions used to create the distribution models.
"3.0...": This creates the project habitat maps for each species. This script also removes species from the analysis if the species are mostly found on water, have limited climate variability in their current habitat range, and for several other reasons. There are two versions of this script - one works for reptiles, and the other works for birds, mammals, and amphibians.
"4.0...": This tests the accuracy of the projected distribution maps. Accuracy is measured as the similarity between the species' modelled habitat in 2010 and the area of habitat map (birds, mammals, amphibians) and IUCN range map (reptiles).
"5.0...": This gets climate suitable habitat maps for each species.
"5.1...": This removes water from each species' climate suitable habitat map. This does not need to be run, as this is also done in the "6.0..." script in this folder, and also in the calculating biodiversity scripts. 
"6.0...": This estimates each species' climate accessible habitat maps.
Mitigation Scenario Folder: Contains scripts as above, but projects climate suitable and climate accessible habitat areas under the SSP1 climate change scenario.

###
9. Calculating Biodiversity Scripts: This calculates estimated remaining habitat area, population abundance, and habitat connectivity for each species. Each of these outcomes are measured for different combinations of stressors. For each species, outputs are saved as a .csv file, with option to save spatial estimates as raster maps (.tif files).
"0.0...": Contains the functions used to calculated biodiversity outcomes.
"1.0...": Estimates the biodiversity outcomes for each species. There are a series of files with this name - one file per 5 year interval from 2020 to 2050, and then two additional files to estimate biodiversity outcomes when using a species' climate suitable habitat range in 2020 and 2050. These scripts calculate outcomes for each species' total habitat area, and for each habitat patch in each species' remaining habitat area. They do not calculate outcomes for populations. 
"2.0...": This estimates habitat connectivity for each species - i.e. uses the outputs from the "1.0..." scripts to identify populations for each species.
Mitigation Scenario Folder: As above, but scripts that run on the mitigation scenario.

###
10. Analysis Scripts. Contains scripts that compile and aggregate results from the "Calculating Biodiversity Scripts" folder. 
"0.0...": Compiles estimate species-level biodiversity outcomes. 
"1.0...": Compiles patch-level estimates for each species into a format that takes less memory space. There are multiple versions of this script depending on what type of results are being compiled.
"2.0...": Compiles results for all species into a single big data frame. There are multiple versions of this script depending on what type of results are being compiled.
Country Analyses Folder: Contains a series of scripts that takes species-level spatial results and converts them into country-level raster maps indicating change in biodiversity outcomes for all species in each raster cell.
Species By Geographic Region Folder: Calculates non-spatial biodiversity results for species found in different biodiversity hotspots or ecoregions. There are scripts in this folder that perform this calculation when defining whether species are in an ecoregion (or hotspot) based on their area of habitat map or their modelled climate accessible habitat map.
Mitigation Scenario Folder: Contains scripts that compile species-level results and then aggregate results across all species. 

###
11. Figure Scripts: Contains the scripts used to make the figures in the manuscript. All scripts in this folder are named so that they refer back to the figure number in the manuscript text.

###
12. Table Scripts: Contains the scripts used to make the tables in the manuscript. All scripts in this folder are named so that they refer back to the table number in the manuscript text.



