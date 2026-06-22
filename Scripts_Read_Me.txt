##### General
This GitHub folder contains the scripts used in the paper titled: "Cumulative Impact of Multiple Anthropogenic Stressors to Terrestrial Vertebrate Biodiversity".
Each folder contains scripts for different purposes e.g. land cover forecasting, estimating agricultural intensity, calculating biodiversity impacts, constructing distribution models, etc.

##### Contents of this file
# Running Order of Different Parts of the Analyses
# High-level description of the scripts in each folder
# Description of the Demo Dataset.
# Description of the folders contained I the Demo Dataset
# Description of the System Used to Run the Analyses

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
"1.1...": This is an old script - does not need to be run. Calculates mean temperature across the year
"2.0...": Calculates growing degree days and variance in precipitation (Across months). The mean annual temperature calculations have been incorporated into this script.
"3.0...": Reprojects rasters. This version of the script uses the projectRaster function from rgdal and gdalUtils packages. I suggest updating these to the terra package equivalent if rerunning the scripts. This script takes a long time to run. 
"4.0...": Reprojects the minimum temperature rasters. This version of the script uses the projectRaster function from rgdal and gdalUtils packages. I suggest updating these to the terra package equivalent if rerunning the scripts. This script takes a long time to run.
"5.0...": Calculates water balance by cell
"6.0...": Calculates annual precipitation by cell

###
2. Creating GDD Rasters: This creates GDD rasters - this is an old script and does not need to be run.

###
3. Species Richness: This contains scripts that calculates species richness using area of habitat maps (amphibians, birds, mammals) and IUCN range maps (reptiles). The output from this script is a series of maps that contain current estimated species richness for amphibians, birds, mammals, and reptiles.
"1.0...": There are four versions of this script. Each version is used to estimate species richness for a different taxa (amphibians, birds, mammals, reptiles). This script takes a long time to run.

###
4. Population Density: This gets coefficients used to estimate population density of terrestrial vertebrates. It is based off Santini et al (2017). https://onlinelibrary.wiley.com/doi/pdfdirect/10.1111/geb.12758. The output from these scripts are correlation coefficients that are used to estimate population density, measured in the number of mature breeding individuals per sq km, for each species.
"1.0...": There are two versions of this script. One version is for amphibians, birds, and mammals. The other is for reptiles. Both scripts recreate the method from Santini et al (2017), but using the CMIP6 climate data for internal consistency in this analysis.

###
5. Urbanisation Projections: This projects urban land cover from 2010 to 2050 using data from Seto et al (2019). The output from these scripts are projected urban land cover from 2010 to 2050 at 5-year intervals.
"0.0...": This contains the functions used to project urban land cover to 2050, and then runs these functions. This scripts take a while to run. It also uses functions from rgdal and gdaUtils. If you plan to run this, I suggest updating the scripts to run on the terra:: versions of the same functions.

###
6. Land Cover Modelling: This contains scripts that allocate projected changes in agricultural land use from national levels (i.e. non-spatial) to subnational (i.e. spatial). E.g., where in Brazil an additional 100 sq km of land would be allocated.
"0.0...": These files contain the functions that are called in the other land cover modelling scripts.
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
Note that several of these scripts require rgdal and gdalUtils, which may not be available on current versions of R (>V 4.2). If you have a newer version of R, I strongly suggest updating the rgdal and gdalUtils functions with the terra:: equivalents. I have provided outputs of the "5.0" and "6.0" scripts for the 100 randomly selected species that are included in the demo dataset.
These scripts (particularly the 3.0, 5.0, and 6.0 scripts) take a very long time to run.
"0.0...": This contains the functions used to create the distribution models.
"1.0...": This script randomly selects 100 species (25 from each taxa) for the demo dataset.
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
"2.0...": This estimates habitat connectivity for each species - i.e. uses the outputs from the "1.0..." scripts to identify populations for each species. There are two versions of these scripts - each version works on one of the two different SDM thresholds used to identify a species' habitat range.
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

##### Description of the Demo Dataset.
The demo dataset provided as a supplementary material runs the full analysis on a small subset of species (40) that exist in one geographic region (Meso America). It contains all of the data inputs needed to run all scripts used in the analysis. The R scripts used for the full analysis are available on the project GitHub (https://github.com/michaelaclark19/Glob2Loc/). 
The demo dataset contains a small number of species (99) that exist in one geographic region (Meso America) due to the storage space and RAM required to run the full analysis. 
The demo dataset is designed to provide an idea of how all scripts and parts of the analysis work. It was tested with R Version 4.2.2 and Python on a MacBook Pro with 32GB of RAM.

##### Description of the folders contained in the Demo Dataset

Ag Intensity Outputs: Contains four subfolders, that collectively contain (a) R objects (.Rda files) that contain the saved classification models used to predict cropland intensity (Classification Model Outputs), (b) cropland intensity projections in different land cover modelling scenarios (Cropland Intensity Forecasts), (c) summary plots of the multinomial classification models used to predict agricultural intensity (Classification_Model_Outputs), and (d) correlation coefficients describing how different types of species respond to different types of land use (Response to Habitat Intensification). This folder and its subfolders will be created when running the 'Agricultural Intensity' scripts.

Analyses: This folder will be created by running the analysis scripts, and will contain outputs created during these scripts.

CMIP6_Climate_Data: Contains raw and manipulated (e.g. reprojected, summarised, etc) from the CMIP6 model ensemble. This contains three subfolders, with each subfolder containing data for a different climate projection (Historic; SSP1-2.6; SSP2-4.5). Each of these folders in turn contain subfolders that either (a) contain data for different periods of time provided by CMIP6, or (b) global rasters that contain interpolated values between the time periods CMIP6 provides data for (these are in the folder named Interpolated_Rasters).

COO_Data: This contains a series of CSV files that indicate the country(ies) a single species is estimated to be found in, based on each species' area of habitat maps (for amphibians, birds, and mammals) or IUCN range maps (for reptiles).

Ecoregions_Feb2023: This contains a series of raster maps that indicate the spatial boundaries of different ecoregions, realms, biomes, the intersection of realms and biomes, and biodiversity hotspots. It also contains a csv file indicating the boundaries of each realm by biome combination. These boundaries are indicated by the min and max row and column indices each realm/biome combination is found in the realm/biome combination raster.

ESA LandCov Maps: Contains historic raster maps derived from ESA GlobMod. These maps have been reprojected to the CRS and resolution used in this analysis, and then updated so that the sum of all land covers do not exceed 1.

ESH_RCPs: This folder will be created when running the 'SDM Scripts', and contains outputs that are created from running the SDM Scripts. There are five subfolders in this folder. Three subfolders correspond with the different CMIP6 scenarios (Historic; SSP1-2.6; and SSP2-4.5). One of the other folders contain estimated accuracy as measured by AUC (in the folder named 'Mod_Accuracy'. The other folder contains overall estimated accuracy of the weighted SDM maps (folder named 'Weighted_Threshold_Accuracy'). Note that the 'Weighted_Threshold_Accuracy' folder only contains outputs for species which pass quality checks when constructing the SDMs. This folder only contains data for ~100 species due to file size limitations. 

ESH_Tifs_12Oct: This contains area of habitat (amphibians, birds, mammals) or IUCN range maps (reptiles) for ~100 species due to file size limitations. If you would like to rerun the analysis using different base habitat range maps, replace the files in these folders with the habitat range maps you would like to use. There are four subfolders in this folder, with one subfolder for each taxa.

Files_For_Fig2 (and 3-5): These folders contain the data inputs used to make Figs 2-5. The files to make Figs 1 are saved in Outputs/Aggregated_CSV_Outputs

Forecasting_Data: This contains a csv for the region 'MesoAmerica'. This csv file contains the raw data that (a) is used to create the agricultural land use regression models; and (b) that is then used to project agricultural land cover change to 2050. This data is limited to one world region due to file size limitations.

Global Mollweide Maps: This folder contains a series of raster files that are used in different parts of the analysis. File names correspond with data in each raster file. The subfolder (CurrentUrbanGlobMod) contains a geo package .gpkg file of current urban land cover.

Land Forecast Outputs: This folder is created by running the 'Land Cover Modelling' scripts. The structure of this folder is as follows. The first level of subfolders indicates the land use scenario (e.g. BAU, etc). The next level of folders then contain a series of country-level raster maps of agricultural land cover (in the folder named 'tifs'), that are then mosaiced together to create a global map of agricultural land cover (in the folder named 'Global tifs'). The names of the files in the 'tifs' and 'Global tifs' folders are descriptive of the data contained in the rasters. The demo dataset contains only the global tifs due to reduce size of the dataset. This folder will be created when running the 'Land Cover Modelling' scripts.

Land Forecasting Coefs and Plots: This folder contains the regression model coefficients that are used to project agricultural land use under different land use scenarios. This contains three subfolders, which contain (a) level plots showing accuracy of the multinomial models used to select cells for changes in agricultural land cover extent ('Model Plots'); (b) coefficients for multinomial models that select cells for changes in agricultural land cover ('Multinomial Accuracy'); and coefficients for the regression models that estimate the proportion of a cell that would experience an increase/decrease in agricultural land cover ('Model Coefficients'). This folder will be created when running the 'Land Cover Modelling' scripts.

Land Forecasting Targets: This contains a series of csv files that indicate how much more land is projected in each 5 year time period for each country. 

Other Data Inputs: This contains other data inputs that are used for various parts of the analysis, including crop yield forecasts, species richness rasters, country-level raster maps that contain the row index for each raster cell in a country (these are used to mosaic country-level agricultural land use projections back together into a global map), inputs used for the population density models, species-specific habitat preferences, the number of cells in each species' current habitat area, estimated patch dispersal distances for different species, and lookup tables that indicate which countries are contained in different regions.

Outputs: This folder will be created when r running the 'Calculating Biodiversity Scripts'. It contains a series of subfolders that contain (a) outputs for each habitat patch for each species ('Patch_Outputs'), (b) spatial outputs for each species ('Raster_Outputs'), (c) habitat size, population abundance, and population density estimates across each species' entire habitat ('CSV_Outputs'), (d) summary stats for each species (e.g. number of patches, number of populations, total habitat area, etc; 'Managed_CSV_Files' folder), and (e) a series of large CSV files that are compiled versions of the summary stats for each species ('Aggregated_CSV_Files'). There are a series of the 'Aggregated_CSV_Files' folders that correspond with different types of outputs (no migration between habitats, migration between habitats ('_Migration'), and the mitigation scenario ('_EAT_Lancet').

Population Density Estimates: This folder contains a series of csv files that contain correlation coefficients from regression models that are used to estimate population density for different types of species. This folder will be created when running the 'Population Density Estimates' scripts. 

Rarity_Weighted_Rasters: This folder contains a series of rasters that indicate rarity-weighted richness for different types of species. These are not used in this analysis.

Raw Urban 2050 Projections: This contains the raw projections of urban land cover extent to 2050, as from Huang et al 2019. https://iopscience.iop.org/article/10.1088/1748-9326/ab4b71

Scenario_Analysis_Inputs: This contains a series of raster maps that indicate different types of conservation areas. These rasters are used in the mitigation scenario, and specifically to define areas where agriculture cannot expand into, as well as areas where existing agricultural land needs to be removed from. 

Scripts: This contains the scripts used in this analysis. See above for a description of the scripts and the order in which the scripts should be run.

TM_WORLD_BORDERS-0.3: This contains shapefiles that indicate national borders. This folder will be created by running the 'Urbanisation Projection' scripts. 

TMP_SDM_FILE_1: This is a temporary folder that will be created when running the 'SDM Scripts'. Intermediate outputs that are not saved will be written to this folder (and then deleted).

TMP_SDM_FILE_1_MOSAIC: As above with the TMP_SDM_FILE_1 folder. 

Urbanization_Forecasts: This contains projections of urban land cover from 2010 to 2050.  

##### Running the scripts on the Demo Dataset.
To run the scripts on the demo dataset, you should only need to (a) install any R libraries that are not installed on your computer, and (b) change the working directory that is set at the top of each R script so that the working directory refers to the demo dataset.
From there, the scripts need to be run in the order indicated above. Running the scripts will create additional folders that create intermediate outputs (e.g. folders that contain projected habitat range maps for each species) and final outputs (e.g. folders that contain estimated biodiversity trends for each species). 

##### Description of the System Used to Run the Analyses
All analyses were run on Oxford’s Advanced Research Computing (ARC) research cluster. A full description of available on ARC is here (https://arc-user-guide.readthedocs.io/en/latest/arc-systems.html). 
The maximum memory used when running the analyses was 200GB RAM. The space used to store all input and output files for the full analysis exceeds 5TB.  
All analyses were conducted in R Version 4.2.2. Parts of the analysis interfaced Python 3 with R via the R package reticulate. More information on the reticulate package is available here (https://rstudio.github.io/reticulate/).