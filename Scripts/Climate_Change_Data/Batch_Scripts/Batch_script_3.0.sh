#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=24:00:00
#SBATCH --job-name=Reproject_MinTemp
#SBATCH --partition=medium

module purge
module load R/4.1.0-foss-2021a
Rscript /data/pubh-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/Scripts/Climate_Change_Data/3.0_Reprojecting_Rasters.R
