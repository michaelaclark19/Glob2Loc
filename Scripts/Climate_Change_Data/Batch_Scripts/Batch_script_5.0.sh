#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem-per-cpu=20G
#SBATCH --time=24:00:00
#SBATCH --job-name=Reproject_Climate
#SBATCH --partition=medium

module purge
module load R/4.1.0-foss-2021a
Rscript /data/pubh-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/Scripts/Climate_Change_Data/5.0_Water_Balance.R
