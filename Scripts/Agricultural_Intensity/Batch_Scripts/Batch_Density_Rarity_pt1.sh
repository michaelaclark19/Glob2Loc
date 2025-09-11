#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem-per-cpu=30G
#SBATCH --time=23:59:00
#SBATCH --job-name=land_s
#SBATCH --partition=medium
#SBATCH --account=ouce-glob2loc

module purge
module load R/4.2.2-foss-2022a-ARC
Rscript /data/ouce-glob2loc/pubh0329/Multiple_Stresses_of_Biodiversity/Scripts/Agricultural\ Intensity/3.0_Classying_Rarity_Pop_Density_6Sep2023_pt3.R
