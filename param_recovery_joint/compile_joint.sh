#!/bin/bash
#SBATCH --job-name=compile_joint
#SBATCH --partition=Calcul
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=./output/compile_joint.log
#SBATCH --error=./error/compile_joint.err
#SBATCH --mail-user=thibault.gajdos@univ-amu.fr
#SBATCH --mail-type=END,FAIL

## Step 1: compile the model once
Rscript R/compile_joint.r
