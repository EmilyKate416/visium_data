#!/bin/bash

#SBATCH --error=nnSVG_single.err
#SBATCH --time=0100:00:00             
#SBATCH --cpus-per-task=4      
#SBATCH --mem=24G                   
#SBATCH --partition=epyc     

# Get the sample name as an argument
SAMPLE_NAME=$1


#Rscript /mnt/scratchc/fmlab/lythgo02/visium_data/scripts/1_sample_specific_filtering_wholeworkflow.r

# Run the R script for a single sample
Rscript /mnt/scratchc/fmlab/lythgo02/visium_data/scripts/2_end_stage_hvg_nnsvg_singlesamplejob.r $SAMPLE_NAME




