
# Spatial QC 

## Purpose
This repository contains R and Python workflows for spatial transcriptomics (Visium) data:
- **Quality control (QC)** of raw Visium outputs
- **Normalization & HVG selection** for downstream analysis
- **Spatial gene analysis** with nnSVG

---

## Repo Structure

- `20240206_QC_rstudio_version.Rmd` - used to figure out QC and downstream steps of spatial data before processing in ipynb/via cluster
- `R_spatial_instructions.ipynb` - Main master instructions for various options of QC and post-processing of spatial data - go with option #2 
  - #0 manual processing in interactive R session on cluster 
  - #1 to process each sample individually via slurm submission - end up with different dimensions so difficult to integrate later
  - **#2 to process in bulk up to or including nnsvg step  ** - preferred
    - `1_sample_specific_filtering_wholeworkflow.r` -   Loads multiple 10x Visium samples for QC together then splits and writes each sample to individual RDS for downstream HVG/nnSVG processing.
    - `2_end_stage_hvg_nnsvg_singlesamplejob.r` - End section of the workflow to submit single jobs for each sample for hvg and nnsvg.
    - `submit_spatial_scripts.sh` - used by #2 to submit 1_sample_specific_filtering_wholeworkflow.r and 2_end_stage_hvg_nnsvg_singlesamplejob.r
---

