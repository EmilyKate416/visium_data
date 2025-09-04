#!/usr/bin/env Rscript
# -------------------------------------------------------------------------
# Script: Single-sample HVG + nnSVG processing
#' Description:
#' This script performs downstream processing of a single Visium sample that has
#' already undergone initial QC and filtering in 1_sample_specific_filtering_wholeworkflow.r. 
#'
# Purpose: Run normalization, HVG selection, and nnSVG on one Visium sample
# Input:  sample ID (argument), RDS file from QC step
# Output: RDS files with HVGs, variance decomposition, nnSVG results
# Steps:  - Load pre-filtered sample
#         - Remove zero-count genes/spots
#         - Normalize & log-transform
#         - Identify HVGs
#         - Run nnSVG to detect spatially variable genes
# -------------------------------------------------------------------------

# Load required libraries
library(SpatialExperiment)
library(scater)
library(scran)
library(nnSVG)
#library(BiocParallel)

# Set up parallel processing
#bpp <- MulticoreParam(workers = 16, progressbar = TRUE)
#register(bpp)

set.seed(123)

# Capture the sample ID from command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("No sample ID provided!")
}
sample <- args[1]

tryCatch({
  cat("Processing HVG and nnSVG for sample:", sample, "\n")
  
  # Load the sample
  spe_sub <- readRDS(paste0("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/20250321_spe_sub_", sample, ".rds"))
  
  # Remove genes with zero expression
  ix_zero_genes <- rowSums(counts(spe_sub)) == 0
  if (sum(ix_zero_genes) > 0) {
    spe_sub <- spe_sub[!ix_zero_genes, ]
  }

  # Remove spots with zero expression
  ix_zero_spots <- colSums(counts(spe_sub)) == 0
  if (sum(ix_zero_spots) > 0) {
    spe_sub <- spe_sub[, !ix_zero_spots]
  }

  # Use nnSVG function to filter lowly expressed genes
  spe_sub <- filter_genes(spe_sub, filter_mito=FALSE)

  # Normalize counts
  spe_sub <- computeLibraryFactors(spe_sub)
  spe_sub <- logNormCounts(spe_sub, BPPARAM = bpp)

  # HVG Calculation
  dec <- modelGeneVar(spe_sub, assay.type = "logcounts")
  top_hvgs <- getTopHVGs(dec, prop = 0.1, var.field = "bio")

  rowData(spe_sub)$hvg <-rownames(spe_sub)  %in% top_hvgs

  # Save intermediate results
  saveRDS(spe_sub, file = paste0("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/20250321_spe_sub_modelgenevar_", sample, ".rds"))
  #saveRDS(top_hvgs, file = paste0("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/20250321_hvg_", sample, ".rds"))
  saveRDS(dec, file = paste0("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/20250321_dec_", sample, ".rds"))

  # Ensure 'logcounts' exists before running nnSVG
  if (!"logcounts" %in% assayNames(spe_sub)) {
    stop(paste("No 'logcounts' found for sample:", sample))
  }

  # Run nnSVG
  spe_sub_nnSVG <- nnSVG(spe_sub, assay_name = "logcounts")

  # Save the final spe_sub with nnSVG results
  saveRDS(spe_sub_nnSVG, file = paste0("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/20250321_spe_sub_nnSVG_", sample, ".rds"))

 # cat("Successfully processed and saved sample:", sample, "\n")

}, error = function(e) {
  cat("Error processing sample:", sample, "\n")
  cat("Error message:", conditionMessage(e), "\n")
})
