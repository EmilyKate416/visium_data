#!/usr/bin/env Rscript

#' Title: Bulk Visium QC, cohort-specific filtering, and per-sample export
#'
#' Description:
#'   Loads multiple 10x Visium samples into a single SpatialExperiment, performs
#'   annotation, per-spot QC, cohort-specific mitochondrial filtering, and basic
#'   spot/gene filtering. Saves QC summaries/plots, then splits and writes each
#'   sample as an individual RDS for downstream HVG/nnSVG processing.
#'
#' Inputs:
#'   - Visium outputs under `load_all_data` (default: "/mnt/scratchc/fmlab/lythgo02/visium_data/").
#'   - Sample directories whose names start with "SITS".
#'   - Excludes file: "SITSA2.mri.tgz".
#'
#' Key Dependencies (R):
#'   SpatialExperiment, scater, scran, AnnotationHub, ggspavis, tidyverse, DT,
#'   BiocParallel, nnSVG.
#'
#' System Requirements (for package install on cluster if needed):
#'   - BiocManager
#'   - libmagick++-dev (for magick); install via apt: `sudo apt-get update && sudo apt-get install libmagick++-dev`
#'
#' QC & Filtering Logic:
#'   1) Keep only tissue-covered spots (`in_tissue == 1`).
#'   2) Remove genes with zero counts across all spots.
#'   3) Per-spot QC metrics via `addPerCellQC()` with MT genes defined from Ensembl
#'      annotation (AnnotationHub EnsDb v109; "MT" chromosome).
#'   4) Thresholds applied:
#'        - Library size (UMIs): exclude if sum < 300.
#'        - Genes detected: exclude if detected < 300.
#'        - Mitochondrial %:
#'            * Cohort "arbitrary": fixed > 20%.
#'            * Cohort "log": sample-specific robust outliers on log(mito% + 0.1):
#'              flag if value > median + 3 * MAD (upper-tail only).
#'      Cohort membership:
#'        log cohort = c("SITSA3","SITSB2","SITSE4","SITSC3","SITSF4","SITSF2","SITSC1","SITSD3","SITSB4");
#'        all others = "arbitrary".
#'   5) Final spot exclusion if ANY of the above flags are TRUE.
#'   6) Post-filter gene pruning: drop genes detected in < 5 spots.
#'
#' Gene Annotation:
#'   - Adds Ensembl-based chromosome info via AnnotationHub EnsDb (Homo sapiens, v109).
#'   - Row names set to gene symbols where available; otherwise Ensembl IDs retained.
#'
#' Outputs (written to /mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/ and subfolders):
#'   - Intermediate RDS:
#'       * 20250321_post_!rowsWithZeroCounts.rds
#'       * 20250321_post_addpercellQC.rds
#'       * 20250321_spe_postfilter_prelist.rds
#'   - Per-sample RDS:
#'       * 20250321_spe_sub_<SAMPLE>.rds  (after all QC/filters)
#'   - QC plots (PNG):
#'       * spots detected, UMI histograms/log, spot maps for excluded spots,
#'         gene-detected plots, mito% plots.
#'   - Tables (CSV):
#'       * table_report.csv                            (per-sample QC summary)
#'       * combined_stats_arbitrary_report.csv         (arbitrary cohort mito summary)
#'       * combined_stats_arbitrary_report.csv         (log cohort mito summary; file name reused)
#'       * filtersummary1.csv / filtersummary2.csv     (final exclusion tallies per cohort)
#'


library(SpatialExperiment) #in terminal: sudo Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos =                     "https://cran.rstudio.com"); BiocManager::install("SpatialExperiment")'
                                  #if this fails due to magick, first run:  sudo apt-get update
                                                                   # then: sudo apt-get install libmagick++-dev
                                                                   # then sudo R to open R session and install.packages("magick")     
library(scater)
library(AnnotationHub)
library(tidyverse)
library(ggspavis)
library(scran)
library(BiocParallel)
library(DT)
library(nnSVG)
bpp <- MulticoreParam(workers = 16, progressbar = TRUE) #to use with bplapply to spread comp expensive tasks across processors
register(bpp)

bpstart(bpp) #restart cluster
# Check if it is active
bpisup(bpp)

set.seed(123)


projDir = "/mnt/scratchc/fmlab/lythgo02/visium_data/"

load_all_data <- "/mnt/scratchc/fmlab/lythgo02/visium_data/"

# List files in the directory
files <- list.files(path = load_all_data, pattern = "^SITS", full.names = TRUE)

# Subset to exclude this file (which I think is in the wrong location anyway)
files <- files[!grepl("SITSA2\\.mri\\.tgz", files)] 


#assign names so they are loaded into the spe with the correct identifiers
samplesheet <- c("SITSA1", "SITSA2", "SITSA3", "SITSA4", 
                 "SITSB1", "SITSB2", "SITSB3", "SITSB4", 
                 "SITSC1", "SITSC2", "SITSC3", "SITSC4",
                 "SITSD1", "SITSD2", "SITSD3", "SITSD4",
                 "SITSE2", "SITSE4",
                 "SITSF2", "SITSF4",
                 "SITSG2", "SITSH2")

names(files) <- samplesheet

# Initialize a list to store each sample
spe_list <- list()

# Loop through each file and process it
for (sample_name in names(files)) {
  # Read the Visium data
  spe <- read10xVisium(samples = files[sample_name])

  # Add unique barcodes
  spe$barcode <- rownames(colData(spe))
  spe$barcodeid <- gsub("-1$", paste0("-", sample_name), spe$barcode)
  rownames(colData(spe)) <- spe$barcodeid

  # Store the processed sample in the list
  spe_list[[sample_name]] <- spe
}

# Combine all processed samples into one SpatialExperiment object
spe <- do.call(cbind, spe_list)

# This bit gets more annotation than just the gene name that given by spaceranger. For example we need which chromosome its on to check if its mitochondrial.
ah <- AnnotationHub()
# change this to whichever species etc. you require
HumanEnsDb <- query(ah, c("EnsDb", "Homo sapiens", "109"))[[1]]
annotations <- genes(HumanEnsDb, return.type = "data.frame")  %>%
  dplyr::select(gene_id, seq_name) %>%
  dplyr::rename(ID=gene_id)

# adds the new info to a version of the rowdata (gene info) from the spe object
gene_metadata <- as.data.frame(rowData(spe)) %>%
  mutate(ID = rownames(.)) %>%
  left_join(annotations, by = "ID") %>%
  dplyr::rename(Chromosome=seq_name)

# it must be a DataFrame not a data.frame to fit back into the spe object
rowData(spe) <- DataFrame(gene_metadata) 

# this changes the rownames so they are the gene name if we have it but stay as the ensembl id if not
rownames(spe) <- uniquifyFeatureNames(rowData(spe)$ID,
                                      rowData(spe)$symbol)


# This bit just gets rid of any genes that we haven't detected in any spot
rowsWithZeroCounts <- rowSums(counts(spe)) == 0

spe <- spe[!rowsWithZeroCounts, ]

saveRDS(spe, "/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/20250321_post_!rowsWithZeroCounts.rds") #because the steo above can take a long time to run

# subset to keep only spots over tissue
spe <- spe[, colData(spe)$in_tissue == 1]

spots <- plotSpots(spe) +
  facet_wrap(~spe$sample_id) 

ggsave("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/plots/QC/20250321_spotsdetected.png")

# calculate per-spot QC metrics and store in colData (from scRNA approach)
# detected = how many genes in that spot
#sum = how many UMIs in that spot
#subsets_Mito_percent = the percentage of UMIs that come from a mitochondrial gene

spe <- addPerCellQC(spe, 
                    subsets = list(Mito = which(rowData(spe)$Chromosome == "MT")))

# this just makes a table of the stuff for the report
table <- colData(spe) %>%
  as_tibble() %>%
  group_by(sample_id) %>%
  summarize(
    `Number of spots` = n(),
    `Median genes per spot` = round(median(detected)),
    `Median number of UMIs per spot` = round(median(sum)),
    `Total number of UMIs` = sum(sum),
    `Total mito UMIs` = sum(subsets_Mito_sum, na.rm = TRUE),  # Total mito UMIs per sample
    `Mito percent` = (`Total mito UMIs` / `Total number of UMIs`) * 100  # Compute % mito reads
  ) %>%
  ungroup() #%>%
 # datatable(rownames = FALSE, options = list(dom = 't'))

# Write the table to CSV
write.csv(table, "table_report.csv", row.names = FALSE)

saveRDS(spe, "/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/20250321_post_addpercellQC.rds") #because the steo above can take a long time to run


# Define the sample IDs for each cohort
cohort_log <- c("SITSA3", "SITSB2", "SITSE4","SITSC3", "SITSF4", "SITSF2", "SITSC1", "SITSD3", "SITSB4")
cohort_arbitrary <- setdiff(unique(colData(spe)$sample_id), cohort_log)  # Remaining samples

# Split `spe` into two cohorts
spe_cohort1 <- spe[, colData(spe)$sample_id %in% cohort_log]
spe_cohort2 <- spe[, colData(spe)$sample_id %in% cohort_arbitrary]

#library size = total UMI counts per spot (in sum column)
# histogram of library sizes across spots
hist(colData(spe_cohort1)$sum, breaks = 20)
hist(colData(spe_cohort2)$sum, breaks = 20)
#check the distribution is smooth and there are no obvious issues like spike at low lib sizes
# Filter the spe object where the sum is greater than 0

#set a relatively arbitrary threshold of 600 UMI counts per spot (as per tutorial) 
#and then check the number of spots below this thresholdI change this looking at the plots
qc_lib_size1 <- colData(spe_cohort1)$sum < 300
colData(spe_cohort1)$qc_lib_size1 <- colData(spe_cohort1)$sum < 300
table(qc_lib_size1)

qc_lib_size2 <- colData(spe_cohort2)$sum < 300
colData(spe_cohort2)$qc_lib_size2 <- colData(spe_cohort2)$sum < 300
table(qc_lib_size2)

umi <- plotColData(spe_cohort1, y = "sum", x = "sample_id", colour_by = "qc_lib_size1") +
  labs(title = "Total number of UMIs for each spot across all genes",
       y = "Total UMI count (library size)", x = "Sample Name") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/plots/QC/20250321_umi_qc_1.png", plot = umi,width= 8, height = 6, units = "in")

umi <- plotColData(spe_cohort2, y = "sum", x = "sample_id", colour_by = "qc_lib_size2") +
  labs(title = "Total number of UMIs for each spot across all genes",
       y = "Total UMI count (library size)", x = "Sample Name") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/plots/QC/20250321_umi_qc_2.png", plot = umi,width= 8, height = 6, units = "in")


spe_plot1 <- spe_cohort1
colData(spe_plot1)$sum <- colData(spe_plot1)$sum + 0.1

spe_plot2 <- spe_cohort2
colData(spe_plot2)$sum <- colData(spe_plot2)$sum + 0.1

umi_log <- plotColData(spe_plot1, y = "sum", x = "sample_id", colour_by = "qc_lib_size1") +
  labs(
    title = "Total number of UMIs for each spot across all genes",
    y = "Log of total UMI count (library size)", x = "Sample Name"
  ) +
theme(legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_log10()

ggsave("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/plots/QC/20250321_umi_qc_log_1.png", plot = umi_log, width = 8, height = 6, units = "in")

umi_log <- plotColData(spe_plot2, y = "sum", x = "sample_id", colour_by = "qc_lib_size2") +
  labs(
    title = "Total number of UMIs for each spot across all genes",
    y = "Log of total UMI count (library size)", x = "Sample Name"
  ) +
theme(legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_log10()

ggsave("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/plots/QC/20250321_umi_qc_log_2.png", plot = umi_log, width = 8, height = 6, units = "in")

# plots all the umi level info as gradients (useful to see if its spatially affected)
plotSpots(spe_cohort1,  
       annotate = "qc_lib_size1") + facet_wrap(~spe_cohort1$sample_id)

plotSpots(spe_cohort2,  
       annotate = "qc_lib_size2") + facet_wrap(~spe_cohort2$sample_id)
#check that the discarded spots do not have any obvious spatial pattern that correlates 
#with known biological features. Otherwise, removing these spots could indicate that we 
#have set the threshold too high, and are removing biologically informative spots

# plots your final decisions
plotSpots(spe_cohort1, 
       annotate = "qc_lib_size1") + 
        scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
       facet_wrap(~spe_cohort1$sample_id)  
ggsave("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/plots/QC/20250321_spotsexcluded_librarysize_1.png")


# plots your final decisions
plotSpots(spe_cohort2,
       annotate = "qc_lib_size2") + 
        scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey"))+ 
       facet_wrap(~spe_cohort2$sample_id)

ggsave("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/plots/QC/20250321_spotsexcluded_librarysize_2.png")

# histogram of numbers of expressed genes
hist(colData(spe_cohort1)$detected, breaks = 20)

hist(colData(spe_cohort2)$detected, breaks = 20)
#arbitrary threshold of < 300 genes detected
qc_genes_detected1 <- colData(spe_cohort1)$detected < 300
colData(spe_cohort1)$qc_genes_detected1 <- colData(spe_cohort1)$detected < 300

qc_genes_detected2 <- colData(spe_cohort2)$detected < 300
colData(spe_cohort2)$qc_genes_detected2 <- colData(spe_cohort2)$detected < 300
     
plotColData(spe_cohort1, y = "detected", x = "sample_id", colour_by="qc_genes_detected1") +
  labs(y = "Number of detected genes", x = "Sample Name", title = "Genes detected per spot") +
    theme(legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/plots/QC/20250321_col_genes_detected_1.png", width = 8, height = 6, units = "in")

plotSpots(spe_cohort1, annotate = "qc_genes_detected1") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey"))+
  facet_wrap(~spe_cohort1$sample_id)

ggsave("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/plots/QC/20250321_spotexclusion_genesdetected_1.png")
     
plotColData(spe_cohort2, y = "detected", x = "sample_id", colour_by="qc_genes_detected2") +
  labs(y = "Number of detected genes", x = "Sample Name", title = "Genes detected per spot") +
    theme(legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/plots/QC/20250321_col_genes_detected_2.png", width = 8, height = 6, units = "in")

plotSpots(spe_cohort2,  annotate = "qc_genes_detected2") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey"))+
  facet_wrap(~spe_cohort2$sample_id)

ggsave("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/plots/QC/20250321_spotexclusion_genesdetected_2.png")

# its expected to see some coorelation here as the more UMIs you have from a spot in theory the more genes you will sample.
colData(spe_cohort1) %>%
  as_tibble() %>%
  ggplot(mapping = aes(x = detected, y = sum)) +
  geom_point(alpha = 0.4, size = 0.5) +
  labs(
    x = "Number of genes detected",
    y = "Library size"
  ) +
  theme_bw() + 
  facet_wrap(~spe_cohort1$sample_id)

colData(spe_cohort2) %>%
  as_tibble() %>%
  ggplot(mapping = aes(x = detected, y = sum)) +
  geom_point(alpha = 0.4, size = 0.5) +
  labs(
    x = "Number of genes detected",
    y = "Library size"
  ) +
  theme_bw() + 
  facet_wrap(~spe_cohort2$sample_id)

is_mito1 <- grepl("(^MT-)|(^mt-)", rowData(spe_cohort1)$symbol)
table(is_mito1)


mito_outlier_stats_arbitrary1 <- as.data.frame(colData(spe_cohort1)) %>%
  group_by(sample_id) %>%
  summarise(
    total_spots = n(),  # Total spots per sample
    mito_outlier_spots = sum(subsets_Mito_percent > 20, na.rm = TRUE),  # Spots with >20% mito
    mito_outlier_percentage = (mito_outlier_spots / total_spots)*100  # Proportion of outliers
  )

is_mito2 <- grepl("(^MT-)|(^mt-)", rowData(spe_cohort2)$symbol)
table(is_mito2)


mito_outlier_stats_arbitrary2 <- as.data.frame(colData(spe_cohort2)) %>%
  group_by(sample_id) %>%
  summarise(
    total_spots = n(),  # Total spots per sample
    mito_outlier_spots = sum(subsets_Mito_percent > 20, na.rm = TRUE),  # Spots with >20% mito
    mito_outlier_percentage = (mito_outlier_spots / total_spots)*100  # Proportion of outliers
  )

# This plot is more useful in single cell I think because alot of visiums just look like blob
colData(spe_cohort1) %>%
  as_tibble() %>%
  ggplot(mapping = aes(x = detected, y = subsets_Mito_percent)) +
  geom_point(alpha = 0.4, size = 0.5) +
  labs(
    x = "Number of genes detected",
    y = "% UMIs from mitochonial genes"
  ) +
  theme_bw()  +facet_wrap(~spe_cohort1$sample_id)

# This plot is more useful in single cell I think because alot of visiums just look like blob
colData(spe_cohort2) %>%
  as_tibble() %>%
  ggplot(mapping = aes(x = detected, y = subsets_Mito_percent)) +
  geom_point(alpha = 0.4, size = 0.5) +
  labs(
    x = "Number of genes detected",
    y = "% UMIs from mitochonial genes"
  ) +
  theme_bw()  +facet_wrap(~spe_cohort2$sample_id)

#arbitrary
qc_mito <- colData(spe_cohort2)$subsets_Mito_percent > 20
colData(spe_cohort2)$qc_mito <- colData(spe_cohort2)$subsets_Mito_percent >20

plotColData(spe_cohort2, y = "subsets_Mito_percent", x = "sample_id", colour_by = "qc_mito") +
  labs(y = "% UMIs from mitochonial genes", x = "Sample Name",
       title= "Mitochondrial content") +
  theme(legend.position = "none",
   axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/plots/QC/20250321_mito_content_qc_20.png", width = 8, height = 6, units = "in")


plotSpots(spe_cohort2, annotate = "qc_mito") +
 scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  facet_wrap(~spe_cohort2$sample_id) 

ggsave("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/plots/QC/20250321_spotexclusion_mito_20.png")


#sample specific thresholds without log transformation excludes too many


colData(spe_cohort1)$log_mito_percent <- log(colData(spe_cohort1)$subsets_Mito_percent + 0.1)

# Step 2: Compute the sample-wise statistics (median and MAD on log-transformed data)
log_sample_stats <- as.data.frame(colData(spe_cohort1)) %>%
  group_by(sample_id) %>%
  summarise(
    log_median_mito_percent = median(log_mito_percent, na.rm = TRUE),
    log_mad_mito_percent = mad(log_mito_percent, na.rm = TRUE))


# Step 3: Add the median and MAD to `colData(spe)`
colData(spe_cohort1)$log_median_mito_percent <- log_sample_stats$log_median_mito_percent[match(colData(spe_cohort1)$sample_id, log_sample_stats$sample_id)]
colData(spe_cohort1)$log_mad_mito_percent <- log_sample_stats$log_mad_mito_percent[match(colData(spe_cohort1)$sample_id, log_sample_stats$sample_id)]

# Step 4: Calculate the deviation from the median (only for values above the median)
colData(spe_cohort1)$log_deviation_above_median <- colData(spe_cohort1)$log_mito_percent - colData(spe_cohort1)$log_median_mito_percent

# Step 5: Flag outliers as those where the deviation is greater than 3 times the MAD (only for positive deviations)
colData(spe_cohort1)$log_is_outlier_mad <- colData(spe_cohort1)$log_deviation_above_median > (3 * colData(spe_cohort1)$log_mad_mito_percent)
log_is_outlier_mad <- colData(spe_cohort1)$log_deviation_above_median > (3 * colData(spe_cohort1)$log_mad_mito_percent)

plotColData(spe_cohort1, y = "subsets_Mito_percent", x = "sample_id", colour_by = "log_is_outlier_mad") +
  labs(y = "% UMIs from mitochonial genes", x = "Sample Name",
       title= "Mitochondrial content") +
  theme(legend.position = "none",
   axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/plots/QC/20250321_logmitocontent_qc_samplespecific_upperoutliers.png", width = 8, height = 6, units = "in")

plotSpots(spe_cohort1,  annotate = "log_is_outlier_mad") +
 scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  facet_wrap(~spe_cohort1$sample_id)

ggsave("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/plots/QC/20250321_spotexclusion_logmito_samplespecific_upperoutliers.png")

# Combine the non-log and log sample statistics by matching on 'sample_id'


combined_stats_arbitrary <- as.data.frame(colData(spe_cohort2)) %>%
  group_by(sample_id) %>%
  summarise(
    # Arbitrary threshold of 20% mitochondrial content
    mito_outlier_arbitrary_spots = sum(subsets_Mito_percent > 20, na.rm = TRUE),  # Spots with >20% mito
    mito_outlier_arbitrary_percentage = (mito_outlier_arbitrary_spots / n()) * 100,  # Proportion of outliers
    mito_outlier_arbitrary_threshold = 20  # Threshold for arbitrary outlier definition
 )

# View the result

write.csv(combined_stats_arbitrary, "combined_stats_arbitrary_report.csv") 


combined_stats_log <- as.data.frame(colData(spe_cohort1)) %>%
  group_by(sample_id) %>%
  summarise(log_mad_outlier_spots = sum(log_is_outlier_mad == TRUE, na.rm = TRUE),  # Number of log-transformed outliers
    log_mad_outlier_percentage = (log_mad_outlier_spots / n()) * 100,  # Percentage of log-transformed outliers
    log_mad_outlier_threshold = median(log_mito_percent, na.rm = TRUE) + 3 * mad(log_mito_percent, na.rm = TRUE)  # MAD threshold for log-transformed data
  )

write.csv(combined_stats_log, "combined_stats_log_report.csv") 



#sets up logical vector for final filterin, filtered will be TRUE for any spot where either qc_lib_size is TRUE or qc_genes_detected is TRUE (or both).
filtered_arb <- qc_lib_size2 | qc_genes_detected2 | qc_mito  
colData(spe_cohort2)$filtered <- qc_lib_size2 | qc_genes_detected2 | qc_mito
# this makes a table of the final filtering (remember to chage this code where indicated to the filters you want)
filter_summary <- tibble(
  Group = spe_cohort2$sample_id,
  `Number of spots` = 1,
  `Library size filter` = qc_lib_size2, # here
  `Genes detected filter` = qc_genes_detected2, # here
  `Mitochondrial filter` = qc_mito, # here
  `Filtered spots` = filtered_arb,
  `Retained spots` = !filtered_arb
) %>%
  mutate(Group = spe_cohort2$sample_id) %>%
  group_by(Group) %>%
  summarize_all(sum)

filter_summary %>%
  datatable(rownames = FALSE, options = list(dom = 't'))


write.csv(filter_summary, "filtersummary2.csv") 

#sets up logical vector for final filterin, filtered will be TRUE for any spot where either qc_lib_size is TRUE or qc_genes_detected is TRUE (or both).
filtered_log <- qc_lib_size1 | qc_genes_detected1 | log_is_outlier_mad    
colData(spe_cohort1)$filtered <- qc_lib_size1 | qc_genes_detected1 | log_is_outlier_mad
# this makes a table of the final filtering (remember to chage this code where indicated to the filters you want)
filter_summary <- tibble(
  Group = spe_cohort1$sample_id,
  `Number of spots` = 1,
  `Library size filter` = qc_lib_size1, # here
  `Genes detected filter` = qc_genes_detected1, # here
  `Mitochondrial filter` = log_is_outlier_mad, # here
  `Filtered spots` = filtered_log,
  `Retained spots` = !filtered_log
) %>%
  mutate(Group = spe_cohort1$sample_id) %>%
  group_by(Group) %>%
  summarize_all(sum)

filter_summary %>%
  datatable(rownames = FALSE, options = list(dom = 't'))

write.csv(filter_summary, "filtersummary1.csv") 


# plots your final decisions
plotSpots(spe_cohort1, 
       annotate = "filtered") + 
 scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey"))+
 facet_wrap(~spe_cohort1$sample_id)

ggsave("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/plots/QC/20250321_spotslost_cohort1.png")
# plots your final decisions
plotSpots(spe_cohort2, 
       annotate = "filtered") + 
       scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey"))+
       facet_wrap(~spe_cohort2$sample_id)

ggsave("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/plots/QC/20250321_spotslost_cohort2.png")

#actually does the filtering and makes a new object (keep the old one so you can go back if needed)
spe_filt1  <- spe_cohort1[, !filtered_log]
coldata1 <- colData(spe_filt1) %>% as.data.frame() %>% colnames()

colData(spe_filt1)$qc_lib_size <- colData(spe_filt1)$qc_lib_size1
colData(spe_filt1)$qc_genes_detected <- colData(spe_filt1)$qc_genes_detected1
colData(spe_filt1)$qc_mito <- colData(spe_filt1)$log_is_outlier_mad
colData(spe_filt1) <- colData(spe_filt1)[ , !(colnames(colData(spe_filt1)) %in% c("qc_genes_detected1","qc_lib_size1","median_mito_percent","mad_mito_percent","deviation_from_median","is_outlier_mad"           , "log_mito_percent","log_median_mito_percent","log_mad_mito_percent","log_deviation_above_median","log_is_outlier_mad"))]

#actually does the filtering and makes a new object (keep the old one so you can go back if needed)
spe_filt2  <- spe_cohort2[, !filtered_arb]
coldata2 <- colData(spe_filt2) %>% as.data.frame() %>%  colnames()

colData(spe_filt2)$qc_lib_size <- colData(spe_filt2)$qc_lib_size2
colData(spe_filt2)$qc_genes_detected <- colData(spe_filt2)$qc_genes_detected2

colData(spe_filt2) <- colData(spe_filt2)[ , !(colnames(colData(spe_filt2)) %in% c("qc_genes_detected2","qc_lib_size2"))]


spe_filt <- cbind(spe_filt1, spe_filt2)

# here I'm trying to look for genes that aren't 'really' detected, basically so low they aren't useful
number_of_spots_in_which_gene_detected <- rowSums(counts(spe_filt) != 0)


tibble(nspots = number_of_spots_in_which_gene_detected) %>%
  ggplot(aes(x = nspots)) +
  geom_histogram(bins = 50, fill = "grey80", colour = "black") +
  labs(x = "Number of spots in which a gene is detected", y = "Frequency") +
  theme_bw()

# I normally have this at 5 unless i'm dealing with a huge dataset
genes_detected_in_too_few_spots <- number_of_spots_in_which_gene_detected < 5
spe_filt <- spe_filt[!genes_detected_in_too_few_spots, ]
saveRDS(spe_filt, "/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/20250321_spe_postfilter_prelist.rds")


# Extract unique sample names
unique_samples <- unique(colData(spe_filt)$sample_id)

# Loop through each sample and save it separately
for (sample in unique_samples) {
  cat("Saving sample:", sample, "\n")
  
  # Subset data
  spe_sub <- spe_filt[, colData(spe_filt)$sample_id == sample]
  
  # Remove genes and spots with zero expression
  spe_sub <- spe_sub[rowSums(counts(spe_sub)) > 0, ]
  spe_sub <- spe_sub[, colSums(counts(spe_sub)) > 0]
  
  # Save the individual sample
  saveRDS(spe_sub, file = paste0("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/20250321_spe_sub_", sample, ".rds"))
}

cat("All samples saved successfully.\n")




##########################




#check again for any zero rows which may have been introduced 
# Sample-specific filtering and processing in loops

# Extract unique sample names
#unique_samples <- unique(colData(spe_filt)$sample_id)

# Initialize lists to store HVGs, decomposed variance, and processed spe_sub objects
#hvg_list <- list()
#dec_list <- list()
#spe_sub_list <- list()  # Store processed spe_sub objects

# First loop: Process samples for HVG and normalization
#for (sample in unique_samples) {
#  tryCatch({
#    cat("Processing sample:", sample, "\n")

    # Subset the data for the current sample
 #   spe_sub <- spe_filt[, colData(spe_filt)$sample_id == sample]

    # Remove genes with zero expression
 #   ix_zero_genes <- rowSums(counts(spe_sub)) == 0
 #   if (sum(ix_zero_genes) > 0) {
 #     spe_sub <- spe_sub[!ix_zero_genes, ]
 #   }

    # Remove spots with zero expression
 #   ix_zero_spots <- colSums(counts(spe_sub)) == 0
  #  if (sum(ix_zero_spots) > 0) {
   #   spe_sub <- spe_sub[, !ix_zero_spots]
    #}

    # Normalization
    #spe_sub <- computeLibraryFactors(spe_sub)
    #spe_sub <- logNormCounts(spe_sub, BPPARAM = bpp)

    # HVG Calculation
    #dec <- modelGeneVar(spe_sub, assay.type = "logcounts")
    #top_hvgs <- getTopHVGs(dec, prop = 0.1, var.field = "bio")

    # Save results to lists
    #hvg_list[[sample]] <- top_hvgs
    #dec_list[[sample]] <- dec
    #spe_sub_list[[sample]] <- spe_sub  # Store processed spe_sub

# Save the current state for the sample
    #saveRDS(spe_sub, file = paste0("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/20250321_spe_sub_", sample, ".rds"))
    #saveRDS(top_hvgs, file = paste0("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/20250321_hvg_", sample, ".rds"))
    #saveRDS(dec, file = paste0("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/20250321_dec_", sample, ".rds"))


 # }, error = function(e) {
   # cat("Error processing sample:", sample, "\n")
   # cat("Error message:", conditionMessage(e), "\n")
 # })
#}

# Save all processed objects as a single RDS file
#saveRDS(spe_sub_list, "/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/20250321_spe_sub_list.rds")
#saveRDS(hvg_list, "/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/20250321_hvg_results_all_samples.rds")
#saveRDS(dec_list, "/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/20250321_dec_results_all_samples.rds")


# Load the processed spe_sub_list from the first loop
#spe_sub_list <- readRDS("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/20250321_spe_sub_list.rds")

# Initialize a list to store nnSVG results
#svg_list <- list()

# Second loop: nnSVG Calculation
#for (sample in unique_samples) {
 # tryCatch({
 #   cat("Processing nnSVG for sample:", sample, "\n")

    # Load the pre-processed and normalized spe_sub object
   # spe_sub <- spe_sub_list[[sample]]
    
    # Ensure 'logcounts' exists
  #  if (!"logcounts" %in% assayNames(spe_sub)) {
  #    stop(paste("No 'logcounts' found for sample:", sample))
  #  }

    # nnSVG Calculation
   # spe_sub_nnSVG <- nnSVG(spe_sub, assay_name = "logcounts")
   # svg_list[[sample]] <- rowData(spe_sub_nnSVG)

    # Save the updated spe_sub with nnSVG results
   # spe_sub_list[[sample]] <- spe_sub_nnSVG  # Save the updated object

  # Save nnSVG results per sample
    #saveRDS(spe_sub_nnSVG, file = paste0("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/20250321_spe_sub_nnSVG_", sample, ".rds"))
    #saveRDS(rowData(spe_sub_nnSVG), file = paste0("/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/20250321_svg_", sample, ".rds"))
    
 # }, error = function(e) {
 #   cat("Error processing nnSVG for sample:", sample, "\n")
 #   cat("Error message:", conditionMessage(e), "\n")
  #})
#}

# Save updated spe_sub_list and svg_list
#saveRDS(spe_sub_list, "/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/20250321_spe_sub_list_nnSVG.rds")
#saveRDS(svg_list, "/mnt/scratchc/fmlab/lythgo02/visium_data/single_sample_from_filtering/20250321_svg_results_all_samples.rds")
