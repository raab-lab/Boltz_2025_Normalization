# PCA analysis and plots

suppressPackageStartupMessages({
  library(csaw)
  library(plyranges)
  library(tidyverse)
  library(DESeq2)  
  library(sva)
  library(janitor)
  library(readr)
  library(GenomicRanges)
  library(Repitools)
  library(SummarizedExperiment)
  library(ggplot2)
  library(ggrepel)
})

################################################################################

# Load union peaks
k4me3 <- read_bed('/proj/jraablab/users/jraab/menin-mll/data/derived_data/peak_sets/cnr_consensus/k4me3_union.bed') # union_bed made from peaks in all samples

# Load size factors for normalization
combined_data <- read_csv('/proj/jraablab/users/jlboltz/Normalization/combined_data.csv')
sf_barcodes <- as.vector(combined_data$sf_barcodes) # size factors for barcodes
sf_ecoli <- as.vector(combined_data$sf_ecoli) # size factors for e. coli

################################################################################

# Write functions to calculate PCA

# Csaw method for normalization
pca_cnr_csaw <- function(peaks, ab) {
  # get sample_sheet details
  samples <- read.csv('/proj/jraablab/users/jlboltz/Normalization/norm_ss.csv') |> janitor::clean_names()
  samples <- samples |> mutate(sample_number = as.character(sample_number))
  samples <- samples |> mutate(experiment_id = as.character(experiment_id))
  samples <- samples |> dplyr::select(sample_number, experiment_id)
  
  # get bam files
  bams <- list.files(path = '/proj/jraablab/users/jlboltz/Normalization/bams/filtered/', pattern = paste0(ab, '_*.*.bam$'), full.names = T)
  print(bams)
  
  # Sample Sheet
  ss <- data.frame(filename = bams, names = basename(bams))
  ss <- ss |> 
    mutate(short_name = str_replace(names, '.aligned.bam', '') ) |> 
    separate(short_name, into = c('airtable_id', 'NGSID', 'cell_line', 'antibody', 'condition', 'rep'), sep = '_' ) 
  ss$condition <- toupper(ss$condition)
  ss <- ss |> left_join(samples, by = c('airtable_id' = 'sample_number'))
  n_exp <- length(unique(ss$experiment_i_id))
  print(n_exp)
  # CSAW Setup
  param = readParam(pe = 'both', dedup = F, minq = 10) # set this to 10 which matches my bam parameters, higher and we lose some numbers
  peak_counts <- regionCounts(bams, peaks, param = param)
  window_counts <- windowCounts(bams, bin = T, width = 10000) 
  # normalize on the regions we've counted
  
  # Need to set up some sort of sample sheet for comparison
  colData(peak_counts) <- DataFrame(left_join(as.data.frame(colData(peak_counts)), 
                                              ss, by = c('bam.files' ='filename') ) )
  colData(window_counts) <- DataFrame(left_join(as.data.frame(colData(window_counts)), 
                                                ss, by= c('bam.files' = 'filename') ) )
  
  # 
  if (n_exp > 1) {
    dds <- DESeqDataSet(peak_counts, design = ~  experiment_i_id + condition + rep) 
    dds$condition <- relevel(dds$condition, ref = "DMSO")
    norm_dds <- DESeqDataSet(window_counts, design = ~ experiment_i_id + condition + rep)
  }
  else { 
    dds <- DESeqDataSet(peak_counts, design = ~   condition + rep) 
    dds$condition <- relevel(dds$condition, ref = "DMSO")
    norm_dds <- DESeqDataSet(window_counts, design = ~  condition + rep)
  }
  norm_dds <- estimateSizeFactors(norm_dds)
  sf <- sizeFactors(norm_dds)
  sizeFactors(dds) <- sf
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
  pca_data <- plotPCA(vsd, intgroup = c("condition", 'rep'), returnData = TRUE)
  return(pca_data)
}

# Uses spike in method for normalization

pca_cnr_spike_in <- function(peaks, ab, sf) {
  # get sample_sheet details
  samples <- read.csv('/proj/jraablab/users/jlboltz/Normalization/norm_ss.csv') |> janitor::clean_names()
  samples <- samples |> mutate(sample_number = as.character(sample_number))
  samples <- samples |> mutate(experiment_id = as.character(experiment_id))
  samples <- samples |> dplyr::select(sample_number, experiment_id)
  
  # get bam files
  bams <- list.files(path = '/proj/jraablab/users/jlboltz/Normalization/bams/filtered/', pattern = paste0(ab, '_*.*.bam$'), full.names = T)
  print(bams)
  
  # Sample Sheet
  ss <- data.frame(filename = bams, names = basename(bams)) 
  ss <- ss |> 
    mutate(short_name = str_replace(names, '.aligned.bam', '') ) |> 
    separate(short_name, into = c('airtable_id', 'NGSID', 'cell_line', 'antibody', 'condition', 'rep'), sep = '_' ) 
  ss$condition <- toupper(ss$condition)
  ss <- ss |> left_join(samples, by = c('airtable_id' = 'sample_number'))
  n_exp <- length(unique(ss$experiment_i_id))
  print(n_exp)
  # CSAW Setup
  param = readParam(pe = 'both', dedup = F, minq = 10) # set this to 10 which matches my bam parameters, higher and we lose some numbers
  peak_counts <- regionCounts(bams, peaks, param = param)
  # normalize on the regions we've counted
  
  # Need to set up some sort of sample sheet for comparison
  colData(peak_counts) <- DataFrame(left_join(as.data.frame(colData(peak_counts)), 
                                              ss, by = c('bam.files' ='filename') ) )
  
  dds <- DESeqDataSet(peak_counts, design = ~   condition + rep) 
  dds$condition <- relevel(dds$condition, ref = "DMSO")
  
  sizeFactors(dds) <- sf
  dds <- estimateDispersionsGeneEst(dds) 
  dispersions(dds) <- mcols(dds)$dispGeneEst
  dds <- nbinomWaldTest(dds)
  vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
  pca_data <- plotPCA(vsd, intgroup = c("condition", 'rep'), returnData = TRUE)
  return(pca_data)
}

################################################################################

# apply functions using 3 different normalization methods

k4me3_pca_csaw <- pca_cnr_csaw(k4me3, 'H3K4me3')
percentVar_csaw <- round(100 * attr(k4me3_pca_csaw, "percentVar"))

k4me3_pca_barcodes <- pca_cnr_spike_in(k4me3, 'H3K4me3', sf_barcodes) # barcode spike-ins
percentVar_barcodes <- round(100 * attr(k4me3_pca_barcodes, "percentVar"))

k4me3_pca_ecoli <- pca_cnr_spike_in(k4me3, 'H3K4me3', sf_ecoli) # e. coli spike ins
percentVar_ecoli <- round(100 * attr(k4me3_pca_ecoli, "percentVar"))

################################################################################
# Graph PCA data

# csaw
ggplot(k4me3_pca_csaw, aes(PC1, PC2, label = k4me3_pca_csaw$group, color = condition)) +
  geom_point(size = 4) +
  geom_text_repel(max.overlaps = 100) +
  xlab(paste0("PC1: ", percentVar_csaw[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_csaw[2], "% variance")) +
  ggtitle("PCA Plot of Data Normalized with CSAW Method") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# barcodes
ggplot(k4me3_pca_barcodes, aes(PC1, PC2, label = k4me3_pca_barcodes$group, color = condition)) +
  geom_point(size = 4) +
  geom_text_repel(max.overlaps = 100) +
  xlab(paste0("PC1: ", percentVar_barcodes[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_barcodes[2], "% variance")) +
  ggtitle("PCA Plot of Data Normalized with Barcode Spike-In Method") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# e. coli
ggplot(k4me3_pca_ecoli, aes(PC1, PC2, label = k4me3_pca_ecoli$group, color = condition)) +
  geom_point(size = 4) +
  geom_text_repel(max.overlaps = 100) +
  xlab(paste0("PC1: ", percentVar_ecoli[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_ecoli[2], "% variance")) +
  ggtitle("PCA Plot of Data Normalized with Ecoli Spike-In Method") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
