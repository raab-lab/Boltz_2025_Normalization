# cnr differential peaks with 3 different normalization methods
# Author: Jesse Raab
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
})
################################################################################

# load union peaks - these are regions of interested for Differential Occupancy

k4me3 <- read_bed('/proj/jraablab/users/jraab/menin-mll/data/derived_data/peak_sets/cnr_consensus/k4me3_union.bed') # union_bed made from peaks in all samples
combined_data <- read_csv('/proj/jraablab/users/jlboltz/Normalization/combined_data.csv')
sf_barcodes <- as.vector(combined_data$sf_barcodes) # size factors for barcodes
sf_ecoli <- as.vector(combined_data$sf_ecoli) # size factors for e. coli

################################################################################

#Write a function that takes a peak set and antibody and returns the DE object

# Uses csaw method for normalization
run_de_cnr_csaw <- function(peaks, ab) {
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
  param = readParam(pe = 'both', dedup = F, minq = 10) # set this to 10 which matches the bam parameters, higher and we lose some numbers
  peak_counts <- regionCounts(bams, peaks, param = param)
  window_counts <- windowCounts(bams, bin = T, width = 10000) 
  # normalize on the regions counted
  
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
  return(dds)
}

# Uses spike in method for normalization, add in 'sf' which is the vector of calculated size factors

run_de_cnr_spike_in <- function(peaks, ab, sf) {
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
  return(dds)
}

################################################################################
# Run functions
# csaw method
k4me3_des_csaw <- run_de_cnr_csaw(k4me3, 'H3K4me3')
resultsNames(k4me3_des_csaw) 
k4me3_res_csaw <- lfcShrink(k4me3_des_csaw, coef = 2, type = 'apeglm', format = 'GRanges')  

# barcode spike-in method
k4me3_des_barcodes <- run_de_cnr_spike_in(k4me3, 'H3K4me3', sf_barcodes)
resultsNames(k4me3_des_barcodes) 
k4me3_res_barcodes <- lfcShrink(k4me3_des_barcodes, coef = 2, type = 'apeglm', format = 'GRanges') 

# e. coli spike-in method
k4me3_des_ecoli <- run_de_cnr_spike_in(k4me3, 'H3K4me3', sf_ecoli)
resultsNames(k4me3_des_ecoli) 
k4me3_res_ecoli <- lfcShrink(k4me3_des_ecoli, coef = 2, type = 'apeglm', format = 'GRanges') 

################################################################################
# Save data

save(k4me3_des_csaw, k4me3_res_csaw, k4me3_des_barcodes, k4me3_res_barcodes, k4me3_des_ecoli, k4me3_res_ecoli, file = '/proj/jraablab/users/jlboltz/Normalization/cnr_peaks.Rda')
write_tsv(k4me3_res_csaw |> as_tibble(), file = '/proj/jraablab/users/jlboltz/Normalization/k4me3_res_csaw.tsv')
write_tsv(k4me3_res_barcodes |> as_tibble(), file = '/proj/jraablab/users/jlboltz/Normalization/k4me3_res_barcodes.tsv')
write_tsv(k4me3_res_ecoli |> as_tibble(), file = '/proj/jraablab/users/jlboltz/Normalization/k4me3_res_ecoli.tsv')
