## Implements qc from the epicypher excel sheet to get barcode counts for a specific antibody
library(magrittr)
library(dplyr)
library(tidyr)
library(readr)

read_spikeins <- function(r1, r2) {
  r1_reads <- read.table(r1, header = F,
                         col.names = c("barcode", "fastq", "nuc_reads")) %>%
    mutate(NGSID = gsub("_.*", "", basename(fastq)))
  
  r2_reads <- read.table(r2, header = F,
                         col.names = c("barcode", "fastq", "nuc_reads")) %>%
    mutate(NGSID = gsub("_.*", "", basename(fastq)))
  spikein_reads <- inner_join(
    x = r1_reads,
    y = r2_reads,
    by = c('barcode', 'NGSID'),
    suffix = c('_r1', '_r2')
  )
  return(spikein_reads)
}

aligned_reads <- read.table("/proj/jraablab/users/jlboltz/Normalization/multiqc_data/multiqc_bowtie2.txt",
                            header = T) %>%
  mutate(NGSID = gsub("[0-9]+_(JR[0-9]{3})_.*", "\\1", Sample)) %>%
  dplyr::select(NGSID, total_reads, paired_aligned_one)
barcodes <- read.table("/proj/jraablab/users/pkuhlers/seq_resources/spikein_barcodes.txt")
spikein_reads <-
  read_spikeins(
    "/proj/jraablab/users/jlboltz/Normalization/R1_spikein_counts.txt",
    "/proj/jraablab/users/jlboltz/Normalization/R2_spikein_counts.txt"
  )
barcode_identities <- read.csv("/proj/jraablab/users/pkuhlers/seq_resources/spikein_barcode_identities.csv")

barcode_reads <- inner_join(spikein_reads,
                            barcode_identities,
                            by = c('barcode' = 'sequence')) %>%
  dplyr::select(NGSID, sequence = barcode, barcode = barcode.y, target, nuc_reads_r1, nuc_reads_r2) %>%
  mutate(tot_reads = nuc_reads_r1 + nuc_reads_r2)

on_target <- barcode_reads %>%
  group_by(NGSID, target) %>%
  reframe(barcode_counts = sum(tot_reads)) %>%
  group_by(NGSID) %>%
  mutate(on_target_norm = barcode_counts / barcode_counts[target == "H3K4me3"] * 100 )

barcodes <- on_target[on_target$target == 'H3K4me3', ]

barcodes <- subset(barcodes, select = c(NGSID, barcode_counts))

barcodes

write_csv(barcodes |> as_tibble(), file = '/proj/jraablab/users/jlboltz/Normalization/barcodes.csv')