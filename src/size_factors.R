# Calculate size factors from ecoli and barcode spike-ins

suppressPackageStartupMessages({
  library(tidyverse)
})

################################################################################

# load e. coli and hlf alignments

align_ecoli <- list.files(path = '/proj/jraablab/users/jlboltz/Normalization/ecoli_alignments/', pattern = paste0('ecoli.alignment_stats'), full.names = T)
align_hlf <- list.files(path = '/proj/jraablab/users/jlboltz/Normalization/bams', pattern = paste0('.alignment_stats'), full.names = T)

load_csv <- function(bams, column_name) {
  df <- data.frame()
  for (bam in bams){
    csv <- read_lines(bam)
    line_4 <- csv[4]
    line_4_new <- sub(".*?(\\d+).*", "\\1", line_4)
    df <- rbind(df, line_4_new)}
  names(df)[1] <- column_name
  print(df)
  return(df)
}

ecoli <- load_csv(align_ecoli, 'Ecoli')
hlf <- load_csv(align_hlf, 'HLF')
barcodes <- read_csv('/proj/jraablab/users/jlboltz/Normalization/barcodes.csv')
  
# combine into one data frame

combined_data <- cbind(barcodes, hlf, ecoli)
combined_data$HLF <- as.numeric(combined_data$HLF)
combined_data$Ecoli <- as.numeric(combined_data$Ecoli)
combined_data <- combined_data[c(-1, -4), ] # first and fourth sample had low counts and were not used
combined_data <- combined_data[order(combined_data$NGSID), ]

################################################################################

# calculate size factors

combined_data$sf_barcodes <- ((combined_data$barcode_counts)/(combined_data$barcode_counts + combined_data$HLF + combined_data$Ecoli)) * 1000000
combined_data$sf_ecoli <- ((combined_data$Ecoli)/(combined_data$barcode_counts + combined_data$HLF + combined_data$Ecoli)) * 1000000

write_csv(combined_data |> as_tibble(), file = '/proj/jraablab/users/jlboltz/Normalization/combined_data.csv')

