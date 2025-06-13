# Plot Log2FoldChange between samples

library(readr)
library(ggplot2)
library(ggrepel)

################################################################################
# load files and merge them
 
k4me3_res_csaw <- read_tsv('/proj/jraablab/users/jlboltz/Normalization/k4me3_res_csaw.tsv')
k4me3_res_barcodes <- read_tsv('/proj/jraablab/users/jlboltz/Normalization/k4me3_res_barcodes.tsv')
k4me3_res_ecoli <- read_tsv('/proj/jraablab/users/jlboltz/Normalization/k4me3_res_ecoli.tsv')

merged <- merge(x= as.data.frame(k4me3_res_csaw), y = as.data.frame(k4me3_res_barcodes), by = "start")
merged <- merge(x= as.data.frame(merged), y = as.data.frame(k4me3_res_ecoli), by = "start")

write_tsv(merged |> as_tibble(), file = '/proj/jraablab/users/jlboltz/Normalization/merged.tsv')

################################################################################
# plot merged data

# barcodes vs. csaw method
graph <- ggplot(merged, aes(x = log2FoldChange.x, y = log2FoldChange.y, color = padj.x < 0.05 & padj.y < 0.05)) +
  geom_point() +
  xlab("log2FoldChange of csaw normalization") +
  ylab("log2FoldChange of barcode spike-in normalization") +
  scale_color_manual(values = c("black", "red")) +  # Significant adjusted p-values are red
  ggtitle("log2FoldChange barcode vs. csaw methods") +
  theme_minimal()
graph

# e. coli vs. csaw method
graph <- ggplot(merged, aes(x = log2FoldChange.x, y = log2FoldChange, color = padj.x < 0.05 & padj < 0.05)) +
  geom_point() +
  xlab("log2FoldChange of csaw normalization") +
  ylab("log2FoldChange of ecoli spike-in normalization") +
  scale_color_manual(values = c("black", "red")) +  # Significant adjusted p-values are red
  ggtitle("log2FoldChange ecoli vs. csaw methods") +
  theme_minimal()
graph

# e. coli vs. barcodes method
graph <- ggplot(merged, aes(x = log2FoldChange.y, y = log2FoldChange, color = padj.y < 0.05 & padj < 0.05)) +
  geom_point() +
  xlab("log2FoldChange of barcode spike-in normalization") +
  ylab("log2FoldChange of ecoli spike-in normalization") +
  scale_color_manual(values = c("black", "red")) +  # Significant adjusted p-values are red
  ggtitle("log2FoldChange ecoli vs. barcodes methods") +
  theme_minimal()
graph
