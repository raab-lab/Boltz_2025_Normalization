# MA plot
library(DESeq2)

# Load data

k4me3_res_csaw <- read_tsv('/proj/jraablab/users/jlboltz/Normalization/k4me3_res_csaw.tsv')
k4me3_res_barcodes <- read_tsv('/proj/jraablab/users/jlboltz/Normalization/k4me3_res_barcodes.tsv')
k4me3_res_ecoli <- read_tsv('/proj/jraablab/users/jlboltz/Normalization/k4me3_res_ecoli.tsv')

################################################################################

# Standard MA plots

# CSAW 
k4me3_des_csaw_shrink <- lfcShrink(k4me3_des_csaw, coef = 2, type = 'apeglm') 
plotMA(k4me3_des_csaw_shrink, ylim = c(-5, 5), main = "CSAW MA Plot")

# Barcodes
k4me3_des_barcodes_shrink <- lfcShrink(k4me3_des_barcodes, coef = 2, type = 'apeglm') 
plotMA(k4me3_des_barcodes_shrink, ylim = c(-5, 5), main = "Barcodes MA Plot")

# E. coli
k4me3_des_ecoli_shrink <- lfcShrink(k4me3_des_ecoli, coef = 2, type = 'apeglm') 
plotMA(k4me3_des_ecoli_shrink, ylim = c(-5, 5), main = "Ecoli MA Plot")

################################################################################

# Alternative MA plot that highlights significant upregulated peaks in red and downregulated in blue

res_df <- as.data.frame(k4me3_des_csaw_shrink)
res_df$regulation <- "Not significant"
res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange > 0] <- "Upregulated"
res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange < 0] <- "Downregulated"

# Plot using ggplot2
ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = regulation)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "gray")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),               
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    axis.line = element_line(color = "black")   
  ) +scale_y_continuous(limits = c(-10, 10)) +
  scale_x_log10(breaks = c(1, 100, 10000)) + # Log scale for x-axis
  labs(x = "Mean expression", y = "Log2 Fold Change", color = "Regulation")
