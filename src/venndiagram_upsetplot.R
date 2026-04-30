# Venn diagram plot and upset plot
library(VennDiagram)
library(UpSetR)

################################################################################
# Venn diagram plot overlapping significant peaks
merged <- read_tsv('/proj/jraablab/users/jlboltz/Normalization/merged.tsv')

merged_real <- na.omit(merged) # remove rows with NA
set1 <- merged_real[merged_real$padj.x < 0.01, ] # select the significant peaks
set1 <- set1$start 
set2 <-  merged_real[merged_real$padj.y < 0.01, ]
set2 <- set2$start
set3 <-  merged_real[merged_real$padj < 0.01, ]
set3 <- set3$start

venn.plot <- venn.diagram(
  x = list(
    CSAW = set1,
    Barcodes = set2,
    Ecoli = set3
  ),
  filename = NULL,
  fill = c("lightblue", "lightyellow", "maroon"),
  alpha = 0.5,
  cat.cex = 1.5,
  main = "Significant Peaks Shared Between Normalization Methods"
)
grid::grid.newpage()
grid::grid.draw(venn.plot)

################################################################################
# Upset plot for significant peaks

# Create a matrix for significant genes in each comparison
upset_data <- data.frame(
  Gene = merged_real$start,
  CSAW = as.integer(merged_real$padj.x < 0.01),
  Barcodes = as.integer(merged_real$padj.y < 0.01),
  Ecoli = as.integer(merged_real$padj < 0.01)
)

UpSetR::upset(upset_data, 
      sets = c("CSAW", "Barcodes", "Ecoli"))

