# Venn diagram plot and upset plot
library(VennDiagram)
library(UpSetR)
library(tidyverse)
library(grid)

################################################################################
# Venn diagram plot of negative significant peaks
merged_real <- read_tsv(
  '/proj/jraablab/users/jlboltz/Normalization/merged_id.tsv',
  show_col_types = FALSE
)

cleaned_data <- merged_real %>%
  select(ID, padj.x, padj.y, padj,
         log2FoldChange.x, log2FoldChange.y, log2FoldChange) %>%
  mutate(
    # Ensure numeric
    across(starts_with("padj"), as.numeric),
    across(starts_with("log2FoldChange"), as.numeric),
    
    # Replace NA values
    across(starts_with("padj"), ~replace_na(.x, 1)),
    across(starts_with("log2FoldChange"), ~replace_na(.x, 0))
  ) %>%
  distinct(ID, .keep_all = TRUE)

set1 <- cleaned_data %>%
  filter(padj.x < 0.05 & log2FoldChange.x < 0) %>%
  pull(ID)

set2 <- cleaned_data %>%
  filter(padj.y < 0.05 & log2FoldChange.y < 0) %>%
  pull(ID)

set3 <- cleaned_data %>%
  filter(padj < 0.05 & log2FoldChange < 0) %>%
  pull(ID)

png(
  filename = "Negative_Significant_Peaks_Venn.png",
  width = 2000,
  height = 2000,
  res = 300
)

venn.plot <- venn.diagram(
  x = list(
    CSAW = set1,
    Barcodes = set2,
    Ecoli = set3
  ),
  filename = NULL,
  fill = c("#DDCC77", "#EE7733", "#44AA99"),
  alpha = 0.5,
  cat.cex = 2,
  cex = 1.5,
  main = "Significant Negative Peaks",
  main.cex = 2,
  cat.fontfamily = "Helvetica",
  fontfamily = "Helvetica",
  main.fontfamily = "Helvetica",
  cat.pos = c(-20, 20, 180),
  cat.dist = c(0.06, 0.06, 0.06)
)

grid.newpage()
grid.draw(venn.plot)
dev.off()

################################################################################
# Venn diagram plot of positive significant peaks
merged_real <- read_tsv(
  '/proj/jraablab/users/jlboltz/Normalization/merged_id.tsv',
  show_col_types = FALSE
)

cleaned_data <- merged_real %>%
  select(ID, padj.x, padj.y, padj,
         log2FoldChange.x, log2FoldChange.y, log2FoldChange) %>%
  mutate(
    # Ensure numeric
    across(starts_with("padj"), as.numeric),
    across(starts_with("log2FoldChange"), as.numeric),
    
    # Replace NA values
    across(starts_with("padj"), ~replace_na(.x, 1)),
    across(starts_with("log2FoldChange"), ~replace_na(.x, 0))
  ) %>%
  distinct(ID, .keep_all = TRUE)

set1 <- cleaned_data %>%
  filter(padj.x < 0.05 & log2FoldChange.x > 0) %>%
  pull(ID)

set2 <- cleaned_data %>%
  filter(padj.y < 0.05 & log2FoldChange.y > 0) %>%
  pull(ID)

set3 <- cleaned_data %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>%
  pull(ID)

png(
  filename = "Positive_Significant_Peaks_Venn.png",
  width = 2000,
  height = 2000,
  res = 300
)

venn.plot <- venn.diagram(
  x = list(
    CSAW = set1,
    Barcodes = set2,
    Ecoli = set3
  ),
  filename = NULL,
  fill = c("#DDCC77", "#EE7733", "#44AA99"),
  alpha = 0.5,
  cat.cex = 2,
  cex = 1.5,
  main = "Significant Positive Peaks",
  main.cex = 2,
  cat.fontfamily = "Helvetica",
  fontfamily = "Helvetica",
  main.fontfamily = "Helvetica",
  cat.pos = c(-20, 20, 180),
  cat.dist = c(0.06, 0.06, 0.06)
)

grid.newpage()
grid.draw(venn.plot)
dev.off()

################################################################################
# Upset plot for negative significant peaks

# Load data
merged_real <- read_tsv(
  '/proj/jraablab/users/jlboltz/Normalization/merged_id.tsv',
  show_col_types = FALSE
)

# Clean and prepare data
cleaned_data <- merged_real %>%
  select(ID, padj.x, padj.y, padj,
         log2FoldChange.x, log2FoldChange.y, log2FoldChange) %>%
  mutate(
    # Ensure numeric
    across(starts_with("padj"), as.numeric),
    across(starts_with("log2FoldChange"), as.numeric),
    
    # Handle missing values
    across(starts_with("padj"), ~replace_na(.x, 1)),
    across(starts_with("log2FoldChange"), ~replace_na(.x, 0))
  ) %>%
  distinct(ID, .keep_all = TRUE)

# Create UpSet input (ONLY binary columns)
upset_data <- cleaned_data %>%
  transmute(
    CSAW = as.integer(padj.x < 0.05 & log2FoldChange.x < 0),
    Barcodes = as.integer(padj.y < 0.05 & log2FoldChange.y < 0),
    Ecoli = as.integer(padj < 0.05 & log2FoldChange < 0)
  ) %>%
  as.data.frame()

# Sanity check (should all be 0)
print(colSums(is.na(upset_data)))

# Save plot
png(
  filename = "Negative_Upset_Significant_Peaks.png",
  width = 2000,
  height = 1500,
  res = 300
)

UpSetR::upset(
  upset_data,
  text.scale = c(2.5, 2, 2.5, 1.5, 3)
)

dev.off()

################################################################################
# Upset plot for positive significant peaks

# Load data
merged_real <- read_tsv(
  '/proj/jraablab/users/jlboltz/Normalization/merged_id.tsv',
  show_col_types = FALSE
)

# Clean and prepare data
cleaned_data <- merged_real %>%
  select(ID, padj.x, padj.y, padj,
         log2FoldChange.x, log2FoldChange.y, log2FoldChange) %>%
  mutate(
    # Ensure numeric
    across(starts_with("padj"), as.numeric),
    across(starts_with("log2FoldChange"), as.numeric),
    
    # Handle missing values
    across(starts_with("padj"), ~replace_na(.x, 1)),
    across(starts_with("log2FoldChange"), ~replace_na(.x, 0))
  ) %>%
  distinct(ID, .keep_all = TRUE)

# Create UpSet input (ONLY binary columns)
upset_data <- cleaned_data %>%
  transmute(
    CSAW = as.integer(padj.x < 0.05 & log2FoldChange.x > 0),
    Barcodes = as.integer(padj.y < 0.05 & log2FoldChange.y > 0),
    Ecoli = as.integer(padj < 0.05 & log2FoldChange > 0)
  ) %>%
  as.data.frame()

# Sanity check (should all be 0)
print(colSums(is.na(upset_data)))

# Save plot
png(
  filename = "Positive_Upset_Significant_Peaks.png",
  width = 2000,
  height = 1500,
  res = 300
)

UpSetR::upset(
  upset_data,
  text.scale = c(2.5, 2, 2.5, 1.5, 3)
)
