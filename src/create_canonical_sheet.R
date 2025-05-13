## This example script reads in the samplesheet compiled from airtable
## and adds the group_norm, group_avg, and params column to be used for the 'canonical' run of the pipeline
library(dplyr)
library(magrittr)
sheet <- read.csv("data/source_data/samplesheet.csv", check.names = F)

## EXAMPLE
full_sheet <- sheet %>%
  mutate(group_norm = Antibody,
         group_avg = paste0(Antibody, "_", Treatment)) %>%
  mutate(params = case_when(
    Antibody == 'menin' ~ '--bins --binSize=10000 --threshold=3',
    Antibody == 'MLL1' ~ '--peaks --winSize=50 --threshold=3',
    Antibody == 'H3K4me3' ~ '--peaks --winSize=50 --threshold=5',
    Antibody == 'H3K27me3' ~ '--peaks --winSize=100 --threshold=4',
    Antibody == 'IgG' ~ '--peaks'
  ))
write.csv(full_sheet, "data/source_data/canonical_samplesheet.csv", quote = F, row.names = F)
