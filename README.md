## Normalization project
1. CNR pipeline: 
Ran samples through Raab Lab CUT&RUN pipeline (cnr.sh)
2. Spike-in: E. coli counts
Counted the number of E. coli spike-in reads in each sample (Snakefile)
3. Spike-in: barcode counts

Counted number of synthetic DNA barcode spike-in reads in each sample (count_spikeint.sh followed by barcode_counts.R)
5. Calculating Size Factors:

Manually calculated the size factors for the barcode spike-ins and E. coli spike-ins (size_factors.R)
7. DESeq2 analysis:
   
Differential analysis done using each method of normalization (differential_binding.R)
9. Plot data:

(log2foldchange_plots.R, ma_plots.R, pca_plots.R, and venndiagram_upsetplot.R)
```
project-template/
├── data
│   ├── derived_data		<- 'Final' analysis outputs, like differential binding/expression tables.
│   │   ├── cnr
│   ├── external		<- Location for potential public or non-pipeline processed data
│   ├── processed_data		<- Intermediate datasets, like summarizedExperiment RDS files or consensus peak calls
│   │   ├── cnr
│   └── source_data		<- Here is where you should link in data from pipeline runs
│       ├── cnr
├── figures			<- Any graphical outputs -- boxplots, scatterplots, etc
│   ├── cnr
├── project-template.Rproj	<- Rename this to your project name
├── README.md
├── renv
├── renv.lock			<- For use with the `renv` package. Helps everyone stay in the same R environment
├── Snakefile			<- See above
└── src				<- All code necessary to generate figures and derived data
    ├── cnr
    │   └── clean_cnr.R
    ├── create_canonical_sheet.R
```
