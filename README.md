## Normalization project
1. CNR pipeline: cnr.sh
2. CNR normalization: normalization.sh
3. Spike-in: ecoli counts
4. Spike-in: barcode counts
5. Calculating Size Factors
6. DESeq2 analysis
7. Plots
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
    ├── atac
    │   └── clean_atac.R
    ├── cnr
    │   └── clean_cnr.R
    ├── create_canonical_sheet.R
    ├── integrate
    │   └── integrate.R
    └── rna
        └── clean_rna.R
```
