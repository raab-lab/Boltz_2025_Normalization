## Normalization project
1. CNR pipeline: Ran samples through [Raab Lab CUT&RUN pipeline](https://github.com/raab-lab/cut-n-run)
   (cnr.sh then normalize.sh)
3. Spike-in E. coli counts:  
Counted the number of E. coli spike-in reads in each sample (Snakefile)
4. Spike-in barcode counts:
Counted number of synthetic DNA barcode spike-in reads in each sample (count_spikeint.sh followed by barcode_counts.R)
5. Calculating Size Factors:
Manually calculated the size factors for the barcode spike-ins and E. coli spike-ins (size_factors.R)
6. DESeq2 analysis:
   Differential analysis done using each method of normalization (differential_binding.R)
7. Plot data:
   (log2foldchange_plots.R, ma_plots.R, pca_plots.R, and venndiagram_upsetplot.R)
   
