## Normalization code
Process:
1. CUT&RUN data processing: Ran samples through [Raab Lab CUT&RUN pipeline](https://github.com/raab-lab/cut-n-run)<br>([cnr.sh](https://github.com/raab-lab/Boltz_2025_Normalization/blob/main/src/cnr.sh) then [normalize.sh](https://github.com/raab-lab/Boltz_2025_Normalization/blob/main/src/normalize.sh))
3. Spike-in E. coli counts:  <br>
Counted the number of E. coli spike-in reads in each sample ([Snakefile](https://github.com/raab-lab/Boltz_2025_Normalization/blob/main/src/Snakefile))
4. Spike-in barcode counts:<br>
Counted number of synthetic DNA barcode spike-in reads in each sample ([count_spikeint.sh](https://github.com/raab-lab/Boltz_2025_Normalization/blob/main/src/count_spikeins.sh) followed by [barcode_counts.R](https://github.com/raab-lab/Boltz_2025_Normalization/blob/main/src/barcode_counts.R))
5. Calculating Size Factors:<br>
Manually calculated the size factors for the barcode spike-ins and E. coli spike-ins (size_factors.R)
6. DESeq2 analysis:<br>
   Differential analysis done using each method of normalization (differential_binding.R)
7. Plot data:<br>
   (log2foldchange_plots.R, ma_plots.R, pca_plots.R, and venndiagram_upsetplot.R)
