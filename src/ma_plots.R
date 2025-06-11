# MA plot

# load in data

k4me3_res_csaw <- read_tsv('/proj/jraablab/users/jlboltz/Normalization/k4me3_res_csaw.tsv')
k4me3_res_barcodes <- read_tsv('/proj/jraablab/users/jlboltz/Normalization/k4me3_res_barcodes.tsv')
k4me3_res_ecoli <- read_tsv('/proj/jraablab/users/jlboltz/Normalization/k4me3_res_ecoli.tsv')

################################################################################

# CSAW 
k4me3_des_csaw_shrink <- lfcShrink(k4me3_des_csaw, coef = 2, type = 'apeglm') 
plotMA(k4me3_des_csaw_shrink, ylim = c(-5, 5), main = "CSAW MA Plot")

# Barcodes
k4me3_des_barcodes_shrink <- lfcShrink(k4me3_des_barcodes, coef = 2, type = 'apeglm') 
plotMA(k4me3_des_barcodes_shrink, ylim = c(-5, 5), main = "Barcodes MA Plot")

# E. coli
k4me3_des_ecoli_shrink <- lfcShrink(k4me3_des_ecoli, coef = 2, type = 'apeglm') 
plotMA(k4me3_des_ecoli_shrink, ylim = c(-5, 5), main = "Ecoli MA Plot")
