#!/bin/bash

#SBATCH -c 2

#SBATCH --mem=10G

#SBATCH -t 24:00:00

#SBATCH -J NF

module add nextflow

## MODIFY THE SAMPLESHEET AND EMAIL FOR YOUR RUN
nextflow run raab-lab/cut-n-run \
		--sample_sheet /proj/jraablab/users/jlboltz/Normalization/norm_ss.csv \
		-profile hg38 \ # change genome as needed
		-w /work/users/j/l/jlboltz \
		--outdir /proj/jraablab/users/jlboltz/Normalization \
		-with-report \
		-N jlboltz@unc.edu \ # emails a report
		-latest \
		-ansi-log false \
		-resume
