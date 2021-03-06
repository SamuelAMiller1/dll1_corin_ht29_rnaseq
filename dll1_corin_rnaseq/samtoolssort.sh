#!/bin/bash

# Activate Enviornment

export PATH=/N/u/millesaa/Carbonate/miniconda3/envs/rnaseq/bin:$PATH

# Navigate to parent directory

cd /N/slate/millesaa/dll1_corin_ht29_rnaseq/dll1_corin_rnaseq

# Make directory for sorted files

mkdir -p ./results/aligned/sorted/

# Create array for bam files

BAMS=($(find ./results/aligned -name "*\.bam"))

# Note: STAR sorts reads by coordinates, for PE fragments, the input bam files need to be re-sorted to take into account which read is aligning prior to indexing

# Sort the BAMs

for BAM in ${BAMS[@]}; do
  samtools sort -@ 7 -o results/aligned/sorted/$(basename $BAM .bam)_sorted.bam -O bam $BAM
done

