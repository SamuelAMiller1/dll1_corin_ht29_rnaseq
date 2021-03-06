#!/bin/bash

# Activate Enviornment

export PATH=/N/u/millesaa/Carbonate/miniconda3/envs/rnaseq/bin:$PATH

# Navigate to parent directory

cd /N/slate/millesaa/dll1_corin_ht29_rnaseq/dll1_corin_rnaseq

# Create array for bam files

SBAMS=($(find ./results/aligned/sorted -name "*\.bam"))

# Index the bam files

for SBAM in ${SBAMS[@]}; do
 samtools index $SBAM
done



